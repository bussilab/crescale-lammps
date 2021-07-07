/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_press_crescale.h"
#include <cstring>
#include <cmath>

#include "atom.h"
#include "force.h"
#include "comm.h"
#include "modify.h"
#include "fix_deform.h"
#include "compute.h"
#include "kspace.h"
#include "random_mars.h"
#include "update.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define TILTMAX 1.5

enum{NOBIAS,BIAS};
enum{NONE,XYZ,XY,YZ,XZ};
enum{ISO,ANISO,TRICLINIC};

/* ---------------------------------------------------------------------- */

FixPressCRescale::FixPressCRescale(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  id_temp(nullptr), id_press(nullptr), random(nullptr), tflag(0), pflag(0)
{
  if (narg < 5) error->all(FLERR,"Illegal fix press/crescale command");

  // CRescale barostat applied every step

  nevery = 1;

  // default values

  pcouple = NONE;
  bulkmodulus = 10.0;
  allremap = 1;

  for (int i = 0; i < 6; i++) {
    p_start[i] = p_stop[i] = p_period[i] = 0.0;
    p_flag[i] = 0;
    p_period[i] = 0.0;
  }

  // set fixed-point to default = center of cell

  fixedpoint[0] = 0.5*(domain->boxlo[0]+domain->boxhi[0]);
  fixedpoint[1] = 0.5*(domain->boxlo[1]+domain->boxhi[1]);
  fixedpoint[2] = 0.5*(domain->boxlo[2]+domain->boxhi[2]);

  // process keywords

  dimension = domain->dimension;

  int iarg = 3;
  int seed = 1998;
  t_start = t_target = t_stop = -1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"iso") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      pcouple = XYZ;
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[0] = p_period[1] = p_period[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;
    } else if (strcmp(arg[iarg],"aniso") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      pcouple = NONE;
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[0] = p_period[1] = p_period[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[0] = p_flag[1] = p_flag[2] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
      }
      iarg += 4;
    } else if (strcmp(arg[iarg],"tri") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      pcouple = NONE;
      p_start[0] = p_start[1] = p_start[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[0] = p_stop[1] = p_stop[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_start[3] = p_start[4] = p_start[5] = 0.0;
      p_stop[3] = p_stop[4] = p_stop[5] = 0.0;
      p_period[0] = p_period[1] = p_period[2] 
        = p_period[3] = p_period[4] = p_period[5] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[0] = p_flag[1] = p_flag[2]  
        = p_flag[3] = p_flag[4] = p_flag[5] = 1;
      if (dimension == 2) {
        p_start[2] = p_stop[2] = p_period[2] = 0.0;
        p_flag[2] = 0;
        p_start[3] = p_stop[3] = p_period[3] = 0.0;
        p_flag[3] = 0;
        p_start[4] = p_stop[4] = p_period[4] = 0.0;
        p_flag[4] = 0;
      }
      iarg += 4;

    } else if (strcmp(arg[iarg],"x") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      p_start[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[0] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[0] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[0] = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"y") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      p_start[1] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[1] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[1] = 1;
      iarg += 4;
    } else if (strcmp(arg[iarg],"z") == 0) {
      if (iarg+4 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      p_start[2] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[2] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[2] = 1;
      iarg += 4;
      if (dimension == 2)
        error->all(FLERR,"Invalid fix press/crescale for a 2d simulation");

    } else if (strcmp(arg[iarg],"yz") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix press/crescale command");
      p_start[3] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[3] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[3] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[3] = 1;
      iarg += 4;
      if (dimension == 2)
        error->all(FLERR,"Invalid fix press/crescale command for a 2d simulation");
    } else if (strcmp(arg[iarg],"xz") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal press/crescale command");
      p_start[4] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[4] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[4] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[4] = 1;
      iarg += 4;
      if (dimension == 2)
        error->all(FLERR,"Invalid fix press/crescale command for a 2d simulation");
    } else if (strcmp(arg[iarg],"xy") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix press/crescale command");
      p_start[5] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      p_stop[5] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      p_period[5] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      p_flag[5] = 1;
      iarg += 4;

    } else if (strcmp(arg[iarg],"couple") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      if (strcmp(arg[iarg+1],"xyz") == 0) pcouple = XYZ;
      else if (strcmp(arg[iarg+1],"xy") == 0) pcouple = XY;
      else if (strcmp(arg[iarg+1],"yz") == 0) pcouple = YZ;
      else if (strcmp(arg[iarg+1],"xz") == 0) pcouple = XZ;
      else if (strcmp(arg[iarg+1],"none") == 0) pcouple = NONE;
      else error->all(FLERR,"Illegal fix press/crescale command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix press/crescale command");
      t_start = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      t_target = t_start;
      t_stop = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (t_start <= 0.0 || t_stop <= 0.0)
        error->all(FLERR,"Target temperature for fix press/crescale cannot be 0.0");
      iarg += 3;
    } else if (strcmp(arg[iarg],"modulus") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      bulkmodulus = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (bulkmodulus <= 0.0)
        error->all(FLERR,"Illegal fix press/crescale command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"dilate") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      if (strcmp(arg[iarg+1],"all") == 0) allremap = 1;
      else if (strcmp(arg[iarg+1],"partial") == 0) allremap = 0;
      else error->all(FLERR,"Illegal fix press/crescale command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"seed") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal fix press/crescale command");
      seed = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (seed <= 0) error->all(FLERR,"Illegal fix temp/csvr command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"fixedpoint") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix press/crescale command");
      fixedpoint[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      fixedpoint[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      fixedpoint[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else error->all(FLERR,"Illegal fix press/crescale command");
  }

  if (allremap == 0) restart_pbc = 1;

  // error checks

  if (t_start == -1)
    error->all(FLERR,"Cannot use fix press/crescale without temp specification");

  if (dimension == 2 && (p_flag[2] || p_flag[3] || p_flag[4]))
    error->all(FLERR,"Invalid fix press/crescale for a 2d simulation");
  if (dimension == 2 && (pcouple == YZ || pcouple == XZ))
    error->all(FLERR,"Invalid fix press/crescale for a 2d simulation");

  if (pcouple == XYZ && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix press/crescale pressure settings");
  if (pcouple == XYZ && dimension == 3 && p_flag[2] == 0)
    error->all(FLERR,"Invalid fix press/crescale pressure settings");
  if (pcouple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all(FLERR,"Invalid fix press/crescale pressure settings");
  if (pcouple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix press/crescale pressure settings");
  if (pcouple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all(FLERR,"Invalid fix press/crescale pressure settings");

  random = new RanMars(lmp,seed + comm->me);

  // require periodicity in tensile dimension

  if (p_flag[0] && domain->xperiodic == 0)
    error->all(FLERR,
               "Cannot use fix press/crescale on a non-periodic dimension");
  if (p_flag[1] && domain->yperiodic == 0)
    error->all(FLERR,
               "Cannot use fix press/crescale on a non-periodic dimension");
  if (p_flag[2] && domain->zperiodic == 0)
    error->all(FLERR,
               "Cannot use fix press/crescale on a non-periodic dimension");

  // require periodicity in 2nd dim of off-diagonal tilt component

  if (p_flag[3] && domain->zperiodic == 0)
    error->all(FLERR,
               "Cannot use fix press/crescale on a 2nd non-periodic dimension");
  if (p_flag[4] && domain->zperiodic == 0)
    error->all(FLERR,
               "Cannot use fix press/crescale on a 2nd non-periodic dimension");
  if (p_flag[5] && domain->yperiodic == 0)
    error->all(FLERR,
               "Cannot use fix press/crescale on a 2nd non-periodic dimension");

  if (!domain->triclinic && (p_flag[3] || p_flag[4] || p_flag[5]))
    error->all(FLERR,"Cannot specify Pxy/Pxz/Pyz in "
               "fix press/crescale with non-triclinic box");

  if (pcouple == XYZ && dimension == 3 &&
      (p_start[0] != p_start[1] || p_start[0] != p_start[2] ||
       p_stop[0] != p_stop[1] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[1] || p_period[0] != p_period[2]))
    error->all(FLERR,"Invalid fix press/crescale pressure settings");
  if (pcouple == XYZ && dimension == 2 &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR,"Invalid fix press/crescale pressure settings");
  if (pcouple == XY &&
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1] ||
       p_period[0] != p_period[1]))
    error->all(FLERR,"Invalid fix press/crescale pressure settings");
  if (pcouple == YZ &&
      (p_start[1] != p_start[2] || p_stop[1] != p_stop[2] ||
       p_period[1] != p_period[2]))
    error->all(FLERR,"Invalid fix press/crescale pressure settings");
  if (pcouple == XZ &&
      (p_start[0] != p_start[2] || p_stop[0] != p_stop[2] ||
       p_period[0] != p_period[2]))
    error->all(FLERR,"Invalid fix press/crescale pressure settings");

  if ((p_flag[0] && p_period[0] <= 0.0) ||
      (p_flag[1] && p_period[1] <= 0.0) ||
      (p_flag[2] && p_period[2] <= 0.0) ||
      (p_flag[3] && p_period[3] <= 0.0) ||
      (p_flag[4] && p_period[4] <= 0.0) ||
      (p_flag[5] && p_period[5] <= 0.0))
    error->all(FLERR,"Fix press/crescale damping parameters must be > 0.0");

  if (p_flag[0]) box_change |= BOX_CHANGE_X;
  if (p_flag[1]) box_change |= BOX_CHANGE_Y;
  if (p_flag[2]) box_change |= BOX_CHANGE_Z;
  if (p_flag[3]) box_change |= BOX_CHANGE_YZ;
  if (p_flag[4]) box_change |= BOX_CHANGE_XZ;
  if (p_flag[5]) box_change |= BOX_CHANGE_XY;
  no_change_box = 1;                    //?????????????????????????????????
  if (allremap == 0) restart_pbc = 1;   //?????????????????????????????????

  // pstyle = TRICLINIC if any off-diagonal term is controlled -> 6 dof
  // else pstyle = ISO if XYZ coupling or XY coupling in 2d -> 1 dof
  // else pstyle = ANISO -> 3 dof

  if (p_flag[3] || p_flag[4] || p_flag[5]) pstyle = TRICLINIC;
  else if (pcouple == XYZ || (dimension == 2 && pcouple == XY)) pstyle = ISO;
  else pstyle = ANISO;  

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  std::string tcmd = id + std::string("_temp");
  id_temp = new char[tcmd.size()+1];
  strcpy(id_temp,tcmd.c_str());

  tcmd += " all temp";
  modify->add_compute(tcmd);
  tflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  std::string pcmd = id + std::string("_press");
  id_press = new char[pcmd.size()+1];
  strcpy(id_press,pcmd.c_str());

  pcmd += " all pressure " + std::string(id_temp);
  modify->add_compute(pcmd);
  pflag = 1;

  nrigid = 0;
  rfix = nullptr;
}

/* ---------------------------------------------------------------------- */

FixPressCRescale::~FixPressCRescale()
{
  delete [] rfix;

  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete [] id_temp;
  delete [] id_press;
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixPressCRescale::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPressCRescale::init()
{
  // insure no conflict with fix deform

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      int *dimflag = ((FixDeform *) modify->fix[i])->dimflag;
      if ((p_flag[0] && dimflag[0]) || (p_flag[1] && dimflag[1]) ||
          (p_flag[2] && dimflag[2]) || (p_flag[3] && dimflag[3]) ||
          (p_flag[4] && dimflag[4]) || (p_flag[5] && dimflag[5]))
        error->all(FLERR,"Cannot use fix press/crescale and "
                   "fix deform on same component of stress tensor");
    }

  // set temperature and pressure ptrs

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix press/crescale does not exist");
  temperature = modify->compute[icompute];

  if (temperature->tempbias) which = BIAS;
  else which = NOBIAS;

  icompute = modify->find_compute(id_press);
  if (icompute < 0)
    error->all(FLERR,"Pressure ID for fix press/crescale does not exist");
  pressure = modify->compute[icompute];

  pdim = p_flag[0] + p_flag[1] + p_flag[2];

  // Kspace setting

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  // detect if any rigid fixes exist so rigid bodies move when box is remapped
  // rfix[] = indices to each fix rigid

  delete [] rfix;
  nrigid = 0;
  rfix = nullptr;

  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;
  if (nrigid) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->rigid_flag) rfix[nrigid++] = i;
  }
}

/* ----------------------------------------------------------------------
   compute T,P before integrator starts
------------------------------------------------------------------------- */

void FixPressCRescale::setup(int /*vflag*/)
{
  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixPressCRescale::end_of_step()
{
  // compute new T,P
  
  compute_temp_target();
  if (pstyle == ISO) {
    temperature->compute_scalar();
    pressure->compute_scalar();
  } else {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  couple();

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  /*double volume = (domain->boxhi[0] - domain->boxlo[0]) * (domain->boxhi[1] - domain->boxlo[1]);
  if (dimension == 3) volume *= (domain->boxhi[2] - domain->boxlo[2]);*/

  double volume;
  if (dimension == 3) volume = domain->xprd * domain->yprd * domain->zprd;
  else volume = domain->xprd * domain->yprd; 
 
  kt = force->boltz * t_target;
  double noise_prefactor = sqrt(2.0*kt*update->dt/(pdim*bulkmodulus*volume));

  compute_press_target();
  if (pstyle != TRICLINIC) {
    dilation[0] = dilation[1] = dilation[2] = 1.0;
    dilation[3] = dilation[4] = dilation[5] = 0.0;
    for (int i = 0; i < 3; i++) {
      if (p_flag[i]) {
        //p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
        //dilation[i] =
          //pow(1.0 - update->dt/p_period[i] *
             //(p_target[i]-p_current[i])/bulkmodulus,1.0/3.0);
        dilation[i] += - update->dt/(3.0*p_period[i]*bulkmodulus) * 
             (p_hydro-p_current[i]-kt/volume) + 
             noise_prefactor/sqrt(p_period[i]) * randoms[i];
      }
    }
  } else {
    double *h = domain->h; // Voigt order: xx, yy, zz, yz, xz, xy
    double *h_inv = domain->h_inv;
    double deltah[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    for (int i = 0; i < 3; i++) {
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
      if (p_flag[i]) {
        deltah[i] = - update->dt/(pdim*p_period[i]*bulkmodulus) * h[i] *
             (p_hydro-p_current[i]-kt/volume) +
             noise_prefactor/sqrt(p_period[i]) * h[i] * randoms[i];
      }
    }
    /*for (int i = 3; i < 6; i++) {
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
    }*/
    if (p_flag[3]) 
      deltah[3] = - update->dt/(pdim*p_period[3]*bulkmodulus) * 
           ((p_hydro-p_current[1]-kt/volume)*h[3] - p_current[3]*h[2]) +
           noise_prefactor/sqrt(p_period[3]) * (randoms[1]*h[3] + randoms[3]*h[2]);
    if (p_flag[4])
      deltah[4] = - update->dt/(pdim*p_period[4]*bulkmodulus) * 
           ((p_hydro-p_current[0]-kt/volume)*h[4] - p_current[5]*h[3] - p_current[4]*h[2]) +
           noise_prefactor/sqrt(p_period[4]) * (randoms[0]*h[4] + randoms[5]*h[3] + randoms[4]*h[2]);
    if (p_flag[5])
      deltah[5] = - update->dt/(pdim*p_period[5]*bulkmodulus) * 
           ((p_hydro-p_current[0]-kt/volume)*h[5] - p_current[5]*h[1]) +
           noise_prefactor/sqrt(p_period[5]) * (randoms[0]*h[5] + randoms[5]*h[1]);

    matrix_prod(deltah,h_inv,dilation);
    
    for (int i = 0; i < 3; i++) {
      dilation[i] += 1.0;
    }
  } 

  // remap simulation box and atoms
  // redo KSpace coeffs since volume has changed

  remap();
  if (kspace_flag) force->kspace->setup();

  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);
}

/* ---------------------------------------------------------------------- */

void FixPressCRescale::couple()
{
  double *tensor = pressure->vector;

  for(int i = 0; i < 6;i++) 
    if (p_flag[i]) 
      randoms[i] = random->gaussian(); 

  if (pstyle == ISO) {
    p_current[0] = p_current[1] = p_current[2] = pressure->scalar;
    double ave_rand = 1.0/sqrt(3.0) * (randoms[0] + randoms[1] + randoms[2]);
    randoms[0] = randoms[1] = randoms[2] = ave_rand;
  } else if (pcouple == XYZ) {
    double ave = 1.0/3.0 * (tensor[0] + tensor[1] + tensor[2]);
    p_current[0] = p_current[1] = p_current[2] = ave;
    double ave_rand = 1.0/sqrt(3.0) * (randoms[0] + randoms[1] + randoms[2]);
    randoms[0] = randoms[1] = randoms[2] = ave_rand;
  } else if (pcouple == XY) {
    double ave = 0.5 * (tensor[0] + tensor[1]);
    p_current[0] = p_current[1] = ave;
    p_current[2] = tensor[2];
    double ave_rand = sqrt(0.5) * (randoms[0] + randoms[1]);
    randoms[0] = randoms[1] = ave_rand;
  } else if (pcouple == YZ) {
    double ave = 0.5 * (tensor[1] + tensor[2]);
    p_current[1] = p_current[2] = ave;
    p_current[0] = tensor[0];
    double ave_rand = sqrt(0.5) * (randoms[1] + randoms[2]);
    randoms[1] = randoms[2] = ave_rand;
  } else if (pcouple == XZ) {
    double ave = 0.5 * (tensor[0] + tensor[2]);
    p_current[0] = p_current[2] = ave;
    p_current[1] = tensor[1];
    double ave_rand = sqrt(0.5) * (randoms[0] + randoms[2]);
    randoms[0] = randoms[2] = ave_rand;
  } else {
    p_current[0] = tensor[0];
    p_current[1] = tensor[1];
    p_current[2] = tensor[2];
  }

  if (!std::isfinite(p_current[0]) || !std::isfinite(p_current[1]) || !std::isfinite(p_current[2]))
    error->all(FLERR,"Non-numeric pressure - simulation unstable");

  // switch order from xy-xz-yz to Voigt

  if (pstyle == TRICLINIC) {
    p_current[3] = tensor[5];
    p_current[4] = tensor[4];
    p_current[5] = tensor[3];

  if (!std::isfinite(p_current[3]) || !std::isfinite(p_current[4]) || !std::isfinite(p_current[5]))
    error->all(FLERR,"Non-numeric pressure - simulation unstable");
  }
}

/* ----------------------------------------------------------------------
   change box size
   remap all atoms or fix group atoms depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass
------------------------------------------------------------------------- */

void FixPressCRescale::remap()
{
  int i, j;
  double oldlo,oldhi,ctr;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // convert pertinent atoms and rigid bodies to lamda coords

  if (allremap) domain->x2lamda(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        domain->x2lamda(x[i],x[i]);
  }

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(0);

  // reset global and local box to new size/shape

  if (pstyle != TRICLINIC) { 
    for (i = 0; i < 3; i++) {
      if (p_flag[i]) {
        oldlo = domain->boxlo[i];
        oldhi = domain->boxhi[i];
        domain->boxlo[i] = (oldlo-fixedpoint[i])*dilation[i] + fixedpoint[i];
        domain->boxhi[i] = (oldhi-fixedpoint[i])*dilation[i] + fixedpoint[i];       
      }
    }

    // rescale velocities

    if (which == NOBIAS) {
      for (i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          v[i][0] /= dilation[0];
          v[i][1] /= dilation[1];
          v[i][2] /= dilation[2];
        }
      }
    } else if (which == BIAS) {
      for (i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          temperature->remove_bias(i,v[i]);
          v[i][0] /= dilation[0];
          v[i][1] /= dilation[1];
          v[i][2] /= dilation[2];
          temperature->remove_bias(i,v[i]);
        }
      }
    }
  } 
  else {
    double *h = domain->h;
    double h_rescaled[6] = {};
    matrix_prod(dilation,h,h_rescaled);
    h = h_rescaled;

    for (i = 0; i < 3; i++) {
      domain->boxlo[i] = fixedpoint[i] - 0.5*h[i];
      domain->boxhi[i] = fixedpoint[i] + 0.5*h[i];
    }    
    domain->yz = h[3];
    domain->xz = h[4];
    domain->xy = h[5];

    if (domain->yz < -TILTMAX*domain->yprd ||
      domain->yz > TILTMAX*domain->yprd ||
      domain->xz < -TILTMAX*domain->xprd ||
      domain->xz > TILTMAX*domain->xprd ||
      domain->xy < -TILTMAX*domain->xprd ||
      domain->xy > TILTMAX*domain->xprd)
      error->all(FLERR,"Fix press/crescale has tilted box too far in one step - "
               "periodic cell is too far from equilibrium state");

    // rescale velocities

    double v_rescaled[3] = {};
    double dilation_inv[6] = {};
    inverse_matrix(dilation,dilation_inv);
    if (which == NOBIAS) {
      for (i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          vector_matrix_prod(v[i],dilation_inv,v_rescaled);    
          v[i][0] = v_rescaled[0];
          v[i][1] = v_rescaled[1];
          v[i][2] = v_rescaled[2];
        }
      }
    } else if (which == BIAS) {
      for (i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          temperature->remove_bias(i,v[i]);
          vector_matrix_prod(v[i],dilation_inv,v_rescaled);    
          v[i][0] = v_rescaled[0];
          v[i][1] = v_rescaled[1];
          v[i][2] = v_rescaled[2];
          temperature->remove_bias(i,v[i]);
        }
      }
    }  
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  if (allremap) domain->lamda2x(nlocal);
  else {
    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        domain->lamda2x(x[i],x[i]);
  }

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(1);
}

/* ---------------------------------------------------------------------- */

int FixPressCRescale::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,"Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning(FLERR,"Temperature for NPT is not for group all");

    // reset id_temp of pressure to new temperature ID

    icompute = modify->find_compute(id_press);
    if (icompute < 0)
      error->all(FLERR,"Pressure ID for fix press/crescale does not exist");
    modify->compute[icompute]->reset_extra_compute_fix(id_temp);

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (pflag) {
      modify->delete_compute(id_press);
      pflag = 0;
    }
    delete [] id_press;
    int n = strlen(arg[1]) + 1;
    id_press = new char[n];
    strcpy(id_press,arg[1]);

    int icompute = modify->find_compute(arg[1]);
    if (icompute < 0) error->all(FLERR,"Could not find fix_modify pressure ID");
    pressure = modify->compute[icompute];

    if (pressure->pressflag == 0)
      error->all(FLERR,"Fix_modify pressure ID does not compute pressure");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixPressCRescale::compute_temp_target()
{
  // compute reference temperature (target temperature of thermostat)

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  t_target = t_start + delta * (t_stop-t_start);
}

/* ---------------------------------------------------------------------- */

void FixPressCRescale::matrix_prod(double *in1, double *in2, double *out) 
{
  // compute matrix product in1 * in2 (in1, in2 and output represented in Voigt order)

  out[0] = in1[0] * in2[0];
  out[1] = in1[1] * in2[1];
  out[2] = in1[2] * in2[2];
  out[3] = in1[1] * in2[3] + in1[3] * in2[2];
  out[4] = in1[0] * in2[4] + in1[5] * in2[3] + in1[4] * in2[2];
  out[5] = in1[0] * in2[5] + in1[5] * in2[1];
}

/* ---------------------------------------------------------------------- */

void FixPressCRescale::vector_matrix_prod(double *in_v, double *in_m, double *out_v) 
{
  // compute vector-matrix product in_v * in_m (matrix in_m represented in Voigt order)

  out_v[0] = in_v[0] * in_m[0];
  out_v[1] = in_v[0] * in_m[5] + in_v[1] * in_m[1];
  out_v[2] = in_v[0] * in_m[4] + in_v[1] * in_m[3] + in_v[2] * in_m[2];  
}

/* ---------------------------------------------------------------------- */

void FixPressCRescale::inverse_matrix(double *in, double *out) 
{
  // invert upper triangular matrix in (represented in Voigt order)

  out[0] = 1/in[0];
  out[1] = 1/in[1];
  out[2] = 1/in[2];
  out[3] = -in[3] / (in[1]*in[2]);
  out[4] = (in[3]*in[5] - in[1]*in[4]) / (in[0]*in[1]*in[2]);
  out[5] = -in[5] / (in[0]*in[1]);
}

void FixPressCRescale::compute_press_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  p_hydro = 0.0;
  for (int i = 0; i < 3; i++)
    if (p_flag[i]) {
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
      p_hydro += p_target[i];
    }
  if (pdim > 0) p_hydro /= pdim;

  if (pstyle == TRICLINIC)
    for (int i = 3; i < 6; i++)
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
}
