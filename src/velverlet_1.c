#include "ljmd.h"

/* velocity verlet */
void velverlet_1(mdsys_t *sys)
{
    int i;
    double coef = 0.5*sys->dt/ (sys->redmass);
    /* first part: propagate velocities by half and positions by full step */
    for (i=0; i<sys->natoms; ++i) {
        sys->vx[i] += coef * sys->fx[i] ;
        sys->vy[i] += coef * sys->fy[i] ;
        sys->vz[i] += coef * sys->fz[i] ;
        sys->rx[i] += sys->dt*sys->vx[i];
        sys->ry[i] += sys->dt*sys->vy[i];
        sys->rz[i] += sys->dt*sys->vz[i];
    }

}
