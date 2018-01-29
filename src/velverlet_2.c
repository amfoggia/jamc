#include "ljmd.h"

/* velocity verlet */
void velverlet_2(mdsys_t *sys)
{
    int i;
    double coef = 0.5*sys->dt / sys->redmass;
    int local_size, lower_bound, upper_bound;
    
    varsize(sys->natoms, &local_size, &lower_bound, &upper_bound);

    /* second part: propagate velocities by another half step */
    for (i=0; i<local_size; ++i) {
        sys->vx[i] += coef * sys->fx[i];
        sys->vy[i] += coef * sys->fy[i];
        sys->vz[i] += coef * sys->fz[i];
    }
}
