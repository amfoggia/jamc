#include "ljmd.h"
#include <sys/time.h>

/* helper function: zero out an array */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* compute kinetic energy */
void ekin(mdsys_t *sys)
{   
    int i;
    int local_size, lower_bound, upper_bound;

    varsize(sys->natoms, &local_size, &lower_bound, &upper_bound);
    
    sys->ekin=0.0;
    for (i=0; i<local_size; ++i) {
        sys->ekin += 0.5*sys->redmass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

/* timing */
double stamp(){
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*1e+3 + tv.tv_usec*1e-3;
}

/* helper function: local size variables for mpi */
void varsize(int natoms, int * local_size, int * lower_bound, int * upper_bound) {

  int rank, nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int rest = natoms%nprocs;

  if (rank < rest) {
    *local_size = natoms/nprocs + 1;
    int a = *local_size;
    *lower_bound = a*rank;
    *upper_bound = a*(rank + 1);
  }
  else {
    *local_size = natoms/nprocs;
    int b= *local_size;
    *lower_bound = b*rank + rest;
    *upper_bound = b*(rank + 1) + rest;
  }
}
