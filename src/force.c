#include "ljmd.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/* newton omp aggressive with replicated memory */
void force(mdsys_t *sys)
{
    double rsq,rsq_inv,r6,ffac;
    double rx,ry,rz;
    //printf("%f\n",sys->rx[0]);
    
    
    int n,i,j,offset;
    double * fx, * fy, * fz; // local pointers for openmp
    double epot=0.0; // needed for reduction with openmp
    int niters = (sys->natoms) * (sys->natoms - 1) / 2; // to linearize (and balance) the loop with openmp / mpi
    // accordingly to MPI policies
    
    int nprocs, rank, local_niter, lower_bound, upper_bound;
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rest = niters%nprocs;
    if (rank < rest) {
        local_niter = niters/nprocs + 1;
        lower_bound = local_niter * rank;
        upper_bound = (rank + 1) * local_niter;
    }
    else {
        local_niter = niters/nprocs;
        lower_bound = local_niter * rank + rest;
        upper_bound = (rank + 1) * local_niter + rest;
    }
    // Broadcast the positions to all processes
    MPI_Bcast(sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    nprocs = 1;
    rank = 0;
    local_niter= niters;
    lower_bound = 0;
    upper_bound = niters;
#endif /* USE_MPI */
    
    double boxby2 = 0.5*sys->box;// pre-calculate
    double rcutsq = sys->rcut*sys->rcut;// pre-calculate, take square
    double c6 = sys->epsilon*pow(sys->sigma,6);
    double c12 = sys->epsilon*pow(sys->sigma,12);
    
#ifdef _OPENMP
#pragma omp parallel private(fx,fy,fz,i,j,rx,ry,rz,rsq,rsq_inv,r6,ffac,offset) reduction(+:epot)
#endif
    {
        
        // who is who
        int thid; // thread id
        int nthds; // number of threads
        
#ifdef _OPENMP
        thid = omp_get_thread_num();
        nthds = omp_get_num_threads();
#else
        thid = 0;
        nthds = 1;
#endif
        
        // assign slices of force array to threads
        fx = sys->fx + (thid * sys->natoms);
        fy = sys->fy + (thid * sys->natoms);
        fz = sys->fz + (thid * sys->natoms);
        
        // zero energy and forces
        azzero(fx,sys->natoms);
        azzero(fy,sys->natoms);
        azzero(fz,sys->natoms);
        
        // main loop
#ifdef _OPENMP
#ifndef CHUNKSIZE
#define CHUNKSIZE 1
#endif
#pragma omp for schedule(guided,CHUNKSIZE)
#endif
        for(n=lower_bound; n<upper_bound; ++n) {

            // compute indexes as a function of n
            i = sys->natoms - 2 - (int)(sqrt(-8*n + 4*sys->natoms*(sys->natoms-1)-7)/2.0 - 0.5);
            j = n + i + 1 - sys->natoms*(sys->natoms-1)/2 + (sys->natoms-i)*((sys->natoms-i)-1)/2;

            // get distance between particle i and j
            rx=pbc(sys->rx[i] - sys->rx[j], boxby2);
            rsq = rx*rx;
            if(rsq>rcutsq) continue;
            ry=pbc(sys->ry[i] - sys->ry[j], boxby2);
            rsq += ry*ry;
            if(rsq>rcutsq) continue;
            rz=pbc(sys->rz[i] - sys->rz[j], boxby2);
            rsq += rz*rz;
            if(rsq>rcutsq) continue;
            
            rsq_inv = 1.0/rsq;
            r6 = rsq_inv*rsq_inv*rsq_inv;
            ffac = (48*c12*r6-24*c6)*r6*rsq_inv;
            epot += 4*r6*(c12*r6-c6);
            fx[i] += rx*ffac;
            fx[j] -= rx*ffac;
            fy[i] += ry*ffac;
            fy[j] -= ry*ffac;
            fz[i] += rz*ffac;
            fz[j] -= rz*ffac;
        }
        
        // after a work sharing construct an omp barrier is implied
#ifdef _OPENMP
#pragma omp for
#endif
        for (i=0; i<sys->natoms; ++i) {
            for (offset=sys->natoms; offset<nthds*sys->natoms; offset+=sys->natoms) {
                sys->fx[i] += sys->fx[offset + i];
                sys->fy[i] += sys->fy[offset + i];
                sys->fz[i] += sys->fz[offset + i];
            }
        }
        
    }
    
    sys->epot = epot;
    
#ifdef USE_MPI
    // Communicate forces
    if (rank == 0){
        MPI_Reduce(MPI_IN_PLACE, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Reduce(sys->fx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(sys->fy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(sys->fz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    
    // Communicate epot
    if (rank == 0)
        MPI_Reduce(MPI_IN_PLACE, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    else
        MPI_Reduce(&sys->epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif /* USE_MPI */

}
