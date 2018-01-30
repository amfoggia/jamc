#include "ljmd.h"
#include <stdio.h>
#ifdef _OPENMP
#include "omp.h"
#endif

/* main */
int main(int argc, char **argv)
{
    int nprint, i;
    int nprocs, rank;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    mdsys_t sys;
    FILE *fp;
    cell_t *cel;
    
#ifdef _OPENMP
    int nthds = omp_get_max_threads(); // openmp
#else
    int nthds = 1;
#endif
#ifdef TIMING
    double timer = 0; // for recording time
#endif
    
#ifdef USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    nprocs = 1;
    rank = 0;
#endif /* USE_MPI */
    
    
    // Initialize the struct to zero
    memset(&sys, 0, sizeof(sys));
    
    if (rank == 0) {
        /* read input file */
        if(get_a_line(stdin,line)) return 1;
        sys.natoms=atoi(line);
        if(get_a_line(stdin,line)) return 1;
        sys.redmass = atof(line)*mvsq2e; // reduced mass
        if(get_a_line(stdin,line)) return 1;
        sys.epsilon=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.sigma=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.rcut=atof(line);
        if(get_a_line(stdin,line)) return 1;
        sys.box=atof(line);
        if(get_a_line(stdin,restfile)) return 1;
        if(get_a_line(stdin,trajfile)) return 1;
        if(get_a_line(stdin,ergfile)) return 1;
        if(get_a_line(stdin,line)) return 1;
        sys.nsteps=atoi(line);
        if(get_a_line(stdin,line)) return 1;
        sys.dt=atof(line);
        if(get_a_line(stdin,line)) return 1;
        nprint=atoi(line);
    }

#ifdef CELL
    sys.cl = sys.rcut;
    sys.cn = floor(sys.box/sys.cl);

    // Allocate memory for the array of cells
    cel = new(cell_t[sys.cn*sys.cn*sys.cn]);
    
    // Create the list of cell pairs that are going to be used in the force calculation
    pair(&sys);
#endif /* CELL */
    
#ifdef USE_MPI
    // Communicate the struct to the rest of the processes
    MPI_Bcast(&sys, sizeof(sys), MPI_CHAR, 0, MPI_COMM_WORLD);
#endif /* USE_MPI */
    
    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*nthds*sizeof(double));
    sys.fy=(double *)malloc(sys.natoms*nthds*sizeof(double));
    sys.fz=(double *)malloc(sys.natoms*nthds*sizeof(double));
    // Only the master process will perform the verlet algorithm
    if (rank == 0) {
        sys.vx=(double *)malloc(sys.natoms*sizeof(double));
        sys.vy=(double *)malloc(sys.natoms*sizeof(double));
        sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    }
    
    /* read restart */    
    if (rank == 0) {
        fp=fopen(restfile,"r");
        if(fp) {
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
            }
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
            }
            fclose(fp);
            azzero(sys.fx, sys.natoms);
            azzero(sys.fy, sys.natoms);
            azzero(sys.fz, sys.natoms);
        } else {
            perror("cannot read restart file");
            return 3;
        }
    }
    
#ifdef CELL
    // Sort the particles inside the cells
    sort(&sys, cel);
#endif /* CELL */    
    
    /* initialize forces and energies.*/
    sys.nfi=0;
    
#ifdef CELL
    cell_force(&sys, cel);
#else
    force(&sys);
#endif /* CELL */    
    
    if (rank == 0)
        ekin(&sys);
    
#ifndef TIMING
    FILE *traj,*erg;
    if (rank == 0) {
        erg=fopen(ergfile,"w");
        traj=fopen(trajfile,"w");
        printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        output(&sys, erg, traj);
    }
#endif /* NO TIMING */
    
    
    /**************************************************/
    /* main MD loop */
    
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
        /* propagate system and recompute energies */
#ifdef TIMING
        velverlet_1(&sys);

	if (rank == 0)
	  timer -= stamp();
#ifdef CELL
	// Sort the particles inside the cells
	sort(&sys, cel);

	cell_force(&sys, cel);
#else
	force(&sys);
#endif /* CELL */
        
        if (rank == 0)
            timer += stamp();
        
        velverlet_2(&sys);
        
#else /* TIMING */
        /* write output, if requested */
        if (rank == 0) {
            if ((sys.nfi % nprint) == 0)
                output(&sys, erg, traj);
        }
        
        if (rank == 0)
            velverlet_1(&sys);

#ifdef CELL
	// Sort the particles inside the cells
	sort(&sys, cel);
        cell_force(&sys, cel);
#else
	force(&sys);
#endif /* CELL */
        
        if (rank == 0)
            velverlet_2(&sys);
#endif /* TIMING */
        
        if (rank == 0)
            ekin(&sys);
    }
    /**************************************************/
    
    /* clean up: close files, free memory */
    if (rank == 0) {
#ifdef TIMING
        printf("\n timing \n\t %f ms \n\n",timer);
#else /* TIMING */
        printf("Simulation Done.\n");
        fclose(erg);
        fclose(traj);
#endif /* TIMING */
    }
    
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
    
#ifdef CELL    
    delete [] cel;
#endif    
    
#ifdef USE_MPI
    MPI_Finalize();
#endif /* USE_MPI */
    
    return 0;
}
