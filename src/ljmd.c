#include "ljmd.h"
//#include "cell.h"

/* main */
int main(int argc, char **argv) 
{
    int nprint, i, proc;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double timer=0; // for recording time
    double total_epot, total_temp, total_ekin;

    // Parallel part with MPI ----------------------------------------------

    /* MPI variables*/
    int nprocs, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

    MPI_Bcast(&sys.natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.redmass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.rcut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.box, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Request sendrequest[3];
    MPI_Request recvrequest[3];

    // Initialize requests
    for (i = 0; i < 3; ++i) {
      sendrequest[i] = MPI_REQUEST_NULL;
      recvrequest[i] = MPI_REQUEST_NULL;
    }
    
    /* Size variables */
    int local_size = 0;
    int lower_bound = 0;
    int upper_bound = 0;

    varsize(sys.natoms, &local_size, &lower_bound, &upper_bound);

    printf("FIRST: Process %d localsize= %d\n", rank, local_size);    

    // Necessary arrays for the MPI_Allgatherv
    int recv_size[nprocs];
    int recv_offset[nprocs];

    for (i = 0; i < nprocs; ++i) {
      if (i < sys.natoms%nprocs) {
	recv_size[i] = sys.natoms/nprocs + 1;
	recv_offset[i] = i*recv_size[i];
      }
      else {
	recv_size[i] = sys.natoms/nprocs;
	recv_offset[i] = i*recv_size[i] + sys.natoms%nprocs;
      }
    }

    printf("HELLO\n");
    int sys_size_bites = sys.natoms*sizeof(double);
    int local_size_bites = local_size*sizeof(double);
    
    /* allocate memory */
    sys.rx=(double *)malloc(sys_size_bites);
    sys.ry=(double *)malloc(sys_size_bites);
    sys.rz=(double *)malloc(sys_size_bites);
    sys.vx=(double *)malloc(local_size_bites);
    sys.vy=(double *)malloc(local_size_bites);
    sys.vz=(double *)malloc(local_size_bites);
    sys.fx=(double *)malloc(local_size_bites);
    sys.fy=(double *)malloc(local_size_bites);
    sys.fz=(double *)malloc(local_size_bites);

    // read restart
    if (rank == 0) {
      FILE * fp=fopen(restfile,"r");
      if(fp) {
        for (i=0; i<sys.natoms; ++i) {
      	  fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
      	// to read velocity only process 0 does it
      	// First the velocity of process 0
      	for (i=lower_bound; i<upper_bound; ++i) {
      	  fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
      	}
	
	for (proc = 1; proc < nprocs; ++proc) {
	    
      	  if (proc < sys.natoms%nprocs)
      	    local_size = sys.natoms/nprocs + 1;
      	  else
      	    local_size = sys.natoms/nprocs;

      	  printf("local_size= %d\n", local_size);
	  // Now, the velocities of every process
        double * helper_buffx = (double *) malloc(local_size_bites);
        double * helper_buffy = (double *) malloc(local_size_bites);
        double * helper_buffz = (double *) malloc(local_size_bites);

      	  for (i=0; i<local_size; ++i) {
      	    fscanf(fp,"%lf%lf%lf",helper_buffx+i, helper_buffy+i, helper_buffz+i);
      	  }
 
      	  // Communicate it to the corresponding process
      	  MPI_Isend(helper_buffx, local_size, MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, &sendrequest[0]);
      	  MPI_Isend(helper_buffy, local_size, MPI_DOUBLE, proc, 2, MPI_COMM_WORLD, &sendrequest[1]);
      	  MPI_Isend(helper_buffz, local_size, MPI_DOUBLE, proc, 3, MPI_COMM_WORLD, &sendrequest[2]);

	  free(helper_buffx);
	  free(helper_buffy);
	  free(helper_buffz);	
	}
      }
      else {
      perror("cannot read restart file");
      return 3;
      }
      fclose(fp);
    }
    else {
      MPI_Irecv(sys.vx, local_size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &recvrequest[0]);
      MPI_Irecv(sys.vy, local_size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &recvrequest[1]);
      MPI_Irecv(sys.vz, local_size, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &recvrequest[2]);
    }
    
    varsize(sys.natoms, &local_size, &lower_bound, &upper_bound);
    
    MPI_Bcast(sys.rx, sys.natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys.ry, sys.natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(sys.rz, sys.natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Each process initializes its local force to zero
    azzero(sys.fx, local_size);
    azzero(sys.fy, local_size);
    azzero(sys.fz, local_size);
    
    MPI_Waitall(3, sendrequest, MPI_STATUS_IGNORE);
    MPI_Waitall(3, recvrequest, MPI_STATUS_IGNORE);

    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys);
    
    ekin(&sys);

    // As we want to print the energy and the temperature we need to communicate the partial values of
    // sys.ekin, sys.epot and sys.temp to process 0 and sum all the contributions
    MPI_Reduce(&sys.epot, &total_epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sys.ekin, &total_ekin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sys.temp, &total_temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {

      FILE * erg = fopen(ergfile,"w");
      FILE * traj = fopen(trajfile,"w");

      sys.epot = total_epot;
      sys.ekin = total_ekin;
      sys.temp = total_temp;
      
      printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
      printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
      output(&sys, erg, traj);

      fclose(erg);
      fclose(traj);
    }

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

      /* write output, if requested, ONLY PROCESS 0 */
      /* if (rank == 0) { */
      /*   if ((sys.nfi % nprint) == 0) { */
      
      /* 	  FILE * erg=fopen(ergfile,"w"); */
      /*     FILE * traj=fopen(trajfile,"w"); */
	  
      /* 	  output(&sys, erg, traj); */
      /* 	} */
      /* } */

      
      if (rank == 0) {
	if ((sys.nfi % nprint) == 0){
	  sys.epot = total_epot;
	  sys.ekin = total_ekin;
	  sys.temp = total_temp;
	  printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys.nfi, sys.temp, sys.ekin, sys.epot, sys.ekin+sys.epot);
	}
      }
        
      /* propagate system and recompute energies */
      velverlet_1(&sys);

      // We need the barrier because the positions on every process need to be already
      // calculated in order to be communicated to the other processes
      MPI_Barrier(MPI_COMM_WORLD);

      // Communicate positions before calculating forces
      MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sys.rx, recv_size, recv_offset, MPI_DOUBLE, MPI_COMM_WORLD);
      /* MPI_Allgatherv(&sys.ry[lower_bound], local_size, MPI_DOUBLE, sys.rx, recv_size, recv_offset, MPI_DOUBLE, MPI_COMM_WORLD); */
      /* MPI_Allgatherv(&sys.rz[lower_bound], local_size, MPI_DOUBLE, sys.rx, recv_size, recv_offset, MPI_DOUBLE, MPI_COMM_WORLD); */
	
      force(&sys);

      velverlet_2(&sys);
	
      ekin(&sys);

      // Bring from every process sys.epot, sys.ekin and sys.temp and sum the contributions
      MPI_Reduce(&sys.epot, &total_epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&sys.ekin, &total_ekin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(&sys.temp, &total_temp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    /* if (rank == 0) { */
    /*   //      FILE * erg = fopen(ergfile,"w"); */
    /*   //FILE * traj = fopen(trajfile,"w"); */
    /*   printf("Simulation Done.\n"); */
    /*   printf("\n timing: %.2f\n\n",timer); */
    /*   //fclose(erg); */
    /*   // fclose(traj); */
    /* } */

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);
    
    MPI_Finalize();
    return 0;
}
