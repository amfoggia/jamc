// reference output of this test should be test_force.dat

#include "ljmd.h"

int main(int argc, char * argv[]) {

  char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
  FILE *fp, *fforce;
  mdsys_t sys;
  int i;

  // Read input file: Initialize the sys object
  if(get_a_line(stdin,line)) return 1;
  sys.natoms=atoi(line);
  if(get_a_line(stdin,line)) return 1;
  sys.mass=atof(line);
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
  
  /* allocate memory */
  sys.rx=(double *)malloc(sys.natoms*sizeof(double));
  sys.ry=(double *)malloc(sys.natoms*sizeof(double));
  sys.rz=(double *)malloc(sys.natoms*sizeof(double));
  sys.vx=(double *)malloc(sys.natoms*sizeof(double));
  sys.vy=(double *)malloc(sys.natoms*sizeof(double));
  sys.vz=(double *)malloc(sys.natoms*sizeof(double));
  sys.fx=(double *)malloc(sys.natoms*sizeof(double));
  sys.fy=(double *)malloc(sys.natoms*sizeof(double));
  sys.fz=(double *)malloc(sys.natoms*sizeof(double));

  /* read restart */
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

  // Call the force function
  force(&sys);

  fforce =  fopen("test_force.dat", "a");

  for (i = 0; i < (sys.natoms); ++i)
    fprintf(fforce, "%d%20.8f%20.8f%20.8f\n", i, sys.fx[i], sys.fy[i], sys.fz[i]);

  fclose(fforce);
  
  free(sys.rx);
  free(sys.ry);
  free(sys.rz);
  free(sys.vx);
  free(sys.vy);
  free(sys.vz);
  free(sys.fx);
  free(sys.fy);
  free(sys.fz);
  
  return 0;
}
