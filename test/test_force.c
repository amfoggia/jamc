#include "ljmd.h"
#include <stdlib.h>

int main(int argc, char * argv[]) {

  mdsys_t sys3, sys4 ;
  int i;
  FILE * fp;

  fp = fopen("test_force.dat", "a");

  // Initialize the system
  sys3.natoms = 3;
  sys3.redmass = 39.948 * mvsq2e;
  sys3.epsilon = 0.2379;
  sys3.sigma = 3.405;
  sys3.rcut = 8.5;
  sys3.box = 20.0;
  
  // Allocate memory
  sys3.rx=(double *)malloc(sys3.natoms*sizeof(double));
  sys3.ry=(double *)malloc(sys3.natoms*sizeof(double));
  sys3.rz=(double *)malloc(sys3.natoms*sizeof(double));
  sys3.fx=(double *)malloc(sys3.natoms*sizeof(double));
  sys3.fy=(double *)malloc(sys3.natoms*sizeof(double));
  sys3.fz=(double *)malloc(sys3.natoms*sizeof(double));

  // FIRST CASE: all far apart
  // Initialize the positions and function
  // Particle 1
  sys3.rx[0] = 8.51;
  sys3.ry[0] = 0.0;
  sys3.rz[0] = 0.0;
 
  // Particle 2
  sys3.rx[1] = 0.0;
  sys3.ry[1] = 0.0;
  sys3.rz[1] = 0.0;
  
  // Particle 3
  sys3.rx[2] = 0.0;
  sys3.ry[2] = 9.1;
  sys3.rz[2] = 0.0;

  // Call the force function
  force(&sys3);

  for (i = 0; i < sys3.natoms; ++i) {
    fprintf(fp, "%d%20.8f%20.8f%20.8f\n", i, sys3.fx[i], sys3.fy[i], sys3.fz[i]);
  }

  // SECOND CASE: two inside the cutoff, one outside
  // Initialize the positions and function
  // Particle 1
  sys3.rx[0] = 8.51;
  sys3.ry[0] = 0.0;
  sys3.rz[0] = 0.0;

  // Particle 2
  sys3.rx[1] = 0.0;
  sys3.ry[1] = 4.2;
  sys3.rz[1] = 0.0;

  // Particle 3
  sys3.rx[2] = 0.0;
  sys3.ry[2] = 0.0;
  sys3.rz[2] = 0.0;
  
  // Call the force function
  force(&sys3);

  for (i = 0; i < sys3.natoms; ++i) {
    fprintf(fp, "%d%20.8f%20.8f%20.8f\n", i, sys3.fx[i], sys3.fy[i], sys3.fz[i]);
  }

  // THIRD CASE: all inside the cutoff
  // Initialize the positions and function
  // Particle 1
  sys3.rx[0] = 1.0;
  sys3.ry[0] = 1.0;
  sys3.rz[0] = 1.0;

  // Particle 2
  sys3.rx[1] = -1.5;
  sys3.ry[1] = -1.5;
  sys3.rz[1] = -1.5;

  // Particle 3
  sys3.rx[2] = 0.0;
  sys3.ry[2] = 0.0;
  sys3.rz[2] = 0.0;
  
  // Call the force function
  force(&sys3);

  for (i = 0; i < sys3.natoms; ++i) {
    fprintf(fp, "%d%20.8f%20.8f%20.8f\n", i, sys3.fx[i], sys3.fy[i], sys3.fz[i]);
  }
  
  free(sys3.rx);
  free(sys3.ry);
  free(sys3.rz);
  free(sys3.fx);
  free(sys3.fy);
  free(sys3.fz);

  // FOURTH CASE: four particles, three inside the cutoff, one outside

  // Initialize the system
  sys4.natoms = 4;
  sys3.redmass = 39.948 * mvsq2e;
  sys4.epsilon = 0.2379;
  sys4.sigma = 3.405;
  sys4.rcut = 8.5;
  sys4.box = 17.1580;
  
  // Allocate memory
  sys4.rx=(double *)malloc(sys4.natoms*sizeof(double));
  sys4.ry=(double *)malloc(sys4.natoms*sizeof(double));
  sys4.rz=(double *)malloc(sys4.natoms*sizeof(double));
  sys4.fx=(double *)malloc(sys4.natoms*sizeof(double));
  sys4.fy=(double *)malloc(sys4.natoms*sizeof(double));
  sys4.fz=(double *)malloc(sys4.natoms*sizeof(double));
  
  // Initialize the positions and function
  // Particle 1
  sys4.rx[0] = 8.51;
  sys4.ry[0] = 0.0;
  sys4.rz[0] = 0.0;
  
  // Particle 2
  sys4.rx[1] = 0.0;
  sys4.ry[1] = 9.3;
  sys4.rz[1] = 0.0;

  // Particle 3
  sys4.rx[2] = 0.0;
  sys4.ry[2] = 0.0;
  sys4.rz[2] = 8.98;

  // Particle 4
  sys4.rx[3] = 0.0;
  sys4.ry[3] = 0.0;
  sys4.rz[3] = 0.0;

  // Call the force function
  force(&sys4);

  for (i = 0; i < sys4.natoms; ++i) {
    fprintf(fp, "%d%20.8f%20.8f%20.8f\n", i, sys4.fx[i], sys4.fy[i], sys4.fz[i]);
  }
  
  free(sys4.rx);
  free(sys4.ry);
  free(sys4.rz);
  free(sys4.fx);
  free(sys4.fy);
  free(sys4.fz);
  
  return 0;
}
