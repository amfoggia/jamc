#include "ljmd.h"



/* main */
int main(int argc, char **argv) {

    mdsys_t sys;
    int max_test_sample = 4; // max number of atoms used for testing
    int i;
    FILE *ofp;

    /* allocate memory */
    sys.rx=(double *)malloc(max_test_sample * sizeof(double));
    sys.ry=(double *)malloc(max_test_sample * sizeof(double));
    sys.rz=(double *)malloc(max_test_sample * sizeof(double));
    sys.vx=(double *)malloc(max_test_sample * sizeof(double));
    sys.vy=(double *)malloc(max_test_sample * sizeof(double));
    sys.vz=(double *)malloc(max_test_sample * sizeof(double));
    sys.fx=(double *)malloc(max_test_sample * sizeof(double));
    sys.fy=(double *)malloc(max_test_sample * sizeof(double));
    sys.fz=(double *)malloc(max_test_sample * sizeof(double));

    /* open output file */
    ofp = fopen("test_velverlet_2.dat","w");

    /* TEST # 1a: 3 particles far away form each other; v = 0 */

    /* initialize system */
    sys.natoms = 3;
    sys.redmass = 39.948 * mvsq2e;
    sys.epsilon = 0.2379;
    sys.sigma = 3.405;
    sys.rcut = 8.5;
    sys.box = 17.1580;
    sys.nsteps = 1;
    sys.dt = 5.0;

    /* initialize other members (not needed actually) */
    sys.nfi=0;
    sys.ekin=1.;
    sys.epot=0.;
    sys.temp=1.;

    /* initialize system positions */
    sys.rx[0] = 0.;
    sys.ry[0] = 0.;
    sys.rz[0] = 0.;

    sys.rx[1] = 8.54;
    sys.ry[1] = 0.;
    sys.rz[1] = 0.;

    sys.rx[2] = 0.;
    sys.ry[2] = 8.54;
    sys.rz[2] = 0.;

    /* initialize system velocities to 0 */
    azzero(sys.vx, sys.natoms);
    azzero(sys.vy, sys.natoms);
    azzero(sys.vz, sys.natoms);

    /* initialize system forces */
    sys.fx[0] = 0.;
    sys.fy[0] = 0.;
    sys.fz[0] = 0.;

    sys.fx[1] = 0.;
    sys.fy[1] = 0.;
    sys.fz[1] = 0.;

    sys.fx[2] = 0.;
    sys.fy[2] = 0.;
    sys.fz[2] = 0.;

    /* propagate system with velverlet and check */
    velverlet_2(&sys);

    fprintf(ofp, "TEST # 1a: 3 particles far away form each other; v = 0\n");
    for (i=0; i < sys.natoms; ++i) {
      fprintf(ofp, "r: %20.8f%20.8f%20.8f\n", sys.rx[i], sys.ry[i], sys.rz[i]);
      fprintf(ofp, "v: %20.8f%20.8f%20.8f\n", sys.vx[i], sys.vy[i], sys.vz[i]);
    }

    /* ************************************************* */



    /* TEST # 1b: 3 particles far away form each other; v != 0 */

    /* initialize system positions */
    sys.rx[0] = 0.;
    sys.ry[0] = 0.;
    sys.rz[0] = 0.;

    sys.rx[1] = 8.54;
    sys.ry[1] = 0.;
    sys.rz[1] = 0.;

    sys.rx[2] = 0.;
    sys.ry[2] = 8.54;
    sys.rz[2] = 0.;

    /* system velocities are not reinitialized */

    /* system forces do not need to be reinitialized */

    /* propagate system with velverlet and check */
    velverlet_2(&sys);

    fprintf(ofp, "TEST # 1b: 3 particles far away form each other; v != 0\n");
    for (i=0; i < sys.natoms; ++i) {
      fprintf(ofp, "r: %20.8f%20.8f%20.8f\n", sys.rx[i], sys.ry[i], sys.rz[i]);
      fprintf(ofp, "v: %20.8f%20.8f%20.8f\n", sys.vx[i], sys.vy[i], sys.vz[i]);
    }

    /* ************************************************* */



    /* TEST # 2a: 3 particles, 2 close to each other; v = 0 */

    /* initialize system */

    /* initialize system positions */
    sys.rx[0] = 0.;
    sys.ry[0] = 0.;
    sys.rz[0] = 0.;

    sys.rx[1] = 2.;
    sys.ry[1] = 2.;
    sys.rz[1] = 0.;

    sys.rx[2] = 0.;
    sys.ry[2] = 0.;
    sys.rz[2] = 8.54;

    /* initialize system velocities to 0 */
    azzero(sys.vx, sys.natoms);
    azzero(sys.vy, sys.natoms);
    azzero(sys.vz, sys.natoms);

    /* initialize system forces */
    sys.fx[0] = -22.10605666;
    sys.fy[0] = -22.10605666;
    sys.fz[0] = 0.;

    sys.fx[1] = 22.10605666;
    sys.fy[1] = 22.10605666;
    sys.fz[1] = 0.;

    sys.fx[2] = 0.;
    sys.fy[2] = 0.;
    sys.fz[2] = 0.;

    /* propagate system with velverlet and check */
    velverlet_2(&sys);

    fprintf(ofp, "TEST # 2a: 3 particles, 2 close to each other; v = 0\n");
    for (i=0; i < sys.natoms; ++i) {
      fprintf(ofp, "r: %20.8f%20.8f%20.8f\n", sys.rx[i], sys.ry[i], sys.rz[i]);
      fprintf(ofp, "v: %20.8f%20.8f%20.8f\n", sys.vx[i], sys.vy[i], sys.vz[i]);
    }

    /* ************************************************* */



    /* TEST # 2b: 3 particles, 2 close to each other; v != 0 */

    /* initialize system */

    /* initialize system positions */
    sys.rx[0] = 0.;
    sys.ry[0] = 0.;
    sys.rz[0] = 0.;

    sys.rx[1] = 2.;
    sys.ry[1] = 2.;
    sys.rz[1] = 0.;

    sys.rx[2] = 0.;
    sys.ry[2] = 0.;
    sys.rz[2] = 8.54;

    /* system velocities are not reinitialized */

    /* system forces do not need to be reinitialized */

    /* propagate system with velverlet and check */
    velverlet_2(&sys);

    fprintf(ofp, "TEST # 2b: 3 particles, 2 close to each other; v != 0\n");
    for (i=0; i < sys.natoms; ++i) {
      fprintf(ofp, "r: %20.8f%20.8f%20.8f\n", sys.rx[i], sys.ry[i], sys.rz[i]);
      fprintf(ofp, "v: %20.8f%20.8f%20.8f\n", sys.vx[i], sys.vy[i], sys.vz[i]);
    }

    /* ************************************************* */



    /* TEST # 3a: 3 particles, 3 close to each other; v = 0 */

    /* initialize system */

    /* initialize system positions */
    sys.rx[0] = 1.5;
    sys.ry[0] = 0.;
    sys.rz[0] = 0.;

    sys.rx[1] = 0.;
    sys.ry[1] = 1.5;
    sys.rz[1] = 0.;

    sys.rx[2] = 0.;
    sys.ry[2] = 0.;
    sys.rz[2] = 1.5;

    /* initialize system velocities to 0 */
    azzero(sys.vx, sys.natoms);
    azzero(sys.vy, sys.natoms);
    azzero(sys.vz, sys.natoms);

    /* initialize system forces */
    sys.fx[0] = 2161.66698602;
    sys.fy[0] = -1080.83349301;
    sys.fz[0] = -1080.83349301;

    sys.fx[1] = -1080.83349301;
    sys.fy[1] = 2161.66698602;
    sys.fz[1] = -1080.83349301;

    sys.fx[2] = -1080.83349301;
    sys.fy[2] = -1080.83349301;
    sys.fz[2] = 2161.66698602;

    /* propagate system with velverlet and check */
    velverlet_2(&sys);

    fprintf(ofp, "TEST # 3a: 3 particles, 3 close to each other; v = 0\n");
    for (i=0; i < sys.natoms; ++i) {
      fprintf(ofp, "r: %20.8f%20.8f%20.8f\n", sys.rx[i], sys.ry[i], sys.rz[i]);
      fprintf(ofp, "v: %20.8f%20.8f%20.8f\n", sys.vx[i], sys.vy[i], sys.vz[i]);
    }

    /* ************************************************* */


    
    /* TEST # 3b: 3 particles, 3 close to each other; v != 0 */

    /* initialize system */

    /* initialize system positions */
    sys.rx[0] = 1.5;
    sys.ry[0] = 0.;
    sys.rz[0] = 0.;

    sys.rx[1] = 0.;
    sys.ry[1] = 1.5;
    sys.rz[1] = 0.;

    sys.rx[2] = 0.;
    sys.ry[2] = 0.;
    sys.rz[2] = 1.5;

    /* system velocities are not reinitialized */

    /* system forces do not need to be reinitialized */

    /* propagate system with velverlet and check */
    velverlet_2(&sys);

    fprintf(ofp, "TEST # 3a: 3 particles, 3 close to each other; v = 0\n");
    for (i=0; i < sys.natoms; ++i) {
      fprintf(ofp, "r: %20.8f%20.8f%20.8f\n", sys.rx[i], sys.ry[i], sys.rz[i]);
      fprintf(ofp, "v: %20.8f%20.8f%20.8f\n", sys.vx[i], sys.vy[i], sys.vz[i]);
    }

    /* ************************************************* */



    /* clean up: free memory */
    fclose(ofp);

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
