// reference output of this test should be test_ekin.dat
#include "ljmd.h"
#include <stdlib.h>

/* test main */
int main(int argc, char **argv){
    int i;
    FILE *fn;
    mdsys_t sys;
    
    sys.natoms=10;
    //sys.redmass = 39.948*2390.05736153349; // reduced mass
    sys.mass = 39.948;
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    
    for (i=0; i<sys.natoms; ++i) {
        *(sys.vx+i) = 1e-3*(rand()%10);
        *(sys.vy+i) = 1e-3*(rand()%10);
        *(sys.vz+i) = 1e-3*(rand()%10);
    }
    
    ekin(&sys);
    //fn = fopen("../reference/test_ekin.dat","w"); // only for first run
    fn = fopen("test_ekin.dat","w");
    // print out to reference file
    fprintf(fn,"%.8f\n",sys.ekin);
    
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    return 0;
}
