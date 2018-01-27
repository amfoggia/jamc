// reference output of this test should be test_ekin.dat
#include "ljmd.h"
#include <stdlib.h>

/* test main */
int main(int argc, char **argv){
    int i;
    mdsys_t sys;
    
    sys.natoms=10;
    sys.redmass = 39.948*2390.05736153349; // reduced mass
    //sys.mass = 39.948;
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    
    for (i=0; i<sys.natoms; ++i) {
        *(sys.vx+i) = 1e-3;
        *(sys.vy+i) = 1e-3;
        *(sys.vz+i) = 1e-3;
    }
    
    ekin(&sys);
    
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    
    double expected = 10*0.5*sys.redmass*(3.e-6);
    if(sys.ekin==expected) {
        return 0;
    }
    else {
        return 1;
    }
}
