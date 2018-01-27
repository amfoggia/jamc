#include "ljmd.h"

/* main */
int main(int argc, char **argv) 
{
    FILE *traj,*erg;
    mdsys_t sys;
    const char * trajfile="test_output.xyz";
    const char * ergfile="test_output.dat";
    
    sys.natoms=1;
    sys.nfi=10;
    sys.temp=20;
    sys.ekin=30;
    sys.epot=40;
    
    

    /* allocate memory*/
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    
    sys.rx[0]=2;
    sys.ry[0]=3;
    sys.rz[0]=4;


    
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    output(&sys, erg, traj);

    fclose(erg);
    fclose(traj);
     
    
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);


    return 0;
}
