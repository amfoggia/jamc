#include "../src/ljmd.h"

int main(int argc, char **argv)
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    
    sys.nfi=4;
    sys.temp=4.;
    sys.ekin=4.;
    sys.epot=4.;
    
    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    sys.natoms=atoi(line);
    //this fpritntf is checking the get_a_line function
    fprintf("simple input check/n", line)
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
    if(get_a_line(stdin,line)) return 1;
    nprint=atoi(line);
    
    check_1=fopen("test_inout.dat", "w")
    fprintf(")
    
    //allocate memory
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
   
    
    //read restart
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        
    
        fclose(fp);
        
    }
    
        else {
        perror("cannot read restart file");
        return 3;
            
        }
    
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");


    output(&sys, erg,traj);
    
    /* clean up: close files, free memory */
    printf("Simulation Done.\n");
    
    fclose(erg);
    fclose(traj);
    
    free(sys.rx);
    free(sys.ry);
    free(sys.rz);

    
    return 0;
}

