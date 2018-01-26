#include "../src/ljmd.h"
#include <string.h>


int main(int argc, char **argv)
{
    int nprint;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE * fn;
    mdsys_t sys;
    
    char * atom= "108";
    char * rest= "argon_108.rest";
    
    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    if (strcmp (line, atom)!= 0) return 1;
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
    if (strcmp (restfile, rest)!= 0) return 1;
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    nprint=atoi(line);
    
    int test_atom, test_steps, test_nprint;
    double test_mass=sys.mass, test_epsilon, test_sigma, test_rcut, test_box, test_dt;
    char * test_rest , * test_traject=trajfile, * test_erg;
    
    fn= fopen("test_in.dat", "w");
    
     test_atom= 108;
     if (sys.natoms != test_atom) return 1;
     else fprintf(fn, "%d\n",sys.natoms );
     test_mass= 39.948;
     if (sys.mass != test_mass) return 1;
     else fprintf(fn, "%g\n", sys.mass);
     test_epsilon= 0.2379;
     if (sys.epsilon != test_epsilon) return 1;
     else fprintf(fn, "%g\n", sys.epsilon);
     test_sigma= 3.405;
     if (sys.sigma != test_sigma) return 1;
     else fprintf(fn, "%g\n",sys.sigma);
     test_rcut= 8.5;
     if (sys.rcut != test_rcut) return 1;
     else fprintf(fn, "%g\n",sys.rcut);
     test_box= 17.1580;
     if (sys.box != test_box) return 1;
     else fprintf(fn, "%g\n",sys.box);
     test_rest= "argon_108.rest";
     if (strcmp (restfile, test_rest)!= 0) return 1;
     else fprintf(fn, "%s\n",restfile);
     test_traject= "argon_108.xyz";
     if (strcmp (trajfile, test_traject) !=0) return 1;
     else fprintf(fn, "%s\n",trajfile);
     test_erg= "argon_108.dat";
     if (strcmp (ergfile, test_erg) !=0) return 1;
     else fprintf(fn, "%s\n",ergfile);
     test_steps= 10000;
     if (sys.nsteps != test_steps) return 1;
     else fprintf(fn, "%d\n",sys.nsteps);
     test_dt= 5.0;
     if (sys.dt != test_dt) return 1;
     else fprintf(fn, "%g\n",sys.dt);
     test_nprint= 100;
     if (nprint != test_nprint) return 1;
     else fprintf(fn, "%d\n",nprint);
    
     fclose(fn);
    
    return 0;
}




