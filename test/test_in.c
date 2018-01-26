#include "../src/ljmd.h"
#include <string.h>

int main(int argc, char **argv)
{
    int nprint;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    
    mdsys_t sys;
    
    char * atom= "2";
    char * rest= "rs";
    
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
    double test_mass, test_epsilon, test_sigma, test_rcut, test_box, test_dt;
    char * test_rest , * test_traject, * test_erg;
    
    test_atom= 2;
    if (sys.natoms != test_atom) return 1;
    test_mass= 3;
    if (sys.mass != test_mass) return 1;
    test_epsilon= 4;
    if (sys.epsilon != test_epsilon) return 1;
    test_sigma= 5;
    if (sys.sigma != test_sigma) return 1;
    test_rcut= 6;
    if (sys.rcut != test_rcut) return 1;
    test_box= 7;
    if (sys.box != test_box) return 1;
    test_rest= "rs";
    if (strcmp (restfile, test_rest)!= 0) return 1;
    test_traject= "tr";
    if (strcmp (trajfile, test_traject) !=0) return 1;
    test_erg= "er";
    if (strcmp (ergfile, test_erg) !=0) return 1;
    test_steps= 8;
    if (sys.nsteps != test_steps) return 1;
    test_dt= 9;
    if (sys.dt != test_dt) return 1;
    test_nprint= 10;
    if (nprint != test_nprint) return 1;
    
    
    return 0;
}

