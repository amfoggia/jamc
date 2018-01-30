/*
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * baseline c version.
 */

#ifndef JLMD_H
#define JLMD_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <ctype.h>
#ifdef USE_MPI
#include <mpi.h>
#endif

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
#define kboltz 0.0019872067     /* boltzman constant in kcal/mol/K */
#define mvsq2e 2390.05736153349 /* m*v^2 in kcal/mol */

/* structure to hold the complete information 
 * about the MD system */

extern "C"{

struct _mdsys {
    int natoms,nfi,nsteps;
    double dt, redmass, epsilon, sigma, box, rcut;
    double ekin, epot, temp;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    
    // for cells *****
    int cn;
    double cl;
    std::vector<int> pair;
    //*******
};

typedef struct _mdsys mdsys_t;

/* helper function: read a line and then return
   the first string with whitespace stripped off */
int get_a_line(FILE *fp, char *buf);

/* helper function: zero out an array */
void azzero(double *d, const int n);

/* helper function: apply minimum image convention */
/*
static inline double pbc(double x, const double boxby2) {
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}
 */
double pbc(double x, const double boxby2); // for profiling

/* compute kinetic energy */
void ekin(mdsys_t *sys);

/* compute forces */
void force(mdsys_t *sys);
void force_index_array(mdsys_t *sys);
void force_atomic(mdsys_t *sys);
void force_basic(mdsys_t *sys);

/* velocity verlet */
void velverlet_1(mdsys_t *sys);
void velverlet_2(mdsys_t *sys);

/* append data to output. */
void output(mdsys_t *sys, FILE *erg, FILE *traj);

/* timer in ms */
double stamp();


//*********************** CELL *********

// cell structure
struct _mdcell {
    std::vector<int> idx;
};

typedef struct _mdcell cell_t;

/* pair up cell indecies */
void pair(mdsys_t *sys);

/* 3D->1D index translator */
int index3d(mdsys_t *sys, int i, int j, int k);

/* place elements in cells */
void sort(mdsys_t *sys, cell_t *cel);

/* calculate force with cell list */
void cell_force(mdsys_t *sys, cell_t *cel);

#endif
}
