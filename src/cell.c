#include <stdlib.h>
#include <math.h>
#include <vector>
#include "ljmd.h"

// sort particles into cells
void sort(mdsys_t *sys,cell_t *cel) {
    for(int it=0;it<sys->cn*sys->cn*sys->cn;++it) {
        cel[it].idx.clear();// clean idx list of each cell
    }
    int pos[3];
    for(int i=0;i<sys->natoms;++i) {// loop through atoms
        pos[0] = floor((sys->rx[i]+0.5*sys->box)/sys->cl);
        pos[1] = floor((sys->ry[i]+0.5*sys->box)/sys->cl);
        pos[2] = floor((sys->rz[i]+0.5*sys->box)/sys->cl);
        // apply pierodic boundary condition
        int cidx = index3d(sys,pos[0],pos[1],pos[2]);// which cell
        cel[cidx].idx.push_back(i);
    }
}

/* aggressive cell list */
void cell_force(mdsys_t *sys,cell_t *cell){
    double rsq,rsq_inv,r6,ffac;
    double rx,ry,rz;
    // zero energy and forces
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    double boxby2 = 0.5*sys->box;// pre-calculate
    double rcutsq = sys->rcut*sys->rcut;// pre-calculate, take square
    double c6 = sys->epsilon*pow(sys->sigma,6);
    double c12 = sys->epsilon*pow(sys->sigma,12);

    // LOOP THROUGH CELL PAIRS
    for(unsigned it=0;it<sys->pair.size();it+=2) {
        // loop through cell pairs
        int c1 = sys->pair[it];
        int c2 = sys->pair[it+1];
        for(unsigned e1=0;e1<cell[c1].idx.size();++e1){
            // element in cell 1
            int i = cell[c1].idx[e1];
	    unsigned e2 = (c1==c2)? e1+1 : 0;
            for(;e2<cell[c2].idx.size();++e2){
                // element in cell 2
                int j = cell[c2].idx[e2];
                if(i==j) continue;
                // get distance between particle i and j
                rx=pbc(sys->rx[i] - sys->rx[j], boxby2);
                rsq = rx*rx;
                if(rsq>rcutsq) continue;
                ry=pbc(sys->ry[i] - sys->ry[j], boxby2);
                rsq += ry*ry;
                if(rsq>rcutsq) continue;
                rz=pbc(sys->rz[i] - sys->rz[j], boxby2);
                rsq += rz*rz;
                if(rsq>rcutsq) continue;
                
                rsq_inv = 1.0/rsq;
                r6 = rsq_inv*rsq_inv*rsq_inv;
                ffac = (48*c12*r6-24*c6)*r6*rsq_inv;
                sys->epot += 4*r6*(c12*r6-c6);
                sys->fx[i] += rx*ffac; sys->fx[j] -= rx*ffac;
                sys->fy[i] += ry*ffac; sys->fy[j] -= ry*ffac;
                sys->fz[i] += rz*ffac; sys->fz[j] -= rz*ffac;
            }// c1
        }// c2
    }// cell pairs
}
