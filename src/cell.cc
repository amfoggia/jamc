#include <stdlib.h>
#include <math.h>
#include "ljmd.h"
#include <vector>
#include "cell.h"

// translate cell (i,j,k) to cel position
int index3d(mdsys_t *sys, int i, int j, int k) {
    return k+sys->cn*(j+i*sys->cn);
}

// sort particles into cells
void sort(mdsys_t *sys,cell_t *cel) {
    int i,j;
    int pos[3];
    for(i=0;i<sys->natoms;++i) {// loop through atoms
        pos[0] = floor((sys->rx[i]+0.5*sys->box)/sys->cl);
        pos[1] = floor((sys->ry[i]+0.5*sys->box)/sys->cl);
        pos[2] = floor((sys->rz[i]+0.5*sys->box)/sys->cl);
        // apply pierodic boundary condition
        for(j=0;j<3;++j) {
            while(pos[j]>(sys->cn-1)) {
                pos[j] -= sys->cn;
            }
            while(pos[j]<0) {
                pos[j] += sys->cn;
            }
        }
        int cidx = index3d(sys,pos[0],pos[1],pos[2]);// which cell
        cel[cidx].idx.push_back(i);
    }
}

