#include <stdlib.h>
#include "ljmd.h"
#include "cell.h"

int main(int argc, char **argv) {
    mdsys_t sys;
    
    sys.natoms=100;
    sys.box=10;
    sys.rcut=3;
    // get number of cells in 1D
    sys.cn = floor(sys.box/sys.rcut);
    if(sys.rcut!=3) {printf("wrong sys.cn");return 1;}
    sys.cl = sys.box/sys.cn;
    
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    
    azzero(sys.rx,sys.natoms);
    azzero(sys.ry,sys.natoms);
    azzero(sys.rz,sys.natoms);
    
    sys.rx[0] = 2; sys.ry[0] = 2; sys.rz[0] = 2;
    sys.rx[1] = 6; sys.ry[1] = 0; sys.rz[1] = 1;
    sys.rx[2] = -6; sys.ry[1] = 0; sys.rz[1] = 0;
    
    // get number of cells in 1D
    sys.cn = floor(sys.box/sys.rcut);
    sys.cl = sys.box/sys.cn;
    cell_t* cel = new cell_t[(sys.cn)*(sys.cn)*(sys.cn)];
    sort(&sys,cel);
    
    if(cel[4].idx[0]==1 && cel[22].idx[0]==2 && cel[26].idx[0]==0){
        free(sys.rx);free(sys.ry);free(sys.rz);delete [] cel;
        return 0;
    }
    free(sys.rx);free(sys.ry);free(sys.rz);delete [] cel;
    return 1;
}
