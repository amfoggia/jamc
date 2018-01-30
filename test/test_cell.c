#include <stdlib.h>
#include "ljmd.h"

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
    pair(&sys);
    
    int pair_count=0;
    for(int it=0;it<sys.pair.size();it+=2){
        pair_count+=1;
    }
    // non-repeatable Permutation of 27 to 2 should be 351
    if(pair_count!=351+27) exit(1);
    
    sort(&sys,cel);
    
    sys.rx[0] = -6; sys.ry[0] = 0; sys.rz[0] = 0;
    
    sort(&sys,cel);
    
    int count = 0;
    for(int it=0;it<sys.cn*sys.cn*sys.cn;++it){
        count += cel[it].idx.size();
        //printf("cell id %d, load: %lu\n",it,cel[it].idx.size());
    }
    if(count != sys.natoms) {
        printf("wrong cell load\n");
        exit(1);
    }
    
    if(cel[4].idx[0]==1 && cel[22].idx[0]==0 && cel[22].idx[1]==2 && cel[26].idx.size()==0){
        free(sys.rx);free(sys.ry);free(sys.rz); delete [] cel;
        return 0;
    }
    
    free(sys.rx);free(sys.ry);free(sys.rz); delete [] cel;
    printf("wrong allocation\n");
    exit(1);
}
