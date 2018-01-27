#include "ljmd.h"

/* compute forces */
void force(mdsys_t *sys)
{
    double rsq,rsq_inv,r6,ffac;
    double rx,ry,rz;
    int i,j;
    
    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    double boxby2 = 0.5*sys->box;// pre-calculate
    double rcutsq = sys->rcut*sys->rcut;// pre-calculate, take square
    double c6 = sys->epsilon*pow(sys->sigma,6);
    double c12 = sys->epsilon*pow(sys->sigma,12);
    
    for(i=0; i < (sys->natoms); ++i) {
        for(j=0; j < i; ++j) {
            /* particles have no interactions with themselves */
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], boxby2);
            rsq = rx*rx;
            if(rsq>rcutsq) continue; // 1D pre truncate
            ry=pbc(sys->ry[i] - sys->ry[j], boxby2);
            rsq += ry*ry;
            if(rsq>rcutsq) continue; // 1D pre truncate
            rz=pbc(sys->rz[i] - sys->rz[j], boxby2);
            rsq += rz*rz;
            if(rsq>rcutsq) continue; // 1D pre truncate
            
            rsq = rx*rx + ry*ry + rz*rz;
            rsq_inv = 1.0/rsq;
            r6 = rsq_inv*rsq_inv*rsq_inv;
            ffac = (48*c12*r6-24*c6)*r6*rsq_inv;
            sys->epot += 2*r6*(c12*r6-c6);
            sys->fx[i] += rx*ffac;
            sys->fy[i] += ry*ffac;
            sys->fz[i] += rz*ffac;
            
        }
        /* particles have no interactions with themselves */
        // rm (i==j) decision with separated loop
        for(j=i+1; j < sys->natoms; ++j) {
            /* particles have no interactions with themselves */
            
            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], boxby2);
            rsq = rx*rx;
            if(rsq>rcutsq) continue; // 1D pre truncate
            ry=pbc(sys->ry[i] - sys->ry[j], boxby2);
            rsq += ry*ry;
            if(rsq>rcutsq) continue; // 1D pre truncate
            rz=pbc(sys->rz[i] - sys->rz[j], boxby2);
            rsq += rz*rz;
            if(rsq>rcutsq) continue; // 1D pre truncate
            
            rsq = rx*rx + ry*ry + rz*rz;
            rsq_inv = 1.0/rsq;
            r6 = rsq_inv*rsq_inv*rsq_inv;
            ffac = (48*c12*r6-24*c6)*r6*rsq_inv;
            sys->epot += 2*r6*(c12*r6-c6);
            sys->fx[i] += rx*ffac;
            sys->fy[i] += ry*ffac;
            sys->fz[i] += rz*ffac;
            
        }
    }
}
