#include <sys/time.h>
#include <vector>
#include "ljmd.h"

/* helper function: zero out an array */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* compute kinetic energy */
void ekin(mdsys_t *sys)
{   
    int i;
    
    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*sys->redmass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}

double pbc(double x, const double boxby2) {
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

/* timing */
double stamp(){
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*1e+3 + tv.tv_usec*1e-3;
}

void pair(mdsys_t *sys) {
    int n = sys->cn;
    int total_cell = n*n*n;
    for(int c1=0;c1<total_cell;++c1){
      sys->pair.push_back(c1);
      sys->pair.push_back(c1);
        int k1 = c1%(n);
        int j1 = ((c1-k1)/n)%n;
        int i1 = (c1-j1*n-k1)/(n*n);
        for(int c2=c1+1;c2<total_cell;++c2){
            // calculate distance between two cells
            int k2 = c2%(n);
            int j2 = ((c2-k2)/n)%n;
            int i2 = (c2-j2*n-k2)/(n*n);
            // exclude rule, 26 different periodic mirroring cases
            int dk = k1-k2; int dj = j1-j2; int di = i1-i2;
            int dist = dk*dk+dj*dj+di*di;
            if(dist<4){
                sys->pair.push_back(c1);
                sys->pair.push_back(c2);
            }
            // pierodic in i, 2 cases
            else if (dist+n*n+2*n*di<4 || dist+n*n-2*n*di<4){
                sys->pair.push_back(c1);
                sys->pair.push_back(c2);
            }
            // pierodic in j, 2 cases
            else if (dist+n*n+2*n*dj<4 || dist+n*n-2*n*dj<4){
                sys->pair.push_back(c1);
                sys->pair.push_back(c2);
            }
            // pierodic in k, 2 cases
            else if (dist+n*n+2*n*dk<4 || dist+n*n-2*n*dk<4){
                sys->pair.push_back(c1);
                sys->pair.push_back(c2);
            }
            // pierodic in i,j, 4 cases
            else if (dist+2*n*n+2*n*(di+dj)<4 || dist+2*n*n+2*n*(di-dj)<4 || dist+2*n*n+2*n*(dj-di)<4|| dist+2*n*n-2*n*(di+dj)<4){
                sys->pair.push_back(c1);
                sys->pair.push_back(c2);
            }
            // pierodic in i,k, 4 cases
            else if (dist+2*n*n+2*n*(di+dk)<4 || dist+2*n*n+2*n*(di-dk)<4 || dist+2*n*n+2*n*(dk-di)<4|| dist+2*n*n-2*n*(di+dk)<4){
                sys->pair.push_back(c1);
                sys->pair.push_back(c2);
            }
            // pierodic in j,k, 4 cases
            else if (dist+2*n*n+2*n*(dk+dj)<4 || dist+2*n*n+2*n*(dk-dj)<4 || dist+2*n*n+2*n*(dj-dk)<4|| dist+2*n*n-2*n*(dk+dj)<4){
                sys->pair.push_back(c1);
                sys->pair.push_back(c2);
            }
            // peirodic in ijk, 8 cases
            else if (dist+3*n*n+2*n*(dk+dj+di)<4 || dist+3*n*n+2*n*(dk+dj-di)<4 || dist+3*n*n+2*n*(-dk+dj+di)<4 || dist+3*n*n+2*n*(dk-dj+di)<4 || dist+3*n*n+2*n*(-dk-dj+di)<4 || dist+3*n*n+2*n*(-dk+dj-di)<4 || dist+3*n*n+2*n*(dk-dj-di)<4 || dist+3*n*n-2*n*(dk+dj+di)<4){
                sys->pair.push_back(c1);
                sys->pair.push_back(c2);
            }
        }
    }
}


// translate cell (i,j,k) to cel position but with boundary fixing
int index3d(mdsys_t *sys, int i, int j, int k) {
  i %= sys->cn;
  j %= sys->cn;
  k %= sys->cn;
  if(i<0) i += sys->cn;
  if(j<0) j += sys->cn;
  if(k<0) k += sys->cn;
  return k+sys->cn*(j+i*sys->cn);
}
