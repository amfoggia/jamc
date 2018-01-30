from ctypes import *
from math import *
import sys

class mdsys_t(Structure):
    _fields_ = [ ("natoms", c_int), ("nfi", c_int), ("nsteps", c_int), ("dt", c_double), ("redmass", c_double), ("epsilon", c_double), ("sigma", c_double), ("box", c_double), ("rcut", c_double), ("ekin", c_double), ("epot", c_double), ("temp", c_double), ("rx", POINTER(c_double)), ("ry", POINTER(c_double)), ("rz", POINTER(c_double)),("vx", POINTER(c_double)),("vy", POINTER(c_double)), ("vz", POINTER(c_double)), ("fx", POINTER(c_double)),("fy", POINTER(c_double)), ("fz", POINTER(c_double))]



def main():
    
#import DSO
    dso = CDLL("../ljmd.so")

    mdsys=mdsys_t()
    
    mvsq2e=2390.05736153349
    
#pass values for system with 3 atoms
    mdsys.natoms = 3
    mdsys.redmass = 39.948 * mvsq2e
    mdsys.epsilon = 0.2379
    mdsys.sigma = 3.405
    mdsys.rcut = 8.5
    mdsys.box = 17.1580
    mdsys.nsteps = 1
    mdsys.dt = 5.0
#other members
    mdsys.nfi=0
    mdsys.ekin=1.0
    mdsys.epot=0.0
    mdsys.temp=1.0

#allocate memory for sys.vx sys.vy sys.vz

    mdsys.rx= (c_double * mdsys.natoms)()
    mdsys.ry= (c_double * mdsys.natoms)()
    mdsys.rz= (c_double * mdsys.natoms)()
    mdsys.vx= (c_double * mdsys.natoms)()
    mdsys.vy= (c_double * mdsys.natoms)()
    mdsys.vz= (c_double * mdsys.natoms)()
    mdsys.fx= (c_double * mdsys.natoms)()
    mdsys.fy= (c_double * mdsys.natoms)()
    mdsys.fz= (c_double * mdsys.natoms)()
    
#1st CASE
#Initialize the positions and function

    mdsys.rx[0] = 0.0
    mdsys.ry[0] = 0.0
    mdsys.rz[0] = 0.0

    mdsys.rx[1] = 2.0
    mdsys.ry[1] = 2.0
    mdsys.rz[1] = 0.0

    mdsys.rx[2] = 0.0
    mdsys.ry[2] = 0.0
    mdsys.rz[2] = 8.54
    
# azzero velocity
    mdsys.vx[0]= 0.0
    mdsys.vy[0]= 0.0
    mdsys.vz[0]= 0.0
    
    mdsys.vx[1]= 0.0
    mdsys.vy[1]= 0.0
    mdsys.vz[1]= 0.0
    
    mdsys.vx[2]= 0.0
    mdsys.vy[2]= 0.0
    mdsys.vz[2]= 0.0

#forces
    mdsys.fx[0] = -22.10605666
    mdsys.fy[0] = -22.10605666
    mdsys.fz[0] = 0.0
    
    mdsys.fx[1] = -22.10605666
    mdsys.fy[1] = -22.10605666
    mdsys.fz[1] = 0.0
    
    mdsys.fx[2] = 0.0
    mdsys.fy[2] = 0.0
    mdsys.fz[2] = 0.0
    
#propagate system with velverlet
#function's call
    dso.velverlet_1.argtypes = [POINTER(mdsys_t)]
#do velverlet_1
    dso.velverlet_1(byref(mdsys))
    
#print results

    with open ("test_velverlet_1.dat", "w") as f:
        for i in range(mdsys.natoms):
            f.write ('{0} {1: >#020.8f} {2: >#020.8f} {3: >#020.8f}\n'.format(i, mdsys.rx[i], mdsys.ry[i], mdsys.rz[i]))
    with open ("test_velverlet_1.dat", "a") as f:
        for i in range(mdsys.natoms):
            f.write ('{0} {1: >#020.8f} {2: >#020.8f} {3: >#020.8f}\n'.format(i, mdsys.vx[i], mdsys.vy[i], mdsys.vz[i]))

#########################################################################################################################

#2nd CASE
#Initialize the positions and function

    mdsys.rx[0] = 0.0
    mdsys.ry[0] = 0.0
    mdsys.rz[0] = 0.0

    mdsys.rx[1] = 2.0
    mdsys.ry[1] = 2.0
    mdsys.rz[1] = 0.0

    mdsys.rx[2] = 0.0
    mdsys.ry[2] = 0.0
    mdsys.rz[2] = 8.54
    
    #propagate system with velverlet
    #function's call
    dso.velverlet_1.argtypes = [POINTER(mdsys_t)]

    #do velverlet_1
    dso.velverlet_1(byref(mdsys))
    
    #print results
    
    with open ("test_velverlet_1.dat", "w") as f:
        for i in range(mdsys.natoms):
            f.write ('{0} {1: >#020.8f} {2: >#020.8f} {3: >#020.8f}\n'.format(i, mdsys.rx[i], mdsys.ry[i], mdsys.rz[i]))
    with open ("test_velverlet_1.dat", "a") as f:
        for i in range(mdsys.natoms):
            f.write ('{0} {1: >#020.8f} {2: >#020.8f} {3: >#020.8f}\n'.format(i, mdsys.vx[i], mdsys.vy[i], mdsys.vz[i]))

if __name__ == "__main__":
    main()





