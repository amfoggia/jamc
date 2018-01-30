from ctypes import *
from math import *
import sys

class mdsys_t(Structure):
    _fields_ = [ ("natoms", c_int), ("nfi", c_int), ("nsteps", c_int), ("dt", c_double), ("redmass", c_double), ("epsilon", c_double), ("sigma", c_double), ("box", c_double), ("rcut", c_double), ("ekin", c_double), ("epot", c_double), ("temp", c_double), ("rx", POINTER(c_double)), ("ry", POINTER(c_double)), ("rz", POINTER(c_double)),("vx", POINTER(c_double)),("vy", POINTER(c_double)), ("vz", POINTER(c_double)), ("fx", POINTER(c_double)),("fy", POINTER(c_double)), ("fz", POINTER(c_double))]

def main():

    mdsys= mdsys_t()


#import DSO
    dso = CDLL("../ljmd.so")

#pass size
    mdsys.natoms=10
    
#pass constants
    mdsys.redmass = 39.948*2390.05736153349

#pass value for verify
    expected = 10*0.5*mdsys.redmass*(3.e-6)

#allocate memory for sys.vx sys.vy sys.vz

    mdsys.vx= (c_double * mdsys.natoms)()
    mdsys.vy= (c_double * mdsys.natoms)()
    mdsys.vz= (c_double * mdsys.natoms)()

#initializze sys.vx sys.vy sys.vz

    for i in range(mdsys.natoms):
        mdsys.vx[i] = 1e-3
        mdsys.vy[i] = 1e-3
        mdsys.vz[i] = 1e-3

#function's call
    dso.ekin.argtypes = [POINTER(mdsys_t)]

#do ekin
    dso.ekin(byref(mdsys))

    if(mdsys.ekin==expected):
        
        return 0
    else:
        exit
if __name__ == "__main__":
    main()

