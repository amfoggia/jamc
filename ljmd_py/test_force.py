from ctypes import *
from math import *
import sys

class mdsys_t (Structure):
    _fields_ = [ ("natoms", c_int), ("nfi", c_int), ("nsteps", c_int), ("dt", c_double), ("redmass", c_double), ("epsilon", c_double), ("sigma", c_double), ("box", c_double), ("rcut", c_double), ("ekin", c_double), ("epot", c_double), ("temp", c_double), ("rx", POINTER(c_double)), ("ry", POINTER(c_double)), ("rz", POINTER(c_double)),("vx", POINTER(c_double)),("vy", POINTER(c_double)), ("vz", POINTER(c_double)), ("fx", POINTER(c_double)),("fy", POINTER(c_double)), ("fz", POINTER(c_double))]



def main():
    
#import DSO
    dso = CDLL("../ljmd.so")

    mdsys= mdsys_t()
    
    mvsq2e=2390.05736153349

#pass values for system with 3 atoms
    mdsys.natoms = 3
    mdsys.redmass = 39.948 * mvsq2e
    mdsys.epsilon = 0.2379
    mdsys.sigma = 3.405
    mdsys.rcut = 8.5
    mdsys.box = 20.0

#allocate memory for sys.vx sys.vy sys.vz

    mdsys.rx= (c_double * mdsys.natoms)()
    mdsys.ry= (c_double * mdsys.natoms)()
    mdsys.rz= (c_double * mdsys.natoms)()
    mdsys.fx= (c_double * mdsys.natoms)()
    mdsys.fy= (c_double * mdsys.natoms)()
    mdsys.fz= (c_double * mdsys.natoms)()

#1st CASE: all far apart
#Initialize the positions and function
#Particle 1
    mdsys.rx[0] = 8.51
    mdsys.ry[0] = 0.0
    mdsys.rz[0] = 0.0

#Particle 2
    mdsys.rx[1] = 0.0
    mdsys.ry[1] = 0.0
    mdsys.rz[1] = 0.0

#Particle 3
    mdsys.rx[2] = 0.0
    mdsys.ry[2] = 9.1
    mdsys.rz[2] = 0.0
    
#function's call
    dso.force.argtypes = [POINTER(mdsys_t)]

#do force
    dso.force(byref(mdsys))
    
#print results

    with open ("test_force.dat", "w") as f:
        for i in range(mdsys.natoms):
            f.write ('{0} {1: >#019.8f}{2: >#020.8f}{3: >#020.8f}\n'.format(i, mdsys.fx[i], mdsys.fy[i], mdsys.fz[i]))

############################################################################################################

#2nd CASE: two within cutoff, one outside
#Initialize the positions and function
#Particle 1
    mdsys.rx[0] = 8.51
    mdsys.ry[0] = 0.0
    mdsys.rz[0] = 0.0
    
    #Particle 2
    mdsys.rx[1] = 0.0
    mdsys.ry[1] = 4.2
    mdsys.rz[1] = 0.0
    
    #Particle 3
    mdsys.rx[2] = 0.0
    mdsys.ry[2] = 0.0
    mdsys.rz[2] = 0.0
    
    #function's call
    dso.force.argtypes = [POINTER(mdsys_t)]
   
    #do force
    dso.force(byref(mdsys))

#print results

    with open ("test_force.dat", "a") as f:
        for i in range(mdsys.natoms):
            f.write ('{0} {1: >#019.8f}{2: >#020.8f}{3: >#020.8f}\n'.format(i, mdsys.fx[i], mdsys.fy[i], mdsys.fz[i]))

##########################################################################################################
#3rd CASE: all inside cutoff
#Initialize the positions and function
#Particle 1

    mdsys.rx[0] = 1.0
    mdsys.ry[0] = 1.0
    mdsys.rz[0] = 1.0

#Particle 2
    mdsys.rx[1] = -1.5
    mdsys.ry[1] = -1.5
    mdsys.rz[1] = -1.5
    
    #Particle 3
    mdsys.rx[2] = 0.0
    mdsys.ry[2] = 0.0
    mdsys.rz[2] = 0.0
    
    #function's call
    dso.force.argtypes = [POINTER(mdsys_t)]
    
    #do force
    dso.force(byref(mdsys))

#print results

    with open ("test_force.dat", "a") as f:
        for i in range(mdsys.natoms):
            f.write ('{0} {1: >#019.8f}{2: >#020.8f}{3: >#020.8f}\n'.format(i, mdsys.fx[i], mdsys.fy[i], mdsys.fz[i]))

########################################################################################################################

#pass values for system with 4 atoms
    mdsys.natoms = 4
    mdsys.redmass = 39.948 * mvsq2e
    mdsys.epsilon = 0.2379
    mdsys.sigma = 3.405
    mdsys.rcut = 8.5
    mdsys.box = 17.1580

#allocate memory for sys.vx sys.vy sys.vz - 4 atmos

    mdsys.rx= (c_double * mdsys.natoms)()
    mdsys.ry= (c_double * mdsys.natoms)()
    mdsys.rz= (c_double * mdsys.natoms)()
    mdsys.fx= (c_double * mdsys.natoms)()
    mdsys.fy= (c_double * mdsys.natoms)()
    mdsys.fz= (c_double * mdsys.natoms)()

#Initialize the positions and function
#Particle 1
    mdsys.rx[0] = 8.51
    mdsys.ry[0] = 0.0
    mdsys.rz[0] = 0.0

#Particle 2
    mdsys.rx[1] = 0.0
    mdsys.ry[1] = 9.3
    mdsys.rz[1] = 0.0

#Particle 3
    mdsys.rx[2] = 0.0
    mdsys.ry[2] = 0.0
    mdsys.rz[2] = 8.98

#Particle 4
    mdsys.rx[3] = 0.0
    mdsys.ry[3] = 0.0
    mdsys.rz[3] = 0.0

#function's call
    dso.force.argtypes = [POINTER(mdsys_t)]
    
    #do force
    py_force = dso.force(byref(mdsys))

#print results

    with open ("test_force.dat", "a") as f:
        for i in range(mdsys.natoms):
            f.write ('{0} {1: >#019.8f}{2: >#020.8f}{3: >#020.8f}\n'.format(i, mdsys.fx[i], mdsys.fy[i], mdsys.fz[i]))

##############################################################################################


if __name__ == "__main__":
    main()





