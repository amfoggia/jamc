# jamc report
P1.6 group assignment (GROUP 1): Lennard-Jones Molecular Dynamics

## Collaborators:

             - Jiaxin Wang         ---> GitHub: ricphy 
             
             - Marco Bettiol       ---> GitHub: marcob3t
             
             - Carolina Bonivento  ---> GitHub: carolinabonivento
             
             - Alejandra Foggia    ---> GitHub: amfoggia

## Contributions:

#### Jiaxin Wang:
* unit test for kinetic energy
* calculation optimization

#### Marco Bettiol:
* unit test for integration (velvervelt)
* multi-threading

#### Carolina Bonivento:
* unit test for input/output
* python interface

#### Alejandra Foggia:
* unit test for force calculation
* MPI

## Remarks:

#### optimization testing log (clang, i7, 3.1GHz):
* force function, serial, testing with argon_2916

|timing|feature|
|------|-------|
|120-130 ms|original|
|115-120 ms|+ 1D pre truncate|
|110-120 ms|+ use r-square|
|95-105  ms|+ inline pbc|
|93-100  ms|+ replace i==j|
|93-100  ms|+ rsq_inv|
(in old gcc std, we can only use "static inline")

* velverlet function, serial, testing with argon_2916

|time(consider looping)|feature|
|----------------------|-------|
|~95 s|original|
|~94 s|+ reduced mass|
|~80 s|+ tearing down spherical truncate into 1D aggressive truncate|

(the sys.mass is always used as sys.mass*mvq2e, 1s would be negligible in single simulation,
but would be meaningful in, ie., Bayesian analysis)

1D aggressive truncate means:

```
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
```

This package contains simplified MD code with multi-threading
parallelization for simulating atoms with a Lennard-Jones potential.

The bundled makefiles are set up to compile the executable once
with OpenMP disabled and once with OpenMP enabled with each build
placing the various object files in separate directories.

The examples directory contains 3 sets of example input decks
and the reference directory the corresponding outputs.

Type: make
to compile everything and: make clean
to remove all compiled objects
