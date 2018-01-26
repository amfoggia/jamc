# jamc
P1.6 group assignment (GROUP 1): Lennard-Jones Molecular Dynamics

Collaborators:

             - Jiaxin Wang         ---> GitHub: ricphy 
             
             - Marco Bettiol       ---> GitHub: marcob3t
             
             - Carolina Bonivento  ---> GitHub: carolinabonivento
             
             - Alejandra Foggia    ---> GitHub: amfoggia

Contributions:

Jiaxin Wang:
    unit test for kinetic energy
    calculation optimization

Marco Bettiol:
    unit test for integration (velvervelt)
    multi-threading

Carolina Bonivento:
    unit test for input/output
    python interface

Alejandra Foggia:
    unit test for force calculation
    MPI

Remarks:





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
