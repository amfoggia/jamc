# Script to plot scaling
reset
## omp agg --------------------------------
set title 'LJMD with OpenMP - argon 2916 (aggressive truncation)'
set xlabel 'omp threads'
set ylabel 'speedup'
set grid ytics
set grid xtics
set key left

plot 'omp_agg.txt' u 3:2 w lp pt 7 ps 1.4 title '', x lt 7 title 'reference'

set terminal png
set output 'omp_agg.png'
rep
unset terminal
reset
## omp aggNewtonAtomic --------------------------
set title 'LJMD with OpenMP - argon 2916 (Newton + Atomic)'
set xlabel 'omp threads'
set ylabel 'speedup'
set grid ytics
set grid xtics
set key left

plot 'omp_aggAtomic.txt' u 3:2 w lp pt 7 ps 1.4 title '', x lt 7 title 'reference'

set terminal png
set output 'omp_aggAtomic.png'
rep
unset terminal
reset
## omp aggNewtonReplicatedMemLINLOOPSIndexArray --------------------------
set title 'LJMD with OpenMP - argon 2916 (Newton + Linearized Loops + IndexArray)'
set xlabel 'omp threads'
set ylabel 'speedup'
set grid ytics
set grid xtics
set key left

plot 'omp_indexArray.txt' u 3:2 w lp pt 7 ps 1.4 title '', x lt 7 title 'reference'

set terminal png
set output 'omp_IndexArray.png'
rep
unset terminal
reset
## omp aggNewtonReplicatedMemLinearLoopInPlaceIndex --------------------------
set title 'LJMD with OpenMP - argon 2916 (Newton + Linearized loops + In-Place index array)'
set xlabel 'omp threads'
set ylabel 'speedup'
set grid ytics
set grid xtics
set key left

plot 'omp_inPlace.txt' u 3:2 w lp pt 7 ps 1.4 title '', x lt 7 title 'reference'

set terminal png
set output 'omp_InPlace.png'
rep
unset terminal
reset
## omp all together ----------------------------------------------------------
set title 'LJMD with OpenMP - argon 2916'
set xlabel 'omp threads'
set ylabel 'speedup'
set key left
set grid ytics
set grid xtics

plot 'omp_agg.txt' u 3:2 w lp pt 7 ps 1.4 title 'agg trunc', \
     'omp_aggAtomic.txt' u 3:2 w lp pt 7 ps 1.4 title 'Newt+Atomic', \
     'omp_indexArray.txt' u 3:2 w lp pt 7 ps 1.4 title 'Newt+Lin Loops+Index', \
     'omp_inPlace.txt' u 3:2 w lp pt 7 ps 1.4 title 'Newt+Lin loops+In-Place index', \
     x lt 7 title 'reference'

set terminal png
set output 'omp_all.png'
rep
unset terminal
reset
## mpi ------------------------------------
set title 'LJMD with MPI - argon 2916'
set xlabel 'MPI processes'
set ylabel 'speedup'
set grid ytics
set grid xtics
set key left

plot 'mpi.txt' u 3:2 w lp pt 7 ps 1.4 title '', x lt 7 title 'reference'

set terminal png
set output 'mpi.png'
rep
unset terminal
reset

## omp all together ----------------------------------------------------------
set title 'LJMD with OpenMP - argon 2916'
set xlabel 'omp threads'
set ylabel 'speedup'
set key left
set grid ytics
set grid xtics

plot 'omp_agg.txt' u 3:2 w lp pt 7 ps 1.4 title 'agg trunc', \
     'omp_aggAtomic.txt' u 3:2 w lp pt 7 ps 1.4 title 'Newt+Atomic', \
     'omp_indexArray.txt' u 3:2 w lp pt 7 ps 1.4 title 'Newt+Lin Loops+Index', \
     'omp_inPlace.txt' u 3:2 w lp pt 7 ps 1.4 title 'Newt+Lin loops+In-Place index', \
     'mpi.txt' u 3:2 w lp pt 7 ps 1.4 title 'mpi', \
     x lt 7 title 'reference'

set terminal png
set output 'ompMpi_all.png'
rep
unset terminal
reset