# jamc report
P1.6 group assignment (GROUP 1): Lennard-Jones Molecular Dynamics

## Collaborators:

             - Jiaxin Wang         ---> GitHub: ricphy 
             
             - Marco Bettiol       ---> GitHub: marcob3t
             
             - Carolina Bonivento  ---> GitHub: carolinabonivento
             
             - Alejandra Foggia    ---> GitHub: amfoggia

## Contributions:

#### Jiaxin Wang:
* unit test for kinetic energy: [test_ekin](./test/test_ekin.c)
* optimization: report below in "OPT log" section
* cell list module: [cell](./src/cell.c)
* unit test for cell list: [test_cell](./test/test_cell.c)

#### Marco Bettiol:
* unit test for integration: [test_velverlet_1](./test/test_velverlet_1.c) [test_velverlet_2](./test/test_velverlet_2.c)
* multi-threading (including benchmarking): [force](./src/force.c)

#### Carolina Bonivento:
* unit test for input/output: [test_in](./test/test_in.c) [test_out](./test/test_out.c)
* python interface:

#### Alejandra Foggia:
* unit test for force calculation: [test_force](./test/test_force.c)
* MPI: [ljmd_mpi](./src/ljmd_mpi.c.old) (already integrated in ljdm.c)
* apply cell list: [cell](./src/cell.c) [cell_aux](./src/utilities.c)


## OPT log:
(related source files in [here](./opt/) )

* profiling original code

```
Each sample counts as 0.01 seconds.
%   cumulative   self              self     total
time   seconds   seconds    calls  us/call  us/call  name
68.75      4.85     4.85    10001   485.34   684.99  force(_mdsys*)
28.00      6.83     1.98 346714668     0.01     0.01  pbc(double, double)
0.57      6.87     0.04    10001     4.00     4.00  ekin(_mdsys*)
0.57      6.91     0.04    10000     4.00     4.00  velverlet_1(_mdsys*)
0.28      6.93     0.02    30006     0.67     0.67  azzero(double*, int)
0.07      6.94     0.01    10000     0.50     0.50  velverlet_2(_mdsys*)
0.00      6.94     0.00      101     0.00     0.00  output(_mdsys*, _IO_FILE*, _IO_FILE*)
0.00      6.94     0.00       12     0.00     0.00  get_a_line(_IO_FILE*, char*)
```

* serial code optimization without cell list

aggressive truncation is better than spherical truncation:

```
/* get distance between particle i and j */
rx=pbc(sys->rx[i] - sys->rx[j], boxby2);
rsq = rx*rx;
if(rsq>rcutsq) continue; 
ry=pbc(sys->ry[i] - sys->ry[j], boxby2);
rsq += ry*ry;
if(rsq>rcutsq) continue; 
rz=pbc(sys->rz[i] - sys->rz[j], boxby2);
rsq += rz*rz;
if(rsq>rcutsq) continue; 
```
which is a temporary solution before implementing cell-list

profiling with using aggressive truncatioin

```
Each sample counts as 0.01 seconds.
%   cumulative   self              self     total
time   seconds   seconds    calls  us/call  us/call  name
73.49      5.50     5.50    10001   550.36   737.48  force(_mdsys*)
24.99      7.38     1.87 319175220     0.01     0.01  pbc(double, double)
1.00      7.45     0.08    10001     7.50     7.50  ekin(_mdsys*)
0.33      7.48     0.03                             stamp()
0.13      7.49     0.01    10000     1.00     1.00  velverlet_1(_mdsys*)
0.13      7.50     0.01    10000     1.00     1.00  velverlet_2(_mdsys*)
0.00      7.50     0.00    30006     0.00     0.00  azzero(double*, int)
0.00      7.50     0.00      101     0.00     0.00  output(_mdsys*, _IO_FILE*, _IO_FILE*)
0.00      7.50     0.00       12     0.00     0.00  get_a_line(_IO_FILE*, char*)
```

Newton 3rd law:

```
for(i=0; i < (sys->natoms); ++i) {
    for(j=i+1; j < (sys->natoms); ++j) {
        // get distance between particle i and j
        rx=pbc(sys->rx[i] - sys->rx[j], boxby2);
        rsq = rx*rx;
        if(rsq>rcutsq) continue;
        ry=pbc(sys->ry[i] - sys->ry[j], boxby2);
        rsq += ry*ry;
        if(rsq>rcutsq) continue;
        rz=pbc(sys->rz[i] - sys->rz[j], boxby2);
        rsq += rz*rz;
        if(rsq>rcutsq) continue;

        rsq_inv = 1.0/rsq;
        r6 = rsq_inv*rsq_inv*rsq_inv;
        ffac = (48*c12*r6-24*c6)*r6*rsq_inv;
        sys->epot += 4*r6*(c12*r6-c6);
        sys->fx[i] += rx*ffac; sys->fx[j] -= rx*ffac;
        sys->fy[i] += ry*ffac; sys->fy[j] -= ry*ffac;
        sys->fz[i] += rz*ffac; sys->fz[j] -= rz*ffac;
    }
}
```

profiling with Newton 3rd law and aggressive truncation

```
Each sample counts as 0.01 seconds.
%   cumulative   self              self     total
time   seconds   seconds    calls  us/call  us/call  name
73.59      3.00     3.00    10001   300.21   401.78  force(_mdsys*)
24.65      4.01     1.01 159587610     0.01     0.01  pbc(double, double)
1.23      4.06     0.05    10001     5.00     5.00  ekin(_mdsys*)
0.25      4.07     0.01    30006     0.33     0.33  azzero(double*, int)
0.25      4.08     0.01    10000     1.00     1.00  velverlet_1(_mdsys*)
0.12      4.08     0.01                             stamp()
0.00      4.08     0.00    10000     0.00     0.00  velverlet_2(_mdsys*)
0.00      4.08     0.00      101     0.00     0.00  output(_mdsys*, _IO_FILE*, _IO_FILE*)
0.00      4.08     0.00       12     0.00     0.00  get_a_line(_IO_FILE*, char*)
```

force function accumulated timing:

|argon_108 (ms)|argon_2916 (ms)|feature|
|--------------|---------------|-------|
|2999|99288|original|
|2477|76481|+agg. trunc|
|1384|42099|+Newton|
|1198|37016|+Newton +agg. trunc|

comments: aggressive truncation is benefitial, so we keep it as standard *until* cell list

* compiling flags

```
-ffast-math -O3
```
are the flags can give significant boost to computing time, '-O3' gives almost 2x

others like
```
-msse3 -fexpensive-optimizations
```
give only a slight/negligible improvement


* OMP and MPI performances

force function accumulated timing:

|argon_2916 (ms)|speedup|OMP_threads|feature|
|--------------|--------|---------------------|-------|
|73217|1.00|1|+agg. trunc|
|37895|1.93|2|+agg. trunc|
|26387|2.77|3|+agg. trunc|
|20478|3.58|4|+agg. trunc|
|16802|4.36|5|+agg. trunc|
|13997|5.23|6|+agg. trunc|
|12071|6.07|7|+agg. trunc|
|10585|6.92|8|+agg. trunc|
|9361|7.82|9|+agg. trunc|
|8516|8.60|10|+agg. trunc|

ulysses cluster, full node, no binding to socket, 5 meas. sample, err O(10), first taken as best serial

|argon_2916 (ms)|speedup|OMP_threads|feature|
|--------------|--------|---------------------|-------|
|48112|1.00|1|+agg. trunc newt(atomic)|
|36787|1.31|2|+agg. trunc newt(atomic)|
|29179|1.65|3|+agg. trunc newt(atomic)|
|25873|1.86|4|+agg. trunc newt(atomic)|
|21927|2.19|5|+agg. trunc newt(atomic)|
|19927|2.41|6|+agg. trunc newt(atomic)|
|17925|2.68|7|+agg. trunc newt(atomic)|
|16813|2.86|8|+agg. trunc newt(atomic)|
|15510|3.10|9|+agg. trunc newt(atomic)|
|14796|3.25|10|+agg. trunc newt(atomic)|

comments: if apply Newton's law with atomic patch of updating shared memory, we have 2x speedup with 1 thread only,
the advantage disappeared with multi-threading.

|argon_2916 (ms)|speedup|MPI_procs/OMP_threads|feature|
|--------------|--------|---------------------|-------|
|64287|1.0|1/1|+Newton +agg.|
|39563|1.62|2/1|+Newton +agg.|
|30437|2.11|3/1|+Newton +agg.|
|24761|2.60|4/1|+Newton +agg.|


* applying cell-list

profiling with original code (argon_2916)

```
Each sample counts as 0.01 seconds.
%   cumulative   self              self     total
time   seconds   seconds    calls  ms/call  ms/call  name
88.67     95.83    95.83     1001    95.74   107.17  force(_mdsys*)
10.59    107.28    11.45 6961373659     0.00     0.00  pbc(double, double)
0.82    108.17     0.89     1001     0.89     0.89  ekin(_mdsys*)
0.02    108.19     0.02     1000     0.02     0.02  velverlet_2(_mdsys*)
0.01    108.20     0.01     1000     0.01     0.01  velverlet_1(_mdsys*)
0.00    108.20     0.00     3006     0.00     0.00  azzero(double*, int)
0.00    108.20     0.00     2000     0.00     0.00  stamp()
0.00    108.20     0.00       12     0.00     0.00  get_a_line(_IO_FILE*, char*)
```

profiling with cell list (argon_2916)

```
Each sample counts as 0.01 seconds.
%   cumulative   self              self     total
time   seconds   seconds    calls  ms/call  ms/call  name
75.32     18.41    18.41     1001    18.39    24.18  cell_force(_mdsys*, _mdcell*)
23.69     24.20     5.79 3463801270     0.00     0.00  pbc(double, double)
0.59     24.35     0.15     1001     0.15     0.15  ekin(_mdsys*)
0.20     24.40     0.05     1001     0.05     0.09  sort(_mdsys*, _mdcell*)
0.16     24.44     0.04  2918916     0.00     0.00  index3d(_mdsys*, int, int, int)
0.08     24.46     0.02     1000     0.02     0.02  velverlet_2(_mdsys*)
0.04     24.47     0.01     1000     0.01     0.01  velverlet_1(_mdsys*)
0.00     24.47     0.00     3006     0.00     0.00  azzero(double*, int)
0.00     24.47     0.00     2000     0.00     0.00  stamp()
0.00     24.47     0.00      469     0.00     0.00  std::vector<int, std::allocator<int> >::_M_insert_aux(__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, int const&)
0.00     24.47     0.00       12     0.00     0.00  get_a_line(_IO_FILE*, char*)
0.00     24.47     0.00        8     0.00     0.00  std::vector<int, std::allocator<int> >::push_back(int const&)
0.00     24.47     0.00        1     0.00     0.00  pair(_mdsys*)

```

force function (+ sorting atoms into cells) accumulated timing:

|argon_108 (ms)|argon_2916 (ms)|feature|
|--------------|---------------|-------|
|2999|99288|original|
|1543|22633|serial cell-list|





