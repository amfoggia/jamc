# -*- Makefile -*-

default: check

clean:
	rm -f *.xyz *.dat

check: ../test_force.x ../test_out.x ../test_in.x ../test_ekin.x ../test_cell.x ../test_velverlet_1.x ../test_velverlet_2.x ../ljmd-serial.x ../ljmd-openmp.x ../ljmd-mpi.x 
	../test_out.x
	cmp test_output.xyz ../reference/test_output.xyz || exit 1
	cmp test_output.dat ../reference/test_output.dat || exit 1
	../test_in.x <argon_108.inp
	cmp test_in.dat ../reference/test_in.dat || exit 1
	../test_ekin.x
	../test_cell.x
	../test_velverlet_1.x
	cmp test_velverlet_1.dat ../reference/test_velverlet_1.dat || exit 1
	../test_velverlet_2.x
	cmp test_velverlet_2.dat ../reference/test_velverlet_2.dat || exit 1
	../test_force.x
	cmp test_force.dat ../reference/test_force.dat || exit 1
	../ljmd-serial.x < argon_108.inp
	head -10 argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	head -10 ../reference/argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	cmp a.dat b.dat || exit 1
	../ljmd-serial.x < argon_2916.inp
	head -10 argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	head -10 ../reference/argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	cmp a.dat b.dat || exit 1
	rm -f a.dat b.dat
	export OMP_NUM_THREADS=4
	../ljmd-openmp.x < argon_108.inp
	head -10 argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	head -10 ../reference/argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	cmp a.dat b.dat || exit 1
	../ljmd-openmp.x < argon_2916.inp
	head -10 argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	head -10 ../reference/argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	cmp a.dat b.dat || exit 1
	rm -f a.dat b.dat
	mpirun -np 4 ../ljmd-mpi.x < argon_108.inp
	head -10 argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	head -10 ../reference/argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	cmp a.dat b.dat || exit 
	rm -f a.dat b.dat
	mpirun -np 4 ../ljmd-mpi.x < argon_2916.inp
	head -10 argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	head -10 ../reference/argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	cmp a.dat b.dat || exit 1
	rm -f a.dat b.dat
