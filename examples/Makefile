# -*- Makefile -*-

default: check

clean:
	rm -f *.xyz *.dat

check: ../ljmd-mpi.x #../test_ekin.x ../test_cell.x ../test_velverlet_1.x ../test_velverlet_2.x ../ljmd-serial.x
	#../test_ekin.x
	#../test_force.x
	#../test_cell.x
	#../test_velverlet_1.x
	#cmp test_velverlet_1.dat ../reference/test_velverlet_1.dat || exit 1
	#../test_velverlet_2.x
	#cmp test_velverlet_2.dat ../reference/test_velverlet_2.dat || exit 1

	mpirun -np 2 ../ljmd-mpi.x < argon_108.inp
	#head -10 argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	#head -10 ../reference/argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	#cmp a.dat b.dat || exit 1
	#rm -f a.dat b.dat

	#../ljmd-serial.x < argon_108.inp
	#head -10 argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	#head -10 ../reference/argon_108.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	#cmp a.dat b.dat || exit 1
	#../ljmd-serial.x < argon_2916.inp
	#head -10 argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> a.dat
	#head -10 ../reference/argon_2916.dat | awk '{printf("%d %.6f %.6f %.6f\n",$$1,$$2,$$3,$$4);}'> b.dat
	#cmp a.dat b.dat || exit 1
	#rm -f a.dat b.dat
