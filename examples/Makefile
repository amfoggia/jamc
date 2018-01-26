# -*- Makefile -*-

default: check

clean:
	rm -f *.xyz *.dat

check:  ../test_in.x ../test_velverlet.x ../test_ekin.x ../test_force.x ../ljmd-serial.x
	../test_in.x <argon_108.inp
	cmp test_in.dat ../reference/test_in.dat || exit 1
	../test_velverlet.x ../reference/test_force.dat < argon_108.inp
	cmp test_velverlet.dat ../reference/test_velverlet.dat || exit 1
	../test_ekin.x < argon_108.inp
	cmp test_ekin.dat ../reference/test_ekin.dat || exit 1
	../test_force.x < argon_108.inp
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
