# -*- Makefile -*-
SHELL=/bin/sh
############################################
# derived makefile variables
OBJ_SERIAL=$(SRC:src/%.f90=objects_serial/%.o)
############################################

default: serial

serial:
	$(MAKE) $(MFLAGS) -C objects_$@

mpi:
	$(MAKE) $(MFLAGS) -C objects_$@

unitest:
	$(MAKE) $(MFLAGS) -C test default

clean:
	$(MAKE) $(MFLAGS) -C objects_serial clean
	$(MAKE) $(MFLAGS) -C objects_mpi clean
	$(MAKE) $(MFLAGS) -C examples clean
	$(MAKE) $(MFLAGS) -C test clean

check: mpi unitest #serial unitest 
	$(MAKE) $(MFLAGS) -C examples check

