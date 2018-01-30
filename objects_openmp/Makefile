# -*- Makefile -*-
CC=g++
CFLAGS=-Wall -g -O3 -ffast-math -fexpensive-optimizations -fomit-frame-pointer -fopenmp -DCHUNKSIZE=1 #-DTIMING #-pg 
LDLIBS=-lm

# list of source files
SRCDIR=../src
OBJDIR=./
SOURCES= $(shell find $(SRCDIR) -type f -name *.c)
OBJ= $(patsubst $(SRCDIR)/%,./%,$(SOURCES:.c=.o))

default: ../ljmd-openmp.x

clean:
	rm -f $(OBJDIR)/*.o ../ljmd-openmp.x

# compilation pattern rule for objects
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

# linker rule
../ljmd-openmp.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)


