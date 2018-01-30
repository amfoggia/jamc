# -*- Makefile -*-
CC=g++
CFLAGS=-Wall -fPIC -g -fomit-frame-pointer -O3 -ffast-math -fexpensive-optimizations -DCELL #-DTIMING #-pg
LDLIBS=-lm

# list of source files
SRCDIR=../src
OBJDIR=./
SOURCES= $(shell find $(SRCDIR) -type f -name *.c)
OBJ= $(patsubst $(SRCDIR)/%,./%,$(SOURCES:.c=.o))

default: ../ljmd-serial.x ../ljmd.so


# compilation pattern rule for objects
$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

# linker rule
../ljmd-serial.x: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# dynamic lib
../ljmd.so: $(filter-out $(OBJDIR)/ljmd.o,$(OBJ))
	$(CC) -o $@ -shared -fPIC $(LDLIBS) $^

clean:
	rm -f $(OBJDIR)/*.o ../*.so ../ljmd-serial.x
