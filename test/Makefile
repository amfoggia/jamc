# -*- Makefile for unit test -*-
SHELL=/bin/sh
CC=g++
CFLAGS=-Wall -g -O3 -ffast-math -fomit-frame-pointer -fexpensive-optimizations -DCELL
LDLIBS=-lm
INC= -I$(INCDIR) -I.

OBJDIR=../objects_serial
SRCDIR=../src
INCDIR=../src
EXEDIR=..
TESTSRC= $(wildcard *.c)
EXE= $(addprefix $(EXEDIR)/,$(TESTSRC:.c=.x))

SOURCES= $(wildcard $(SRCDIR)/*.c)
OBJ= $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SOURCES:.c=.o))
FOBJ= $(filter-out $(OBJDIR)/ljmd.o,$(OBJ))
HEADS=$(wildcard $(INCDIR)*.h)

default: $(EXE)

# linker rule
$(EXEDIR)/%.x: %.o $(FOBJ)
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)

# compile true source files
$(OBJDIR)/%.o: $(SRCDIR)/%.c $(HEADS)
	$(CC) -c $(CFLAGS) $< -o $@

# compile testing source file
%.o: %.c $(HEADS)
	$(CC) $(INC) -c $(CFLAGS) $< -o $@

help:
	echo $(EXE)

clean:
	rm -f $(EXE) *.o $(OBJDIR)/*.o
