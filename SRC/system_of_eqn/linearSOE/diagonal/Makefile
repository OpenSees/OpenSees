# Makefile for fe objects

include ../../../../Makefile.def

PROGRAM = test

OBJS       = DiagonalSOE.o \
	DiagonalSolver.o \
	DiagonalDirectSolver.o \
	DistributedDiagonalSOE.o \
	DistributedDiagonalSolver.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
