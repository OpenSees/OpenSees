include ../../../../Makefile.def

PROGRAM = test

OBJS       = DiagonalSOE.o \
	DiagonalSolver.o \
	DiagonalDirectSolver.o \
	DistributedDiagonalSOE.o \
	DistributedDiagonalSolver.o \
	MPIDiagonalSOE.o \
	MPIDiagonalSolver.o

ifeq ($(PROGRAMMING_MODE), SEQUENTIAL)
OBJS       = DiagonalSOE.o \
	DiagonalSolver.o \
	DiagonalDirectSolver.o \
	DistributedDiagonalSOE.o \
	DistributedDiagonalSolver.o
endif

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
