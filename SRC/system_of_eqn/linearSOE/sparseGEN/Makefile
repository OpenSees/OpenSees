# Makefile for fe objects

include ../../../../Makefile.def

PROGRAM         = test

OBJS       = SparseGenColLinSOE.o \
	SparseGenColLinSolver.o \
	SparseGenRowLinSOE.o \
	SparseGenRowLinSolver.o \
	SuperLU.o 

ifeq ($(PROGRAMMING_MODE), PARALLEL)

OBJS       = SparseGenColLinSOE.o \
	SparseGenColLinSolver.o \
	SparseGenRowLinSOE.o \
	SparseGenRowLinSolver.o \
	SuperLU.o \
	DistributedSuperLU.o \
	DistributedSparseGenColLinSOE.o \
	DistributedSparseGenRowLinSOE.o

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
