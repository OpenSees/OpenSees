# Makefile for fe objects

include ../../../../Makefile.def

PROGRAM         = test

ifndef CULA
CULA_SOLVER = 
endif

ifeq ($(PROGRAMMING_MODE), PARALLEL)

OBJS       = SparseGenColLinSOE.o \
	SparseGenColLinSolver.o \
	SparseGenRowLinSOE.o \
	SparseGenRowLinSolver.o \
	SuperLU.o \
	DistributedSuperLU.o \
	DistributedSparseGenColLinSOE.o \
	DistributedSparseGenRowLinSOE.o \
	PFEMSolver.o \
	PFEMLinSOE.o \
	PFEMCompressibleSolver.o \
	PFEMCompressibleLinSOE.o 
else

ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)

OBJS       = SparseGenColLinSOE.o \
	SparseGenColLinSolver.o \
	SparseGenRowLinSOE.o \
	SparseGenRowLinSolver.o \
	SuperLU.o \
	DistributedSuperLU.o \
	DistributedSparseGenColLinSOE.o \
	DistributedSparseGenRowLinSOE.o \
	PFEMSolver.o \
	PFEMLinSOE.o \
	PFEMCompressibleSolver.o \
	PFEMCompressibleLinSOE.o 
else

OBJS       = $(CULA_SOLVER) SparseGenColLinSOE.o \
	SparseGenColLinSolver.o \
	SparseGenRowLinSOE.o \
	SparseGenRowLinSolver.o \
	SuperLU.o \
	PFEMSolver.o \
	PFEMLinSOE.o \
	PFEMCompressibleSolver.o \
	PFEMCompressibleLinSOE.o 

endif
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
