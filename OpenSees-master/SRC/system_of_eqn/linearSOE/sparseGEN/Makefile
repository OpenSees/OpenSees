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
	PFEMCompressibleLinSOE.o \
	PFEMSolver_Mumps.o \
	PFEMSolver_LumpM.o \
	PFEMSolver_Umfpack.o \
	PFEMCompressibleSolver_Mumps.o \
	PFEMQuasiLinSOE.o \
	PFEMQuasiSolver.o \
	PFEMDiaLinSOE.o \
	PFEMDiaSolver.o \
	PFEMSolver_Laplace.o \
	PFEMSolver_LumpM.o
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
	PFEMCompressibleLinSOE.o \
	PFEMSolver_Mumps.o \
	PFEMSolver_Umfpack.o \
	PFEMCompressibleSolver_Mumps.o \
	PFEMQuasiLinSOE.o \
	PFEMQuasiSolver.o \
	PFEMDiaLinSOE.o \
	PFEMDiaSolver.o \
	PFEMSolver_Laplace.o \
	PFEMSolver_LumpM.o
else

OBJS       = $(CULA_SOLVER) \
	SparseGenColLinSOE.o \
	SparseGenColLinSolver.o \
	SparseGenRowLinSOE.o \
	SparseGenRowLinSolver.o \
	SuperLU.o \
	PFEMSolver.o \
	PFEMSolver_Umfpack.o \
	PFEMSolver_Mumps.o \
	PFEMLinSOE.o \
	PFEMCompressibleSolver.o \
	PFEMCompressibleLinSOE.o \
	PFEMQuasiLinSOE.o \
	PFEMQuasiSolver.o \
	PFEMDiaLinSOE.o \
	PFEMDiaSolver.o \
	PFEMSolver_Laplace.o \
	PFEMSolver_LumpM.o

endif
endif


all:         $(OBJS)

ifdef NVCC

# PFEMUnifiedSolver_Hybrid.o : PFEMUnifiedSolver_Hybrid.h PFEMUnifiedSolver_Hybrid.cpp
# 	$(NVCC) $(INCLUDES) -c -O3 -x cu PFEMUnifiedSolver_Hybrid.cpp
endif

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
