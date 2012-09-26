# Makefile for fe objects

include ../../../../Makefile.def

PROGRAM = test

OBJS       = ProfileSPDLinSOE.o \
	ProfileSPDLinSolver.o \
	ProfileSPDLinDirectSolver.o \
	ProfileSPDLinSubstrSolver.o \
	ProfileSPDLinDirectBlockSolver.o \
	ProfileSPDLinDirectSkypackSolver.o \
	ProfileSPDLinSolverGather.o \
	ProfileSPDLinSOEGather.o \
	DistributedProfileSPDLinSOE.o \
	SProfileSPDLinSOE.o \
	SProfileSPDLinSolver.o 

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
