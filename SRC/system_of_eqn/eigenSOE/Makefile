include ../../../Makefile.def

OBJS       = EigenSOE.o \
	EigenSolver.o \
	ArpackSOE.o \
	ArpackSolver.o \
	SymBandEigenSOE.o \
	SymBandEigenSolver.o \
	FullGenEigenSOE.o \
	FullGenEigenSolver.o \
	SymmGeneralizedEigenSOE.o \
	SymmGeneralizedEigenSolver.o 

all:    $(OBJS)
	@$(CD) $(FE)/system_of_eqn/eigenSOE/sparsePython; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS)
	@$(CD) $(FE)/system_of_eqn/eigenSOE/sparsePython; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
