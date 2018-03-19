include ../../../../Makefile.def

OBJS       = EquiSolnAlgo.o Linear.o NewtonRaphson.o \
	ModifiedNewton.o NewtonLineSearch.o Broyden.o BFGS.o \
        KrylovNewton.o PeriodicNewton.o AcceleratedNewton.o \
        LineSearch.o InitialInterpolatedLineSearch.o NewtonHallM.o \
	SecantLineSearch.o RegulaFalsiLineSearch.o BisectionLineSearch.o

# Compilation control

all:         $(OBJS)
	@$(CD) $(FE)/analysis/algorithm/equiSolnAlgo/accelerator; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) $(FE)/analysis/algorithm/equiSolnAlgo/accelerator; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
