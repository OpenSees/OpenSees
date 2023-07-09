include ../../../Makefile.def

OBJS       = UpdatedLagrangianBeam2D.o Elastic2DGNL.o \
	Inelastic2DYS02.o InelasticYS2DGNL.o Inelastic2DYS01.o \
	Inelastic2DYS03.o TclElement2dGNL.o TclElement2dYS.o \
	CyclicModel.o LinearCyclic.o BilinearCyclic.o QuadraticCyclic.o \
	TclCyclicModelCommands.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
