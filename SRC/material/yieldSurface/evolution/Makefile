include ../../../../Makefile.def

OBJS       = BkStressLimSurface2D.o  CombinedIsoKin2D01.o  Isotropic2D01.o \
	Kinematic2D02.o  PeakOriented2D01.o  PlasticHardening2D.o \
	YS_Evolution.o BoundingSurface2D.o CombinedIsoKin2D02.o \
	Kinematic2D01.o  NullEvolution.o  PeakOriented2D02.o  YS_Evolution2D.o \
	TclModelBuilderYS_EvolutionCommand.o 

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
