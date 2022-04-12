include ../../../../Makefile.def

OBJS       = HystereticBackbone.o \
	ArctangentBackbone.o \
	TrilinearBackbone.o \
	ReeseSoftClayBackbone.o \
	ReeseStiffClayBelowWS.o \
	ReeseStiffClayAboveWS.o \
	VuggyLimestone.o \
	CementedSoil.o \
	WeakRock.o \
	LiquefiedSand.o \
	ReeseSandBackbone.o \
	ManderBackbone.o \
	RaynorBackbone.o \
	MultilinearBackbone.o \
	CappedBackbone.o \
	LinearCappedBackbone.o \
	MaterialBackbone.o \
	TclModelBuilderBackboneCommand.o


all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
