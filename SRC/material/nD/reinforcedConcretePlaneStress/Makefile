include ../../../../Makefile.def

OBJS       = ConcreteL01.o \
	ConcreteZ01.o \
	TendonL01.o \
	SteelZ01.o \
	ReinforcedConcretePlaneStress.o \
	FAReinforcedConcretePlaneStress.o \
	FAFourSteelRCPlaneStress.o \
	RAFourSteelRCPlaneStress.o \
	PrestressedConcretePlaneStress.o \
	FAPrestressedConcretePlaneStress.o \
	FAFourSteelPCPlaneStress.o \
	RAFourSteelPCPlaneStress.o 

# Compilation control

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.









