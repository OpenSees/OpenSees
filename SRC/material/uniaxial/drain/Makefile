include ../../../../Makefile.def

OBJS       = DrainHardeningMaterial.o anal00.o common00.o \
	DrainBilinearMaterial.o \
	DrainClough1Material.o \
	DrainClough2Material.o \
	DrainPinch1Material.o \
	TclDrainMaterialCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
