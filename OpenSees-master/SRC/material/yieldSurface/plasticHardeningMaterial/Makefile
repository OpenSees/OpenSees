include ../../../../Makefile.def

OBJS       = ExponReducing.o  MultiLinearKp.o  NullPlasticMaterial.o \
	PlasticHardeningMaterial.o \
	TclModelBuilderYSPlasticMaterialCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
