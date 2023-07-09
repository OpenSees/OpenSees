include ../../../../Makefile.def

OBJS       = Coulomb.o \
	FrictionModel.o \
	FrictionResponse.o \
	TclModelBuilderFrictionModelCommand.o \
	VelDependent.o \
	VelDepMultiLinear.o \
    VelNormalFrcDep.o \
	VelPressureDep.o

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
