include ../../../Makefile.def

OBJS       = NDMaterial.o ElasticIsotropicMaterial.o \
	ElasticIsotropicPlaneStrain2D.o ElasticIsotropicPlaneStress2D.o \
	TclModelBuilderNDMaterialCommand.o

all:	$(OBJS)

# Miscellaneous

tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
