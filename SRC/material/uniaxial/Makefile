include ../../../Makefile.def

 OBJS       = UniaxialMaterial.o ElasticMaterial.o \
        ElasticPPMaterial.o ParallelMaterial.o \
        HardeningMaterial.o Concrete01.o Steel01.o \
        HystereticMaterial.o EPPGapMaterial.o \
        ViscousMaterial.o SeriesMaterial.o MinMaxMaterial.o \
        PathIndependentMaterial.o CableMaterial.o \
	ENTMaterial.o FedeasMaterial.o DrainMaterial.o \
	BoucWenMaterial.o \
	Pinching4Material.o BarSlipMaterial.o \
	FatigueMaterial.o \
        TclModelBuilderUniaxialMaterialCommand.o

all:         $(OBJS)
	@$(CD) $(FE)/material/uniaxial/fedeas; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/drain; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/snap; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/PY; $(MAKE);

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/material/uniaxial/fedeas; $(MAKE) wipe;
	@$(CD) $(FE)/material/uniaxial/drain; $(MAKE) wipe;
	@$(CD) $(FE)/material/uniaxial/snap; $(MAKE) wipe;
	@$(CD) $(FE)/material/uniaxial/PY; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
