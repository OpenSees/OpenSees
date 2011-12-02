include ../../../Makefile.def

OBJS       = NDMaterial.o \
             ElasticIsotropic3D.o \
             ElasticIsotropicMaterial.o \
             ElasticIsotropicPlaneStrain2D.o \
             ElasticIsotropicPlaneStress2D.o \
             ElasticIsotropicPlateFiber.o \
             ElasticIsotropicAxiSymm.o\
             TclModelBuilderNDMaterialCommand.o \
             J2PlaneStrain.o \
             PlateFiberMaterial.o \
             PlaneStressMaterial.o \
             J2PlaneStress.o \
             J2Plasticity.o \
             J2AxiSymm.o \
             J2ThreeDimensional.o \
             J2PlateFiber.o \
             PressureDependentElastic3D.o \
             ElasticIsotropicBeamFiber.o \
             BeamFiberMaterial.o \
             FeapMaterial.o \
             ElasticCrossAnisotropic.o 

all: $(OBJS) 
	@$(CD) $(FE)/material/nD/soil; $(MAKE);
	@$(CD) $(FE)/material/nD/feap; $(MAKE);
	@$(CD) $(FE)/material/nD/cyclicSoil; $(MAKE);
	@$(CD) $(FE)/material/nD/Template3Dep; $(MAKE);  
	@$(CD) $(FE)/material/nD/FiniteDeformation; $(MAKE);
	@$(CD) $(FE)/material/nD/NewTemplate3Dep; $(MAKE); 

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o
	@$(RM) $(RMFLAGS) .inslog*
	@$(RM) $(RMFLAGS) tca.map

spotless: clean
	@$(CD) $(FE)/material/nD/soil; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/feap; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/cyclicSoil; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/Template3Dep; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/FiniteDeformation; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/NewTemplate3Dep; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



