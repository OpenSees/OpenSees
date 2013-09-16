include ../../../Makefile.def

#	ElasticCrossAnisotropic.o \

ifdef DAMAGE2P_FLAG
	DAMAGE_2P = Damage2p.o Damage2p3D.o Damage2ppstress.o
else
	DAMAGE_2P =
endif

ifdef PLAINCONCRETE_FLAG
	PLAINCONCRETE = PlainConcrete.o forumat.o
else
	PLAINCONCRETE =
endif

OBJS       = NDMaterial.o \
	$(DAMAGE_2P) $(PLAINCONCRETE) BeamFiberMaterial.o \
	BeamFiberMaterial2d.o \
	ElasticIsotropicAxiSymm.o\
	ElasticIsotropicBeamFiber2d.o \
	ElasticIsotropicBeamFiber.o \
	ElasticIsotropicMaterial.o \
	ElasticIsotropicPlaneStrain2D.o \
	ElasticIsotropicPlaneStress2D.o \
	ElasticIsotropicPlateFiber.o \
	ElasticIsotropicThreeDimensional.o \
	ElasticOrthotropicMaterial.o \
	ElasticOrthotropicThreeDimensional.o \
	PressureDependentElastic3D.o \
	FeapMaterial.o \
	J2AxiSymm.o \
	J2PlaneStrain.o \
	J2PlaneStress.o \
	J2Plasticity.o \
	J2PlateFiber.o \
	J2ThreeDimensional.o \
	PlaneStrainMaterial.o \
	PlaneStressMaterial.o \
	PlaneStressUserMaterial.o \
	PlateFiberMaterial.o \
	PlateFromPlaneStressMaterial.o \
	PlateRebarMaterial.o \
	CycLiqCP.o \
	CycLiqCP3D.o \
	CycLiqCPPlaneStrain.o \
	CycLiqCPSP.o \
	CycLiqCPSP3D.o \
	CycLiqCPSPPlaneStrain.o \
	TclModelBuilderNDMaterialCommand.o \
	WrapperNDMaterial.o \
	InitStressNDMaterial.o \
	ConcreteMcftNonLinear5.o \
	ConcreteMcftNonLinear7.o \
	CapPlasticity.o \
	ConcreteS.o \
	SimplifiedJ2.o \
	PlaneStressSimplifiedJ2.o


all: $(OBJS) 
	@$(CD) $(FE)/material/nD/soil; $(MAKE);
	@$(CD) $(FE)/material/nD/feap; $(MAKE);
	@$(CD) $(FE)/material/nD/cyclicSoil; $(MAKE);
	@$(CD) $(FE)/material/nD/reinforcedConcretePlaneStress; $(MAKE); 
	@$(CD) $(FE)/material/nD/UWmaterials; $(MAKE);

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
	@$(CD) $(FE)/material/nD/reinforcedConcretePlaneStress; $(MAKE) wipe; 
	@$(CD) $(FE)/material/nD/UWmaterials; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



