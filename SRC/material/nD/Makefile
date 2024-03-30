include ../../../Makefile.def

#	ElasticCrossAnisotropic.o \

ifdef DAMAGE2P_FLAG
	DAMAGE_2P = Damage2p.o Damage2p3D.o Damage2ppstress.o \
	Concrete.o CPlaneStress.o CPlaneStrain.o CPlaneStress2d.o \
	CThreeDimensional.o
else
	DAMAGE_2P =
endif

ifdef PLAINCONCRETE_FLAG
	PLAINCONCRETE = PlainConcrete.o forumat.o
else
	PLAINCONCRETE =
endif

#	PlasticDamageConcretePlaneStress.o \

OBJS       = NDMaterial.o \
	$(DAMAGE_2P) $(PLAINCONCRETE) BeamFiberMaterial.o \
	PlasticDamageConcretePlaneStress.o \
	BeamFiberMaterial2d.o \
	BeamFiberMaterial2dPS.o \
	PlasticDamageConcrete3d.o \
	LinearCap.o \
	eig3.o \
	Faria1998.o \
	Faria1998PlaneStrain.o \
	Faria1998PlaneStress.o \
	Faria1998PlaneStress2d.o \
	Faria1998ThreeDimensional.o \
	FSAM.o \
	AcousticMedium.o \
	ElasticIsotropicAxiSymm.o\
	ElasticIsotropicBeamFiber2d.o \
	ElasticIsotropicBeamFiber.o \
	ElasticIsotropicMaterial.o \
	ElasticIsotropicPlaneStrain2D.o \
	ElasticIsotropicPlaneStress2D.o \
	ElasticIsotropicPlateFiber.o \
	ElasticIsotropicThreeDimensional.o \
	IncrementalElasticIsotropicThreeDimensional.o \
	ElasticOrthotropicMaterial.o \
	ElasticOrthotropicThreeDimensional.o \
	PressureDependentElastic3D.o \
	FeapMaterial.o \
	J2AxiSymm.o \
	J2PlaneStrain.o \
	J2PlaneStress.o \
	J2Plasticity.o \
	J2PlateFiber.o \
	J2PlateFibre.o \
	J2BeamFiber2d.o \
	J2BeamFiber3d.o \
	J2ThreeDimensional.o \
	OrthotropicMaterial.o \
	Series3DMaterial.o \
	PlaneStrainMaterial.o \
	PlaneStressMaterial.o \
	PlaneStressUserMaterial.o \
	PlateFiberMaterial.o \
	PlateFromPlaneStressMaterial.o \
	PlaneStressLayeredMaterial.o \
	PlateRebarMaterial.o \
	PlaneStressRebarMaterial.o \
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
	LinearCap.o \
	AcousticMedium.o \
	PlaneStressSimplifiedJ2.o \
	DruckerPrager3DThermal.o \
	DruckerPragerThermal.o \
	ElasticIsotropic3DThermal.o \
	ElasticIsotropicMaterialThermal.o \
	J2PlasticityThermal.o \
	J2ThreeDimensionalThermal.o \
	PlateFiberMaterialThermal.o \
	PlateFromPlaneStressMaterialThermal.o \
	UVCmultiaxial.o	\
	UVCplanestress.o \
	PlateRebarMaterialThermal.o \
	ElasticPlaneStress.o \
	ElasticOrthotropicPlaneStress.o \
	VonPapaDamage.o \
	ASDConcrete3DMaterial.o


all: $(OBJS) 
	@$(CD) $(FE)/material/nD/soil; $(MAKE);
	@$(CD) $(FE)/material/nD/feap; $(MAKE);
	@$(CD) $(FE)/material/nD/cyclicSoil; $(MAKE);
	@$(CD) $(FE)/material/nD/reinforcedConcretePlaneStress; $(MAKE); 
	@$(CD) $(FE)/material/nD/UWmaterials; $(MAKE);
	@$(CD) $(FE)/material/nD/matCMM; $(MAKE);
	@$(CD) $(FE)/material/nD/stressDensityModel; $(MAKE);
	@$(CD) $(FE)/material/nD/UANDESmaterials; $(MAKE);
	@$(CD) $(FE)/material/nD/OrthotropicRotatingAngleConcreteT2DMaterial01; $(MAKE);
	@$(CD) $(FE)/material/nD/SmearedSteelDoubleLayerT2DMaterial01; $(MAKE);

test: $(OBJS) main.o
	$(LINKER) $(LINKFLAGS) main.o ../../../SRC/api/elementAPI_Dummy.o $(FE_LIBRARY) \
	$(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) \
	-o material_test

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o *.mod
	@$(RM) $(RMFLAGS) .inslog*
	@$(RM) $(RMFLAGS) tca.map

spotless: clean
	@$(CD) $(FE)/material/nD/soil; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/feap; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/cyclicSoil; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/reinforcedConcretePlaneStress; $(MAKE) wipe; 
	@$(CD) $(FE)/material/nD/UWmaterials; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/matCMM; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/stressDensityModel; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/UANDESmaterials; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/OrthotropicRotatingAngleConcreteT2DMaterial01; $(MAKE) wipe;
	@$(CD) $(FE)/material/nD/SmearedSteelDoubleLayerT2DMaterial01; $(MAKE) wipe;
	

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



