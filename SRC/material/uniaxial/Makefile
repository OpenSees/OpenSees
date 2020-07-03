include ../../../Makefile.def

OBJS       = UniaxialMaterial.o \
	ElasticMaterial.o \
	Elastic2Material.o \
	ElasticPowerFunc.o \
	UVCuniaxial.o \
	Bilin.o \
	Bilin02.o \
	IMKBilin.o \
	IMKPeakOriented.o \
	IMKPinching.o \
	Steel2.o \
	OriginCentered.o \
	ConcretewBeta.o \
	ConcreteD.o \
	ConcreteSakaiKawashima.o \
	ConcreteCM.o \
	SteelMPF.o \
	c14-SK-M.o \
	ElasticMultiLinear.o \
        ElasticPPMaterial.o \
        ElasticBilin.o \
	ParallelMaterial.o \
	ResilienceLow.o \
	AxialSp.o \
	BWBN.o \
	AxialSpHD.o \
	AxialSp.o \
	Neoprene.o \
	Maxwell.o \
	ViscousDamper.o \
	BilinearOilDamper.o \
	Cast.o \
        HardeningMaterial.o \
        HardeningMaterial2.o \
	Concrete01.o \
	Concrete02Thermal.o \
	Steel02Thermal.o \
	Steel01Thermal.o \
	Concrete04.o \
	Concrete06.o \
	Concrete07.o \
	ResilienceMaterialHR.o \
	RambergOsgoodSteel.o \
	Steel01.o \
        HystereticMaterial.o \
	EPPGapMaterial.o \
        ViscousMaterial.o \
	SeriesMaterial.o \
	InitStrainMaterial.o \
	InitStressMaterial.o \
	SimpleFractureMaterial.o \
	MinMaxMaterial.o \
	TensionOnlyMaterial.o \
        PathIndependentMaterial.o \
	CableMaterial.o \
	ENTMaterial.o \
	FedeasMaterial.o \
	DrainMaterial.o \
	BoucWenOriginal.o \
    BoucWenMaterial.o \
	SecantConcrete.o \
	HardeningMaterial2.o \
	Pinching4Material.o \
	BarSlipMaterial.o \
	FatigueMaterial.o \
	SteelMP.o \
	SmoothPSConcrete.o \
	Steel02.o \
	Steel02Fatigue.o \
	GNGMaterial.o \
	Concrete02.o \
	Concrete02IS.o \
	Steel03.o \
	MultiLinear.o \
	HookGap.o \
	Bond_SP01.o \
	UniaxialJ2Plasticity.o \
	HyperbolicGapMaterial.o \
	ImpactMaterial.o \
	KinematicHardening.o \
	NewUniaxialMaterial.o \
	TclNewUnixialMaterial.o \
	ReinforcingSteel.o \
	PenaltyMaterial.o \
	ContinuumUniaxial.o \
	WrappedMaterial.o \
	SecantMaterial.o \
	Concrete01WithSITC.o \
	ConcretewBeta.o \
	FRPConfinedConcrete.o \
	FRPConfinedConcrete02.o \
	KikuchiAikenHDR.o \
	KikuchiAikenLRB.o \
	Concrete05.o \
	ECC01.o \
	PrestressedSteelMaterial.o \
	ElasticBDMaterial.o \
	ShearPanelMaterial.o \
	SelfCenteringMaterial.o \
	BackboneMaterial.o \
	WrapperUniaxialMaterial.o \
	SAWSMaterial.o \
	ConfinedConcrete01.o \
	ConfinedConcrete02.o \
	ModIMKPeakOriented.o \
	ModIMKPeakOriented02.o \
	ModIMKPinching.o \
	ModIMKPinching02.o \
	pyUCLA.o \
	Dodd_Restrepo.o DoddRestrepo.o \
	STEELDR.o DoddRestr.o \
	SteelBRB.o \
	CubicSpline.o \
	TriMatrix.o \
	CFSSSWP.o \
	CFSWSWP.o \
	Steel4.o \
	ConcreteECThermal.o \
	ElasticMaterialThermal.o \
	StainlessECThermal.o \
	SteelECThermal.o \
	DamperMaterial.o \
	MaterialState.o \
	OOHystereticMaterial.o \
	SPSW02.o \
	TDConcrete.o \
	TDConcreteEXP.o \
	TDConcreteMC10.o \
	TDConcreteMC10NL.o \
	PlateBearingConnectionThermal.o \
	DegradingPinchedBW.o \
        TclModelBuilderUniaxialMaterialCommand.o

all:         $(OBJS)
	@$(CD) $(FE)/material/uniaxial/fedeas; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/drain; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/snap; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/PY; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/backbone; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/stiffness; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/strength; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/unloading; $(MAKE);
	@$(CD) $(FE)/material/uniaxial/limitState; $(MAKE);

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
	@$(CD) $(FE)/material/uniaxial/backbone; $(MAKE) wipe;
	@$(CD) $(FE)/material/uniaxial/stiffness; $(MAKE) wipe;
	@$(CD) $(FE)/material/uniaxial/strength; $(MAKE) wipe;
	@$(CD) $(FE)/material/uniaxial/unloading; $(MAKE) wipe;
	@$(CD) $(FE)/material/uniaxial/limitState; $(MAKE) wipe;

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
