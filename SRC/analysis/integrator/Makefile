# Makefile for element objects

include ../../../Makefile.def

OBJS       = AlphaOS.o \
    AlphaOS_TP.o \
    EQPath.o \
    AlphaOSGeneralized.o \
    ExplicitDifference.o \
	AlphaOSGeneralized_TP.o \
	ArcLength.o \
	ArcLength1.o \
	BackwardEuler.o \
	CentralDifference.o \
	CentralDifferenceAlternative.o \
	CentralDifferenceNoDamping.o \
	Collocation.o \
	CollocationHSFixedNumIter.o \
	CollocationHSIncrLimit.o \
	CollocationHSIncrReduct.o \
	DisplacementControl.o \
	DistributedDisplacementControl.o \
	EigenIntegrator.o \
	GeneralizedAlpha.o \
	HHT.o \
	HHT_TP.o \
	HHTExplicit.o \
	HHTExplicit_TP.o \
	HHTGeneralized.o \
	HHTGeneralized_TP.o \
	HHTGeneralizedExplicit.o \
	HHTGeneralizedExplicit_TP.o \
	HHTHSFixedNumIter.o \
	HHTHSFixedNumIter_TP.o \
	HHTHSIncrLimit.o \
	HHTHSIncrLimit_TP.o \
	HHTHSIncrReduct.o \
	HHTHSIncrReduct_TP.o \
	Houbolt.o \
	HSConstraint.o \
    IncrementalIntegrator.o \
	Integrator.o \
	KRAlphaExplicit.o \
	KRAlphaExplicit_TP.o \
	LoadControl.o \
	StagedLoadControl.o \
	StagedNewmark.o \
	LoadPath.o \
	MinUnbalDispNorm.o \
	Newmark.o \
	Newmark1.o \
	NewmarkExplicit.o \
	NewmarkHSFixedNumIter.o \
	NewmarkHSIncrLimit.o \
	NewmarkHSIncrReduct.o \
	ParkLMS3.o \
	PFEMIntegrator.o \
	StaticIntegrator.o \
	TransientIntegrator.o \
	TRBDF2.o \
    TRBDF3.o \
	WilsonTheta.o \
	GimmeMCK.o 

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
