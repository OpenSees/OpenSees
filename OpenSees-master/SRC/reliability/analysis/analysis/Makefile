include ../../../../Makefile.def

OBJS       = 	FORMAnalysis.o \
	FOSMAnalysis.o \
	BivariateDecomposition.o \
	GFunVisualizationAnalysis.o \
	OutCrossingAnalysis.o \
	SamplingAnalysis.o \
	ReliabilityAnalysis.o \
	SORMAnalysis.o \
	SystemAnalysis.o \
	DP_RSM_Sim.o \
	DP_RSM_Sim_TimeVariant.o \
	ExperimentalPointRule1D.o \
	GridPlane.o \
	MonteCarloResponseAnalysis.o \
	MultiDimVisPrincPlane.o \
	OrthogonalPlaneSamplingAnalysis.o \
	PrincipalAxis.o \
	RespSurfaceSimulation.o \
	SurfaceDesign.o \
	UnivariateDecomposition.o \
	UniformExperimentalPointRule1D.o \
	ImportanceSamplingAnalysis.o 


# Compilation control
all:    $(OBJS)
	@$(CD) $(FE)/reliability/analysis/analysis/system; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(CD) $(FE)/reliability/analysis/analysis/system; $(MAKE) wipe;

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.
