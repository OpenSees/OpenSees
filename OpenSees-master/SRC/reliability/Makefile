include ../../Makefile.def

ifndef RELIABILITY_LIBRARY
  RELIABILITY_LIBRARY = $(HOME)/lib/libReliability.a
endif

RELIABILITY_OBJS = $(FE)/domain/pattern/DiscretizedRandomProcessSeries.o \
		$(FE)/domain/pattern/SimulatedRandomProcessSeries.o \
		$(FE)/reliability/analysis/analysis/FORMAnalysis.o \
		$(FE)/reliability/analysis/analysis/GFunVisualizationAnalysis.o \
		$(FE)/reliability/analysis/analysis/FOSMAnalysis.o \
		$(FE)/reliability/analysis/analysis/OutCrossingAnalysis.o \
		$(FE)/reliability/analysis/analysis/SamplingAnalysis.o \
		$(FE)/reliability/analysis/analysis/ReliabilityAnalysis.o \
		$(FE)/reliability/analysis/analysis/SORMAnalysis.o \
		$(FE)/reliability/analysis/analysis/SystemAnalysis.o \
		$(FE)/reliability/analysis/analysis/system/IPCM.o \
		$(FE)/reliability/analysis/analysis/system/PCM.o \
		$(FE)/reliability/analysis/analysis/system/SCIS.o \
		$(FE)/reliability/analysis/analysis/system/MVNcdf.o \
		$(FE)/reliability/analysis/analysis/BivariateDecomposition.o \
		$(FE)/reliability/analysis/analysis/DP_RSM_Sim.o \
		$(FE)/reliability/analysis/analysis/DP_RSM_Sim_TimeVariant.o \
		$(FE)/reliability/analysis/analysis/ExperimentalPointRule1D.o \
		$(FE)/reliability/analysis/analysis/GridPlane.o \
		$(FE)/reliability/analysis/analysis/ImportanceSamplingAnalysis.o \
		$(FE)/reliability/analysis/analysis/MonteCarloResponseAnalysis.o \
		$(FE)/reliability/analysis/analysis/MultiDimVisPrincPlane.o \
		$(FE)/reliability/analysis/analysis/OrthogonalPlaneSamplingAnalysis.o \
		$(FE)/reliability/analysis/analysis/PrincipalAxis.o \
		$(FE)/reliability/analysis/analysis/RespSurfaceSimulation.o \
		$(FE)/reliability/analysis/analysis/SurfaceDesign.o \
		$(FE)/reliability/analysis/analysis/UniformExperimentalPointRule1D.o \
		$(FE)/reliability/analysis/analysis/UnivariateDecomposition.o \
		$(FE)/reliability/analysis/convergenceCheck/OptimalityConditionReliabilityConvergenceCheck.o \
		$(FE)/reliability/analysis/convergenceCheck/ReliabilityConvergenceCheck.o \
		$(FE)/reliability/analysis/convergenceCheck/StandardReliabilityConvergenceCheck.o \
		$(FE)/reliability/analysis/curvature/CurvaturesBySearchAlgorithm.o \
		$(FE)/reliability/analysis/curvature/FindCurvatures.o \
		$(FE)/reliability/analysis/curvature/FirstPrincipalCurvature.o \
        $(FE)/reliability/analysis/curvature/CurvatureFitting.o \
		$(FE)/reliability/analysis/designPoint/SearchWithStepSizeAndStepDirection.o \
		$(FE)/reliability/analysis/designPoint/FindDesignPointAlgorithm.o \
		$(FE)/reliability/analysis/direction/HLRFSearchDirection.o \
		$(FE)/reliability/analysis/direction/SearchDirection.o \
		$(FE)/reliability/analysis/direction/PolakHeSearchDirectionAndMeritFunction.o \
		$(FE)/reliability/analysis/direction/SQPsearchDirectionMeritFunctionAndHessian.o \
		$(FE)/reliability/analysis/direction/GradientProjectionSearchDirection.o \
		$(FE)/reliability/analysis/gradient/GradientEvaluator.o \
		$(FE)/reliability/analysis/gradient/FiniteDifferenceGradient.o \
		$(FE)/reliability/analysis/gradient/ImplicitGradient.o \
		$(FE)/reliability/analysis/hessian/HessianEvaluator.o \
		$(FE)/reliability/analysis/hessian/FiniteDifferenceHessian.o \
		$(FE)/reliability/analysis/meritFunction/AdkZhangMeritFunctionCheck.o \
		$(FE)/reliability/analysis/meritFunction/CriteriaReductionMeritFunctionCheck.o \
		$(FE)/reliability/analysis/meritFunction/MeritFunctionCheck.o \
		$(FE)/reliability/analysis/misc/MatrixOperations.o \
		$(FE)/reliability/analysis/misc/CorrelatedStandardNormal.o \
		$(FE)/reliability/analysis/randomNumber/CStdLibRandGenerator.o \
		$(FE)/reliability/analysis/randomNumber/RandomNumberGenerator.o \
		$(FE)/reliability/analysis/rootFinding/RootFinding.o \
		$(FE)/reliability/analysis/rootFinding/SecantRootFinding.o \
		$(FE)/reliability/analysis/rootFinding/SafeGuardedZeroFindingAlgorithm.o \
		$(FE)/reliability/analysis/rootFinding/ZeroFindingAlgorithm.o \
		$(FE)/reliability/analysis/stepSize/ArmijoStepSizeRule.o \
		$(FE)/reliability/analysis/stepSize/FixedStepSizeRule.o \
		$(FE)/reliability/analysis/stepSize/StepSizeRule.o \
		$(FE)/reliability/analysis/transformation/NatafProbabilityTransformation.o \
		$(FE)/reliability/analysis/transformation/ProbabilityTransformation.o \
		$(FE)/reliability/analysis/transformation/AllIndependentTransformation.o \
		$(FE)/reliability/domain/components/CorrelationCoefficient.o \
		$(FE)/reliability/domain/components/CorrelationCoefficientIter.o \
		$(FE)/reliability/domain/components/RandomVariable.o \
		$(FE)/reliability/domain/components/RandomVariableIter.o \
		$(FE)/reliability/domain/components/ReliabilityDomainComponent.o \
		$(FE)/reliability/domain/components/ReliabilityDomain.o \
		$(FE)/reliability/domain/components/PrincipalPlane.o \
		$(FE)/reliability/domain/components/Cutset.o \
		$(FE)/reliability/domain/components/CutsetIter.o \
		$(FE)/reliability/domain/distributions/BetaRV.o \
		$(FE)/reliability/domain/distributions/GumbelRV.o \
		$(FE)/reliability/domain/distributions/ParetoRV.o \
		$(FE)/reliability/domain/distributions/Type1LargestValueRV.o \
		$(FE)/reliability/domain/distributions/UniformRV.o \
		$(FE)/reliability/domain/distributions/ChiSquareRV.o \
		$(FE)/reliability/domain/distributions/LaplaceRV.o \
		$(FE)/reliability/domain/distributions/RayleighRV.o \
		$(FE)/reliability/domain/distributions/Type1SmallestValueRV.o \
		$(FE)/reliability/domain/distributions/WeibullRV.o \
		$(FE)/reliability/domain/distributions/ExponentialRV.o \
		$(FE)/reliability/domain/distributions/LognormalRV.o \
		$(FE)/reliability/domain/distributions/ShiftedExponentialRV.o \
		$(FE)/reliability/domain/distributions/Type2LargestValueRV.o \
		$(FE)/reliability/domain/distributions/GammaRV.o \
		$(FE)/reliability/domain/distributions/NormalRV.o \
		$(FE)/reliability/domain/distributions/ShiftedRayleighRV.o \
		$(FE)/reliability/domain/distributions/Type3SmallestValueRV.o \
		$(FE)/reliability/domain/distributions/UserDefinedRV.o \
		$(FE)/reliability/domain/distributions/PythonRV.o \
		$(FE)/reliability/domain/functionEvaluator/FunctionEvaluator.o \
		$(FE)/reliability/domain/functionEvaluator/TclEvaluator.o \
		$(FE)/reliability/domain/performanceFunction/PerformanceFunction.o \
		$(FE)/reliability/domain/performanceFunction/PerformanceFunctionIter.o \
		$(FE)/reliability/domain/performanceFunction/LimitStateFunction.o \
		$(FE)/reliability/domain/performanceFunction/LimitStateFunctionIter.o \
		$(FE)/reliability/domain/modulatingFunction/ModulatingFunction.o \
		$(FE)/reliability/domain/modulatingFunction/ModulatingFunctionIter.o \
		$(FE)/reliability/domain/modulatingFunction/ConstantModulatingFunction.o \
		$(FE)/reliability/domain/modulatingFunction/GammaModulatingFunction.o \
		$(FE)/reliability/domain/modulatingFunction/TrapezoidalModulatingFunction.o \
		$(FE)/reliability/domain/modulatingFunction/KooModulatingFunction.o \
		$(FE)/reliability/domain/spectrum/JonswapSpectrum.o \
		$(FE)/reliability/domain/spectrum/NarrowBandSpectrum.o \
		$(FE)/reliability/domain/spectrum/PointsSpectrum.o \
		$(FE)/reliability/domain/spectrum/Spectrum.o \
		$(FE)/reliability/domain/spectrum/SpectrumIter.o \
		$(FE)/reliability/domain/filter/Filter.o \
		$(FE)/reliability/domain/filter/FilterIter.o \
		$(FE)/reliability/domain/filter/StandardLinearOscillatorAccelerationFilter.o \
		$(FE)/reliability/domain/filter/StandardLinearOscillatorDisplacementFilter.o \
		$(FE)/reliability/domain/filter/StandardLinearOscillatorVelocityFilter.o \
		$(FE)/reliability/domain/filter/KooFilter.o \
		$(FE)/reliability/domain/filter/DeltaFilter.o \
        $(FE)/reliability/domain/storage/ReliabilityStorage.o \
        $(FE)/reliability/domain/storage/FORMStorage.o \
		$(FE)/reliability/FEsensitivity/SensitivityAlgorithm.o \
		$(FE)/reliability/tcl/TclReliabilityBuilder.o 

# Compilation control
all:         
	     @$(CD) $(FE)/reliability/analysis; $(MAKE);
	     @$(CD) $(FE)/reliability/domain; $(MAKE);
	     @$(CD) $(FE)/reliability/FEsensitivity; $(MAKE);
	     @$(CD) $(FE)/reliability/tcl; $(MAKE);
	     $(AR) $(ARFLAGS) $(RELIABILITY_LIBRARY) $(RELIABILITY_OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) $(FE)/reliability/analysis; $(MAKE) wipe;
	@$(CD) $(FE)/reliability/domain; $(MAKE) wipe;
	@$(CD) $(FE)/reliability/FEsensitivity; $(MAKE) wipe;
	@$(CD) $(FE)/reliability/tcl; $(MAKE) wipe;
	@$(RM) $(RMFLAGS) $(RELIABILITY_LIBRARY)

# DO NOT DELETE THIS LINE -- make depend depends on it.
