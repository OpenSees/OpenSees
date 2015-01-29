include ../../../../Makefile.def

OBJS       = RandomVibrationSimulatorResult.o \
	Analyzer.o \
	NewSearchWithStepSizeAndStepDirection.o \
	ranrotb.o \
	AnalyzerGFunEvaluator.o \
	ranrotw.o \
	AnalyzerGradGEvaluator.o \
	NewStandardLinearOscillatorAccelerationFilter.o \
	ReliabilityDirectIntegrationAnalysis.o \
	CrossingRateAnalyzer.o \
	ReliabilityStaticAnalysis.o \
	DynamciAnalyzer.o \
	NewWhitenoiseFilter.o \
	SelectLoadInitialStaticAnalysis.o \
	FirstPassageAnalyzer.o \
	NonStatFirstPassageAnalyzer.o \
	SmoothConcrete02.o \
	FOSeriesSimulation.o \
	NonStatRandomVibrationSimulation.o \
	StatFirstPassageAnalyzer.o \
	GeneralRandGenerator.o \
	OutCrossingResults.o \
	StaticAnalyzer.o \
	GFunEachStepEvaluator.o \
	PerformanceFunctionCoeff.o \
	StatRandomVibrationSimulation.o \
	InitialPointBuilder.o \
	PerformanceFunctionCoefficientIter.o \
	ThresholdIncInitialPointBuilder.o \
	InitialStaticAnalysis.o \
	rancombi.o \
	TimePoints.o \
	mersenne.o \
	RandomProcess.o \
	UniformGenerator.o \
	mother.o \
	RandomVibrationAnalysis.o \
	NewDiscretizedRandomProcessSeries.o \
	RandomVibrationSimulation.o \
	localMath.o

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
