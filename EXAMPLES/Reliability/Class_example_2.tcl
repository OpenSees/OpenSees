# Create the reliability model builder ------------------------------------------------------------
reliability
model basic -ndm 2

# CREATE RANDOM VARIABLES
randomVariable 1 gamma 		-mean 500.0  -stdv 100.0
randomVariable 3 lognormal	-mean 600.0  -stdv 200.0
randomVariable 5 weibull	-mean 550.0  -stdv 150.0
randomVariable 7 lognormal	-mean 500.0  -stdv 290.5
randomVariable 9 normal		-mean 500.0  -stdv 100.0

# SPECIFY CORRELATION
#correlate   1 3  0.3
#correlate   1 5  0.2
#correlate   3 5  0.2

# assign random variables to parameters (do not need same tags)
parameter 1 randomVariable 1
parameter 3 randomVariable 3
parameter 5 randomVariable 5
parameter 7 randomVariable 7
parameter 9 randomVariable 9
#parameter 9 constant theta

# DEFINE LIMIT-STATE FUNCTION(S)
performanceFunction 5 "650 - \$par(9)"
performanceFunction 3 "750 - \$par(3)"
performanceFunction 1 "100 + \$par(1) - 2.0*\$par(3) + 3.0*\$par(5) - 4.0*\$par(7) + \$par(9)"

# Check analytic gradients (it will work with these or without)
gradPerformanceFunction 5 1 "0"
gradPerformanceFunction 5 3 "0"
gradPerformanceFunction 5 5 "0"
gradPerformanceFunction 5 7 "0"
gradPerformanceFunction 5 9 "-1"

# CREATE NECESSARY RELIABILITY ANALYSIS TOOLS
probabilityTransformation    Nataf           -print 1
reliabilityConvergenceCheck  Standard         -e1 1.0e-3    -e2 1.0e-3  -print 1
#reliabilityConvergenceCheck  OptimalityCondition         -e1 1.0e-3    -e2 1.0e-3  -print 1
functionEvaluator            Tcl
gradientEvaluator            FiniteDifference -pert 1000
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 10.0   -factor 0.0  
#stepSizeRule                 Armijo           -maxNum 1    -base 0.5   -initial 1.0 2  -print 0
stepSizeRule                 Fixed           -stepSize 1.0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 15
randomNumberGenerator        CStdLib

# Run some sequences of FORM and FOSM to check maintenance of parameter spaces
runFORMAnalysis		Class_example_2_output/FORMresults1.out
runFOSMAnalysis		Class_example_2_output/FOSMresults1.out
runFORMAnalysis		Class_example_2_output/FORMresults2.out
runFOSMAnalysis		Class_example_2_output/FOSMresults2.out

# create cutsets for generalized system reliability analysis
cutset 1	1 3
cutset 2	3 5

runSystemAnalysis	Class_example_2_output/SYSTEMresults_IPCM.out	IPCM 	allInParallel
runSystemAnalysis	Class_example_2_output/SYSTEMresults_PCM.out	PCM 	allInParallel
runSystemAnalysis	Class_example_2_output/SYSTEMresults_MVN.out	MVN 	allInParallel
runSystemAnalysis	Class_example_2_output/SYSTEMresults_SCIS.out	SCIS 	allInParallel
runSystemAnalysis	Class_example_2_output/SYSTEMresults_general_IPCM.out	IPCM	cutsets
runSystemAnalysis	Class_example_2_output/SYSTEMresults_general_PCM.out	PCM		cutsets
runSystemAnalysis	Class_example_2_output/SYSTEMresults_general_MVN.out	MVN		cutsets
runSystemAnalysis	Class_example_2_output/SYSTEMresults_general_SCIS.out	SCIS	cutsets

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0

runImportanceSamplingAnalysis	Class_example_2_output/SIM1results.out	-type failureProbability -variance 1.0 -maxNum 15000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis	Class_example_2_output/SIM2results.out	-type responseStatistics -variance 1.0 -maxNum 15000 -targetCOV 0.01 -print 0
