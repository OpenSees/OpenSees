# CREATE THE RELIABILITY MODEL BUILDER. NO FINITE ELEMENT MODEL BUILDER IS NECESSARY IN THIS EXAMPLE
reliability
model basic -ndm 2

# Units = kN, m

# CREATE RANDOM VARIABLES
for {set i 1} {$i <= 5} {incr i} {
    randomVariable $i lognormal -mean 134.9 -stdv 13.49
}
randomVariable 6 lognormal -mean 50.0 -stdv 15.0
randomVariable 7 lognormal -mean 40.0 -stdv 12.0

# SPECIFY CORRELATION
#correlate   1 2  0.0
#correlate   1 3  0.0
#correlate   2 3  0.0

# PARAMETERS
parameter 1 randomVariable 1
parameter 2 randomVariable 2
parameter 3 randomVariable 3
parameter 4 randomVariable 4
parameter 5 randomVariable 5
parameter 6 randomVariable 6
parameter 7 randomVariable 7

# Column height
set h 5.0

# DEFINE A LIMIT-STATE FUNCTION
performanceFunction 1 "\$par(1) + \$par(2) + \$par(4) + \$par(5) - $h*\$par(6)"
performanceFunction 2 "\$par(1) + 2*\$par(3) + 2*\$par(4) + \$par(5) - $h*\$par(6) - $h*\$par(7)"
performanceFunction 3 "\$par(2) + 2*\$par(3) + \$par(4) - $h*\$par(7)"

# CREATE NECESSARY RELIABILITY ANALYSIS TOOLS
probabilityTransformation    Nataf            -print 1
reliabilityConvergenceCheck  Standard         -e1 1.0e-3    -e2 1.0e-3  -print 1
#reliabilityConvergenceCheck  OptimalityCondition         -e1 1.0e-3    -e2 1.0e-3  -print 1
functionEvaluator            Tcl
gradientEvaluator            FiniteDifference -pert 1000
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 10.0   -factor 0.0  
#stepSizeRule                 Armijo           -maxNum 1    -base 0.5   -initial 1.0 2  -print 0
stepSizeRule                 Fixed           -stepSize 1.0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 40

# RUN THE FORM ANALYSIS FOLLOWED BY SYSTEM ANALYSIS WITH ALL g-FUNs IN SERIES
runFORMAnalysis		Madsen_eg_5_9_output/FORMresults.out
runFOSMAnalysis		Madsen_eg_5_9_output/FOSMresults.out

hessianEvaluator			FiniteDifference
findCurvatures				curvatureFitting
probabilityTransformation    Nataf           -print 0
runSORMAnalysis		Madsen_eg_5_9_output/SORMresults.out
runSystemAnalysis	Madsen_eg_5_9_output/SYSTEMSresults.out  IPCM 	allInSeries

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0
randomNumberGenerator        CStdLib
runImportanceSamplingAnalysis	Madsen_eg_5_9_output/SIM1results.out	-type failureProbability -variance 1.0 -maxNum 500000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis	Madsen_eg_5_9_output/SIM2results.out	-type responseStatistics -variance 1.0 -maxNum 500000 -targetCOV 0.01 -print 0
