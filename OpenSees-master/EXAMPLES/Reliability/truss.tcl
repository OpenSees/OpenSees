model basic -ndm 2 -ndf 3

set L 144.0

node 1 0.0 0.0
node 2  $L 0.0

fix 1 1 1 0
fix 2 0 1 0

set E 30000.0
set A 25.0
set I 1500.0

section Elastic 1 $E $A $I

geomTransf Linear 1

element forceBeamColumn 1 1 2 1 Lobatto 1 3

set P 25.0

pattern Plain 1 Constant {
    load 2 $P 0 0
}

analysis Static

reliability

randomVariable 62 normal -mean $E -stdv [expr 0.1*$E]
randomVariable 25 normal -mean $A -stdv [expr 0.1*$A]

#randomVariable 32 normal -mean $P -stdv [expr 0.2*$P]
#randomVariable 89 normal -mean 0 -stdv 1
#randomVariable 41 normal -mean -$w -stdv [expr abs(0.2*$w)]
#randomVariable 42 normal -mean -$P -stdv [expr abs(0.2*$P)]

parameter 12 randomVariable 62 element 1 E
parameter 13 randomVariable 25 element 1 A
#addToParameter 12 element 3 E
#parameter 25 randomVariable 32 loadPattern 1 loadAtNode 2 1
#parameter  3 randomVariable 89 node 1 coord 1
#parameter 45 randomVariable 41 loadPattern 1 elementLoad 2 wy
#parameter 46 randomVariable 42 loadPattern 2 elementLoad 2 P

parameter 23 node 2 disp 1

performanceFunction 76 "0.005-\$par(23)"
#performanceFunction 81 "0.005-\[nodeDisp 2 1\]"
#performanceFunction 23 "1500.0+\[sectionForce 1 1 2\]"

#gradPerformanceFunction 76 12 "-\[sensNodeDisp 2 1 12\]"
#gradPerformanceFunction 76 13 "-\[sensNodeDisp 2 1 13\]"
#gradPerformanceFunction 81 62 "-\[sensNodeDisp 2 1 12\]"

sensitivityIntegrator -static
sensitivityAlgorithm -computeAtEachStep

randomNumberGenerator        CStdLib
probabilityTransformation    Nataf            -print 1
reliabilityConvergenceCheck  Standard         -e1 1.0e-2    -e2 1.0e-2  -print 1
functionEvaluator            Tcl -file "analyze 1"
gradientEvaluator            FiniteDifference
gradientEvaluator			 Implicit
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 10.0    -factor 0.5  
#stepSizeRule                 Armijo           -maxNum 10   -base 0.5   -initial 0.3 5
stepSizeRule 				 Fixed -stepSize 1.0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 15  -printDesignPointX designPointX.out

runFORMAnalysis  truss_FORM.out   -relSens 1

print

#foreach perf [getLSFTags] {
#    puts "Performance Function $perf"
#    puts "beta = [format %.4f [betaFORM $perf]]"
#    foreach rv [getRVTags] {
		#puts "\t x*($rv) = [format %7.4f $par($rv,$perf)], alpha($rv) = [format %.4f [alphaFORM $perf $rv]], gamma($rv) = [format %.4f [gammaFORM $perf $rv]]"
		#puts "\t alpha($rv) = [format %.4f [alphaFORM $perf $rv]], gamma($rv) = [format %.4f [gammaFORM $perf $rv]]"
#    }
#}

#analyze 1

#puts "Computed u: [expr [nodeDisp 2 1]]"
#puts "Exact u: [expr $P*$L/($E*$A)]"


#puts "Computed du/dE: [expr [sensNodeDisp 2 1 12]]"
#puts "Exact du/dE: [expr -$P*$L/($A*$E*$E)]"

#puts "Computed du/dA: [expr [sensNodeDisp 2 1 13]]"
#puts "Exact du/dA: [expr -$P*$L/($A*$A*$E)]"

#print

#puts $par(23)

runFOSMAnalysis  truss_FOSM.out

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0

runImportanceSamplingAnalysis	truss_SIM1.out	-type failureProbability -variance 1.0 -maxNum 5000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis	truss_SIM2.out	-type responseStatistics -variance 1.0 -maxNum 5000 -targetCOV 0.01 -print 0
