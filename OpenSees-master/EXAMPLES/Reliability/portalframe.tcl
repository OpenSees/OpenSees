model basic -ndm 2 -ndf 3

node 1 0 0
node 2 0 144
node 3 240 144
node 4 240 0

fix 1 1 1 1
fix 4 1 1 1

set E 30000
set Ag 25
set Ig 1500
set Ac 29
set Ic 2000

section Elastic 1 $E $Ag $Ig
section Elastic 2 $E $Ac $Ic

geomTransf Linear 1

element forceBeamColumn 1 1 2 1 Lobatto 2 3
element forceBeamColumn 2 2 3 1 Lobatto 1 3
element forceBeamColumn 3 3 4 1 Lobatto 2 3

set P 25.0
set w 1.0e-1
pattern Plain 1 Constant {
    load 2 $P 0 0
    eleLoad -ele 2 -type beamUniform -$w 
}

analysis Static


reliability

randomVariable 62 lognormal -mean $E -stdv [expr 0.1*$E]
randomVariable 32 normal -mean $P -stdv [expr 0.2*$P]
randomVariable 89 normal -mean 0 -stdv 1
randomVariable 41 normal -mean -$w -stdv [expr abs(0.2*$w)]

parameter 12 randomVariable 62 element 1 E
addToParameter 12 element 3 E
parameter 25 randomVariable 32 loadPattern 1 loadAtNode 2 1
parameter  3 randomVariable 89 node 1 coord 1
parameter 45 randomVariable 41 loadPattern 1 elementLoad 2 wy

parameter 23 node 2 disp 1

performanceFunction 76 "0.15-\$par(23)"

sensitivityIntegrator -static
sensitivityAlgorithm -computeAtEachStep

randomNumberGenerator        CStdLib
probabilityTransformation    Nataf            -print 3
reliabilityConvergenceCheck  Standard         -e1 1.0e-2    -e2 1.0e-2  -print 1
functionEvaluator                Tcl -file "analyze 1"
gradientEvaluator            Implicit
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 50    -factor 0.5  
stepSizeRule                 Armijo           -maxNum 10   -base 0.5   -initial 0.3 5
stepSizeRule Fixed -stepSize 1.0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 30;# -printDesignPointX designPointX.out


runFORMAnalysis  portalframe.out

foreach perf [getLSFTags] {
    puts "Performance Function $perf"
    puts "beta = [format %.7f $betaFORM($perf)]"
    foreach rv [getRVTags] {
	puts "\t x*($rv) = [format %7.4f $designPointXFORM($perf,$rv)], alpha($rv) = [format %.7f $alphaFORM($perf,$rv)], gamma($rv) = [format %.7f $gammaFORM($perf,$rv)]"
    }
}
