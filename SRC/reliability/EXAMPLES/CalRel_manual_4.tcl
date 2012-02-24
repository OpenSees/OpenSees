# CalRel Manual Example 4
# by Kevin Mackie 2012/02/12
#
# >>>> FIRST-ORDER RELIABILITY ANALYSIS <<<<  
# limit-state function     4
# ----------------------------------------------------------------------------- 
# iteration number ..............iter=         7 
# value of limit-state function..g(x)=    2.8200E-07 
# reliability index .............beta=      1.621329 
# probability ....................Pf1= 5.2473578E-02 
# var          design point                     sensitivity vectors 
#            x*            u*            alpha     gamma     delta      eta 
# x1      5.698E+02     7.586E-01        0.4677    0.4480   -0.4323   -0.2640 
# x2      2.189E+03     3.409E-01        0.2101    0.2801   -0.2822   -0.1052 
# x3      4.635E+00    -7.891E-01       -0.4870   -0.4699    0.4252   -0.3100 
# x4      5.128E+02     7.586E-01        0.4677    0.4480   -0.4323   -0.2640 
# x5      1.970E+03     3.409E-01        0.2101    0.2801   -0.2822   -0.1052 
# x6      4.172E+00    -7.891E-01       -0.4870   -0.4699    0.4252   -0.3100
# -----------------------------------------------------------------------------
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS -- POINT FITTING <<<< 
# limit-state function    4 
# ----------------------------------------------------------------------------- 
# coordinates and  ave. main curvatures of fitting points in rotated space 
# axis  u'i    u'n    G(u)       a'+i         u'i    u'n    G(u)       a'-i 
#    1 1.621  1.479  4.527E-07 -1.080E-01   -1.546  1.693 -6.963E-09  5.983E-02 
#    2 1.621  1.556  5.999E-08 -4.990E-02   -1.621  1.571  1.217E-08 -3.822E-02 
#    3 1.621  1.621  4.687E-12 -5.108E-04   -1.322  1.873  3.314E-08  2.885E-01 
#    4 1.621  1.524  3.055E-07 -7.421E-02   -1.621  1.564  5.066E-09 -4.334E-02 
#    5 1.621  1.563  4.947E-08 -4.456E-02   -1.621  1.578  9.745E-09 -3.325E-02 
# 
#                                       improved Breitung      Tvedt's EI 
# generalized reliability index betag =      1.584027            1.586346 
# probability                     Pf2 = 5.6593818E-02       5.6330469E-02 
# ----------------------------------------------------------------------------- 
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS  -- CURVATURE FITTING <<<< 
# limit-state function    4
# ----------------------------------------------------------------------------- 
# main curvatures in (n-1)x(n-1) space 
#        1          2          3          4          5 
#  1  1.645E-01 -5.180E-02  1.900E-01 -2.100E-01 -6.858E-02
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      1.570532            1.589283 
# probability                     Pf2 = 5.8145721E-02       5.5998269E-02 
# ----------------------------------------------------------------------------- 
# 
# >>>> MONTE CARLO SIMULATION <<<< 
# limit-state function     4
# ---------------------------------------------------------------------- 
#  trials       Pf-mean       betag-mean     coef of var of Pf 
#   1000    4.6000000E-02    1.6849407E+00    1.4408293E-01 
#   2000    6.3000000E-02    1.5300676E+00    8.6256757E-02 
#   3000    5.9333333E-02    1.5603946E+00    7.2707677E-02 
#   4000    5.9750000E-02    1.5568757E+00    6.2730249E-02 
#   5000    5.9200000E-02    1.5615247E+00    5.6382747E-02 
#   6000    6.0666667E-02    1.5492014E+00    5.0803703E-02 
#   7000    6.1428571E-02    1.5428914E+00    4.6722972E-02

# CREATE THE RELIABILITY MODEL BUILDER ------------------------------------------------------------
reliability
model basic -ndm 2

# CREATE RANDOM VARIABLES
randomVariable  1	lognormal		-mean 500.0		-stdv 100.0
randomVariable  2	lognormal		-mean 2000.0	-stdv 400.0
randomVariable	3	uniform			-mean 5.0 		-stdv 0.5
randomVariable  4	lognormal		-mean 450.0 	-stdv 90.0
randomVariable  5	lognormal		-mean 1800.0	-stdv 360.0
randomVariable	6	uniform			-mean 4.5 		-stdv 0.45

# SPECIFY CORRELATION
correlate 1 2	0.3
correlate 1 3	0.2
correlate 2 3	0.2
correlate 4 5	0.3
correlate 4 6	0.2
correlate 5 6	0.2

# PARAMETERS
parameter 1 randomVariable 1
parameter 2 randomVariable 2
parameter 3 randomVariable 3
parameter 4 randomVariable 4
parameter 5 randomVariable 5
parameter 6 randomVariable 6
parameter 7

updateParameter 7  1.7

# DEFINE LIMIT-STATE FUNCTION(s)
performanceFunction 1 "\$par(7) - \$par(2)/1000.0/\$par(3) - pow( \$par(1)/200.0/\$par(3) ,2) - \$par(5)/1000.0/\$par(6) - pow( \$par(4)/200.0/\$par(6) ,2)"

# CREATE NECESSARY RELIABILITY ANALYSIS TOOLS
probabilityTransformation    Nataf            -print 1
reliabilityConvergenceCheck  Standard         -e1 1.0e-3    -e2 1.0e-3  -print 1
functionEvaluator            Tcl
gradientEvaluator            FiniteDifference -pert 1000
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 10.0   -factor 0.0  
#stepSizeRule                 Armijo           -maxNum 5    -base 0.5   -initial 1.0 2  -print 0
stepSizeRule                 Fixed           -stepSize 1.0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 15
randomNumberGenerator        CStdLib
findCurvatures				 bySearchAlgorithm 5

# RUN THE FORM ANALYSIS FOLLOWED BY SYSTEM ANALYSIS WITH ALL g-FUNs IN SERIES
runFORMAnalysis		CalRel_manual_4_output/4_FORM.out   -relSens 1
# SORM does not appear to give correct answer (not same as CalRel)
runSORMAnalysis		CalRel_manual_4_output/4_SORM.out

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0

runImportanceSamplingAnalysis CalRel_manual_4_output/4_SAMPLEa.out -type responseStatistics -maxNum 50000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis CalRel_manual_4_output/4_SAMPLEb.out -type failureProbability -maxNum 50000 -targetCOV 0.01 -print 0
