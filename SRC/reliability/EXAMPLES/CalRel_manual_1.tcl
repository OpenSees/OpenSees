# CalRel Manual Example 1
# by Kevin Mackie 2012/02/12
#
# >>>> FIRST-ORDER RELIABILITY ANALYSIS <<<<  
# limit-state function     1 
# ----------------------------------------------------------------------------- 
# iteration number ..............iter=         3 
# value of limit-state function..g(x)=   -5.9577E-07 
# reliability index .............beta=      0.877644 
# probability ....................Pf1= 1.9006858E-01 
# var          design point                     sensitivity vectors 
#            x*            u*            alpha     gamma     delta      eta 
# x1      2.407E+01    -7.041E-01       -0.8022   -0.8022    0.7438   -0.6681 
# x2      2.407E+01     5.240E-01        0.5970    0.5970   -0.5153   -0.1048 
# -----------------------------------------------------------------------------
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS -- POINT FITTING <<<< 
# limit-state function    1 
# ----------------------------------------------------------------------------- 
# coordinates and  ave. main curvatures of fitting points in rotated space 
# axis  u'i    u'n    G(u)       a'+i         u'i    u'n    G(u)       a'-i 
#    1 1.033  0.839  9.277E-06 -7.251E-02   -1.010  0.866  1.732E-10 -2.223E-02 
# 
#                                       improved Breitung      Tvedt's EI 
# generalized reliability index betag =      0.852718            0.853295 
# probability                     Pf2 = 1.9690791E-01       1.9674786E-01 
# ----------------------------------------------------------------------------- 
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS  -- CURVATURE FITTING <<<< 
# limit-state function    1 
# ----------------------------------------------------------------------------- 
# main curvatures in (n-1)x(n-1) space 
#        1 
#  1 -4.253E-02 
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      0.855571            0.855995 
# probability                     Pf2 = 1.9611764E-01       1.9600027E-01 
# ----------------------------------------------------------------------------- 
# 
# >>>> MONTE CARLO SIMULATION <<<< 
# limit-state function     1 
# ---------------------------------------------------------------------- 
#  trials       Pf-mean       betag-mean     coef of var of Pf 
#   1000    2.0800000E-01    8.1338038E-01    6.1737378E-02 
#   2000    2.0800000E-01    8.1338038E-01    4.3643998E-02

# CREATE THE RELIABILITY MODEL BUILDER ------------------------------------------------------------
reliability
model basic -ndm 2

# CREATE RANDOM VARIABLES
randomVariable  1	uniform				-parameters 0.0 100.0
randomVariable  2	shiftedExponential	-parameters 0.05 0.0

# SPECIFY CORRELATION
# none

# PARAMETERS
parameter 1 randomVariable 1
parameter 2 randomVariable 2
parameter 3
parameter 4
parameter 5

updateParameter 3  1.0
updateParameter 4  1.0
updateParameter 5  0.0

# DEFINE LIMIT-STATE FUNCTION(s)
performanceFunction 1 "pow( \$par(1),\$par(3) ) - \$par(4)*\$par(2) + \$par(5)"

# CREATE NECESSARY RELIABILITY ANALYSIS TOOLS
probabilityTransformation    Nataf            -print 1
reliabilityConvergenceCheck  Standard         -e1 1.0e-3    -e2 1.0e-3  -print 1
functionEvaluator            Tcl
gradientEvaluator            FiniteDifference -pert 1000
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 10.0   -factor 0.5
#stepSizeRule                 Armijo           -maxNum 5    -base 0.5   -initial 1.0 2  -print 0
stepSizeRule                 Fixed           -stepSize 1.0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 15   -printDesignPointX CalRel_manual_1_output/1_designX.out
randomNumberGenerator        CStdLib
findCurvatures				 firstPrincipal

# RUN THE FORM ANALYSIS FOLLOWED BY SYSTEM ANALYSIS WITH ALL g-FUNs IN SERIES
runFORMAnalysis		CalRel_manual_1_output/1_FORM.out   -relSens 1
runSORMAnalysis		CalRel_manual_1_output/1_SORM.out

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0

runImportanceSamplingAnalysis CalRel_manual_1_output/1_SAMPLEa.out -type responseStatistics -maxNum 50000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis CalRel_manual_1_output/1_SAMPLEb.out -type failureProbability -maxNum 50000 -targetCOV 0.01 -print 0
