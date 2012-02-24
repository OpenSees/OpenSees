# CalRel Manual Example 2 (and 3 with restart and convergence)
# by Kevin Mackie 2012/02/12
#
# >>>> FIRST-ORDER RELIABILITY ANALYSIS <<<<  
# limit-state function     2 
# ----------------------------------------------------------------------------- 
# iteration number ..............iter=        14 
# value of limit-state function..g(x)=    3.7757E-07 
# reliability index .............beta=      1.772397 
# probability ....................Pf1= 3.8164339E-02 
# var          design point                     sensitivity vectors 
#            x*            u*            alpha     gamma     delta      eta 
#         6.319E+02     1.281E+00        0.7233    0.6962   -0.5905   -0.7863 
#         2.320E+03     4.814E-01        0.2719    0.3662   -0.3436   -0.2483 
#         4.526E+00    -1.126E+00       -0.6347   -0.6174    0.6303   -0.5981 
# -----------------------------------------------------------------------------
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS -- POINT FITTING <<<< 
# limit-state function    2 
# ----------------------------------------------------------------------------- 
# coordinates and  ave. main curvatures of fitting points in rotated space 
# axis  u'i    u'n    G(u)       a'+i         u'i    u'n    G(u)       a'-i 
#    1 1.682  1.859 -2.062E-08  6.100E-02   -1.308  2.138 -6.503E-12  4.283E-01 
#    2 1.772  1.685  5.294E-07 -5.562E-02   -1.772  1.714  7.873E-08 -3.727E-02
# 
#                                       improved Breitung      Tvedt's EI 
# generalized reliability index betag =      1.834630            1.835973 
# probability                     Pf2 = 3.3280233E-02       3.3180819E-02
# ----------------------------------------------------------------------------- 
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS  -- CURVATURE FITTING <<<< 
# limit-state function    2 
# ----------------------------------------------------------------------------- 
# main curvatures in (n-1)x(n-1) space 
#        1          2 
#  1  3.419E-01 -8.982E-02
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      1.849092            1.854139 
# probability                     Pf2 = 3.2222280E-02       3.1859603E-02
# ----------------------------------------------------------------------------- 
# 
# >>>> DIRECTIONAL SIMULATION <<<< 
# limit-state function     2 
# ---------------------------------------------------------------------- 
#  trials       Pf-mean       betag-mean     coef of var of pf 
#   1000    3.6186313E-02    1.7967668E+00    4.7026543E-02 
#   2000    3.5088238E-02    1.8107698E+00    3.3765775E-02 
#   3000    3.5571511E-02    1.8045633E+00    2.7483781E-02 
#
# >>>> MONTE CARLO SIMULATION <<<< 
# limit-state function     2
# ---------------------------------------------------------------------- 
#  trials       Pf-mean       betag-mean     coef of var of Pf 
#   2000    2.8500000E-02    1.9033107E+00    1.3058478E-01 
#   6000    3.4666667E-02    1.8162414E+00    6.8130751E-02
#  10000    3.5500000E-02    1.8054773E+00    5.2126511E-02 
#  14000    3.3857143E-02    1.8269034E+00    4.5148883E-02 
#  18000    3.3833333E-02    1.8272201E+00    3.9831755E-02 
#  22000    3.3500000E-02    1.8316739E+00    3.6214047E-02 
#  26000    3.4153846E-02    1.8229715E+00    3.2980395E-02 
#  30000    3.4033333E-02    1.8245651E+00    3.0759228E-02 
#  32000    3.3562500E-02    1.8308361E+00    2.9997947E-02

# CREATE THE RELIABILITY MODEL BUILDER ------------------------------------------------------------
reliability
model basic -ndm 2

# CREATE RANDOM VARIABLES
randomVariable  1	lognormal	-mean 500.0		-stdv 100.0
randomVariable  2	lognormal	-mean 2000.0 	-stdv 400.0
randomVariable	3	uniform		-mean 5.0 		-stdv 0.5

# SPECIFY CORRELATION
correlate 1 2	0.3
correlate 1 3	0.2
correlate 2 3	0.2

# PARAMETERS
parameter 1 randomVariable 1
parameter 2 randomVariable 2
parameter 3 randomVariable 3
parameter 4

updateParameter 4  1.0

# DEFINE LIMIT-STATE FUNCTION(s)
performanceFunction 1 "\$par(4) - \$par(2)/1000.0/\$par(3) - pow( \$par(1)/200.0/\$par(3),2)"

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
findDesignPoint              StepSearch       -maxNumIter 15   -printDesignPointX CalRel_manual_2_output/2_designX.out
randomNumberGenerator        CStdLib
findCurvatures				 bySearchAlgorithm 2
#findCurvatures				 firstPrincipal   -exe

# RUN THE FORM ANALYSIS FOLLOWED BY SYSTEM ANALYSIS WITH ALL g-FUNs IN SERIES
runFORMAnalysis		CalRel_manual_2_output/2_FORM.out   -relSens 1
# SORM gives a different answer regardless of firstPrincipal or using 2 curvatures
runSORMAnalysis		CalRel_manual_2_output/2_SORM.out

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0

runImportanceSamplingAnalysis CalRel_manual_2_output/2_SAMPLEa.out -type responseStatistics -maxNum 50000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis CalRel_manual_2_output/2_SAMPLEb.out -type failureProbability -maxNum 50000 -targetCOV 0.01 -print 0
