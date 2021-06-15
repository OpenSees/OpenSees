# CalRel Manual Example 1
# by Kevin Mackie 2012/02/12
#
# >>>> FIRST-ORDER RELIABILITY ANALYSIS <<<<  
# limit-state function     1 
# -----------------------------------------------------------------------------
# iteration number ..............iter=         4
# value of limit-state function..g(x)=    8.6977E-04
# reliability index .............beta=       .877621
# probability ....................Pf1= 1.9007467E-01
# var          design point                     sensitivity vectors
#            x*            u*            alpha     gamma     delta      eta
# x1      2.407E+01    -7.040E-01        -.8022    -.8022     .7437    -.6680
# x2      2.407E+01     5.240E-01         .5970     .5970    -.5153    -.1049
# -----------------------------------------------------------------------------
#
#>>>> SENSITIVITY ANALYSIS AT COMPONENT LEVEL <<<<
# sensitivity with respect to distribution parameters
# limit-state function    1
# ----------------------------------------------------------------------
# d(beta)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1    2.576E-02 -2.314E-02  1.956E-02  6.202E-03                       
# x2   -2.576E-02 -5.243E-03  1.240E+01 -2.576E-02                       
# 
# d(Pf1)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1   -6.993E-03  6.281E-03 -5.310E-03 -1.683E-03                       
# x2    6.993E-03  1.423E-03 -3.366E+00  6.993E-03                       
# ----------------------------------------------------------------------
# 
# sensitivity with respect to limit-state function parameters
# limit-state function    1
# ----------------------------------------------------------------------
# par   d(beta)/d(parameter)     d(Pf1)/d(parameter)
#  1          1.973E+00               -5.355E-01
#  2         -6.201E-01                1.683E-01
#  3          2.576E-02               -6.993E-03
# ----------------------------------------------------------------------
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS -- POINT FITTING <<<< 
# limit-state function    1 
# -----------------------------------------------------------------------------
# coordinates and  ave. main curvatures of fitting points in rotated space
# axis  u'i    u'n    G(u)       a'+i         u'i    u'n    G(u)       a'-i
#    1 1.033   .839  9.273E-06 -7.248E-02   -1.010   .866 -3.684E-05 -2.216E-02
# 
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =       .852719             .853295
# probability                     Pf2 = 1.9690761E-01       1.9674785E-01
# -----------------------------------------------------------------------------
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS  -- CURVATURE FITTING <<<< 
# limit-state function    1 
# -----------------------------------------------------------------------------
# main curvatures in (n-1)x(n-1) space
#        1
#  1 -4.254E-02
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =       .855544             .855969
# probability                     Pf2 = 1.9612501E-01       1.9600759E-01
# -----------------------------------------------------------------------------
# 
# >>>> SECOND-ORDER DIRECTIONAL SIMULATION <<<<
# limit-state function     1
# ----------------------------------------------------------------------
# trials       Pf-mean       betag-mean     coef of var of pf
#   1000    1.9774327E-01    8.4970963E-01    2.0497537E-02
#   2000    1.9693374E-01    8.5262464E-01    1.4607498E-02
#   3000    1.9682853E-01    8.5300404E-01    1.1971354E-02
#   4000    1.9686546E-01    8.5287083E-01    1.0330301E-02
#   5000    1.9641100E-01    8.5451084E-01    9.2536119E-03
#
# >>>> MONTE CARLO SIMULATION <<<< 
# limit-state function     1 
# ----------------------------------------------------------------------
#      trials       Pf-mean       betag-mean     coef of var of Pf
#        2000    2.0150000E-01    8.3627537E-01    4.4523935E-02
#        4000    1.9825000E-01    8.4788861E-01    3.1800737E-02
#        6000    1.9550000E-01    8.5780525E-01    2.6190887E-02
#        8000    1.9512500E-01    8.5916406E-01    2.2708577E-02
#       10000    1.9280000E-01    8.6762435E-01    2.0462504E-02
#       12000    1.9333333E-01    8.6567817E-01    1.8647523E-02
#       14000    1.9385714E-01    8.6376992E-01    1.7235211E-02
#       16000    1.9543750E-01    8.5803161E-01    1.6040929E-02
#       18000    1.9550000E-01    8.5780525E-01    1.5120476E-02
#       20000    1.9615000E-01    8.5545373E-01    1.4314929E-02
#       22000    1.9631818E-01    8.5484606E-01    1.3641444E-02
#       24000    1.9683333E-01    8.5298671E-01    1.3039380E-02
#       26000    1.9642308E-01    8.5446722E-01    1.2544090E-02
#       28000    1.9735714E-01    8.5109912E-01    1.2052121E-02
#       30000    1.9710000E-01    8.5202537E-01    1.1652900E-02
#       32000    1.9662500E-01    8.5373830E-01    1.1299822E-02
#       34000    1.9776471E-01    8.4963253E-01    1.0923037E-02
#       36000    1.9841667E-01    8.4729028E-01    1.0593514E-02
#       38000    1.9773684E-01    8.4973274E-01    1.0333050E-02
#       40000    1.9812500E-01    8.4833755E-01    1.0059100E-02
#       42000    1.9861905E-01    8.4656415E-01    9.8014332E-03

# Create the reliability model builder ------------------------------------------------------------
reliability
model basic -ndm 2

# Create the random variables
randomVariable  1	uniform				-parameters 0.0 100.0
randomVariable  2	shiftedExponential	-parameters 0.05 0.0

# Specify correlation
# none

# Parameters
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
stepSizeRule                 Armijo           -maxNum 5    -base 0.5   -initial 1.0 2  -print 0
stepSizeRule                 Fixed           -stepSize 1.0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 30   -printDesignPointX CalRel_manual_1_output/1_designX.out

# Run the FORM analysis
runFORMAnalysis		CalRel_manual_1_output/1_FORM.out   -relSens 1

# Run the SORM analysis
findCurvatures		firstPrincipal
runSORMAnalysis		CalRel_manual_1_output/1_SORM_pc.out

hessianEvaluator			FiniteDifference -pert 1000
probabilityTransformation   Nataf           -print 0
findCurvatures				curvatureFitting
runSORMAnalysis				CalRel_manual_1_output/1_SORM_cf.out

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0
randomNumberGenerator        CStdLib
runImportanceSamplingAnalysis CalRel_manual_1_output/1_SAMPLEa.out -type responseStatistics -maxNum 50000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis CalRel_manual_1_output/1_SAMPLEb.out -type failureProbability -maxNum 50000 -targetCOV 0.01 -print 0
