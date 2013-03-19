# Class (CES6010) Example 1
# by Kevin Mackie 2012/04/12
#
# >>>> FIRST-ORDER RELIABILITY ANALYSIS <<<<
# limit-state function    44
# -----------------------------------------------------------------------------
# iteration number ..............iter=         6
# value of limit-state function..g(x)=   -3.3135E-04
# reliability index .............beta=      1.596433
# probability ....................Pf1= 5.5196131E-02
# var          design point                     sensitivity vectors
#            x*            u*            alpha     gamma     delta      eta
# x1      5.248E-01    -7.455E-01        -.4669    -.4669     .5864    -.4493
# x2      1.299E+02     1.412E+00         .8843     .8843    -.5838    -.8739
# -----------------------------------------------------------------------------
# 
# >>>> SENSITIVITY ANALYSIS AT COMPONENT LEVEL <<<<
# sensitivity with respect to distribution parameters
# 
# limit-state function   44
# ----------------------------------------------------------------------
# d(beta)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1    3.609E+00 -2.765E+00  1.896E+00 -1.414E+00                       
# x2   -2.919E-02 -4.370E-02 -2.919E-02  1.772E+01                       
# 
# d(Pf1)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1   -4.026E-01  3.084E-01 -2.116E-01  1.577E-01                       
# x2    3.256E-03  4.874E-03  3.256E-03 -1.977E+00                       
# ----------------------------------------------------------------------
# 
# sensitivity with respect to limit-state function parameters
# limit-state function   44
# ----------------------------------------------------------------------
# par   d(beta)/d(parameter)     d(Pf1)/d(parameter)
#  1          1.896E-03               -2.116E-04
#  2          5.895E-02               -6.576E-03
#  3         -9.230E+00                1.030E+00
#  4          1.123E-04               -1.253E-05
# ----------------------------------------------------------------------
# 
# >>>> SECOND-ORDER RELIABILITY ANALYSIS -- POINT FITTING <<<<
# limit-state function   44
# -----------------------------------------------------------------------------
# coordinates and  ave. main curvatures of fitting points in rotated space
# axis  u'i    u'n    G(u)       a'+i         u'i    u'n    G(u)       a'-i
#    1 1.596  1.568  7.754E-03 -2.253E-02   -1.596  1.557  3.731E-03 -3.057E-02
# 
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      1.582737            1.582863
# probability                     Pf2 = 5.6740722E-02       5.6726315E-02
# -----------------------------------------------------------------------------
# 
# >>>> SECOND-ORDER RELIABILITY ANALYSIS  -- CURVATURE FITTING <<<<
# limit-state function   44
# -----------------------------------------------------------------------------
# main curvatures in (n-1)x(n-1) space
#        1
#  1 -2.661E-02
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      1.582717            1.582843
# probability                     Pf2 = 5.6742976E-02       5.6728589E-02
# -----------------------------------------------------------------------------
# 
# >>>> SECOND-ORDER DIRECTIONAL SIMULATION <<<<
# limit-state function    44
# ----------------------------------------------------------------------
# trials       Pf-mean       betag-mean     coef of var of pf
#   1000    5.7020903E-02    1.5802842E+00    3.0588704E-02
#   2000    5.7919958E-02    1.5724773E+00    2.1015992E-02
#   3000    5.7686747E-02    1.5744931E+00    1.7266431E-02
#   4000    5.8026405E-02    1.5715593E+00    1.4837642E-02
#   5000    5.7512039E-02    1.5760075E+00    1.3335932E-02
#   6000    5.6842022E-02    1.5818490E+00    1.2241487E-02
#   7000    5.6761499E-02    1.5825547E+00    1.1320950E-02
#   8000    5.6603092E-02    1.5839452E+00    1.0631125E-02
#   9000    5.6639877E-02    1.5836221E+00    1.0016061E-02
#  10000    5.6809349E-02    1.5821353E+00    9.4794801E-03
# 
# >>>> MONTE CARLO SIMULATION <<<<
# limit-state function    44
# ----------------------------------------------------------------------
#      trials       Pf-mean       betag-mean     coef of var of Pf
#        2000    6.3500000E-02    1.5260397E+00    8.5893566E-02
#        4000    6.2250000E-02    1.5361566E+00    6.1375940E-02
#        6000    6.0500000E-02    1.5505899E+00    5.0878146E-02
#        8000    5.9500000E-02    1.5589847E+00    4.4453204E-02
#       10000    5.8800000E-02    1.5649271E+00    4.0010503E-02
#       12000    5.9833333E-02    1.5561743E+00    3.6187476E-02
#       14000    5.9428571E-02    1.5595886E+00    3.3624028E-02
#       16000    5.8875000E-02    1.5642878E+00    3.1609074E-02
#       18000    5.8888889E-02    1.5641695E+00    2.9797483E-02
#       20000    5.9150000E-02    1.5619491E+00    2.8201920E-02
#       22000    5.8181818E-02    1.5702213E+00    2.7126166E-02
#       24000    5.7916667E-02    1.5725057E+00    2.6034324E-02
#       26000    5.8230769E-02    1.5698005E+00    2.4941218E-02
#       28000    5.8107143E-02    1.5708638E+00    2.4061059E-02
#       30000    5.7733333E-02    1.5740899E+00    2.3324922E-02
#       32000    5.7656250E-02    1.5747572E+00    2.2600252E-02
#       34000    5.7705882E-02    1.5743275E+00    2.1915437E-02
#       36000    5.7277778E-02    1.5780438E+00    2.1382258E-02
#       38000    5.7421053E-02    1.5767976E+00    2.0784386E-02
#       40000    5.7475000E-02    1.5763290E+00    2.0248012E-02
#       42000    5.7261905E-02    1.5781820E+00    1.9798996E-02
#       44000    5.7022727E-02    1.5802682E+00    1.9386759E-02
#       46000    5.6608696E-02    1.5838960E+00    1.9034004E-02
#       48000    5.6666667E-02    1.5833868E+00    1.8623129E-02
#       50000    5.6400000E-02    1.5857325E+00    1.8292529E-02
#       52000    5.6384615E-02    1.5858681E+00    1.7939887E-02
#       54000    5.6537037E-02    1.5845260E+00    1.7579359E-02
#       56000    5.6321429E-02    1.5864254E+00    1.7297569E-02
#       58000    5.6396552E-02    1.5857629E+00    1.6984714E-02
#       60000    5.6200000E-02    1.5874976E+00    1.6730150E-02
#       62000    5.6435484E-02    1.5854199E+00    1.6421670E-02
#       64000    5.6640625E-02    1.5836155E+00    1.6131991E-02
#       66000    5.6560606E-02    1.5843187E+00    1.5897590E-02
#       68000    5.6514706E-02    1.5847225E+00    1.5668794E-02
#       70000    5.6157143E-02    1.5878764E+00    1.5495351E-02
#       72000    5.5916667E-02    1.5900065E+00    1.5313388E-02
#       74000    5.6148649E-02    1.5879516E+00    1.5071942E-02
#       76000    5.6078947E-02    1.5885683E+00    1.4882092E-02
#       78000    5.6192308E-02    1.5875656E+00    1.4674348E-02
#       80000    5.6275000E-02    1.5868351E+00    1.4478471E-02

# Create the reliability model builder ------------------------------------------------------------
reliability
model basic -ndm 2

# CREATE RANDOM VARIABLES
randomVariable  1	lognormal			-mean 0.65		-stdv [expr 0.25*0.65]
randomVariable  2	gumbel				-mean 100.0		-stdv [expr 0.2*100.0]

set gbar 32.17
set rbar 1000

# SPECIFY CORRELATION
# none

# PARAMETERS
parameter 1 randomVariable 1
parameter 2 randomVariable 2
parameter 3
parameter 4

updateParameter 3  2.0
updateParameter 4  0.0

# DEFINE LIMIT-STATE FUNCTION(s)
performanceFunction 1 "\$par(1)*$rbar*$gbar - pow($\par(2),\$par(3)) + \$par(4)"

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
hessianEvaluator			 FiniteDifference
findDesignPoint              StepSearch       -maxNumIter 15   -printDesignPointX CalRel_manual_1_output/1_designX.out

# RUN THE FORM ANALYSIS
runFORMAnalysis		Class_example_1_output/1_FORM.out   -relSens 1

# RUN THE SORM ANALYSIS
findCurvatures				firstPrincipal
runSORMAnalysis				Class_example_1_output/1_SORM_pc.out

probabilityTransformation    Nataf           -print 0
hessianEvaluator			FiniteDifference
findCurvatures				curvatureFitting
runSORMAnalysis				Class_example_1_output/1_SORM_cf.out

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0
randomNumberGenerator        CStdLib
runImportanceSamplingAnalysis Class_example_1_output/1_SAMPLEa.out -type responseStatistics -maxNum 50000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis Class_example_1_output/1_SAMPLEb.out -type failureProbability -maxNum 50000 -targetCOV 0.01 -print 0
