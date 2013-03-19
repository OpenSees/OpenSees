# CalRel Manual Example 4
# by Kevin Mackie 2012/02/12
#
# >>>> FIRST-ORDER RELIABILITY ANALYSIS <<<<  
# limit-state function     4
# ----------------------------------------------------------------------------- 
# iteration number ..............iter=         7
# value of limit-state function..g(x)=    1.4119E-07
# reliability index .............beta=      1.621329
# probability ....................Pf1= 5.2473539E-02
# var          design point                     sensitivity vectors
#            x*            u*            alpha     gamma     delta      eta
# x1      5.697E+02     7.581E-01         .4677     .4481    -.4323    -.2638
# x2      2.189E+03     3.406E-01         .2101     .2802    -.2823    -.1051
# x3      4.635E+00    -7.897E-01        -.4869    -.4698     .4252    -.3105
# x4      5.127E+02     7.581E-01         .4677     .4481    -.4323    -.2638
# x5      1.970E+03     3.406E-01         .2101     .2802    -.2823    -.1051
# x6      4.171E+00    -7.897E-01        -.4869    -.4698     .4252    -.3105
# -----------------------------------------------------------------------------
#
# >>>> SENSITIVITY ANALYSIS AT COMPONENT LEVEL <<<<
# sensitivity with respect to distribution parameters
# 
# limit-state function    4
# ----------------------------------------------------------------------
# d(beta)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1   -4.323E-03 -2.638E-03 -2.426E+00 -1.839E+00                       
# x2   -7.057E-04 -2.628E-04 -1.516E+00 -8.416E-01                       
# x3    8.505E-01 -6.210E-01  6.045E-01  2.460E-01                       
# x4   -4.804E-03 -2.931E-03 -2.426E+00 -1.839E+00                       
# x5   -7.841E-04 -2.920E-04 -1.516E+00 -8.416E-01                       
# x6    9.450E-01 -6.900E-01  6.717E-01  2.733E-01                       
# 
# d(Pf1)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1    4.634E-04  2.828E-04  2.600E-01  1.971E-01                       
# x2    7.563E-05  2.817E-05  1.625E-01  9.020E-02                       
# x3   -9.115E-02  6.655E-02 -6.479E-02 -2.636E-02                       
# x4    5.149E-04  3.142E-04  2.600E-01  1.971E-01                       
# x5    8.403E-05  3.130E-05  1.625E-01  9.020E-02                       
# x6   -1.013E-01  7.395E-02 -7.199E-02 -2.929E-02                       
# ----------------------------------------------------------------------
# 
# sensitivity with respect to limit-state function parameters
# limit-state function    4
# ----------------------------------------------------------------------
# par   d(beta)/d(parameter)     d(Pf1)/d(parameter)
#  1          3.211E+00               -3.441E-01
# ----------------------------------------------------------------------
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS -- POINT FITTING <<<< 
# limit-state function    4 
# ----------------------------------------------------------------------------- 
# coordinates and  ave. main curvatures of fitting points in rotated space
# axis  u'i    u'n    G(u)       a'+i         u'i    u'n    G(u)       a'-i
#    1 1.621  1.480  4.518E-07 -1.078E-01   -1.547  1.693 -6.953E-09  5.964E-02
#    2 1.621  1.556  5.987E-08 -4.983E-02   -1.621  1.571  1.225E-08 -3.831E-02
#    3 1.621  1.621  4.507E-12 -5.006E-04   -1.322  1.873  3.481E-08  2.883E-01
#    4 1.621  1.524  3.047E-07 -7.409E-02   -1.621  1.564  5.122E-09 -4.353E-02
#    5 1.621  1.563  4.936E-08 -4.449E-02   -1.621  1.578  9.811E-09 -3.334E-02
# 
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      1.583981            1.586300
# probability                     Pf2 = 5.6598969E-02       5.6335643E-02
# ----------------------------------------------------------------------------- 
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS  -- CURVATURE FITTING <<<< 
# limit-state function    4
# ----------------------------------------------------------------------------- 
# main curvatures in (n-1)x(n-1) space
# 
#        1          2          3          4          5
#  1  1.645E-01 -5.181E-02  1.899E-01 -2.100E-01 -6.859E-02
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      1.570465            1.589225
# probability                     Pf2 = 5.8153465E-02       5.6004824E-02
# ----------------------------------------------------------------------------- 
# 
# >>>> SECOND-ORDER DIRECTIONAL SIMULATION <<<<
# limit-state function     4
# ----------------------------------------------------------------------
# 
# trials       Pf-mean       betag-mean     coef of var of pf
#   2000    5.5750437E-02    1.5914832E+00    3.9090965E-02
#   4000    5.5802224E-02    1.5910228E+00    2.7697814E-02
#   6000    5.5463253E-02    1.5940426E+00    2.2643525E-02
#   8000    5.5484941E-02    1.5938489E+00    1.9685152E-02
#  10000    5.5439453E-02    1.5942551E+00    1.7683081E-02
#  12000    5.5227766E-02    1.5961490E+00    1.6186748E-02
#  14000    5.5131744E-02    1.5970100E+00    1.5042016E-02
#  16000    5.5576403E-02    1.5930329E+00    1.4004901E-02
#  18000    5.5789019E-02    1.5911402E+00    1.3170522E-02
#  20000    5.5811508E-02    1.5909403E+00    1.2472705E-02
#
# >>>> MONTE CARLO SIMULATION <<<< 
# limit-state function     4
# ---------------------------------------------------------------------- 
#      trials       Pf-mean       betag-mean     coef of var of Pf
#        2000    5.7500000E-02    1.5761120E+00    9.0552482E-02
#        4000    5.7250000E-02    1.5782857E+00    6.4170411E-02
#        6000    5.7166667E-02    1.5790119E+00    5.2433228E-02
#        8000    5.7875000E-02    1.5728654E+00    4.5111879E-02
#       10000    5.9100000E-02    1.5623736E+00    3.9902464E-02
#       12000    5.9833333E-02    1.5561743E+00    3.6187476E-02
#       14000    5.9642857E-02    1.5577787E+00    3.3559748E-02
#       16000    6.0437500E-02    1.5511114E+00    3.1171893E-02
#       18000    6.1777778E-02    1.5400197E+00    2.9047752E-02
#       20000    6.2250000E-02    1.5361566E+00    2.7445410E-02
#       22000    6.2000000E-02    1.5381989E+00    2.6224312E-02
#       24000    6.2250000E-02    1.5361566E+00    2.5054012E-02
#       26000    6.1961538E-02    1.5385137E+00    2.4130755E-02
#       28000    6.2464286E-02    1.5344111E+00    2.3152973E-02
#       30000    6.3033333E-02    1.5297983E+00    2.2259917E-02
#       32000    6.3343750E-02    1.5272957E+00    2.1496613E-02
#       34000    6.3088235E-02    1.5293550E+00    2.0899798E-02
#       36000    6.2500000E-02    1.5341206E+00    2.0412698E-02
#       38000    6.2789474E-02    1.5317711E+00    1.9819335E-02
#       40000    6.2400000E-02    1.5349342E+00    1.9381703E-02
#       42000    6.2047619E-02    1.5378094E+00    1.8971792E-02
#       44000    6.2181818E-02    1.5367129E+00    1.8514254E-02
#       46000    6.2000000E-02    1.5381989E+00    1.8135577E-02
#       48000    6.1708333E-02    1.5405898E+00    1.7798399E-02
#       50000    6.2240000E-02    1.5362381E+00    1.7359228E-02
#       52000    6.2269231E-02    1.5359997E+00    1.7017856E-02
#       54000    6.2518519E-02    1.5339700E+00    1.6664188E-02
#       56000    6.2642857E-02    1.5329600E+00    1.6346570E-02
#       58000    6.2724138E-02    1.5323006E+00    1.6051149E-02
#       60000    6.2750000E-02    1.5320910E+00    1.5777887E-02
#       62000    6.2822581E-02    1.5315029E+00    1.5511746E-02
#       64000    6.2906250E-02    1.5308256E+00    1.5256609E-02
#       66000    6.2818182E-02    1.5315385E+00    1.5034901E-02
#       68000    6.2823529E-02    1.5314952E+00    1.4811473E-02
#       70000    6.2714286E-02    1.5323805E+00    1.4611905E-02
#       72000    6.2500000E-02    1.5341206E+00    1.4433857E-02
#       74000    6.2486486E-02    1.5342305E+00    1.4239108E-02
#       76000    6.2434211E-02    1.5346558E+00    1.4056773E-02
#       78000    6.2461538E-02    1.5344334E+00    1.3872147E-02
#       80000    6.2375000E-02    1.5351378E+00    1.3707777E-02
#       82000    6.2256098E-02    1.5361068E+00    1.3553358E-02
#       84000    6.2357143E-02    1.5352832E+00    1.3379459E-02
#       86000    6.2267442E-02    1.5360143E+00    1.3233121E-02
#       88000    6.2443182E-02    1.5345827E+00    1.3062233E-02
#       90000    6.2633333E-02    1.5330373E+00    1.2895350E-02
#       92000    6.2608696E-02    1.5332374E+00    1.2757089E-02
#       94000    6.2691489E-02    1.5325654E+00    1.2611750E-02
#       96000    6.2562500E-02    1.5336126E+00    1.2493403E-02
#       98000    6.2744898E-02    1.5321323E+00    1.2346074E-02
#      100000    6.2830000E-02    1.5314428E+00    1.2213154E-02
#      102000    6.2862745E-02    1.5311777E+00    1.2089462E-02
#      104000    6.2913462E-02    1.5307673E+00    1.1967501E-02
#      106000    6.3075472E-02    1.5294580E+00    1.1837805E-02
#      108000    6.2962963E-02    1.5303670E+00    1.1738860E-02
#      110000    6.2954545E-02    1.5304350E+00    1.1632483E-02
#      112000    6.2964286E-02    1.5303563E+00    1.1527201E-02
#      114000    6.2991228E-02    1.5301385E+00    1.1423029E-02
#      116000    6.3112069E-02    1.5291627E+00    1.1312549E-02
#      118000    6.2974576E-02    1.5302731E+00    1.1229331E-02
#      120000    6.2925000E-02    1.5306740E+00    1.1140040E-02
#      122000    6.2893443E-02    1.5309293E+00    1.1051308E-02
#      124000    6.2822581E-02    1.5315029E+00    1.0968416E-02
#      126000    6.2968254E-02    1.5303242E+00    1.0867578E-02
#      128000    6.3132812E-02    1.5289953E+00    1.0767333E-02
#      130000    6.3323077E-02    1.5274621E+00    1.0667039E-02
#      132000    6.3401515E-02    1.5268311E+00    1.0578926E-02
#      134000    6.3380597E-02    1.5269993E+00    1.0501531E-02
#      136000    6.3272059E-02    1.5278729E+00    1.0433569E-02
#      138000    6.3217391E-02    1.5283133E+00    1.0362467E-02
#      140000    6.3228571E-02    1.5282232E+00    1.0287211E-02
#      142000    6.3274648E-02    1.5278520E+00    1.0210538E-02
#      144000    6.3298611E-02    1.5276591E+00    1.0137334E-02
#      146000    6.3445205E-02    1.5264799E+00    1.0055236E-02
#      148000    6.3418919E-02    1.5266912E+00    9.9892730E-03

# Create the reliability model builder ------------------------------------------------------------
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
findDesignPoint              StepSearch       -maxNumIter 30

# Run the FORM analysis
runFORMAnalysis		CalRel_manual_4_output/4_FORM.out   -relSens 1

# Run the SORM analysis
#findCurvatures				 bySearchAlgorithm 5
findCurvatures		firstPrincipal
runSORMAnalysis		CalRel_manual_4_output/4_SORM_pc.out

hessianEvaluator			FiniteDifference -pert 1000
probabilityTransformation   Nataf           -print 0
findCurvatures				curvatureFitting
runSORMAnalysis				CalRel_manual_4_output/4_SORM_cf.out

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0
randomNumberGenerator        CStdLib
runImportanceSamplingAnalysis CalRel_manual_4_output/4_SAMPLEa.out -type responseStatistics -maxNum 50000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis CalRel_manual_4_output/4_SAMPLEb.out -type failureProbability -maxNum 50000 -targetCOV 0.01 -print 0
