# CalRel Manual Example 7
# by Kevin Mackie 2012/02/12
#
# >>>> FIRST-ORDER RELIABILITY ANALYSIS <<<<  
# limit-state function     7
# -----------------------------------------------------------------------------
# iteration number ..............iter=        32
# value of limit-state function..g(x)=    6.2394E-04
# reliability index .............beta=      1.933953
# probability ....................Pf1= 2.6559427E-02
# var          design point                     sensitivity vectors
#            x*            u*            alpha     gamma     delta      eta
# x1      1.126E+02    -9.391E-01        -.4856    -.2947     .2191    -.2165
# x2      1.126E+02    -6.109E-01        -.3159    -.2947     .2191    -.2165
# x3      1.370E+02    -2.965E-01        -.1534     .0000     .0000     .0000
# x4      1.169E+02    -7.807E-01        -.4038    -.3811     .2848    -.3235
# x5      1.169E+02    -5.685E-01        -.2940    -.3811     .2848    -.3235
# x6      9.180E+01     1.210E+00         .6256     .7321    -.9414   -1.2617
# x7      6.001E+01     0.000E+00         .0000     .0000     .0000     .0000
# -----------------------------------------------------------------------------
# 
# limit-state function     8
# -----------------------------------------------------------------------------
# iteration number ..............iter=         9
# value of limit-state function..g(x)=    9.9720E-05
# reliability index .............beta=      1.514151
# probability ....................Pf1= 6.4993803E-02
# var          design point                     sensitivity vectors
#            x*            u*            alpha     gamma     delta      eta
# x1      1.250E+02    -4.587E-01        -.3031    -.1537     .1164    -.0568
# x2      1.248E+02    -3.069E-01        -.2028     .0000     .0000     .0000
# x3      1.291E+02    -7.446E-01        -.4920    -.4841     .3541    -.3360
# x4      1.236E+02    -5.458E-01        -.3606    -.4210     .3029    -.2866
# x5      1.328E+02    -2.235E-01        -.1477    -.2040     .1524    -.1024
# x6      8.715E+01     9.140E-01         .6033     .5638    -.4969    -.5502
# x7      6.549E+01     5.004E-01         .3305     .4531    -.3726    -.3849
# -----------------------------------------------------------------------------
# 
# limit-state function     9
# -----------------------------------------------------------------------------
# iteration number ..............iter=         9
# value of limit-state function..g(x)=    4.2855E-05
# reliability index .............beta=      2.698649
# probability ....................Pf1= 3.4810809E-03
# var          design point                     sensitivity vectors
#            x*            u*            alpha     gamma     delta      eta
# x1      1.175E+02    -7.529E-01        -.2789     .0000     .0000     .0000
# x2      9.807E+01    -1.270E+00        -.4703    -.2780     .1880    -.2861
# x3      7.560E+01    -1.797E+00        -.6658    -.7337     .5039   -1.0563
# x4      9.903E+01    -7.154E-01        -.2650    -.3469     .2409    -.3876
# x5      1.271E+02     0.000E+00         .0000     .0000     .0000     .0000
# x6      7.790E+01     4.852E-01         .1800     .0000     .0000     .0000
# x7      6.966E+01     1.061E+00         .3937     .5139    -.6189    -.8113
# -----------------------------------------------------------------------------
#
# >>>> SENSITIVITY ANALYSIS AT COMPONENT LEVEL <<<<
# sensitivity with respect to distribution parameters
# 
# limit-state function    7
# ----------------------------------------------------------------------
# d(beta)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1    9.526E-03 -9.413E-03  7.393E-03  3.741E-02                       
# x2    9.526E-03 -9.413E-03  7.393E-03  3.741E-02                       
# x3   -5.701E-19  4.282E-19 -4.388E-19 -3.439E-18                       
# x4    9.492E-03 -1.078E-02  6.793E-03  6.191E-02                       
# x5    9.492E-03 -1.078E-02  6.793E-03  6.191E-02                       
# x6   -4.707E-02 -6.309E-02 -5.324E-03 -4.175E-02                       
# x7    0.000E+00  0.000E+00  0.000E+00  0.000E+00                       
# 
# d(Pf1)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1   -5.857E-04  5.787E-04 -4.545E-04 -2.300E-03                       
# x2   -5.857E-04  5.787E-04 -4.545E-04 -2.300E-03                       
# x3    3.505E-20 -2.633E-20  2.698E-20  2.114E-19                       
# x4   -5.836E-04  6.630E-04 -4.177E-04 -3.806E-03                       
# x5   -5.836E-04  6.630E-04 -4.177E-04 -3.806E-03                       
# x6    2.894E-03  3.879E-03  3.273E-04  2.567E-03                       
# x7    0.000E+00  0.000E+00  0.000E+00  0.000E+00                       
# ----------------------------------------------------------------------
# limit-state function    8
# ----------------------------------------------------------------------
# d(beta)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1    5.059E-03 -2.472E-03  4.331E-03  1.245E-02                       
# x2   -1.067E-18  5.298E-19 -9.119E-19 -2.651E-18                       
# x3    1.012E-02 -9.599E-03  7.383E-03  7.245E-02                       
# x4    1.010E-02 -9.554E-03  7.581E-03  5.729E-02                       
# x5    5.081E-03 -3.413E-03  4.073E-03  2.260E-02                       
# x6   -2.484E-02 -2.751E-02 -4.480E-03 -2.036E-02                       
# x7   -2.484E-02 -2.566E-02 -5.014E-03 -1.983E-02                       
# 
# d(Pf1)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1   -6.414E-04  3.134E-04 -5.491E-04 -1.578E-03                       
# x2    1.352E-19 -6.717E-20  1.156E-19  3.361E-19                       
# x3   -1.283E-03  1.217E-03 -9.361E-04 -9.185E-03                       
# x4   -1.280E-03  1.211E-03 -9.612E-04 -7.264E-03                       
# x5   -6.442E-04  4.327E-04 -5.164E-04 -2.865E-03                       
# x6    3.150E-03  3.488E-03  5.681E-04  2.582E-03                       
# x7    3.150E-03  3.253E-03  6.357E-04  2.514E-03                       
# ----------------------------------------------------------------------
# limit-state function    9
# ----------------------------------------------------------------------
# d(beta)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1    8.124E-19 -6.463E-19  6.555E-19  2.731E-18                       
# x2    8.172E-03 -1.244E-02  5.643E-03  4.490E-02                       
# x3    1.440E-02 -3.018E-02  7.181E-03  1.974E-01                       
# x4    8.031E-03 -1.292E-02  5.045E-03  6.934E-02                       
# x5    0.000E+00  0.000E+00  0.000E+00  0.000E+00                       
# x6    0.000E+00  0.000E+00  0.000E+00  0.000E+00                       
# x7   -4.126E-02 -5.408E-02 -5.018E-03 -3.624E-02                       
# 
# d(Pf1)/d(parameter) :
# var     mean      std dev     par 1      par 2      par 3      par 4
# x1   -8.496E-21  6.759E-21 -6.856E-21 -2.856E-20                       
# x2   -8.547E-05  1.301E-04 -5.902E-05 -4.696E-04                       
# x3   -1.506E-04  3.156E-04 -7.511E-05 -2.065E-03                       
# x4   -8.400E-05  1.351E-04 -5.276E-05 -7.252E-04                       
# x5    0.000E+00  0.000E+00  0.000E+00  0.000E+00                       
# x6    0.000E+00  0.000E+00  0.000E+00  0.000E+00                       
# x7    4.316E-04  5.657E-04  5.248E-05  3.791E-04                       
# ----------------------------------------------------------------------
# 
# sensitivity with respect to limit-state function parameters
# 
# limit-state function    7
# ----------------------------------------------------------------------
# par   d(beta)/d(parameter)     d(Pf1)/d(parameter)
#  1         -8.643E-01                5.314E-02
# ----------------------------------------------------------------------
# limit-state function    8
# ----------------------------------------------------------------------
# par   d(beta)/d(parameter)     d(Pf1)/d(parameter)
#  1         -7.584E-01                9.615E-02
# ----------------------------------------------------------------------
# limit-state function    9
# ----------------------------------------------------------------------
# par   d(beta)/d(parameter)     d(Pf1)/d(parameter)
#  1         -5.749E-01                6.013E-03
# ----------------------------------------------------------------------
#
# >>>> SECOND-ORDER RELIABILITY ANALYSIS -- POINT FITTING <<<< 
# limit-state function    7
# -----------------------------------------------------------------------------
# coordinates and  ave. main curvatures of fitting points in rotated space
# axis  u'i    u'n    G(u)       a'+i         u'i    u'n    G(u)       a'-i
#    1 1.499  2.287  1.233E-07  3.145E-01   -1.652  2.180 -9.831E-05  1.800E-01
#    2 1.853  2.012 -2.052E-05  4.536E-02   -1.860  2.005 -8.491E-04  4.132E-02
#    3 1.917  1.951 -1.090E-05  9.291E-03   -1.917  1.951 -6.249E-05  9.313E-03
#    4 1.842  2.021 -1.708E-05  5.151E-02   -1.841  2.022 -9.630E-04  5.218E-02
#    5 1.930  1.938 -3.881E-07  2.232E-03   -1.894  1.973 -2.792E-04  2.189E-02
#    6 1.934  1.934 -9.132E-05  3.070E-06   -1.934  1.934 -9.132E-05  3.070E-06
# 
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      2.081729            2.086987
# probability                     Pf2 = 1.8683599E-02       1.8444667E-02
# -----------------------------------------------------------------------------
# 
# limit-state function    8
# -----------------------------------------------------------------------------
# coordinates and  ave. main curvatures of fitting points in rotated space
# axis  u'i    u'n    G(u)       a'+i         u'i    u'n    G(u)       a'-i
#    1 1.391  1.628 -8.719E-05  1.179E-01   -1.441  1.584 -2.179E-04  6.751E-02
#    2 1.490  1.538 -4.173E-05  2.107E-02   -1.495  1.533 -6.791E-05  1.673E-02
#    3 1.458  1.569 -9.134E-05  5.138E-02   -1.473  1.554 -1.898E-04  3.678E-02
#    4 1.514  1.503  1.263E-06 -9.519E-03   -1.514  1.510  9.554E-08 -3.992E-03
#    5 1.514  1.497  2.491E-06 -1.516E-02   -1.514  1.504  6.346E-07 -9.237E-03
#    6 1.232  1.751  1.748E-05  3.125E-01   -1.356  1.657 -7.147E-05  1.555E-01
# 
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      1.666790            1.673479
# probability                     Pf2 = 4.7778128E-02       4.7116539E-02
# -----------------------------------------------------------------------------
# 
# limit-state function    9
# -----------------------------------------------------------------------------
# coordinates and  ave. main curvatures of fitting points in rotated space
# axis  u'i    u'n    G(u)       a'+i         u'i    u'n    G(u)       a'-i
#    1 2.312  3.037 -1.484E-05  1.265E-01   -2.356  3.003 -6.252E-05  1.096E-01
#    2 2.330  3.023 -1.988E-05  1.195E-01   -2.361  2.998 -7.158E-05  1.075E-01
#    3 2.443  2.932 -4.861E-05  7.818E-02   -2.430  2.942 -2.347E-04  8.256E-02
#    4 2.699  2.686  2.631E-06 -3.406E-03   -2.600  2.793 -3.186E-04  2.801E-02
#    5 2.699  2.699 -7.154E-06  1.114E-07   -2.699  2.699 -7.154E-06  1.114E-07
#    6 2.499  2.885 -1.393E-04  5.957E-02   -2.517  2.869 -7.403E-04  5.374E-02
# 
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      2.861782            2.866446
# probability                     Pf2 = 2.1063360E-03       2.0755464E-03
# ----------------------------------------------------------------------------- 
# 
# >>>> SECOND-ORDER RELIABILITY ANALYSIS  -- CURVATURE FITTING <<<<
# limit-state function    7
# -----------------------------------------------------------------------------
# main curvatures in (n-1)x(n-1) space
# 
#        1          2          3          4          5          6
#  1  4.396E-01 -2.092E-02 -2.052E-10 -3.216E-02 -1.571E-02  0.000E+00
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      2.047672            2.052785
# probability                     Pf2 = 2.0296067E-02       2.0046705E-02
# -----------------------------------------------------------------------------
# 
# limit-state function    8
# -----------------------------------------------------------------------------
# main curvatures in (n-1)x(n-1) space
# 
#        1          2          3          4          5          6
#  1  1.953E-01  7.606E-11 -1.918E-02 -2.081E-02 -1.404E-02  2.547E-01
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      1.667588            1.675708
# probability                     Pf2 = 4.7698732E-02       4.6897717E-02
# -----------------------------------------------------------------------------
# 
# limit-state function    9
# -----------------------------------------------------------------------------
# main curvatures in (n-1)x(n-1) space
# 
#        1          2          3          4          5          6
#  1  4.259E-01  1.689E-10 -1.082E-02  1.210E-02  1.023E-09  0.000E+00
#                                       improved Breitung      Tvedt's EI
# generalized reliability index betag =      2.833430            2.837097
# probability                     Pf2 = 2.3025720E-03       2.2762907E-03
# -----------------------------------------------------------------------------

# Create the reliability model builder ------------------------------------------------------------
reliability
model basic -ndm 2

# CREATE RANDOM VARIABLES
randomVariable	1	weibull		-mean 134.0 	-stdv 23.0
randomVariable	2	weibull		-mean 134.0 	-stdv 23.0
randomVariable	3	weibull		-mean 160.0 	-stdv 35.0
randomVariable	4	weibull		-mean 150.0		-stdv 30.0
randomVariable	5	weibull		-mean 150.0		-stdv 30.0

randomVariable	6	uniform		-mean 65.0		-stdv 20.0
randomVariable	7	uniform		-mean 50.0		-stdv 15.0

# SPECIFY CORRELATION
correlate   1 2  0.4
correlate   1 3  0.2
correlate   1 4  0.2
correlate	1 5  0.2
correlate	2 3  0.4
correlate   2 4  0.2
correlate	2 5  0.2
correlate	3 4  0.4
correlate	3 5  0.2
correlate	4 5  0.4

correlate	6 7  0.4

# PARAMETERS
parameter 1 randomVariable 1
parameter 2 randomVariable 2
parameter 3 randomVariable 3
parameter 4 randomVariable 4
parameter 5 randomVariable 5
parameter 6 randomVariable 6
parameter 7 randomVariable 7
parameter 8

updateParameter 8  5.0

# DEFINE LIMIT-STATE FUNCTION(s)
performanceFunction 1 "\$par(1) + \$par(2) + \$par(4) + \$par(5) - \$par(6) * \$par(8)"
performanceFunction 2 "\$par(1) + 2.0 * \$par(3) + 2.0 * \$par(4) + \$par(5) - \$par(6) * \$par(8) - \$par(7) * \$par(8)"
performanceFunction 3 "\$par(2) + 2.0 * \$par(3) + \$par(4) - \$par(7) * \$par(8)"

# CREATE NECESSARY RELIABILITY ANALYSIS TOOLS
probabilityTransformation    Nataf            -print 1
reliabilityConvergenceCheck  Standard         -e1 1.0e-3    -e2 1.0e-3  -print 1
functionEvaluator            Tcl
gradientEvaluator            FiniteDifference -pert 1000
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 10.0   -factor 0.5
stepSizeRule                 Fixed            -stepSize 0.75
#stepSizeRule                 Armijo           -maxNum 10    -base 0.5   -initial 1.0 2  -print 1
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 40

# Run the FORM analysis
runFORMAnalysis		CalRel_manual_7_output/7_FORM.out   -relSens 1

# Run the SORM analysis
#findCurvatures				 bySearchAlgorithm 6
findCurvatures		firstPrincipal
runSORMAnalysis		CalRel_manual_7_output/7_SORM_pc.out

hessianEvaluator			FiniteDifference -pert 1000
probabilityTransformation   Nataf           -print 0
findCurvatures				curvatureFitting
runSORMAnalysis				CalRel_manual_7_output/7_SORM_cf.out

# Now reset probability transformation so we don't get as much output for simulation
randomNumberGenerator        CStdLib
probabilityTransformation    Nataf           -print 0
runImportanceSamplingAnalysis CalRel_manual_7_output/7_SAMPLEa.out -type responseStatistics -maxNum 50000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis CalRel_manual_7_output/7_SAMPLEb.out -type failureProbability -maxNum 50000 -targetCOV 0.01 -print 0

# Run some system analyses
runSystemAnalysis	CalRel_manual_7_output/7_SERIES_pcm.out   PCM   allInSeries
runSystemAnalysis	CalRel_manual_7_output/7_SERIES_ipcm.out  IPCM  allInSeries
runSystemAnalysis	CalRel_manual_7_output/7_SERIES_mvn.out   MVN   allInSeries  -Nmax 50000 -tol 1.0e-8
runSystemAnalysis	CalRel_manual_7_output/7_SERIES_scis.out  SCIS  allInSeries  -Nmax 50000 -tol 1.0e-8

runSystemAnalysis	CalRel_manual_7_output/7_PARALLEL_pcm.out   PCM   allInParallel
runSystemAnalysis	CalRel_manual_7_output/7_PARALLEL_ipcm.out  IPCM  allInParallel
runSystemAnalysis	CalRel_manual_7_output/7_PARALLEL_mvn.out   MVN   allInParallel  -Nmax 50000 -tol 1.0e-8
runSystemAnalysis	CalRel_manual_7_output/7_PARALLEL_scis.out  SCIS  allInParallel  -Nmax 50000 -tol 1.0e-8

printReliability 1

# test the ability of cutset reliability routines
cutset 1	1
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET1_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET1_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET1_mvn.out   MVN   cutsets  -Nmax 50000 -tol 1.0e-8
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET1_scis.out  SCIS  cutsets  -Nmax 50000 -tol 1.0e-8
remove cutset 1

cutset 1	1 2
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET2_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET2_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET2_mvn.out   MVN   cutsets  -Nmax 50000 -tol 1.0e-8
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET2_scis.out  SCIS  cutsets  -Nmax 50000 -tol 1.0e-8
remove cutset 1

cutset 1	1 2 3
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET3_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET3_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET3_mvn.out   MVN   cutsets  -Nmax 50000 -tol 1.0e-8
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET3_scis.out  SCIS  cutsets  -Nmax 50000 -tol 1.0e-8
remove cutset 1

cutset 1	1
cutset 2	2
cutset 3	3
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET4_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET4_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET4_mvn.out   MVN   cutsets  -Nmax 50000 -tol 1.0e-8
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET4_scis.out  SCIS  cutsets  -Nmax 50000 -tol 1.0e-8
remove cutset 1
remove cutset 2
remove cutset 3

cutset 1	-2
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET5_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET5_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET5_mvn.out   MVN   cutsets  -Nmax 50000 -tol 1.0e-8
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET5_scis.out  SCIS  cutsets  -Nmax 50000 -tol 1.0e-8
remove cutset 1

cutset 1	-1 -2 -3
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET6_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET6_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET6_mvn.out   MVN   cutsets  -Nmax 50000 -tol 1.0e-8
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET6_scis.out  SCIS  cutsets  -Nmax 50000 -tol 1.0e-8
remove cutset 1
