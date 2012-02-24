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
# x1      1.126E+02    -9.391E-01       -0.4856   -0.2947    0.2191   -0.2165 
# x2      1.126E+02    -6.109E-01       -0.3159   -0.2947    0.2191   -0.2165 
# x3      1.370E+02    -2.965E-01       -0.1534    0.0000    0.0000    0.0000 
# x4      1.169E+02    -7.807E-01       -0.4038   -0.3811    0.2848   -0.3235 
# x5      1.169E+02    -5.685E-01       -0.2940   -0.3811    0.2848   -0.3235 
# x6      9.180E+01     1.210E+00        0.6256    0.7321   -0.9414   -1.2617 
# x7      6.001E+01     0.000E+00        0.0000    0.0000    0.0000    0.0000 
# ----------------------------------------------------------------------------- 
# limit-state function     8 
# ----------------------------------------------------------------------------- 
# iteration number ..............iter=         9 
# value of limit-state function..g(x)=    9.9720E-05 
# reliability index .............beta=      1.514151 
# probability ....................Pf1= 6.4993803E-02 
# var          design point                     sensitivity vectors 
#            x*            u*            alpha     gamma     delta      eta 
# x1      1.250E+02    -4.587E-01       -0.3031   -0.1537    0.1164   -0.0568 
# x2      1.248E+02    -3.069E-01       -0.2028    0.0000    0.0000    0.0000 
# x3      1.291E+02    -7.446E-01       -0.4920   -0.4841    0.3541   -0.3360 
# x4      1.236E+02    -5.458E-01       -0.3606   -0.4210    0.3029   -0.2866 
# x5      1.328E+02    -2.235E-01       -0.1477   -0.2040    0.1524   -0.1024 
# x6      8.715E+01     9.140E-01        0.6033    0.5638   -0.4969   -0.5502 
# x7      6.549E+01     5.004E-01        0.3305    0.4531   -0.3726   -0.3849 
# ----------------------------------------------------------------------------- 
# limit-state function     9 
# ----------------------------------------------------------------------------- 
# iteration number ..............iter=         9 
# value of limit-state function..g(x)=    4.2855E-05 
# reliability index .............beta=      2.698649 
# probability ....................Pf1= 3.4810809E-03 
# var          design point                     sensitivity vectors 
#            x*            u*            alpha     gamma     delta      eta 
# x1      1.175E+02    -7.529E-01       -0.2789    0.0000    0.0000    0.0000 
# x2      9.807E+01    -1.270E+00       -0.4703   -0.2780    0.1880   -0.2861 
# x3      7.560E+01    -1.797E+00       -0.6658   -0.7337    0.5039   -1.0563 
# x4      9.903E+01    -7.154E-01       -0.2650   -0.3469    0.2409   -0.3876 
# x5      1.271E+02     0.000E+00        0.0000    0.0000    0.0000    0.0000 
# x6      7.790E+01     4.852E-01        0.1800    0.0000    0.0000    0.0000 
# x7      6.966E+01     1.061E+00        0.3937    0.5139   -0.6189   -0.8113
# -----------------------------------------------------------------------------
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
# probability                     Pf2 = 4.7778128E-02       4.7116535E-02 
# ----------------------------------------------------------------------------- 
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
# >>>> DIRECTIONAL SIMULATION <<<< 
# ---------------------------------------------------------------------- 
#  trials       Pf-mean       betag-mean     coef of var of pf 
#    500    5.0143914E-02    1.6434598E+00    9.2961329E-02 
#   1000    5.0685641E-02    1.6382417E+00    6.4194953E-02 
#
# >>>> MONTE CARLO SIMULATION <<<< 
# ---------------------------------------------------------------------- 
#  trials       Pf-mean       betag-mean     coef of var of Pf 
#   2000    5.7000000E-02    1.5804668E+00    9.0972896E-02 
#   4000    5.4500000E-02    1.6027041E+00    6.5865321E-02 

# CREATE THE RELIABILITY MODEL BUILDER. ------------------------------------------------------------
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
reliabilityConvergenceCheck  Standard         -e1 1.0e-4    -e2 1.0e-4  -print 1
functionEvaluator            Tcl
gradientEvaluator            FiniteDifference -pert 1000
searchDirection              iHLRF
meritFunctionCheck           AdkZhang         -multi 2.0    -add 10.0   -factor 0.5
stepSizeRule                 Fixed            -stepSize 1.0
#stepSizeRule                 Armijo           -maxNum 5    -base 0.5   -initial 1.0 2  -print 0
startPoint                   Mean
findDesignPoint              StepSearch       -maxNumIter 30
randomNumberGenerator        CStdLib
findCurvatures				 bySearchAlgorithm 6

# RUN THE FORM ANALYSIS FOLLOWED BY SYSTEM ANALYSIS WITH ALL g-FUNs IN SERIES
runFORMAnalysis		CalRel_manual_7_output/7_FORM.out   -relSens 1
# SORM does not appear to give correct answer (not same as CalRel)
runSORMAnalysis		CalRel_manual_7_output/7_SORM.out

# Now reset probability transformation so we don't get as much output for simulation
probabilityTransformation    Nataf           -print 0

runImportanceSamplingAnalysis CalRel_manual_7_output/7_SAMPLEa.out -type responseStatistics -maxNum 50000 -targetCOV 0.01 -print 0
runImportanceSamplingAnalysis CalRel_manual_7_output/7_SAMPLEb.out -type failureProbability -maxNum 50000 -targetCOV 0.01 -print 0

runSystemAnalysis	CalRel_manual_7_output/7_SERIES_pcm.out   PCM   allInSeries
runSystemAnalysis	CalRel_manual_7_output/7_SERIES_ipcm.out  IPCM  allInSeries
runSystemAnalysis	CalRel_manual_7_output/7_SERIES_mvn.out   MVN   allInSeries  -Nmax 1e10 -tol 1.0e-12
runSystemAnalysis	CalRel_manual_7_output/7_SERIES_scis.out  SCIS  allInSeries  -Nmax 1e10 -tol 1.0e-12

runSystemAnalysis	CalRel_manual_7_output/7_PARALLEL_pcm.out   PCM   allInParallel
runSystemAnalysis	CalRel_manual_7_output/7_PARALLEL_ipcm.out  IPCM  allInParallel
runSystemAnalysis	CalRel_manual_7_output/7_PARALLEL_mvn.out   MVN   allInParallel  -Nmax 1e10 -tol 1.0e-12
runSystemAnalysis	CalRel_manual_7_output/7_PARALLEL_scis.out  SCIS  allInParallel  -Nmax 1e10 -tol 1.0e-12

printReliability 1

# test the ability of cutset reliability routines
cutset 1	1
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET1_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET1_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET1_mvn.out   MVN   cutsets  -Nmax 1e10 -tol 1.0e-12
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET1_scis.out  SCIS  cutsets  -Nmax 1e10 -tol 1.0e-12
remove cutset 1

cutset 1	1 2
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET2_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET2_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET2_mvn.out   MVN   cutsets  -Nmax 1e10 -tol 1.0e-12
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET2_scis.out  SCIS  cutsets  -Nmax 1e10 -tol 1.0e-12
remove cutset 1

cutset 1	1 2 3
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET3_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET3_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET3_mvn.out   MVN   cutsets  -Nmax 1e10 -tol 1.0e-12
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET3_scis.out  SCIS  cutsets  -Nmax 1e10 -tol 1.0e-12
remove cutset 1

cutset 1	1
cutset 2	2
cutset 3	3
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET4_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET4_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET4_mvn.out   MVN   cutsets  -Nmax 1e10 -tol 1.0e-12
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET4_scis.out  SCIS  cutsets  -Nmax 1e10 -tol 1.0e-12
remove cutset 1
remove cutset 2
remove cutset 3

cutset 1	-2
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET5_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET5_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET5_mvn.out   MVN   cutsets  -Nmax 1e10 -tol 1.0e-12
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET5_scis.out  SCIS  cutsets  -Nmax 1e10 -tol 1.0e-12
remove cutset 1

cutset 1	-1 -2 -3
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET6_pcm.out   PCM   cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET6_ipcm.out  IPCM  cutsets
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET6_mvn.out   MVN   cutsets  -Nmax 1e10 -tol 1.0e-12
runSystemAnalysis	CalRel_manual_7_output/7_CUTSET6_scis.out  SCIS  cutsets  -Nmax 1e10 -tol 1.0e-12
remove cutset 1
