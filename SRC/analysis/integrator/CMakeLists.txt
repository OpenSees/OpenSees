#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
       ExplicitDifference.cpp
       ArcLength.cpp
       ArcLength1.cpp
       BackwardEuler.cpp
       CentralDifference.cpp
       IncrementalIntegrator.cpp
       Integrator.cpp
       LoadControl.cpp
       HarmonicSteadyState.cpp
       LoadPath.cpp
       MinUnbalDispNorm.cpp
       Newmark.cpp
       Newmark1.cpp
       StaticIntegrator.cpp
       StagedLoadControl.cpp
       TransientIntegrator.cpp
    PUBLIC
       ExplicitDifference.h
       ArcLength.h
       ArcLength1.h
       BackwardEuler.h
       CentralDifference.h
       IncrementalIntegrator.h
       Integrator.h
       LoadControl.h
       HarmonicSteadyState.h
       LoadPath.h
       MinUnbalDispNorm.h
       Newmark.h
       Newmark1.h
       StaticIntegrator.h
       StagedLoadControl.h
       TransientIntegrator.h
)
target_sources(OPS_Analysis
    PRIVATE
       AlphaOS.cpp
       AlphaOS_TP.cpp
       AlphaOSGeneralized.cpp
       AlphaOSGeneralized_TP.cpp            # Andreas Schellenberg
       CentralDifferenceAlternative.cpp     # fmk
       CentralDifferenceNoDamping.cpp
       Collocation.cpp                      # Andreas Schellenberg
       CollocationHSFixedNumIter.cpp
       CollocationHSIncrLimit.cpp
       CollocationHSIncrReduct.cpp
       DisplacementControl.cpp
       DistributedDisplacementControl.cpp
       EigenIntegrator.cpp                  # Jun Peng
       EQPath.cpp                           # M. Salehi
       GeneralizedAlpha.cpp
       HHT.cpp
       HHT_TP.cpp
       HHTExplicit.cpp
       HHTExplicit_TP.cpp
       HHTGeneralized.cpp
       HHTGeneralized_TP.cpp
       HHTGeneralizedExplicit.cpp
       HHTGeneralizedExplicit_TP.cpp
       HHTHSFixedNumIter.cpp
       HHTHSFixedNumIter_TP.cpp
       HHTHSIncrLimit.cpp
       HHTHSIncrLimit_TP.cpp
       HHTHSIncrReduct.cpp
       HHTHSIncrReduct_TP.cpp
       Houbolt.cpp                      # krm
       HSConstraint.cpp
       KRAlphaExplicit.cpp
       KRAlphaExplicit_TP.cpp
       NewmarkExplicit.cpp              # Andreas Schellenberg
       NewmarkHSFixedNumIter.cpp        # Andreas Schellenberg
       NewmarkHSIncrLimit.cpp           # Andreas Schellenberg
       NewmarkHSIncrReduct.cpp          # Andreas Schellenberg
       ParkLMS3.cpp
       TRBDF2.cpp
       TRBDF3.cpp
       WilsonTheta.cpp
       StagedNewmark.cpp
       GimmeMCK.cpp                     # MHS
    PUBLIC
       AlphaOS.h
       AlphaOS_TP.h
       AlphaOSGeneralized.h
       AlphaOSGeneralized_TP.h
       CentralDifferenceAlternative.h
       CentralDifferenceNoDamping.h
       Collocation.h
       CollocationHSFixedNumIter.h
       CollocationHSIncrLimit.h
       CollocationHSIncrReduct.h
       DisplacementControl.h
       DistributedDisplacementControl.h
       EigenIntegrator.h
       EQPath.h
       GeneralizedAlpha.h
       HHT.h
       HHT_TP.h
       HHTExplicit.h
       HHTExplicit_TP.h
       HHTGeneralized.h
       HHTGeneralized_TP.h
       HHTGeneralizedExplicit.h
       HHTGeneralizedExplicit_TP.h
       HHTHSFixedNumIter.h
       HHTHSFixedNumIter_TP.h
       HHTHSIncrLimit.h
       HHTHSIncrLimit_TP.h
       HHTHSIncrReduct.h
       HHTHSIncrReduct_TP.h
       Houbolt.h
       HSConstraint.h
       KRAlphaExplicit.h
       KRAlphaExplicit_TP.h
       NewmarkExplicit.h
       NewmarkHSFixedNumIter.h
       NewmarkHSIncrLimit.h
       NewmarkHSIncrReduct.h
       ParkLMS3.h
       TRBDF2.h
       TRBDF3.h
       WilsonTheta.h
       StagedNewmark.h
       GimmeMCK.h
)

target_sources(OPS_PFEM
    PRIVATE
       PFEMIntegrator.cpp
    PUBLIC
       PFEMIntegrator.h  # Uses ops_Dt
)

target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

