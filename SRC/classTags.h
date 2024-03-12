/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written: fmk
// Revision: A
//
// Purpose: This file contains the declarations of all class tags used.

#ifndef classTags_h
#define classTags_h

#define svnTest  0

#define intType    1
#define doubleType 2
#define idType     3
#define vectorType 4
#define matrixType 5



#define ACTOR_TAGS_SUBDOMAIN 1

#define DMG_TAG_HystereticEnergy 1
#define DMG_TAG_ParkAng          2
#define DMG_TAG_Kratzig          3
#define DMG_TAG_Mehanny          4
#define DMG_TAG_NormalizedPeak   5


#define EigenSOE_TAGS_BandArpackSOE 	1
#define EigenSOE_TAGS_SymArpackSOE 	2
#define EigenSOE_TAGS_SymBandEigenSOE   3
#define EigenSOE_TAGS_FullGenEigenSOE   4
#define EigenSOE_TAGS_ArpackSOE 	5
#define EigenSOE_TAGS_GeneralArpackSOE 	6
#define EigenSOLVER_TAGS_BandArpackSolver 	1
#define EigenSOLVER_TAGS_SymArpackSolver 	2
#define EigenSOLVER_TAGS_SymBandEigenSolver     3
#define EigenSOLVER_TAGS_FullGenEigenSolver  4
#define EigenSOLVER_TAGS_ArpackSolver  5
#define EigenSOLVER_TAGS_GeneralArpackSolver  6

#define EigenALGORITHM_TAGS_Frequency 1
#define EigenALGORITHM_TAGS_Standard  2

#define EigenINTEGRATOR_TAGS_Eigen 1

#define CONVERGENCE_TEST_CTestNormUnbalance                 1
#define CONVERGENCE_TEST_CTestNormDispIncr                  2
#define CONVERGENCE_TEST_CTestEnergyIncr                    3
#define CONVERGENCE_TEST_CTestRelativeNormUnbalance         4
#define CONVERGENCE_TEST_CTestRelativeNormDispIncr          5
#define CONVERGENCE_TEST_CTestRelativeEnergyIncr            6
#define CONVERGENCE_TEST_CTestRelativeTotalNormDispIncr     7
#define CONVERGENCE_TEST_CTestFixedNumIter                  8
#define CONVERGENCE_TEST_NormDispAndUnbalance               9
#define CONVERGENCE_TEST_NormDispOrUnbalance               10
#define CONVERGENCE_TEST_CTestPFEM                         11


#define GRND_TAG_ElCentroGroundMotion                 1
#define GROUND_MOTION_TAG_GroundMotionRecord          2
#define GROUND_MOTION_TAG_InterpolatedGroundMotion    3
#define GROUND_MOTION_TAG_GroundMotion                4

#define REGION_TAG_MeshRegion      1

#define TIMESERIES_INTEGRATOR_TAG_Trapezoidal 1
#define TIMESERIES_INTEGRATOR_TAG_Simpson     2

#define SECT_TAG_Section         1

#define TSERIES_TAG_LinearSeries         1
#define TSERIES_TAG_RectangularSeries          2
#define TSERIES_TAG_PathTimeSeries       3
#define TSERIES_TAG_PathSeries       4
#define TSERIES_TAG_ConstantSeries       5
#define TSERIES_TAG_TrigSeries       6
#define TSERIES_TAG_DiscretizedRandomProcessSeries 7
#define TSERIES_TAG_SimulatedRandomProcessSeries 8
#define TSERIES_TAG_PulseSeries       9
#define TSERIES_TAG_TriangleSeries       10
#define TSERIES_TAG_PeerMotion       11
#define TSERIES_TAG_PeerNGAMotion       12
#define TSERIES_TAG_PathTimeSeriesThermal  13  //L.Jiang [ SIF ]
#define TSERIES_TAG_RampSeries  14  //CDM
#define TSERIES_TAG_MPAccSeries       15 //Tang.S[SEU]

#define PARAMETER_TAG_Parameter			   1
#define PARAMETER_TAG_MaterialStageParameter       2
#define PARAMETER_TAG_MatParameter                 3
#define PARAMETER_TAG_InitialStateParameter        4
#define PARAMETER_TAG_ElementStateParameter        5
#define PARAMETER_TAG_ElementParameter             6



#define MAT_TAG_ElasticMaterial			 1
#define MAT_TAG_ElasticPPMaterial		 2
#define MAT_TAG_ParallelMaterial		 3
#define MAT_TAG_Concrete01			 4
#define MAT_TAG_Concrete01A			 444
#define MAT_TAG_Steel01				 5
#define MAT_TAG_Hardening			 6
#define MAT_TAG_Hysteretic			 7
#define MAT_TAG_HystereticSM		1969  // Silvia Mazzoni, 2022
#define MAT_TAG_EPPGap				 8
#define MAT_TAG_Viscous				 9
#define MAT_TAG_Backbone			10
#define MAT_TAG_PathIndependent			11
#define MAT_TAG_Multiplier			111
#define MAT_TAG_SeriesMaterial			12
#define MAT_TAG_CableMaterial                   13
#define MAT_TAG_ENTMaterial			14
#define MAT_TAG_Penalty				15
#define MAT_TAG_MinMax				16
#define MAT_TAG_TensionOnly			1601
#define MAT_TAG_BoucWen				17
#define MAT_TAG_Pinching4			18
#define MAT_TAG_BarSlip				19
#define MAT_TAG_Fatigue			        20
#define MAT_TAG_SmoothSteel01			21
#define MAT_TAG_SmoothConcrete01		22
#define MAT_TAG_Steel03				23
#define MAT_TAG_ReinforcingSteel		24
#define MAT_TAG_Concrete02			25
#define MAT_TAG_Steel02				26
#define MAT_TAG_Bond_SP01                       27
#define MAT_TAG_Hysteretic2			28
#define MAT_TAG_SteelDRC			29
#define MAT_TAG_Concrete04                      30
#define MAT_TAG_SecantConcrete                  31
#define MAT_TAG_ContinuumUniaxial               32
#define MAT_TAG_Concrete05                      33
#define MAT_TAG_Concrete06                      34
#define MAT_TAG_Concrete07                      37
#define MAT_TAG_HyperbolicGapMaterial           38
#define MAT_TAG_ImpactMaterial                  39
#define MAT_TAG_Hertzdamp                  420
#define MAT_TAG_JankowskiImpact            421
#define MAT_TAG_ViscoelasticGap            422
#define MAT_TAG_ShearPanelMaterial		40
#define MAT_TAG_SAWSMaterial			41
#define MAT_TAG_ConcreteL01			42
#define MAT_TAG_ConcreteZ01			43
#define MAT_TAG_TendonL01			44
#define MAT_TAG_SteelZ01			45
#define MAT_TAG_ElasticMultiLinear		46
#define MAT_TAG_InitStrain			47
#define MAT_TAG_InitStress			48
#define MAT_TAG_pyUCLA  			49
#define MAT_TAG_Maxwell			        50
#define MAT_TAG_Cast			        51
#define MAT_TAG_MultiLinear			52
#define MAT_TAG_ElasticBilin			53
#define MAT_TAG_SMA                             54
#define MAT_TAG_SelfCentering                   55
#define MAT_TAG_Clough1	                        56
#define MAT_TAG_Clough2	                        57
#define MAT_TAG_Pinch1	                        58
#define MAT_TAG_BiLinear                        59
#define MAT_TAG_Pinching                        60
#define MAT_TAG_HookGap 			61
#define MAT_TAG_FRPConfinedConcrete             62
#define MAT_TAG_Steel01Thermal		        63
#define MAT_TAG_Steel02Thermal			64
#define MAT_TAG_Concrete02Thermal		65
#define MAT_TAG_ModIMKPinching                  66
#define MAT_TAG_ModIMKPeakOriented              67
#define MAT_TAG_RambergOsgoodSteel              68
#define MAT_TAG_PinchingLimitStateMaterial      69	
#define MAT_TAG_BraceMaterial                   70	
#define MAT_TAG_ViscousDamper                   71
#define MAT_TAG_ConcretewBeta                   72
#define MAT_TAG_WrapperUniaxialMaterial         73
#define MAT_TAG_UniaxialJ2Plasticity            74
#define MAT_TAG_BWBN                            75
#define MAT_TAG_OriginCentered                  76
#define MAT_TAG_Steel2                          77
#define MAT_TAG_DoddRestr                       78
#define MAT_TAG_ConcreteSakaiKawashima          79
#define MAT_TAG_ResilienceMaterialHR            80
#define MAT_TAG_CFSSSWP                         81
#define MAT_TAG_CFSWSWP                         82
#define MAT_TAG_ResilienceLow                   83
#define MAT_TAG_Bilin02                         84
#define MAT_TAG_ModIMKPeakOriented02            85
#define MAT_TAG_ModIMKPinching02                86
#define MAT_TAG_Steel4                          87
#define MAT_TAG_SimpleFractureMaterial          88
#define MAT_TAG_BilinearOilDamper               89
#define MAT_TAG_ConcreteCM                      90
#define MAT_TAG_SteelMPF                        91
#define MAT_TAG_ElasticMaterialThermal          92   //L.Jiang [SIF]
#define MAT_TAG_SteelECThermal                  93   //L.Jiang [SIF]
#define MAT_TAG_StainlessECThermal              94   //L.Jiang [SIF]
#define MAT_TAG_ConcreteECThermal               95   //L.Jiang [SIF]
#define MAT_TAG_BoucWenOriginal                 96
#define MAT_TAG_DamperMaterial                  97
#define MAT_TAG_SPSW02                          98   //SAJalali
#define MAT_TAG_Steel02Fatigue                  99 //nassermarafi
#define MAT_TAG_Concrete02IS                    100 //nassermarafi
#define MAT_TAG_ConfinedConcrete01              101
#define MAT_TAG_ElasticPowerFunc                102
#define MAT_TAG_UVCuniaxial                     103
#define MAT_TAG_IMKBilin                        104
#define MAT_TAG_IMKPeakOriented                 105
#define MAT_TAG_IMKPinching                     106
#define MAT_TAG_SLModel                         107
#define MAT_TAG_PySimple1                    205
#define MAT_TAG_TzSimple1                    206
#define MAT_TAG_QzSimple1                    207
#define MAT_TAG_PyLiq1                       208
#define MAT_TAG_TzLiq1                       209
#define MAT_TAG_QzLiq1                       210
#define MAT_TAG_PySimple2                    211
#define MAT_TAG_TzSimple2                    212
#define MAT_TAG_QzSimple2                    213
#define MAT_TAG_SteelBRB                     214
#define MAT_TAG_PySimple3                    215
#define MAT_TAG_PlateBearingConnectionThermal 216
#define MAT_TAG_ASD_SMA_3K                    217
#define MAT_TAG_SteelFractureDI			218 // galvisf
#define MAT_TAG_Masonry 219
#define MAT_TAG_Masonryt 220
#define MAT_TAG_Trilinwp 221
#define MAT_TAG_Trilinwp2 222
#define MAT_TAG_Trilinwpd 223
#define MAT_TAG_TDConcrete 224
#define MAT_TAG_TDConcreteNL 2240
#define MAT_TAG_TDConcreteEXP 225
#define MAT_TAG_TDConcreteMC10 226
#define MAT_TAG_TDConcreteMC10NL 227
#define MAT_TAG_CoulombDamperMaterial 228
#define MAT_TAG_FlagShapeMaterial 229
#define MAT_TAG_CreepMaterial 230

#define MAT_TAG_FedeasMaterial    1000
#define MAT_TAG_FedeasBond1       1001
#define MAT_TAG_FedeasBond2       1002
#define MAT_TAG_FedeasConcrete1   1003
#define MAT_TAG_FedeasConcrete2   1004
#define MAT_TAG_FedeasConcrete3   1005
#define MAT_TAG_FedeasHardening   1006
#define MAT_TAG_FedeasHysteretic1 1007
#define MAT_TAG_FedeasHysteretic2 1008
#define MAT_TAG_FedeasSteel1      1009
#define MAT_TAG_FedeasSteel2      1010
#define MAT_TAG_PlasticDamage	  1011

#define MAT_TAG_LimitState	   1972
#define MAT_TAG_Elastic2Material   1973
#define MAT_TAG_FRCC		   1980		   		
#define MAT_TAG_DrainMaterial		2000
#define MAT_TAG_DrainHardening		2001
#define MAT_TAG_DrainBilinear		2002
#define MAT_TAG_DrainClough1		2003
#define MAT_TAG_DrainClough2		2004
#define MAT_TAG_DrainPinch1			2005
#define MAT_TAG_DrainPinch2			2006
#define MAT_TAG_Bilin		2007

#define MAT_TAG_SnapMaterial		3000
#define MAT_TAG_SnapBilinear		3001
#define MAT_TAG_SnapClough		3002
#define MAT_TAG_SnapPinch		3003
#define MAT_TAG_SnapCloughDamage	3004
#define MAT_TAG_SnapPinchingDamage	3005

#define MAT_TAG_ECC01 3010
#define MAT_TAG_Concrete01WithSITC 3011

#define MAT_TAG_KikuchiAikenHDR 6102
#define MAT_TAG_KikuchiAikenLRB 6105
#define MAT_TAG_AxialSp   6111
#define MAT_TAG_AxialSpHD 6112

#define MAT_TAG_HystereticPoly 6113			// Salvatore Sessa 14-Jan-2021 Mail: salvatore.sessa2@unina.it
#define MAT_TAG_DowelType  6114

#define MAT_TAG_DuctileFracture 6115 // Kuanshi Zhong

#define MAT_TAG_HystereticSmooth 6116			// Salvatore Sessa 19-Apr-2022 Mail: salvatore.sessa2@unina.it
#define MAT_TAG_HystereticAsym 6117			// Salvatore Sessa 19-Apr-2022 Mail: salvatore.sessa2@unina.it

#define ND_TAG_ExternalNDMaterial 999901
#define MAT_TAG_ExternalUniaxialMaterial 999901

#define MAT_TAG_BoucWenInfill  6666    // Stefano Sirotti 09-Feb-2022 stefano.sirotti@unimore.it

#define MAT_TAG_GMG_CyclicReinforcedConcrete    9999    // Rasool Ghorbani


// GNG material - J.Cook UCanterbury
#define MAT_TAG_GNG 7001

#define SEC_TAG_Elastic2d                        3
#define SEC_TAG_Elastic3d                        4
#define SEC_TAG_Generic1d	                 5
#define SEC_TAG_GenericNd	                 6
#define SEC_TAG_Aggregator	                 7
#define SEC_TAG_Parallel	                 77
#define SEC_TAG_Fiber		                 8
#define SEC_TAG_FiberSection2d		         9
#define SEC_TAG_NDFiberSection2d		         900
#define SEC_TAG_FiberSection3d		        10
#define SEC_TAG_FiberSectionWarping3d		        1010
#define SEC_TAG_FiberSectionAsym3d		        1011
#define SEC_TAG_NDFiberSection3d		         1000
#define SEC_TAG_FiberSectionGJ		        11
#define SEC_TAG_BeamFiberSection	        12
#define SEC_TAG_ElasticPlateSection	        13
#define SEC_TAG_ElasticMembranePlateSection	14
#define SEC_TAG_MembranePlateFiberSection	15
#define SEC_TAG_Bidirectional	                16
#define SEC_TAG_WSection2d	                17
#define SEC_TAG_Isolator2spring                 18
#define SEC_TAG_SoilFooting2d                   19
#define SEC_TAG_YieldSurface2d                  20
#define SEC_TAG_YieldSurface2D02                21
#define SEC_TAG_YieldSurface2D01                22
#define SEC_TAG_ElasticShear2d                  23
#define SEC_TAG_ElasticBDShear2d 2376
#define SEC_TAG_ElasticShear3d                  24
#define SEC_TAG_FiberSection2dInt		25
#define SEC_TAG_FiberSection2dThermal		26
#define SEC_TAG_LayeredShellFiberSection        27
#define SEC_TAG_ElasticWarpingShear2d           28
#define SEC_TAG_DoubleMembranePlateFiberSection 29
#define SEC_TAG_NDFiberSectionWarping2d         30
#define SEC_TAG_Elliptical2                     31
#define SEC_TAG_FiberSection3dThermal           32   // L.Jiang[SIF]
#define SEC_TAG_FiberSectionGJThermal           33   // L.Jiang[SIF]
#define SEC_TAG_MembranePlateFiberSectionThermal 34  // L.Jiang[SIF]
#define SEC_TAG_LayeredShellFiberSectionThermal 35     //L.Jiang[SIF]
#define SEC_TAG_BiaxialHysteretic 36
#define SEC_TAG_ElasticTube3d 37
#define SEC_TAG_CreepSection 38

#define SEC_TAG_MCFTFiberSection2d 7601

#define SEC_TAG_ReinforcedConcreteLayeredMembraneSection 7701 // M. J. Nunez - UChile
#define SEC_TAG_LayeredMembraneSection 7702 // M. J. Nunez - UChile
#define SEC_TAG_ElasticMembraneSection 7703 // M. J. Nunez - UChile

#define SECTION_INTEGRATION_TAG_WideFlange 1
#define SECTION_INTEGRATION_TAG_RC 2
#define SECTION_INTEGRATION_TAG_RCT 3
#define SECTION_INTEGRATION_TAG_RCTUM 4
#define SECTION_INTEGRATION_TAG_RCCIRCULAR 5
#define SECTION_INTEGRATION_TAG_RCTUNNEL 6
#define SECTION_INTEGRATION_TAG_Tube 7
#define SECTION_INTEGRATION_TAG_HSS 8

#define ND_TAG_WrapperNDMaterial		9
#define ND_TAG_ElasticIsotropic			10
#define ND_TAG_ElasticIsotropicPlaneStrain2d	11
#define ND_TAG_ElasticIsotropicPlaneStress2d	12
#define ND_TAG_ElasticIsotropicAxiSymm          13
#define ND_TAG_ElasticIsotropicPlateFiber	14
#define ND_TAG_ElasticIsotropicBeamFiber	15
#define ND_TAG_ElasticIsotropicThreeDimensional 16
#define ND_TAG_ElasticCrossAnisotropic3D        17
#define ND_TAG_ElasticIsotropicBeamFiber2d	18
#define ND_TAG_CycLiqCP3D                       19
#define ND_TAG_CycLiqCPPlaneStrain              20
#define ND_TAG_PressureDependentElastic3D       22
#define ND_TAG_Damage2p 			23
#define ND_TAG_Damage2p3D 			24
#define ND_TAG_Damage2ppstress 			25
#define ND_TAG_SimplifiedJ2                     26
#define ND_TAG_CapPlasticity                    27
#define ND_TAG_PlaneStressUserMaterial          28
#define ND_TAG_PlateFromPlaneStressMaterial     29
#define ND_TAG_PlateRebarMaterial               30
#define ND_TAG_ElasticOrthotropic		  31
#define ND_TAG_ElasticOrthotropicPlaneStrain2d	  32
#define ND_TAG_ElasticOrthotropicPlaneStress2d	  33
#define ND_TAG_ElasticOrthotropicAxiSymm          34
#define ND_TAG_ElasticOrthotropicPlateFiber	  35
#define ND_TAG_ElasticOrthotropicBeamFiber	  36
#define ND_TAG_ElasticOrthotropicThreeDimensional 37
#define ND_TAG_ElasticOrthotropicBeamFiber2d	  38
#define ND_TAG_CycLiqCPSP3D                       39
#define ND_TAG_CycLiqCPSPPlaneStrain              40
#define ND_TAG_ConcreteS                          41
#define ND_TAG_MaterialCMM                        42
#define ND_TAG_FSAM                               43
#define ND_TAG_PlasticDamageConcrete3d            44
#define ND_TAG_PlaneStressLayeredMaterial         45
#define ND_TAG_PlaneStressRebarMaterial           46
#define ND_TAG_Faria1998PlaneStrain               48
#define ND_TAG_Faria1998PlaneStress               49
#define ND_TAG_Faria1998PlaneStress2d               50
#define ND_TAG_Faria1998               51
#define ND_TAG_Faria1998ThreeDimensional               52
#define ND_TAG_CPlaneStress   53
#define ND_TAG_CPlaneStrain   54
#define ND_TAG_CPlaneStress2d   55
#define ND_TAG_CThreeDimensional   55
#define ND_TAG_StressDensityModel2D 56
#define ND_TAG_StressDensityModel3D 57
#define ND_TAG_UVCmultiaxial  58
#define ND_TAG_UVCplanestress 59
#define ND_TAG_LinearCap 60

#define ND_TAG_LowTension 65
#define ND_TAG_LowTensionPlaneStress 66
#define ND_TAG_ExponentialTS 67
#define ND_TAG_ExponentialTS2D 68
#define ND_TAG_ElasticTS 69
#define ND_TAG_ElasticTS2D 70
#define ND_TAG_BilinearTS 71
#define ND_TAG_BilinearTS2D 72

#define ND_TAG_FluidSolidPorousMaterial        100
#define ND_TAG_PressureDependMultiYield		101
#define ND_TAG_PressureIndependMultiYield		102
#define ND_TAG_PressureDependMultiYield02		103
#define ND_TAG_ReinforcedConcretePlaneStress  104
#define ND_TAG_FAReinforcedConcretePlaneStress  105
#define ND_TAG_FAFourSteelRCPlaneStress  106
#define ND_TAG_RAFourSteelRCPlaneStress  107
#define ND_TAG_PrestressedConcretePlaneStress  108
#define ND_TAG_FAPrestressedConcretePlaneStress  109
#define ND_TAG_FAFourSteelPCPlaneStress  110
#define ND_TAG_RAFourSteelPCPlaneStress  111
#define ND_TAG_PressureDependMultiYield03		112

#define ND_TAG_J2PlaneStrain                  3005
#define ND_TAG_J2PlaneStress                  3006
#define ND_TAG_J2AxiSymm                      3007
#define ND_TAG_J2ThreeDimensional             3009
#define ND_TAG_J2PlateFiber		      3010
#define ND_TAG_J2BeamFiber		      3011
#define ND_TAG_J2BeamFiber2d 91625
#define ND_TAG_J2BeamFiber3d 92516

#define ND_TAG_FeapMaterial                 1000
#define ND_TAG_FeapMaterial01                 1001
#define ND_TAG_FeapMaterial02                 1002
#define ND_TAG_FeapMaterial03                 1003
#define ND_TAG_PlaneStressMaterial          2000
#define ND_TAG_PlateFiberMaterial          2001
#define ND_TAG_PlaneStrainMaterial          2003
#define ND_TAG_BeamFiberMaterial		2002
#define ND_TAG_BeamFiberMaterial2d		2004
#define ND_TAG_BeamFiberMaterial2dPS		2005
#define ND_TAG_OrthotropicMaterial		2006
#define ND_TAG_Series3DMaterial		2007
#define ND_TAG_CompressibleFluid		3001
#define ND_TAG_GeneralizedPlasticity 3002
#define ND_TAG_J2Plasticity02  3003
#define ND_TAG_FiniteDeformationElastic3D	8002
#define ND_TAG_NeoHookeanCompressible3D	        8003
#define ND_TAG_FDdecoupledElastic3D	        8004
#define ND_TAG_FiniteDeformationEP3D	        8005
// Contact Material - P.Arduino
#define ND_TAG_ContactMaterial2D				14001
#define ND_TAG_ContactMaterial3D				14002
// Drucker-Prager - P.Arduino
#define ND_TAG_DruckerPrager					14003
#define ND_TAG_DruckerPragerThreeDimensional	14004
#define ND_TAG_DruckerPragerTensionCutoff		14005
#define ND_TAG_DruckerPrager3D	                14006
#define ND_TAG_DruckerPragerPlaneStrain         14007
// CamClay with Bounding Surface - C.McGann
#define ND_TAG_BoundingCamClay                  14008
#define ND_TAG_BoundingCamClay3D                14009
#define ND_TAG_BoundingCamClayPlaneStrain       14010
// Initial state analysis material wrapper - C.McGann
#define ND_TAG_InitialStateAnalysisWrapper      14011
// Manzari Dafalias material - P. Arduino
#define ND_TAG_ManzariDafalias                  14012
#define ND_TAG_ManzariDafalias3D                14013
#define ND_TAG_ManzariDafaliasPlaneStrain       14014
// Manzari Dafalias material - A. Ghofrani
#define ND_TAG_ManzariDafaliasRO                14015
#define ND_TAG_ManzariDafalias3DRO              14016
#define ND_TAG_ManzariDafaliasPlaneStrainRO     14017
// Stress Density material - C.McGann
#define ND_TAG_stressDensity                  14018
// PM4Sand material - L.Chen
#define ND_TAG_PM4Sand                        14021
// PM4Silt material - L.Chen
#define ND_TAG_PM4Silt                        14022
// J2CyclicBoundingSurface material - P. Arduino,  D.Turello
#define ND_TAG_J2CyclicBoundingSurface            14023
#define ND_TAG_J2CyclicBoundingSurface3D          14024
#define ND_TAG_J2CyclicBoundingSurfacePlaneStrain 14025

// MultiaxialCyclicPlasticity, add by Gang Wang
#define ND_TAG_MultiaxialCyclicPlasticity             10031
#define ND_TAG_MultiaxialCyclicPlasticity3D           10032
#define ND_TAG_MultiaxialCyclicPlasticityAxiSymm      10033
#define ND_TAG_MultiaxialCyclicPlasticityPlaneStrain  10034

#define ND_TAG_ConcreteMcftNonLinear5 7601
#define ND_TAG_ConcreteMcftNonLinear7 7602

#define ND_TAG_ElasticIsotropicThermal	      7000   //L.Jiang[SIF]
#define ND_TAG_ElasticIsotropic3DThermal      7001   //L.Jiang[SIF]
#define ND_TAG_J2ThreeDimensionalThermal      7002   //L.Jiang[SIF]
#define ND_TAG_DruckerPragerThermal	      7003   //L.Jiang[SIF]
#define ND_TAG_DruckerPrager3DThermal         7004   //L.Jiang[SIF]
#define ND_TAG_PlasticDamageConcrete3dThermal 7005   //L.Jiang[SIF]
#define ND_TAG_PlateFiberMaterialThermal      7006   //L.Jiang[SIF]
#define ND_TAG_PlateRebarMaterialThermal      7007   //L.Jiang[SIF]
#define ND_TAG_PlateFromPlaneStressMaterialThermal 7008   //L.Jiang[SIF]

#define ND_TAG_InitStressNDMaterial 7009

#define ND_TAG_IncrementalElasticIsotropicThreeDimensional 7010 //Chile

#define ND_TAG_SAniSandMS 7011 //UANDES - Chile
#define ND_TAG_SAniSandMSPlaneStrain 7012 //UANDES - Chile
#define ND_TAG_SAniSandMS3D 7013 //UANDES - Chile

#define ND_TAG_ElasticPlaneStress 7014
#define ND_TAG_ElasticOrthotropicPlaneStress 7015
#define ND_TAG_VonPapaDamage 7016

#define ND_TAG_ASDConcrete3DMaterial 7017 // Massimo Petracca ASDEA Software

#define ND_TAG_OrthotropicRotatingAngleConcreteT2DMaterial01 7018 // M. J. Nunez - UChile
#define ND_TAG_SmearedSteelDoubleLayerT2DMaterial01 7019		  // M. J. Nunez - UChile

#define FIBER_TAG_Uniaxial2d	1
#define FIBER_TAG_Uniaxial3d	2
#define FIBER_TAG_ND2d	3
#define FIBER_TAG_ND3d	4

#define BACKBONE_TAG_Capped		1
#define BACKBONE_TAG_LinearCapped	2
#define BACKBONE_TAG_Material		3
#define BACKBONE_TAG_Arctangent		4
#define BACKBONE_TAG_Trilinear		5
#define BACKBONE_TAG_Multilinear	6
#define BACKBONE_TAG_Mander		7
#define BACKBONE_TAG_KentPark		8
#define BACKBONE_TAG_Raynor		9
#define BACKBONE_TAG_ReeseStiffClayBelowWS 10
#define BACKBONE_TAG_ReeseSoftClay      11
#define BACKBONE_TAG_ReeseSand          12
#define BACKBONE_TAG_ReeseStiffClayAboveWS 13
#define BACKBONE_TAG_VuggyLimestone 14
#define BACKBONE_TAG_CementedSoil 15
#define BACKBONE_TAG_WeakRock 16
#define BACKBONE_TAG_LiquefiedSand 17


#define DEG_TAG_STIFF_Constant		1
#define DEG_TAG_STIFF_Ductility		2
#define DEG_TAG_STIFF_Energy		3
#define DEG_TAG_STIFF_Pincheira		4

#define DEG_TAG_UNLOAD_Constant		1
#define DEG_TAG_UNLOAD_Takeda		2
#define DEG_TAG_UNLOAD_Energy		3
#define DEG_TAG_UNLOAD_Karsan		4

#define DEG_TAG_STRENGTH_ACI		1
#define DEG_TAG_STRENGTH_Constant	2
#define DEG_TAG_STRENGTH_Ductility	3
#define DEG_TAG_STRENGTH_Petrangeli	4
#define DEG_TAG_STRENGTH_Energy		5
#define DEG_TAG_STRENGTH_Section	6

#define PATTERN_TAG_LoadPattern           1
#define PATTERN_TAG_MultiSupportPattern	  3
#define PATTERN_TAG_UniformExcitation     2
#define PATTERN_TAG_FirePattern           3
#define PATTERN_TAG_PBowlLoading          4
#define PATTERN_TAG_DRMLoadPattern        5
#define PATTERN_TAG_H5DRM                 6

#define LOAD_TAG_Beam2dUniformLoad        3
#define LOAD_TAG_Beam2dPointLoad          4
#define LOAD_TAG_Beam3dUniformLoad        5
#define LOAD_TAG_Beam3dPointLoad          6
#define LOAD_TAG_BrickSelfWeight          7
#define LOAD_TAG_Beam2dTempLoad           8
#define LOAD_TAG_SurfaceLoader            9 // C.McGann, U.W.
#define LOAD_TAG_SelfWeight              10 // C.McGann, U.W.
#define LOAD_TAG_Beam2dThermalAction      11
#define LOAD_TAG_Beam2dPartialUniformLoad 12
#define LOAD_TAG_Beam3dPartialUniformLoad 121
#define LOAD_TAG_Beam3dThermalAction      13 // L.Jiang [ SIF ]
#define LOAD_TAG_ShellThermalAction       14 // L.Jiang [ SIF ]
#define LOAD_TAG_NodalThermalAction       15 //L.Jiang [ SIF ]
#define LOAD_TAG_ThermalActionWrapper     16 //L.Jiang [ SIF ]
#define LOAD_TAG_LysmerVelocityLoader      17  //Jose Abell (UANDES)
#define LOAD_TAG_IGAFollowerLoad      18  //Jose Abell (UANDES)


#define MAT_TAG_IsotropicLinElastic         1001
#define MAT_TAG_IsotropicLinElasticPoint    1002
#define MAT_TAG_OrthotropicLinElastic       1003
#define MAT_TAG_OrthotropicLinElasticPoint  1004

#define ELE_TAG_Subdomain     	         1
#define ELEMENT_TAGS_WrapperElement      2
#define ELE_TAG_ElasticBeam2d            3
#define ELE_TAG_ModElasticBeam2d         4
#define ELE_TAG_ModElasticBeam3d         41234
#define ELE_TAG_ElasticBeam3d            5
#define ELE_TAG_ElasticBeamWarping3d            5001
#define ELE_TAG_Beam2d    	         6
#define ELE_TAG_beam2d02    	         7
#define ELE_TAG_beam2d03    	         8
#define ELE_TAG_beam2d04    	         9
#define ELE_TAG_beam3d01    	        10
#define ELE_TAG_beam3d02    	        11
#define ELE_TAG_Truss    	        12
#define ELE_TAG_TrussSection            13
#define ELE_TAG_CorotTruss    	        14
#define ELE_TAG_CorotTrussSection    	15
#define ELE_TAG_fElmt05	                16
#define ELE_TAG_fElmt02	                17
#define ELE_TAG_MyTruss    	        18
#define ELE_TAG_ZeroLength	        19
#define ELE_TAG_ZeroLengthSection	20
#define ELE_TAG_ZeroLengthND	        21
#define ELE_TAG_ZeroLengthContact2D	22
#define ELE_TAG_ZeroLengthContact3D	23
#define ELE_TAG_ZeroLengthContactNTS2D	24
#define ELE_TAG_ZeroLengthInterface2D	25
#define ELE_TAG_CoupledZeroLength	26
#define ELE_TAG_BiaxialZeroLength	2626
#define ELE_TAG_ZeroLengthRocking       27
#define ELE_TAG_NLBeamColumn2d	        28
#define ELE_TAG_NLBeamColumn3d	        29
#define ELE_TAG_LargeDispBeamColumn3d	30
#define ELE_TAG_FourNodeQuad	        31
#define ELE_TAG_FourNodeQuad3d	        32
#define ELE_TAG_Tri31	                33    //Added by Roozbeh Geraili Mikola
#define ELE_TAG_BeamWithHinges2d        34
#define ELE_TAG_BeamWithHinges3d        35
#define ELE_TAG_EightNodeBrick          36
#define ELE_TAG_TwentyNodeBrick         37
#define ELE_TAG_EightNodeBrick_u_p_U    38
#define ELE_TAG_TwentyNodeBrick_u_p_U   39
#define ELE_TAG_FourNodeQuadUP          40
#define ELE_TAG_TotalLagrangianFD20NodeBrick 41
#define ELE_TAG_TotalLagrangianFD8NodeBrick  42
#define ELE_TAG_EightNode_LDBrick_u_p        43
#define ELE_TAG_EightNode_Brick_u_p     44
#define ELE_TAG_TwentySevenNodeBrick    45
#define ELE_TAG_BrickUP                 46
#define ELE_TAG_Nine_Four_Node_QuadUP   47
#define ELE_TAG_Twenty_Eight_Node_BrickUP    48
#define ELE_TAG_Twenty_Node_Brick       49
#define ELE_TAG_BBarFourNodeQuadUP      50
#define ELE_TAG_BBarBrickUP             51
#define ELE_TAG_PlateMITC4              52
#define ELE_TAG_ShellMITC4              53
#define ELE_TAG_ShellMITC9              54 //Tesser
#define ELE_TAG_Plate1                  55
#define ELE_TAG_Brick                   56
#define ELE_TAG_BbarBrick               57
#define ELE_TAG_FLBrick                 58
#define ELE_TAG_EnhancedQuad            59
#define ELE_TAG_ConstantPressureVolumeQuad 60
#define ELE_TAG_NineNodeMixedQuad          61
#define ELE_TAG_DispBeamColumn2d        62
#define ELE_TAG_DispBeamColumnNL2d        621
#define ELE_TAG_TimoshenkoBeamColumn2d  63
#define ELE_TAG_TimoshenkoBeamColumn3d  631
#define ELE_TAG_DispBeamColumn3d        64
#define ELE_TAG_DispBeamColumnNL3d        640
#define ELE_TAG_DispBeamColumnWarping3d        641
#define ELE_TAG_DispBeamColumnAsym3d           642
#define ELE_TAG_HingedBeam2d            65
#define ELE_TAG_HingedBeam3d            66
#define ELE_TAG_TwoPointHingedBeam2d    67
#define ELE_TAG_TwoPointHingedBeam3d    68
#define ELE_TAG_OnePointHingedBeam2d    69
#define ELE_TAG_OnePointHingedBeam3d    70
#define ELE_TAG_BeamColumnJoint2d       71
#define ELE_TAG_BeamColumnJoint3d       72
#define ELE_TAG_ForceBeamColumn2d       73
#define ELE_TAG_ForceBeamColumnWarping2d 731
#define ELE_TAG_ForceBeamColumn3d       74
#define ELE_TAG_ElasticForceBeamColumn2d 75
#define ELE_TAG_ElasticForceBeamColumnWarping2d 751
#define ELE_TAG_ElasticForceBeamColumn3d 76
#define ELE_TAG_ForceBeamColumnCBDI2d   77
#define ELE_TAG_ForceBeamColumnCBDI3d   78
#define ELE_TAG_MixedBeamColumn2d 30766
#define ELE_TAG_MixedBeamColumn3d 30765
#define ELE_TAG_MixedBeamColumnAsym3d 30767
#define ELE_TAG_DispBeamColumn2dInt     79
#define ELE_TAG_InternalSpring          80
#define ELE_TAG_SimpleJoint2D           81
#define ELE_TAG_LehighJoint2d           8181
#define ELE_TAG_Joint2D                 82
#define ELE_TAG_Joint3D                 83
#define ELE_TAG_ElastomericBearingPlasticity3d 84
#define ELE_TAG_ElastomericBearingPlasticity2d 85
#define ELE_TAG_TwoNodeLink             86
#define ELE_TAG_ActuatorCorot           87
#define ELE_TAG_Actuator                88
#define ELE_TAG_Adapter                 89
#define ELE_TAG_ElastomericBearingBoucWen2d 90
#define ELE_TAG_ElastomericBearingBoucWen3d 91
#define ELE_TAG_FlatSliderSimple2d      92
#define ELE_TAG_FlatSliderSimple3d      93
#define ELE_TAG_FlatSlider2d            94
#define ELE_TAG_FlatSlider3d            95
#define ELE_TAG_SingleFPSimple2d        96
#define ELE_TAG_SingleFPSimple3d        97
#define ELE_TAG_SingleFP2d              98
#define ELE_TAG_SingleFP3d              99
#define ELE_TAG_DoubleFPSimple2d       100
#define ELE_TAG_DoubleFPSimple3d       101
#define ELE_TAG_DoubleFP2d             102
#define ELE_TAG_DoubleFP3d             103
#define ELE_TAG_TripleFPSimple2d       104
#define ELE_TAG_TripleFPSimple3d       105
#define ELE_TAG_TripleFP2d             106
#define ELE_TAG_TripleFP3d             107
#define ELE_TAG_MultiFP2d              108
#define ELE_TAG_MultiFP3d              109
#define ELE_TAG_GenericClient          110
#define ELE_TAG_GenericCopy            111
#define ELE_TAG_PY_MACRO2D             112
// elements added by U.W. - P.Arduino
#define ELE_TAG_SimpleContact2D        113
#define ELE_TAG_SimpleContact3D        114
#define ELE_TAG_BeamContact3D          115
#define ELE_TAG_SurfaceLoad            116
#define ELE_TAG_BeamContact2D          117
#define ELE_TAG_BeamEndContact3D       118
#define ELE_TAG_SSPquad                119
#define ELE_TAG_SSPquadUP              120
#define ELE_TAG_SSPbrick               121
#define ELE_TAG_SSPbrickUP             122
#define ELE_TAG_BeamContact2Dp         123
#define ELE_TAG_BeamContact3Dp         124
#define ELE_TAG_BeamEndContact3Dp      125
#define ELE_TAG_Quad4FiberOverlay      126
#define ELE_TAG_Brick8FiberOverlay     127
#define ELE_TAG_DispBeamColumn2dThermal 128
#define ELE_TAG_TPB1D                  129
#define ELE_TAG_TFP_Bearing            130
#define ELE_TAG_TFP_Bearing2d          131
#define ELE_TAG_TripleFrictionPendulum 132
#define ELE_TAG_PFEMElement2D          133
#define ELE_TAG_FourNodeQuad02         134
#define ELE_TAG_cont2d01    	       135	// provisional
#define ELE_TAG_cont2d02    	       136 	// provisional
#define ELE_TAG_CST	    	           137
#define ELE_TAG_Truss2                 138
#define ELE_TAG_CorotTruss2            139
#define ELE_Tag_ZeroLengthImpact3D     140
#define ELE_TAG_PFEMElement3D          141
#define ELE_TAG_PFEMElement2DCompressible 142
#define ELE_TAG_PFEMElement2DBubble       143
#define ELE_TAG_PFEMElement2Dmini         144
#define ELE_TAG_ElasticTimoshenkoBeam2d   145
#define ELE_TAG_ElasticTimoshenkoBeam3d   146
#define ELE_TAG_ElastomericBearingUFRP2d  147
#define ELE_TAG_ElastomericBearingUFRP3d  148
#define ELE_TAG_RJWatsonEQS2d             149
#define ELE_TAG_RJWatsonEQS3d             150
#define ELE_TAG_HDR                       151
#define ELE_TAG_ElastomericX              152
#define ELE_TAG_LeadRubberX               153
#define ELE_TAG_PileToe3D                 154
#define ELE_TAG_N4BiaxialTruss            155
#define ELE_TAG_ShellDKGQ                 156
#define ELE_TAG_ShellNLDKGQ               157
#define ELE_TAG_MultipleShearSpring       158
#define ELE_TAG_MultipleNormalSpring      159
#define ELE_TAG_KikuchiBearing            160
#define ELE_TAG_YamamotoBiaxialHDR        161
#define ELE_TAG_MVLEM                     162 // Kristijan Kolozvari
#define ELE_TAG_SFI_MVLEM                 163 // Kristijan Kolozvari
#define ELE_TAG_PFEMElement2DFIC          164
#define ELE_TAG_ElastomericBearingBoucWenMod3d 165
#define ELE_TAG_FPBearingPTV              166
#define ELE_TAG_ShellDKGT                 167
#define ELE_TAG_ShellNLDKGT               168
#define ELE_TAG_CatenaryCable             169
#define ELE_TAG_DispBeamColumn3dThermal   170  // L.Jiang [SIF]
#define ELE_TAG_ForceBeamColumn2dThermal  171  //L.Jiang [SIF]
#define ELE_TAG_ForceBeamColumn3dThermal  172  //L.Jiang [SIF] //Still testing
#define ELE_TAG_ShellMITC4Thermal         173   //L.Jiang [SIF]
#define ELE_TAG_ShellNLDKGQThermal        174   //L.Jiang [SIF]
#define ELE_TAG_ShellANDeS                175 //by jaabell (UANDES)
#define ELE_TAG_AxEqDispBeamColumn2d      178
#define ELE_TAG_FourNodeTetrahedron       179 //by jaabell (UANDES)
#define ELE_TAG_TriSurfaceLoad            180 //by jaabell (UANDES) 
#define ELE_TAG_QuadBeamEmbedContact      181
#define ELE_TAG_EmbeddedBeamInterfaceL    182
#define ELE_TAG_EmbeddedBeamInterfaceP    183
#define ELE_TAG_EmbeddedEPBeamInterface   184
#define ELE_TAG_LysmerTriangle            185
#define ELE_TAG_TaylorHood2D              186
#define ELE_TAG_PFEMElement2DQuasi        187
#define ELE_TAG_MINI                      188
#define ELE_TAG_PFEMElement3DBubble       189
#define ELE_TAG_LinearElasticSpring       190
#define ELE_TAG_Inerter                   191
#define ELE_TAG_GradientInelasticBeamColumn2d	192
#define ELE_TAG_GradientInelasticBeamColumn3d	193
#define ELE_TAG_CohesiveZoneQuad 194
#define ELE_TAG_ComponentElement2d       195
#define ELE_TAG_ComponentElement3d       195195
#define ELE_TAG_InerterElement 196
#define ELE_TAG_BeamColumn2DwLHNMYS 197
#define ELE_TAG_BeamColumn3DwLHNMYS 198
#define ELE_TAG_PFEMLink                  199
#define ELE_TAG_PFEMContact2D             200
#define ELE_TAG_PML3D                     201
#define ELE_TAG_PML2D                     202
#define ELE_TAG_ASDShellQ4                203  // Massimo Petracca (ASDEA)
#define ELE_TAG_ASDShellT3                204  // Massimo Petracca (ASDEA)
#define ELE_TAG_WheelRail                 205
#define ELE_TAG_DispBeamColumn3dID        206 // Jose Abell the Chileno added 
#define ELE_TAG_NineNodeQuad              207
#define ELE_TAG_EightNodeQuad             208
#define ELE_TAG_SixNodeTri                209
#define ELE_TAG_RockingBC	          210
#define ELE_TAG_BeamColumn2DwLHNMYS_Damage 211
#define ELE_TAG_MVLEM_3D	          212 // Kristijan Kolozvari
#define ELE_TAG_SFI_MVLEM_3D	          213 // Kristijan Kolozvari
#define ELE_TAG_BeamGT                    214
#define ELE_TAG_MasonPan12                    215
#define ELE_TAG_MasonPan3D                    216
#define ELE_TAG_ASDEmbeddedNodeElement             217  // Massimo Petracca (ASDEA)
#define ELE_TAG_InertiaTruss              218	//Added by Xiaodong Ji, Yuhao Cheng, Yue Yu
#define ELE_TAG_ASDAbsorbingBoundary2D    219  // Massimo Petracca (ASDEA)
#define ELE_TAG_ASDAbsorbingBoundary3D    220  // Massimo Petracca (ASDEA)
#define ELE_TAG_ZeroLengthContactASDimplex  221  // Onur Deniz Akan (IUSS), Massimo Petracca (ASDEA)
#define ELE_TAG_IGALinePatch       	  250 // IGA Shell by Felipe Elgueta and jaabell (UANDES)
#define ELE_TAG_IGASurfacePatch       	  251 // IGA Shell by Felipe Elgueta and jaabell (UANDES)
#define ELE_TAG_IGAVolumePatch       	  252 // IGA Shell by Felipe Elgueta and jaabell (UANDES)
#define ELE_TAG_IGAKLShell       	  253 // IGA Shell by Felipe Elgueta and jaabell (UANDES)
#define ELE_TAG_IGAKLShell_BendingStrip   254 // IGA Shell by Felipe Elgueta and jaabell (UANDES) 216 because 208 was taken
#define ELE_TAG_PFEMContact3D             255
#define ELE_TAG_TenNodeTetrahedron        256 //by jaabell and j0selarenas (UANDES)
#define ELE_TAG_E_SFI        			257 // C. N. Lopez
#define ELE_TAG_TripleFrictionPendulumX               258
#define ELE_TAG_E_SFI_MVLEM_3D	          259259 // Kristijan Kolozvari
#define ELE_TAG_ExternalElement           99990
#define ELE_TAG_PML2D_3                   259
#define ELE_TAG_PML2D_5                   260
#define ELE_TAG_PML2D_12                  261
#define ELE_TAG_PML2DVISCOUS              262
#define ELE_TAG_MEFI        			  270 // C. N. Lopez


#define FRN_TAG_Coulomb            1
#define FRN_TAG_VelDependent       2
#define FRN_TAG_VelPressureDep     3
#define FRN_TAG_VelDepMultiLinear  4
#define FRN_TAG_VelNormalFrcDep    5

#define BEAM_INTEGRATION_TAG_Lobatto         1
#define BEAM_INTEGRATION_TAG_Legendre        2
#define BEAM_INTEGRATION_TAG_Radau           3
#define BEAM_INTEGRATION_TAG_NewtonCotes           4
#define BEAM_INTEGRATION_TAG_Trapezoidal           5
#define BEAM_INTEGRATION_TAG_CompositeSimpson           55
#define BEAM_INTEGRATION_TAG_Simpson           56
#define BEAM_INTEGRATION_TAG_Midpoint           6
#define BEAM_INTEGRATION_TAG_UserDefined     7
#define BEAM_INTEGRATION_TAG_FixedLocation     8
#define BEAM_INTEGRATION_TAG_LowOrder     9
#define BEAM_INTEGRATION_TAG_MidDistance     40
#define BEAM_INTEGRATION_TAG_Chebyshev     401

#define BEAM_INTEGRATION_TAG_ConcentratedPlasticity     41
#define BEAM_INTEGRATION_TAG_ConcentratedCurvature     42

#define BEAM_INTEGRATION_TAG_HingeMidpoint 10
#define BEAM_INTEGRATION_TAG_HingeEndpoint 11
#define BEAM_INTEGRATION_TAG_HingeRadau    12
#define BEAM_INTEGRATION_TAG_HingeRadauTwo    13
#define BEAM_INTEGRATION_TAG_UserHinge     14
#define BEAM_INTEGRATION_TAG_DistHinge     15
#define BEAM_INTEGRATION_TAG_RegularizedHinge     16

#define BEAM_INTEGRATION_TAG_HingeMidpoint2d 20
#define BEAM_INTEGRATION_TAG_HingeEndpoint2d 21
#define BEAM_INTEGRATION_TAG_HingeRadau2d    22
#define BEAM_INTEGRATION_TAG_HingeRadauTwo2d    23
#define BEAM_INTEGRATION_TAG_UserHinge2d     24
#define BEAM_INTEGRATION_TAG_DistHinge2d     25

#define BEAM_INTEGRATION_TAG_HingeMidpoint3d 30
#define BEAM_INTEGRATION_TAG_HingeEndpoint3d 31
#define BEAM_INTEGRATION_TAG_HingeRadau3d    32
#define BEAM_INTEGRATION_TAG_HingeRadauTwo3d    33
#define BEAM_INTEGRATION_TAG_UserHinge3d     34
#define BEAM_INTEGRATION_TAG_DistHinge3d     35


#define CRDTR_TAG_LinearCrdTransf2d 1
#define CRDTR_TAG_PDeltaCrdTransf2d 2
#define CRDTR_TAG_ModerateDispCrdTransf2d 8
#define CRDTR_TAG_CorotCrdTransf2d  3
#define CRDTR_TAG_CorotCrdTransfWarping2d 31
#define CRDTR_TAG_LinearCrdTransf3d 4
#define CRDTR_TAG_PDeltaCrdTransf3d 5
#define CRDTR_TAG_ModerateDispCrdTransf3d 9
#define CRDTR_TAG_CorotCrdTransf3d  6
#define CRDTR_TAG_CorotCrdTransfWarping3d  61
#define CRDTR_TAG_LinearCrdTransf2dInt 7

#define DMP_TAG_UniformDamping 1
#define DMP_TAG_SecStifDamping 2
#define DMP_TAG_URDDamping 3
#define DMP_TAG_URDDampingbeta 4

#define NOD_TAG_Node      	1
#define NOD_TAG_DummyNode 	2

#define LOAD_TAG_LoadCase  	0
#define LOAD_TAG_NodalLoad 	1
#define LOAD_TAG_EarthquakeNodalLoad 	2
#define LOAD_TAG_SingleExcitation 	3
#define LOAD_TAG_RectPulseNodalLoad 	4

#define CNSTRNT_TAG_SP_Constraint 	1
#define CNSTRNT_TAG_MP_Constraint 	2
#define CNSTRNT_TAG_ImposedMotionSP	3
#define CNSTRNT_TAG_ImposedMotionSP1	4
#define CNSTRNT_TAG_MP_Joint2D          5
#define CNSTRNT_TAG_MP_SimpleJoint2D    6
#define CNSTRNT_TAG_MP_Joint3D          7
#define CNSTRNT_TAG_Pressure_Constraint 8


#define MATRIX_TAG_Matrix 	1

#define VECTOR_TAG_Vector 	1

#define ID_TAG_ID 		1

#define HANDLER_TAG_PlainHandler 			1
#define HANDLER_TAG_LagrangeConstraintHandler   	2
#define HANDLER_TAG_PenaltyConstraintHandler    	3
#define HANDLER_TAG_TransformationConstraintHandler    	4
#define HANDLER_TAG_PenaltyHandlerNoHomoSPMultipliers   5

#define NUMBERER_TAG_DOF_Numberer      	1
#define NUMBERER_TAG_PlainNumberer 	2
#define NUMBERER_TAG_ParallelNumberer 	3

#define GraphNUMBERER_TAG_RCM   		1

#define GraphNUMBERER_TAG_SimpleNumberer   	2
#define GraphNUMBERER_TAG_MyRCM   		3
#define GraphNUMBERER_TAG_Metis   		4
#define GraphNUMBERER_TAG_AMD   		5


#define AnaMODEL_TAGS_AnalysisModel 	1

#define EquiALGORITHM_TAGS_Linear 		1
#define EquiALGORITHM_TAGS_NewtonRaphson       	2
#define EquiALGORITHM_TAGS_ModifiedNewton 	3
#define EquiALGORITHM_TAGS_Broyden 		4
#define EquiALGORITHM_TAGS_BFGS 		5
#define EquiALGORITHM_TAGS_SplitNewton 		6
#define EquiALGORITHM_TAGS_KrylovNewton         7
#define EquiALGORITHM_TAGS_NewtonLineSearch     8
#define EquiALGORITHM_TAGS_PeriodicNewton       9
#define EquiALGORITHM_TAGS_SecantNewton         10
#define EquiALGORITHM_TAGS_AcceleratedNewton          11
#define EquiALGORITHM_TAGS_AcceleratedNewtonLineSearch          12
#define EquiALGORITHM_TAGS_InitialNewton          13
#define EquiALGORITHM_TAGS_ElasticAlgorithm 14
#define EquiALGORITHM_TAGS_NewtonHallM 15
#define EquiALGORITHM_TAGS_ExpressNewton 16

#define ACCELERATOR_TAGS_Krylov		1
#define ACCELERATOR_TAGS_Secant		2
#define ACCELERATOR_TAGS_Miller         3
#define ACCELERATOR_TAGS_Monitored      4
#define ACCELERATOR_TAGS_Raphson        5
#define ACCELERATOR_TAGS_Periodic       6
#define ACCELERATOR_TAGS_Difference     7

#define LINESEARCH_TAGS_InitialInterpolatedLineSearch 1
#define LINESEARCH_TAGS_BisectionLineSearch           2
#define LINESEARCH_TAGS_RegulaFalsiLineSearch         3
#define LINESEARCH_TAGS_SecantLineSearch              4


#define INTEGRATOR_TAGS_Newmark                          1
#define INTEGRATOR_TAGS_HHT                              2
#define INTEGRATOR_TAGS_HHT_TP                           3
#define INTEGRATOR_TAGS_WilsonTheta                      4
#define INTEGRATOR_TAGS_CentralDifference                5
#define INTEGRATOR_TAGS_LoadControl                      6
#define INTEGRATOR_TAGS_DisplacementControl              7
#define INTEGRATOR_TAGS_ArcLength                        8
#define INTEGRATOR_TAGS_LoadPath                         9
#define INTEGRATOR_TAGS_Newmark1                        10
#define INTEGRATOR_TAGS_HHT1                            11
#define INTEGRATOR_TAGS_MinUnbalDispNorm                12
#define INTEGRATOR_TAGS_ArcLength1                      13
#define INTEGRATOR_TAGS_StaticSensitivity               14
#define INTEGRATOR_TAGS_HSConstraint                    15
#define INTEGRATOR_TAGS_DistributedDisplacementControl  16
#define INTEGRATOR_TAGS_CentralDifferenceAlternative    17
#define INTEGRATOR_TAGS_CentralDifferenceNoDamping      18
#define INTEGRATOR_TAGS_NewmarkExplicit                 19
#define INTEGRATOR_TAGS_NewmarkHSIncrReduct             20
#define INTEGRATOR_TAGS_NewmarkHSIncrLimit              21
#define INTEGRATOR_TAGS_NewmarkHSFixedNumIter           22
#define INTEGRATOR_TAGS_HHTExplicit                     23
#define INTEGRATOR_TAGS_HHTExplicit_TP                  24
#define INTEGRATOR_TAGS_HHTGeneralized                  25
#define INTEGRATOR_TAGS_HHTGeneralized_TP               26
#define INTEGRATOR_TAGS_HHTGeneralizedExplicit          27
#define INTEGRATOR_TAGS_HHTGeneralizedExplicit_TP       28
#define INTEGRATOR_TAGS_HHTHSIncrReduct                 29
#define INTEGRATOR_TAGS_HHTHSIncrReduct_TP              30
#define INTEGRATOR_TAGS_HHTHSIncrLimit                  31
#define INTEGRATOR_TAGS_HHTHSIncrLimit_TP               32
#define INTEGRATOR_TAGS_HHTHSFixedNumIter               33
#define INTEGRATOR_TAGS_HHTHSFixedNumIter_TP            34
#define INTEGRATOR_TAGS_AlphaOS                         35
#define INTEGRATOR_TAGS_AlphaOS_TP                      36
#define INTEGRATOR_TAGS_AlphaOSGeneralized              37
#define INTEGRATOR_TAGS_AlphaOSGeneralized_TP           38
#define INTEGRATOR_TAGS_Collocation                     39
#define INTEGRATOR_TAGS_CollocationHSIncrReduct         40
#define INTEGRATOR_TAGS_CollocationHSIncrLimit          41
#define INTEGRATOR_TAGS_CollocationHSFixedNumIter       42
#define INTEGRATOR_TAGS_TRBDF2                          43
#define INTEGRATOR_TAGS_GeneralizedAlpha                44
#define INTEGRATOR_TAGS_DisplacementPath                45
#define INTEGRATOR_TAGS_FSI                             46
#define INTEGRATOR_TAGS_TRBDF3                          47
#define INTEGRATOR_TAGS_Houbolt                         48
#define INTEGRATOR_TAGS_ParkLMS3                        49
#define INTEGRATOR_TAGS_BackwardEuler                   50
#define INTEGRATOR_TAGS_EnergyConserved                 51
#define INTEGRATOR_TAGS_PFEMIntegrator                  52
#define INTEGRATOR_TAGS_KRAlphaExplicit                 53
#define INTEGRATOR_TAGS_KRAlphaExplicit_TP              54
#define INTEGRATOR_TAGS_ExplicitDifference              55
#define INTEGRATOR_TAGS_EQPath                          56
#define INTEGRATOR_TAGS_GimmeMCK       	                57
#define INTEGRATOR_TAGS_StagedLoadControl               58
#define INTEGRATOR_TAGS_StagedNewmark                   59
#define INTEGRATOR_TAGS_HarmonicSteadyState             60


#define LinSOE_TAGS_FullGenLinSOE		1
#define LinSOE_TAGS_BandGenLinSOE		2
#define LinSOE_TAGS_BandSPDLinSOE		3
#define LinSOE_TAGS_ProfileSPDLinSOE		4
#define LinSOE_TAGS_SlowLinearSOE		5
#define LinSOE_TAGS_SparseGenColLinSOE		6
#define LinSOE_TAGS_PetscSOE       		7
#define LinSOE_TAGS_ShadowPetscSOE		8
#define LinSOE_TAGS_ActorPetscSOE		9
#define LinSOE_TAGS_UmfpackGenLinSOE		10
#define LinSOE_TAGS_SymSparseLinSOE         11
#define LinSOE_TAGS_DiagonalLinSOE         12
#define LinSOE_TAGS_ItpackLinSOE           13
#define LinSOE_TAGS_ProfileSPDLinSOEGather	14
#define LinSOE_TAGS_DistributedBandGenLinSOE		15
#define LinSOE_TAGS_DistributedBandSPDLinSOE		16
#define LinSOE_TAGS_DistributedProfileSPDLinSOE		17
#define LinSOE_TAGS_DistributedSparseGenColLinSOE       18
#define LinSOE_TAGS_DiagonalSOE       19
#define LinSOE_TAGS_SparseGenRowLinSOE		20
#define LinSOE_TAGS_DistributedSparseGenRowLinSOE       21
#define LinSOE_TAGS_DistributedDiagonalSOE 22
#define LinSOE_TAGS_MumpsSOE 23
#define LinSOE_TAGS_MumpsParallelSOE 24
#define LinSOE_TAGS_MPIDiagonalSOE 25
#define LinSOE_TAGS_PFEMLinSOE 26
#define LinSOE_TAGS_SProfileSPDLinSOE		27
#define LinSOE_TAGS_PFEMCompressibleLinSOE 28
#define LinSOE_TAGS_PFEMQuasiLinSOE 29
#define LinSOE_TAGS_PFEMDiaLinSOE 30
#define LinSOE_TAGS_PARDISOGenLinSOE 99990


#define SOLVER_TAGS_FullGenLinLapackSolver  	1
#define SOLVER_TAGS_BandGenLinLapackSolver  	2
#define SOLVER_TAGS_BandSPDLinLapackSolver  	3
#define SOLVER_TAGS_ProfileSPDLinDirectSolver  	4
#define SOLVER_TAGS_ProfileSPDLinSubstrSolver  	5
#define SOLVER_TAGS_SlowLinearSOESolver  	6
#define SOLVER_TAGS_BandSPDLinThreadSolver  	7
#define SOLVER_TAGS_ProfileSPDLinDirectThreadSolver  	8
#define SOLVER_TAGS_ProfileSPDLinDirectBlockSolver  	9
#define SOLVER_TAGS_ProfileSPDLinDirectSkypackSolver  	10
#define SOLVER_TAGS_SuperLU			      	11
#define SOLVER_TAGS_ThreadedSuperLU		      	12
#define SOLVER_TAGS_PetscSolver      			13
#define SOLVER_TAGS_UmfpackGenLinSolver      		14
#define SOLVER_TAGS_SymSparseLinSolver 15
#define SOLVER_TAGS_DiagonalLinSolver 16
#define SOLVER_TAGS_Itpack            17
#define SOLVER_TAGS_ProfileSPDLinSolverGather  	18
#define SOLVER_TAGS_DistributedSuperLU		      	19
#define SOLVER_TAGS_DiagonalDirectSolver 20
#define SOLVER_TAGS_PetscSparseSeqSolver 21
#define SOLVER_TAGS_DistributedDiagonalSolver 22
#define SOLVER_TAGS_MumpsSolver			      	23
#define SOLVER_TAGS_MumpsParallelSolver			24
#define SOLVER_TAGS_MPIDiagonalSolver                   25
#define SOLVER_TAGS_PFEMSolver                          26
#define SOLVER_TAGS_SProfileSPDLinSolver  	        27
#define SOLVER_TAGS_PFEMCompressibleSolver              28
#define SOLVER_TAGS_CulaSparseS4                        29
#define SOLVER_TAGS_CulaSparseS5                        30
#define SOLVER_TAGS_CuSP                                31
#define SOLVER_TAGS_PFEMQuasiSolver                     32
#define SOLVER_TAGS_PFEMDiaSolver                       33

#define RECORDER_TAGS_ElementRecorder		1
#define RECORDER_TAGS_NodeRecorder		2
#define RECORDER_TAGS_EnvelopeNodeRecorder	3
#define RECORDER_TAGS_EnvelopeElementRecorder	4
#define RECORDER_TAGS_DatastoreRecorder		5
#define RECORDER_TAGS_MaxNodeDispRecorder	6
#define RECORDER_TAGS_FilePlotter		7
#define RECORDER_TAGS_AlgorithmIncrements	8
#define RECORDER_TAGS_DriftRecorder		9
#define RECORDER_TAGS_EnvelopeDriftRecorder	15
#define RECORDER_TAGS_GSA_Recorder		10
#define RECORDER_TAGS_YsVisual                  11
#define RECORDER_TAGS_DamageRecorder		12
#define RECORDER_TAGS_PatternRecorder		13
#define RECORDER_TAGS_TclFeViewer		14
#define RECORDER_TAGS_NormElementRecorder	16
#define RECORDER_TAGS_NormNodeRecorder	        17
#define RECORDER_TAGS_NormEnvelopeElementRecorder	18
#define RECORDER_TAGS_PVDRecorder               19
#define RECORDER_TAGS_MPCORecorder               20
#define RECORDER_TAGS_GmshRecorder               21
#define RECORDER_TAGS_VTK_Recorder               22
#define RECORDER_TAGS_NodeRecorderRMS               23
#define RECORDER_TAGS_ElementRecorderRMS               24

#define OPS_STREAM_TAGS_FileStream		1
#define OPS_STREAM_TAGS_StandardStream		2
#define OPS_STREAM_TAGS_XmlFileStream		3
#define OPS_STREAM_TAGS_DataFileStream		4
#define OPS_STREAM_TAGS_DatabaseStream		5
#define OPS_STREAM_TAGS_DummyStream		6
#define OPS_STREAM_TAGS_BinaryFileStream        7
#define OPS_STREAM_TAGS_TCP_Stream              8
#define OPS_STREAM_TAGS_ChannelStream           9
#define OPS_STREAM_TAGS_DataTurbineStream      10
#define OPS_STREAM_TAGS_DataFileStreamAdd      11


#define DomDecompALGORITHM_TAGS_DomainDecompAlgo 1

#define DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis 1
#define ANALYSIS_TAGS_StaticDomainDecompositionAnalysis 2
#define ANALYSIS_TAGS_TransientDomainDecompositionAnalysis 3

#define PartitionedModelBuilder_TAGS_PartitionedQuick2dFrameModel 1

#define RANDOM_VARIABLE_beta				1
#define RANDOM_VARIABLE_chisquare			2
#define RANDOM_VARIABLE_exponential			3
#define RANDOM_VARIABLE_gamma				4
#define RANDOM_VARIABLE_gumbel				5
#define RANDOM_VARIABLE_laplace				6
#define RANDOM_VARIABLE_lognormal			7
#define RANDOM_VARIABLE_normal				8
#define RANDOM_VARIABLE_pareto				9
#define RANDOM_VARIABLE_rayleigh			10
#define RANDOM_VARIABLE_shiftedexponential	11
#define RANDOM_VARIABLE_shiftedrayleigh		12
#define RANDOM_VARIABLE_type1largestvalue	13
#define RANDOM_VARIABLE_type1smallestvalue	14
#define RANDOM_VARIABLE_type2largestvalue	15
#define RANDOM_VARIABLE_type3smallestvalue	16
#define RANDOM_VARIABLE_uniform				17
#define RANDOM_VARIABLE_weibull				18
#define RANDOM_VARIABLE_userdefined             19
#define RANDOM_VARIABLE_python             20

#define RANDOM_VARIABLE_POSITIONER        1
#define PARAMETER_POSITIONER              2

#define CORRELATION_COEFFICIENT           1

#define LIMIT_STATE_FUNCTION		  1
#define LIMCRV_TAG_WrapperLimitCurve      1

#define CUTSET			  1

#define MODULATING_FUNCTION_gamma         1
#define MODULATING_FUNCTION_constant      2
#define MODULATING_FUNCTION_trapezoidal   3

#define FILTER_standardLinearOscillator   1

#define SPECTRUM_jonswap                  1
#define SPECTRUM_constant                 2
#define SPECTRUM_points                   3


#define CHANNEL_TAGS_FileDatastore	  1

#endif
