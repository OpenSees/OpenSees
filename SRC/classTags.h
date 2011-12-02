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
                                                                        
// $Revision: 1.73 $
// $Date: 2006-01-17 20:48:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/classTags.h,v $
                                                                        
// Written: fmk
// Revision: A
//
// Purpose: This file contains the declarations of all class tags used.
//
// What: "@(#) classTags.h, revA"

#ifndef classTags_h
#define classTags_h

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

#define EigenSOLVER_TAGS_BandArpackSolver 	1
#define EigenSOLVER_TAGS_SymArpackSolver 	2
#define EigenSOLVER_TAGS_SymBandEigenSolver     3

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


#define GRND_TAG_ElCentroGroundMotion                 1
#define GROUND_MOTION_TAG_GroundMotionRecord          2
#define GROUND_MOTION_TAG_InterpolatedGroundMotion    3
#define GROUND_MOTION_TAG_GroundMotion                4

#define REGION_TAG_MeshRegion      1

#define TIMESERIES_INTEGRATOR_TAG_Trapezoidal 1

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


#define MAT_TAG_ElasticMaterial			1
#define MAT_TAG_ElasticPPMaterial		2
#define MAT_TAG_ParallelMaterial		3
#define MAT_TAG_Concrete01				4
#define MAT_TAG_Steel01					5

#define MAT_TAG_Hardening				6
#define MAT_TAG_Hysteretic				7
#define MAT_TAG_EPPGap					8
#define MAT_TAG_Viscous					9
#define MAT_TAG_Backbone				10
#define MAT_TAG_PathIndependent			11
#define MAT_TAG_SeriesMaterial			12
#define MAT_TAG_CableMaterial          13
#define MAT_TAG_ENTMaterial				14
#define MAT_TAG_Penalty					15
#define MAT_TAG_MinMax					16
#define MAT_TAG_BoucWen					17
#define MAT_TAG_Pinching4					18
#define MAT_TAG_BarSlip					19
#define MAT_TAG_Fatigue			20
#define MAT_TAG_Steel03					23
#define MAT_TAG_ReinforcingSteel					24

// Terje:
#define MAT_TAG_SmoothSteel01			21
#define MAT_TAG_SmoothConcrete01		22


#define MAT_TAG_Concrete04 30

//B Jeremic
#define MAT_TAG_PySimple1        205
#define MAT_TAG_TzSimple1        206
#define MAT_TAG_QzSimple1        207
#define MAT_TAG_PyLiq1           208
#define MAT_TAG_TzLiq1           209

#define MAT_TAG_FedeasMaterial         1000
#define MAT_TAG_FedeasBond1       1001
#define MAT_TAG_FedeasBond2       1002
#define MAT_TAG_FedeasConcrete1       1003
#define MAT_TAG_FedeasConcrete2       1004
#define MAT_TAG_FedeasConcrete3       1005
#define MAT_TAG_FedeasHardening       1006
#define MAT_TAG_FedeasHysteretic1       1007
#define MAT_TAG_FedeasHysteretic2       1008
#define MAT_TAG_FedeasSteel1       1009
#define MAT_TAG_FedeasSteel2       1010

#define MAT_TAG_DrainMaterial		2000
#define MAT_TAG_DrainHardening		2001
#define MAT_TAG_DrainBilinear		2002
#define MAT_TAG_DrainClough1		2003
#define MAT_TAG_DrainClough2		2004
#define MAT_TAG_DrainPinch1			2005
#define MAT_TAG_DrainPinch2			2006

#define MAT_TAG_SnapMaterial		3000
#define MAT_TAG_SnapBilinear		3001
#define MAT_TAG_SnapClough		3002
#define MAT_TAG_SnapPinch		3003
#define MAT_TAG_SnapCloughDamage	3004
#define MAT_TAG_SnapPinchingDamage	3005

#define MAT_TAG_Clough1	201
#define MAT_TAG_Clough2	202
#define MAT_TAG_Pinch1	203
#define MAT_TAG_BiLinear	204
#define MAT_TAG_Pinching	205


#define SEC_TAG_Elastic2d   3
#define SEC_TAG_Elastic3d   4
#define SEC_TAG_Generic1d	5
#define SEC_TAG_GenericNd	6
#define SEC_TAG_Aggregator	7
#define SEC_TAG_Fiber		8
#define SEC_TAG_FiberSection2d		9
#define SEC_TAG_FiberSection3d		10
#define SEC_TAG_FiberSectionGJ		11
#define SEC_TAG_BeamFiberSection	12
#define SEC_TAG_ElasticPlateSection	13
#define SEC_TAG_ElasticMembranePlateSection	14
#define SEC_TAG_MembranePlateFiberSection	15 
#define SEC_TAG_Bidirectional	16
#define SEC_TAG_Isolator2spring	17

#define ND_TAG_ElasticIsotropic					10
#define ND_TAG_ElasticIsotropicPlaneStrain2d	11
#define ND_TAG_ElasticIsotropicPlaneStress2d	12
#define ND_TAG_ElasticIsotropicAxiSymm          13
#define ND_TAG_ElasticIsotropicPlateFiber		14
#define ND_TAG_ElasticIsotropicBeamFiber		15 
#define ND_TAG_ElasticIsotropic3D               16
#define ND_TAG_ElasticCrossAnisotropic3D        17
#define ND_TAG_J2PlaneStrain                  3005 
#define ND_TAG_J2PlaneStress                  3006 
#define ND_TAG_J2AxiSymm                      3007 
#define ND_TAG_J2ThreeDimensional             3009 
#define ND_TAG_J2PlateFiber					3010
#define ND_TAG_J2BeamFiber					3011
#define ND_TAG_PressureDependentElastic3D       22
#define ND_TAG_Template3Dep 			        31
#define ND_TAG_FluidSolidPorousMaterial        100
#define ND_TAG_PressureDependMultiYield		101
#define ND_TAG_PressureIndependMultiYield		102
#define ND_TAG_PressureDependMultiYield02		103
#define ND_TAG_FeapMaterial                 1000
#define ND_TAG_FeapMaterial01                 1001
#define ND_TAG_FeapMaterial02                 1002
#define ND_TAG_FeapMaterial03                 1003
#define ND_TAG_PlaneStressMaterial          2000
#define ND_TAG_PlateFiberMaterial          2001
#define ND_TAG_BeamFiberMaterial		2002
#define ND_TAG_CompressibleFluid		3001
#define ND_TAG_GeneralizedPlasticity 3002
#define ND_TAG_J2Plasticity02  3003
#define ND_TAG_FiniteDeformationElastic3D	8002
#define ND_TAG_NeoHookeanCompressible3D	        8003
#define ND_TAG_FDdecoupledElastic3D	        8004
#define ND_TAG_FiniteDeformationEP3D	        8005

// MultiaxialCyclicPlasticity, add by Gang Wang
#define ND_TAG_MultiaxialCyclicPlasticity             10031
#define ND_TAG_MultiaxialCyclicPlasticity3D           10032 
#define ND_TAG_MultiaxialCyclicPlasticityAxiSymm      10033
#define ND_TAG_MultiaxialCyclicPlasticityPlaneStrain  10034

#define FIBER_TAG_Uniaxial2d	1
#define FIBER_TAG_Uniaxial3d	2

#define BACKBONE_TAG_Capped		1
#define BACKBONE_TAG_LinearCapped	2
#define BACKBONE_TAG_Material		3
#define BACKBONE_TAG_Arctangent		4
#define BACKBONE_TAG_Trilinear		5
#define BACKBONE_TAG_Multilinear	6
#define BACKBONE_TAG_Mander		7
#define BACKBONE_TAG_KentPark		8

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

#define PATTERN_TAG_LoadPattern		  1
#define PATTERN_TAG_MultiSupportPattern	  3
#define PATTERN_TAG_UniformExcitation     2
#define LOAD_TAG_Beam2dUniformLoad        3
#define LOAD_TAG_Beam2dPointLoad          4
#define LOAD_TAG_Beam3dUniformLoad        5
#define LOAD_TAG_Beam3dPointLoad          6
#define LOAD_TAG_BrickSelfWeight          7
#define LOAD_TAG_Beam2dTempLoad           8
#define PATTERN_TAG_PBowlLoading	  10


#define MAT_TAG_IsotropicLinElastic         1001
#define MAT_TAG_IsotropicLinElasticPoint    1002
#define MAT_TAG_OrthotropicLinElastic       1003
#define MAT_TAG_OrthotropicLinElasticPoint  1004

#define ELE_TAG_cont2d01    	2101	// provisional
#define ELE_TAG_cont2d02    	2102	// provisional
#define ELE_TAG_CST	    	4050

#define ELE_TAG_Subdomain     	1
#define ELE_TAG_ElasticBeam2d   2000
#define ELE_TAG_ElasticBeam3d   3000
#define ELE_TAG_Beam2d    	2001
#define ELE_TAG_beam2d02    	2002
#define ELE_TAG_beam2d03    	2003
#define ELE_TAG_beam2d04    	2004
#define ELE_TAG_beam3d01    	3001
#define ELE_TAG_beam3d02    	3002
#define ELE_TAG_Truss    	4001
#define ELE_TAG_TrussSection    4005
#define ELE_TAG_CorotTruss    	4003
#define ELE_TAG_CorotTrussSection    	4004
#define ELE_TAG_fElmt05	           5
#define ELE_TAG_fElmt02	           2
#define ELE_TAG_MyTruss    	 4002
#define ELE_TAG_ZeroLength	 5000
#define ELE_TAG_ZeroLengthSection	 5001
#define ELE_TAG_ZeroLengthND	 5002
#define ELE_TAG_NLBeamColumn2d	 6000
#define ELE_TAG_NLBeamColumn3d	 6001
#define ELE_TAG_LargeDispBeamColumn3d	 6002
#define ELE_TAG_FourNodeQuad	 1010
#define ELE_TAG_BeamWithHinges2d  401  
#define ELE_TAG_BeamWithHinges3d  402
#define ELE_TAG_EightNodeBrick   7001
#define ELE_TAG_TwentyNodeBrick   7002
#define ELE_TAG_EightNodeBrick_u_p_U            7003
#define ELE_TAG_TwentyNodeBrick_u_p_U            7004
#define ELE_TAG_FourNodeQuadUP  7005
#define ELE_TAG_TotalLagrangianFD20NodeBrick 7006 // ZC added
#define ELE_TAG_TwentySevenNodeBrick   7010
#define ELE_TAG_BrickUP         7011
#define ELE_TAG_Nine_Four_Node_QuadUP         7012
#define ELE_TAG_Twenty_Eight_Node_BrickUP         7013
#define ELE_TAG_Twenty_Node_Brick         7014
#define ELE_TAG_PlateMITC4      2023 
#define ELE_TAG_ShellMITC4      2024 
#define ELE_TAG_Plate1          2022 
#define ELE_TAG_Brick                      3458 
#define ELE_TAG_BbarBrick                  3457 
#define ELE_TAG_FLBrick                  3455 
#define ELE_TAG_EnhancedQuad               3459
#define ELE_TAG_ConstantPressureVolumeQuad 3456 
#define ELE_TAG_NineNodeMixedQuad          3359 
#define ELE_TAG_DispBeamColumn2d 9870
#define ELE_TAG_DispBeamColumn3d 9871
#define ELE_TAG_HingedBeam2d     9872
#define ELE_TAG_HingedBeam3d     9873
#define ELE_TAG_TwoPointHingedBeam2d     9874
#define ELE_TAG_TwoPointHingedBeam3d     9875
#define ELE_TAG_OnePointHingedBeam2d     9876
#define ELE_TAG_OnePointHingedBeam3d     9877
#define ELE_TAG_BeamColumnJoint2d 4444
#define ELE_TAG_BeamColumnJoint3d 4445

#define ELE_TAG_ForceBeamColumn2d 9878
#define ELE_TAG_ForceBeamColumn3d 9879

#define ELE_TAG_InternalSpring   9900
#define ELE_TAG_SimpleJoint2D    9901
#define ELE_TAG_Joint2D    9902
#define ELE_TAG_Joint3D    9903


// Terje:
#define ELE_TAG_FourNodeQuad02  9903

#define BEAM_INTEGRATION_TAG_Lobatto         1
#define BEAM_INTEGRATION_TAG_Legendre        2
#define BEAM_INTEGRATION_TAG_Radau           3
#define BEAM_INTEGRATION_TAG_NewtonCotes           4
#define BEAM_INTEGRATION_TAG_Trapezoidal           5
#define BEAM_INTEGRATION_TAG_Midpoint           6
#define BEAM_INTEGRATION_TAG_UserDefined     7
#define BEAM_INTEGRATION_TAG_FixedLocation     8

#define BEAM_INTEGRATION_TAG_HingeMidpoint 10
#define BEAM_INTEGRATION_TAG_HingeEndpoint 11
#define BEAM_INTEGRATION_TAG_HingeRadau    12
#define BEAM_INTEGRATION_TAG_HingeRadauTwo    13
#define BEAM_INTEGRATION_TAG_UserHinge     14
#define BEAM_INTEGRATION_TAG_DistHinge     15

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
#define CRDTR_TAG_CorotCrdTransf2d  3
#define CRDTR_TAG_LinearCrdTransf3d 4
#define CRDTR_TAG_PDeltaCrdTransf3d 5
#define CRDTR_TAG_CorotCrdTransf3d  6

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
#define EquiALGORITHM_TAGS_AccelNewton          11

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


#define INTEGRATOR_TAGS_Newmark          	             1
#define INTEGRATOR_TAGS_HHT 	          	             2
#define INTEGRATOR_TAGS_WilsonTheta                 	 3
#define INTEGRATOR_TAGS_CentralDifference                4
#define INTEGRATOR_TAGS_LoadControl          	         5
#define INTEGRATOR_TAGS_DisplacementControl              6
#define INTEGRATOR_TAGS_ArcLength	     	             7
#define INTEGRATOR_TAGS_LoadPath          	             8
#define INTEGRATOR_TAGS_Newmark1          	             9
#define INTEGRATOR_TAGS_HHT1 	          	            10
#define INTEGRATOR_TAGS_MinUnbalDispNorm 	            11
#define INTEGRATOR_TAGS_ArcLength1	     	            12
#define INTEGRATOR_TAGS_StaticSensitivity				13
#define INTEGRATOR_TAGS_HSConstraint					14
#define INTEGRATOR_TAGS_DistributedDisplacementControl  15
#define INTEGRATOR_TAGS_CentralDifferenceAlternative	16
#define INTEGRATOR_TAGS_CentralDifferenceNoDamping		17
#define INTEGRATOR_TAGS_NewmarkExplicit		            18
#define INTEGRATOR_TAGS_NewmarkHybridSimulation		    19
#define INTEGRATOR_TAGS_HHTExplicit 	          	    20
#define INTEGRATOR_TAGS_HHTGeneralized 	          	    21
#define INTEGRATOR_TAGS_HHTGeneralizedExplicit          22
#define INTEGRATOR_TAGS_HHTHybridSimulation 	        23
#define INTEGRATOR_TAGS_AlphaOS                         24
#define INTEGRATOR_TAGS_AlphaOSGeneralized              25
#define INTEGRATOR_TAGS_Collocation 	          	    26
#define INTEGRATOR_TAGS_CollocationHybridSimulation 	27



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


#define RECORDER_TAGS_ElementRecorder		1
#define RECORDER_TAGS_NodeRecorder		2
#define RECORDER_TAGS_EnvelopeNodeRecorder	3
#define RECORDER_TAGS_EnvelopeElementRecorder	4
#define RECORDER_TAGS_DatastoreRecorder		5
#define RECORDER_TAGS_MaxNodeDispRecorder	6
#define RECORDER_TAGS_FilePlotter		7
#define RECORDER_TAGS_AlgorithmIncrements	8
#define RECORDER_TAGS_DriftRecorder		9
#define RECORDER_TAGS_GSA_Recorder		10
#define RECORDER_TAGS_YsVisual                  11
#define RECORDER_TAGS_DamageRecorder		12
#define RECORDER_TAGS_PatternRecorder		13
#define RECORDER_TAGS_TclFeViewer		14

#define DATAHANDLER_TAGS_DataOutputStreamHandler		1
#define DATAHANDLER_TAGS_DataOutputFileHandler		2
#define DATAHANDLER_TAGS_DataOutputDatabaseHandler		3

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

#define RANDOM_VARIABLE_POSITIONER        1
#define PARAMETER_POSITIONER              2

#define CORRELATION_COEFFICIENT           1

#define LIMIT_STATE_FUNCTION			  1

#define MODULATING_FUNCTION_gamma         1
#define MODULATING_FUNCTION_constant      2
#define MODULATING_FUNCTION_trapezoidal   3

#define FILTER_standardLinearOscillator   1

#define SPECTRUM_jonswap                  1
#define SPECTRUM_constant                 2
#define SPECTRUM_points                   3


#endif


