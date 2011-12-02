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
                                                                        
// $Revision: 1.4 $
// $Date: 2000-12-13 08:30:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/classTags.h,v $
                                                                        
                                                                        
// File: ~classTags.h
//
// Written: fmk
// Created: Fri Sept 20 12:27:47: 1996
// Revision: A
//
// Purpose: This file contains the declarations of all class tags used.
//
// What: "@(#) classTags.h, revA"

#ifndef classTags_h
#define classTags_h

#define intType 1
#define doubleType 2
#define idType 3
#define vectorType 4
#define matrixType 5

#define EigenSOE_TAGS_BandArpackSOE 	1
#define EigenSOE_TAGS_SymArpackSOE 	2
#define EigenSOLVER_TAGS_BandArpackSolver 	1
#define EigenSOLVER_TAGS_SymArpackSolver 	2

#define EigenALGORITHM_TAGS_Frequency 1

#define EigenINTEGRATOR_TAGS_Eigen 1

#define CONVERGENCE_TEST_CTestNormUnbalance      1
#define CONVERGENCE_TEST_CTestNormDispIncr       2
#define CONVERGENCE_TEST_CTestEnergyIncr         3

#define GRND_TAG_ElCentroGroundMotion            1
#define GROUND_MOTION_TAG_GroundMotionRecord     2
#define GROUND_MOTION_TAG_InterpolatedGroundMotion    3
#define GROUND_MOTION_TAG_GroundMotion     4


#define TIMESERIES_INTEGRATOR_TAG_Trapezoidal 1

#define SECT_TAG_Section         1

#define TSERIES_TAG_LinearSeries         1
#define TSERIES_TAG_RectangularSeries          2
#define TSERIES_TAG_PathTimeSeries       3
#define TSERIES_TAG_PathSeries       4
#define TSERIES_TAG_ConstantSeries       5
#define TSERIES_TAG_TrigSeries       6

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
#define MAT_TAG_Clough1	202
#define MAT_TAG_Clough2	202
#define MAT_TAG_Pinch1	203
#define MAT_TAG_BiLinear	204


#define SEC_TAG_Elastic2d   3
#define SEC_TAG_Elastic3d   4
#define SEC_TAG_Generic1d	5
#define SEC_TAG_GenericNd	6
#define SEC_TAG_Aggregator	7
#define SEC_TAG_Fiber		8

#define ND_TAG_ElasticIsotropic					10
#define ND_TAG_ElasticIsotropicPlaneStrain2d	11
#define ND_TAG_ElasticIsotropicPlaneStress2d	12
#define ND_TAG_ElasticIsotropic3D               21
#define ND_TAG_Template3Dep 			31
#define ND_TAG_J2PlaneStrain                  3005 
#define ND_TAG_J2PlaneStress                  3006 
#define ND_TAG_J2AxiSymm                      3007 
#define ND_TAG_J2ThreeDimensional             3009 
#define MAT_TAG_FluidSolidPorousMaterial       102

#define FIBER_TAG_Uniaxial2d	1
#define FIBER_TAG_Uniaxial3d	2

#define BACKBONE_TAG_Capped			1
#define BACKBONE_TAG_LinearCapped	2
#define BACKBONE_TAG_Material		3
#define BACKBONE_TAG_Petrangeli		4
#define BACKBONE_TAG_Trilinear		5

#define DEG_TAG_STIFF_Constant		1
#define DEG_TAG_STIFF_Ductility		2
#define DEG_TAG_STIFF_Stanford		3

#define DEG_TAG_DEF_Constant		1
#define DEG_TAG_DEF_Ductility		2
#define DEG_TAG_DEF_Stanford		3

#define DEG_TAG_STRENGTH_ACI		1
#define DEG_TAG_STRENGTH_Constant	2
#define DEG_TAG_STRENGTH_Petrangeli	3
#define DEG_TAG_STRENGTH_Stanford	4

#define PATTERN_TAG_LoadPattern		  1
#define PATTERN_TAG_MultiSupportPattern	  3
#define LOAD_TAG_UniformExcitation        2

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
#define ELE_TAG_fElmt05	           5
#define ELE_TAG_fElmt02	           2
// #define ELE_TAG_MyTruss    	 4002
#define ELE_TAG_ZeroLength	 5000
#define ELE_TAG_ZeroLengthSection	 5001
#define ELE_TAG_ZeroLengthND	 5002
#define ELE_TAG_NLBeamColumn2d	 6000
#define ELE_TAG_NLBeamColumn3d	 6001
#define ELE_TAG_FourNodeQuad	 1010
#define ELE_TAG_BeamWithHinges2d  401  
#define ELE_TAG_BeamWithHinges3d  402

#define CRDTR_TAG_LinearCrdTransf2d 1
#define CRDTR_TAG_CorotCrdTransf2d  2
#define CRDTR_TAG_LinearCrdTransf3d 3
#define CRDTR_TAG_CorotCrdTransf3d  4

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

#define MATRIX_TAG_Matrix 	1

#define VECTOR_TAG_Vector 	1

#define ID_TAG_ID 		1

#define HANDLER_TAG_PlainHandler 			1
#define HANDLER_TAG_LagrangeConstraintHandler   	2
#define HANDLER_TAG_PenaltyConstraintHandler    	3
#define HANDLER_TAG_TransformationConstraintHandler    	4

#define NUMBERER_TAG_DOF_Numberer      	1
#define NUMBERER_TAG_PlainNumberer 	2

#define GraphNUMBERER_TAG_RCM   		1
#define GraphNUMBERER_TAG_SimpleNumberer   	2
#define GraphNUMBERER_TAG_MyRCM   		3
#define GraphNUMBERER_TAG_Metis   		4


#define AnaMODEL_TAGS_AnalysisModel 	1

#define EquiALGORITHM_TAGS_Linear 		1
#define EquiALGORITHM_TAGS_NewtonRaphson       	2
#define EquiALGORITHM_TAGS_ModifiedNewton 	3
#define EquiALGORITHM_TAGS_BFGS 		4
#define EquiALGORITHM_TAGS_SplitNewton 		5

#define INTEGRATOR_TAGS_Newmark          	2
#define INTEGRATOR_TAGS_HHT 	          	3
#define INTEGRATOR_TAGS_WilsonTheta          	4
#define INTEGRATOR_TAGS_CentralDifference       5
#define INTEGRATOR_TAGS_LoadControl          	6
#define INTEGRATOR_TAGS_DisplacementControl     7
#define INTEGRATOR_TAGS_ArcLength	     	5
#define INTEGRATOR_TAGS_LoadPath          	8
#define INTEGRATOR_TAGS_Newmark1          	9
#define INTEGRATOR_TAGS_HHT1 	          	10
#define INTEGRATOR_TAGS_MinUnbalDispNorm 	11
#define INTEGRATOR_TAGS_ArcLength1	     	12

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

#define DomDecompALGORITHM_TAGS_DomainDecompAlgo 1

#define DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis 1

#define PartitionedModelBuilder_TAGS_PartitionedQuick2dFrameModel 1

#endif


