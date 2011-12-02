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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/objectBroker/FEM_ObjectBroker.cpp,v $
                                                                        
                                                                        
// File: ~/actor/broker/FEM_ObjectBroker.C
//
// Written: fmk
// Created: 10/96
// Revision: A
//
// Purpose: This file contains the class definition for FEM_ObjectBroker.
// FEM_ObjectBroker is is an object broker class for the finite element
// method. All methods are virtual to allow for subclasses; which can be
// used by programmers when introducing new subclasses of the main objects.
//
// What: "@(#) FEM_ObjectBroker.C, revA"

#include <FEM_ObjectBroker.h>

// Convergence tests
#include <CTestNormUnbalance.h>
#include <CTestNormDispIncr.h>
#include <CTestEnergyIncr.h>

// graph numbering schemes
#include <RCM.h>
#include <MyRCM.h>
#include <SimpleNumberer.h>

// uniaxial material model header files
#include <ElasticMaterial.h>
#include <ElasticPPMaterial.h>
#include <ParallelMaterial.h>
#include <Concrete01.h>
#include <Steel01.h>
#include <HardeningMaterial.h>
#include <HystereticMaterial.h>
#include <EPPGapMaterial.h>
#include <ViscousMaterial.h>
#include <PathIndependentMaterial.h>
#include <SeriesMaterial.h>

// Sections
#include <ElasticSection2d.h>
#include <ElasticSection3d.h>
#include <GenericSection1d.h>
#include <GenericSectionNd.h>
#include <SectionAggregator.h>
#include <FiberSection.h>

// NDMaterials
#include <ElasticIsotropicPlaneStrain2D.h>
#include <ElasticIsotropicPlaneStress2D.h>

// Fibers
#include <UniaxialFiber2d.h>
#include <UniaxialFiber3d.h>

// element header files
#include <Element.h>
#include <beam2d02.h>
#include <beam2d03.h>
#include <beam2d04.h>
#include <beam3d01.h>
#include <beam3d02.h>
#include <Truss.h>
#include <TrussSection.h>
#include <ZeroLength.h>
#include <FourNodeQuad.h>
#include <ElasticBeam2d.h>
#include <ElasticBeam3d.h>
#include <BeamWithHinges2d.h>
#include <BeamWithHinges3d.h>
#include <NLBeamColumn2d.h>
#include <NLBeamColumn3d.h>

#include <LinearCrdTransf2d.h>
#include <LinearCrdTransf3d.h>

// node header files
#include <Node.h>

// mp_constraint header files
#include <MP_Constraint.h>

// sp_constraint header files
#include <SP_Constraint.h>

// nodal load header files
#include <NodalLoad.h>

// elemental load header files
#include <ElementalLoad.h>

// matrix, vector & id header files
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

// subdomain header files
#include <Subdomain.h>

// constraint handler header files
#include <ConstraintHandler.h>
#include <PlainHandler.h>

// dof numberer header files
#include <DOF_Numberer.h>   
#include <PlainNumberer.h>

// analysis model header files
#include <AnalysisModel.h>    

// equi soln algo header files
#include <EquiSolnAlgo.h>
#include <Linear.h>
#include <NewtonRaphson.h>
#include <ModifiedNewton.h>

// domain decomp soln algo header files
#include <DomainDecompAlgo.h>

// integrator header files
#include <LoadControl.h>
#include <ArcLength.h>
#include <TransientIntegrator.h>
#include <Newmark.h>

// system of eqn header files
#include <LinearSOE.h>
#include <DomainSolver.h>
#include <FullGenLinSOE.h>
#include <FullGenLinLapackSolver.h>
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>
#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <ProfileSPDLinSubstrSolver.h>


#include <DomainDecompositionAnalysis.h>

// load patterns
#include <LoadPattern.h>
#include <GroundMotion.h>

// time series
#include <LinearSeries.h>
#include <PathSeries.h>
#include <PathTimeSeries.h>
#include <RectangularSeries.h>
#include <ConstantSeries.h>
#include <TrigSeries.h>

FEM_ObjectBroker::FEM_ObjectBroker()
:lastLinearSolver(0),lastDomainSolver(0)
{

}


FEM_ObjectBroker::~FEM_ObjectBroker()
{

}



PartitionedModelBuilder          *
FEM_ObjectBroker::getPtrNewPartitionedModelBuilder(Subdomain &theSubdomain,
						   int classTag)
{
    switch(classTag) {
	/*
	case PartitionedModelBuilder_TAGS_PartitionedQuick2dFrameModel:  
	     return new PartitionedQuick2dFrame(theSubdomain);
	     */

	default:
	     cerr << "FEM_ObjectBroker::getPtrNewPartitionedModelBuilder - ";
	     cerr << " - no PartitionedModelBuilder type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }    
}


GraphNumberer *
FEM_ObjectBroker::getPtrNewGraphNumberer(int classTag)
{
    switch(classTag) {
	case GraphNUMBERER_TAG_RCM:  
 	     return new RCM();
	     
	     
	case GraphNUMBERER_TAG_MyRCM:  
	     return new MyRCM();
	     	     
	     
	case GraphNUMBERER_TAG_SimpleNumberer:  
	     return new SimpleNumberer();				
	     
	     
	default:
	     cerr << "ObjectBroker::getPtrNewGraphNumberer - ";
	     cerr << " - no GraphNumberer type exists for class tag " ;
	     cerr << classTag << endl;
	     return 0;
	     
	 }
}

/*****************************************
 *
 * METHODS TO GET NEW MODELLING CLASSES
 *
 *****************************************/








Element       *
FEM_ObjectBroker::getNewElement(int classTag)
{
    switch(classTag) {
	     
	case ELE_TAG_Truss:  
	     return new Truss(); 
	     
	case ELE_TAG_TrussSection:  
	     return new TrussSection(); 	     
	     
	case ELE_TAG_ZeroLength:  
	     return new ZeroLength(); 	     
	     
	case ELE_TAG_FourNodeQuad:  
	     return new FourNodeQuad(); 	     
	     
	case ELE_TAG_ElasticBeam2d:
		return new ElasticBeam2d();

	case ELE_TAG_ElasticBeam3d:
		return new ElasticBeam3d();

	case ELE_TAG_BeamWithHinges2d:  
	     return new BeamWithHinges2d(); 	     	     

	case ELE_TAG_BeamWithHinges3d:  
	     return new BeamWithHinges3d(); 	     	     

	case ELE_TAG_NLBeamColumn2d:  
	     return new NLBeamColumn2d();					     

	case ELE_TAG_NLBeamColumn3d:  
	     return new NLBeamColumn3d();  
				
	default:
	     cerr << "FEM_ObjectBroker::getNewElement - ";
	     cerr << " - no Element type exists for class tag " ;
	     cerr << classTag << endl;
	     return 0;
	     
	 }
}
				
Node          *
FEM_ObjectBroker::getNewNode(int classTag)
{
    switch(classTag) {
	case NOD_TAG_Node:  
	     return new Node(classTag);
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewNode - ";
	     cerr << " - no Node type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }    
}


MP_Constraint *
FEM_ObjectBroker::getNewMP(int classTag)
{
    switch(classTag) {
	case CNSTRNT_TAG_MP_Constraint:  
	     return new MP_Constraint(classTag);
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewMP - ";
	     cerr << " - no SP_Constraint type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }    
}


SP_Constraint *
FEM_ObjectBroker::getNewSP(int classTag)
{
    switch(classTag) {
	case CNSTRNT_TAG_SP_Constraint:  
	     return new SP_Constraint(classTag);
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewSP - ";
	     cerr << " - no SP_Constraint type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }    
}

NodalLoad     *
FEM_ObjectBroker::getNewNodalLoad(int classTag)
{
    switch(classTag) {
	case LOAD_TAG_NodalLoad:  
	     return new NodalLoad(classTag);
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewNodalLoad - ";
	     cerr << " - no NodalLoad type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }    
}


ElementalLoad *
FEM_ObjectBroker::getNewElementalLoad(int classTag)
{
    cerr << "FEM_ObjectBroker: NOT IMPLEMENTED YET";
    return 0;
}

CrdTransf2d*
FEM_ObjectBroker::getNewCrdTransf2d(int classTag)
{
	switch(classTag) {
	case CRDTR_TAG_LinearCrdTransf2d:
		return new LinearCrdTransf2d();
	//case CRDTR_TAG_CorotCrdTransf2d :
	//	return new CorotCrdTransf2d();
	default:
		cerr << "FEM_ObjectBroker::getPtrNewCrdTransf2d - ";
	    cerr << " - no CrdTransf2d type exists for class tag ";
	    cerr << classTag << endl;
	    return 0;
	}

}

CrdTransf3d*
FEM_ObjectBroker::getNewCrdTransf3d(int classTag)
{
	switch(classTag) {
	case CRDTR_TAG_LinearCrdTransf3d:
		return new LinearCrdTransf3d();
	//case CRDTR_TAG_CorotCrdTransf3d :
	//	return new CorotCrdTransf3d();
	default:
		cerr << "FEM_ObjectBroker::getPtrNewCrdTransf3d - ";
	    cerr << " - no CrdTransf3d type exists for class tag ";
	    cerr << classTag << endl;
	    return 0;
	}

}

UniaxialMaterial *
FEM_ObjectBroker::getNewUniaxialMaterial(int classTag)
{
    switch(classTag) {
	case MAT_TAG_ElasticMaterial:  
	     return new ElasticMaterial(); // values set in recvSelf
	     
	case MAT_TAG_ElasticPPMaterial:  
	     return new ElasticPPMaterial(); // values set in recvSelf
	     	     
	case MAT_TAG_ParallelMaterial:  
	     return new ParallelMaterial();

	case MAT_TAG_Concrete01:  
	     return new Concrete01();

	case MAT_TAG_Steel01:  
	     return new Steel01();

	case MAT_TAG_Hardening:
		return new HardeningMaterial();

	case MAT_TAG_Hysteretic:
		return new HystereticMaterial();

	case MAT_TAG_EPPGap:
		return new EPPGapMaterial();

	case MAT_TAG_Viscous:
		return new ViscousMaterial();

	case MAT_TAG_PathIndependent:
		return new PathIndependentMaterial();

	case MAT_TAG_SeriesMaterial:
		return new SeriesMaterial();
	     
	default:
	     cerr << "FEM_ObjectBroker::getPtrNewUniaxialMaterial - ";
	     cerr << " - no UniaxialMaterial type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}

SectionForceDeformation *
FEM_ObjectBroker::getNewSection(int classTag)
{
    switch(classTag) {
	case SEC_TAG_Elastic2d:
	     return new ElasticSection2d();
	     
	case SEC_TAG_Elastic3d:
	     return new ElasticSection3d();	     
	     
	case SEC_TAG_Generic1d:
	     return new GenericSection1d();
	     
	case SEC_TAG_GenericNd:
	     return new GenericSectionNd();	     

	case SEC_TAG_Aggregator:
	     return new SectionAggregator();

	case SEC_TAG_Fiber:
		return new FiberSection();
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewSection - ";
	     cerr << " - no section type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}

NDMaterial*
FEM_ObjectBroker::getNewNDMaterial(int classTag)
{
    switch(classTag) {
		case ND_TAG_ElasticIsotropicPlaneStrain2d:
			return new ElasticIsotropicPlaneStrain2D();
		
		case ND_TAG_ElasticIsotropicPlaneStress2d:
			return new ElasticIsotropicPlaneStress2D();
		
		default:
			cerr << "FEM_ObjectBroker::getNewNDMaterial - ";
			cerr << " - no NDMaterial type exists for class tag ";
			cerr << classTag << endl;
			return 0;   
	 }
}

Fiber*
FEM_ObjectBroker::getNewFiber(int classTag)
{
	switch(classTag) {
	case FIBER_TAG_Uniaxial2d:
		return new UniaxialFiber2d();

	case FIBER_TAG_Uniaxial3d:
		return new UniaxialFiber3d();

	default:
		cerr << "FEM_ObjectBroker::getNewFiber - ";
		cerr << " - no Fiber type exists for class tag ";
		cerr << classTag << endl;
		return 0;
	}
}

ConvergenceTest *
FEM_ObjectBroker::getNewConvergenceTest(int classTag)
{
    switch(classTag) {
	case CONVERGENCE_TEST_CTestNormUnbalance:  
	     return new CTestNormUnbalance();
	     
	case CONVERGENCE_TEST_CTestNormDispIncr:  
	     return new CTestNormDispIncr();
	     	     
	case CONVERGENCE_TEST_CTestEnergyIncr:  
	     return new CTestEnergyIncr();
	     	     	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getPtrNewConvergenceTest - ";
	     cerr << " - no ConvergenceTest type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}


LoadPattern *
FEM_ObjectBroker::getNewLoadPattern(int classTag)
{
    switch(classTag) {
	case PATTERN_TAG_LoadPattern:
	     return new LoadPattern();

	default:
	     cerr << "FEM_ObjectBroker::getPtrLoadPattern - ";
	     cerr << " - no Load type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}


GroundMotion *
FEM_ObjectBroker::getNewGroundMotion(int classTag)
{
    switch(classTag) {
	default:
	     cerr << "FEM_ObjectBroker::getPtrGroundMotion - ";
	     cerr << " - no Load type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}

TimeSeries *
FEM_ObjectBroker::getNewTimeSeries(int classTag)
{
    switch(classTag) {
        case TSERIES_TAG_LinearSeries:
	  return new LinearSeries;
      
        case TSERIES_TAG_RectangularSeries:
	  return new RectangularSeries;

        case TSERIES_TAG_PathTimeSeries:
	  return new PathTimeSeries;

        case TSERIES_TAG_PathSeries:
	  return new PathSeries;

        case TSERIES_TAG_ConstantSeries:
	  return new ConstantSeries;

		case TSERIES_TAG_TrigSeries:
	  return new TrigSeries;

	default:
	     cerr << "FEM_ObjectBroker::getPtrTimeSeries - ";
	     cerr << " - no Load type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}


Matrix	  *
FEM_ObjectBroker::getPtrNewMatrix(int classTag, int noRows, int noCols)
{
    switch(classTag) {
	case MATRIX_TAG_Matrix:  
	     return new Matrix(noRows,noCols);
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getPtrNewMatrix - ";
	     cerr << " - no NodalLoad type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}


Vector	  *
FEM_ObjectBroker::getPtrNewVector(int classTag, int size)
{
    switch(classTag) {
	case VECTOR_TAG_Vector:  
	     return new Vector(size);
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getPtrNewVector - ";
	     cerr << " - no Vector type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}


ID	          *
FEM_ObjectBroker::getPtrNewID(int classTag, int size)
{
    switch(classTag) {
	case ID_TAG_ID:  
	     return new ID(size);
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getPtrNewID - ";
	     cerr << " - no ID type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}



/*****************************************
 *
 * METHODS TO GET NEW ANALYSIS CLASSES
 *
 *****************************************/

ConstraintHandler   *
FEM_ObjectBroker::getNewConstraintHandler(int classTag)
{
    switch(classTag) {
	case HANDLER_TAG_PlainHandler:  
	     return new PlainHandler();
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewConstraintHandler - ";
	     cerr << " - no ConstraintHandler type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}


DOF_Numberer        *
FEM_ObjectBroker::getNewNumberer(int classTag)
{
    switch(classTag) {
	case NUMBERER_TAG_DOF_Numberer:  
	     return new DOF_Numberer();
	     
	     
	case NUMBERER_TAG_PlainNumberer:  
	     return new PlainNumberer();
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewConstraintHandler - ";
	     cerr << " - no ConstraintHandler type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }
}


AnalysisModel       *
FEM_ObjectBroker::getNewAnalysisModel(int classTag)
{
    switch(classTag) {
	case AnaMODEL_TAGS_AnalysisModel:  
	     return new AnalysisModel();
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewAnalysisModel - ";
	     cerr << " - no AnalysisModel type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}


EquiSolnAlgo        *
FEM_ObjectBroker::getNewEquiSolnAlgo(int classTag)
{
    switch(classTag) {
	case EquiALGORITHM_TAGS_Linear:  
	     return new Linear();
	     
	     
	case EquiALGORITHM_TAGS_NewtonRaphson:  
	     return new NewtonRaphson();
	     
	     
	case EquiALGORITHM_TAGS_ModifiedNewton:  
	     return new ModifiedNewton();
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewEquiSolnAlgo - ";
	     cerr << " - no EquiSolnAlgo type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }        
}


DomainDecompAlgo    *
FEM_ObjectBroker::getNewDomainDecompAlgo(int classTag)
{
    switch(classTag) {
	case DomDecompALGORITHM_TAGS_DomainDecompAlgo:  
	     return new DomainDecompAlgo();
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewDomainDecompAlgo - ";
	     cerr << " - no DomainDecompAlgo type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }
}


StaticIntegrator    *
FEM_ObjectBroker::getNewStaticIntegrator(int classTag)
{
    switch(classTag) {
	case INTEGRATOR_TAGS_LoadControl:  
	     return new LoadControl(1.0,1.0,1.0,.10); // must recvSelf
	     
	     
	case INTEGRATOR_TAGS_ArcLength:  
	     return new ArcLength(1.0);      // must recvSelf
	     	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewStaticIntegrator - ";
	     cerr << " - no StaticIntegrator type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }
}


TransientIntegrator *
FEM_ObjectBroker::getNewTransientIntegrator(int classTag)
{
    switch(classTag) {
	case INTEGRATOR_TAGS_Newmark:  
	     return new Newmark();
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewTransientIntegrator - ";
	     cerr << " - no TransientIntegrator type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }
}


IncrementalIntegrator *
FEM_ObjectBroker::getNewIncrementalIntegrator(int classTag)
{
    switch(classTag) {
	case INTEGRATOR_TAGS_LoadControl:  
	     return new LoadControl(1.0,1.0,1.0,1.0); // must recvSelf
	    
	     
	case INTEGRATOR_TAGS_ArcLength:  
	     return new ArcLength(1.0);      // must recvSelf
	     	     
	     
	case INTEGRATOR_TAGS_Newmark:  
	     return new Newmark();
	     
	     
	default:
	     cerr << "FEM_ObjectBroker::getNewIncrementalIntegrator - ";
	     cerr << " - no IncrementalIntegrator type exists for class tag ";
	     cerr << classTag << endl;
	     return 0;
	     
	 }
}


LinearSOESolver *
FEM_ObjectBroker::getNewLinearSolver(void)
{
    return lastLinearSolver;
}

LinearSOE *
FEM_ObjectBroker::getNewLinearSOE(int classTagSOE, 
				     int classTagSolver)
{
    LinearSOE *theSOE =0;
//    SlowLinearSOESolver *theSlowSolver =0;	    
    FullGenLinSolver *theGenSolver =0;
    BandGenLinSolver *theGenBandSolver =0;
    BandSPDLinSolver *theBandSPDSolver =0;
    ProfileSPDLinSolver *theProfileSPDSolver =0;    

    /*
      case LinSOE_TAGS_SlowLinearSOE:  
	if (classTagSolver == SOLVER_TAGS_SlowLinearSOESolver) {
	    theSlowSolver = new SlowLinearSOESolver();
	    theSOE = new SlowLinearSOE(*theSlowSolver);
	    lastLinearSolver = theSlowSolver;
	    return theSOE;
	} else {
	    cerr << "FEM_ObjectBroker::getNewLinearSOE - ";
	    cerr << " - no SlowLinearSOESolver type exists for class tag ";
	    cerr << classTagSolver << endl;
	    return 0;		 
	}
	
	*/

    
    switch(classTagSOE) {
      case LinSOE_TAGS_FullGenLinSOE:  

	if (classTagSolver == SOLVER_TAGS_FullGenLinLapackSolver) {
	    theGenSolver = new FullGenLinLapackSolver();
	    theSOE = new FullGenLinSOE(*theGenSolver);
	    lastLinearSolver = theGenSolver;
	    return theSOE;
	} else {
	    cerr << "FEM_ObjectBroker::getNewLinearSOE - ";
	    cerr << " - no FullGenLinSOESolver type exists for class tag ";
	    cerr << classTagSolver << endl;
	    return 0;		 
	}	     
	
					 
      case LinSOE_TAGS_BandGenLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_BandGenLinLapackSolver) {
	      theGenBandSolver = new BandGenLinLapackSolver();
	      theSOE = new BandGenLinSOE(*theGenBandSolver);
	      lastLinearSolver = theGenBandSolver;
	      return theSOE;
	  } else {
	      cerr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      cerr << " - no BandGenLinSolver type exists for class tag ";
	      cerr << classTagSolver << endl;
	      return 0;
	  }		     
	  
					 
	case LinSOE_TAGS_BandSPDLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_BandSPDLinLapackSolver) {
	      theBandSPDSolver = new BandSPDLinLapackSolver();
	      theSOE = new BandSPDLinSOE(*theBandSPDSolver);
	      lastLinearSolver = theBandSPDSolver;
	      return theSOE;
	  } else {
	      cerr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      cerr << " - no BandSPDLinSOESolver type exists for class tag ";
	      cerr << classTagSolver << endl;
	      return 0;
	  }	     
	  
					 
	case LinSOE_TAGS_ProfileSPDLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_ProfileSPDLinDirectSolver) {
	      theProfileSPDSolver = new ProfileSPDLinDirectSolver();
	      theSOE = new ProfileSPDLinSOE(*theProfileSPDSolver);
	      lastLinearSolver = theProfileSPDSolver;
	      return theSOE;
	  } else if (classTagSolver == SOLVER_TAGS_ProfileSPDLinSubstrSolver) {
	      theProfileSPDSolver = new ProfileSPDLinSubstrSolver();
	      theSOE = new ProfileSPDLinSOE(*theProfileSPDSolver);
	      lastLinearSolver = theProfileSPDSolver;
	      return 0;		 
	  }	
	  else {
	      cerr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      cerr << " - no ProfileSPD_LinSolver type exists for class tag ";
	      cerr << classTagSolver << endl;
	      return 0;		 
	  }	     
	  

	default:
	  cerr << "FEM_ObjectBroker::getNewLinearSOE - ";
	  cerr << " - no LinearSOE type exists for class tag ";
	  cerr << classTagSOE << endl;
	  return 0;
	  
      
    }
}




DomainSolver *
FEM_ObjectBroker::getNewDomainSolver(void)
{
    return lastDomainSolver;
}
    
LinearSOE *
FEM_ObjectBroker::getPtrNewDDLinearSOE(int classTagSOE, 
				       int classTagDDSolver)
{
    ProfileSPDLinSubstrSolver *theProfileSPDSolver =0;    

    switch(classTagSOE) {
      case LinSOE_TAGS_ProfileSPDLinSOE:  

	if (classTagDDSolver == SOLVER_TAGS_ProfileSPDLinSubstrSolver) {
	    theProfileSPDSolver = new ProfileSPDLinSubstrSolver();
	    LinearSOE *theSOE = new ProfileSPDLinSOE(*theProfileSPDSolver);
	    lastDomainSolver = theProfileSPDSolver;
	    return theSOE;		 
	}
	else {
	    cerr << "FEM_ObjectBroker::getNewLinearSOE - ";
	    cerr << " - no ProfileSPD Domain Solver type exists for class tag ";
	    cerr << classTagDDSolver << endl;
	    return 0;		 
	}	     
	
					    
      default:
	cerr << "FEM_ObjectBroker::getNewLinearSOE - ";
	cerr << " - no LinearSOE type exists for class tag ";
	cerr << classTagSOE << endl;
	return 0;
	
    }
}


DomainDecompositionAnalysis *
FEM_ObjectBroker::getNewDomainDecompAnalysis(int classTag, 
						Subdomain &theSubdomain)
{
    switch(classTag) {
      case DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis:  
	return new DomainDecompositionAnalysis(theSubdomain);
	
	
      default:
	cerr << "ObjectBroker::getNewDomainDecompAnalysis ";
	cerr << " - no DomainDecompAnalysis type exists for class tag " ;
	cerr << classTag << endl;
	return 0;
	
    }
}


Subdomain 	  *
FEM_ObjectBroker::getSubdomainPtr(int classTag)
{
    cerr << "FEM_ObjectBroker: NOT IMPLEMENTED YET";
    return 0;
}




