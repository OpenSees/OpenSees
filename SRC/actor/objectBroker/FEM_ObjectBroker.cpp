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
                                                                        
// $Revision: 1.30 $
// $Date: 2006-01-13 01:09:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/objectBroker/FEM_ObjectBroker.cpp,v $
                                                                        
                                                                        
// Written: fmk
// Revision: A
//
// Purpose: This file contains the class definition for FEM_ObjectBroker.
// FEM_ObjectBroker is is an object broker class for the finite element
// method. All methods are virtual to allow for subclasses; which can be
// used by programmers when introducing new subclasses of the main objects.
//
// What: "@(#) FEM_ObjectBroker.C, revA"


#include <FEM_ObjectBroker.h>

// ActorTypes
#include <ActorSubdomain.h>

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
#include <CableMaterial.h>
#include <ENTMaterial.h>
#include <MinMaxMaterial.h>

//PY springs: RWBoulanger and BJeremic
#include <PySimple1.h>
#include <TzSimple1.h>
#include <QzSimple1.h>
#include <PyLiq1.h>
#include <TzLiq1.h>


#include <FedeasBond1Material.h>
#include <FedeasBond2Material.h>
#include <FedeasConcr1Material.h>
#include <FedeasConcr2Material.h>
#include <FedeasConcr3Material.h>
#include <FedeasHardeningMaterial.h>
#include <FedeasHyster1Material.h>
#include <FedeasHyster2Material.h>
#include <FedeasSteel1Material.h>
#include <FedeasSteel2Material.h>

#include <DrainBilinearMaterial.h>
#include <DrainClough1Material.h>
#include <DrainClough2Material.h>
#include <DrainPinch1Material.h>

// Sections
#include <ElasticSection2d.h>
#include <ElasticSection3d.h>
#include <GenericSection1d.h>
//#include <GenericSectionNd.h>
#include <SectionAggregator.h>
//#include <FiberSection.h>
#include <FiberSection2d.h>
#include <FiberSection3d.h>
#include <ElasticPlateSection.h>
#include <ElasticMembranePlateSection.h>
#include <MembranePlateFiberSection.h>
#include <Bidirectional.h>

// NDMaterials
#include <ElasticIsotropicPlaneStrain2D.h>
#include <ElasticIsotropicPlaneStress2D.h>
#include <ElasticIsotropicPlateFiber.h>
#include <ElasticIsotropicAxiSymm.h>
#include <ElasticIsotropic3D.h>
#include <J2PlaneStrain.h>
#include <J2PlaneStress.h>
#include <J2PlateFiber.h>
#include <J2AxiSymm.h>
#include <J2ThreeDimensional.h>
#include <PlaneStressMaterial.h>
#include <PlateFiberMaterial.h>
#include <FeapMaterial03.h>

#include <FluidSolidPorousMaterial.h>
#include <PressureDependMultiYield.h>
#include <PressureIndependMultiYield.h>

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
#include <CorotTruss.h>
#include <CorotTrussSection.h>
#include <ZeroLength.h>
#include <ZeroLengthSection.h>
//#include <ZeroLengthND.h>
#include <FourNodeQuad.h>
#include <EnhancedQuad.h>
#include <NineNodeMixedQuad.h>
#include <ConstantPressureVolumeQuad.h>
#include <ElasticBeam2d.h>
#include <ElasticBeam3d.h>
#include <BeamWithHinges2d.h>
#include <BeamWithHinges3d.h>
#include <NLBeamColumn2d.h>
#include <NLBeamColumn3d.h>
#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn3d.h>

#include <DispBeamColumn2d.h>
#include <DispBeamColumn3d.h>
#include <ShellMITC4.h>
#include <Brick.h>
#include <BbarBrick.h>
#include <Joint2D.h>		// Arash


#include <LinearCrdTransf2d.h>
#include <LinearCrdTransf3d.h>
#include <PDeltaCrdTransf2d.h>
#include <PDeltaCrdTransf3d.h>
#include <CorotCrdTransf2d.h>
#include <CorotCrdTransf3d.h>

#include <HingeMidpointBeamIntegration2d.h>
#include <HingeMidpointBeamIntegration3d.h>
#include <HingeRadauBeamIntegration2d.h>
#include <HingeRadauBeamIntegration3d.h>
#include <HingeRadauTwoBeamIntegration2d.h>
#include <HingeRadauTwoBeamIntegration3d.h>
#include <LobattoBeamIntegration.h>

// node header files
#include <Node.h>


#include <DataOutputStreamHandler.h>
#include <DataOutputFileHandler.h>
#include <DataOutputDatabaseHandler.h>

#include <NodeRecorder.h>
#include <ElementRecorder.h>
#include <EnvelopeNodeRecorder.h>
#include <EnvelopeElementRecorder.h>

// mp_constraint header files
#include <MP_Constraint.h>
#include <MP_Joint2D.h>

// sp_constraint header files
#include <SP_Constraint.h>
#include <SP_Constraint.h>
#include <ImposedMotionSP.h>
#include <ImposedMotionSP1.h>

// nodal load header files
#include <NodalLoad.h>

// elemental load header files
#include <ElementalLoad.h>
#include<Beam2dUniformLoad.h>
#include<Beam2dPointLoad.h>
#include<Beam3dUniformLoad.h>
#include<Beam3dPointLoad.h>
#include<BrickSelfWeight.h>

// matrix, vector & id header files
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

// subdomain header files
#include <Subdomain.h>

// constraint handler header files
#include <ConstraintHandler.h>
#include <PlainHandler.h>
#include <PenaltyConstraintHandler.h>
#include <LagrangeConstraintHandler.h>
#include <TransformationConstraintHandler.h>

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
#include <DisplacementControl.h>
#include <CentralDifferenceNoDamping.h>
#include <CentralDifferenceAlternative.h>

#ifdef _PARALLEL_PROCESSING
#include <DistributedDisplacementControl.h>
#endif
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

#include <SparseGenColLinSOE.h>
#include <SuperLU.h>

#include <DomainDecompositionAnalysis.h>

// load patterns
#include <LoadPattern.h>
#include <UniformExcitation.h>
#include <MultiSupportPattern.h>
#include <GroundMotion.h>
#include <InterpolatedGroundMotion.h>

// time series
#include <LinearSeries.h>
#include <PathSeries.h>
#include <PathTimeSeries.h>
#include <RectangularSeries.h>
#include <ConstantSeries.h>
#include <TrigSeries.h>

// time series integrators
#include <TrapezoidalTimeSeriesIntegrator.h>

#ifdef _PETSC
#include <PetscSOE.h>
#include <PetscSolver.h>
#include <SparseGenColLinSOE.h>
#include <PetscSparseSeqSolver.h>
#endif

#ifdef _PARALLEL_PROCESSING
#include <DistributedBandSPDLinSOE.h>
#include <DistributedProfileSPDLinSOE.h>
#include <DistributedSparseGenColLinSOE.h>
#include <DistributedSparseGenRowLinSOE.h>
#include <DistributedSparseGenRowLinSolver.h>
#include <DistributedBandGenLinSOE.h>
#include <DistributedSuperLU.h>
#include <ParallelNumberer.h>
#include <StaticDomainDecompositionAnalysis.h>
#include <TransientDomainDecompositionAnalysis.h>
#include <DistributedDiagonalSOE.h>
#include <DistributedDiagonalSolver.h>
#endif

#include <packages.h>

typedef struct uniaxialPackage {
  int classTag;
  char *libName;
  char *funcName;
  UniaxialMaterial *(*funcPtr)(void);
  struct uniaxialPackage *next;
} UniaxialPackage;

static UniaxialPackage *theUniaxialPackage = NULL;



FEM_ObjectBroker::FEM_ObjectBroker()
:lastLinearSolver(0),lastDomainSolver(0)
{

}


FEM_ObjectBroker::~FEM_ObjectBroker()
{

}


Actor *
FEM_ObjectBroker::getNewActor(int classTag, Channel *theChannel)
{
  switch(classTag) {

#ifdef _PARALLEL_PROCESSING
  case ACTOR_TAGS_SUBDOMAIN:  
    return new ActorSubdomain(*theChannel, *this);
#endif

  default:
    opserr << "FEM_ObjectBroker::getNewActor - ";
    opserr << " - no ActorType type exists for class tag ";
    opserr << classTag << endln;
    return 0;
  }
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
	     opserr << "FEM_ObjectBroker::getPtrNewPartitionedModelBuilder - ";
	     opserr << " - no PartitionedModelBuilder type exists for class tag ";
	     opserr << classTag << endln;
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
	     opserr << "ObjectBroker::getPtrNewGraphNumberer - ";
	     opserr << " - no GraphNumberer type exists for class tag " ;
	     opserr << classTag << endln;
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
	     
	case ELE_TAG_CorotTruss:  
	     return new CorotTruss(); 
	     
	case ELE_TAG_CorotTrussSection:  
	     return new CorotTrussSection(); 	     

	case ELE_TAG_ZeroLength:  
	     return new ZeroLength(); 	     

	case ELE_TAG_ZeroLengthSection:  
	     return new ZeroLengthSection(); 	     

	     //case ELE_TAG_ZeroLengthND:  
	     //return new ZeroLengthND(); 	     
	     
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

	case ELE_TAG_ForceBeamColumn2d:  
	     return new ForceBeamColumn2d();					     

	case ELE_TAG_ForceBeamColumn3d:  
	     return new ForceBeamColumn3d();  
				
	case ELE_TAG_DispBeamColumn2d:  
	     return new DispBeamColumn2d();					     

	case ELE_TAG_DispBeamColumn3d:  
	     return new DispBeamColumn3d(); 
		 
	case ELE_TAG_EnhancedQuad:
		return new EnhancedQuad();

        case ELE_TAG_NineNodeMixedQuad:
		return new NineNodeMixedQuad();

	case ELE_TAG_ConstantPressureVolumeQuad:
		return new ConstantPressureVolumeQuad();

	case ELE_TAG_Brick:
		return new Brick();

	case ELE_TAG_ShellMITC4:
		return new ShellMITC4();

	case ELE_TAG_BbarBrick:
		return new BbarBrick();
			
	case ELE_TAG_Joint2D:				// Arash
		return new Joint2D();			// Arash
	
	default:
	     opserr << "FEM_ObjectBroker::getNewElement - ";
	     opserr << " - no Element type exists for class tag " ;
	     opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getNewNode - ";
	     opserr << " - no Node type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }    
}


MP_Constraint *
FEM_ObjectBroker::getNewMP(int classTag)
{
    switch(classTag) {
	case CNSTRNT_TAG_MP_Constraint:  
	     return new MP_Constraint( 0 , classTag);

 	case CNSTRNT_TAG_MP_Joint2D:			// Arash
	     return new MP_Joint2D();			// Arash
	
	default:
	     opserr << "FEM_ObjectBroker::getNewMP - ";
	     opserr << " - no MP_Constraint type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }    
}


SP_Constraint *
FEM_ObjectBroker::getNewSP(int classTag)
{
    switch(classTag) {
	case CNSTRNT_TAG_SP_Constraint:  
	     return new SP_Constraint(classTag);

	case CNSTRNT_TAG_ImposedMotionSP:  
	     return new ImposedMotionSP();

	case CNSTRNT_TAG_ImposedMotionSP1:  
	     return new ImposedMotionSP1();
	     
	default:
	     opserr << "FEM_ObjectBroker::getNewSP - ";
	     opserr << " - no SP_Constraint type exists for class tag ";
	     opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getNewNodalLoad - ";
	     opserr << " - no NodalLoad type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }    
}


ElementalLoad *
FEM_ObjectBroker::getNewElementalLoad(int classTag)
{
  switch(classTag) {
    
    case LOAD_TAG_Beam2dUniformLoad:
      return new Beam2dUniformLoad();
    
    case LOAD_TAG_Beam2dPointLoad:
      return new Beam2dPointLoad();
    
    case LOAD_TAG_Beam3dUniformLoad:
      return new Beam3dUniformLoad();
    
    case LOAD_TAG_Beam3dPointLoad:
      return new Beam3dPointLoad();
    
    case LOAD_TAG_BrickSelfWeight:
      return new BrickSelfWeight();	     
	     
  default:
    opserr << "FEM_ObjectBroker::getNewNodalLoad - ";
    opserr << " - no NodalLoad type exists for class tag ";
    opserr << classTag << endln;
    return 0;
    
  }    
  
  return 0;
}

CrdTransf2d*
FEM_ObjectBroker::getNewCrdTransf2d(int classTag)
{
	switch(classTag) {
	case CRDTR_TAG_LinearCrdTransf2d:
		return new LinearCrdTransf2d();
	case CRDTR_TAG_PDeltaCrdTransf2d:
		return new PDeltaCrdTransf2d();
#ifdef _COROTATIONAL
	case CRDTR_TAG_CorotCrdTransf2d:
		return new CorotCrdTransf2d();
#endif
	default:
		opserr << "FEM_ObjectBroker::getCrdTransf2d - ";
	    opserr << " - no CrdTransf2d type exists for class tag ";
	    opserr << classTag << endln;
	    return 0;
	}

}

CrdTransf3d*
FEM_ObjectBroker::getNewCrdTransf3d(int classTag)
{
	switch(classTag) {
	case CRDTR_TAG_LinearCrdTransf3d:
		return new LinearCrdTransf3d();
	case CRDTR_TAG_PDeltaCrdTransf3d:
		return new PDeltaCrdTransf3d();
#ifdef _COROTATIONAL
	case CRDTR_TAG_CorotCrdTransf3d:
		return new CorotCrdTransf3d();
#endif
	default:
		opserr << "FEM_ObjectBroker::getCrdTransf3d - ";
	    opserr << " - no CrdTransf3d type exists for class tag ";
	    opserr << classTag << endln;
	    return 0;
	}

}

BeamIntegration *
FEM_ObjectBroker::getNewBeamIntegration(int classTag)
{
  switch(classTag) {
    case BEAM_INTEGRATION_TAG_Lobatto:        
      return new LobattoBeamIntegration();

    case BEAM_INTEGRATION_TAG_HingeMidpoint2d:
      return new HingeMidpointBeamIntegration2d();

    case BEAM_INTEGRATION_TAG_HingeRadau2d:
      return new HingeRadauBeamIntegration2d();

    case BEAM_INTEGRATION_TAG_HingeRadauTwo2d:
      return new HingeRadauTwoBeamIntegration2d();

    case BEAM_INTEGRATION_TAG_HingeMidpoint3d:
      return new HingeMidpointBeamIntegration3d();

    case BEAM_INTEGRATION_TAG_HingeRadau3d:   
      return new HingeRadauBeamIntegration3d();

    case BEAM_INTEGRATION_TAG_HingeRadauTwo3d:   
      return new HingeRadauTwoBeamIntegration3d();

    default:
      opserr << "FEM_ObjectBroker::getBeamIntegration - ";
      opserr << " - no BeamIntegration type exists for class tag ";
      opserr << classTag << endln;
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

	//PY springs: RWBoulanger and BJeremic
	case MAT_TAG_PySimple1:
		return new PySimple1();

	case MAT_TAG_PyLiq1:
		return new PyLiq1();

	case MAT_TAG_TzSimple1:
		return new TzSimple1();

	case MAT_TAG_TzLiq1:
		return new TzLiq1();

	case MAT_TAG_QzSimple1:
		return new QzSimple1();


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

	case MAT_TAG_CableMaterial:
		return new CableMaterial();
	     
	case MAT_TAG_ENTMaterial:
		return new ENTMaterial();

	case MAT_TAG_FedeasBond1:
		return new FedeasBond1Material();

	case MAT_TAG_FedeasBond2:
		return new FedeasBond2Material();

	case MAT_TAG_FedeasConcrete1:
		return new FedeasConcr1Material();

	case MAT_TAG_FedeasConcrete2:
		return new FedeasConcr2Material();

	case MAT_TAG_FedeasConcrete3:
		return new FedeasConcr3Material();

	case MAT_TAG_FedeasHardening:
		return new FedeasHardeningMaterial();

	case MAT_TAG_FedeasHysteretic1:
		return new FedeasHyster1Material();

	case MAT_TAG_FedeasHysteretic2:
		return new FedeasHyster2Material();

	case MAT_TAG_FedeasSteel1:
		return new FedeasSteel1Material();

	case MAT_TAG_FedeasSteel2:
		return new FedeasSteel2Material();

	case MAT_TAG_DrainBilinear:
		return new DrainBilinearMaterial();

	case MAT_TAG_DrainClough1:
		return new DrainClough1Material();

	case MAT_TAG_DrainClough2:
		return new DrainClough2Material();

	case MAT_TAG_DrainPinch1:
		return new DrainPinch1Material();

        case MAT_TAG_MinMax:
	  return new MinMaxMaterial();

	default:

	  UniaxialPackage *matCommands = theUniaxialPackage;
	  bool found = false;
	  while (matCommands != NULL && found == false) {
	    if ((matCommands->classTag == classTag) && (matCommands->funcPtr != 0)){
	      UniaxialMaterial *result = (*(matCommands->funcPtr))();
	      return result;
	    } 
	    matCommands = matCommands->next;
	  }	  

	  opserr << "FEM_ObjectBroker::getNewUniaxialMaterial - ";
	  opserr << " - no UniaxialMaterial type exists for class tag ";
	  opserr << classTag << endln;
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
	     
	     //case SEC_TAG_GenericNd:
	     //return new GenericSectionNd();	     

	case SEC_TAG_Aggregator:
	     return new SectionAggregator();

	     //case SEC_TAG_Fiber:
	     //return new FiberSection();
	
	case SEC_TAG_FiberSection2d:
		return new FiberSection2d();
      
	case SEC_TAG_FiberSection3d:
		return new FiberSection3d();

	case SEC_TAG_ElasticPlateSection:
		return new ElasticPlateSection();

	case SEC_TAG_ElasticMembranePlateSection:
		return new ElasticMembranePlateSection();

	case SEC_TAG_MembranePlateFiberSection:
		return new MembranePlateFiberSection();

	case SEC_TAG_Bidirectional:
		return new Bidirectional();

	default:
	     opserr << "FEM_ObjectBroker::getNewSection - ";
	     opserr << " - no section type exists for class tag ";
	     opserr << classTag << endln;
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
		
  case ND_TAG_ElasticIsotropicAxiSymm:
    return new ElasticIsotropicAxiSymm();
    
  case ND_TAG_ElasticIsotropicPlateFiber:
    return new ElasticIsotropicPlateFiber();
    
  case ND_TAG_ElasticIsotropic3D:
    return new ElasticIsotropic3D();
		  
  case ND_TAG_J2PlaneStrain:
    return new J2PlaneStrain();
    
  case ND_TAG_J2PlaneStress:
    return new J2PlaneStress();
    
  case ND_TAG_J2AxiSymm:
    return new J2AxiSymm();
    
  case ND_TAG_J2PlateFiber:
    return new J2PlateFiber();
    
  case ND_TAG_J2ThreeDimensional:
    return new J2ThreeDimensional();
    
  case ND_TAG_PlaneStressMaterial:
    return new PlaneStressMaterial();
		  
  case ND_TAG_PlateFiberMaterial:
    return new PlateFiberMaterial();
    
  case ND_TAG_FluidSolidPorousMaterial:
    return new FluidSolidPorousMaterial();

  case ND_TAG_PressureDependMultiYield:
    return new PressureDependMultiYield();

  case ND_TAG_PressureIndependMultiYield:
    return new PressureIndependMultiYield();

  case ND_TAG_FeapMaterial03:
    return new FeapMaterial03();
    
  default:
    opserr << "FEM_ObjectBroker::getNewNDMaterial - ";
    opserr << " - no NDMaterial type exists for class tag ";
    opserr << classTag << endln;
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
		opserr << "FEM_ObjectBroker::getNewFiber - ";
		opserr << " - no Fiber type exists for class tag ";
		opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getNewConvergenceTest - ";
	     opserr << " - no ConvergenceTest type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}


LoadPattern *
FEM_ObjectBroker::getNewLoadPattern(int classTag)
{
    switch(classTag) {
	case PATTERN_TAG_LoadPattern:
	     return new LoadPattern();

	case PATTERN_TAG_UniformExcitation:
	     return new UniformExcitation();

	case PATTERN_TAG_MultiSupportPattern:
	     return new MultiSupportPattern();

	default:
	     opserr << "FEM_ObjectBroker::getPtrLoadPattern - ";
	     opserr << " - no Load type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}


GroundMotion *
FEM_ObjectBroker::getNewGroundMotion(int classTag)
{
    switch(classTag) {

        case GROUND_MOTION_TAG_GroundMotion:
	  return new GroundMotion(GROUND_MOTION_TAG_GroundMotion);

        case GROUND_MOTION_TAG_InterpolatedGroundMotion:
	  return new GroundMotion(GROUND_MOTION_TAG_InterpolatedGroundMotion);

	default:
	     opserr << "FEM_ObjectBroker::getPtrGroundMotion - ";
	     opserr << " - no Load type exists for class tag ";
	     opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getPtrTimeSeries - ";
	     opserr << " - no Load type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}

TimeSeriesIntegrator *
FEM_ObjectBroker::getNewTimeSeriesIntegrator(int classTag)
{
    switch(classTag) {
    case TIMESERIES_INTEGRATOR_TAG_Trapezoidal:
	  return new TrapezoidalTimeSeriesIntegrator();

	default:
	     opserr << "FEM_ObjectBroker::getPtrTimeSeriesIntegrator - ";
	     opserr << " - no Load type exists for class tag ";
	     opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getPtrNewMatrix - ";
	     opserr << " - no NodalLoad type exists for class tag ";
	     opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getPtrNewVector - ";
	     opserr << " - no Vector type exists for class tag ";
	     opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getPtrNewID - ";
	     opserr << " - no ID type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}

/*****************************************
 *
 * METHODS TO GET NEW OUTPUT CLASS OBJECTS
 *
 *****************************************/

DataOutputHandler *
FEM_ObjectBroker::getPtrNewDataOutputHandler(int classTag)
{
    switch(classTag) {
	case DATAHANDLER_TAGS_DataOutputStreamHandler:  
	     return new DataOutputStreamHandler();

	case DATAHANDLER_TAGS_DataOutputFileHandler:  
	     return new DataOutputFileHandler();

	case DATAHANDLER_TAGS_DataOutputDatabaseHandler:  
	     return new DataOutputDatabaseHandler();
	     
	default:
	     opserr << "FEM_ObjectBroker::getPtrNewDataOutputHandler - ";
	     opserr << " - no DataOutputHandler type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }        
}

Recorder *
FEM_ObjectBroker::getPtrNewRecorder(int classTag)
{
    switch(classTag) {
	case RECORDER_TAGS_ElementRecorder:  
	     return new ElementRecorder();

	case RECORDER_TAGS_NodeRecorder:  
	     return new NodeRecorder();

	case RECORDER_TAGS_EnvelopeNodeRecorder:  
	     return new EnvelopeNodeRecorder();

	case RECORDER_TAGS_EnvelopeElementRecorder:  
	     return new EnvelopeElementRecorder();
	     
	default:
	     opserr << "FEM_ObjectBroker::getNewConstraintHandler - ";
	     opserr << " - no ConstraintHandler type exists for class tag ";
	     opserr << classTag << endln;
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
	     
	case HANDLER_TAG_PenaltyConstraintHandler:  
	     return new PenaltyConstraintHandler(1.0e12, 1.0e12);

	case HANDLER_TAG_LagrangeConstraintHandler:  
	     return new LagrangeConstraintHandler(1.0, 1.0);

	case HANDLER_TAG_TransformationConstraintHandler:  
	     return new TransformationConstraintHandler();
	     
	default:
	     opserr << "FEM_ObjectBroker::getNewConstraintHandler - ";
	     opserr << " - no ConstraintHandler type exists for class tag ";
	     opserr << classTag << endln;
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


#ifdef _PARALLEL_PROCESSING
	case NUMBERER_TAG_ParallelNumberer:  
	     return new ParallelNumberer();
#endif
	     
	default:
	     opserr << "FEM_ObjectBroker::getNewConstraintHandler - ";
	     opserr << " - no ConstraintHandler type exists for class tag ";
	     opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getNewAnalysisModel - ";
	     opserr << " - no AnalysisModel type exists for class tag ";
	     opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getNewEquiSolnAlgo - ";
	     opserr << " - no EquiSolnAlgo type exists for class tag ";
	     opserr << classTag << endln;
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
	     opserr << "FEM_ObjectBroker::getNewDomainDecompAlgo - ";
	     opserr << " - no DomainDecompAlgo type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}


StaticIntegrator    *
FEM_ObjectBroker::getNewStaticIntegrator(int classTag)
{
    switch(classTag) {
	case INTEGRATOR_TAGS_LoadControl:  
	     return new LoadControl(1.0,1,1.0,.10); // must recvSelf

#ifdef _PARALLEL_PROCESSING
	case INTEGRATOR_TAGS_DistributedDisplacementControl:  
	     return new DistributedDisplacementControl(); // must recvSelf
#endif	     
	     
	case INTEGRATOR_TAGS_ArcLength:  
	     return new ArcLength(1.0);      // must recvSelf

	     
	default:
	     opserr << "FEM_ObjectBroker::getNewStaticIntegrator - ";
	     opserr << " - no StaticIntegrator type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}


TransientIntegrator *
FEM_ObjectBroker::getNewTransientIntegrator(int classTag)
{
    switch(classTag) {
	case INTEGRATOR_TAGS_Newmark:  
	     return new Newmark();

	case INTEGRATOR_TAGS_CentralDifferenceNoDamping:  
	     return new CentralDifferenceNoDamping();      // must recvSelf

	case INTEGRATOR_TAGS_CentralDifferenceAlternative:  
	     return new CentralDifferenceAlternative();      // must recvSelf
	     	     
	default:
	     opserr << "FEM_ObjectBroker::getNewTransientIntegrator - ";
	     opserr << " - no TransientIntegrator type exists for class tag ";
	     opserr << classTag << endln;
	     return 0;
	     
	 }
}


IncrementalIntegrator *
FEM_ObjectBroker::getNewIncrementalIntegrator(int classTag)
{
    switch(classTag) {
	case INTEGRATOR_TAGS_LoadControl:  
	     return new LoadControl(1.0,1,1.0,1.0); // must recvSelf
	    
	     
	case INTEGRATOR_TAGS_ArcLength:  
	     return new ArcLength(1.0);      // must recvSelf
	     	     
	     
	case INTEGRATOR_TAGS_Newmark:  
	     return new Newmark();

#ifdef _PARALLEL_PROCESSING	     
	case INTEGRATOR_TAGS_DistributedDisplacementControl:  
	     return new DistributedDisplacementControl(); // must recvSelf
#endif
	     
	default:
	     opserr << "FEM_ObjectBroker::getNewIncrementalIntegrator - ";
	     opserr << " - no IncrementalIntegrator type exists for class tag ";
	     opserr << classTag << endln;
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
    SuperLU *theSparseGenLinSolver =0;

#ifdef _PETSC
    PetscSolver *thePetscSolver = 0;
#endif

#ifdef _PARALLEL_PROCESSING
    DistributedSuperLU *theDistributedSparseGenLinSolver =0;
    DistributedDiagonalSolver *theDistributedDiagonalSolver =0;
#endif    

    /*
      case LinSOE_TAGS_SlowLinearSOE:  
	if (classTagSolver == SOLVER_TAGS_SlowLinearSOESolver) {
	    theSlowSolver = new SlowLinearSOESolver();
	    theSOE = new SlowLinearSOE(*theSlowSolver);
	    lastLinearSolver = theSlowSolver;
	    return theSOE;
	} else {
	    opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	    opserr << " - no SlowLinearSOESolver type exists for class tag ";
	    opserr << classTagSolver << endln;
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
	    opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	    opserr << " - no FullGenLinSOESolver type exists for class tag ";
	    opserr << classTagSolver << endln;
	    return 0;		 
	}	     
	
					 
      case LinSOE_TAGS_BandGenLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_BandGenLinLapackSolver) {
	      theGenBandSolver = new BandGenLinLapackSolver();
	      theSOE = new BandGenLinSOE(*theGenBandSolver);
	      lastLinearSolver = theGenBandSolver;
	      return theSOE;
	  } else {
	      opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      opserr << " - no BandGenLinSolver type exists for class tag ";
	      opserr << classTagSolver << endln;
	      return 0;
	  }		     

	case LinSOE_TAGS_BandSPDLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_BandSPDLinLapackSolver) {
	      theBandSPDSolver = new BandSPDLinLapackSolver();
	      theSOE = new BandSPDLinSOE(*theBandSPDSolver);
	      lastLinearSolver = theBandSPDSolver;
	      return theSOE;
	  } else {
	      opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      opserr << " - no BandSPDLinSOESolver type exists for class tag ";
	      opserr << classTagSolver << endln;
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
	      opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      opserr << " - no ProfileSPD_LinSolver type exists for class tag ";
	      opserr << classTagSolver << endln;
	      return 0;		 
	  }	     


#ifdef _PETSC
      case LinSOE_TAGS_PetscSOE:  
	
	  if (classTagSolver == SOLVER_TAGS_PetscSolver) {
	      thePetscSolver = new PetscSolver();
	      theSOE = new PetscSOE(*thePetscSolver);
	      return theSOE;
	  } else {
	      opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      opserr << " - no PetscSolver type exists for class tag ";
	      opserr << classTagSolver << endln;
	      return 0;
	  }		     
#endif

#ifdef _PARALLEL_PROCESSING
      case LinSOE_TAGS_DistributedBandGenLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_BandGenLinLapackSolver) {
	      theGenBandSolver = new BandGenLinLapackSolver();
	      theSOE = new DistributedBandGenLinSOE(*theGenBandSolver);
	      lastLinearSolver = theGenBandSolver;
	      return theSOE;
	  } else {
	      opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      opserr << " - no DistributedBandGenLinSolver type exists for class tag ";
	      opserr << classTagSolver << endln;
	      return 0;
	  }		     

        case LinSOE_TAGS_DistributedBandSPDLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_BandSPDLinLapackSolver) {
	      theBandSPDSolver = new BandSPDLinLapackSolver();
	      theSOE = new DistributedBandSPDLinSOE(*theBandSPDSolver);
	      lastLinearSolver = theBandSPDSolver;
	      return theSOE;
	  } else {
	      opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      opserr << " - no DistributedBandSPDLinSolver type exists for class tag ";
	      opserr << classTagSolver << endln;
	      return 0;
	  }		     
	  

	case LinSOE_TAGS_DistributedProfileSPDLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_ProfileSPDLinDirectSolver) {
	      theProfileSPDSolver = new ProfileSPDLinDirectSolver();
	      theSOE = new DistributedProfileSPDLinSOE(*theProfileSPDSolver);
	      lastLinearSolver = theProfileSPDSolver;
	      return theSOE;
	  } else if (classTagSolver == SOLVER_TAGS_ProfileSPDLinSubstrSolver) {
	      theProfileSPDSolver = new ProfileSPDLinSubstrSolver();
	      theSOE = new DistributedProfileSPDLinSOE(*theProfileSPDSolver);
	      lastLinearSolver = theProfileSPDSolver;
	      return 0;		 
	  }	
	  else {
	      opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      opserr << " - no ProfileSPD_LinSolver type exists for class tag ";
	      opserr << classTagSolver << endln;
	      return 0;		 
	  }	     
	  
	case LinSOE_TAGS_DistributedDiagonalSOE:  

	  if (classTagSolver == SOLVER_TAGS_DistributedDiagonalSolver) {
	    theDistributedDiagonalSolver = new DistributedDiagonalSolver();
	    theSOE = new DistributedDiagonalSOE(*theDistributedDiagonalSolver);
	    lastLinearSolver = theDistributedDiagonalSolver;
	    return theSOE;
	  } else {
	    opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	    opserr << " - no DistributedSparseGenLinSolverSolver type exists for class tag ";
	    opserr << classTagSolver << endln;
	    return 0;
	  }	     

	case LinSOE_TAGS_DistributedSparseGenColLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_SuperLU) {
	      theSparseGenLinSolver = new SuperLU();
	      theSOE = new DistributedSparseGenColLinSOE(*theSparseGenLinSolver);
	      lastLinearSolver = theSparseGenLinSolver;
	      return theSOE;
	  } else if (classTagSolver == SOLVER_TAGS_DistributedSuperLU) {
	      theDistributedSparseGenLinSolver = new DistributedSuperLU();
	      theSOE = new DistributedSparseGenColLinSOE(*theDistributedSparseGenLinSolver);
	      lastLinearSolver = theSparseGenLinSolver;
	      return theSOE;
	  } else {
	      opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      opserr << " - no DistributedSparseGenLinSolverSolver type exists for class tag ";
	      opserr << classTagSolver << endln;
	      return 0;
	  }	     

#else
	case LinSOE_TAGS_SparseGenColLinSOE:  

	  if (classTagSolver == SOLVER_TAGS_SuperLU) {
	      theSparseGenLinSolver = new SuperLU();
	      theSOE = new SparseGenColLinSOE(*theSparseGenLinSolver);
	      lastLinearSolver = theSparseGenLinSolver;
	      return theSOE;
	  } else {
	      opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	      opserr << " - no SparseGenLinSolverSolver type exists for class tag ";
	      opserr << classTagSolver << endln;
	      return 0;
	  }	     


#endif

	default:
	  opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	  opserr << " - no LinearSOE type exists for class tag ";
	  opserr << classTagSOE << endln;
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
	    opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	    opserr << " - no ProfileSPD Domain Solver type exists for class tag ";
	    opserr << classTagDDSolver << endln;
	    return 0;		 
	}	     
	
					    
      default:
	opserr << "FEM_ObjectBroker::getNewLinearSOE - ";
	opserr << " - no LinearSOE type exists for class tag ";
	opserr << classTagSOE << endln;
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

#ifdef _PARALLEL_PROCESSING
      case ANALYSIS_TAGS_StaticDomainDecompositionAnalysis:
	return new StaticDomainDecompositionAnalysis(theSubdomain);      

      case ANALYSIS_TAGS_TransientDomainDecompositionAnalysis:
	return new TransientDomainDecompositionAnalysis(theSubdomain);      
#endif
	
      default:
	opserr << "ObjectBroker::getNewDomainDecompAnalysis ";
	opserr << " - no DomainDecompAnalysis type exists for class tag " ;
	opserr << classTag << endln;
	return 0;
	
    }
}


Subdomain 	  *
FEM_ObjectBroker::getSubdomainPtr(int classTag)
{
    opserr << "FEM_ObjectBroker: NOT IMPLEMENTED YET";
    return 0;
}


int 
FEM_ObjectBroker::addUniaxialMaterial(int classTag, 
				      const char *lib, 
				      const char *funcName, 
				      UniaxialMaterial *(*funcPtr)(void))
{
  // check to see if it's already added

  UniaxialPackage *matCommands = theUniaxialPackage;
  bool found = false;
  while (matCommands != NULL && found == false) {
    if ((strcmp(lib, matCommands->libName) == 0) && (strcmp(funcName, matCommands->funcName) == 0)) {
      return 0;
    }
  }

  //
  // if funPtr == 0; go get the handle
  //

  void *libHandle;
  if (funcPtr == 0) {
    if (getLibraryFunction(lib, funcName, &libHandle, (void **)&funcPtr) != 0) {
      opserr << "FEM_ObjectBroker::addUniaxialMaterial - could not find function\n";
      return -1;
    }
  } 
  
  //
  // add the new funcPtr
  //
  
  char *libNameCopy = new char[strlen(lib)+1];
  char *funcNameCopy = new char[strlen(funcName)+1];
  UniaxialPackage *theMat = new UniaxialPackage;
  if (libNameCopy == 0 || funcNameCopy == 0 || theMat == 0) {
      opserr << "FEM_ObjectBroker::addUniaxialMaterial - could not add lib, out of memory\n";
      return -1;
  }
  strcpy(libNameCopy, lib);
  strcpy(funcNameCopy, funcName);

  theMat->classTag = classTag;	
  theMat->funcName = funcNameCopy;	
  theMat->libName = libNameCopy;	
  theMat->funcPtr = funcPtr;
  theMat->next = theUniaxialPackage;
  theUniaxialPackage = theMat;

  return 0;

}
