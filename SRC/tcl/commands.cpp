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
                                                                        
// $Revision: 1.77 $
// $Date: 2006-05-26 18:30:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/commands.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 04/98
//
// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified,
// see tkAppInit.C for command names.
//
// What: "@(#) commands.C, revA"

extern "C" {
#include <tcl.h>
EXTERN int      Tcl_SetObjCmd _ANSI_ARGS_((ClientData clientData,
					   Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[]));

}

#include <OPS_Globals.h>
#include <SimulationInformation.h>
extern SimulationInformation simulationInfo;

// the following is a little kludgy but it works!
#ifdef _USING_STL_STREAMS

#include <iomanip>
using std::ios;
#include <iostream>
using std::ofstream;

#else

#include <StandardStream.h>
#include <FileStream.h>
StandardStream sserr;
//OPS_Stream &opserr = sserr;
OPS_Stream *opserrPtr = &sserr;

#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <packages.h>

#include <FEM_ObjectBroker.h>
#include <RigidRod.h>
#include <RigidBeam.h>
#include <RigidDiaphragm.h>

#include <Timer.h>
#include <ModelBuilder.h>
#include "commands.h"

// domain
 #ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
#else
#include <Domain.h>
#endif

#include <Element.h>
#include <Node.h>
#include <ElementIter.h>
#include <NodeIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>
#include <SP_Constraint.h> //Joey UC Davis
#include <SP_ConstraintIter.h> //Joey UC Davis


// analysis model
#include <AnalysisModel.h>

// convergence tests
#include <CTestNormUnbalance.h>
#include <CTestNormDispIncr.h>
#include <CTestEnergyIncr.h>
#include <CTestRelativeNormUnbalance.h>
#include <CTestRelativeNormDispIncr.h>
#include <CTestRelativeEnergyIncr.h>
#include <CTestRelativeTotalNormDispIncr.h>
#include <CTestFixedNumIter.h>

// soln algorithms
#include <Linear.h>
#include <NewtonRaphson.h>
#include <NewtonLineSearch.h>
#include <ModifiedNewton.h>
#include <FrequencyAlgo.h>
#include <StandardEigenAlgo.h>
#include <Broyden.h>
#include <BFGS.h>
#include <KrylovNewton.h>
#include <PeriodicNewton.h>


// line searches
#include <BisectionLineSearch.h>
#include <InitialInterpolatedLineSearch.h>
#include <RegulaFalsiLineSearch.h>
#include <SecantLineSearch.h>

// constraint handlers
#include <PlainHandler.h>
#include <PenaltyConstraintHandler.h>
//#include <PenaltyHandlerNoHomoSPMultipliers.h>
#include <LagrangeConstraintHandler.h>
#include <TransformationConstraintHandler.h>

// numberers
#include <PlainNumberer.h>
#include <DOF_Numberer.h>

// integrators
#include <LoadControl.h>
#include <ArcLength.h>
#include <ArcLength1.h>
/******************************/
#include <HSConstraint.h>
/******************************/
#include <MinUnbalDispNorm.h>
#include <DisplacementControl.h>
#include <Newmark.h>
#include <WilsonTheta.h>
#include <HHT.h>
#include <HHT1.h>
#include <Newmark1.h> 
#include <EigenIntegrator.h>
#include <CentralDifferenceAlternative.h>
#include <CentralDifferenceNoDamping.h>
#include <CentralDifference.h>
#include <NewmarkExplicit.h>
#include <NewmarkHybridSimulation.h>
#include <HHTExplicit.h>
#include <HHTGeneralized.h>
#include <HHTGeneralizedExplicit.h>
#include <HHTHybridSimulation.h>
#include <AlphaOS.h>
#include <AlphaOSGeneralized.h>
#include <Collocation.h>
#include <CollocationHybridSimulation.h>


// analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>
#include <EigenAnalysis.h>

// system of eqn and solvers
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>

#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>

#include <ConjugateGradientSolver.h>

#include <FullGenLinSOE.h>
#include <FullGenLinLapackSolver.h>

#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <DiagonalSOE.h>
#include <DiagonalDirectSolver.h>

// #include <ProfileSPDLinDirectBlockSolver.h>
// #include <ProfileSPDLinDirectThreadSolver.h>
// #include <ProfileSPDLinDirectSkypackSolver.h>
// #include <BandSPDLinThreadSolver.h>

#include <SparseGenColLinSOE.h>
#ifdef _THREADS
#include <ThreadedSuperLU.h>
#else
#include <SuperLU.h>
#endif


#ifdef _MUMPS
#ifdef _PARALLEL_PROCESSING
#include <MumpsParallelSOE.h>
#include <MumpsParallelSolver.h>
#else
#include <MumpsSOE.h>
#include <MumpsSolver.h>
#endif
#endif

#ifdef _PETSC
#include <PetscSOE.h>
#include <PetscSolver.h>
#include <SparseGenRowLinSOE.h>
#include <PetscSparseSeqSolver.h>
#endif

#include <SymSparseLinSOE.h>
#include <SymSparseLinSolver.h>
#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>

#include <EigenSOE.h>
#include <EigenSolver.h>
#include <SymArpackSOE.h>
#include <SymArpackSolver.h>
#include <BandArpackSOE.h>
#include <BandArpackSolver.h>
#include <SymBandEigenSOE.h>
#include <SymBandEigenSolver.h>


// graph
#include <RCM.h>

#include <ErrorHandler.h>
#include <ConsoleErrorHandler.h>

#ifdef _NOGRAPHICS

#else
#include <TclVideoPlayer.h>
#endif

#include <FE_Datastore.h>

#ifdef _RELIABILITY
// AddingSensitivity:BEGIN /////////////////////////////////////////////////
#include <ReliabilityDomain.h>
#include <SensitivityAlgorithm.h>
#include <SensitivityIntegrator.h>
#include <StaticSensitivityIntegrator.h>
//#include <DynamicSensitivityIntegrator.h>
#include <NewmarkSensitivityIntegrator.h>
#include <RandomVariablePositioner.h>
// AddingSensitivity:END /////////////////////////////////////////////////
#include <TclReliabilityBuilder.h>
static TclReliabilityBuilder *theReliabilityBuilder = 0;
int reliability(ClientData, Tcl_Interp *, int, TCL_Char **);
int wipeReliability(ClientData, Tcl_Interp *, int, TCL_Char **);
#endif

ModelBuilder *theBuilder =0;

// some global variables 
#ifdef _PARALLEL_PROCESSING
#include <DistributedDisplacementControl.h>
#include <ShadowSubdomain.h>
#include <Metis.h>
#include <ShedHeaviest.h>
#include <DomainPartitioner.h>
#include <GraphPartitioner.h>
#include <FEM_ObjectBroker.h>
#include <Subdomain.h>
#include <SubdomainIter.h>
#include <MachineBroker.h>

// parallel analysis
#include <StaticDomainDecompositionAnalysis.h>
#include <TransientDomainDecompositionAnalysis.h>
#include <ParallelNumberer.h>

//  parallel soe & solvers
#include <DistributedBandSPDLinSOE.h>
#include <DistributedSparseGenColLinSOE.h>
#include <DistributedSparseGenRowLinSOE.h>
#include <DistributedBandGenLinSOE.h>
#include <DistributedDiagonalSOE.h>
#include <DistributedDiagonalSolver.h>

#define MPIPP_H
#include <DistributedSuperLU.h>
#include <DistributedProfileSPDLinSOE.h>

PartitionedDomain theDomain;
int OPS_PARALLEL_PROCESSING =0;
int OPS_NUM_SUBDOMAINS      =0;
bool OPS_PARTITIONED        =false;
bool OPS_USING_MAIN_DOMAIN  = false;
int OPS_MAIN_DOMAIN_PARTITION_ID =0;

DomainPartitioner *OPS_DOMAIN_PARTITIONER =0;
GraphPartitioner  *OPS_GRAPH_PARTITIONER =0;
LoadBalancer      *OPS_BALANCER = 0;
FEM_ObjectBroker  *OPS_OBJECT_BROKER;
MachineBroker     *OPS_MACHINE;
Channel          **OPS_theChannels = 0;

#else

Domain theDomain;

#endif

static AnalysisModel *theAnalysisModel =0;
static EquiSolnAlgo *theAlgorithm =0;
static ConstraintHandler *theHandler =0;
static DOF_Numberer *theNumberer =0;
static LinearSOE *theSOE =0;
static StaticAnalysis *theStaticAnalysis = 0;
static DirectIntegrationAnalysis *theTransientAnalysis = 0;
static VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis = 0;

// AddingSensitivity:BEGIN /////////////////////////////////////////////
#ifdef _RELIABILITY
SensitivityAlgorithm *theSensitivityAlgorithm = 0;
static SensitivityIntegrator *theSensitivityIntegrator = 0;
static NewmarkSensitivityIntegrator *theNSI = 0;
#endif
// AddingSensitivity:END ///////////////////////////////////////////////

static StaticIntegrator *theStaticIntegrator =0;
static TransientIntegrator *theTransientIntegrator =0;
static ConvergenceTest *theTest =0;
static bool builtModel = false;

static EigenAnalysis *theEigenAnalysis = 0;

static char *resDataPtr = 0;
static int resDataSize = 0;
static Timer *theTimer = 0;

FE_Datastore *theDatabase  =0;
FEM_ObjectBroker theBroker;

// init the global variabled defined in OPS_Globals.h
double        ops_Dt = 1.0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;

#ifdef _NOGRAPHICS

#else
TclVideoPlayer *theTclVideoPlayer =0;
#endif

// g3AppInit() is the method called by tkAppInit() when the
// interpreter is being set up .. this is where all the
// commands defined in this file are registered with the interpreter.

int 
logFile(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

extern int OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

extern int myCommands(Tcl_Interp *interp);

int g3AppInit(Tcl_Interp *interp) {

#ifndef _LINUX  
    opserr.setFloatField(SCIENTIFIC);
    opserr.setFloatField(FIXEDD);
#endif
    Tcl_CreateObjCommand(interp, "pset", &OPS_SetObjCmd,
			 (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
    Tcl_CreateCommand(interp, "source", &OPS_SourceCmd,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
    Tcl_CreateCommand(interp, "wipe", &wipeModel,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "wipeAnalysis", &wipeAnalysis,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "reset", &resetModel,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);        
    Tcl_CreateCommand(interp, "initialize", &initializeAnalysis,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);        
    Tcl_CreateCommand(interp, "loadConst", &setLoadConst,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);            
    Tcl_CreateCommand(interp, "setTime", &setTime,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);     
    Tcl_CreateCommand(interp, "getTime", &getTime,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);     
    Tcl_CreateCommand(interp, "build", &buildModel,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "analyze", &analyzeModel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "print", &printModel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "analysis", &specifyAnalysis, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "system", &specifySOE, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "numberer", &specifyNumberer, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "constraints", &specifyConstraintHandler, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "algorithm", &specifyAlgorithm, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "test", &specifyCTest, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "integrator", &specifyIntegrator, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "recorder", &addRecorder, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "algorithmRecorder", &addAlgoRecorder, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "database", &addDatabase, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "rigidLink", &rigidLink, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);                
    Tcl_CreateCommand(interp, "rigidDiaphragm", &rigidDiaphragm, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);   
    Tcl_CreateCommand(interp, "eigen", &eigenAnalysis, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "video", &videoPlayer, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "remove", &removeObject, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeDisp", &nodeDisp, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeVel", &nodeVel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeCoord", &nodeCoord, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeBounds", &nodeBounds, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "start", &startTimer, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "stop", &stopTimer, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "rayleigh", &rayleighDamping, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "region", &addRegion, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "logFile", &logFile, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "exit", &OpenSeesExit, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

#ifdef _RELIABILITY
    Tcl_CreateCommand(interp, "wipeReliability", wipeReliability, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
    Tcl_CreateCommand(interp, "reliability", reliability, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
    theReliabilityBuilder = 0;
// AddingSensitivity:BEGIN //////////////////////////////////
    Tcl_CreateCommand(interp, "computeGradients", &computeGradients, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "sensitivityAlgorithm", &sensitivityAlgorithm, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "sensitivityIntegrator", &sensitivityIntegrator, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "sensNodeDisp", &sensNodeDisp, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "sensNodeVel", &sensNodeVel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       

	theSensitivityAlgorithm =0;
	theSensitivityIntegrator =0;
// AddingSensitivity:END //////////////////////////////////

#endif

    theAlgorithm =0;
    theHandler =0;
    theNumberer =0;
    theAnalysisModel =0;  
    theSOE =0;
    theStaticIntegrator =0;
    theTransientIntegrator =0;
    theStaticAnalysis =0;
    theTransientAnalysis =0;    
    theVariableTimeStepTransientAnalysis =0;    
    theTest = 0;

    // create an error handler

#ifdef _NOGRAPHICS

#else
    theTclVideoPlayer = 0;
#endif

    return myCommands(interp);
}

int 
OPS_SetObjCmd(ClientData clientData, Tcl_Interp *interp, int argc, Tcl_Obj * const *argv)
{
  
  if (argc > 2)
    simulationInfo.addParameter(Tcl_GetString(argv[1]), Tcl_GetString(argv[2]));


  Tcl_SetObjCmd(clientData, interp, argc, argv);
  return 0;
}

int 
OPS_SourceCmd(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int ok = TCL_OK;
  if (argc > 1) {
    simulationInfo.addReadFile(argv[1]);

    ok = Tcl_EvalFile(interp, argv[1]);
  }
  return ok;
}



#ifdef _RELIABILITY

int 
reliability(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theReliabilityBuilder == 0) {

    theReliabilityBuilder = new TclReliabilityBuilder(theDomain,interp);
    return TCL_OK;
  }
  else
    return TCL_ERROR;
}

int 
wipeReliability(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theReliabilityBuilder != 0) {
    delete theReliabilityBuilder;
    theReliabilityBuilder = 0;
  }
  return TCL_OK;

}

int 
sensitivityAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	bool withRespectToRVs = true;
	int analysisTypeTag;  // 1: compute at each step wrt. random variables
						  // 2: compute at each step wrt. parameters
						  // 3: compute by command wrt. random variables
						  // 4: compute by command wrt. parameters

	if (argc < 2) {
		opserr << "ERROR: Wrong number of parameters to sensitivity algorithm." << endln;
		return TCL_ERROR;
	}
	if (theReliabilityBuilder == 0) {
		opserr << "The command 'reliability' needs to be issued before " << endln
		<< " the sensitivity algorithm can be created." << endln;
		return TCL_ERROR;
	}
	else if (theSensitivityIntegrator == 0) {
		opserr << "The sensitivity integrator needs to be instantiated before " << endln
		<< " the sensitivity algorithm can be created." << endln;
		return TCL_ERROR;
	}
	else {
		if (argc == 3) {
			if (strcmp(argv[2],"-parameters") == 0) {
				withRespectToRVs = false;
			}
			else if (strcmp(argv[2],"-randomVariables") == 0) {
				withRespectToRVs = true;
			}
			else {
				opserr << "ERROR: Invalid argument in sensitivity algorithm command." << endln;
				return TCL_ERROR;
			}
		}
		if (strcmp(argv[1],"-computeAtEachStep") == 0) {
			if (withRespectToRVs) {
				analysisTypeTag = 1;
			}
			else {
				analysisTypeTag = 2;
			}
		}
		else if (strcmp(argv[1],"-computeByCommand") == 0) {
			if (withRespectToRVs) {
				analysisTypeTag = 3;
			}
			else {
				analysisTypeTag = 4;
			}
		}
		else {
			opserr << "WARNING: Invalid type of sensitivity algorithm." << endln;
			return TCL_ERROR;
		}
		ReliabilityDomain *theReliabilityDomain;
		theReliabilityDomain = theReliabilityBuilder->getReliabilityDomain();
		theSensitivityAlgorithm = new SensitivityAlgorithm(theReliabilityDomain,
									theAlgorithm,
									theSensitivityIntegrator,
									analysisTypeTag);
		if (theSensitivityAlgorithm == 0) {
			opserr << "ERROR: Could not create theSensitivityAlgorithm. " << endln;
			return TCL_ERROR;
		}
	}
	return TCL_OK;
}
int 
sensitivityIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	if (strcmp(argv[1],"-static") == 0) {

		if (theAnalysisModel == 0) {
			theAnalysisModel = new AnalysisModel();
		}
		theSensitivityIntegrator = new StaticSensitivityIntegrator(theAnalysisModel, theSOE);
		return TCL_OK;
	}
	else if (strcmp(argv[1],"-definedAbove") == 0) {  

		if (theNSI == 0) {
			opserr << "ERROR: No sensitivity integrator has been specified. " << endln;
			return TCL_ERROR;
		}
		else {
			theSensitivityIntegrator = theNSI;
			return TCL_OK;
		}
	}
//	else if (strcmp(argv[1],"Dynamic") == 0) {  
//
//		if (theAnalysisModel == 0) {
//			theAnalysisModel = new AnalysisModel();
//		}
//
//		theSensitivityIntegrator = new DynamicSensitivityIntegrator(theAnalysisModel, theSOE, theTransientIntegrator);
//		return TCL_OK;
//	}
	else {
		opserr << "WARNING: Invalid type of sensitivity integrator." << endln;
		return TCL_ERROR;
	}
}
// AddingSensitivity:END /////////////////////////////////////////////////

#endif





int 
wipeModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // to build the model make sure the ModelBuilder has been constructed
  // and that the model has not already been constructed
  if (theBuilder != 0) {
    delete theBuilder;
    builtModel = false;
    theBuilder = 0;
  }

  if (theStaticAnalysis != 0) {
      theStaticAnalysis->clearAll();
      delete theStaticAnalysis;
  }
  
  if (theTransientAnalysis != 0) {
      theTransientAnalysis->clearAll();
      delete theTransientAnalysis;  
  }


  // NOTE : DON'T do the above on theVariableTimeStepAnalysis
  // as it and theTansientAnalysis are one in the same

  /*
  if (theEigenAnalysis != 0) {
    delete theEigenAnalysis;
    theEigenAnalysis = 0;
  }
  */

  if (theDatabase != 0)
    delete theDatabase;
  
  theDomain.clearAll();

#ifdef _NOGRAPHICS

#else

  if (theTclVideoPlayer != 0) {
	  delete theTclVideoPlayer;
	  theTclVideoPlayer = 0;
  }
#endif

  theAlgorithm =0;
  theHandler =0;
  theNumberer =0;
  theAnalysisModel =0;  
  theSOE =0;
  theStaticIntegrator =0;
  theTransientIntegrator =0;
  theStaticAnalysis =0;
  theTransientAnalysis =0;    
  theVariableTimeStepTransientAnalysis =0;    

  theTest = 0;
  theDatabase = 0;

// AddingSensitivity:BEGIN /////////////////////////////////////////////////
#ifdef _RELIABILITY
  theSensitivityAlgorithm =0;
  theSensitivityIntegrator =0;
#endif
// AddingSensitivity:END /////////////////////////////////////////////////

  // the domain deletes the record objects, 
  // just have to delete the private array

  return TCL_OK;  
}

int 
wipeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

#ifdef _PARALLEL_PROCESSING
  if (OPS_PARTITIONED == true && OPS_NUM_SUBDOMAINS > 1) {
    SubdomainIter &theSubdomains = theDomain.getSubdomains();
    Subdomain *theSub =0;
    
    // create the appropriate domain decomposition analysis
    while ((theSub = theSubdomains()) != 0) 
      theSub->wipeAnalysis();
  }
#endif

  if (theStaticAnalysis != 0) {
      theStaticAnalysis->clearAll();
      delete theStaticAnalysis;
  }
  
  if (theTransientAnalysis != 0) {
      theTransientAnalysis->clearAll();
      delete theTransientAnalysis;  
  }

  // NOTE : DON'T do the above on theVariableTimeStepAnalysis
  // as it and theTansientAnalysis are one in the same

  theAlgorithm =0;
  theHandler =0;
  theNumberer =0;
  theAnalysisModel =0;  
  theSOE =0;
  theStaticIntegrator =0;
  theTransientIntegrator =0;
  theStaticAnalysis =0;
  theTransientAnalysis =0;    
  theVariableTimeStepTransientAnalysis =0;    
  theTest = 0;

// AddingSensitivity:BEGIN /////////////////////////////////////////////////
#ifdef _RELIABILITY
  theSensitivityAlgorithm =0;
  theSensitivityIntegrator =0;
#endif
// AddingSensitivity:END /////////////////////////////////////////////////
  // the domain deletes the record objects, 
  // just have to delete the private array

  return TCL_OK;  
}


int 
resetModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	theDomain.revertToStart();


// AddingSensitivity:BEGIN ////////////////////////////////////
	// (Calling domainChanged() don't work because the ::getX()
	// of the linear SOE is being called without the X having
	// been instantitated.)
	if (theTransientIntegrator != 0) {
		theTransientIntegrator->revertToStart();
	}
// AddingSensitivity:END //////////////////////////////////////

	
	return TCL_OK;
}

int
initializeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theTransientAnalysis != 0)
    theTransientAnalysis->initialize();
  else if (theStaticAnalysis != 0)
    theStaticAnalysis->initialize();
  
  theDomain.initialize();

  return TCL_OK;
}


int 
setLoadConst(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  theDomain.setLoadConstant();
  if (argc == 3) {
      if( strcmp(argv[1],"-time") == 0) {
	  double newTime;
	  if (Tcl_GetDouble(interp, argv[2], &newTime) != TCL_OK) {
	      opserr << "WARNING readingvalue - loadConst -time value \n";
	      return TCL_ERROR;
	  } else {
	      theDomain.setCurrentTime(newTime);
	      theDomain.setCommittedTime(newTime);
	  }
      }    	  
  }
	  
  return TCL_OK;
}


int 
setTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
      opserr << "WARNING illegal command - time pseudoTime? \n";
      return TCL_ERROR;
  }
  double newTime;
  if (Tcl_GetDouble(interp, argv[1], &newTime) != TCL_OK) {
      opserr << "WARNING reading time value - time pseudoTime? \n";
      return TCL_ERROR;
  } else {
      theDomain.setCurrentTime(newTime);
      theDomain.setCommittedTime(newTime);
  }
  return TCL_OK;
}

int 
getTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  double time = theDomain.getCurrentTime();
  // now we copy the value to the tcl string that is returned
  sprintf(interp->result,"%f",time);
  return TCL_OK;
}


// command invoked to build the model, i.e. to invoke buildFE_Model() 
// on the ModelBuilder
int 
buildModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // to build the model make sure the ModelBuilder has been constructed
  // and that the model has not already been constructed
  if (theBuilder != 0 && builtModel == false) {
    builtModel = true;
    return theBuilder->buildFE_Model();
  }  else if (theBuilder != 0 && builtModel == true) {
      opserr << "WARNING Model has already been built - not built again \n";
      return TCL_ERROR;
  }
  else {
      opserr << "WARNING No ModelBuilder type has been specified \n";
      return TCL_ERROR;
  }    
}


#ifdef _PARALLEL_PROCESSING

int 
partitionModel(void)
{
  int result = 0;
  
  if (OPS_theChannels != 0)
    delete [] OPS_theChannels;

  OPS_theChannels = new Channel *[OPS_NUM_SUBDOMAINS];
  
  // create some subdomains
  for (int i=1; i<=OPS_NUM_SUBDOMAINS; i++) {
    if (i != OPS_MAIN_DOMAIN_PARTITION_ID) {
      opserr << "Creating SUBDOMAIN: " << i << endln;
      ShadowSubdomain *theSubdomain = new ShadowSubdomain(i, *OPS_MACHINE, *OPS_OBJECT_BROKER);
      theDomain.addSubdomain(theSubdomain);
      OPS_theChannels[i-1] = theSubdomain->getChannelPtr();
    }
  }

  // create a partitioner & partition the domain
  if (OPS_DOMAIN_PARTITIONER == 0) {
    //      OPS_BALANCER = new ShedHeaviest();
    OPS_GRAPH_PARTITIONER  = new Metis;
    //OPS_DOMAIN_PARTITIONER = new DomainPartitioner(*OPS_GRAPH_PARTITIONER, *OPS_BALANCER);
    OPS_DOMAIN_PARTITIONER = new DomainPartitioner(*OPS_GRAPH_PARTITIONER);
    theDomain.setPartitioner(OPS_DOMAIN_PARTITIONER);
  }

  theDomain.partition(OPS_NUM_SUBDOMAINS, OPS_USING_MAIN_DOMAIN, OPS_MAIN_DOMAIN_PARTITION_ID);
  
  OPS_PARTITIONED = true;
  
  DomainDecompositionAnalysis *theSubAnalysis;
  SubdomainIter &theSubdomains = theDomain.getSubdomains();
  Subdomain *theSub =0;
  
  // create the appropriate domain decomposition analysis
  while ((theSub = theSubdomains()) != 0) {
    if (theStaticAnalysis != 0) {      
      theSubAnalysis = new StaticDomainDecompositionAnalysis(*theSub,
							     *theHandler,
							     *theNumberer,
							     *theAnalysisModel,
							     *theAlgorithm,
							     *theSOE,
							     *theStaticIntegrator,
							     theTest,
							     false);
      
    } else {
      theSubAnalysis = new TransientDomainDecompositionAnalysis(*theSub,
								*theHandler,
								*theNumberer,
								*theAnalysisModel,
								*theAlgorithm,
								*theSOE,
								*theTransientIntegrator,
								theTest,
								false);
    }       
    theSub->setDomainDecompAnalysis(*theSubAnalysis);
    //  delete theSubAnalysis;
  }
  return result;
}

#endif






//
// command invoked to build the model, i.e. to invoke analyze() 
// on the Analysis object
//
int 
analyzeModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int result = 0;

#ifdef _PARALLEL_PROCESSING
  if (OPS_PARTITIONED == false && OPS_NUM_SUBDOMAINS > 1) 
    if (partitionModel() < 0) {
      opserr << "WARNING before analysis; partition failed\n";
      return TCL_ERROR;
    }
#endif

  if (theStaticAnalysis != 0) {
    if (argc < 2) {
      opserr << "WARNING static analysis: analysis numIncr?\n";
      return TCL_ERROR;
    }
    int numIncr;

    if (Tcl_GetInt(interp, argv[1], &numIncr) != TCL_OK)	
      return TCL_ERROR;	      
    result = theStaticAnalysis->analyze(numIncr);

  } else if (theTransientAnalysis != 0) {
    if (argc < 3) {
      opserr << "WARNING transient analysis: analysis numIncr? deltaT?\n";
      return TCL_ERROR;
    }
    int numIncr;
    if (Tcl_GetInt(interp, argv[1], &numIncr) != TCL_OK)	
      return TCL_ERROR;
    double dT;
    if (Tcl_GetDouble(interp, argv[2], &dT) != TCL_OK)	
      return TCL_ERROR;

    // Set global timestep variable
    ops_Dt = dT;

    if (argc == 6) {
      int Jd;
      double dtMin, dtMax;
      if (Tcl_GetDouble(interp, argv[3], &dtMin) != TCL_OK)	
	return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[4], &dtMax) != TCL_OK)	
	return TCL_ERROR;
      if (Tcl_GetInt(interp, argv[5], &Jd) != TCL_OK)	
	return TCL_ERROR;

      if (theVariableTimeStepTransientAnalysis != 0)
	result =  theVariableTimeStepTransientAnalysis->analyze(numIncr, dT, dtMin, dtMax, Jd);
      else {
	opserr << "WARNING analyze - no variable time step transient analysis object constructed\n";
	return TCL_ERROR;
      }

    } else {
      result = theTransientAnalysis->analyze(numIncr, dT);
    }

  } else {
    opserr << "WARNING No Analysis type has been specified \n";
    return TCL_ERROR;
  }    

  if (result < 0) {
    opserr << "OpenSees > analyze failed, returned: " << result << " error flag\n";
  }

  sprintf(interp->result,"%d",result);    

  return TCL_OK;

}

int 
printElement(ClientData clientData, Tcl_Interp *interp, int argc, 
	     TCL_Char **argv, int nodeArg, OPS_Stream &output);


int 
printNode(ClientData clientData, Tcl_Interp *interp, int argc, 
	  TCL_Char **argv, int nodeArg, OPS_Stream &output);
	  
int 
printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, 
		TCL_Char **argv, int nodeArg, OPS_Stream &output);	  
		
int 
printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, 
	       TCL_Char **argv, int nodeArg, OPS_Stream &output);	  		


int 
printModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // if just 'print' then print out the entire domain
  if (argc == 1) {
    opserr << theDomain;
    return TCL_OK;
  }    

  // if 'print ele i j k..' print out some elements
  if ((strcmp(argv[1],"-ele") == 0) || (strcmp(argv[1],"ele") == 0))
    return printElement(clientData, interp, argc, argv, 3, opserr);    

  // if 'print node i j k ..' print out some nodes
  else if ((strcmp(argv[1],"-node") == 0) || (strcmp(argv[1],"node") == 0)) 
      return printNode(clientData, interp, argc, argv, 3, opserr);
  
  // if 'print integrator flag' print out the integrator
  else if ((strcmp(argv[1],"integrator") == 0) || 
	   (strcmp(argv[1],"-integrator") == 0)) 
    return printIntegrator(clientData, interp, argc, argv, 3, opserr);  
  
  // if 'print algorithm flag' print out the algorithm
  else if ((strcmp(argv[1],"algorithm") == 0) || 
	   (strcmp(argv[1],"-algorithm") == 0))
    return printAlgorithm(clientData, interp, argc, argv, 3, opserr);    

  else { // it must be a file we are going to print to
    
#ifdef _USING_OpenSees_STREAMS
    FileStream output;
    if (output.setFile(argv[1], APPEND) != 0) {
      opserr << "print <filename> .. - failed to open file: " << argv[1] << endln;
      return TCL_ERROR;
    }
#else
    ofstream output(argv[1],ios::app); // open for appending to
    if (!output) {
      opserr << "print <filename> .. - failed to open file: " << argv[1] << endln;
      return TCL_ERROR;
    }
#endif

    // if just 'print <filename>' then print out the entire domain to eof
    if (argc == 2) {
      output << theDomain;
      return TCL_OK;
    }    

    int pos = 2;
    if ((strcmp(argv[pos],"string") == 0) || 
	(strcmp(argv[pos],"-string") == 0)) {
	output << argv[3] << endln;
	pos +=2;
    }
    int res = TCL_OK;    

    // if 'print <filename> ele i j k..' print out some elements
    if ((strcmp(argv[pos],"ele") == 0) || 
	(strcmp(argv[pos],"-ele") == 0))
      res = printElement(clientData, interp, argc, argv, pos+2, output);    

    // if 'print <filename> node i j k ..' print out some nodes
    else if ((strcmp(argv[pos],"node") == 0) || (strcmp(argv[pos],"-node") == 0))
      res = printNode(clientData, interp, argc, argv, pos+2, output);
    
    // if 'print integrator flag' print out the integrator
    else if ((strcmp(argv[pos],"integrator") == 0) 
	     || (strcmp(argv[pos],"-integrator") == 0))
	return printIntegrator(clientData, interp, argc, argv, pos+2, opserr);  
  
    // if 'print algorithm flag' print out the algorithm
    else if ((strcmp(argv[pos],"-algorithm") == 0)|| 
	     (strcmp(argv[pos],"algorithm") == 0))
	return printAlgorithm(clientData, interp, argc, argv, pos+2, opserr);    


    // close the output file
    output.close();
    return res;
  }
  
}


// printNode():
// function to print out the nodal information conatined in line
//     print <filename> node <flag int> <int int int>
// input: nodeArg: integer equal to arg count to node plus 1
//        output: output stream to which the results are sent
// 
int 
printNode(ClientData clientData, Tcl_Interp *interp, int argc, 
	  TCL_Char **argv, int nodeArg, OPS_Stream &output)
{
  int flag = 0; // default flag sent to a nodes Print() method

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc < nodeArg) { 
      NodeIter &theNodes = theDomain.getNodes();
      Node *theNode;
      while ((theNode = theNodes()) != 0)
	theNode->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> node flag int <int int ..>' get the flag
  if ((strcmp(argv[nodeArg-1],"flag") == 0) ||
      (strcmp(argv[nodeArg-1],"-flag") == 0)) { 
      // get the specified flag
    if (argc <= nodeArg) {
      opserr << "WARNING print <filename> node <flag int> no int specified \n";
      return TCL_ERROR;
    }    
    if (Tcl_GetInt(interp, argv[nodeArg], &flag) != TCL_OK) {
      opserr << "WARNING print node failed to get integer flag: \n";
      opserr << argv[nodeArg] << endln; 
      return TCL_ERROR;
    }    
    nodeArg += 2;
  }

  // now print the nodes with the specified flag, 0 by default

  // if 'print <filename> node flag' 
  //     print out all the nodes in the domain with flag
  if (argc < nodeArg) { 
    NodeIter &theNodes = theDomain.getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0)
      theNode->Print(output, flag);
    return TCL_OK;
  } else { 
    // otherwise print out the specified nodes i j k .. with flag
    for (int i= nodeArg-1; i<argc; i++) {
      int nodeTag;
      if (Tcl_GetInt(interp, argv[i], &nodeTag) != TCL_OK) {
	opserr << "WARNING print node failed to get integer: " << argv[i] << endln;
	return TCL_ERROR;
      }
      Node *theNode = theDomain.getNode(nodeTag);
      if (theNode != 0)
	theNode->Print(output,flag);
    }
    return TCL_OK;
  }
}


int 
printElement(ClientData clientData, Tcl_Interp *interp, int argc, 
	  TCL_Char **argv, int eleArg, OPS_Stream &output)
{
  int flag = 0; // default flag sent to a nodes Print() method

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc < eleArg) { 
      ElementIter &theElements = theDomain.getElements();
      Element *theElement;
      while ((theElement = theElements()) != 0)
	theElement->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> Element flag int <int int ..>' get the flag
  if ((strcmp(argv[eleArg-1],"flag") == 0) ||
      (strcmp(argv[eleArg-1],"-flag")) == 0) { // get the specified flag
    if (argc <= eleArg) {
      opserr << "WARNING print <filename> ele <flag int> no int specified \n";
      return TCL_ERROR;
    }    
    if (Tcl_GetInt(interp, argv[eleArg], &flag) != TCL_OK) {
      opserr << "WARNING print ele failed to get integer flag: \n";
      opserr << argv[eleArg] << endln; 
      return TCL_ERROR;
    }    
    eleArg += 2;
  }

  // now print the Elements with the specified flag, 0 by default

  // if 'print <filename> Element flag' 
  //     print out all the Elements in the domain with flag
  if (flag == 2)
    output << "#FRAME\n";
  
  if (argc < eleArg) { 
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)      
      theElement->Print(output, flag);
    return TCL_OK;
  } else { 
    // otherwise print out the specified Elements i j k .. with flag
    for (int i= eleArg-1; i<argc; i++) {
      int ElementTag;
      if (Tcl_GetInt(interp, argv[i], &ElementTag) != TCL_OK) {
	opserr << "WARNING print ele failed to get integer: " << argv[i] << endln;
	return TCL_ERROR;
      }
      Element *theElement = theDomain.getElement(ElementTag);
      if (theElement != 0)
	theElement->Print(output,flag);
    }
    return TCL_OK;
  }
}


int 
printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, 
	       TCL_Char **argv, int eleArg, OPS_Stream &output)
{

  if (theAlgorithm == 0)
      return TCL_OK;

  // if just 'print <filename> algorithm'- no flag
  if (argc < eleArg) { 
      theAlgorithm->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> Algorithm flag' get the flag
  int flag;  
  if (Tcl_GetInt(interp, argv[eleArg-1], &flag) != TCL_OK) {  
      opserr << "WARNING print algorithm failed to get integer flag: \n";
      opserr << argv[eleArg] << endln; 
      return TCL_ERROR;
  }    
  theAlgorithm->Print(output,flag);
  return TCL_OK;  
}


int 
printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, 
		TCL_Char **argv, int eleArg, OPS_Stream &output)
{

  if (theStaticIntegrator == 0 && theTransientIntegrator == 0)
      return TCL_OK;
  
  IncrementalIntegrator *theIntegrator;
  if (theStaticIntegrator != 0)
      theIntegrator = theStaticIntegrator;
  else
      theIntegrator = theTransientIntegrator;

  // if just 'print <filename> algorithm'- no flag
  if (argc < eleArg) { 
      theIntegrator->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> Algorithm flag' get the flag
  int flag;  
  if (Tcl_GetInt(interp, argv[eleArg-1], &flag) != TCL_OK) {  
      opserr << "WARNING print algorithm failed to get integer flag: \n";
      opserr << argv[eleArg] << endln; 
      return TCL_ERROR;
  }    
  theIntegrator->Print(output,flag);
  return TCL_OK;  
}


//
// command invoked to allow the Analysis object to be built
//
int 
specifyAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, 
		TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING need to specify an analysis type (Static, Transient)\n";
	return TCL_ERROR;
    }    

    // delete the old analysis
    if (theStaticAnalysis != 0) {
	delete theStaticAnalysis;
	theStaticAnalysis = 0;
    }
    if (theTransientAnalysis != 0) {
	delete theTransientAnalysis;
	theTransientAnalysis = 0;
	theVariableTimeStepTransientAnalysis = 0;
    }
    
    // check argv[1] for type of SOE and create it
    if (strcmp(argv[1],"Static") == 0) {
	// make sure all the components have been built,
	// otherwise print a warning and use some defaults
	if (theAnalysisModel == 0) 
	    theAnalysisModel = new AnalysisModel();

	if (theTest == 0) 
	  theTest = new CTestNormUnbalance(1.0e-6,25,0);       
	
	if (theAlgorithm == 0) {
	    opserr << "WARNING analysis Static - no Algorithm yet specified, \n";
	    opserr << " NewtonRaphson default will be used\n";	    

	    theAlgorithm = new NewtonRaphson(*theTest); 
	}
	if (theHandler == 0) {
	    opserr << "WARNING analysis Static - no ConstraintHandler yet specified, \n";
	    opserr << " PlainHandler default will be used\n";
	    theHandler = new PlainHandler();       
	}
	if (theNumberer == 0) {
	    opserr << "WARNING analysis Static - no Numberer specified, \n";
	    opserr << " RCM default will be used\n";
	    RCM *theRCM = new RCM(false);	
	    theNumberer = new DOF_Numberer(*theRCM);    	
	}
	if (theStaticIntegrator == 0) {
	    opserr << "WARNING analysis Static - no Integrator specified, \n";
	    opserr << " StaticIntegrator default will be used\n";
	    theStaticIntegrator = new LoadControl(1, 1, 1, 1);       
	}
	if (theSOE == 0) {
	    opserr << "WARNING analysis Static - no LinearSOE specified, \n";
	    opserr << " ProfileSPDLinSOE default will be used\n";
	    ProfileSPDLinSolver *theSolver;
	    theSolver = new ProfileSPDLinDirectSolver(); 	
	    theSOE = new ProfileSPDLinSOE(*theSolver);      
	}
    
	theStaticAnalysis = new StaticAnalysis(theDomain,
					       *theHandler,
					       *theNumberer,
					       *theAnalysisModel,
					       *theAlgorithm,
					       *theSOE,
					       *theStaticIntegrator,
					       theTest);

// AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
	if (theSensitivityAlgorithm != 0 && theSensitivityAlgorithm->shouldComputeAtEachStep()) {
		theStaticAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
	}
#endif
// AddingSensitivity:END /////////////////////////////////

    } else if (strcmp(argv[1],"Transient") == 0) {
	// make sure all the components have been built,
	// otherwise print a warning and use some defaults
	if (theAnalysisModel == 0) 
	    theAnalysisModel = new AnalysisModel();

	if (theTest == 0) 
	  theTest = new CTestNormUnbalance(1.0e-6,25,0);       
	
	if (theAlgorithm == 0) {
	    opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
	    opserr << " NewtonRaphson default will be used\n";	    

	    theAlgorithm = new NewtonRaphson(*theTest); 
	}
	if (theHandler == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
	    opserr << " yet specified, PlainHandler default will be used\n";
	    theHandler = new PlainHandler();       
	}
	if (theNumberer == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no Numberer specified, \n";
	    opserr << " RCM default will be used\n";
	    RCM *theRCM = new RCM(false);	
	    theNumberer = new DOF_Numberer(*theRCM);    	
	}
	if (theTransientIntegrator == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no Integrator specified, \n";
	    opserr << " Newmark(.5,.25) default will be used\n";
	    theTransientIntegrator = new Newmark(0.5,0.25);       
	}
	if (theSOE == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no LinearSOE specified, \n";
	    opserr << " ProfileSPDLinSOE default will be used\n";
	    ProfileSPDLinSolver *theSolver;
	    theSolver = new ProfileSPDLinDirectSolver(); 	
	    theSOE = new ProfileSPDLinSOE(*theSolver);      
	}
    
	theTransientAnalysis = new DirectIntegrationAnalysis(theDomain,
							     *theHandler,
							     *theNumberer,
							     *theAnalysisModel,
							     *theAlgorithm,
							     *theSOE,
							     *theTransientIntegrator,
							     theTest);


// AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
	if (theSensitivityAlgorithm != 0 && theSensitivityAlgorithm->shouldComputeAtEachStep()) {
		theTransientAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
	}
#endif
// AddingSensitivity:END /////////////////////////////////

    } else if ((strcmp(argv[1],"VariableTimeStepTransient") == 0) ||
	       (strcmp(argv[1],"TransientWithVariableTimeStep") == 0) ||
	       (strcmp(argv[1],"VariableTransient") == 0)) {
	// make sure all the components have been built,
	// otherwise print a warning and use some defaults
	if (theAnalysisModel == 0) 
	    theAnalysisModel = new AnalysisModel();

	if (theTest == 0) 
	  theTest = new CTestNormUnbalance(1.0e-6,25,0);       
	
	if (theAlgorithm == 0) {
	    opserr << "WARNING analysis Transient - no Algorithm yet specified, \n";
	    opserr << " NewtonRaphson default will be used\n";	    
	    theAlgorithm = new NewtonRaphson(*theTest); 
	}

	if (theHandler == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no ConstraintHandler\n";
	    opserr << " yet specified, PlainHandler default will be used\n";
	    theHandler = new PlainHandler();       
	}

	if (theNumberer == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no Numberer specified, \n";
	    opserr << " RCM default will be used\n";
	    RCM *theRCM = new RCM(false);	
	    theNumberer = new DOF_Numberer(*theRCM);    	
	}

	if (theTransientIntegrator == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no Integrator specified, \n";
	    opserr << " Newmark(.5,.25) default will be used\n";
	    theTransientIntegrator = new Newmark(0.5,0.25);       
	}

	if (theSOE == 0) {
	    opserr << "WARNING analysis Transient dt tFinal - no LinearSOE specified, \n";
	    opserr << " ProfileSPDLinSOE default will be used\n";
	    ProfileSPDLinSolver *theSolver;
	    theSolver = new ProfileSPDLinDirectSolver(); 	
	    theSOE = new ProfileSPDLinSOE(*theSolver);      
	}
    
	theVariableTimeStepTransientAnalysis = new VariableTimeStepDirectIntegrationAnalysis
	  (theDomain,
	   *theHandler,
	   *theNumberer,
	   *theAnalysisModel,
	   *theAlgorithm,
	   *theSOE,
	   *theTransientIntegrator,
	   theTest);

	// set the pointer for variabble time step analysis
	theTransientAnalysis = theVariableTimeStepTransientAnalysis;
	
    } else {
	opserr << "WARNING No Analysis type exists (Static Transient only) \n";
	return TCL_ERROR;
    }

#ifdef _PARALLEL_PROCESSING
    if (OPS_PARTITIONED == true && OPS_NUM_SUBDOMAINS > 1) {
      DomainDecompositionAnalysis *theSubAnalysis;
      SubdomainIter &theSubdomains = theDomain.getSubdomains();
      Subdomain *theSub =0;
      // create the appropriate domain decomposition analysis
      while ((theSub = theSubdomains()) != 0) {
	if (theStaticAnalysis != 0) {      
	  theSubAnalysis = new StaticDomainDecompositionAnalysis(*theSub,
								 *theHandler,
								 *theNumberer,
								 *theAnalysisModel,
								 *theAlgorithm,
								 *theSOE,
								 *theStaticIntegrator,
								 theTest,
								 false);
	  
	} else {
	  theSubAnalysis = new TransientDomainDecompositionAnalysis(*theSub,
								    *theHandler,
								    *theNumberer,
								    *theAnalysisModel,
								    *theAlgorithm,
								    *theSOE,
								    *theTransientIntegrator,
								    theTest,
								  false);
	}       
	
	theSub->setDomainDecompAnalysis(*theSubAnalysis);
	//	delete theSubAnalysis;
      }
    }
#endif
    return TCL_OK;
}


//
// command invoked to allow the SystemOfEqn and Solver objects to be built
//


typedef struct linearSOE_PackageCommand {
  char *funcName;
  int (*funcPtr)(ClientData clientData, 
		 Tcl_Interp *interp,  
		 int argc, 
		 TCL_Char **argv, 
		 FEM_ObjectBroker *,
		 LinearSOE **); 
  struct linearSOE_PackageCommand *next;
} LinearSOE_PackageCommand;


// static variables
static LinearSOE_PackageCommand *theLinearSOE_PackageCommands = NULL;

int 
specifySOE(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
      opserr << "WARNING need to specify a model type \n";
      return TCL_ERROR;
  }    



  // check argv[1] for type of SOE and create it
  // BAND GENERAL SOE & SOLVER
  if ((strcmp(argv[1],"BandGeneral") == 0) || (strcmp(argv[1],"BandGEN") == 0)
      || (strcmp(argv[1],"BandGen") == 0)){
    BandGenLinSolver    *theSolver = new BandGenLinLapackSolver();
#ifdef _PARALLEL_PROCESSING
    theSOE = new DistributedBandGenLinSOE(*theSolver);      
#else
    theSOE = new BandGenLinSOE(*theSolver);      
#endif
  } 



  // BAND SPD SOE & SOLVER
  else if (strcmp(argv[1],"BandSPD") == 0) {
      BandSPDLinSolver    *theSolver = new BandSPDLinLapackSolver();   
#ifdef _PARALLEL_PROCESSING
      theSOE = new DistributedBandSPDLinSOE(*theSolver);        
#else
      theSOE = new BandSPDLinSOE(*theSolver);        
#endif

  } 

  // Diagonal SOE & SOLVER
  else if (strcmp(argv[1],"Diagonal") == 0) {
#ifdef _PARALLEL_PROCESSING
      DistributedDiagonalSolver    *theSolver = new DistributedDiagonalSolver();   
      theSOE = new DistributedDiagonalSOE(*theSolver);
#else
      DiagonalSolver    *theSolver = new DiagonalDirectSolver();   
      theSOE = new DiagonalSOE(*theSolver);
#endif


  } 

  // PROFILE SPD SOE * SOLVER
  else if (strcmp(argv[1],"ProfileSPD") == 0) {
    // now must determine the type of solver to create from rest of args
    ProfileSPDLinSolver *theSolver = new ProfileSPDLinDirectSolver(); 	

    /* *********** Some misc solvers i play with ******************
    else if (strcmp(argv[2],"Normal") == 0) {
      theSolver = new ProfileSPDLinDirectSolver(); 	
    } 

    else if (strcmp(argv[2],"Block") == 0) {  
      int blockSize = 4;
      if (argc == 4) {
	if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
	  return TCL_ERROR;
      }
      theSolver = theSolver = new ProfileSPDLinDirectBlockSolver(1.0e-12,blockSize); 
    }

    
      int blockSize = 4;
      int numThreads = 1;
      if (argc == 5) {
	if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[4], &numThreads) != TCL_OK)
	  return TCL_ERROR;
      }
      theSolver = new ProfileSPDLinDirectThreadSolver(numThreads,blockSize,1.0e-12); 
      } else if (strcmp(argv[2],"Thread") == 0) {  
      int blockSize = 4;
      int numThreads = 1;
      if (argc == 5) {
	if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[4], &numThreads) != TCL_OK)
	  return TCL_ERROR;
      }
      theSolver = new ProfileSPDLinDirectThreadSolver(numThreads,blockSize,1.0e-12); 
    } 
    else if (strcmp(argv[2],"Skypack") == 0) {  
      if (argc == 5) {
	int mCols, mRows;
	if (Tcl_GetInt(interp, argv[3], &mCols) != TCL_OK)
	  return TCL_ERROR;
	if (Tcl_GetInt(interp, argv[4], &mRows) != TCL_OK)
	  return TCL_ERROR;
	theSolver = new ProfileSPDLinDirectSkypackSolver(mCols, mRows); 
      } else 
	theSolver = new ProfileSPDLinDirectSkypackSolver(); 	
    }
    else 
      theSolver = new ProfileSPDLinDirectSolver(); 	
    ***************************************************************  */

#ifdef _PARALLEL_PROCESSING
    theSOE = new DistributedProfileSPDLinSOE(*theSolver);
#else
    theSOE = new ProfileSPDLinSOE(*theSolver);      
#endif
  }

  // SPARSE GENERAL SOE * SOLVER
  else if ((strcmp(argv[1],"SparseGeneral") == 0) || (strcmp(argv[1],"SuperLU") == 0) ||
	   (strcmp(argv[1],"SparseGEN") == 0)) {
    
    SparseGenColLinSolver *theSolver =0;    

    double thresh = 0.0;
    int npRow = 1;
    int npCol = 1;
    int np = 1;

    // defaults for threaded SuperLU
    int count = 2;
    int permSpec = 0;
    int panelSize = 6;
    int relax = 6;




    while (count < argc) {

      if ((strcmp(argv[count],"p") == 0) || (strcmp(argv[count],"piv") == 0)||
	  (strcmp(argv[count],"-piv") == 0)) {
	thresh = 1.0;
      }
      else if ((strcmp(argv[count],"-np") == 0) || (strcmp(argv[count],"np") == 0)) {
	count++;
	if (count < argc)
	  if (Tcl_GetInt(interp, argv[count], &np) != TCL_OK)
	    return TCL_ERROR;		     
      }
      else if ((strcmp(argv[count],"npRow") == 0) || (strcmp(argv[count],"-npRow") ==0)) {
	count++;
	if (count < argc)
	  if (Tcl_GetInt(interp, argv[count], &npRow) != TCL_OK)
	    return TCL_ERROR;		     
      } else if ((strcmp(argv[count],"npCol") == 0) || (strcmp(argv[count],"-npCol") ==0)) {
	count++;
	if (count < argc)
	  if (Tcl_GetInt(interp, argv[count], &npCol) != TCL_OK)
	    return TCL_ERROR;		     
      }
      count++;
    }
  

#ifdef _THREADS
    if (np != 0)
      theSolver = new ThreadedSuperLU(np, permSpec, panelSize, relax, thresh); 	
#endif

#ifdef _PARALLEL_PROCESSING
    if (theSolver != 0)
      delete theSolver;
    theSolver = 0;

    if (npRow != 0 && npCol != 0) {
      theSolver = new DistributedSuperLU(npRow, npCol);
    }
#else

    char symmetric = 'N';
    double drop_tol = 0.0;

    while (count < argc) {
      if (strcmp(argv[count],"s") == 0 || strcmp(argv[count],"symmetric") ||
	  strcmp(argv[count],"-symm")) {
	symmetric = 'Y';
      }
      count++;
    }
    theSolver = new SuperLU(permSpec, drop_tol, panelSize, relax, symmetric); 	

#endif

#ifdef _PARALLEL_PROCESSING
    theSOE = new DistributedSparseGenColLinSOE(*theSolver);      
#else
    theSOE = new SparseGenColLinSOE(*theSolver);      
#endif
  }

  
  else if (strcmp(argv[1],"SparseSPD") == 0) {
    // now must determine the type of solver to create from rest of args

    // now determine ordering scheme
    //   1 -- MMD
    //   2 -- ND
    //   3 -- RCM
    int lSparse = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[3], &lSparse) != TCL_OK)
	return TCL_ERROR;
    }

    SymSparseLinSolver *theSolver = new SymSparseLinSolver();
    theSOE = new SymSparseLinSOE(*theSolver, lSparse);      
  }	
  
  else if (strcmp(argv[1],"UmfPack") == 0) {
    // now must determine the type of solver to create from rest of args
      UmfpackGenLinSolver *theSolver = new UmfpackGenLinSolver();
      theSOE = new UmfpackGenLinSOE(*theSolver);      
  }	  

  else if (strcmp(argv[1],"FullGeneral") == 0) {
    // now must determine the type of solver to create from rest of args
    FullGenLinLapackSolver *theSolver = new FullGenLinLapackSolver();
    theSOE = new FullGenLinSOE(*theSolver);
  }

#ifdef _PETSC

  else if (strcmp(argv[1],"Petsc") == 0) {
    // now must determine the type of solver to create from rest of args
    KSPType method = KSPCG;            // KSPCG KSPGMRES
    PCType preconditioner = PCJACOBI; // PCJACOBI PCILU PCBJACOBI
    int matType = 0;
    
    double rTol = 1.0e-5;
    double aTol = 1.0e-50;
    double dTol = 1.0e5;
    int maxIts = 100000;
    int count = 2;
    while (count < argc-1) {
      if (strcmp(argv[count],"-matrixType") == 0 || strcmp(argv[count],"-matrix")){	
	if (strcmp(argv[count+1],"sparse") == 0)
	  matType = 1;
      }
      else if (strcmp(argv[count],"-rTol") == 0 || strcmp(argv[count],"-relTol") ||
	       strcmp(argv[count],"-relativeTolerance")) {
	if (Tcl_GetDouble(interp, argv[count+1], &rTol) != TCL_OK)
	  return TCL_ERROR;		     
      } else if (strcmp(argv[count],"-aTol") == 0 || strcmp(argv[count],"-absTol") ||
		 strcmp(argv[count],"-absoluteTolerance")) {
	if (Tcl_GetDouble(interp, argv[count+1], &aTol) != TCL_OK)
	  return TCL_ERROR;		     
      } else if (strcmp(argv[count],"-dTol") == 0 || strcmp(argv[count],"-divTol") ||
		 strcmp(argv[count],"-divergenceTolerance")) {
	if (Tcl_GetDouble(interp, argv[count+1], &dTol) != TCL_OK)
	  return TCL_ERROR;		     
      } else if (strcmp(argv[count],"-mIts") == 0 || strcmp(argv[count],"-maxIts") ||
		 strcmp(argv[count],"-maxIterations")) {
	if (Tcl_GetInt(interp, argv[count+1], &maxIts) != TCL_OK)
	  return TCL_ERROR;		     
      } else if (strcmp(argv[count],"-KSP") == 0 || strcmp(argv[count],"-KSPType")){	
	if (strcmp(argv[count+1],"KSPCG") == 0)
	  method = KSPCG;
	else if (strcmp(argv[count+1],"KSPBICG") == 0)
	  method = KSPBICG;
	else if (strcmp(argv[count+1],"KSPRICHARDSON") == 0)
	  method = KSPRICHARDSON;
	else if (strcmp(argv[count+1],"KSPCHEBYCHEV") == 0)
	  method = KSPCHEBYCHEV;
	else if (strcmp(argv[count+1],"KSPGMRES") == 0)
	  method = KSPGMRES;
      } else if (strcmp(argv[count],"-PC") == 0 || strcmp(argv[count],"-PCType")){	
	if ((strcmp(argv[count+1],"PCJACOBI") == 0) || (strcmp(argv[count+1],"JACOBI") == 0))
	  preconditioner = PCJACOBI;
	else if ((strcmp(argv[count+1],"PCILU") == 0) || (strcmp(argv[count+1],"ILU") == 0))
	  preconditioner = PCILU;
	else if ((strcmp(argv[count+1],"PCICC") == 0) || (strcmp(argv[count+1],"ICC") == 0)) 
	  preconditioner = PCICC;
	else if ((strcmp(argv[count+1],"PCBJACOBI") == 0) || (strcmp(argv[count+1],"BIJACOBI") == 0))
	  preconditioner = PCBJACOBI;
	else if ((strcmp(argv[count+1],"PCNONE") == 0) || (strcmp(argv[count+1],"NONE") == 0))
	  preconditioner = PCNONE;
      }
      count+=2;
    }

    if (matType == 0) {
      PetscSolver *theSolver = new PetscSolver(method, preconditioner, rTol, aTol, dTol, maxIts);
      theSOE = new PetscSOE(*theSolver);
    } else {
      PetscSparseSeqSolver *theSolver = new PetscSparseSeqSolver(method, preconditioner, rTol, aTol, dTol, maxIts);
      theSOE = new SparseGenRowLinSOE(*theSolver);
    }
  }


#endif


#ifdef _MUMPS

  else if (strcmp(argv[1],"Mumps") == 0) {
    
#ifdef _PARALLEL_PROCESSING
    MumpsParallelSolver *theSolver = new MumpsParallelSolver();
    theSOE = new MumpsParallelSOE(*theSolver);
#else
    MumpsSolver *theSolver = new MumpsSolver();
    theSOE = new MumpsSOE(*theSolver);

#endif

  }

#endif

  
  else {

    //
    // maybe a package
    //
    
    // try existing loaded packages
    
    LinearSOE_PackageCommand *soeCommands = theLinearSOE_PackageCommands;
    bool found = false;
    while (soeCommands != NULL && found == false) {
      if (strcmp(argv[1], soeCommands->funcName) == 0) {
	int result = (*(soeCommands->funcPtr))(clientData, interp, argc, argv, &theBroker, &theSOE);
	found = true;
      } else
	soeCommands = soeCommands->next;
    }
    
    if (found == false) {
      // load new package
    
      void *libHandle;
      int (*funcPtr)(ClientData clientData, Tcl_Interp *interp,  int argc, 
		     TCL_Char **argv, FEM_ObjectBroker *, LinearSOE **);       
      int linearSOE_NameLength = strlen(argv[1]);
      char *tclFuncName = new char[linearSOE_NameLength+12];
      strcpy(tclFuncName, "TclCommand_");
      strcpy(&tclFuncName[11], argv[1]);    
      
      int res = getLibraryFunction(argv[1], tclFuncName, &libHandle, (void **)&funcPtr);
      
      if (res == 0) {
	char *linearSOE_Name = new char[linearSOE_NameLength+1];
	strcpy(linearSOE_Name, argv[1]);
	LinearSOE_PackageCommand *theSOE_Command = new LinearSOE_PackageCommand;
	theSOE_Command->funcPtr = funcPtr;
	theSOE_Command->funcName = linearSOE_Name;	
	theSOE_Command->next = theLinearSOE_PackageCommands;
	theLinearSOE_PackageCommands = theSOE_Command;
	
	int result = (*funcPtr)(clientData, interp,
				argc, 
				argv,
				&theBroker,
				&theSOE);	
      }
    }
  }
    
  // if the analysis exists - we want to change the SOEif

  if (theSOE != 0) {
    if (theStaticAnalysis != 0)
      theStaticAnalysis->setLinearSOE(*theSOE);
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setLinearSOE(*theSOE);
    
    return TCL_OK;
  }


  // if the analysis exists - we want to change the SOE
  if (theStaticAnalysis != 0)
    theStaticAnalysis->setLinearSOE(*theSOE);
  else if (theTransientAnalysis != 0)
    theTransientAnalysis->setLinearSOE(*theSOE);  

#ifdef _PARALLEL_PROCESSING
    if (theStaticAnalysis != 0 || theTransientAnalysis != 0) {
      SubdomainIter &theSubdomains = theDomain.getSubdomains();
      Subdomain *theSub;
      while ((theSub = theSubdomains()) != 0) {
	theSub->setAnalysisLinearSOE(*theSOE);
      }
    }
#endif

  return TCL_ERROR;
}



//
// command invoked to allow the Numberer objects to be built
//
int 
specifyNumberer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain numberer
  if (argc < 2) {
      opserr << "WARNING need to specify a Nemberer type \n";
      return TCL_ERROR;
  }    

#ifdef _PARALLEL_PROCESSING

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1],"Plain") == 0) {
    theNumberer = new ParallelNumberer();       
  } else if (strcmp(argv[1],"RCM") == 0) {
    RCM *theRCM = new RCM(false);	
    theNumberer = new ParallelNumberer(*theRCM);    	
  } else {
    opserr << "WARNING No Numberer type exists (Plain, RCM only) \n";
    return TCL_ERROR;
  }    

#else

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1],"Plain") == 0) {
    theNumberer = new PlainNumberer();       
  } else if (strcmp(argv[1],"RCM") == 0) {
    RCM *theRCM = new RCM(false);	
    theNumberer = new DOF_Numberer(*theRCM);    	
  } else {
    opserr << "WARNING No Numberer type exists (Plain, RCM only) \n";
    return TCL_ERROR;
  }    
#endif

  return TCL_OK;
}




//
// command invoked to allow the ConstraintHandler object to be built
//
int 
specifyConstraintHandler(ClientData clientData, Tcl_Interp *interp, int argc, 
			 TCL_Char **argv)
{
  // make sure at least one other argument to contain numberer
  if (argc < 2) {
      opserr << "WARNING need to specify a Nemberer type \n";
      return TCL_ERROR;
  }    

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1],"Plain") == 0) 
    theHandler = new PlainHandler();       

  else if (strcmp(argv[1],"Penalty") == 0) {
    if (argc < 4) {
      opserr << "WARNING: need to specify alpha: handler Penalty alpha \n";
      return TCL_ERROR;
    }    
    double alpha1, alpha2;
    if (Tcl_GetDouble(interp, argv[2], &alpha1) != TCL_OK)	
      return TCL_ERROR;	
    if (Tcl_GetDouble(interp, argv[3], &alpha2) != TCL_OK)	
      return TCL_ERROR;	
    theHandler = new PenaltyConstraintHandler(alpha1, alpha2);
  }

  /****** adding later
  else if (strcmp(argv[1],"PenaltyNoHomoSPMultipliers") == 0) {
    if (argc < 4) {
      opserr << "WARNING: need to specify alpha: handler Penalty alpha \n";
      return TCL_ERROR;
    }    
    double alpha1, alpha2;
    if (Tcl_GetDouble(interp, argv[2], &alpha1) != TCL_OK)	
      return TCL_ERROR;	
    if (Tcl_GetDouble(interp, argv[3], &alpha2) != TCL_OK)	
      return TCL_ERROR;	
    theHandler = new PenaltyHandlerNoHomoSPMultipliers(alpha1, alpha2);
  }
  ***********************/
  else if (strcmp(argv[1],"Lagrange") == 0) {
    double alpha1 = 1.0;
    double alpha2 = 1.0;
    if (argc == 4) {
      if (Tcl_GetDouble(interp, argv[2], &alpha1) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[3], &alpha2) != TCL_OK)	
	return TCL_ERROR;	
    }
    theHandler = new LagrangeConstraintHandler(alpha1, alpha2);
  }  
  
  else if (strcmp(argv[1],"Transformation") == 0) {
    theHandler = new TransformationConstraintHandler();
  }    

  else {
    opserr << "WARNING No ConstraintHandler type exists (Plain, Penalty,\n";
    opserr << " Lagrange, Transformation) only\n";
    return TCL_ERROR;
  }    
  return TCL_OK;
}



//
// command invoked to allow the SolnAlgorithm object to be built
//
int
specifyAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, 
		 TCL_Char **argv)
{
  // make sure at least one other argument to contain numberer
  if (argc < 2) {
      opserr << "WARNING need to specify an Algorithm type \n";
      return TCL_ERROR;
  }    
  EquiSolnAlgo *theNewAlgo = 0;

  // check argv[1] for type of Algorithm and create the object
  if (strcmp(argv[1],"Linear") == 0) 
    theNewAlgo = new Linear();       

  else if (strcmp(argv[1],"Newton") == 0) {
    int formTangent = CURRENT_TANGENT;
    if (argc > 2) {
      if (strcmp(argv[2],"-secant") == 0) {
	formTangent = CURRENT_SECANT;
      } else if (strcmp(argv[2],"-initial") == 0) {
	formTangent = INITIAL_TANGENT;
      } else if ((strcmp(argv[2],"-initialThenCurrent") == 0) || 
		 (strcmp(argv[2],"-initialCurrent") == 0))  {
	formTangent = INITIAL_THEN_CURRENT_TANGENT;
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }
    theNewAlgo = new NewtonRaphson(*theTest, formTangent); 
  }

  else if (strcmp(argv[1],"KrylovNewton") == 0) {
    int formTangent = CURRENT_TANGENT;
    int maxDim = -1;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i],"-secant") == 0) {
	formTangent = CURRENT_SECANT;
      } else if (strcmp(argv[i],"-initial") == 0) {
	formTangent = INITIAL_TANGENT;
      } else if (strcmp(argv[i++],"-maxDim") == 0 && i < argc) {
	maxDim = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }
    if (maxDim == -1)
      theNewAlgo = new KrylovNewton(*theTest, formTangent); 
    else
      theNewAlgo = new KrylovNewton(*theTest, formTangent, maxDim); 
  }

  else if (strcmp(argv[1],"PeriodicNewton") == 0) {
    int formTangent = CURRENT_TANGENT;
    int maxDim = -1;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i],"-secant") == 0) {
	formTangent = CURRENT_SECANT;
      } else if (strcmp(argv[i],"-initial") == 0) {
	formTangent = INITIAL_TANGENT;
      } else if (strcmp(argv[i++],"-maxCount") == 0 && i < argc) {
	maxDim = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }
    if (maxDim == -1)
      theNewAlgo = new PeriodicNewton(*theTest, formTangent); 
    else
      theNewAlgo = new PeriodicNewton(*theTest, formTangent, maxDim); 
  }

  else if (strcmp(argv[1],"Broyden") == 0) {
    int formTangent = CURRENT_TANGENT;
    int count = -1;

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i],"-secant") == 0) {
	formTangent = CURRENT_SECANT;
      } else if (strcmp(argv[i],"-initial") == 0) {
	formTangent = INITIAL_TANGENT;
      } else if (strcmp(argv[i++],"-count") == 0 && i < argc) {
	count = atoi(argv[i]);
      }
    }

    if (count == -1)
      theNewAlgo = new Broyden(*theTest, formTangent); 
    else
      theNewAlgo = new Broyden(*theTest, formTangent, count); 
  }

  else if (strcmp(argv[1],"BFGS") == 0) {
    int formTangent = CURRENT_TANGENT;
    int count = -1;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i],"-secant") == 0) {
	formTangent = CURRENT_SECANT;
      } else if (strcmp(argv[i],"-initial") == 0) {
	formTangent = INITIAL_TANGENT;
      } else if (strcmp(argv[i++],"-count") == 0 && i < argc) {
	count = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }

    if (count == -1)
      theNewAlgo = new BFGS(*theTest, formTangent); 
    else
      theNewAlgo = new BFGS(*theTest, formTangent, count); 
  }
  
  else if (strcmp(argv[1],"ModifiedNewton") == 0) {
    int formTangent = CURRENT_TANGENT;
    if (argc > 2) {
      if (strcmp(argv[2],"-secant") == 0) {
	formTangent = CURRENT_SECANT;
      } else if (strcmp(argv[2],"-initial") == 0) {
	formTangent = INITIAL_TANGENT;
      }
    }
    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }
      
    theNewAlgo = new ModifiedNewton(*theTest, formTangent); 
  }  
  
  else if (strcmp(argv[1],"NewtonLineSearch") == 0) {
      if (theTest == 0) {
	  opserr << "ERROR: No ConvergenceTest yet specified\n";
	  return TCL_ERROR;	  
      }

      int    count = 2;
      
      // set some default variable
      double tol        = 0.8;
      int    maxIter    = 10;
      double maxEta     = 10.0;
      double minEta     = 0.1;
      int    pFlag      = 1;
      int    typeSearch = 0;
      
      while (count < argc) {
	if (strcmp(argv[count], "-tol") == 0) {
	  count++;
	  if (Tcl_GetDouble(interp, argv[count], &tol) != TCL_OK)	
	    return TCL_ERROR;	      	  
	  count++;
	} else if (strcmp(argv[count], "-maxIter") == 0) {
	  count++;
	  if (Tcl_GetInt(interp, argv[count], &maxIter) != TCL_OK)	
	    return TCL_ERROR;	      	  
	  count++;	  
	} else if (strcmp(argv[count], "-pFlag") == 0) {
	  count++;
	  if (Tcl_GetInt(interp, argv[count], &pFlag) != TCL_OK)	
	    return TCL_ERROR;	      	  
	  count++;
	} else if (strcmp(argv[count], "-minEta") == 0) {
	  count++;
	  if (Tcl_GetDouble(interp, argv[count], &minEta) != TCL_OK)	
	    return TCL_ERROR;	      	  
	  count++;
	} else if (strcmp(argv[count], "-maxEta") == 0) {
	  count++;
	  if (Tcl_GetDouble(interp, argv[count], &maxEta) != TCL_OK)	
	    return TCL_ERROR;	      	  
	  count++;
	} else if (strcmp(argv[count], "-type") == 0) {
	  count++;
	  if (strcmp(argv[count], "Bisection") == 0) 
	    typeSearch = 1;
	  else if (strcmp(argv[count], "Secant") == 0) 
	    typeSearch = 2;
	  else if (strcmp(argv[count], "RegulaFalsi") == 0) 
	    typeSearch = 3;
	  else if (strcmp(argv[count], "LinearInterpolated") == 0) 
	    typeSearch = 3;
	  else if (strcmp(argv[count], "InitialInterpolated") == 0) 
	    typeSearch = 0;
	  count++;
	} else
	  count++;
      }
      
      LineSearch *theLineSearch = 0;      
      if (typeSearch == 0)
	theLineSearch = new InitialInterpolatedLineSearch(tol, maxIter, minEta, maxEta, pFlag);
							  
      else if (typeSearch == 1)
	theLineSearch = new BisectionLineSearch(tol, maxIter, minEta, maxEta, pFlag);
      else if (typeSearch == 2)
	theLineSearch = new SecantLineSearch(tol, maxIter, minEta, maxEta, pFlag);
      else if (typeSearch == 3)
	theLineSearch = new RegulaFalsiLineSearch(tol, maxIter, minEta, maxEta, pFlag);

      theNewAlgo = new NewtonLineSearch(*theTest, theLineSearch); 
  }

  else {
      opserr << "WARNING No EquiSolnAlgo type exists (Linear, Newton only) \n";
      return TCL_ERROR;
  }    


  if (theNewAlgo != 0) {
    theAlgorithm = theNewAlgo;
    
    // if the analysis exists - we want to change the SOE
    if (theStaticAnalysis != 0)
      theStaticAnalysis->setAlgorithm(*theAlgorithm);
    else if (theTransientAnalysis != 0)
      theTransientAnalysis->setAlgorithm(*theAlgorithm);  

#ifdef _PARALLEL_PROCESSING
    if (theStaticAnalysis != 0 || theTransientAnalysis != 0) {
      SubdomainIter &theSubdomains = theDomain.getSubdomains();
      Subdomain *theSub;
      while ((theSub = theSubdomains()) != 0) {
	theSub->setAnalysisAlgorithm(*theAlgorithm);
      }
    }
#endif
  }

  return TCL_OK;
}


//
// command invoked to allow the SolnAlgorithm object to be built
//
int
specifyCTest(ClientData clientData, Tcl_Interp *interp, int argc, 
	     TCL_Char **argv)
{
  // make sure at least one other argument to contain numberer
  if (argc < 2) {
      opserr << "WARNING need to specify a ConvergenceTest Type type \n";
      return TCL_ERROR;
  }    

  // get the tolerence first
  double tol = 0.0;
  int numIter = 0;
  int printIt = 0;
  int normType = 2;

  if (strcmp(argv[1],"FixedNumIter") != 0) {
    if (argc == 4) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
    } else if (argc == 5) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &printIt) != TCL_OK)	
	return TCL_ERROR;			  
    } else if (argc == 6) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &printIt) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[5], &normType) != TCL_OK)	
	return TCL_ERROR;			  
    }
  } else {
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
    } else if (argc == 4) {
      if (Tcl_GetInt(interp, argv[2], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[3], &printIt) != TCL_OK)	
	return TCL_ERROR;			  
    } else if (argc == 5) {
      if (Tcl_GetInt(interp, argv[2], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[3], &printIt) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &normType) != TCL_OK)	
	return TCL_ERROR;			  
    }
  }

  ConvergenceTest *theNewTest = 0;

  if (numIter == 0) {
    opserr << "ERROR: no numIter specified in test command\n";
    return TCL_ERROR;
  }

  if (strcmp(argv[1],"FixedNumIter") == 0)
    theNewTest = new CTestFixedNumIter(numIter,printIt,normType);             
  else {
    if (tol == 0.0) {
      opserr << "ERROR: no tolerance specified in test command\n";
      return TCL_ERROR;
    }
    if (strcmp(argv[1],"NormUnbalance") == 0) 
      theNewTest = new CTestNormUnbalance(tol,numIter,printIt,normType);       
    else if (strcmp(argv[1],"NormDispIncr") == 0) 
      theNewTest = new CTestNormDispIncr(tol,numIter,printIt,normType);             
    else if (strcmp(argv[1],"EnergyIncr") == 0) 
      theNewTest = new CTestEnergyIncr(tol,numIter,printIt,normType);             
    else if (strcmp(argv[1],"RelativeNormUnbalance") == 0) 
      theNewTest = new CTestRelativeNormUnbalance(tol,numIter,printIt,normType);       
    else if (strcmp(argv[1],"RelativeNormDispIncr") == 0) 
      theNewTest = new CTestRelativeNormDispIncr(tol,numIter,printIt,normType);             
    else if (strcmp(argv[1],"RelativeEnergyIncr") == 0) 
      theNewTest = new CTestRelativeEnergyIncr(tol,numIter,printIt,normType);             
    else if (strcmp(argv[1],"RelativeTotalNormDispIncr") == 0) 
      theNewTest = new CTestRelativeTotalNormDispIncr(tol,numIter,printIt,normType);             
    else {
      opserr << "WARNING No ConvergenceTest type (NormUnbalance, NormDispIncr, EnergyIncr, \n";
      opserr << "RelativeNormUnbalance, RelativeNormDispIncr, RelativeEnergyIncr, \n";
      opserr << "RelativeTotalNormDispIncr, FixedNumIter)\n";
      return TCL_ERROR;
    }    
  }

  if (theNewTest != 0) {
    theTest = theNewTest;
  

  // if the analysis exists - we want to change the Test
  if (theStaticAnalysis != 0)
    theStaticAnalysis->setConvergenceTest(*theTest);

  else if (theTransientAnalysis != 0)
    theTransientAnalysis->setConvergenceTest(*theTest); 

#ifdef _PARALLEL_PROCESSING
    if (theStaticAnalysis != 0 || theTransientAnalysis != 0) {
      SubdomainIter &theSubdomains = theDomain.getSubdomains();
      Subdomain *theSub;
      while ((theSub = theSubdomains()) != 0) {
	theSub->setAnalysisConvergenceTest(*theTest);;
      }
    }
#endif
  }
  
  return TCL_OK;
}





//
// command invoked to allow the Integrator object to be built
//
int 
specifyIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, 
		  TCL_Char **argv)
{
  // make sure at least one other argument to contain integrator
  if (argc < 2) {
      opserr << "WARNING need to specify an Integrator type \n";
      return TCL_ERROR;
  }    

  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1],"LoadControl") == 0) {
      double dLambda;
      double minIncr, maxIncr;
      int numIter;
      if (argc < 3) {
	opserr << "WARNING incorrect # args - integrator LoadControl dlam <Jd dlamMin dlamMax>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &dLambda) != TCL_OK)	
	return TCL_ERROR;	
      if (argc > 5) {
	if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	  return TCL_ERROR;	
	if (Tcl_GetDouble(interp, argv[4], &minIncr) != TCL_OK)	
	  return TCL_ERROR;	
	if (Tcl_GetDouble(interp, argv[5], &maxIncr) != TCL_OK)	
	  return TCL_ERROR;	  
      }
      else {
	minIncr = dLambda;
	maxIncr = dLambda;
	numIter = 1;
      }
      theStaticIntegrator = new LoadControl(dLambda, numIter, minIncr, maxIncr);       

  // if the analysis exists - we want to change the Integrator
  if (theStaticAnalysis != 0)
    theStaticAnalysis->setIntegrator(*theStaticIntegrator);
  }

  
  else if (strcmp(argv[1],"ArcLength") == 0) {
      double arcLength;
      double alpha;
      if (argc != 4) {
	opserr << "WARNING integrator ArcLength arcLength alpha \n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK)	
	return TCL_ERROR;	
      theStaticIntegrator = new ArcLength(arcLength,alpha);       

  // if the analysis exists - we want to change the Integrator
  if (theStaticAnalysis != 0)
    theStaticAnalysis->setIntegrator(*theStaticIntegrator);
  }

  else if (strcmp(argv[1],"ArcLength1") == 0) {
      double arcLength;
      double alpha;
      if (argc != 4) {
	opserr << "WARNING integrator ArcLength1 arcLength alpha \n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK)	
	return TCL_ERROR;	
      theStaticIntegrator = new ArcLength1(arcLength,alpha);       

  // if the analysis exists - we want to change the Integrator
  if (theStaticAnalysis != 0)
    theStaticAnalysis->setIntegrator(*theStaticIntegrator);
  }
  /************************added for HSConstraint*************************************/
  
  else if (strcmp(argv[1],"HSConstraint") == 0) {
      double arcLength;
      double psi_u;
      double psi_f;
      double u_ref;
      if (argc < 3) {
	opserr << "WARNING integrator HSConstraint <arcLength> <psi_u> <psi_f> <u_ref> \n";
	return TCL_ERROR;
      }    
      if (argc >= 3 && Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)	
	return TCL_ERROR;	
      if (argc>=4 && Tcl_GetDouble(interp, argv[3], &psi_u) != TCL_OK)	
	return TCL_ERROR;	
      if (argc>=5 && Tcl_GetDouble(interp, argv[4], &psi_f) != TCL_OK)	
	return TCL_ERROR;	
      if (argc==6 && Tcl_GetDouble(interp, argv[5], &u_ref) != TCL_OK)	
	return TCL_ERROR;	
	switch(argc)
	{
		case 3:
		    	theStaticIntegrator = new HSConstraint(arcLength);       
		case 4:
		      	theStaticIntegrator = new HSConstraint(arcLength, psi_u);       
		case 5:
		      	theStaticIntegrator = new HSConstraint(arcLength, psi_u, psi_f);       
		case 6:
		      	theStaticIntegrator = new HSConstraint(arcLength, psi_u, psi_f, u_ref);       
	}
    // if the analysis exists - we want to change the Integrator
    if (theStaticAnalysis != 0)
    	theStaticAnalysis->setIntegrator(*theStaticIntegrator);
  }
  /*********************************************************************************/
  
  else if (strcmp(argv[1],"MinUnbalDispNorm") == 0) {
      double lambda11, minlambda, maxlambda;
      int numIter;
      if (argc < 3) {
	opserr << "WARNING integrator MinUnbalDispNorm lambda11 <Jd minLambda1j maxLambda1j>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &lambda11) != TCL_OK)	
	return TCL_ERROR;	
      if (argc > 5) {
	if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	  return TCL_ERROR;	
	if (Tcl_GetDouble(interp, argv[4], &minlambda) != TCL_OK)	
	  return TCL_ERROR;	
	if (Tcl_GetDouble(interp, argv[5], &maxlambda) != TCL_OK)	
	  return TCL_ERROR;	
      }
      else {
	minlambda = lambda11;
	maxlambda = lambda11;
	numIter = 1;
	argc += 3;
      }

      int signFirstStepMethod = SIGN_LAST_STEP;
      if (argc == 7)
	if ((strcmp(argv[argc-1],"-determinant") == 0) ||
	    (strcmp(argv[argc-1],"-det") == 0))
	    signFirstStepMethod = CHANGE_DETERMINANT;	    

      theStaticIntegrator = new MinUnbalDispNorm(lambda11,numIter,minlambda,maxlambda,signFirstStepMethod);

      // if the analysis exists - we want to change the Integrator
      if (theStaticAnalysis != 0)
	theStaticAnalysis->setIntegrator(*theStaticIntegrator);
  }
  
  else if (strcmp(argv[1],"DisplacementControl") == 0) {
      int node;
      int dof;
      double increment, minIncr, maxIncr;
      int numIter;
      if (argc < 5) {
	opserr << "WARNING integrator DisplacementControl node dof dU \n";
	opserr << "<Jd minIncrement maxIncrement>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetInt(interp, argv[2], &node) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[4], &increment) != TCL_OK)	
	return TCL_ERROR;	      
      if (argc > 7) {
	if (Tcl_GetInt(interp, argv[5], &numIter) != TCL_OK)	
	  return TCL_ERROR;	
	if (Tcl_GetDouble(interp, argv[6], &minIncr) != TCL_OK)	
	  return TCL_ERROR;	
	if (Tcl_GetDouble(interp, argv[7], &maxIncr) != TCL_OK)	
	  return TCL_ERROR;	  
      }
      else {
	minIncr = increment;
	maxIncr = increment;
	numIter = 1;
      }

      Node *theNode = theDomain.getNode(node);
      if (theNode == 0) {
	opserr << "WARNING integrator DisplacementControl node dof dU : Node does not exist\n";
	return TCL_ERROR;	  
      }


#ifdef _PARALLEL_PROCESSING
      theStaticIntegrator = new DistributedDisplacementControl(node,dof-1,increment,
							       numIter, minIncr, maxIncr);
#else
      int numDOF = theNode->getNumberDOF();
      if (dof <= 0 || dof > numDOF) {
	opserr << "WARNING integrator DisplacementControl node dof dU : invalid dof given\n";
	return TCL_ERROR;	  
      }

      theStaticIntegrator = new DisplacementControl(node, dof-1, increment, &theDomain,
						    numIter, minIncr, maxIncr);
#endif

      // if the analysis exists - we want to change the Integrator
      if (theStaticAnalysis != 0)
	theStaticAnalysis->setIntegrator(*theStaticIntegrator);
  }  
  
  else if (strcmp(argv[1],"Newmark") == 0) {
      double gamma;
      double beta;
      double alphaM = 0.0;
      double betaK  = 0.0;
      double betaKi = 0.0;
      double betaKc = 0.0;
      if (argc != 4 && argc != 8) {
	opserr << "WARNING integrator Newmark gamma beta <alphaM betaK betaK0 betaKc>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator Newmark gamma beta - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
	  opserr << "WARNING integrator Newmark gamma beta - undefined beta\n";
	  return TCL_ERROR;	
      }
      if (argc == 7 || argc == 8) {
	  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator Newmark gamma beta alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator Newmark gamma beta alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator Newmark gamma beta alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator Newmark gamma beta alphaM betaK betaKi betaKc - betaKc\n";
	    return TCL_ERROR;	
	  }
      }
      if (argc == 4)
	  theTransientIntegrator = new Newmark(gamma,beta);       
      else
	  theTransientIntegrator = new Newmark(gamma,beta,alphaM,betaK,betaKi,betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }  

  else if (strcmp(argv[1],"NewmarkExplicit") == 0) {
      double gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 7) {
	opserr << "WARNING integrator NewmarkExplicit gamma <alphaM betaK betaKi betaKc>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator NewmarkExplicit gamma - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator NewmarkExplicit gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator NewmarkExplicit gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator NewmarkExplicit gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator NewmarkExplicit gamma alphaM betaK betaKi betaKc - betaKc\n";
	    return TCL_ERROR;	
	  }
      }
      if (argc == 3)
	  theTransientIntegrator = new NewmarkExplicit(gamma);       
      else
	  theTransientIntegrator = new NewmarkExplicit(gamma, alphaM, betaK, betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }  

  else if (strcmp(argv[1],"NewmarkHybridSimulation") == 0) {
      double beta, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 4 && argc != 8) {
	opserr << "WARNING integrator NewmarkHybridSimulation gamma beta <alphaM betaK betaKi betaKc>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator NewmarkHybridSimulation gamma beta - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
	  opserr << "WARNING integrator NewmarkHybridSimulation gamma beta - undefined beta\n";
	  return TCL_ERROR;	
      }
      if (argc == 8) {
	  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator NewmarkHybridSimulation gamma beta alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator NewmarkHybridSimulation gamma beta alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator NewmarkHybridSimulation gamma beta alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator NewmarkHybridSimulation gamma beta alphaM betaK betaKi betaKc - betaKc\n";
	    return TCL_ERROR;	
	  }
      }
      if (theTest == 0) {
        opserr << "WARNING no ConvergenceTest specified yet\n";
        return TCL_ERROR;	  
      }
      if (argc == 4)
	  theTransientIntegrator = new NewmarkHybridSimulation(gamma, beta, *theTest);       
      else
	  theTransientIntegrator = new NewmarkHybridSimulation(gamma, beta, *theTest, alphaM, betaK, betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

#ifdef _RELIABILITY
  else if (strcmp(argv[1],"NewmarkWithSensitivity") == 0) {
	  int assemblyFlag = 0;
      double gamma;
      double beta;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 4 && argc != 6 && argc != 8 && argc != 10) {
	     interp->result = "WARNING integrator Newmark gamma beta <alphaM?  betaKcurrent?  betaKi? betaKlastCommitted?> <-assemble tag?> ";
	     return TCL_ERROR;
      }
	  
	  // Take care of argc == 4, the basic case
      if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
		  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;	  
		  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
		  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
		  return TCL_ERROR;	
      }

	  // If only assembly flag is given extra
	  if (argc == 6) {
		  if (strcmp(argv[4],"-assemble") != 0) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
		  }
		  if (Tcl_GetInt(interp, argv[5], &assemblyFlag) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
	  }
	  // If only extra integrator (damping) parameters are given extra
      if (argc == 8) {
		  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
		  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
		  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
		  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
      }
	  // If everything is given extra
	  if (argc == 10) {
		  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
		  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
		  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
		  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
		  if (strcmp(argv[8],"-assemble") != 0) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
		  }
		  if (Tcl_GetInt(interp, argv[9], &assemblyFlag) != TCL_OK) {
			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
			  return TCL_ERROR;	
		  }
	  }

      if (argc == 4 || argc == 6) {
	theNSI = new NewmarkSensitivityIntegrator(assemblyFlag, gamma,beta);       
      }
      else {
	theNSI = new NewmarkSensitivityIntegrator(assemblyFlag, gamma,beta,alphaM,betaK,betaKi,betaKc);
      }
      theTransientIntegrator = theNSI;


      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }  

#endif
  
  else if (strcmp(argv[1],"HHT") == 0) {
      double alpha, beta, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 5 && argc != 7 && argc != 9) {
	opserr << "WARNING integrator HHT alpha <alphaM betaK betaKi betaKc>\n";
	opserr << "     or integrator HHT alpha beta gamma <alphaM betaK betaKi betaKc>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
	  opserr << "WARNING integrator HHT alpha - undefined alpha\n";	  
	  return TCL_ERROR;	
      }
      if (argc == 5 || argc == 9) {
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
	  opserr << "WARNING integrator HHT alpha beta gamma - undefined beta\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator HHT alpha beta gamma - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHT alpha alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHT alpha alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHT alpha alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHT alpha alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 9) {
	  if (Tcl_GetDouble(interp, argv[5], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHT alpha beta gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHT alpha beta gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHT alpha beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[8], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHT alpha beta gamma alphaM betaK betaKi betaKc - betaKc\n";
	    return TCL_ERROR;	
	  }
      }      
      if (argc == 3)
	  theTransientIntegrator = new HHT(alpha);       
      else if (argc == 5)
      theTransientIntegrator = new HHT(alpha, beta, gamma);
      else if (argc == 7)
	  theTransientIntegrator = new HHT(alpha, alphaM, betaK, betaKi, betaKc);
      else if (argc == 9)
	  theTransientIntegrator = new HHT(alpha, beta, gamma, alphaM, betaK, betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }    

  else if (strcmp(argv[1],"HHTGeneralized") == 0) {
      double rhoInf, alphaI, alphaF, beta, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 6 && argc != 7 && argc != 10) {
	opserr << "WARNING integrator HHTGeneralized rhoInf <alphaM betaK betaKi betaKc>\n";
    opserr << "     or integrator HHTGeneralized alphaI alphaF beta gamma <alphaM betaK betaKi betaKc>\n";
	return TCL_ERROR;
      }
      if (argc == 3 || argc == 7) {
      if (Tcl_GetDouble(interp, argv[2], &rhoInf) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralized rhoInf - undefined rhoInf\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 6 || argc == 10) {
      if (Tcl_GetDouble(interp, argv[2], &alphaI) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralized alphaI alphaF beta gamma - undefined alphaI\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &alphaF) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralized alphaI alphaF beta gamma - undefined alphaF\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &beta) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralized alphaI alphaF beta gamma - undefined beta\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[5], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralized alphaI alphaF beta gamma - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralized rhoInf alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralized rhoInf alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralized rhoInf alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHTGeneralized rhoInf alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 10) {
	  if (Tcl_GetDouble(interp, argv[6], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralized alphaI alphaF beta gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralized alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[8], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralized alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[9], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHTGeneralized alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaKc\n";
	    return TCL_ERROR;	
	  }
      }      
      if (argc == 3)
	  theTransientIntegrator = new HHTGeneralized(rhoInf);       
      else if (argc == 6)
      theTransientIntegrator = new HHTGeneralized(alphaI, alphaF, beta, gamma);
      else if (argc == 7)
	  theTransientIntegrator = new HHTGeneralized(rhoInf, alphaM, betaK, betaKi, betaKc);
      else if (argc == 10)
	  theTransientIntegrator = new HHTGeneralized(alphaI, alphaF, beta, gamma, alphaM, betaK, betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1],"HHTExplicit") == 0) {
      double alpha, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 4 && argc != 7 && argc != 8) {
	      opserr << "WARNING integrator HHTExplicit alpha <alphaM betaK betaKi betaKc>\n";
          opserr << "     or integrator HHTExplicit alpha gamma <alphaM betaK betaKi betaKc>\n";
	      return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
	      opserr << "WARNING integrator HHTExplicit alpha - undefined alpha\n";	  
	      return TCL_ERROR;	
      }
      if (argc == 4 || argc == 8) {
      if (Tcl_GetDouble(interp, argv[3], &gamma) != TCL_OK) {
	      opserr << "WARNING integrator HHTExplicit alpha gamma - undefined gamma\n";	  
	      return TCL_ERROR;	
      }
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHTExplicit alpha alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHTExplicit alpha alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHTExplicit alpha alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHTExplicit alpha alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 8) {
	  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHTExplicit alpha gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHTExplicit alpha gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHTExplicit alpha gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHTExplicit alpha gamma alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 3)
	  theTransientIntegrator = new HHTExplicit(alpha);       
      else if (argc == 4)
      theTransientIntegrator = new HHTExplicit(alpha, gamma);
      else if (argc == 7)
	  theTransientIntegrator = new HHTExplicit(alpha, alphaM, betaK, betaKi, betaKc);
      else if (argc == 8)
      theTransientIntegrator = new HHTExplicit(alpha, gamma, alphaM, betaK, betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1],"HHTGeneralizedExplicit") == 0) {
      double rhoB, alphaI, alphaF, beta, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 4 && argc != 6 && argc != 8 && argc != 10) {
	opserr << "WARNING integrator HHTGeneralizedExplicit rhoB alphaF <alphaM betaK betaKi betaKc>\n";
    opserr << "     or integrator HHTGeneralizedExplicit alphaI alphaF beta gamma <alphaM betaK betaKi betaKc>\n";
	return TCL_ERROR;
      }
      if (argc == 4 || argc == 8) {
      if (Tcl_GetDouble(interp, argv[2], &rhoB) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralizedExplicit rhoB alphaF - undefined rhoB\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &alphaF) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralizedExplicit rhoB alphaF - undefined alphaF\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 6 || argc == 10) {
      if (Tcl_GetDouble(interp, argv[2], &alphaI) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralizedExplicit alphaI alphaF beta gamma - undefined alphaI\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &alphaF) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralizedExplicit alphaI alphaF beta gamma - undefined alphaF\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &beta) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralizedExplicit alphaI alphaF beta gamma - undefined beta\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[5], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator HHTGeneralizedExplicit alphaI alphaF beta gamma - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 8) {
	  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralizedExplicit rhoInf alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralizedExplicit rhoInf alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralizedExplicit rhoInf alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHTGeneralizedExplicit rhoInf alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 10) {
	  if (Tcl_GetDouble(interp, argv[6], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralizedExplicit alphaI alphaF beta gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralizedExplicit alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[8], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHTGeneralizedExplicit alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[9], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHTGeneralizedExplicit alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaKc\n";
	    return TCL_ERROR;	
	  }
      }      
      if (argc == 4)
	  theTransientIntegrator = new HHTGeneralizedExplicit(rhoB, alphaF);       
      else if (argc == 6)
      theTransientIntegrator = new HHTGeneralizedExplicit(alphaI, alphaF, beta, gamma);
      else if (argc == 8)
	  theTransientIntegrator = new HHTGeneralizedExplicit(rhoB, alphaF, alphaM, betaK, betaKi, betaKc);
      else if (argc == 10)
	  theTransientIntegrator = new HHTGeneralizedExplicit(alphaI, alphaF, beta, gamma, alphaM, betaK, betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTHybridSimulation") == 0) {
      double rhoInf, alphaI, alphaF, beta, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 6 && argc != 7 && argc != 10) {
	opserr << "WARNING integrator HHTHybridSimulation rhoInf <alphaM betaK betaKi betaKc>\n";
    opserr << "     or integrator HHTHybridSimulation alphaI alphaF beta gamma <alphaM betaK betaKi betaKc>\n";
	return TCL_ERROR;
      }
      if (argc == 3 || argc == 7) {
      if (Tcl_GetDouble(interp, argv[2], &rhoInf) != TCL_OK) {
	  opserr << "WARNING integrator HHTHybridSimulation rhoInf - undefined rhoInf\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 6 || argc == 10) {
      if (Tcl_GetDouble(interp, argv[2], &alphaI) != TCL_OK) {
	  opserr << "WARNING integrator HHTHybridSimulation alphaI alphaF beta gamma - undefined alphaI\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &alphaF) != TCL_OK) {
	  opserr << "WARNING integrator HHTHybridSimulation alphaI alphaF beta gamma - undefined alphaF\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &beta) != TCL_OK) {
	  opserr << "WARNING integrator HHTHybridSimulation alphaI alphaF beta gamma - undefined beta\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[5], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator HHTHybridSimulation alphaI alphaF beta gamma - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHTHybridSimulation rhoInf alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHTHybridSimulation rhoInf alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHTHybridSimulation rhoInf alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHTHybridSimulation rhoInf alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 10) {
	  if (Tcl_GetDouble(interp, argv[6], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHTHybridSimulation alphaI alphaF beta gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHTHybridSimulation alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[8], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHTHybridSimulation alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[9], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHTHybridSimulation alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaKc\n";
	    return TCL_ERROR;	
	  }
      }      
      if (theTest == 0) {
        opserr << "WARNING no ConvergenceTest specified yet\n";
        return TCL_ERROR;	  
      }
      if (argc == 3)
	  theTransientIntegrator = new HHTHybridSimulation(rhoInf, *theTest);       
      else if (argc == 6)
      theTransientIntegrator = new HHTHybridSimulation(alphaI, alphaF, beta, gamma, *theTest);
      else if (argc == 7)
	  theTransientIntegrator = new HHTHybridSimulation(rhoInf, *theTest, alphaM, betaK, betaKi, betaKc);
      else if (argc == 10)
	  theTransientIntegrator = new HHTHybridSimulation(alphaI, alphaF, beta, gamma, *theTest, alphaM, betaK, betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1],"AlphaOS") == 0) {
      double alpha, beta, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 5 && argc != 7 && argc != 9) {
	      opserr << "WARNING integrator AlphaOS alpha <alphaM betaK betaKi betaKc>\n";
	      opserr << "     or integrator AlphaOS alpha beta gamma <alphaM betaK betaKi betaKc>\n";
		  return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOS alpha - undefined alpha\n";	  
	      return TCL_ERROR;	
      }
      if (argc == 5 || argc == 9) {
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
	  opserr << "WARNING integrator AlphaOS alpha beta gamma - undefined beta\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator AlphaOS alpha beta gamma - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOS alpha alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOS alpha alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOS alpha alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator AlphaOS alpha alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 9) {
	  if (Tcl_GetDouble(interp, argv[5], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOS alpha beta gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOS alpha beta gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOS alpha beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[8], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator AlphaOS alpha beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 3)
	  theTransientIntegrator = new AlphaOS(alpha);       
      else if (argc == 5)
	  theTransientIntegrator = new AlphaOS(alpha, beta, gamma);
      else if (argc == 7)
      theTransientIntegrator = new AlphaOS(alpha, alphaM, betaK, betaKi, betaKc);
      else if (argc == 9)
      theTransientIntegrator = new AlphaOS(alpha, beta, gamma, alphaM, betaK, betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }    

  else if (strcmp(argv[1],"AlphaOSGeneralized") == 0) {
      double rhoInf, alphaI, alphaF, beta, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 6 && argc != 7 && argc != 10) {
	      opserr << "WARNING integrator AlphaOSGeneralized rhoInf <alphaM betaK betaKi betaKc>\n";
	      opserr << "     or integrator AlphaOSGeneralized alphaI alphaF beta gamma <alphaM betaK betaKi betaKc>\n";
		  return TCL_ERROR;
      }
      if (argc == 3 || argc == 7) {
      if (Tcl_GetDouble(interp, argv[2], &rhoInf) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOSGeneralized rhoInf - undefined rhoInf\n";	  
	      return TCL_ERROR;	
      }
      }
      if (argc == 6 || argc == 10) {
      if (Tcl_GetDouble(interp, argv[2], &alphaI) != TCL_OK) {
	  opserr << "WARNING integrator AlphaOSGeneralized alphaI alphaF beta gamma - undefined alphaI\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &alphaF) != TCL_OK) {
	  opserr << "WARNING integrator AlphaOSGeneralized alphaI alphaF beta gamma - undefined alphaF\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &beta) != TCL_OK) {
	  opserr << "WARNING integrator AlphaOSGeneralized alphaI alphaF beta gamma - undefined beta\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[5], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator AlphaOSGeneralized alphaI alphaF beta gamma - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOSGeneralized rhoInf alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOSGeneralized rhoInf alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOSGeneralized rhoInf alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator AlphaOSGeneralized rhoInf alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 10) {
	  if (Tcl_GetDouble(interp, argv[6], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOSGeneralized alphaI alphaF beta gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOSGeneralized alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[8], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator AlphaOSGeneralized alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[9], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator AlphaOSGeneralized alphaI alphaF beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 3)
	  theTransientIntegrator = new AlphaOSGeneralized(rhoInf);       
      else if (argc == 6)
	  theTransientIntegrator = new AlphaOSGeneralized(alphaI, alphaF, beta, gamma);
      else if (argc == 7)
      theTransientIntegrator = new AlphaOSGeneralized(rhoInf, alphaM, betaK, betaKi, betaKc);
      else if (argc == 10)
      theTransientIntegrator = new AlphaOSGeneralized(alphaI, alphaF, beta, gamma, alphaM, betaK, betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }    

  else if (strcmp(argv[1],"Collocation") == 0) {
      double theta, beta, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 5 && argc != 7 && argc != 9) {
	opserr << "WARNING integrator Collocation theta <alphaM betaK betaKi betaKc>\n";
	opserr << "     or integrator Collocation theta beta gamma <alphaM betaK betaKi betaKc>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &theta) != TCL_OK) {
	  opserr << "WARNING integrator Collocation theta - undefined theta\n";	  
	  return TCL_ERROR;	
      }
      if (argc == 5 || argc == 9) {
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
	  opserr << "WARNING integrator Collocation theta beta gamma - undefined beta\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator Collocation theta beta gamma - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator Collocation theta alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator Collocation theta alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator Collocation theta alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator Collocation theta alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 9) {
	  if (Tcl_GetDouble(interp, argv[5], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator Collocation theta beta gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator Collocation theta beta gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator Collocation theta beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[8], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator Collocation theta beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 3)
	  theTransientIntegrator = new Collocation(theta);       
      else if (argc == 5)
	  theTransientIntegrator = new Collocation(theta, beta, gamma);       
      else if (argc == 7)
	  theTransientIntegrator = new Collocation(theta, alphaM, betaK, betaKi, betaKc);       
      else if (argc == 9)
	  theTransientIntegrator = new Collocation(theta, beta, gamma, alphaM, betaK, betaKi, betaKc);       

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }    
 
  else if (strcmp(argv[1],"CollocationHybridSimulation") == 0) {
      double theta, beta, gamma;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 5 && argc != 7 && argc != 9) {
	opserr << "WARNING integrator CollocationHybridSimulation theta <alphaM betaK betaKi betaKc>\n";
	opserr << "     or integrator CollocationHybridSimulation theta beta gamma <alphaM betaK betaKi betaKc>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &theta) != TCL_OK) {
	  opserr << "WARNING integrator CollocationHybridSimulation theta - undefined theta\n";	  
	  return TCL_ERROR;	
      }
      if (argc == 5 || argc == 9) {
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
	  opserr << "WARNING integrator CollocationHybridSimulation theta beta gamma - undefined beta\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[4], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator CollocationHybridSimulation theta beta gamma - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator CollocationHybridSimulation theta alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator CollocationHybridSimulation theta alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator CollocationHybridSimulation theta alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator CollocationHybridSimulation theta alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (argc == 9) {
	  if (Tcl_GetDouble(interp, argv[5], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator CollocationHybridSimulation theta beta gamma alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator CollocationHybridSimulation theta beta gamma alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator CollocationHybridSimulation theta beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[8], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator CollocationHybridSimulation theta beta gamma alphaM betaK betaKi betaKc - betaKi\n";
	    return TCL_ERROR;	
	  }
      } 
      if (theTest == 0) {
        opserr << "WARNING no ConvergenceTest specified yet\n";
        return TCL_ERROR;	  
      }
      if (argc == 3)
	  theTransientIntegrator = new CollocationHybridSimulation(theta, *theTest);       
      else if (argc == 5)
	  theTransientIntegrator = new CollocationHybridSimulation(theta, beta, gamma, *theTest);       
      else if (argc == 7)
	  theTransientIntegrator = new CollocationHybridSimulation(theta, *theTest, alphaM, betaK, betaKi, betaKc);       
      else if (argc == 9)
	  theTransientIntegrator = new CollocationHybridSimulation(theta, beta, gamma, *theTest, alphaM, betaK, betaKi, betaKc);       

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }    

  else if (strcmp(argv[1],"Newmark1") == 0) {
      double gamma;
      double beta;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 4 && argc != 8) {
	opserr << "WARNING integrator Newmark1 gamma beta <alphaM> <betaKcurrent> <betaKi> <betaKlastCommitted>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
	  opserr << "WARNING integrator Newmark1 gamma beta - undefined gamma\n";	  
	  return TCL_ERROR;	
      }
      if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
	  opserr << "WARNING integrator Newmark1 gamma beta - undefined beta\n";
	  return TCL_ERROR;	
      }

      if (argc == 8 || argc == 7) {
	  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator Newmark1 gamma beta alphaM betaK betaKi betaKc - betaKc\n";
	    return TCL_ERROR;	
	  }
      }
      if (argc == 4)
	  theTransientIntegrator = new Newmark1(gamma,beta);       
      else
	  theTransientIntegrator = new Newmark1(gamma,beta,alphaM,betaK,betaKi, betaKc);

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1],"HHT1") == 0) {
      double alpha;
      double alphaM, betaK, betaKi, betaKc;
      if (argc != 3 && argc != 7) {
	opserr << "WARNING integrator HHT1 alpha <alphaM> <betaKcurrent> <betaKi> <betaKlastCommitted>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &alpha) != TCL_OK) {
	  opserr << "WARNING integrator HHT alpha - undefined alpha\n";	  
	  return TCL_ERROR;	
      }
      if (argc == 7 || argc == 6) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator HHT1 alpha alphaM betaK betaKi betaKc - alphaM\n";	  
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator HHT1 alpha alphaM betaK betaKi betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator HHT1 alpha alphaM betaK betaKi betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	    opserr << "WARNING integrator HHT1 alpha alphaM betaK betaKi betaKc - betaKc\n";
	    return TCL_ERROR;	
	  }
      }      
      if (argc == 3)
	  theTransientIntegrator = new HHT1(alpha);       
      else
	  theTransientIntegrator = new HHT1(alpha,alphaM,betaK, betaKi, betaKc);       

      // if the analysis exists - we want to change the Integrator
	  if (theTransientAnalysis != 0)
		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }    

  
  else if (strcmp(argv[1],"WilsonTheta") == 0) {
    double theta, alphaM,betaK, betaKi, betaKc;
      if (argc != 3 && argc != 7) {
	opserr << "WARNING integrator WilsonTheta theta <alphaM> <betaK> <betaK0> <betaKc>\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetDouble(interp, argv[2], &theta) != TCL_OK) {
	  opserr << "WARNING integrator WilsonTheta theta - undefined theta\n";	  
	  return TCL_ERROR;	
      }
      if (argc == 7) {
	  if (Tcl_GetDouble(interp, argv[3], &alphaM) != TCL_OK) {
	      opserr << "WARNING integrator WilsonTheta gamma beta alphaM betaK betaK0 betaKc- alphaM\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[4], &betaK) != TCL_OK) {
	      opserr << "WARNING integrator WilsonTheta gamma beta alphaM betaK betaK0 betaKc - betaK\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[5], &betaKi) != TCL_OK) {
	      opserr << "WARNING integrator WilsonTheta gamma beta alphaM betaK betaK0 betaKc - betaKi\n";
	      return TCL_ERROR;	
	  }
	  if (Tcl_GetDouble(interp, argv[6], &betaKc) != TCL_OK) {
	      opserr << "WARNING integrator WilsonTheta gamma beta alphaM betaK betaK0 betaKc - betaKc\n";
	      return TCL_ERROR;	
	  }
      }            
      if (argc == 3)
	  theTransientIntegrator = new WilsonTheta(theta);       
      else
	theTransientIntegrator = new WilsonTheta(theta, alphaM, betaK, betaKi, betaKc);
      
      // if the analysis exists - we want to change the Integrator
      if (theTransientAnalysis != 0)
	theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }      

  else if (strcmp(argv[1],"CentralDifference") == 0) {
    theTransientIntegrator = new CentralDifference();       

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }      

  else if (strcmp(argv[1],"CentralDifferenceAlternative") == 0) {
    theTransientIntegrator = new CentralDifferenceAlternative();       

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }      

  else if (strcmp(argv[1],"CentralDifferenceNoDamping") == 0) {
    theTransientIntegrator = new CentralDifferenceNoDamping();       

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }      

  else {
    opserr << "WARNING No Integrator type exists \n";
    return TCL_ERROR;
  }    

#ifdef _PARALLEL_PROCESSING
  if (theStaticAnalysis != 0 && theStaticIntegrator != 0) {
    IncrementalIntegrator *theIntegrator;
    theIntegrator = theStaticIntegrator;

    SubdomainIter &theSubdomains = theDomain.getSubdomains();
    Subdomain *theSub;
    while ((theSub = theSubdomains()) != 0) {
      theSub->setAnalysisIntegrator(*theIntegrator);
    }
  } else if (theTransientAnalysis != 0 && theTransientIntegrator != 0) {
    IncrementalIntegrator *theIntegrator;
    theIntegrator = theStaticIntegrator;
    
    SubdomainIter &theSubdomains = theDomain.getSubdomains();
    Subdomain *theSub;
    while ((theSub = theSubdomains()) != 0) {
      theSub->setAnalysisIntegrator(*theIntegrator);
    }
  }
#endif

  return TCL_OK;
}


extern int
TclAddRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
	       TCL_Char **argv, Domain &theDomain);

int 
addRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
	    TCL_Char **argv)
{
  return TclAddRecorder(clientData, interp, argc, argv, theDomain);
}

extern int
TclAddAlgorithmRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
			TCL_Char **argv, Domain &theDomain, EquiSolnAlgo *theAlgorithm);

int 
addAlgoRecorder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	if (theAlgorithm != 0)
		return TclAddAlgorithmRecorder(clientData, interp, argc, argv,
			theDomain, theAlgorithm);

	else
		return 0;
}

extern int
TclAddDatabase(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv, 
	       Domain &theDomain, 
	       FEM_ObjectBroker &theBroker);

int 
addDatabase(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  return TclAddDatabase(clientData, interp, argc, argv, theDomain, theBroker);
}


/*
int 
groundExcitation(ClientData clientData, Tcl_Interp *interp, int argc, 
		  TCL_Char **argv)
{
  // make sure at least one other argument to contain integrator
  if (argc < 2) {
      opserr << "WARNING need to specify the commitTag \n";
      return TCL_ERROR;
  }    

  if (strcmp(argv[1],"Single") == 0) {
      if (argc < 4) {
	opserr << "WARNING quake single dof motion\n";
	return TCL_ERROR;
      }    

      int dof;
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK)	
	  return TCL_ERROR;	      
      
      // read in the ground motion
      GroundMotion *theMotion;
      if (strcmp(argv[3],"ElCentro") == 0) {
	  double fact = 1.0;
	  if (argc == 5) {
	      if (Tcl_GetDouble(interp, argv[4], &fact) != TCL_OK)	
		  return TCL_ERROR;	
	  }
	  theMotion = new ElCentroGroundMotion(fact);
      } else {
	  opserr << "WARNING quake Single motion - no motion type exists \n";
	  return TCL_ERROR;      
      }

      Load *theLoad = new SingleExcitation(*theMotion, dof, nextTag++);
      theDomain.addOtherLoad(theLoad);
      return TCL_OK;
  }  
  
  else {
    opserr << "WARNING No quake type exists \n";
    return TCL_ERROR;
  }    
}
*/

int 
rigidLink(ClientData clientData, Tcl_Interp *interp, int argc, 
	  TCL_Char **argv)
{
  if (argc < 4) {
      opserr << "WARNING rigidLink linkType? rNode? cNode?\n";
      return TCL_ERROR;
  }    

  int numMPs = theDomain.getNumMPs();
  int rNode, cNode;
  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read rNode \n";
      return TCL_ERROR;	        
  }
  if (Tcl_GetInt(interp, argv[3], &cNode) != TCL_OK) {
      opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read CNode \n";
      return TCL_ERROR;	        
  }

  // construct a rigid rod or beam depending on 1st arg
  if ((strcmp(argv[1],"-bar") == 0) || (strcmp(argv[1],"bar") == 0)) {
    RigidRod theLink(theDomain, rNode, cNode, numMPs);
  } else if ((strcmp(argv[1],"-beam") == 0) || (strcmp(argv[1],"beam") == 0)) {
    RigidBeam theLink(theDomain, rNode, cNode, numMPs);
  } else {
      opserr << "WARNING rigidLink linkType? rNode? cNode? - unrecognised link type (-bar, -beam) \n";
      return TCL_ERROR;	        
  }

  return TCL_OK;
}



int 
rigidDiaphragm(ClientData clientData, Tcl_Interp *interp, int argc, 
	   TCL_Char **argv)
{
  if (argc < 3) {
      opserr << "WARNING rigidLink perpDirn? rNode? <cNodes?>\n";
      return TCL_ERROR;
  }    

  int rNode, perpDirn;
  if (Tcl_GetInt(interp, argv[1], &perpDirn) != TCL_OK) {
      opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read perpDirn? \n";
      return TCL_ERROR;	        
  }

  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read rNode \n";
      return TCL_ERROR;	        
  }
  
  // read in the constrained Nodes
  int numConstrainedNodes = argc - 3;
  ID constrainedNodes(numConstrainedNodes);
  for (int i=0; i<numConstrainedNodes; i++) {
      int cNode;
      if (Tcl_GetInt(interp, argv[3+i], &cNode) != TCL_OK) {
	  opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read a cNode\n";
	  return TCL_ERROR;	        
      }
      constrainedNodes(i) = cNode;
  }
  int numMPs = theDomain.getNumMPs();
  RigidDiaphragm theLink(theDomain, rNode, constrainedNodes, 
	perpDirn-1, numMPs);

  return TCL_OK;
}


int 
eigenAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, 
	      TCL_Char **argv)
{
     // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - eigen <type> numModes?\n";
	return TCL_ERROR;
    }    

    int typeAlgo = 0; // 0 - frequency/generalized (default), 2 - standard, 2 - buckling
    int typeSolver = 2; // 0 - SymmBandLapack, 1 - SymmSparseArpack, 2 - GenBandArpack (default)
    int loc = 1;

    // Check type of eigenvalue analysis
    while (loc < (argc-1)) {
      
      if ((strcmp(argv[loc],"frequency") == 0) || 
	  (strcmp(argv[loc],"-frequency") == 0) ||
	  (strcmp(argv[loc],"generalized") == 0) ||
	  (strcmp(argv[loc],"-generalized") == 0))
	typeAlgo = 0;

      else if ((strcmp(argv[loc],"standard") == 0) || 
	  (strcmp(argv[loc],"-standard") == 0))
	typeAlgo = 1;

      else if ((strcmp(argv[loc],"symmBandLapack") == 0) || 
	  (strcmp(argv[loc],"-symmBandLapack") == 0))
	typeSolver = 0;

      else if ((strcmp(argv[loc],"symmSparseArpack") == 0) || 
	  (strcmp(argv[loc],"-symmSparseArpack") == 0))
	typeSolver = 1;

      else if ((strcmp(argv[loc],"genBandArpack") == 0) || 
	  (strcmp(argv[loc],"-genBandArpack") == 0))
	typeSolver = 2;

      else {
	opserr << "eigen - unknown option specified " << argv[loc] << endln;
	return TCL_ERROR;
      }

      loc++;
    }


    // check argv[loc] for number of modes
    int numEigen;
    if ((Tcl_GetInt(interp, argv[loc], &numEigen) != TCL_OK) || numEigen < 0) {
      opserr << "WARNING eigen numModes?  - illegal numModes\n";    
      return TCL_ERROR;	
    }

    EigenAlgorithm *theEigenAlgo = 0;
    EigenSOE       *theEigenSOE = 0;
    AnalysisModel *theEigenModel = new AnalysisModel();

    // create the algorithm
    if (typeAlgo == 0) 
      theEigenAlgo = new FrequencyAlgo();
    else if (typeAlgo == 1) {
      theEigenAlgo = new StandardEigenAlgo();

      // temporarily will place here .. only solver that will work with standard
      SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver(); 
      theEigenSOE = new SymBandEigenSOE(*theEigenSolver, *theEigenModel);    

    }

    // again temporary until i rewrite these solvers.
    if (typeAlgo == 0) {

      // create the eigen system and solver
      if (typeSolver == 0) {
	SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver(); 
	theEigenSOE = new SymBandEigenSOE(*theEigenSolver, *theEigenModel);    
      } else if (typeSolver == 1) {
	SymArpackSolver *theEigenSolver = new SymArpackSolver(numEigen); 
	theEigenSOE = new SymArpackSOE(*theEigenSolver, *theEigenModel);    
      } else if (typeSolver == 2) {  
	BandArpackSolver *theEigenSolver = new BandArpackSolver(numEigen); 
	theEigenSOE = new BandArpackSOE(*theEigenSolver, *theEigenModel);    
      }      
    }
    // create the rest of components of an eigen analysis
    EigenIntegrator  *theEigenIntegrator = new EigenIntegrator();    
    RCM *theRCM = new RCM(false);	
    DOF_Numberer *theEigenNumberer = new DOF_Numberer(*theRCM);    	
    ConstraintHandler *theEigenHandler = new TransformationConstraintHandler();


    // create the eigen analysis
    theEigenAnalysis = new EigenAnalysis(theDomain,
					 *theEigenHandler,
					 *theEigenNumberer,
					 *theEigenModel,
					 *theEigenAlgo,
					 *theEigenSOE,
					 *theEigenIntegrator);

    int requiredDataSize = 20*numEigen;
    if (requiredDataSize > resDataSize) {
      if (resDataPtr != 0) {
	delete [] resDataPtr;
      }
      resDataPtr = new char[requiredDataSize];
      resDataSize = requiredDataSize;
    }

    for (int i=0; i<requiredDataSize; i++)
      resDataPtr[i] = '\n';

    // perfrom the eigen analysis & store the results with the interpreter

    if (theEigenAnalysis->analyze(numEigen) == 0) {
      //      char *eigenvalueS = new char[15 * numEigen];    
      const Vector &eigenvalues = theDomain.getEigenvalues();
      int cnt = 0;

      for (int i=0; i<numEigen; i++) {
	cnt += sprintf(&resDataPtr[cnt], "%.6e  ", eigenvalues[i]);
      }
      
      Tcl_SetResult(interp, resDataPtr, TCL_STATIC);
    }

    // finally invoke the destructor on the eigen analysis
    delete theEigenAnalysis;

    theEigenAnalysis = 0;
    
    return TCL_OK;
}








int 
videoPlayer(ClientData clientData, Tcl_Interp *interp, int argc, 
	    TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 5) {
	opserr << "WARNING want - video -window windowTitle? -file fileName?\n";
	return TCL_ERROR;
    }    

    TCL_Char *wTitle =0;
    TCL_Char *fName = 0;
    TCL_Char *imageName = 0;
    TCL_Char *offsetName = 0;

    int endMarker = 1;
    while (endMarker < (argc-1)) {
      if (strcmp(argv[endMarker],"-window") == 0) {
	wTitle = argv[endMarker+1];
	endMarker+=2;
      } else if (strcmp(argv[endMarker],"-file") == 0) {
	fName = argv[endMarker+1];
	endMarker+=2;
      } else if (strcmp(argv[endMarker],"-image") == 0) {
	imageName = argv[endMarker+1];
	endMarker += 2;
      } else if (strcmp(argv[endMarker],"-offset") == 0) {
	offsetName = argv[endMarker+1];
	endMarker += 2;
      }
      else {
	opserr << "WARNING unknown " << argv[endMarker] << 
	  " want - video -window windowTitle? -file fileName?\n";
				
	return TCL_ERROR;
      }
    }


#ifdef _NOGRAPHICS

#else    
    if (wTitle != 0 && fName != 0) {
      // delete the old video player if one exists
      if (theTclVideoPlayer != 0)
	delete theTclVideoPlayer;

      // create a new player
      theTclVideoPlayer = new TclVideoPlayer(wTitle, fName, imageName, interp, offsetName);
    }
    else
      return TCL_ERROR;
#endif
    return TCL_OK;
}



int 
removeObject(ClientData clientData, Tcl_Interp *interp, int argc, 
	     TCL_Char **argv)
{

    // make sure at least one other argument to contain type of system
      if (argc < 2) {
	opserr << "WARNING want - remove objectType?\n";
	return TCL_ERROR;
      }    

    int tag;
    if ((strcmp(argv[1],"element") == 0) || (strcmp(argv[1],"ele") == 0)) {
      if (argc < 3) {
	opserr << "WARNING want - remove element eleTag?\n";
	return TCL_ERROR;
      }    

      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING remove element tag? failed to read tag: " << argv[2] << endln;
	return TCL_ERROR;
      }      
      Element *theEle = theDomain.removeElement(tag);
      if (theEle != 0) {
	// we also have to remove any elemental loads from the domain
	LoadPatternIter &theLoadPatterns = theDomain.getLoadPatterns();
	LoadPattern *thePattern;
	
	// go through all load patterns
	while ((thePattern = theLoadPatterns()) != 0) {
	  ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
	  ElementalLoad *theLoad;

	  // go through all elemental loads in the pattern
	  while ((theLoad = theEleLoads()) != 0) {

	    // remove & destroy elemental from elemental load if there
	    // note - if last element in load, remove the load and delete it
	    
	    /* *****************
	    int numLoadsLeft = theLoad->removeElement(tag);
	    if (numLoadsLeft == 0) {
	      thePattern->removeElementalLoad(theLoad->getTag());
	      delete theLoad;
	    }
	    *********************/
	  }
	}

	// finally invoke the destructor on the element
	delete theEle;
      }
    }

    else if (strcmp(argv[1],"loadPattern") == 0) {
      if (argc < 3) {
	opserr << "WARNING want - remove loadPattern patternTag?\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING remove loadPattern tag? failed to read tag: " << argv[2] << endln;
	return TCL_ERROR;
      }      
      LoadPattern *thePattern = theDomain.removeLoadPattern(tag);
      if (thePattern != 0) {
	delete thePattern;
      }
    }


    else if (strcmp(argv[1],"node") == 0) {
      if (argc < 3) {
	opserr << "WARNING want - remove node nodeTag?\n";
	return TCL_ERROR;
      }    
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING remove node tag? failed to read tag: " << argv[2] << endln;
	return TCL_ERROR;
      }      
      Node *theNode = theDomain.removeNode(tag);
      if (theNode != 0) {
	delete theNode;
      }
    }


    else if (strcmp(argv[1],"recorders") == 0) {
      theDomain.removeRecorders();
    }

    //Boris Jeremic and Joey Yang -- UC Davis
    else if ((strcmp(argv[1],"SPconstraint") == 0) || (strcmp(argv[1],"sp") == 0)) {
      if (argc < 3) {
      	opserr << "WARNING want - remove SPconstraint spTag?\n";
      	return TCL_ERROR;
      }    
      if (argc == 3) {
	if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	  opserr << "WARNING remove loadPattern tag? failed to read tag: " << argv[2] << endln;
	  return TCL_ERROR;
	}      
      
	SP_Constraint *theSPconstraint = theDomain.removeSP_Constraint(tag);
	if (theSPconstraint != 0) {
	  delete theSPconstraint;
	}
      } else {
	int nodeTag, dofTag;
	if (Tcl_GetInt(interp, argv[2], &nodeTag) != TCL_OK) {
	  opserr << "WARNING remove loadPattern tag? failed to read node tag: " << argv[2] << endln;
	  return TCL_ERROR;
	}      
	if (Tcl_GetInt(interp, argv[3], &dofTag) != TCL_OK) {
	  opserr << "WARNING remove loadPattern tag? failed to read node tag: " << argv[2] << endln;
	  return TCL_ERROR;
	}      
	dofTag--;  // one for C++ indexing of dof

	SP_ConstraintIter &theSPs = theDomain.getSPs();
	SP_Constraint *theSP;

	ID theSPTags(0,12); int cnt=0;
	while ((theSP = theSPs()) != 0) {

	  int spNode = theSP->getNodeTag();
	  if (spNode == nodeTag) {
	    int spDofTag = theSP->getDOF_Number();
	    if (spDofTag == dofTag) {
	      theSPTags(cnt) = theSP->getTag();
	      cnt++;
	    }
	  }
	}

	for (int i=0; i<cnt; i++) {
	  SP_Constraint *theSPconstraint = theDomain.removeSP_Constraint(theSPTags(i));
	  if (theSPconstraint != 0) {
	    delete theSPconstraint;
	  }	
	}
      }
    }

#ifdef _RELIABILITY
// AddingSensitivity:BEGIN ///////////////////////////////////////
    else if (strcmp(argv[1],"randomVariablePositioner") == 0) {
		int rvPosTag;
		if (Tcl_GetInt(interp, argv[2], &rvPosTag) != TCL_OK) {
			opserr << "WARNING invalid input: rvPositionerTag \n";
			return TCL_ERROR;
		}
		ReliabilityDomain *theReliabilityDomain = theReliabilityBuilder->getReliabilityDomain();
		theReliabilityDomain->removeRandomVariablePositioner(rvPosTag);
	}
    else if (strcmp(argv[1],"performanceFunction") == 0) {
		int lsfTag;
		if (Tcl_GetInt(interp, argv[2], &lsfTag) != TCL_OK) {
			opserr << "WARNING invalid input: lsfTag \n";
			return TCL_ERROR;
		}
		ReliabilityDomain *theReliabilityDomain = theReliabilityBuilder->getReliabilityDomain();
		theReliabilityDomain->removePerformanceFunction(lsfTag);
	}
    else if (strcmp(argv[1],"sensitivityAlgorithm") == 0) {
		if (theSensitivityAlgorithm != 0) {
			theStaticAnalysis->setSensitivityAlgorithm(0);
			theSensitivityAlgorithm = 0;
		}
	}
// AddingSensitivity:END ///////////////////////////////////////
#endif

    else
      opserr << "WARNING remove element, loadPattern - only commands  available at the moment: " << endln;

    return TCL_OK;
}


int 
nodeDisp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    
    // make sure at least one other argument to contain type of system
    if (argc < 3) {
	opserr << "WARNING want - nodeDisp nodeTag? dof?\n";
	return TCL_ERROR;
   }    

    int tag, dof;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
    }        
    
    Node *theNode = theDomain.getNode(tag);
    double value = 0.0;
    if (theNode != 0) {
	const Vector &disp = theNode->getTrialDisp();
	if (disp.Size() >= dof && dof > 0) {
	    value = disp(dof-1); // -1 for OpenSees vs C indexing
	}
    }
    
    // now we copy the value to the tcl string that is returned
    sprintf(interp->result,"%35.20f",value);
    
	
    return TCL_OK;
}

int 
nodeCoord(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - nodeCoord nodeTag? dim?\n";
    return TCL_ERROR;
  }    
  
  int tag, dim;
  
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeCoord nodeTag? dim? - could not read nodeTag? \n";
    return TCL_ERROR;	        
  }    
  if (strcmp(argv[2],"X") == 0 || strcmp(argv[2],"x") == 0 ||
      strcmp(argv[2],"1") == 0)
    dim = 1;
  else if (strcmp(argv[2],"Y") == 0 || strcmp(argv[2],"y") == 0 ||
      strcmp(argv[2],"2") == 0)
    dim = 2;
  else if (strcmp(argv[2],"Z") == 0 || strcmp(argv[2],"z") == 0 ||
      strcmp(argv[2],"3") == 0)
    dim = 3;
  else {
    opserr << "WARNING nodeCoord nodeTag? dim? - could not read dim? \n";
    return TCL_ERROR;	        
  }        
  
  Node *theNode = theDomain.getNode(tag);
  double value = 0.0;
  if (theNode != 0) {
    const Vector &coords = theNode->getCrds();
    if (coords.Size() >= dim && dim > 0) {
      value = coords(dim-1); // -1 for OpenSees vs C indexing
    }
  }
  
  // now we copy the value to the tcl string that is returned
  sprintf(interp->result,"%35.20f",value);
  
  return TCL_OK;
}

int 
nodeBounds(ClientData clientData, Tcl_Interp *interp, int argc, 
	   TCL_Char **argv)
{
  int requiredDataSize = 20*6;
  if (requiredDataSize > resDataSize) {
    if (resDataPtr != 0) {
      delete [] resDataPtr;
    }
    resDataPtr = new char[requiredDataSize];
    resDataSize = requiredDataSize;
  }
  
  for (int i=0; i<requiredDataSize; i++)
    resDataPtr[i] = '\n';
  
  const Vector &bounds = theDomain.getPhysicalBounds();
  
  int cnt = 0;
  for (int i=0; i<6; i++) {
    cnt += sprintf(&resDataPtr[cnt], "%.6e  ", bounds(i));
  }
  
  Tcl_SetResult(interp, resDataPtr, TCL_STATIC);
  
  return TCL_OK;
}


// AddingSensitivity:BEGIN ////////////////////////////////////
int 
nodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    
    // make sure at least one other argument to contain type of system
    if (argc < 3) {
	interp->result = "WARNING want - nodeVel nodeTag? dof?\n";
	return TCL_ERROR;
   }    

    int tag, dof;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? - could not read nodeTag? ";
	return TCL_ERROR;	        
    }    
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? - could not read dof? ";
	return TCL_ERROR;	        
    }        
    
    Node *theNode = theDomain.getNode(tag);
    double value = 0.0;
    if (theNode != 0) {
	const Vector &vel = theNode->getTrialVel();
	if (vel.Size() >= dof && dof > 0) {
	    value = vel(dof-1); // -1 for OpenSees vs C indexing
	}
    }
    
    // now we copy the value to the tcl string that is returned
    sprintf(interp->result,"%35.20f",value);
    
	
    return TCL_OK;
}


int 
sensNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

    // make sure at least one other argument to contain type of system
    if (argc < 4) {
	interp->result = "WARNING want - sensNodeDisp nodeTag? dof?\n";
	return TCL_ERROR;
   }    

    int tag, dof, gradNum;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? gradNum?- could not read nodeTag? ";
	return TCL_ERROR;	        
    }    
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? gradNum?- could not read dof? ";
	return TCL_ERROR;	        
    }        
    if (Tcl_GetInt(interp, argv[3], &gradNum) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? gradNum?- could not read dof? ";
	return TCL_ERROR;	        
    }        
    
    Node *theNode = theDomain.getNode(tag);
	double value = theNode->getDispSensitivity(dof,gradNum);
    
    // copy the value to the tcl string that is returned
    sprintf(interp->result,"%35.20f",value);
	
    return TCL_OK;
}
int 
sensNodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

    // make sure at least one other argument to contain type of system
    if (argc < 4) {
	interp->result = "WARNING want - sensNodeDisp nodeTag? dof?\n";
	return TCL_ERROR;
   }    

    int tag, dof, gradNum;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? gradNum?- could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? gradNum?- could not read dof? \n";
	return TCL_ERROR;	        
    }        
    if (Tcl_GetInt(interp, argv[3], &gradNum) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? gradNum?- could not read dof? \n";
	return TCL_ERROR;	        
    }        
    
    Node *theNode = theDomain.getNode(tag);
	double value = theNode->getVelSensitivity(dof,gradNum);
    
    // copy the value to the tcl string that is returned
    sprintf(interp->result,"%35.20f",value);
	
    return TCL_OK;
}



int 
computeGradients(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// Comment to the developer:
	// This Tcl command is meant to be called only for
	// path-independent problems.  In such cases the 
	// gradients can be computed AFTER the complete
	// structural analysis is completed.
	// No error messages is returned if the user tries to invoke this command
	// during a path-dependent analysis.  The reason is that the reliability
	// analysis will call this method BOTH for path dependent and independent problems. 

#ifdef _RELIABILITY
	if (theSensitivityAlgorithm->shouldComputeAtEachStep()) {
	}
	else {

		theSensitivityAlgorithm->computeSensitivities();
	}
#endif
    return TCL_OK;
}
// AddingSensitivity:END //////////////////////////////////////


int 
startTimer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theTimer == 0)
    theTimer = new Timer();
  
  theTimer->start();
  return TCL_OK;
}

int 
stopTimer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theTimer == 0)
    return TCL_OK;
  
  theTimer->pause();
  opserr << *theTimer;
  return TCL_OK;
}

int 
rayleighDamping(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 5) { 
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - not enough arguments to command\n";
    return TCL_ERROR;
  }
  double alphaM, betaK, betaK0, betaKc;
  if (Tcl_GetDouble(interp, argv[1], &alphaM) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read alphaM? \n";
    return TCL_ERROR;	        
  }    
  if (Tcl_GetDouble(interp, argv[2], &betaK) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaK? \n";
    return TCL_ERROR;	        
  }        
  if (Tcl_GetDouble(interp, argv[3], &betaK0) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaK0? \n";
    return TCL_ERROR;	        
  }        
  if (Tcl_GetDouble(interp, argv[4], &betaKc) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaKc? \n";
    return TCL_ERROR;	        
  }        
  
  theDomain.setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  return TCL_OK;
}


extern int
TclAddMeshRegion(ClientData clientData, Tcl_Interp *interp, int argc, 
	     TCL_Char **argv, Domain &theDomain);

int 
addRegion(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  return TclAddMeshRegion(clientData, interp, argc, argv, theDomain);
}

int 
logFile(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

  if (argc < 2) { 
    opserr << "WARNING logFile fileName? - no filename supplied\n";
    return TCL_ERROR;
  }
  openMode mode = OVERWRITE;
  if (argc >= 3) 
    if (strcmp(argv[2],"-append") == 0) 
      mode = APPEND;

  if (opserr.setFile(argv[1], mode) < 0) 
    opserr << "WARNING logFile " << argv[1] << " failed to set the file\n";

  
  simulationInfo.addWriteFile(argv[1]);

  return TCL_OK;
}

int 
exit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  Tcl_Finalize();
  return TCL_OK;
}


