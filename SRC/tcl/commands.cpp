/* ****************************************************************** **
**    OpenSees System for Earthquake Engineering Simulation           **
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
                                                                        
// $Revision$
// $Date$
// $URL$
                                                                        
                                                                        
// Written: fmk 
// Created: 04/98
//
// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified,
// see tkAppInit.C for command names.

#include <classTags.h>

#include <DOF_Group.h>

#ifdef _PARALLEL_PROCESSING
#include <mpi.h>
#elif _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif

extern "C" {
#include <tcl.h>
}

#include <OPS_Globals.h>
#include <TclModelBuilder.h>
#include <Matrix.h>
#include <iostream>
#include <set>
#include <algorithm>

extern void OPS_clearAllUniaxialMaterial(void);
extern void OPS_clearAllNDMaterial(void);
extern void OPS_clearAllSectionForceDeformation(void);

extern void OPS_clearAllHystereticBackbone(void);
extern void OPS_clearAllStiffnessDegradation(void);
extern void OPS_clearAllStrengthDegradation(void);
extern void OPS_clearAllUnloadingRule(void);

// the following is a little kludgy but it works!
#ifdef _USING_STL_STREAMS

#include <iomanip>
using std::ios;
#include <iostream>
using std::ofstream;
 
#else

#include <StandardStream.h>
#include <FileStream.h>
#include <DummyStream.h>

bool OPS_suppressOpenSeesOutput = false;
StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <elementAPI.h>
extern "C" int         OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp * interp, int cArg, int mArg, TCL_Char * *argv, Domain * domain);

#include <packages.h>

#include <FEM_ObjectBrokerAllClasses.h>

#include <Timer.h>
#include <ModelBuilder.h>
#include "commands.h"

// domain
 #ifdef _PARALLEL_PROCESSING
#include <PartitionedDomain.h>
#else
#include <Domain.h>
#endif

#include <Information.h>
#include <Element.h>
#include <Node.h>
#include <ElementIter.h>
#include <NodeIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>
#include <ParameterIter.h>
#include <SP_Constraint.h> //Joey UC Davis
#include <SP_ConstraintIter.h> //Joey UC Davis
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <InitialStateParameter.h>
#include <ElementStateParameter.h>
#include <Pressure_Constraint.h>

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
#include <NormDispAndUnbalance.h>
#include <NormDispOrUnbalance.h>
#include <CTestPFEM.h>

// soln algorithms
#include <Linear.h>
#include <NewtonRaphson.h>
#include <NewtonLineSearch.h>
#include <ModifiedNewton.h>
#include <Broyden.h>
#include <BFGS.h>
#include <KrylovNewton.h>
#include <PeriodicNewton.h>
#include <AcceleratedNewton.h>
#include <ExpressNewton.h>

// accelerators
#include <RaphsonAccelerator.h>
#include <PeriodicAccelerator.h>
#include <KrylovAccelerator.h>
#include <SecantAccelerator1.h>
#include <SecantAccelerator2.h>
#include <SecantAccelerator3.h>
//#include <MillerAccelerator.h>

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
#include <StagedLoadControl.h>
#include <ArcLength.h>
#include <ArcLength1.h>
#include <HSConstraint.h>
#include <MinUnbalDispNorm.h>
#include <DisplacementControl.h>
#include <EQPath.h>

#include <PFEMIntegrator.h>
#include <Integrator.h>//Abbas

//  recorders
#include <Recorder.h> //SAJalali

extern void *OPS_NewtonRaphsonAlgorithm(void);
extern void *OPS_ModifiedNewton(void);
extern void *OPS_NewtonHallM(void);

extern void *OPS_Newmark(void);
extern void *OPS_StagedNewmark(void);
extern void *OPS_GimmeMCK(void);
extern void *OPS_AlphaOS(void);
extern void *OPS_AlphaOS_TP(void);
extern void *OPS_AlphaOSGeneralized(void);
extern void *OPS_AlphaOSGeneralized_TP(void);
extern void *OPS_CentralDifference(void);
extern void *OPS_CentralDifferenceAlternative(void);
extern void *OPS_CentralDifferenceNoDamping(void);
extern void *OPS_Collocation(void);
extern void *OPS_CollocationHSFixedNumIter(void);
extern void *OPS_CollocationHSIncrLimit(void);
extern void *OPS_CollocationHSIncrReduct(void);
extern void *OPS_GeneralizedAlpha(void);
extern void *OPS_HHT(void);
extern void *OPS_HHT_TP(void);
extern void *OPS_HHTExplicit(void);
extern void *OPS_HHTExplicit_TP(void);
extern void *OPS_HHTGeneralized(void);
extern void *OPS_HHTGeneralized_TP(void);
extern void *OPS_HHTGeneralizedExplicit(void);
extern void *OPS_HHTGeneralizedExplicit_TP(void);
extern void *OPS_HHTHSFixedNumIter(void);
extern void *OPS_HHTHSFixedNumIter_TP(void);
extern void *OPS_HHTHSIncrLimit(void);
extern void *OPS_HHTHSIncrLimit_TP(void);
extern void *OPS_HHTHSIncrReduct(void);
extern void *OPS_HHTHSIncrReduct_TP(void);
extern void *OPS_KRAlphaExplicit(void);
extern void *OPS_KRAlphaExplicit_TP(void);
extern void *OPS_NewmarkExplicit(void);
extern void *OPS_NewmarkHSFixedNumIter(void);
extern void *OPS_NewmarkHSIncrLimit(void);
extern void *OPS_NewmarkHSIncrReduct(void);
extern void *OPS_WilsonTheta(void);

// for response spectrum analysis
extern void OPS_DomainModalProperties(void);
extern void OPS_ResponseSpectrumAnalysis(void);

#include <Newmark.h>
#include <StagedNewmark.h>
#include <TRBDF2.h>
#include <TRBDF3.h>
#include <Newmark1.h>
#include <Houbolt.h>
#include <ParkLMS3.h>
#include <BackwardEuler.h>

// analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>
#include <PFEMAnalysis.h>

// system of eqn and solvers
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>

#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>

#include <ConjugateGradientSolver.h>

#ifdef _ITPACK
//#include <ItpackLinSOE.h>
//#include <ItpackLinSolver.h>
#endif

#include <FullGenLinSOE.h>
#include <FullGenLinLapackSolver.h>

#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectSolver.h>
#include <DiagonalSOE.h>
#include <DiagonalDirectSolver.h>

#include <SProfileSPDLinSolver.h>
#include <SProfileSPDLinSOE.h>

// #include <ProfileSPDLinDirectBlockSolver.h>
// #include <ProfileSPDLinDirectThreadSolver.h>
// #include <ProfileSPDLinDirectSkypackSolver.h>
// #include <BandSPDLinThreadSolver.h>

#include <SparseGenColLinSOE.h>
#include <PFEMSolver.h>
#include <PFEMSolver_Umfpack.h>
#include <PFEMLinSOE.h>
#include <PFEMCompressibleSolver.h>
#include <PFEMCompressibleLinSOE.h>
#ifdef _MUMPS
#include <PFEMSolver_Mumps.h>
#include <PFEMCompressibleSolver_Mumps.h>
#endif

#ifdef _THREADS
#include <ThreadedSuperLU.h>
#else
#include <SuperLU.h>
#endif

#ifdef _CUSP
#include <CuSPSolver.h>
#endif

#ifdef _CULAS4
#include <CulaSparseSolverS4.h>
#endif

#ifdef _CULAS5
#include <CulaSparseSolverS5.h>
#endif

#ifdef _MUMPS
#ifdef _PARALLEL_PROCESSING
#include <MumpsParallelSOE.h>
#include <MumpsParallelSolver.h>
#elif _PARALLEL_INTERPRETERS
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

#include <SparseGenRowLinSOE.h>
#include <SymSparseLinSOE.h>
#include <SymSparseLinSolver.h>
#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <EigenSOE.h>
#include <EigenSolver.h>
#include <ArpackSOE.h>
#include <ArpackSolver.h>
#include <SymArpackSOE.h>
#include <SymArpackSolver.h>
#include <BandArpackSOE.h>
#include <BandArpackSolver.h>
#include <SymBandEigenSOE.h>
#include <SymBandEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <FullGenEigenSolver.h>

#ifdef _CUDA
#include <BandGenLinSOE_Single.h>
#include <BandGenLinLapackSolver_Single.h>
#endif


// graph
#include <RCM.h>
#include <AMDNumberer.h>

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
// #include <SensitivityIntegrator.h>
// #include <StaticSensitivityIntegrator.h>
//#include <DynamicSensitivityIntegrator.h>
// #include <NewmarkSensitivityIntegrator.h>
// #include <NewNewmarkSensitivityIntegrator.h>
// #include <NewStaticSensitivityIntegrator.h>
// #include <PFEMSensitivityIntegrator.h>
//#include <OrigSensitivityAlgorithm.h>
//#include <NewSensitivityAlgorithm.h>
//#include <ReliabilityStaticAnalysis.h>
//#include <ReliabilityDirectIntegrationAnalysis.h>
// AddingSensitivity:END /////////////////////////////////////////////////
#include <TclReliabilityBuilder.h>

int reliability(ClientData, Tcl_Interp *, int, TCL_Char **);
int wipeReliability(ClientData, Tcl_Interp *, int, TCL_Char **);
int optimization(ClientData, Tcl_Interp *, int, TCL_Char **);  //Quan  (2)

#endif

const char * getInterpPWD(Tcl_Interp *interp);

#include <XmlFileStream.h>


#include <Response.h>

ModelBuilder *theBuilder =0;

// some global variables 
#ifdef _PARALLEL_PROCESSING

#include <DistributedDisplacementControl.h>
#include <ShadowSubdomain.h>
#include <Metis.h>
#include <ShedHeaviest.h>
#include <DomainPartitioner.h>
#include <GraphPartitioner.h>
#include <FEM_ObjectBrokerAllClasses.h>
#include <Subdomain.h>
#include <SubdomainIter.h>
#include <MachineBroker.h>
#include <MPIDiagonalSOE.h>
#include <MPIDiagonalSolver.h>

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

//MachineBroker *theMachineBroker = 0;
PartitionedDomain theDomain;
int OPS_PARALLEL_PROCESSING =0;
int OPS_NUM_SUBDOMAINS      =0;
bool OPS_PARTITIONED        =false;
bool OPS_USING_MAIN_DOMAIN  = false;
int OPS_MAIN_DOMAIN_PARTITION_ID =0;

DomainPartitioner *OPS_DOMAIN_PARTITIONER =0;
GraphPartitioner  *OPS_GRAPH_PARTITIONER =0;
LoadBalancer      *OPS_BALANCER = 0;
FEM_ObjectBroker  *OPS_OBJECT_BROKER =0;
MachineBroker     *OPS_MACHINE =0;
Channel          **OPS_theChannels = 0;

bool setMPIDSOEFlag = false;

#elif _PARALLEL_INTERPRETERS

bool setMPIDSOEFlag = false;

// parallel analysis
#include <ParallelNumberer.h>
#include <DistributedDisplacementControl.h>

//  parallel soe & solvers
#include <DistributedBandSPDLinSOE.h>
#include <DistributedSparseGenColLinSOE.h>
#include <DistributedSparseGenRowLinSOE.h>


#include <DistributedBandGenLinSOE.h>
#include <DistributedDiagonalSOE.h>
#include <DistributedDiagonalSolver.h>
#include <MPIDiagonalSOE.h>
#include <MPIDiagonalSolver.h>
#define MPIPP_H
#include <DistributedSuperLU.h>
#include <DistributedProfileSPDLinSOE.h>

Domain theDomain;

#else

Domain theDomain;

#endif

#include <MachineBroker.h>

MachineBroker *theMachineBroker =0;
Channel **theChannels =0;
int numChannels =0;
int OPS_rank =0;
int OPS_np =0;


typedef struct parameterValues {
  char *value;
  struct parameterValues *next;
} OpenSeesTcl_ParameterValues;


typedef struct parameter {
  char *name;
  OpenSeesTcl_ParameterValues *values;
  struct parameter *next;
} OpenSeesTcl_Parameter;


static OpenSeesTcl_Parameter *theParameters = NULL;
static OpenSeesTcl_Parameter *endParameters = NULL;

static int numParam = 0;
static char **paramNames =0;
static char **paramValues =0;

AnalysisModel *theAnalysisModel =0;
EquiSolnAlgo *theAlgorithm =0;
ConstraintHandler *theHandler =0;
DOF_Numberer *theNumberer =0;
LinearSOE *theSOE =0;
EigenSOE *theEigenSOE =0;
StaticAnalysis *theStaticAnalysis = 0;
DirectIntegrationAnalysis *theTransientAnalysis = 0;
VariableTimeStepDirectIntegrationAnalysis *theVariableTimeStepTransientAnalysis = 0;
int numEigen = 0;

static PFEMAnalysis* thePFEMAnalysis = 0;

// AddingSensitivity:BEGIN /////////////////////////////////////////////
#ifdef _RELIABILITY
static TclReliabilityBuilder *theReliabilityBuilder = 0;

Integrator *theSensitivityAlgorithm = 0;
Integrator *theSensitivityIntegrator = 0;
//FMK RELIABILITY ReliabilityStaticAnalysis *theReliabilityStaticAnalysis = 0;
//FMK RELIABILITY ReliabilityDirectIntegrationAnalysis *theReliabilityTransientAnalysis = 0;

// static NewmarkSensitivityIntegrator *theNSI = 0;
// static NewNewmarkSensitivityIntegrator *theNNSI = 0;
// static PFEMSensitivityIntegrator* thePFEMSI = 0;

//static SensitivityIntegrator *theSensitivityIntegrator = 0;
//static NewmarkSensitivityIntegrator *theNSI = 0;

#include <TclOptimizationBuilder.h>
static TclOptimizationBuilder *theOptimizationBuilder = 0;   // Quan March 2010 (3)

#endif
// AddingSensitivity:END ///////////////////////////////////////////////


StaticIntegrator *theStaticIntegrator =0;
TransientIntegrator *theTransientIntegrator =0;
ConvergenceTest *theTest =0;
bool builtModel = false;

static char *resDataPtr = 0;
static int resDataSize = 0;
static Timer *theTimer = 0;

#include <FileStream.h>
#include <SimulationInformation.h>
SimulationInformation simulationInfo;
SimulationInformation *theSimulationInfoPtr = 0;

char *simulationInfoOutputFilename = 0;


FE_Datastore *theDatabase  =0;
FEM_ObjectBrokerAllClasses theBroker;

// init the global variabled defined in OPS_Globals.h
// double        ops_Dt = 1.0;
// Element    *ops_TheActiveElement = 0;
// bool          ops_InitialStateAnalysis = false; // McGann, U.Washington

#ifdef _NOGRAPHICS

#else
TclVideoPlayer *theTclVideoPlayer =0;
#endif

// g3AppInit() is the method called by tkAppInit() when the
// interpreter is being set up .. this is where all the
// commands defined in this file are registered with the interpreter.


int 
printModelGID(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
printA(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
printB(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
setPrecision(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
logFile(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
version(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
getPID(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
getNP(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
opsBarrier(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
domainChange(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
record(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
opsSend(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
opsRecv(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
opsPartition(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);


int
peerNGA(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
defaultUnits(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
stripOpenSeesXML(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
setParameter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

//extern 
int OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

extern int myCommands(Tcl_Interp *interp);

//extern "C" int Tcl_InterpObjCmd(ClientData clientData,  
//			Tcl_Interp *interp, 
//		int objc, 
//	Tcl_Obj *const objv[]);

int
convertBinaryToText(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
convertTextToBinary(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
maxOpenFiles(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);



// pointer for old putsCommand

static Tcl_ObjCmdProc *Tcl_putsCommand = 0;

//
// revised puts command to send to cerr!
//

int OpenSees_putsCommand(ClientData dummy,  Tcl_Interp *interp, int objc, Tcl_Obj *const objv[])
{
    Tcl_Channel chan;           /* The channel to puts on. */
    Tcl_Obj *string;            /* String to write. */
    Tcl_Obj *chanObjPtr = NULL; /* channel object. */
    int newline;                /* Add a newline at end? */

    switch (objc) {
    case 2: /* [puts $x] */
        string = objv[1];
        newline = 1;
        break;

    case 3: /* [puts -nonewline $x] or [puts $chan $x] */
        if (strcmp(Tcl_GetString(objv[1]), "-nonewline") == 0) {
            newline = 0;
        } else {
            newline = 1;
            chanObjPtr = objv[1];
        }
        string = objv[2];
        break;

    case 4: /* [puts -nonewline $chan $x] or [puts $chan $x nonewline] */
        newline = 0;
        if (strcmp(Tcl_GetString(objv[1]), "-nonewline") == 0) {
            chanObjPtr = objv[2];
            string = objv[3];
            break;
        } else if (strcmp(Tcl_GetString(objv[3]), "nonewline") == 0) {
	  /*
             * The code below provides backwards compatibility with an old
             * form of the command that is no longer recommended or
             * documented. See also [Bug #3151675]. Will be removed in Tcl 9,
             * maybe even earlier.
             */

            chanObjPtr = objv[1];
            string = objv[2];
            break;
        }
        /* Fall through */
    default:
        /* [puts] or [puts some bad number of arguments...] */
        Tcl_WrongNumArgs(interp, 1, objv, "?-nonewline? ?channelId? string");
        return TCL_ERROR;
    }

    if (chanObjPtr == NULL) {
        if (newline == 0)
            opserr << Tcl_GetString(string);
        else
            opserr << Tcl_GetString(string) << endln;
        return TCL_OK;
    } else {
       if (Tcl_putsCommand != 0) {
           return Tcl_putsCommand(dummy, interp, objc, objv);
       } else {
           std::cerr << "MEARD!  commands.cpp .. old puts command not found or set!\n";
           return TCL_ERROR;
    }
    return TCL_OK;
  }
}


int Tcl_InterpOpenSeesObjCmd(ClientData clientData,  Tcl_Interp *interp, int objc, Tcl_Obj *CONST objv[])
{
  int index;
  static TCL_Char *options[] = {
    "alias",	"aliases",	"create",	"delete", 
    "eval",		"exists",	"expose",	"hide", 
    "hidden",	"issafe",	"invokehidden",	"marktrusted", 
    "recursionlimit",		"slaves",	"share",
    "target",	"transfer",
    NULL
  };
  enum option {
    OPT_ALIAS,	OPT_ALIASES,	OPT_CREATE,	OPT_DELETE,
    OPT_EVAL,	OPT_EXISTS,	OPT_EXPOSE,	OPT_HIDE,
    OPT_HIDDEN,	OPT_ISSAFE,	OPT_INVOKEHID,	OPT_MARKTRUSTED,
    OPT_RECLIMIT,			OPT_SLAVES,	OPT_SHARE,
    OPT_TARGET,	OPT_TRANSFER
  };

  int ok = TCL_OK;
  //int ok = Tcl_InterpObjCmd(clientData, interp, objc, objv);
  //if (ok != TCL_OK) 
  //return ok;
  
  if (Tcl_GetIndexFromObj(interp, objv[1], options, "option", 0, &index) != TCL_OK) {
    return TCL_ERROR;
  }
  
  switch ((enum option) index) {
  case OPT_CREATE: {
    TCL_Char *theInterpreterName = Tcl_GetStringResult(interp);
    Tcl_Interp *secondaryInterp = Tcl_GetSlave(interp, theInterpreterName);
    ok = OpenSeesAppInit(secondaryInterp);
    return ok;
    break;
  }
  default:
    return ok;
  }
  
  return ok;
}


int OpenSeesAppInit(Tcl_Interp *interp) {

  ops_TheActiveDomain = &theDomain;

  //
  // redo puts command so we can capture puts into std:cerr
  //

  if (OPS_suppressOpenSeesOutput == false) {
    // get a handle on puts procedure
    Tcl_CmdInfo putsCommandInfo;
    Tcl_GetCommandInfo(interp, "puts", &putsCommandInfo);
    Tcl_putsCommand = putsCommandInfo.objProc;
    // if handle, use ouur procedure as opposed to theirs
    if (Tcl_putsCommand != 0) {
      Tcl_CreateObjCommand(interp, "oldputs", Tcl_putsCommand, NULL, NULL);
      Tcl_CreateObjCommand(interp, "puts", OpenSees_putsCommand, NULL, NULL);
    }
  }

  theSimulationInfoPtr = &simulationInfo;
    
#ifndef _LINUX  
    opserr.setFloatField(SCIENTIFIC);
    opserr.setFloatField(FIXEDD);
#endif
	
    //Tcl_CreateObjCommand(interp, "interp", Tcl_InterpOpenSeesObjCmd, NULL, NULL);
	
	Tcl_CreateCommand(interp, "recorderValue", &OPS_recorderValue,
		(ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); //by SAJalali

	Tcl_CreateObjCommand(interp, "pset", &OPS_SetObjCmd,
			 (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
	
    Tcl_CreateObjCommand(interp, "source", &OPS_SourceCmd,
			 (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 

    Tcl_CreateCommand(interp, "getNDM", &getNDM,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
    Tcl_CreateCommand(interp, "getNDF", &getNDF,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

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

    Tcl_CreateCommand(interp, "setCreep", &setCreep,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "setTime", &setTime,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);     
    Tcl_CreateCommand(interp, "getTime", &getTime,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "getLoadFactor", &getLoadFactor,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
		
    Tcl_CreateCommand(interp, "build", &buildModel,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "analyze", &analyzeModel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "print", &printModel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "printModel", &printModel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "printA", &printA, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "printB", &printB, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    // Talledo Start 
    Tcl_CreateCommand(interp, "printGID", &printModelGID,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    // Talledo End
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
    Tcl_CreateCommand(interp, "testNorms", &getCTestNorms, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "testIter", &getCTestIter, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    

    Tcl_CreateCommand(interp, "integrator", &specifyIntegrator, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "recorder", &addRecorder, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "algorithmRecorder", &addAlgoRecorder, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "database", &addDatabase, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "eigen", &eigenAnalysis, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "modalProperties", &modalProperties,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
    Tcl_CreateCommand(interp, "responseSpectrum", &responseSpectrum,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
    Tcl_CreateCommand(interp, "video", &videoPlayer, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "remove", &removeObject, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       

    Tcl_CreateCommand(interp, "eleForce", &eleForce, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "localForce", &localForce, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);          
    Tcl_CreateCommand(interp, "eleDynamicalForce", &eleDynamicalForce, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "eleResponse", &eleResponse, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeDisp", &nodeDisp, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "setNodeDisp", &setNodeDisp, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);            
    Tcl_CreateCommand(interp, "nodeReaction", &nodeReaction, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeUnbalance", &nodeUnbalance, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeEigenvector", &nodeEigenvector, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeVel", &nodeVel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "setNodeVel", &setNodeVel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeAccel", &nodeAccel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
    Tcl_CreateCommand(interp, "setNodeAccel", &setNodeAccel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);             
    Tcl_CreateCommand(interp, "nodeResponse", &nodeResponse, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "reactions", &calculateNodalReactions, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "nodeDOFs", &nodeDOFs, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "nodeCoord", &nodeCoord, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);      
    Tcl_CreateCommand(interp, "setNodeCoord", &setNodeCoord, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "updateElementDomain", &updateElementDomain, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "eleType", &eleType,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
    Tcl_CreateCommand(interp, "eleNodes", &eleNodes, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);            
    Tcl_CreateCommand(interp, "nodeMass", &nodeMass, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);            
    Tcl_CreateCommand(interp, "nodePressure", &nodePressure, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);   
    Tcl_CreateCommand(interp, "nodeBounds", &nodeBounds, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "start", &startTimer, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "stop", &stopTimer, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "rayleigh", &rayleighDamping, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "modalDamping", &modalDamping, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "modalDampingQ", &modalDampingQ, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "setElementRayleighDampingFactors", 
		      &setElementRayleighDampingFactors, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "region", &addRegion, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "logFile", &logFile, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
    Tcl_CreateCommand(interp, "setPrecision", &setPrecision, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
    Tcl_CreateCommand(interp, "exit", &OpenSeesExit, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "quit", &OpenSeesExit, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "findNodeWithID", &findID, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       

    Tcl_CreateCommand(interp, "getNP", &getNP, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "getPID", &getPID, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "barrier", &opsBarrier, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "send", &opsSend, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "recv", &opsRecv, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "partition", &opsPartition, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "searchPeerNGA", &peerNGA, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    Tcl_CreateCommand(interp, "domainChange",  &domainChange,(ClientData)NULL, NULL);

    Tcl_CreateCommand(interp, "record",  &record,(ClientData)NULL, NULL);

    Tcl_CreateCommand(interp, "defaultUnits", &defaultUnits,(ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "stripXML", &stripOpenSeesXML,(ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "convertBinaryToText", &convertBinaryToText,(ClientData)NULL, NULL);
    Tcl_CreateCommand(interp, "convertTextToBinary", &convertTextToBinary,(ClientData)NULL, NULL);

    Tcl_CreateCommand(interp, "getEleTags", &getEleTags, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "getNodeTags", &getNodeTags, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "getParamTags", &getParamTags, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "getParamValue", &getParamValue, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
              
    Tcl_CreateCommand(interp, "fixedNodes", &fixedNodes,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
    Tcl_CreateCommand(interp, "fixedDOFs", &fixedDOFs,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
    Tcl_CreateCommand(interp, "constrainedNodes", &constrainedNodes,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
    Tcl_CreateCommand(interp, "constrainedDOFs", &constrainedDOFs,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
    Tcl_CreateCommand(interp, "retainedNodes", &retainedNodes,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);
    Tcl_CreateCommand(interp, "retainedDOFs", &retainedDOFs,
        (ClientData)NULL, (Tcl_CmdDeleteProc*)NULL);

    Tcl_CreateCommand(interp, "sdfResponse", &sdfResponse, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  

    Tcl_CreateCommand(interp, "sectionForce", &sectionForce, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "sectionDeformation", &sectionDeformation, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "sectionStiffness", &sectionStiffness, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "sectionFlexibility", &sectionFlexibility, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "sectionLocation", &sectionLocation, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "sectionWeight", &sectionWeight, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "basicDeformation", &basicDeformation, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "basicForce", &basicForce, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "basicStiffness", &basicStiffness, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  

    // command added for initial state analysis for nDMaterials
    // Chris McGann, U.Washington
    Tcl_CreateCommand(interp, "InitialStateAnalysis", &InitialStateAnalysis,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);

    Tcl_CreateCommand(interp, "totalCPU", &totalCPU, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "solveCPU", &solveCPU, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "accelCPU", &accelCPU, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "numFact", &numFact, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "numIter", &numIter, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "systemSize", &systemSize, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  
    Tcl_CreateCommand(interp, "version", &version, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  

    Tcl_CreateCommand(interp, "setParameter", &setParameter, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);  

    Tcl_CreateCommand(interp, "setMaxOpenFiles", &maxOpenFiles, 
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
     Tcl_CreateCommand(interp, "sensLambda", &sensLambda, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); //Abbas
     Tcl_CreateCommand(interp, "sensNodeVel", &sensNodeVel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "sensNodeAccel", &sensNodeAccel, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "sensSectionForce", &sensSectionForce, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);       
    Tcl_CreateCommand(interp, "sensNodePressure", &sensNodePressure, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);


    theSensitivityAlgorithm =0;
    theSensitivityIntegrator =0;
	//FMK RELIABILITY theReliabilityStaticAnalysis =0;
    //FMK RELIABILITY theReliabilityTransientAnalysis =0;    
    // AddingSensitivity:END //////////////////////////////////

    theOptimizationBuilder = 0;
    
    // --- Quan March 2010  (4)
    Tcl_CreateCommand(interp, "optimization", &optimization, 
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL); 
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
OPS_SetObjCmd(ClientData clientData, Tcl_Interp *interp, int objc, Tcl_Obj * const objv[])
{

    if (objc > 2)
    simulationInfo.addParameter(Tcl_GetString(objv[1]), Tcl_GetString(objv[2]));

    Tcl_Obj *varValueObj;
    
    if (objc == 2) {
      varValueObj = Tcl_ObjGetVar2(interp, objv[1], NULL,TCL_LEAVE_ERR_MSG);
      if (varValueObj == NULL) {
	return TCL_ERROR;
      }
      Tcl_SetObjResult(interp, varValueObj);
      return TCL_OK;
    } else if (objc == 3) {
      varValueObj = Tcl_ObjSetVar2(interp, objv[1], NULL, objv[2],
				   TCL_LEAVE_ERR_MSG);
      if (varValueObj == NULL) {
	return TCL_ERROR;
      }
      Tcl_SetObjResult(interp, varValueObj);
      return TCL_OK;
    } else {
      Tcl_WrongNumArgs(interp, 1, objv, "varName ?newValue?");
      return TCL_ERROR;
    }

    //    Tcl_SetObjCmd(clientData, interp, objc, objv);
    return 0;
}

int
OPS_SourceCmd(
    ClientData dummy,		/* Not used. */
    Tcl_Interp *interp,		/* Current interpreter. */
    int objc,			/* Number of arguments. */
    Tcl_Obj *CONST objv[])	/* Argument objects. */
{
    CONST char *encodingName = NULL;
    Tcl_Obj *fileName;
    
    if (objc != 2 && objc !=4) {
	Tcl_WrongNumArgs(interp, 1, objv, "?-encoding name? fileName");
	return TCL_ERROR;
    }

    fileName = objv[objc-1];

    if (objc == 4) {
	static CONST char *options[] = {
	    "-encoding", NULL
	};
	int index;

	if (TCL_ERROR == Tcl_GetIndexFromObj(interp, objv[1], options,
		"option", TCL_EXACT, &index)) {
	    return TCL_ERROR;
	}
	encodingName = Tcl_GetString(objv[2]);
    }

    const char *pwd = getInterpPWD(interp);
    const char *fileN = Tcl_GetString(fileName);

    simulationInfo.addInputFile(fileN, pwd);

#ifndef _TCL85
    return Tcl_EvalFile(interp, fileN);
#else
    return Tcl_FSEvalFileEx(interp, fileName, encodingName);
#endif
}

#ifdef _RELIABILITY   

// -- optimization Quan March 2010  (5)
int 
optimization(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  
  if (theOptimizationBuilder == 0) {

  theOptimizationBuilder = new TclOptimizationBuilder(theDomain , interp);
      

    return TCL_OK;
  }
  else  
    return TCL_ERROR;
}


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
sensitivityIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // Does nothing, but keeping command for backward compatibility
  return TCL_OK;
}

int 
sensitivityAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	bool withRespectToRVs = true;
	bool newalgorithm = false;
	int analysisTypeTag = 1;
	if(theStaticIntegrator != 0) {
	    theSensitivityIntegrator = theStaticIntegrator;
	} else if(theTransientIntegrator != 0) {
	    theSensitivityIntegrator = theTransientIntegrator;
	}
	// 1: compute at each step (default); 2: compute by command

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

	if (strcmp(argv[1],"-computeAtEachStep") == 0)
	  analysisTypeTag = 1;
	else if (strcmp(argv[1],"-computeByCommand") == 0)
	  analysisTypeTag = 2;
	else {
	  opserr << "Unknown sensitivity algorithm option: " << argv[1] << endln;
	  return TCL_ERROR;
	}
	
	ReliabilityDomain *theReliabilityDomain;
	theReliabilityDomain = theReliabilityBuilder->getReliabilityDomain();
	if(newalgorithm){
	  // theSensitivityAlgorithm = new 
	  //   NewSensitivityAlgorithm(theReliabilityDomain, 
	  // 			    &theDomain,
	  // 			    theAlgorithm,
	  // 			    theSensitivityIntegrator,
	  // 			    analysisTypeTag);
	} else {
	//  theSensitivityAlgorithm = new 
	 //   SensitivityAlgorithm(&theDomain,
	//			 theAlgorithm,
	//			 theSensitivityIntegrator,
	//			 analysisTypeTag);
	
	
    IncrementalIntegrator *theIntegrator = 0;
      
	   if (theStaticAnalysis != 0 && theStaticIntegrator != 0) {
               theIntegrator = theStaticIntegrator;
             
        theIntegrator->setComputeType(analysisTypeTag);
       	theIntegrator->activateSensitivityKey();
    
	   } else if (theTransientAnalysis != 0 && theTransientIntegrator != 0) {
  
    theIntegrator = theTransientIntegrator;
	theIntegrator->setComputeType(analysisTypeTag);
    theIntegrator->activateSensitivityKey();

	}
		    
	if (theIntegrator == 0) {
	  opserr << "ERROR: Could not create theSensitivityAlgorithm. " << endln;
	  return TCL_ERROR;
	}	
	// ---- by Quan 2009 for recover the previous framework ---

	if (theIntegrator->shouldComputeAtEachStep()) {
	
	    //if (theStaticAnalysis !=0)
		    //theStaticAnalysis->setSensitivityAlgorithm(theIntegrator);
	    //else if (theTransientAnalysis !=0)
		    //theTransientAnalysis->setSensitivityAlgorithm(theIntegrator);
	    //else if (theVariableTimeStepTransientAnalysis !=0)
		    //theVariableTimeStepTransientAnalysis->setSensitivityAlgorithm(theIntegrator);
		// else {
		// 	// do nothing		
		// }
	
	
	}
        }

	return TCL_OK;
}

// AddingSensitivity:END /////////////////////////////////////////////////

#endif


int 
wipeModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  wipeAnalysis(clientData, interp, argc, argv);

  /*
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
  */

  // NOTE : DON'T do the above on theVariableTimeStepAnalysis
  // as it and theTansientAnalysis are one in the same
  if (theDatabase != 0)
    delete theDatabase;

  theDomain.clearAll();
  OPS_clearAllUniaxialMaterial();
  OPS_clearAllNDMaterial();
  OPS_clearAllSectionForceDeformation();

  OPS_clearAllHystereticBackbone();
  OPS_clearAllStiffnessDegradation();
  OPS_clearAllStrengthDegradation();
  OPS_clearAllUnloadingRule();

  ops_Dt = 0.0;


#ifdef _PARALLEL_PROCESSING
  OPS_PARTITIONED = false;
#endif

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
  //theSensitivityAlgorithm =0;
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
  theEigenSOE =0;
  theStaticIntegrator =0;
  theTransientIntegrator =0;
  theStaticAnalysis =0;
  theTransientAnalysis =0;    
  theVariableTimeStepTransientAnalysis =0;   
  //  theSensitivityAlgorithm=0; 
  thePFEMAnalysis = 0;
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

// by SAJalali
int OPS_recorderValue(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	// make sure at least one other argument to contain type of system

	// clmnID starts from 1
	if (argc < 3) {
		opserr << "WARNING want - recorderValue recorderTag clmnID <rowOffset> <-reset>\n";
		return TCL_ERROR;
	}

	int tag, rowOffset;
	int dof = -1;

	if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
		opserr << "WARNING recorderValue recorderTag? clmnID <rowOffset> <-reset> could not read recorderTag \n";
		return TCL_ERROR;
	}

	if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
		opserr << "WARNING recorderValue recorderTag? clmnID - could not read clmnID \n";
		return TCL_ERROR;
	}
	dof--;
	rowOffset = 0;
	int curArg = 3;
	if (argc > curArg)
	{
		if (Tcl_GetInt(interp, argv[curArg], &rowOffset) != TCL_OK) {
			opserr << "WARNING recorderValue recorderTag? clmnID <rowOffset> <-reset> could not read rowOffset \n";
			return TCL_ERROR;
		}
		curArg++;
	}
	bool reset = false;
	if (argc > curArg)
	{
		if (strcmp(argv[curArg], "-reset") == 0)
			reset = true;
		curArg++;
	}
	Recorder* theRecorder = theDomain.getRecorder(tag);
	double res = theRecorder->getRecordedValue(dof, rowOffset, reset);
	// now we copy the value to the tcl string that is returned
	//sprintf(interp->result, "%35.8f ", res);
	char buffer [40];
	sprintf(buffer,"%35.8f", res);	
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);
	
	return TCL_OK;
}

int 
resetModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	theDomain.revertToStart();

	if (theTransientIntegrator != 0) {
		theTransientIntegrator->revertToStart();
	}
	
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
setCreep(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
      opserr << "WARNING illegal command - setCreep value? \n";
      return TCL_ERROR;
  }
  int newFlag;
  if (Tcl_GetInt(interp, argv[1], &newFlag) != TCL_OK) {
      opserr << "WARNING reading creep value - setCreep newFlag? \n";
      return TCL_ERROR;
  } else {
      theDomain.setCreep(newFlag);
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
  
  // get the display format
  char format[80];
  if (argc == 1) {
    //      strcpy(format,"%f");
    sprintf(format,"%f",time);
  } else if (argc == 2) {
    //      strcpy(format,argv[1]);
    sprintf(format,argv[1],time);    
  }
  
  // now we copy the value to the tcl string that is returned
  //  sprintf(interp->result,format,time);
  Tcl_SetResult(interp, format, TCL_VOLATILE);
  return TCL_OK;
}

int 
getLoadFactor(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "WARNING no load pattern supplied -- getLoadFactor\n";
    return TCL_ERROR;
  }

  int pattern;
  if (Tcl_GetInt(interp, argv[1], &pattern) != TCL_OK) {
    opserr << "ERROR reading load pattern tag -- getLoadFactor\n";
    return TCL_ERROR;
  }

  LoadPattern *thePattern = theDomain.getLoadPattern(pattern);
  if (thePattern == 0) {
    opserr << "ERROR load pattern with tag " << pattern << " not found in domain -- getLoadFactor\n";
    return TCL_ERROR;
  }

  double factor = thePattern->getLoadFactor();

  //  sprintf(interp->result,"%f",factor);

  char buffer [40];
  sprintf(buffer,"%35.20f", factor);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

////////////////////////////////////////////////Abbas//////////////////////////////

int 
sensLambda(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv )
{
  if (argc < 3) {
    opserr << "WARNING no load pattern supplied -- getLoadFactor\n";
    return TCL_ERROR;
  }
  
  int pattern, paramTag;
  if (Tcl_GetInt(interp, argv[1], &pattern) != TCL_OK) {
    opserr << "ERROR reading load pattern tag -- getLoadFactor\n";
    return TCL_ERROR;
  }
  
  LoadPattern *thePattern = theDomain.getLoadPattern(pattern);
  if (thePattern == 0) {
    opserr << "ERROR load pattern with tag " << pattern << " not found in domain -- getLoadFactor\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &paramTag) != TCL_OK) {
    opserr << "WARNING sensLambda patternTag?  paramTag?- could not read paramTag? ";
    return TCL_ERROR;	        
  }        
  Parameter *theParam = theDomain.getParameter(paramTag);
  if (theParam == 0) {
    opserr << "sensLambda: parameter " << paramTag << " not found" << endln;
    return TCL_ERROR;
  }
  
  IncrementalIntegrator *theIntegrator = 0;
  
  if (theStaticAnalysis != 0 && theStaticIntegrator != 0) {
    theIntegrator = theStaticIntegrator;
    //   opserr<<" commands.cpp: calling static integrator"<<endln;
    
  } else if (theTransientAnalysis != 0 && theTransientIntegrator != 0) {
    
    theIntegrator = theTransientIntegrator;
    //   opserr<<"commands.cpp:calling transientIntegrator"<<endln;
  } 
  
  int gradIndex = theParam->getGradIndex();
  // double   factor = thePattern->getSensLambda(theIntegrator);
  //double factor = theIntegrator->dLambdadh(gradIndex);
  double factor = thePattern->getLoadFactorSensitivity(gradIndex);

  char buffer[40];
  sprintf(buffer,"%35.20f",factor);
  
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}



//////////////////////////////////////////////Abbas///////////////////////////////////////


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
partitionModel(int eleTag)
{
  if (OPS_PARTITIONED == true)
    return 0;

  int result = 0;
  
  if (OPS_theChannels != 0)
    delete [] OPS_theChannels;

  OPS_theChannels = new Channel *[OPS_NUM_SUBDOMAINS];

  // create some subdomains
  for (int i=1; i<=OPS_NUM_SUBDOMAINS; i++) {
    if (i != OPS_MAIN_DOMAIN_PARTITION_ID) {
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

 // opserr << "commands.cpp - partition numPartitions: " << OPS_NUM_SUBDOMAINS << endln;

  result = theDomain.partition(OPS_NUM_SUBDOMAINS, OPS_USING_MAIN_DOMAIN, OPS_MAIN_DOMAIN_PARTITION_ID, eleTag);

  if (result < 0) 
    return result;

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


int 
opsPartition(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
#ifdef _PARALLEL_PROCESSING
  int eleTag;
  if (argc == 2) {
    if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
      ;
    }
  }
  partitionModel(eleTag);
#endif
  return TCL_OK;
}

//
// command invoked to build the model, i.e. to invoke analyze() 
// on the Analysis object
//
int 
analyzeModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int result = 0;

#ifdef _PARALLEL_PROCESSING
  if (OPS_PARTITIONED == false && OPS_NUM_SUBDOMAINS > 1) {
    if (partitionModel(0) < 0) {
      opserr << "WARNING before analysis; partition failed - too few elements\n";
      OpenSeesExit(clientData, interp, argc, argv);
      return TCL_ERROR;
    }
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

  } else if(thePFEMAnalysis != 0) {
      result = thePFEMAnalysis->analyze();

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

  char buffer [10];
  sprintf(buffer,"%d", result);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);
  
  //  sprintf(interp->result,"%d",result);    

  return TCL_OK;

}

int 
printElement(ClientData clientData, Tcl_Interp *interp, int argc, 
	     TCL_Char **argv, OPS_Stream &output);


int 
printNode(ClientData clientData, Tcl_Interp *interp, int argc,  
	  TCL_Char **argv, OPS_Stream &output);
	  
int 
printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, 
		TCL_Char **argv, OPS_Stream &output);	  
		
int 
printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, 
	       TCL_Char **argv, OPS_Stream &output);	  		


int 
printModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int currentArg = 1;
  int res = 0;

  int flag = OPS_PRINT_CURRENTSTATE;
  
  FileStream outputFile;
  OPS_Stream *output = &opserr;
  bool done = false;

  // if just 'print' then print out the entire domain
  if (argc == currentArg) {
    opserr << theDomain;
    return TCL_OK;
  }    

  while(done == false) {
    // if 'print ele i j k..' print out some elements
    if ((strcmp(argv[currentArg],"-ele") == 0) || (strcmp(argv[currentArg],"ele") == 0)) {
      currentArg++;
      res = printElement(clientData, interp, argc-currentArg, argv+currentArg, *output);    
      done = true;
    }
    // if 'print node i j k ..' print out some nodes
    else if ((strcmp(argv[currentArg],"-node") == 0) || (strcmp(argv[currentArg],"node") == 0)) {
      currentArg++;      
      res = printNode(clientData, interp, argc-currentArg, argv+currentArg, *output);
      done = true;
    }
  
    // if 'print integrator flag' print out the integrator
    else if ((strcmp(argv[currentArg],"integrator") == 0) || 
	     (strcmp(argv[currentArg],"-integrator") == 0)) {
      currentArg++;
      res = printIntegrator(clientData, interp, argc-currentArg, argv+currentArg, *output);  
      done = true;
    }

    // if 'print algorithm flag' print out the algorithm
    else if ((strcmp(argv[currentArg],"algorithm") == 0) || 
	     (strcmp(argv[currentArg],"-algorithm") == 0)) {
      currentArg++;
      res = printAlgorithm(clientData, interp, argc-currentArg, argv+currentArg, *output);    
      done = true;
    }

    else if ((strcmp(argv[currentArg],"-JSON") == 0)) {
      currentArg++;
      flag = OPS_PRINT_PRINTMODEL_JSON;
    }

    else {

      if ((strcmp(argv[currentArg],"file") == 0) || 
	  (strcmp(argv[currentArg],"-file") == 0)) 
	currentArg++;
	
      openMode mode = APPEND;
      if (flag == OPS_PRINT_PRINTMODEL_JSON)
          mode = OVERWRITE;
      if (outputFile.setFile(argv[currentArg], mode) != 0) {
          opserr << "print <filename> .. - failed to open file: " << argv[currentArg] << endln;
          return TCL_ERROR;
      }
      currentArg++;

      // if just 'print <filename>' then print out the entire domain to eof
      if (argc == currentArg) {
          if (flag == OPS_PRINT_PRINTMODEL_JSON)
              simulationInfo.Print(outputFile, flag);
          
          theDomain.Print(outputFile, flag);
          return TCL_OK;
      }

      output = &outputFile;

    }
  }

  // close the output file
  outputFile.close();
  return res;
}

								   

// printNode():
// function to print out the nodal information conatined in line
//     print <filename> node <flag int> <int int int>
// input: nodeArg: integer equal to arg count to node plus 1
//        output: output stream to which the results are sent
// 
int 
printNode(ClientData clientData, Tcl_Interp *interp, int argc, 
	  TCL_Char **argv, OPS_Stream &output)
{
  int flag = 0; // default flag sent to a nodes Print() method
  int nodeArg = 0;

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc == 0) { 
    NodeIter &theNodes = theDomain.getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0)
      theNode->Print(output);
    return TCL_OK;
  }    

  // if 'print <filename> node flag int <int int ..>' get the flag
  if ((strcmp(argv[0],"flag") == 0) ||
      (strcmp(argv[0],"-flag") == 0)) { 
      // get the specified flag
    if (argc <= nodeArg) {
      opserr << "WARNING print <filename> node <flag int> no int specified \n";
      return TCL_ERROR;
    }    
    if (Tcl_GetInt(interp, argv[1], &flag) != TCL_OK) {
      opserr << "WARNING print node failed to get integer flag: \n";
      opserr << argv[nodeArg] << endln; 
      return TCL_ERROR;
    }    
    nodeArg += 2;
  }

  // now print the nodes with the specified flag, 0 by default

  // if 'print <filename> node flag' 
  //     print out all the nodes in the domain with flag
  if (nodeArg == argc) { 
    NodeIter &theNodes = theDomain.getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0)
      theNode->Print(output, flag);
    return TCL_OK;
  } else { 
    // otherwise print out the specified nodes i j k .. with flag
    int numNodes = argc-nodeArg;
    ID *theNodes = new ID(numNodes);
    for (int i= 0; i<numNodes; i++) {
      int nodeTag;
      if (Tcl_GetInt(interp, argv[nodeArg], &nodeTag) != TCL_OK) {
	opserr << "WARNING print node failed to get integer: " << argv[nodeArg] << endln;
	return TCL_ERROR;
      }
      (*theNodes)(i) = nodeTag;
      nodeArg++;
    }

    theDomain.Print(output, theNodes, 0, flag);
    delete theNodes;
  }    

  return TCL_OK;

}


int 
printElement(ClientData clientData, Tcl_Interp *interp, int argc, 
	  TCL_Char **argv, OPS_Stream &output)
{
  int flag = 0; // default flag sent to a nodes Print() method
  int eleArg = 0;

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc == 0) { 
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)
      theElement->Print(output);
    return TCL_OK;
  }    

  // if 'print <filename> Element flag int <int int ..>' get the flag
  if ((strcmp(argv[0],"flag") == 0) ||
      (strcmp(argv[0],"-flag")) == 0) { // get the specified flag
    if (argc < 2) {
      opserr << "WARNING print <filename> ele <flag int> no int specified \n";
      return TCL_ERROR;
    }    
    if (Tcl_GetInt(interp, argv[1], &flag) != TCL_OK) {
      opserr << "WARNING print ele failed to get integer flag: \n";
      opserr << argv[eleArg] << endln; 
      return TCL_ERROR;
    }    
    eleArg += 2;
  }

  // now print the Elements with the specified flag, 0 by default
  if (argc == eleArg) { 
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)      
      theElement->Print(output, flag);
    return TCL_OK;
  } else { 

    // otherwise print out the specified nodes i j k .. with flag
    int numEle = argc-eleArg;
    ID *theEle = new ID(numEle);
    for (int i= 0; i<numEle; i++) {
      int eleTag;
      if (Tcl_GetInt(interp, argv[i+eleArg], &eleTag) != TCL_OK) {
	opserr << "WARNING print ele failed to get integer: " << argv[i] << endln;
	return TCL_ERROR;
      }
      (*theEle)(i) = eleTag;
    }

    theDomain.Print(output, 0, theEle, flag);
    delete theEle;
  }

  return TCL_OK;
}


int 
printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, 
	       TCL_Char **argv, OPS_Stream &output)
{
  int eleArg = 0;
  if (theAlgorithm == 0)
      return TCL_OK;

  // if just 'print <filename> algorithm'- no flag
  if (argc == 0) { 
      theAlgorithm->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> Algorithm flag' get the flag
  int flag;  
  if (Tcl_GetInt(interp, argv[eleArg], &flag) != TCL_OK) {  
      opserr << "WARNING print algorithm failed to get integer flag: \n";
      opserr << argv[eleArg] << endln; 
      return TCL_ERROR;
  }    
  theAlgorithm->Print(output,flag);
  return TCL_OK;  
}


int 
printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, 
		TCL_Char **argv, OPS_Stream &output)
{
  int eleArg = 0;
  if (theStaticIntegrator == 0 && theTransientIntegrator == 0)
      return TCL_OK;
  
  IncrementalIntegrator *theIntegrator;
  if (theStaticIntegrator != 0)
      theIntegrator = theStaticIntegrator;
  else
      theIntegrator = theTransientIntegrator;

  // if just 'print <filename> algorithm'- no flag
  if (argc == 0) { 
      theIntegrator->Print(output);
      return TCL_OK;
  }    

  // if 'print <filename> Algorithm flag' get the flag
  int flag;  
  if (Tcl_GetInt(interp, argv[eleArg], &flag) != TCL_OK) {  
      opserr << "WARNING print algorithm failed to get integer flag: \n";
      opserr << argv[eleArg] << endln; 
      return TCL_ERROR;
  }    
  theIntegrator->Print(output,flag);
  return TCL_OK;  
}


int 
printA(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int res = 0;

  FileStream outputFile;
  OPS_Stream *output = &opserr;

  int currentArg = 1;

  if (argc > 2) {
    if ((strcmp(argv[currentArg],"file") == 0) || 
	(strcmp(argv[currentArg],"-file") == 0)) {
      currentArg++;
      
      if (outputFile.setFile(argv[currentArg]) != 0) {
	opserr << "print <filename> .. - failed to open file: " << argv[currentArg] << endln;
	return TCL_ERROR;
      }
      output = &outputFile;
    }
  }
  if (theSOE != 0) {
    if (theStaticIntegrator != 0)
      theStaticIntegrator->formTangent();
    else if (theTransientIntegrator != 0)
      theTransientIntegrator->formTangent(0);
      
    const Matrix *A = theSOE->getA();
    if (A != 0) {
      *output << *A;
    }
  }
  
  // close the output file
  outputFile.close();
  
  return res;
}

int 
printB(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int res = 0;

  FileStream outputFile;
  OPS_Stream *output = &opserr;
  //  bool done = false;

  int currentArg = 1;

  if (argc > 2) {
    if ((strcmp(argv[currentArg],"file") == 0) || 
	(strcmp(argv[currentArg],"-file") == 0)) {
      currentArg++;
      
      if (outputFile.setFile(argv[currentArg]) != 0) {
	opserr << "print <filename> .. - failed to open file: " << argv[currentArg] << endln;
	return TCL_ERROR;
      }
      output = &outputFile;
    }
  }
  if (theSOE != 0) {
    if (theStaticIntegrator != 0)
      theStaticIntegrator->formUnbalance();
    else if (theTransientIntegrator != 0)
      theTransientIntegrator->formUnbalance();
      
    const Vector &b = theSOE->getB();
    *output << b;
  }
  
  // close the output file
  outputFile.close();
  
  return res;
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

    //
    // do nothing if request is for the same analysis type!
    //

    if ((strcmp(argv[1],"Static") == 0) && (theStaticAnalysis != 0))
      return TCL_OK;

    if (((strcmp(argv[1],"VariableTimeStepTransient") == 0) ||
	(strcmp(argv[1],"TransientWithVariableTimeStep") == 0) ||
	 (strcmp(argv[1],"VariableTransient") == 0)) && 
	(theVariableTimeStepTransientAnalysis != 0))
      return TCL_OK;

    if ((strcmp(argv[1],"Transient") == 0) && (theTransientAnalysis != 0))
      return TCL_OK;

    //
    // analysis changing .. delete the old analysis
    //

    if (theStaticAnalysis != 0) {
	delete theStaticAnalysis;
	theStaticAnalysis = 0;
	opserr << "WARNING: analysis .. StaticAnalysis already exists => wipeAnalysis not invoked, problems may arise\n";
    }

    if (theTransientAnalysis != 0) {
	delete theTransientAnalysis;
	theTransientAnalysis = 0;
	theVariableTimeStepTransientAnalysis = 0;
	opserr << "WARNING: analysis .. TransientAnalysis already exists => wipeAnalysis not invoked, problems may arise\n";
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
#ifdef _PARALLEL_PROCESSING
	    theSOE = new DistributedProfileSPDLinSOE(*theSolver);
#else
	    theSOE = new ProfileSPDLinSOE(*theSolver);      
#endif
	}
    
	theStaticAnalysis = new StaticAnalysis(theDomain,
					       *theHandler,
					       *theNumberer,
					       *theAnalysisModel,
					       *theAlgorithm,
					       *theSOE,
					       *theStaticIntegrator,
					       theTest);

	#ifdef _PARALLEL_INTERPRETERS
	if (setMPIDSOEFlag) {
	  ((MPIDiagonalSOE*) theSOE)->setAnalysisModel(*theAnalysisModel);
	}
#endif

// AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
	if (theSensitivityAlgorithm != 0 && theSensitivityAlgorithm->shouldComputeAtEachStep()) {
	    //theStaticAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
	}
#endif
// AddingSensitivity:END /////////////////////////////////
    } else if(strcmp(argv[1], "PFEM") == 0) {

        if(argc < 5) {
            opserr<<"WARNING: wrong no of args -- analysis PFEM dtmax dtmin gravity <ratio>\n";
            return TCL_ERROR;
        }
        double dtmax, dtmin, gravity, ratio=0.5;
        if(Tcl_GetDouble(interp, argv[2], &dtmax) != TCL_OK) {  
            opserr<<"WARNING: invalid dtmax "<<argv[2]<<"\n";
            return TCL_ERROR;
        }
        if(Tcl_GetDouble(interp, argv[3], &dtmin) != TCL_OK) {  
            opserr<<"WARNING: invalid dtmin "<<argv[3]<<"\n";
            return TCL_ERROR;
        }
	if(Tcl_GetDouble(interp, argv[4], &gravity) != TCL_OK) {  
            opserr<<"WARNING: invalid gravity "<<argv[4]<<"\n";
            return TCL_ERROR;
        }
        if(argc > 5) {
            if(Tcl_GetDouble(interp, argv[5], &ratio) != TCL_OK) {  
                opserr<<"WARNING: invalid ratio "<<argv[5]<<"\n";
                return TCL_ERROR;
            }
        }

        if(theAnalysisModel == 0) {
            theAnalysisModel = new AnalysisModel();
        }
        if(theTest == 0) {
            //theTest = new CTestNormUnbalance(1e-2,10000,1,2,3);
            theTest = new CTestPFEM(1e-2,1e-2,1e-2,1e-2,1e-4,1e-3,10000,100,1,2);
        }
        if(theAlgorithm == 0) {
            theAlgorithm = new NewtonRaphson(*theTest);
        }
        if(theHandler == 0) {
            theHandler = new TransformationConstraintHandler();
        }
        if(theNumberer == 0) {
            RCM* theRCM = new RCM(false);
            theNumberer = new DOF_Numberer(*theRCM);
        }
        if(theTransientIntegrator == 0) {
            theTransientIntegrator = new PFEMIntegrator();
        }
        if(theSOE == 0) {
            PFEMSolver* theSolver = new PFEMSolver();
            theSOE = new PFEMLinSOE(*theSolver);
        }
        thePFEMAnalysis = new PFEMAnalysis(theDomain,
                                           *theHandler,
                                           *theNumberer,
                                           *theAnalysisModel,
                                           *theAlgorithm,
                                           *theSOE,
                                           *theTransientIntegrator,
                                           theTest,dtmax,dtmin,gravity,ratio);

        theTransientAnalysis = thePFEMAnalysis;

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
#ifdef _PARALLEL_PROCESSING
	    theSOE = new DistributedProfileSPDLinSOE(*theSolver);
#else
	    theSOE = new ProfileSPDLinSOE(*theSolver);      
#endif
	}
	
	int count = 2;
	int numSubLevels = 0;
	int numSubSteps = 10;
	while (count < argc) {
	  if (strcmp(argv[count],"-numSubLevels") == 0) {
	    count++;
	    if (count < argc)
	      if (Tcl_GetInt(interp, argv[count], &numSubLevels) != TCL_OK)
		return TCL_ERROR;		     
	  }
	  else if ((strcmp(argv[count],"-numSubSteps") == 0) ) {
	    count++;
	    if (count < argc)
	      if (Tcl_GetInt(interp, argv[count], &numSubSteps) != TCL_OK)
		return TCL_ERROR;		     
	  }
	  count++;
	}

	theTransientAnalysis = new DirectIntegrationAnalysis(theDomain,
							     *theHandler,
							     *theNumberer,
							     *theAnalysisModel,
							     *theAlgorithm,
							     *theSOE,
							     *theTransientIntegrator,
							     theTest,
							     numSubLevels,
							     numSubSteps);
	  ;
#ifdef _PARALLEL_INTERPRETERS
	if (setMPIDSOEFlag) {
	  ((MPIDiagonalSOE*) theSOE)->setAnalysisModel(*theAnalysisModel);
	}
#endif

// AddingSensitivity:BEGIN ///////////////////////////////
#ifdef _RELIABILITY
	if (theSensitivityAlgorithm != 0 && theSensitivityAlgorithm->shouldComputeAtEachStep()) {

	  /* This if-statement cannot possibly stay in the code -- MHS
	  if(theSensitivityAlgorithm->newAlgorithm()){
	    opserr << "WARNING original sensitivity algorothm needs to be specified \n";
	    opserr << "for static analysis \n";
	    return TCL_ERROR;
	  }
	  */
		
	  //theTransientAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
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
#ifdef _PARALLEL_PROCESSING
	    theSOE = new DistributedProfileSPDLinSOE(*theSolver);
#else
	    theSOE = new ProfileSPDLinSOE(*theSolver);      
#endif
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

	#ifdef _RELIABILITY

	//////////////////////////////////
    ////// added by K Fujimura ///////
	
	//FMK RELIABILITY
	/*
    } else if (strcmp(argv[1],"ReliabilityStatic") == 0) {
		// make sure all the components have been built,
		// otherwise print a warning and use some defaults
		if (theAnalysisModel == 0) 
			theAnalysisModel = new AnalysisModel();
		if (theTest == 0) 
		  theTest = new CTestNormUnbalance(1.0e-6,25,0);       
		if (theAlgorithm == 0) {
		    opserr << "WARNING analysis Static - no Algorithm yet specified, \n";
			opserr << " NewtonRaphson default will be used\n";	    
		    theAlgorithm = new NewtonRaphson(*theTest); 	}
		if (theHandler == 0) {
			opserr << "WARNING analysis Static - no ConstraintHandler yet specified, \n";
			opserr << " PlainHandler default will be used\n";
			theHandler = new PlainHandler();       	}
		if (theNumberer == 0) {
		    opserr << "WARNING analysis Static - no Numberer specified, \n";
			opserr << " RCM default will be used\n";
			RCM *theRCM = new RCM(false);	
			theNumberer = new DOF_Numberer(*theRCM);    		}
		if (theStaticIntegrator == 0) {
			opserr << "Fatal ! theStaticIntegrator must be defined before defining\n";
			opserr << "ReliabilityStaticAnalysis by NewStaticSensitivity\n";
			return TCL_ERROR;
		}
		if (theSOE == 0) {
			opserr << "WARNING analysis Static - no LinearSOE specified, \n";
			opserr << " ProfileSPDLinSOE default will be used\n";
			ProfileSPDLinSolver *theSolver;
			theSolver = new ProfileSPDLinDirectSolver(); 	
			theSOE = new ProfileSPDLinSOE(*theSolver);      	}
    
		theReliabilityStaticAnalysis = new ReliabilityStaticAnalysis(theDomain,
					       *theHandler,
					       *theNumberer,
					       *theAnalysisModel,
					       *theAlgorithm,
					       *theSOE,
					       *theStaticIntegrator,
					       theTest);

		if (theSensitivityAlgorithm != 0 && theSensitivityAlgorithm->shouldComputeAtEachStep()) {

		  //This if-statement cannot stay -- MHS
		  //if(!theSensitivityAlgorithm->newAlgorithm()){
		  //  opserr << "WARNING new sensitivity algorothm needs to be specified \n";
		   // opserr << "for reliability static analysis \n";
		   // return TCL_ERROR;
		  //}
		  

		  //theStaticAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
		} else {
			opserr << "Faltal SensitivityAlgorithm must be definde before defining \n";
			opserr << "ReliabilityStaticAnalysis with computeateachstep\n";
			return TCL_ERROR;
		}

    } else if (strcmp(argv[1],"ReliabilityTransient") == 0) {
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
			opserr << "Fatal ! theTransientIntegrator must be defined before defining\n";
			opserr << "ReliabilityTransientAnalysis by NewNewmarkWithSensitivity\n";
			return TCL_ERROR;
		}
		if (theSOE == 0) {
			opserr << "WARNING analysis Transient dt tFinal - no LinearSOE specified, \n";
			opserr << " ProfileSPDLinSOE default will be used\n";
			ProfileSPDLinSolver *theSolver;
			theSolver = new ProfileSPDLinDirectSolver(); 	
			theSOE = new ProfileSPDLinSOE(*theSolver);      
		}
    
		theReliabilityTransientAnalysis = new ReliabilityDirectIntegrationAnalysis(theDomain,
							     *theHandler,
							     *theNumberer,
							     *theAnalysisModel,
							     *theAlgorithm,
							     *theSOE,
							     *theTransientIntegrator,
							     theTest);

		if (theSensitivityAlgorithm != 0 && theSensitivityAlgorithm->shouldComputeAtEachStep()) {

		  //This if-statement must go -- MHS
		  //if(!theSensitivityAlgorithm->newAlgorithm()){
		   // opserr << "WARNING new sensitivity algorothm needs to be specified \n";
		   // opserr << "for reliability static analysis \n";
		   // return TCL_ERROR;
		  //}
		  

			theReliabilityTransientAnalysis->setSensitivityAlgorithm(theSensitivityAlgorithm);
		}else{
			opserr << "Faltal SensitivityAlgorithm must be definde before defining \n";
			opserr << "ReliabilityStaticAnalysis with computeateachstep\n";
			return TCL_ERROR;
		}
		FMK RELIABILITY
	    *************************/
// AddingSensitivity:END /////////////////////////////////
#endif

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

    if (theEigenSOE != 0) {
      if (theStaticAnalysis != 0) {
	theStaticAnalysis->setEigenSOE(*theEigenSOE);
      } else if (theTransientAnalysis != 0) {
	theTransientAnalysis->setEigenSOE(*theEigenSOE);
      }
    }
	

    return TCL_OK;
}


typedef struct externalClassFunction {
  char *funcName;
  void *(*funcPtr)();
  struct externalClassFunction *next;
} ExternalClassFunction;

static ExternalClassFunction *theExternalSolverCommands = NULL;
static ExternalClassFunction *theExternalStaticIntegratorCommands = NULL;
static ExternalClassFunction *theExternalTransientIntegratorCommands = NULL;
static ExternalClassFunction *theExternalAlgorithmCommands = NULL;

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

#ifdef _CUDA
  else if ((strcmp(argv[1],"BandGeneral_Single") == 0) || (strcmp(argv[1],"BandGEN_Single") == 0)
      || (strcmp(argv[1],"BandGen_Single") == 0)){
    BandGenLinLapackSolver_Single    *theSolver = new BandGenLinLapackSolver_Single();
    theSOE = new BandGenLinSOE_Single(*theSolver);      
  }
#endif

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
  // Diagonal SOE & SOLVER
  else if (strcmp(argv[1],"MPIDiagonal") == 0) {
#ifdef _PARALLEL_INTERPRETERS
      MPIDiagonalSolver    *theSolver = new MPIDiagonalSolver();   
      theSOE = new MPIDiagonalSOE(*theSolver);
      setMPIDSOEFlag = true;
#else
      DiagonalSolver    *theSolver = new DiagonalDirectSolver();   
      theSOE = new DiagonalSOE(*theSolver);
#endif
  } 


  // PROFILE SPD SOE * SOLVER
  else if (strcmp(argv[1],"SProfileSPD") == 0) {
    // now must determine the type of solver to create from rest of args
    SProfileSPDLinSolver *theSolver = new SProfileSPDLinSolver(); 	
    theSOE = new SProfileSPDLinSOE(*theSolver);      
  }

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

#ifdef _PARALLEL_INTERPRETERS
  else if (strcmp(argv[1],"ParallelProfileSPD") == 0) {
    ProfileSPDLinSolver *theSolver = new ProfileSPDLinDirectSolver(); 	
    DistributedProfileSPDLinSOE *theParallelSOE = new DistributedProfileSPDLinSOE(*theSolver);
    theSOE = theParallelSOE;
    theParallelSOE->setProcessID(OPS_rank);
    theParallelSOE->setChannels(numChannels, theChannels);
  }
#endif

  else if(strcmp(argv[1], "PFEM") == 0) {
      if(argc <= 2) {
          PFEMSolver* theSolver = new PFEMSolver();
          theSOE = new PFEMLinSOE(*theSolver);
      } else if(strcmp(argv[2], "-quasi") == 0) {
          PFEMCompressibleSolver* theSolver = new PFEMCompressibleSolver();
          theSOE = new PFEMCompressibleLinSOE(*theSolver);
      } else if (strcmp(argv[2],"-mumps") ==0) {
#ifdef _PARALLEL_INTERPRETERS
	  int relax = 20;
	  if (argc > 3) {
	      if (Tcl_GetInt(interp, argv[3], &relax) != TCL_OK) {
		  opserr<<"WARNING: failed to read relax\n";
		  return TCL_ERROR;
	      }
	  }
	  PFEMSolver_Mumps* theSolver = new PFEMSolver_Mumps(relax,0,0,0);
          theSOE = new PFEMLinSOE(*theSolver);
#endif
      } else if (strcmp(argv[2],"-quasi-mumps")==0) {
#ifdef _PARALLEL_INTERPRETERS
	  int relax = 20;
	  if (argc > 3) {
	      if (Tcl_GetInt(interp, argv[3], &relax) != TCL_OK) {
		  opserr<<"WARNING: failed to read relax\n";
		  return TCL_ERROR;
	      }
	  }
	  PFEMCompressibleSolver_Mumps* theSolver = new PFEMCompressibleSolver_Mumps(relax,0,0);
          theSOE = new PFEMCompressibleLinSOE(*theSolver);
#endif
      }
  }

#ifdef _CUSP
  else if ((_stricmp(argv[1],"CuSP")==0)) {

    
    double relTol=1e-6;
    int maxInteration=100000;
    int preCond=1;			//diagonal
    int solver=0;				//Bicg
    int count = 2;
    
    
    while (count < argc) {
      
      if (_stricmp(argv[count],"-rTol") == 0) {
	count++;
	if (count < argc)
	  if (Tcl_GetDouble(interp, argv[count], &relTol) != TCL_OK)
	    return TCL_ERROR;		     
      }
      else if ((_stricmp(argv[count],"-mInt") == 0) ) {
	count++;
	if (count < argc)
	  if (Tcl_GetInt(interp, argv[count], &maxInteration) != TCL_OK)
	    return TCL_ERROR;		     
      }
      else if ((_stricmp(argv[count],"-pre") == 0) ) {
	count++;
	if (count < argc)
	  if ((_stricmp(argv[count],"none") == 0))
	    preCond=0;
	  else if ((_stricmp(argv[count],"diagonal") == 0))
	    preCond=1;
	  else if ((_stricmp(argv[count],"ainv") == 0))
	    preCond=2;
	  else
	    return TCL_ERROR;
      } else if ((_stricmp(argv[count],"-solver") == 0)) {
	count++;
	if (count < argc)
	  if ((_stricmp(argv[count],"bicg") == 0))
	    solver=0;
	  else if ((_stricmp(argv[count],"bicgstab") == 0))
	    solver=1;
	  else if ((_stricmp(argv[count],"cg") == 0))
	    solver=2;
	  else if ((_stricmp(argv[count],"gmres") == 0))
	    solver=3;
	  else
	    return TCL_ERROR;
      }
      count++;
    }
    
    CuSPSolver* theSolver = new CuSPSolver(maxInteration,relTol,preCond,solver);
    theSOE = new SparseGenRowLinSOE(* theSolver);
    
  }
#endif

#if defined(_CULAS4) || defined(_CULAS5)
  // CULA SPARSE
  else if ((strcmp(argv[1],"CulaSparse")==0))
  {
    double absTol = 1.0e-6;
    double relTol=1e-6;

    int maxInteration=100000;
    
    int preCond=5;		//fainv
#ifdef _CULAS4
    preCond = 1;
#endif
    int solver=0;		//cg
    int count = 2;
    int single=0; 
    int host=0;
    
    while (count < argc) {
      
      if (strcmp(argv[count],"-rTol") == 0) {
	count++;
	if (count < argc)
	  if (Tcl_GetDouble(interp, argv[count], &relTol) != TCL_OK)
	    return TCL_ERROR;		     
      }
      else if ((strcmp(argv[count],"-mInt") == 0) ) {
	count++;
	if (count < argc)
	  if (Tcl_GetInt(interp, argv[count], &maxInteration) != TCL_OK)
	    return TCL_ERROR;		     
      }
      else if ((strcmp(argv[count],"-pre") == 0) ) {
	count++;
	if (count < argc)
	  if ((strcmp(argv[count],"none") == 0))
	    preCond=0;
	  else if ((strcmp(argv[count],"jacobi") == 0))
	    preCond=1;
	  else if ((strcmp(argv[count],"blockjacobi") == 0))
	    preCond=2;
	  else if ((strcmp(argv[count],"ilu0") == 0))
	    preCond=3;
	  else if ((strcmp(argv[count],"ainv") == 0))
	    preCond=4;
	  else if ((strcmp(argv[count],"fainv") == 0))
	    preCond=5;
	  else
	    return TCL_ERROR;
      } else if ((strcmp(argv[count],"-solver") == 0)) {
	count++;
	if (count < argc)
	  if ((strcmp(argv[count],"cg") == 0))
	    solver=0;
	  else if ((strcmp(argv[count],"bicg") == 0))
	    solver=1;
	  else if ((strcmp(argv[count],"blockstab") == 0))
	    solver=2;
	  else if ((strcmp(argv[count],"blockstabl") == 0))
	    solver=3;
	  else if ((strcmp(argv[count],"gmres") == 0))
	    solver=4;
	  else if ((strcmp(argv[count],"minres") == 0))
	    solver=5;	
	  else
	    return TCL_ERROR;
      }
      else if ((strcmp(argv[count],"-single") == 0)) {
	single=1;
      }
      else if ((strcmp(argv[count],"-host") == 0)) {
	host=1;
      }
      count++;
    }
    
#ifdef _CULAS5
    CulaSparseSolverS5* theSolver = new CulaSparseSolverS5(relTol,
							   maxInteration,
							   preCond,
							   solver,
							   single,
							   host);
#else
    CulaSparseSolverS4* theSolver = new CulaSparseSolverS4(relTol,
							   maxInteration,
							   preCond,
							   solver);
#endif

    theSOE = new SparseGenRowLinSOE(*theSolver);

  }
#endif

  // SPARSE GENERAL SOE * SOLVER
  else if ((strcmp(argv[1],"SparseGeneral") == 0) || (strcmp(argv[1],"SuperLU") == 0) ||
	   (strcmp(argv[1],"SparseGEN") == 0)) {
    
    SparseGenColLinSolver *theSolver =0;    
    int count = 2;
    double thresh = 0.0;
    int npRow = 1;
    int npCol = 1;
    int np = 1;

    // defaults for threaded SuperLU

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

    int permSpec = 0;
    int panelSize = 6;
    int relax = 6;


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
      opserr << "commands.cpp: DistributedSuperLU\n";
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
    opserr << "commands.cpp: DistributedSparseGenColLinSOE\n";

    theSOE = new DistributedSparseGenColLinSOE(*theSolver);      
#else
    theSOE = new SparseGenColLinSOE(*theSolver);
#endif
  }

  
  else if ((strcmp(argv[1],"SparseSPD") == 0) || (strcmp(argv[1],"SparseSYM") == 0)) {
    // now must determine the type of solver to create from rest of args

    // now determine ordering scheme
    //   1 -- MMD
    //   2 -- ND
    //   3 -- RCM
    int lSparse = 1;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &lSparse) != TCL_OK)
	return TCL_ERROR;
    }

    SymSparseLinSolver *theSolver = new SymSparseLinSolver();
    theSOE = new SymSparseLinSOE(*theSolver, lSparse);      
  }    
  
  else if ((strcmp(argv[1],"UmfPack") == 0) || (strcmp(argv[1],"Umfpack") == 0)) {
    
    // now must determine the type of solver to create from rest of args
    int factLVALUE = 10;
    int factorOnce=0;
    int printTime = 0;
    int count = 2;

    while (count < argc) {
      if ((strcmp(argv[count],"-lValueFact") == 0) || (strcmp(argv[count],"-lvalueFact") == 0) || (strcmp(argv[count],"-LVALUE") == 0)) {
	if (Tcl_GetInt(interp, argv[count+1], &factLVALUE) != TCL_OK)
	  return TCL_ERROR;
	count++;
      } else if ((strcmp(argv[count],"-factorOnce") == 0) || (strcmp(argv[count],"-FactorOnce") ==0 )) {
	factorOnce = 1;
      } else if ((strcmp(argv[count],"-printTime") == 0) || (strcmp(argv[count],"-time") ==0 )) {
	printTime = 1;
      }
      count++;
    }
    
    UmfpackGenLinSolver *theSolver = new UmfpackGenLinSolver();
    // theSOE = new UmfpackGenLinSOE(*theSolver, factLVALUE, factorOnce, printTime);      
    theSOE = new UmfpackGenLinSOE(*theSolver);      
  }	
  
#ifdef _ITPACK
//  else if (strcmp(argv[1],"Itpack") == 0) {
//    
//    // now must determine the type of solver to create from rest of args
//    int method = 1;
//    if (argc == 3) {
//      if (Tcl_GetInt(interp, argv[2], &method) != TCL_OK)
//	return TCL_ERROR;
//    }
//    ItpackLinSolver *theSolver = new ItpackLinSolver(method);
//    theSOE = new ItpackLinSOE(*theSolver);      
//  }
#endif	 
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
	else if (strcmp(argv[count+1],"KSPCHEBYSHEV") == 0)
	  method = KSPCHEBYSHEV;
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
      // PetscSolver *theSolver = new PetscSolver(method, preconditioner, rTol, aTol, dTol, maxIts);
      PetscSolver *theSolver = new PetscSolver(method, preconditioner, rTol, aTol, dTol, maxIts);
      theSOE = new PetscSOE(*theSolver);
    } else {
      // PetscSparseSeqSolver *theSolver = new PetscSparseSeqSolver(method, preconditioner, rTol, aTol, dTol, maxIts);
      PetscSparseSeqSolver *theSolver = 0;
      theSOE = new SparseGenRowLinSOE(*theSolver);
    }
  }


#endif


#ifdef _MUMPS

  else if (strcmp(argv[1],"Mumps") == 0) {

    int icntl14 = 20;    
    int icntl7 = 7;
    int matType = 0; // 0: unsymmetric, 1: symmetric positive definite, 2: symmetric general

    int currentArg = 2;
    while (currentArg < argc) {
      if (argc > 2) {
	if (strcmp(argv[currentArg],"-ICNTL14") == 0) {
	  if (Tcl_GetInt(interp, argv[currentArg+1], &icntl14) != TCL_OK)	
	    ;
	  currentArg += 2;
	} else  if (strcmp(argv[currentArg],"-ICNTL7") == 0) {
	  if (Tcl_GetInt(interp, argv[currentArg+1], &icntl7) != TCL_OK)	
	    ;
	  currentArg += 2;
	} else  if (strcmp(argv[currentArg],"-matrixType") == 0) {
	  if (Tcl_GetInt(interp, argv[currentArg+1], &matType) != TCL_OK)
		  opserr << "Mumps Warning: failed to get -matrixType. Unsymmetric matrix assumed\n";
	  if (matType < 0 || matType > 2) {
		  opserr << "Mumps Warning: wrong -matrixType value (" << matType << "). Unsymmetric matrix assumed\n";
		  matType = 0;
	  }
	  currentArg += 2;
	} else 
	  currentArg++;
      }    
    }

#ifdef _PARALLEL_PROCESSING
    MumpsParallelSolver *theSolver = new MumpsParallelSolver(icntl7, icntl14);
    theSOE = new MumpsParallelSOE(*theSolver);
#elif _PARALLEL_INTERPRETERS
    MumpsParallelSolver *theSolver = new MumpsParallelSolver(icntl7, icntl14);
    MumpsParallelSOE *theParallelSOE = new MumpsParallelSOE(*theSolver, matType);
    theParallelSOE->setProcessID(OPS_rank);
    theParallelSOE->setChannels(numChannels, theChannels);
    theSOE = theParallelSOE;
#else
    MumpsSolver *theSolver = new MumpsSolver(icntl7, icntl14);
    theSOE = new MumpsSOE(*theSolver, matType);
#endif

  }

#endif

  
  else {

    //
    // maybe a package
    //

    // try existing loaded packages
    ExternalClassFunction  *solverCommands = theExternalSolverCommands;
    bool found = false;
    //    int result = TCL_ERROR;
    while (solverCommands != NULL && found == false) {

      if (strcmp(argv[1], solverCommands->funcName) == 0) {
	
          OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
	void *theRes = (*(solverCommands->funcPtr))();
	if (theRes != 0) {

	  theSOE = (LinearSOE *)theRes;
	  found = true;
	}
      } else
	solverCommands = solverCommands->next;
    }

    //
    // if not there try loading package
    //

    if (found == false) {

      void *libHandle;
      void *(*funcPtr)();
      int solverNameLength = strlen(argv[1]);
      char *tclFuncName = new char[solverNameLength+5];
      strcpy(tclFuncName, "OPS_");
      strcpy(&tclFuncName[4], argv[1]);    

      int res = getLibraryFunction(argv[1], tclFuncName, &libHandle, (void **)&funcPtr);
      
      delete [] tclFuncName;
      
      if (res == 0) {

	char *solverName = new char[solverNameLength+1];
	strcpy(solverName, argv[1]);
	ExternalClassFunction *theSolverCommand = new ExternalClassFunction;
	theSolverCommand->funcPtr = funcPtr;
	theSolverCommand->funcName = solverName;	
	theSolverCommand->next = theExternalSolverCommands;
	theExternalSolverCommands = theSolverCommand;
	
    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);
	
	void *theRes = (*funcPtr)();
	if (theRes != 0) {
	  theSOE = (LinearSOE *)theRes;
	}
      }
    }
  }
    
  // if the analysis exists - we want to change the SOEif

  if (theSOE != 0) {
    if (theStaticAnalysis != 0)
      theStaticAnalysis->setLinearSOE(*theSOE);
    if (theTransientAnalysis != 0)
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
    
    return TCL_OK;
  }

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
  } else if (strcmp(argv[1],"AMD") == 0) {
    AMD *theAMD = new AMD();	
    theNumberer = new DOF_Numberer(*theAMD);    	
  } 

#ifdef _PARALLEL_INTERPRETERS

  else if ((strcmp(argv[1],"ParallelPlain") == 0) || (strcmp(argv[1],"Parallel") == 0)) {
    ParallelNumberer *theParallelNumberer = new ParallelNumberer;
    theNumberer = theParallelNumberer;       
    theParallelNumberer->setProcessID(OPS_rank);
    theParallelNumberer->setChannels(numChannels, theChannels);
  } else if (strcmp(argv[1],"ParallelRCM") == 0) {
    RCM *theRCM = new RCM(false);	
    ParallelNumberer *theParallelNumberer = new ParallelNumberer(*theRCM);    	
    theNumberer = theParallelNumberer;       
    theParallelNumberer->setProcessID(OPS_rank);
    theParallelNumberer->setChannels(numChannels, theChannels);
  }   

#endif

  else {
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
  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);

  // check argv[1] for type of Algorithm and create the object
  if (strcmp(argv[1],"Linear") == 0) {
    int formTangent = CURRENT_TANGENT;
    int factorOnce = 0;
    int count = 2;
    while (count < argc) {
      if ((strcmp(argv[count],"-secant") == 0) || (strcmp(argv[count],"-Secant") == 0)) {
	formTangent = CURRENT_SECANT;
      } else if ((strcmp(argv[count],"-initial") == 0) || (strcmp(argv[count],"-Initial") == 0)) {
	formTangent = INITIAL_TANGENT;
      } else if ((strcmp(argv[count],"-factorOnce") == 0) || (strcmp(argv[count],"-FactorOnce") ==0 )) {
	factorOnce = 1;
      }
      count++;
    }
    theNewAlgo = new Linear(formTangent, factorOnce);
  }

  else if (strcmp(argv[1],"Newton") == 0) {
    void *theNewtonAlgo = OPS_NewtonRaphsonAlgorithm();
    if (theNewtonAlgo == 0)
      return TCL_ERROR;

    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
    if (theTest != 0)
      theNewAlgo->setConvergenceTest(theTest);
  }

  else if ((strcmp(argv[1],"NewtonHallM") == 0) || (strcmp(argv[1],"NewtonHall") == 0)) {
    void *theNewtonAlgo = OPS_NewtonHallM();
    if (theNewtonAlgo == 0)
      return TCL_ERROR;

    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
    if (theTest != 0)
      theNewAlgo->setConvergenceTest(theTest);
  }

  else if (strcmp(argv[1],"ModifiedNewton") == 0) {
    void *theNewtonAlgo = OPS_ModifiedNewton();
    if (theNewtonAlgo == 0)
      return TCL_ERROR;

    theNewAlgo = (EquiSolnAlgo *)theNewtonAlgo;
    if (theTest != 0)
      theNewAlgo->setConvergenceTest(theTest);
  }

  else if (strcmp(argv[1],"KrylovNewton") == 0) {
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i],"-iterate") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  iterateTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  iterateTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  iterateTangent = NO_TANGENT;
      } 
      else if (strcmp(argv[i],"-increment") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  incrementTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  incrementTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  incrementTangent = NO_TANGENT;
      }
      else if (strcmp(argv[i],"-maxDim") == 0 && i+1 < argc) {
	i++;
	maxDim = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }

    Accelerator *theAccel;
    theAccel = new KrylovAccelerator(maxDim, iterateTangent);

    theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  }

  else if (strcmp(argv[1],"RaphsonNewton") == 0) {
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i],"-iterate") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  iterateTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  iterateTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  iterateTangent = NO_TANGENT;
      } 
      else if (strcmp(argv[i],"-increment") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  incrementTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  incrementTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  incrementTangent = NO_TANGENT;
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }

    Accelerator *theAccel;
    theAccel = new RaphsonAccelerator(iterateTangent);

    theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  }

  else if (strcmp(argv[1],"MillerNewton") == 0) {
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;

    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i],"-iterate") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  iterateTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  iterateTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  iterateTangent = NO_TANGENT;
      } 
      else if (strcmp(argv[i],"-increment") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  incrementTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  incrementTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  incrementTangent = NO_TANGENT;
      }
      else if (strcmp(argv[i],"-maxDim") == 0 && i+1 < argc) {
	i++;
	maxDim = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }

    Accelerator *theAccel = 0;
    //theAccel = new MillerAccelerator(maxDim, 0.01, iterateTangent);

    theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  }

  else if (strcmp(argv[1],"SecantNewton") == 0) {
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i],"-iterate") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  iterateTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  iterateTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  iterateTangent = NO_TANGENT;
      } 
      else if (strcmp(argv[i],"-increment") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  incrementTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  incrementTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  incrementTangent = NO_TANGENT;
      }
      else if (strcmp(argv[i],"-maxDim") == 0 && i+1 < argc) {
	i++;
	maxDim = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }

    Accelerator *theAccel;
    theAccel = new SecantAccelerator2(maxDim, iterateTangent); 

    theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
  }

  else if (strcmp(argv[1],"PeriodicNewton") == 0) {
    int incrementTangent = CURRENT_TANGENT;
    int iterateTangent = CURRENT_TANGENT;
    int maxDim = 3;
    for (int i = 2; i < argc; i++) {
      if (strcmp(argv[i],"-iterate") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  iterateTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  iterateTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  iterateTangent = NO_TANGENT;
      } 
      else if (strcmp(argv[i],"-increment") == 0 && i+1 < argc) {
	i++;
	if (strcmp(argv[i],"current") == 0)
	  incrementTangent = CURRENT_TANGENT;
	if (strcmp(argv[i],"initial") == 0)
	  incrementTangent = INITIAL_TANGENT;
	if (strcmp(argv[i],"noTangent") == 0)
	  incrementTangent = NO_TANGENT;
      }
      else if (strcmp(argv[i],"-maxDim") == 0 && i+1 < argc) {
	i++;
	maxDim = atoi(argv[i]);
      }
    }

    if (theTest == 0) {
      opserr << "ERROR: No ConvergenceTest yet specified\n";
      return TCL_ERROR;	  
    }

    Accelerator *theAccel;
    theAccel = new PeriodicAccelerator(maxDim, iterateTangent); 

    theNewAlgo = new AcceleratedNewton(*theTest, theAccel, incrementTangent);
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

  else if (strcmp(argv[1],"ExpressNewton") == 0) {
    int nIter = 2, factorOnce = 0, formTangent = CURRENT_TANGENT;
    double kMultiplier = 1.0;
	  if (argc >= 3 && Tcl_GetInt(interp, argv[2], &nIter) != TCL_OK)
      return TCL_ERROR;
	  if (argc >= 4 && Tcl_GetDouble(interp, argv[3], &kMultiplier) != TCL_OK)
      return TCL_ERROR;
    int count = 4;
    while (argc > count) {
      if ((strcmp(argv[count],"-initialTangent") == 0) || (strcmp(argv[count],"-InitialTangent") == 0)) {
        formTangent = INITIAL_TANGENT;
      } else if ((strcmp(argv[count],"-currentTangent") == 0) || (strcmp(argv[count],"-CurrentTangent") ==0 )) {
        formTangent = CURRENT_TANGENT;
      } else if ((strcmp(argv[count],"-factorOnce") == 0) || (strcmp(argv[count],"-FactorOnce") ==0 )) {
        factorOnce = 1;
      }
      count++;
    }
    theNewAlgo = new ExpressNewton(nIter,kMultiplier,formTangent,factorOnce);
  }

  else {
    opserr << "WARNING No EquiSolnAlgo type " << argv[1] << " exists\n";
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
  double tol2 = 0.0;
  double tolp = 0.0;
  double tolp2 = 0.0;
  double tolrel = 0.0;
  double tolprel = 0.0;
  double maxTol = OPS_MAXTOL;

  int numIter = 0;
  int printIt = 0;
  int normType = 2;
  int maxIncr = -1;

  if ((strcmp(argv[1],"NormDispAndUnbalance") == 0) || 
      (strcmp(argv[1],"NormDispOrUnbalance") == 0)) {
    if (argc == 5) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetDouble(interp, argv[3], &tol2) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
    } else if (argc == 6) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetDouble(interp, argv[3], &tol2) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[5], &printIt) != TCL_OK)	
	return TCL_ERROR;			  
    } else if (argc == 7) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetDouble(interp, argv[3], &tol2) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[5], &printIt) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[6], &normType) != TCL_OK)	
	return TCL_ERROR;			  
    } else if (argc == 8) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetDouble(interp, argv[3], &tol2) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[5], &printIt) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[6], &normType) != TCL_OK)	
	return TCL_ERROR;
      if (Tcl_GetInt(interp, argv[7], &maxIncr) != TCL_OK)	
	return TCL_ERROR;
    }

  } else if (strcmp(argv[1],"PFEM") == 0) {
      if(argc > 8) {
          if(Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
              return TCL_ERROR;
          if(Tcl_GetDouble(interp, argv[3], &tolp) != TCL_OK)	
              return TCL_ERROR;
          if(Tcl_GetDouble(interp, argv[4], &tol2) != TCL_OK)	
              return TCL_ERROR;
          if(Tcl_GetDouble(interp, argv[5], &tolp2) != TCL_OK)	
              return TCL_ERROR;
          if(Tcl_GetDouble(interp, argv[6], &tolrel) != TCL_OK)	
              return TCL_ERROR;
          if(Tcl_GetDouble(interp, argv[7], &tolprel) != TCL_OK)	
              return TCL_ERROR;
          if(Tcl_GetInt(interp, argv[8], &numIter) != TCL_OK)	
              return TCL_ERROR;
      }
      if(argc > 9) {
          if(Tcl_GetInt(interp, argv[9], &maxIncr) != TCL_OK)	
              return TCL_ERROR;
      }
      if(argc > 10) {
          if(Tcl_GetInt(interp, argv[10], &printIt) != TCL_OK)	
              return TCL_ERROR;
      }
      if(argc > 11) {
          if(Tcl_GetInt(interp, argv[11], &normType) != TCL_OK)	
              return TCL_ERROR;
      }

  } else if (strcmp(argv[1],"FixedNumIter") == 0) {

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
    } else if (argc == 6) {
      if (Tcl_GetInt(interp, argv[2], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[3], &printIt) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &normType) != TCL_OK)	
	return TCL_ERROR;			  		  
      if (Tcl_GetDouble(interp, argv[5], &maxTol) != TCL_OK)	
	return TCL_ERROR;			  		  
    }    

  } else {
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
    } else if (argc == 7) {
      if (Tcl_GetDouble(interp, argv[2], &tol) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[4], &printIt) != TCL_OK)	
	return TCL_ERROR;			  
      if (Tcl_GetInt(interp, argv[5], &normType) != TCL_OK)	
	return TCL_ERROR;
      if (Tcl_GetDouble(interp, argv[6], &maxTol) != TCL_OK)	
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
      theNewTest = new CTestNormUnbalance(tol,numIter,printIt,normType,maxIncr, maxTol);       
    else if (strcmp(argv[1],"NormDispIncr") == 0) 
      theNewTest = new CTestNormDispIncr(tol,numIter,printIt,normType, maxTol);             
    else if (strcmp(argv[1],"NormDispAndUnbalance") == 0) 
        theNewTest = new NormDispAndUnbalance(tol,tol2, numIter,printIt,normType,maxIncr);       
    else if (strcmp(argv[1],"NormDispOrUnbalance") == 0) 
        theNewTest = new NormDispOrUnbalance(tol,tol2, numIter,printIt,normType,maxIncr);       
    else if (strcmp(argv[1],"EnergyIncr") == 0) 
      theNewTest = new CTestEnergyIncr(tol,numIter,printIt,normType, maxTol);             
    else if (strcmp(argv[1],"RelativeNormUnbalance") == 0) 
      theNewTest = new CTestRelativeNormUnbalance(tol,numIter,printIt,normType);       
    else if (strcmp(argv[1],"RelativeNormDispIncr") == 0) 
      theNewTest = new CTestRelativeNormDispIncr(tol,numIter,printIt,normType);             
    else if (strcmp(argv[1],"RelativeEnergyIncr") == 0) 
      theNewTest = new CTestRelativeEnergyIncr(tol,numIter,printIt,normType);             
    else if (strcmp(argv[1],"RelativeTotalNormDispIncr") == 0) 
      theNewTest = new CTestRelativeTotalNormDispIncr(tol,numIter,printIt,normType);             
    else if (strcmp(argv[1],"PFEM") == 0) 
        theNewTest = new CTestPFEM(tol,tolp,tol2,tolp2,tolrel,tolprel,numIter,maxIncr,printIt,normType);
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

    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, &theDomain);

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
      } else if (strcmp(argv[1],"StagedLoadControl") == 0) {
      double dLambda;
      double minIncr, maxIncr;
      int numIter;
      if (argc < 3) {
	opserr << "WARNING incorrect # args - integrator StagedLoadControl dlam <Jd dlamMin dlamMax>\n";
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
      theStaticIntegrator = new StagedLoadControl(dLambda, numIter, minIncr, maxIncr);       

  
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

  else if (strcmp(argv[1], "EQPath") == 0) {
		double arcLength;
		int type;
		int numIter;
		if (argc != 4) {
			opserr << "WARNING integrator EQPath $arc_length $type \n";
			opserr << "REFS : \n";
			opserr << " https://doi.org/10.12989/sem.2013.48.6.849	 \n";
			opserr << " https://doi.org/10.12989/sem.2013.48.6.879	 \n";
			return TCL_ERROR;
		}

		if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)
		{
			opserr << "WARNING integrator EQPath $arc_length $type \n";
			opserr << " https://doi.org/10.12989/sem.2013.48.6.849	 \n";
			opserr << " https://doi.org/10.12989/sem.2013.48.6.879	 \n";
			return TCL_ERROR;
			return TCL_ERROR;
		}

		if (Tcl_GetInt(interp, argv[3], &type) != TCL_OK)
		{
			opserr << "WARNING integrator $arc_length $type \n";
			opserr << "$type = 1 Minimum Residual Displacement \n";
			opserr << "$type = 2 Normal Plain \n";
			opserr << "$type = 3 Update Normal Plain \n";
			opserr << "$type = 4 Cylindrical Arc-Length \n";

			return TCL_ERROR;
		}

		theStaticIntegrator = new EQPath(arcLength, type);

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
      int tangFlag = 0;
      
      if (Tcl_GetInt(interp, argv[2], &node) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK)	
	return TCL_ERROR;	
      if (Tcl_GetDouble(interp, argv[4], &increment) != TCL_OK)	
	return TCL_ERROR;

      if (argc == 6 || argc == 9)
	if (argc == 6) {
	  if (strcmp(argv[5],"-initial") == 0)
	    tangFlag = 1;
	} else if (strcmp(argv[8],"-initial") == 0)
	  tangFlag = 1;

      if (argc > 6) {
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



#ifdef _PARALLEL_PROCESSING

      theStaticIntegrator = new DistributedDisplacementControl(node,dof-1,increment,
							       numIter, minIncr, maxIncr);
#else
      Node *theNode = theDomain.getNode(node);
      if (theNode == 0) {
	opserr << "WARNING integrator DisplacementControl node dof dU : Node does not exist\n";
	return TCL_ERROR;	  
      }


      int numDOF = theNode->getNumberDOF();
      if (dof <= 0 || dof > numDOF) {
	opserr << "WARNING integrator DisplacementControl node dof dU : invalid dof given\n";
	return TCL_ERROR;	  
      }

      theStaticIntegrator = new DisplacementControl(node, dof-1, increment, &theDomain,
						    numIter, minIncr, maxIncr, tangFlag);
#endif

      // if the analysis exists - we want to change the Integrator
      if (theStaticAnalysis != 0) 
	theStaticAnalysis->setIntegrator(*theStaticIntegrator);
  }  


#ifdef _PARALLEL_INTERPRETERS

  else if ((strcmp(argv[1],"ParallelDisplacementControl") == 0) || (strcmp(argv[1],"ParallelDisplacementControl") == 0)) {
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


      DistributedDisplacementControl *theDDC  = new DistributedDisplacementControl(node,dof-1,increment,
										   numIter, minIncr, maxIncr);

      theDDC->setProcessID(OPS_rank);
      theDDC->setChannels(numChannels, theChannels);
      theStaticIntegrator = theDDC;

      // if the analysis exists - we want to change the Integrator
      if (theStaticAnalysis != 0)
	theStaticAnalysis->setIntegrator(*theStaticIntegrator);
  }
#endif

  else if ((strcmp(argv[1],"TRBDF2") == 0) || (strcmp(argv[1],"Bathe") == 0)) {
    theTransientIntegrator = new TRBDF2();           
  }

  else if ((strcmp(argv[1],"TRBDF3") == 0) || (strcmp(argv[1],"Bathe3") == 0)) {
      theTransientIntegrator = new TRBDF3();     
  }

  else if (strcmp(argv[1],"Houbolt") == 0) {
      theTransientIntegrator = new Houbolt();     
  }
  
  /*else if (strcmp(argv[1],"ParkLMS3") == 0) {
      theTransientIntegrator = new ParkLMS3();     
  }*/
  
  else if (strcmp(argv[1],"BackwardEuler") == 0) {
      int optn = 0;
      if (argc == 3) {
          if (Tcl_GetInt(interp, argv[2], &optn) != TCL_OK) {
              opserr << "WARNING integrator BackwardEuler <option> - undefined option specified\n";	  
              return TCL_ERROR;	
          }
      }
      theTransientIntegrator = new BackwardEuler(optn);     
  }

  
  else if (strcmp(argv[1],"Newmark") == 0) {
    theTransientIntegrator = (TransientIntegrator*)OPS_Newmark();
    
    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"GimmeMCK") == 0 || strcmp(argv[1],"ZZTop") == 0) {
    theTransientIntegrator = (TransientIntegrator*)OPS_GimmeMCK();
   }
  else if (strcmp(argv[1],"StagedNewmark") == 0) {
    theTransientIntegrator = (TransientIntegrator*)OPS_StagedNewmark();
    
    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"PFEM") == 0) {
    theTransientIntegrator = new PFEMIntegrator();

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  } 
  
  else if (strcmp(argv[1],"NewmarkExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkExplicit();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"NewmarkHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSIncrReduct();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"NewmarkHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSIncrLimit();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"NewmarkHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSFixedNumIter();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
// #ifdef _RELIABILITY
//   else if (strcmp(argv[1],"NewmarkWithSensitivity") == 0) {
// 	  int assemblyFlag = 0;
//       double gamma;
//       double beta;
//       double alphaM, betaK, betaKi, betaKc;
//       if (argc != 4 && argc != 6 && argc != 8 && argc != 10) {
// 	     interp->result = "WARNING integrator Newmark gamma beta <alphaM?  betaKcurrent?  betaKi? betaKlastCommitted?> <-assemble tag?> ";
// 	     return TCL_ERROR;
//       }
	  
// 	  // Take care of argc == 4, the basic case
//       if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
// 		  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;	  
// 		  return TCL_ERROR;	
//       }
//       if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
// 		  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 		  return TCL_ERROR;	
//       }

// 	  // If only assembly flag is given extra
// 	  if (argc == 6) {
// 		  if (strcmp(argv[4],"-assemble") != 0) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 		  }
// 		  if (Tcl_GetInt(interp, argv[5], &assemblyFlag) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 	  }
// 	  // If only extra integrator (damping) parameters are given extra
//       if (argc == 8) {
// 		  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
//       }
// 	  // If everything is given extra
// 	  if (argc == 10) {
// 		  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (strcmp(argv[8],"-assemble") != 0) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 		  }
// 		  if (Tcl_GetInt(interp, argv[9], &assemblyFlag) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 	  }

//       if (argc == 4 || argc == 6) {
// 	theNSI = new NewmarkSensitivityIntegrator(assemblyFlag,gamma,beta);       
//       }
//       else {
// 	theNSI = new NewmarkSensitivityIntegrator(assemblyFlag,gamma,beta,alphaM,betaK,betaKi,betaKc);
//       }
//       theTransientIntegrator = theNSI;


//       // if the analysis exists - we want to change the Integrator
// 	  if (theTransientAnalysis != 0)
// 		theTransientAnalysis->setIntegrator(*theTransientIntegrator);
//   }  

//   else if (strcmp(argv[1],"NewNewmarkWithSensitivity") == 0) {
// 	  int assemblyFlag = 0;
//       double gamma;
//       double beta;
//       double alphaM, betaK, betaKi, betaKc;
//       if (argc != 4 && argc != 6 && argc != 8 && argc != 10) {
// 	     interp->result = "WARNING integrator Newmark gamma beta <alphaM?  betaKcurrent?  betaKi? betaKlastCommitted?> <-assemble tag?> ";
// 	     return TCL_ERROR;
//       }
	  
// 	  // Take care of argc == 4, the basic case
//       if (Tcl_GetDouble(interp, argv[2], &gamma) != TCL_OK) {
// 		  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;	  
// 		  return TCL_ERROR;	
//       }
//       if (Tcl_GetDouble(interp, argv[3], &beta) != TCL_OK) {
// 		  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 		  return TCL_ERROR;	
//       }

// 	  // If only assembly flag is given extra
// 	  if (argc == 6) {
// 		  if (strcmp(argv[4],"-assemble") != 0) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 		  }
// 		  if (Tcl_GetInt(interp, argv[5], &assemblyFlag) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 	  }
// 	  // If only extra integrator (damping) parameters are given extra
//       if (argc == 8) {
// 		  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
//       }
// 	  // If everything is given extra
// 	  if (argc == 10) {
// 		  if (Tcl_GetDouble(interp, argv[4], &alphaM) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[5], &betaK) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[6], &betaKi) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (Tcl_GetDouble(interp, argv[7], &betaKc) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 		  if (strcmp(argv[8],"-assemble") != 0) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 		  }
// 		  if (Tcl_GetInt(interp, argv[9], &assemblyFlag) != TCL_OK) {
// 			  opserr << "WARNING: Error in input to Newmark sensitivity integrator" << endln;
// 			  return TCL_ERROR;	
// 		  }
// 	  }

//       if (argc == 4 || argc == 6) {
// 	theNNSI = new NewNewmarkSensitivityIntegrator(assemblyFlag,gamma,beta);       
//       }
//       else {
// 	theNNSI = new NewNewmarkSensitivityIntegrator(assemblyFlag,gamma,beta,alphaM,betaK,betaKi,betaKc);
//       }
//       theTransientIntegrator = theNNSI;

// 	  //// added by K Fujimura
// 	  if (theTransientAnalysis != 0){
// 	    opserr << "For the TransientAnalysis, the integrator must be \n";
// 	    opserr << "NewmarkSensitivityIntegrator \n";
// 	    return TCL_ERROR;	
// 	  }
//       // if the analysis exists - we want to change the Integrator
// 	  if (theReliabilityTransientAnalysis != 0)
// 		theReliabilityTransientAnalysis->setIntegrator(*theTransientIntegrator);
//   }  
  
//   else if(strcmp(argv[1], "PFEMWithSensitivity") == 0) {
//       int flag = 0;
//       if(argc > 4) {
//           if(strcmp(argv[2],"-assemble") != TCL_OK) {
//               opserr<<"WARNING: Error in input to PFEMSensitivityIntegrator\n";
//               return TCL_ERROR;
//           }
//           if(Tcl_GetInt(interp, argv[3], &flag) != TCL_OK) {
//               opserr<<"WARNING: Error in input to PFEMSensitivityIntegrator\n";
//               return TCL_ERROR;	
//           }
//       }

//       thePFEMSI = new PFEMSensitivityIntegrator(flag);
//       theTransientIntegrator = thePFEMSI;
//       if (theTransientAnalysis != 0) {
//           theTransientAnalysis->setIntegrator(*theTransientIntegrator);
//       }
//   }

// #endif
  
  else if (strcmp(argv[1],"HHT") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHT();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHT_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHT_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTGeneralized") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralized();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTGeneralized_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralized_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTExplicit();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTExplicit_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTGeneralizedExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralizedExplicit();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTGeneralizedExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralizedExplicit_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrLimit();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTHSIncrLimit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrLimit_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrReduct();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTHSIncrReduct_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrReduct_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSFixedNumIter();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"HHTHSFixedNumIter_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSFixedNumIter_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"GeneralizedAlpha") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_GeneralizedAlpha();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"KRAlphaExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_KRAlphaExplicit();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"KRAlphaExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_KRAlphaExplicit_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"AlphaOS") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOS();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"AlphaOS_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOS_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"AlphaOSGeneralized") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOSGeneralized();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"AlphaOSGeneralized_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOSGeneralized_TP();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"Collocation") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_Collocation();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"CollocationHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSIncrReduct();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"CollocationHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSIncrLimit();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"CollocationHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSFixedNumIter();
    
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
  
  else if (strcmp(argv[1],"WilsonTheta") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_WilsonTheta();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"CentralDifference") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifference();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"CentralDifferenceAlternative") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifferenceAlternative();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"CentralDifferenceNoDamping") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifferenceNoDamping();
    
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }
  
  else if (strcmp(argv[1],"Transient") == 0) {

    theTransientIntegrator = 0;
    
    // try existing loaded packages
    ExternalClassFunction  *integratorCommands = theExternalTransientIntegratorCommands;
    bool found = false;
    //    int result = TCL_ERROR;
    while (integratorCommands != NULL && found == false) {

      if (strcmp(argv[2], integratorCommands->funcName) == 0) {
	
          OPS_ResetInputNoBuilder(clientData, interp, 3, argc, argv, &theDomain);
	void *theRes = (*(integratorCommands->funcPtr))();
	if (theRes != 0) {
	  theTransientIntegrator = (TransientIntegrator *)theRes;
	  found = true;
	}
      } else
	integratorCommands = integratorCommands->next;
    }

    //
    // if not there try loading package
    //

    if (found == false) {

      void *libHandle;
      void *(*funcPtr)();
      int integratorNameLength = strlen(argv[2]);
      char *tclFuncName = new char[integratorNameLength+5];
      strcpy(tclFuncName, "OPS_");
      strcpy(&tclFuncName[4], argv[2]);    

      int res = getLibraryFunction(argv[2], tclFuncName, &libHandle, (void **)&funcPtr);
      
      delete [] tclFuncName;
      
      if (res == 0) {

	char *integratorName = new char[integratorNameLength+1];
	strcpy(integratorName, argv[2]);
	ExternalClassFunction *theIntegratorCommand = new ExternalClassFunction;
	theIntegratorCommand->funcPtr = funcPtr;
	theIntegratorCommand->funcName = integratorName;	
	theIntegratorCommand->next = theExternalTransientIntegratorCommands;
	theExternalTransientIntegratorCommands = theIntegratorCommand;
	
    OPS_ResetInputNoBuilder(clientData, interp, 3, argc, argv, &theDomain);
	
	void *theRes = (*funcPtr)();
	if (theRes != 0) {
	  theTransientIntegrator = (TransientIntegrator *)theRes;
	}
      }
    }

    if (theTransientIntegrator == 0) {
      opserr << "Transient Integrator Not Found \n";
      return TCL_ERROR;      
    }

    // if the analysis exists - we want to change the Integrator
    if (theTransientAnalysis != 0)
      theTransientAnalysis->setIntegrator(*theTransientIntegrator);
  }

  else if (strcmp(argv[1],"Static") == 0) {

    theStaticIntegrator = 0;

    // try existing loaded packages
    ExternalClassFunction  *integratorCommands = theExternalStaticIntegratorCommands;
    bool found = false;

    while (integratorCommands != NULL && found == false) {

      if (strcmp(argv[2], integratorCommands->funcName) == 0) {
	
    OPS_ResetInputNoBuilder(clientData, interp, 3, argc, argv, &theDomain);
	void *theRes = (*(integratorCommands->funcPtr))();
	if (theRes != 0) {
	  theStaticIntegrator = (StaticIntegrator *)theRes;
	  found = true;
	}
      } else
	integratorCommands = integratorCommands->next;
    }

    //
    // if not there try loading package
    //

    if (found == false) {

      void *libHandle;
      void *(*funcPtr)();
      int integratorNameLength = strlen(argv[2]);
      char *tclFuncName = new char[integratorNameLength+5];
      strcpy(tclFuncName, "OPS_");
      strcpy(&tclFuncName[4], argv[2]);    

      int res = getLibraryFunction(argv[2], tclFuncName, &libHandle, (void **)&funcPtr);
      
      delete [] tclFuncName;
      
      if (res == 0) {

	char *integratorName = new char[integratorNameLength+1];
	strcpy(integratorName, argv[2]);
	ExternalClassFunction *theIntegratorCommand = new ExternalClassFunction;
	theIntegratorCommand->funcPtr = funcPtr;
	theIntegratorCommand->funcName = integratorName;	
	theIntegratorCommand->next = theExternalStaticIntegratorCommands;
	theExternalStaticIntegratorCommands = theIntegratorCommand;
	
    OPS_ResetInputNoBuilder(clientData, interp, 3, argc, argv, &theDomain);
	
	void *theRes = (*funcPtr)();
	if (theRes != 0) {
	  theStaticIntegrator = (StaticIntegrator *)theRes;
	}
      }
    }

    if (theStaticIntegrator == 0) {
      opserr << "Static Integrator Not Found \n";
      return TCL_ERROR;      
    }

    // if the analysis exists - we want to change the Integrator
    if (theStaticAnalysis != 0)
      theStaticAnalysis->setIntegrator(*theStaticIntegrator);
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
    theIntegrator = theTransientIntegrator;
    
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
eigenAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, 
	      TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - eigen <type> numModes?\n";
    return TCL_ERROR;
  }    
  
  bool generalizedAlgo = true; // 0 - frequency/generalized (default),1 - standard, 2 - buckling
  int typeSolver = EigenSOE_TAGS_ArpackSOE;
  int loc = 1;
  double shift = 0.0;
  bool findSmallest = true;
  
  // Check type of eigenvalue analysis
  while (loc < (argc-1)) {
    if ((strcmp(argv[loc],"frequency") == 0) || 
	(strcmp(argv[loc],"-frequency") == 0) ||
	(strcmp(argv[loc],"generalized") == 0) ||
	(strcmp(argv[loc],"-generalized") == 0))
      generalizedAlgo = true;
    
    else if ((strcmp(argv[loc],"standard") == 0) || 
	     (strcmp(argv[loc],"-standard") == 0))
      generalizedAlgo = false;

    else if ((strcmp(argv[loc],"-findLargest") == 0))
      findSmallest = false;
    
    else if ((strcmp(argv[loc],"genBandArpack") == 0) || 
         (strcmp(argv[loc],"-genBandArpack") == 0) ||
         (strcmp(argv[loc],"genBandArpackEigen") == 0) || 
         (strcmp(argv[loc],"-genBandArpackEigen") == 0))
      typeSolver = EigenSOE_TAGS_ArpackSOE;
    
    else if ((strcmp(argv[loc],"symmBandLapack") == 0) || 
         (strcmp(argv[loc],"-symmBandLapack") == 0) ||
         (strcmp(argv[loc],"symmBandLapackEigen") == 0) || 
         (strcmp(argv[loc],"-symmBandLapackEigen") == 0))
      typeSolver = EigenSOE_TAGS_SymBandEigenSOE;
    
    else if ((strcmp(argv[loc],"fullGenLapack") == 0) || 
         (strcmp(argv[loc],"-fullGenLapack") == 0) ||
         (strcmp(argv[loc],"fullGenLapackEigen") == 0) || 
         (strcmp(argv[loc],"-fullGenLapackEigen") == 0))
      typeSolver = EigenSOE_TAGS_FullGenEigenSOE;
    
    else {
      opserr << "eigen - unknown option specified " << argv[loc] << endln;
    }
    
    loc++;
  }
  
    // check argv[loc] for number of modes
  //    int numEigen;
    if ((Tcl_GetInt(interp, argv[loc], &numEigen) != TCL_OK) || numEigen < 0) {
      opserr << "WARNING eigen numModes?  - illegal numModes\n";    
      return TCL_ERROR;	
    }

    //
    // create a transient analysis if no analysis exists
    //

    if (theStaticAnalysis == 0 && theTransientAnalysis == 0) {

	if (theAnalysisModel == 0) 
	    theAnalysisModel = new AnalysisModel();
	if (theTest == 0) 
	  theTest = new CTestNormUnbalance(1.0e-6,25,0);       
	if (theAlgorithm == 0) {
	  theAlgorithm = new NewtonRaphson(*theTest); 
	}
	if (theHandler == 0) {
	  theHandler = new TransformationConstraintHandler();       
	}
	if (theNumberer == 0) {
	  RCM *theRCM = new RCM(false);	
	  theNumberer = new DOF_Numberer(*theRCM);    	
	}
	if (theTransientIntegrator == 0) {
	  theTransientIntegrator = new Newmark(0.5,0.25);       
	}
	if (theSOE == 0) {
	  ProfileSPDLinSolver *theSolver;
	  theSolver = new ProfileSPDLinDirectSolver(); 	
#ifdef _PARALLEL_PROCESSING
	  theSOE = new DistributedProfileSPDLinSOE(*theSolver);
#else
	  theSOE = new ProfileSPDLinSOE(*theSolver);      
#endif
	}
	
	theTransientAnalysis = new DirectIntegrationAnalysis(theDomain,
							     *theHandler,
							     *theNumberer,
							     *theAnalysisModel,
							     *theAlgorithm,
							     *theSOE,
							     *theTransientIntegrator,
							     theTest);
    }

    //
    // create a new eigen system and solver
    //

    bool setEigen = false;
    if (theEigenSOE != 0) {
      if (theEigenSOE->getClassTag() != typeSolver) {
	//	delete theEigenSOE;
	theEigenSOE = 0;
	setEigen = true;
      }
    }

	
    if (theEigenSOE == 0) {
                       
      if (typeSolver == EigenSOE_TAGS_SymBandEigenSOE) {
	SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver(); 
	theEigenSOE = new SymBandEigenSOE(*theEigenSolver, *theAnalysisModel); 
	
      } else if (typeSolver == EigenSOE_TAGS_FullGenEigenSOE) {
	
	FullGenEigenSolver *theEigenSolver = new FullGenEigenSolver();
	theEigenSOE = new FullGenEigenSOE(*theEigenSolver, *theAnalysisModel);

      } else {

	theEigenSOE = new ArpackSOE(shift);    

      }
      
      //
      // set the eigen soe in the system
      //

      if (theStaticAnalysis != 0) {
	theStaticAnalysis->setEigenSOE(*theEigenSOE);
      } else if (theTransientAnalysis != 0) {
	theTransientAnalysis->setEigenSOE(*theEigenSOE);
      }
 


#ifdef _PARALLEL_PROCESSING

      if (OPS_PARTITIONED == false && OPS_NUM_SUBDOMAINS > 1) {
	if (partitionModel(0) < 0) {
	  opserr << "WARNING before analysis; partition failed - too few elements\n";
	  OpenSeesExit(clientData, interp, argc, argv);
	  return TCL_ERROR;
	}
      }

      if (theStaticAnalysis != 0 || theTransientAnalysis != 0) {
	SubdomainIter &theSubdomains = theDomain.getSubdomains();
	Subdomain *theSub;
	while ((theSub = theSubdomains()) != 0) {
	  theSub->setAnalysisEigenSOE(*theEigenSOE);
	}
      }
#endif

    } // theEigenSOE != 0    


    int requiredDataSize = 40*numEigen;
    if (requiredDataSize > resDataSize) {
      if (resDataPtr != 0) {
	delete [] resDataPtr;
      }
      resDataPtr = new char[requiredDataSize];
      resDataSize = requiredDataSize;
    }
    
    for (int i=0; i<requiredDataSize; i++)
      resDataPtr[i] = '\n';
    
    int result = 0;

    if (theStaticAnalysis != 0) {
      result = theStaticAnalysis->eigen(numEigen, generalizedAlgo, findSmallest);      
    } else if (theTransientAnalysis != 0) {
      result = theTransientAnalysis->eigen(numEigen, generalizedAlgo, findSmallest);      
    }

    if (result == 0) {
      //      char *eigenvalueS = new char[15 * numEigen];    
      const Vector &eigenvalues = theDomain.getEigenvalues();
      int cnt = 0;
      for (int i=0; i<numEigen; i++) {
	cnt += sprintf(&resDataPtr[cnt], "%35.20f  ", eigenvalues[i]);
      }
      
      Tcl_SetResult(interp, resDataPtr, TCL_STATIC);
    }
    
    return TCL_OK;
}

int 
modalProperties(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, &theDomain);
    OPS_DomainModalProperties();
    return TCL_OK;
}

int
responseSpectrum(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, &theDomain);
    OPS_ResponseSpectrumAnalysis();
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


extern bool OPS_removeTimeSeries(int tag);

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
      thePattern->clearAll();
      delete thePattern;
    }
  }

  else if ((strcmp(argv[1],"TimeSeries") == 0) ||
	   (strcmp(argv[1],"timeSeries") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove loadPattern patternTag?\n";
      return TCL_ERROR;
    }    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING remove loadPattern tag? failed to read tag: " << argv[2] << endln;
	return TCL_ERROR;
    }      
    bool ok =  OPS_removeTimeSeries(tag);
    if (ok == true)
      return TCL_OK;
    else
      return TCL_ERROR;
  }
  
  else if (strcmp(argv[1],"parameter") == 0) {
    if (argc < 3) {
      opserr << "WARNING want - remove parameter paramTag?\n";
      return TCL_ERROR;
    }    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove parameter tag? failed to read tag: " << argv[2] << endln;
      return TCL_ERROR;
    }      
    Parameter *theParameter = theDomain.removeParameter(tag);
    if (theParameter != 0) {
      delete theParameter;
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
    Pressure_Constraint* thePC = theDomain.removePressure_Constraint(tag);
    if(thePC != 0) {
        delete thePC;
    }
  }
  
  
  else if (strcmp(argv[1],"recorders") == 0) {
    theDomain.removeRecorders();
  }
  
  else if ((strcmp(argv[1],"recorder") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove recorder recorderTag?\n";
      return TCL_ERROR;
    }    
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove recorder tag? failed to read tag: " << argv[2] << endln;
      return TCL_ERROR;
    }      
    return theDomain.removeRecorder(tag);
  }

  else if ((strcmp(argv[1],"timeSeries") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove timeSeries $tag\n";
      return TCL_ERROR;
    }    
    
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove timeSeries tag? failed to read tag: " << argv[2] << endln;
      return TCL_ERROR;
    }      
    return OPS_removeTimeSeries(tag);
  }
  
  
  else if ((strcmp(argv[1],"SPconstraint") == 0) || (strcmp(argv[1],"sp") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove SPconstraint spTag? -or- remove SPconstraint nodeTag? dofTag? <patternTag?>\n";
      return TCL_ERROR;
    }    
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
	opserr << "WARNING remove sp tag? failed to read tag: " << argv[2] << endln;
	return TCL_ERROR;
      }      

      SP_Constraint *theSPconstraint = theDomain.removeSP_Constraint(tag);
      if (theSPconstraint != 0) {
	delete theSPconstraint;
      }
    } else {
      int nodeTag, dofTag;
      int patternTag = -1;
      
      if (Tcl_GetInt(interp, argv[2], &nodeTag) != TCL_OK) {
	opserr << "WARNING remove sp tag? failed to read node tag: " << argv[2] << endln;
	return TCL_ERROR;
      }      
      if (Tcl_GetInt(interp, argv[3], &dofTag) != TCL_OK) {
	opserr << "WARNING remove sp tag? failed to read dof tag: " << argv[3] << endln;
	return TCL_ERROR;
      }      
      
      if (argc == 5) {
	if (Tcl_GetInt(interp, argv[4], &patternTag) != TCL_OK) {
	  opserr << "WARNING remove sp tag? failed to read pattern tag: " << argv[4] << endln;
	  return TCL_ERROR;
	}      
	
      }
      dofTag--;  // one for C++ indexing of dof
      
      theDomain.removeSP_Constraint(nodeTag, dofTag, patternTag);
      
      return TCL_OK;
    }
  }

  else if ((strcmp(argv[1],"MPconstraint") == 0) || (strcmp(argv[1],"mp") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove MPconstraint nNodeTag? -or- remove MPconstraint -tag mpTag\n";
      return TCL_ERROR;
    }    
    int nodTag = 0;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &nodTag) != TCL_OK) {
	opserr << "WARNING remove mp nodeTag? failed to read nodeTag: " << argv[2] << endln;
	return TCL_ERROR;
      }      
      
      theDomain.removeMP_Constraints(nodTag);
      return TCL_OK;
    }
    if (strcmp(argv[2],"-tag") == 0 && argc > 3) {
      if (Tcl_GetInt(interp, argv[3], &nodTag) != TCL_OK) {
	opserr << "WARNING remove mp -tag mpTag? failed to read mpTag: " << argv[3] << endln;
	return TCL_ERROR;
      }      
      
      theDomain.removeMP_Constraint(nodTag);
      return TCL_OK;
    }
  }

#ifdef _RELIABILITY
// AddingSensitivity:BEGIN ///////////////////////////////////////
    else if (strcmp(argv[1],"randomVariable") == 0) {
		int rvTag;
		if (Tcl_GetInt(interp, argv[2], &rvTag) != TCL_OK) {
			opserr << "WARNING invalid input: rvTag \n";
			return TCL_ERROR;
		}
		ReliabilityDomain *theReliabilityDomain = theReliabilityBuilder->getReliabilityDomain();
		theReliabilityDomain->removeRandomVariable(rvTag);
	}
    else if (strcmp(argv[1],"performanceFunction") == 0) {
		int lsfTag;
		if (Tcl_GetInt(interp, argv[2], &lsfTag) != TCL_OK) {
			opserr << "WARNING invalid input: lsfTag \n";
			return TCL_ERROR;
		}
		ReliabilityDomain *theReliabilityDomain = theReliabilityBuilder->getReliabilityDomain();
		theReliabilityDomain->removeLimitStateFunction(lsfTag);
	}
	else if (strcmp(argv[1],"cutset") == 0) {
		int cutTag;
		if (Tcl_GetInt(interp, argv[2], &cutTag) != TCL_OK) {
			opserr << "WARNING invalid input: cutTag \n";
			return TCL_ERROR;
		}
		ReliabilityDomain *theReliabilityDomain = theReliabilityBuilder->getReliabilityDomain();
		theReliabilityDomain->removeCutset(cutTag);
	}
    else if (strcmp(argv[1],"sensitivityAlgorithm") == 0) {
		if (theSensitivityAlgorithm != 0) {
		    //theStaticAnalysis->setSensitivityAlgorithm(0);
			theSensitivityAlgorithm = 0;
			theSensitivityIntegrator = 0;
		}
	}
// AddingSensitivity:END ///////////////////////////////////////
#endif

    else
      opserr << "WARNING remove " << argv[1] << " not supported" << endln;

    return TCL_OK;
}


int
getCTestNorms(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theTest != 0) {
    const Vector &data = theTest->getNorms();
      
    char buffer [40];
    int size = data.Size();
    for (int i=0; i<size; i++) {
      sprintf(buffer,"%35.20e",data(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
	
    return TCL_OK;
  } 

  opserr << "ERROR testNorms - no convergence test!\n";
  return TCL_ERROR;
}

int
getCTestIter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (theTest != 0) {
    int res = theTest->getNumTests();
  
    char buffer [10];
    sprintf(buffer,"%d",res);
    Tcl_AppendResult(interp, buffer, NULL);

    return TCL_OK;
  }

  opserr << "ERROR testIter - no convergence test!\n";
  return TCL_ERROR;
}

int 
nodeDisp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - nodeDisp nodeTag? <dof?>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    } 

    if (argc > 2) {
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeDisp nodeTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
      }   
    }     
    
    dof--;

    const Vector *nodalResponse = theDomain.getNodeResponse(tag, Disp);

    if (nodalResponse == 0)
      return TCL_ERROR;

    int size = nodalResponse->Size();

    if (dof >= 0) {

      if (dof >= size) {
	opserr << "WARNING nodeDisp nodeTag? dof? - dofTag? too large\n";
	return TCL_ERROR;
      }
      
      double value = (*nodalResponse)(dof);
      
      // now we copy the value to the tcl string that is returned

      char buffer [40];
      sprintf(buffer,"%35.20f", value);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);
      //  sprintf(interp->result,"%35.20f ",value);
    } else {
      char buffer [40];
      for (int i=0; i<size; i++) {
	sprintf(buffer,"%35.20f",(*nodalResponse)(i));
	Tcl_AppendResult(interp, buffer, NULL);
      }
    }	
	
    return TCL_OK;
}

int 
nodeReaction(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - nodeReaction nodeTag? <dof?>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeReaction nodeTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    } 

    if (argc > 2) {
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeReaction nodeTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
      }   
    }     
    
    dof--;

    const Vector *nodalResponse = theDomain.getNodeResponse(tag, Reaction);

    if (nodalResponse == 0)
      return TCL_ERROR;

    int size = nodalResponse->Size();

    if (dof >= 0) {

      if (dof >= size) {
	opserr << "WARNING nodeReaction nodeTag? dof? - dofTag? too large\n";
	return TCL_ERROR;
      }
      
      double value = (*nodalResponse)(dof);
      
      // now we copy the value to the tcl string that is returned

      char buffer [40];
      sprintf(buffer,"%35.20f", value);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);
      //      sprintf(interp->result,"%35.20f ",value);
    } else {
      char buffer [40];
      for (int i=0; i<size; i++) {
	sprintf(buffer,"%35.20f",(*nodalResponse)(i));
	Tcl_AppendResult(interp, buffer, NULL);
      }
    }	
	
    return TCL_OK;
}

int 
nodeUnbalance(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - nodeUnbalance nodeTag? <dof?>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeUnbalance nodeTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    } 

    if (argc > 2) {
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeUnbalance nodeTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
      }   
    }     
    
    dof--;

    const Vector *nodalResponse = theDomain.getNodeResponse(tag, Unbalance);

    if (nodalResponse == 0)
      return TCL_ERROR;

    int size = nodalResponse->Size();

    if (dof >= 0) {

      if (dof >= size) {
	opserr << "WARNING nodeUnbalance nodeTag? dof? - dofTag? too large\n";
	return TCL_ERROR;
      }
      
      double value = (*nodalResponse)(dof);
      
      // now we copy the value to the tcl string that is returned
      //      sprintf(interp->result,"%35.20f ",value);

      char buffer [40];
      sprintf(buffer,"%35.20f", value);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    } else {
      char buffer [40];
      for (int i=0; i<size; i++) {
	sprintf(buffer,"%35.20f",(*nodalResponse)(i));
	Tcl_AppendResult(interp, buffer, NULL);
      }
    }	
	
    return TCL_OK;
}

int 
nodeEigenvector(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 3) {
	opserr << "WARNING want - nodeEigenVector nodeTag? eigenVector? <dof?>\n";
	return TCL_ERROR;
   }    

    int tag;
    int eigenvector = 0;
    int dof = -1;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeEigenvector nodeTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    } 

    if (Tcl_GetInt(interp, argv[2], &eigenvector) != TCL_OK) {
	opserr << "WARNING nodeEigenvector nodeTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
    }        

    if (argc > 3) {
      if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK) {
	opserr << "WARNING nodeEigenvector nodeTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
      }   
    }     
    
    dof--; eigenvector--;
    Node *theNode = theDomain.getNode(tag);
    const Matrix &theEigenvectors = theNode->getEigenvectors();

    int size = theEigenvectors.noRows();
    int numEigen = theEigenvectors.noCols();

    if (eigenvector < 0 || eigenvector >= numEigen) {
	opserr << "WARNING nodeEigenvector nodeTag? dof? - eigenvecor too large\n";
	return TCL_ERROR;
    }

    if (dof >= 0) {
      if (dof >= size) {
	opserr << "WARNING nodeEigenvector nodeTag? dof? - dofTag? too large\n";
	return TCL_ERROR;
      }

      double value = theEigenvectors(dof, eigenvector);      
      // now we copy the value to the tcl string that is returned
      //      sprintf(interp->result,"%35.20f ",value);
      char buffer [40];
      sprintf(buffer,"%35.20f", value);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    } else {

      char buffer [40];
      for (int i=0; i<size; i++) {
	double value = theEigenvectors(i, eigenvector);      
	sprintf(buffer,"%35.20f", value);
	Tcl_AppendResult(interp, buffer, NULL);
      }
    }	

    return TCL_OK;
}



int 
eleForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - eleForce eleTag? <dof?>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    

    if (argc > 2) {
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
      }   
    }     
    
    dof--;

    /*
    Element *theEle = theDomain.getElement(tag);
    if (theEle == 0)
      return TCL_ERROR;
    
    const Vector &force = theEle->getResistingForce();
    */

    const char *myArgv[1];
    char myArgv0[8]; 
    strcpy(myArgv0,"forces");
    myArgv[0] = myArgv0;

    const Vector *force = theDomain.getElementResponse(tag, &myArgv[0], 1);
    if (force != 0) {
      int size = force->Size();
      
      if (dof >= 0) {
	
	if (size < dof)
	  return TCL_ERROR;
	
	double value = (*force)(dof);
	
	// now we copy the value to the tcl string that is returned
	//	sprintf(interp->result,"%35.20f",value);

	char buffer [40];
	sprintf(buffer,"%35.20f", value);
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);	

      } else {
	char buffer[40];
	for (int i=0; i<size; i++) {
	  sprintf(buffer,"%35.20f",(*force)(i));
	  Tcl_AppendResult(interp, buffer, NULL);
	}
      }
    }

    return TCL_OK;
}

int 
localForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - localForce eleTag? <dof?>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING localForce eleTag? dof? - could not read eleTag? \n";
	return TCL_ERROR;	        
    }    

    if (argc > 2) {
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING localForce eleTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
      }   
    }     
    
    dof--;

    /*
    Element *theEle = theDomain.getElement(tag);
    if (theEle == 0)
      return TCL_ERROR;
    
    const Vector &force = theEle->getResistingForce();
    */

    const char *myArgv[1];
    char myArgv0[80]; 
    strcpy(myArgv0,"localForces");
    myArgv[0] = myArgv0;

    const Vector *force = theDomain.getElementResponse(tag, &myArgv[0], 1);
    if (force != 0) {
      int size = force->Size();
      
      if (dof >= 0) {
	
	if (size < dof)
	  return TCL_ERROR;
	
	double value = (*force)(dof);
	
	// now we copy the value to the tcl string that is returned
	//	sprintf(interp->result,"%35.20f",value);

	char buffer [40];
	sprintf(buffer,"%35.20f", value);
	Tcl_SetResult(interp, buffer, TCL_VOLATILE);
	
      } else {
	char buffer[40];
	for (int i=0; i<size; i++) {
	  sprintf(buffer,"%35.20f",(*force)(i));
	  Tcl_AppendResult(interp, buffer, NULL);
	}
      }
    }

    return TCL_OK;
}

int 
eleDynamicalForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - eleForce eleTag? <dof?>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    

    if (argc > 2) {
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
      }   
    }     
    
    dof--;
    Element *theEle = theDomain.getElement(tag);
    if (theEle == 0)
      return TCL_ERROR;
    
    const Vector &force = theEle->getResistingForceIncInertia();
    int size = force.Size();

    if (dof >= 0) {

      if (size < dof)
	return TCL_ERROR;

      double value = force(dof);
      
      // now we copy the value to the tcl string that is returned
      //      sprintf(interp->result,"%35.20f",value);
      char buffer [40];
      sprintf(buffer,"%35.20f", value);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    } else {
      char buffer[40];
      for (int i=0; i<size; i++) {
	sprintf(buffer,"%35.20f",force(i));
	Tcl_AppendResult(interp, buffer, NULL);
      }
    }

    return TCL_OK;
}



int 
eleResponse(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - eleResponse eleTag? eleArgs...\n";
	return TCL_ERROR;
   }    

    int tag;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    

    /*
    Element *theEle = theDomain.getElement(tag);
    if (theEle == 0)
      return TCL_ERROR;
    
    DummyStream dummy;
    Response *theResponse = theEle->setResponse(argv+2, argc-2, dummy);
    if (theResponse == 0) {
      return TCL_ERROR;	  
    }

    if (theResponse->getResponse() < 0) {
      delete theResponse;
      return TCL_ERROR;
    }

    Information &eleInfo = theResponse->getInformation();
    const Vector &data = eleInfo.getData();
    */

    const Vector *data = theDomain.getElementResponse(tag, argv+2, argc-2);
    if (data != 0) {
      int size = data->Size();
      char buffer[40];
      for (int i=0; i<size; i++) {
	sprintf(buffer,"%f ",(*data)(i));
	Tcl_AppendResult(interp, buffer, NULL);
      }
    }

    return TCL_OK;
}



int 
findID(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - findNodesWithID ?id\n";
	return TCL_ERROR;
   }    

    int tag;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    

    NodeIter &theNodes = theDomain.getNodes();
    Node *theNode;
    char buffer[20] ={0};
    

    while ((theNode = theNodes()) != 0) {
      DOF_Group *theGroup = theNode->getDOF_GroupPtr();
      if (theGroup != 0) {
	const ID &nodeID =  theGroup->getID();
	for (int i=0; i<nodeID.Size(); i++) {
	  if (nodeID(i) == tag) {
	    sprintf(buffer, "%d ", theNode->getTag());
	    Tcl_AppendResult(interp, buffer, NULL);
	    break;
	  }
	}
      }
    }
      
    return TCL_OK;
}

int 
nodeCoord(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - nodeCoord nodeTag? <dim?>\n";
    return TCL_ERROR;
  }    
  
  int tag;
  
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeCoord nodeTag? dim? - could not read nodeTag? \n";
    return TCL_ERROR;	        
  }    

  int dim = -1;

  if (argc > 2) {
    if (strcmp(argv[2],"X") == 0 || strcmp(argv[2],"x") == 0 ||
	strcmp(argv[2],"1") == 0)
      dim = 0;
    else if (strcmp(argv[2],"Y") == 0 || strcmp(argv[2],"y") == 0 ||
	     strcmp(argv[2],"2") == 0)
      dim = 1;
    else if (strcmp(argv[2],"Z") == 0 || strcmp(argv[2],"z") == 0 ||
	     strcmp(argv[2],"3") == 0)
      dim = 2;
    else {
      opserr << "WARNING nodeCoord nodeTag? dim? - could not read dim? \n";
      return TCL_ERROR;	        
    }        
  }    


  Node *theNode = theDomain.getNode(tag);

  if (theNode == 0) {
    return TCL_ERROR;
  }

  const Vector &coords = theNode->getCrds();
  
  int size = coords.Size();
  if (dim == -1) {
    char buffer[40];
    for (int i=0; i<size; i++) {
      sprintf(buffer,"%35.20f",coords(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
    return TCL_OK;
  }
  else if (dim < size) {
    double value = coords(dim); // -1 for OpenSees vs C indexing
    //    sprintf(interp->result,"%35.20f",value);
    char buffer [40];
    sprintf(buffer,"%35.20f", value);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
  }

  return TCL_ERROR;
}


int
fixedNodes(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    SP_Constraint* theSP;
    SP_ConstraintIter& spIter = theDomain.getDomainAndLoadPatternSPs();

    // get unique constrained nodes with set
    set<int> tags;
    int tag;
    while ((theSP = spIter()) != 0) {
        tag = theSP->getNodeTag();
        tags.insert(tag);
    }
    // assign set to vector and sort
    vector<int> tagv;
    tagv.assign(tags.begin(), tags.end());
    sort(tagv.begin(), tagv.end());
    // loop through unique, sorted tags, adding to output
    char buffer[20];
    for (int tag : tagv) {
        sprintf(buffer, "%d ", tag);
        Tcl_AppendResult(interp, buffer, NULL);
    }

    return TCL_OK;
}

int
fixedDOFs(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    if (argc < 2) {
        opserr << "WARNING want - fixedDOFs fNode?\n";
        return TCL_ERROR;
    }

    int fNode;
    if (Tcl_GetInt(interp, argv[1], &fNode) != TCL_OK) {
        opserr << "WARNING fixedDOFs fNode? - could not read fNode? \n";
        return TCL_ERROR;
    }
    
    SP_Constraint* theSP;
    SP_ConstraintIter& spIter = theDomain.getDomainAndLoadPatternSPs();

    int tag;
    Vector fixed(6);
    while ((theSP = spIter()) != 0) {
        tag = theSP->getNodeTag();
        if (tag == fNode) {
            fixed(theSP->getDOF_Number()) = 1;
        }
    }

    char buffer[20];
    for (int i = 0; i < 6; i++) {
        if (fixed(i) == 1) {
            sprintf(buffer, "%d ", i + 1);
            Tcl_AppendResult(interp, buffer, NULL);
        }
    }

    return TCL_OK;
}

int
constrainedNodes(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    bool all = 1;
    int rNode;
    if (argc > 1) {
        if (Tcl_GetInt(interp, argv[1], &rNode) != TCL_OK) {
            opserr << "WARNING constrainedNodes <rNode?> - could not read rNode? \n";
            return TCL_ERROR;
        }
        all = 0;
    }

    MP_Constraint* theMP;
    MP_ConstraintIter& mpIter = theDomain.getMPs();

    // get unique constrained nodes with set
    set<int> tags;
    int tag;
    while ((theMP = mpIter()) != 0) {
        tag = theMP->getNodeConstrained();
        if (all || rNode == theMP->getNodeRetained()) {
            tags.insert(tag);
        }
    }
    // assign set to vector and sort
    vector<int> tagv;
    tagv.assign(tags.begin(), tags.end());
    sort(tagv.begin(), tagv.end());
    // loop through unique, sorted tags, adding to output
    char buffer[20];
    for (int tag : tagv) {
        sprintf(buffer, "%d ", tag);
        Tcl_AppendResult(interp, buffer, NULL);
    }

    return TCL_OK;
}

int
constrainedDOFs(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    if (argc < 2) {
        opserr << "WARNING want - constrainedDOFs cNode? <rNode?> <rDOF?>\n";
        return TCL_ERROR;
    }

    int cNode;
    if (Tcl_GetInt(interp, argv[1], &cNode) != TCL_OK) {
        opserr << "WARNING constrainedDOFs cNode? <rNode?> <rDOF?> - could not read cNode? \n";
        return TCL_ERROR;
    }

    int rNode;
    bool allNodes = 1;
    if (argc > 2) {
        if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
            opserr << "WARNING constrainedDOFs cNode? <rNode?> <rDOF?> - could not read rNode? \n";
            return TCL_ERROR;
        }
        allNodes = 0;
    }

    int rDOF;
    bool allDOFs = 1;
    if (argc > 3) {
        if (Tcl_GetInt(interp, argv[3], &rDOF) != TCL_OK) {
            opserr << "WARNING constrainedDOFs cNode? <rNode?> <rDOF?> - could not read rDOF? \n";
            return TCL_ERROR;
        }
        rDOF--;
        allDOFs = 0;
    }

    MP_Constraint* theMP;
    MP_ConstraintIter& mpIter = theDomain.getMPs();

    int tag;
    int i;
    int n;
    Vector constrained(6);
    while ((theMP = mpIter()) != 0) {
        tag = theMP->getNodeConstrained();
        if (tag == cNode) {
            if (allNodes || rNode == theMP->getNodeRetained()) {
                const ID &cDOFs = theMP->getConstrainedDOFs();
                n = cDOFs.Size();
                if (allDOFs) {
                    for (i = 0; i < n; i++) {
                        constrained(cDOFs(i)) = 1;
                    }
                }
                else {
                    const ID &rDOFs = theMP->getRetainedDOFs();
                    for (i = 0; i < n; i++) {
                        if (rDOF == rDOFs(i))
                            constrained(cDOFs(i)) = 1;
                    }
                }
            }
        }
    }
    char buffer[20];
    for (int i = 0; i < 6; i++) {
        if (constrained(i) == 1) {
            sprintf(buffer, "%d ", i + 1);
            Tcl_AppendResult(interp, buffer, NULL);
        }
    }

    return TCL_OK;
}

int
retainedNodes(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    bool all = 1;
    int cNode;
    if (argc > 1) {
        if (Tcl_GetInt(interp, argv[1], &cNode) != TCL_OK) {
            opserr << "WARNING retainedNodes <cNode?> - could not read cNode? \n";
            return TCL_ERROR;
        }
        all = 0;
    }

    MP_Constraint* theMP;
    MP_ConstraintIter& mpIter = theDomain.getMPs();

    // get unique constrained nodes with set
    set<int> tags;
    int tag;
    while ((theMP = mpIter()) != 0) {
        tag = theMP->getNodeRetained();
        if (all || cNode == theMP->getNodeConstrained()) {
            tags.insert(tag);
        }
    }
    // assign set to vector and sort
    vector<int> tagv;
    tagv.assign(tags.begin(), tags.end());
    sort(tagv.begin(), tagv.end());
    // loop through unique, sorted tags, adding to output
    char buffer[20];
    for (int tag : tagv) {
        sprintf(buffer, "%d ", tag);
        Tcl_AppendResult(interp, buffer, NULL);
    }

    return TCL_OK;
}

int
retainedDOFs(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{

    if (argc < 2) {
        opserr << "WARNING want - retainedDOFs rNode? <cNode?> <cDOF?>\n";
        return TCL_ERROR;
    }

    int rNode;
    if (Tcl_GetInt(interp, argv[1], &rNode) != TCL_OK) {
        opserr << "WARNING retainedDOFs rNode? <cNode?> <cDOF?> - could not read rNode? \n";
        return TCL_ERROR;
    }

    int cNode;
    bool allNodes = 1;
    if (argc > 2) {
        if (Tcl_GetInt(interp, argv[2], &cNode) != TCL_OK) {
            opserr << "WARNING retainedDOFs rNode? <cNode?> <cDOF?> - could not read cNode? \n";
            return TCL_ERROR;
        }
        allNodes = 0;
    }

    int cDOF;
    bool allDOFs = 1;
    if (argc > 3) {
        if (Tcl_GetInt(interp, argv[3], &cDOF) != TCL_OK) {
            opserr << "WARNING retainedDOFs rNode? <cNode?> <cDOF?> - could not read cDOF? \n";
            return TCL_ERROR;
        }
        cDOF--;
        allDOFs = 0;
    }

    MP_Constraint* theMP;
    MP_ConstraintIter& mpIter = theDomain.getMPs();

    int tag;
    int i;
    int n;
    Vector retained(6);
    while ((theMP = mpIter()) != 0) {
        tag = theMP->getNodeRetained();
        if (tag == rNode) {
            if (allNodes || cNode == theMP->getNodeConstrained()) {
                const ID& rDOFs = theMP->getRetainedDOFs();
                n = rDOFs.Size();
                if (allDOFs) {
                    for (i = 0; i < n; i++) {
                        retained(rDOFs(i)) = 1;
                    }
                }
                else {
                    const ID& cDOFs = theMP->getConstrainedDOFs();
                    for (i = 0; i < n; i++) {
                        if (cDOF == cDOFs(i))
                            retained(rDOFs(i)) = 1;
                    }
                }
            }
        }
    }
    char buffer[20];
    for (int i = 0; i < 6; i++) {
        if (retained(i) == 1) {
            sprintf(buffer, "%d ", i + 1);
            Tcl_AppendResult(interp, buffer, NULL);
        }
    }

    return TCL_OK;
}

int 
setNodeCoord(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - setNodeCoord nodeTag? dim? value?\n";
    return TCL_ERROR;
  }    
  
  int tag;
  
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read nodeTag? \n";
    return TCL_ERROR;	        
  }    

	int dim; 
	double value;

	if (Tcl_GetInt(interp, argv[2], &dim) != TCL_OK) {
	    opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read dim? \n";
		return TCL_ERROR;
	}
	if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
	    opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read value? \n";
		return TCL_ERROR;
	}

  Node *theNode = theDomain.getNode(tag);

  if (theNode == 0) {
    return TCL_ERROR;
  }

  Vector coords(theNode->getCrds());
  coords(dim-1) = value;
  theNode->setCrds(coords);
  
  return TCL_OK;
}

int 
updateElementDomain(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // Need to "setDomain" to make the change take effect. 
    ElementIter &theElements = theDomain.getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      theElement->setDomain(&theDomain);
    }

	return 0;
}

int
getNDM(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    int ndm;

    if (argc > 1) {
        int tag;
        if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
            opserr << "WARNING ndm nodeTag? \n";
            return TCL_ERROR;
        }
        Node* theNode = theDomain.getNode(tag);
        if (theNode == 0) {
            opserr << "WARNING nodeTag " << tag << " does not exist \n";
            return TCL_ERROR;
        }
        const Vector& coords = theNode->getCrds();
        ndm = coords.Size();
    } else {
        if (theBuilder == 0) {
            return TCL_OK;
        }
        else {
            ndm = OPS_GetNDM();
        }
    }

    char buffer[20];
    sprintf(buffer, "%d", ndm);
    Tcl_AppendResult(interp, buffer, NULL);

    return TCL_OK;
}

int
getNDF(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    int ndf;

    if (argc > 1) {
        int tag;
        if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
            opserr << "WARNING ndf nodeTag? \n";
            return TCL_ERROR;
        }
        Node* theNode = theDomain.getNode(tag);
        if (theNode == 0) {
            opserr << "WARNING nodeTag " << tag << " does not exist \n";
            return TCL_ERROR;
        }
        ndf = theNode->getNumberDOF();
    } else {
        if (theBuilder == 0) {
            return TCL_OK;
        }
        else {
            ndf = OPS_GetNDF();
        }
    }

    char buffer[20];
    sprintf(buffer, "%d", ndf);
    Tcl_AppendResult(interp, buffer, NULL);

    return TCL_OK;
}

int
eleType(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
    if (argc < 2) {
        opserr << "WARNING want - eleType eleTag?\n";
        return TCL_ERROR;
    }

    int tag;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
        opserr << "WARNING eleType eleTag? \n";
        return TCL_ERROR;
    }

    char buffer[20];
    Element* theElement = theDomain.getElement(tag);
    if (theElement == 0) {
        opserr << "WARNING eleType ele " << tag << " not found" << endln;
        return TCL_ERROR;
    }
    const char* type = theElement->getClassType();
    sprintf(buffer, "%s", type);
    Tcl_AppendResult(interp, buffer, NULL);

    return TCL_OK;
}

int 
eleNodes(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "WARNING want - eleNodes eleTag?\n";
    return TCL_ERROR;
  }    
  
  int tag;
  
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING eleNodes eleTag? \n";
    return TCL_ERROR;	        
  }    
  
  char buffer[20];

  const char *myArgv[1];
  char myArgv0[80]; 
  strcpy(myArgv0,"nodeTags");
  myArgv[0] = myArgv0;

  //const Vector *tags = theDomain.getElementResponse(tag, &myArgv[0], 1);
  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
      opserr << "WARNING eleNodes ele " << tag << " not found" << endln;
      return TCL_ERROR;
  }
  int numTags = theElement->getNumExternalNodes();
  const ID &tags = theElement->getExternalNodes();
  for (int i = 0; i < numTags; i++) {
    sprintf(buffer, "%d ", tags(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }
  
  return TCL_OK;
}

int 
nodeDOFs(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "WARNING want - nodeDOFs nodeTag?\n";
    return TCL_ERROR;
  }    
  
  int tag;
  
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
    return TCL_ERROR;	        
  }    
  
  char buffer[40];

  Node *theNode = theDomain.getNode(tag);
  if (theNode == 0) {
    opserr << "WARNING nodeDOFs node " << tag << " not found" << endln;
    return TCL_ERROR;
  }
  int numDOF = theNode->getNumberDOF();

  DOF_Group *theDOFgroup = theNode->getDOF_GroupPtr();
  if (theDOFgroup == 0) {
    opserr << "WARNING nodeDOFs DOF group null" << endln;
    return -1;
  }
  const ID &eqnNumbers = theDOFgroup->getID();
  for (int i = 0; i < numDOF; i++) {
    sprintf(buffer, "%d ", eqnNumbers(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }
  
  return TCL_OK;
}

int 
nodeMass(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 3) {
    opserr << "WARNING want - nodeMass nodeTag? nodeDOF?\n";
    return TCL_ERROR;
  }    
  
  int tag, dof;
  
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
    return TCL_ERROR;	        
  }    
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
    return TCL_ERROR;	        
  }    
  
  char buffer[40];

  Node *theNode = theDomain.getNode(tag);
  if (theNode == 0) {
    opserr << "WARNING nodeMass node " << tag << " not found" << endln;
    return TCL_ERROR;
  }
  int numDOF = theNode->getNumberDOF();
  if (dof < 1 || dof > numDOF) {
    opserr << "WARNING nodeMass dof " << dof << " not in range" << endln;
    return TCL_ERROR;
  }
  else {
    const Matrix &mass = theNode->getMass();
    sprintf(buffer, "%35.20f", mass(dof-1,dof-1));
    Tcl_AppendResult(interp, buffer, NULL);
  }
  
  return TCL_OK;
}

int 
nodePressure(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    if(argc < 2) {
        opserr << "WARNING: want - nodePressure nodeTag?\n";
        return TCL_ERROR;
    }
    int tag;
    if(Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
        opserr << "WARNING: nodePressure "<<argv[1]<<"\n";
        return TCL_ERROR;
    }
    double pressure = 0.0;
    Pressure_Constraint* thePC = theDomain.getPressure_Constraint(tag);
    if(thePC != 0) {
        pressure = thePC->getPressure();
        // int ptag = thePC->getPressureNode();
        // Node* pNode = theDomain.getNode(ptag);
        // if(pNode != 0) {
        //     const Vector& vel = pNode->getVel();
        //     if(vel.Size() > 0)  {
        //         pressure = vel(0);
        //         // opserr<<"pressure = "<<pressure<<"\n";
        //     }
        // }
    }
    char buffer[80];
    sprintf(buffer, "%35.20f", pressure);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

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
  for (int j=0; j<6; j++) {
    cnt += sprintf(&resDataPtr[cnt], "%.6e  ", bounds(j));
  }
  
  Tcl_SetResult(interp, resDataPtr, TCL_STATIC);
  
  return TCL_OK;
}

int 
nodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - nodeVel nodeTag? <dof?>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeVel nodeTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    
    if (argc > 2) {
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeVel nodeTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
      }        
    }    

    dof--;

    const Vector *nodalResponse = theDomain.getNodeResponse(tag, Vel);

    if (nodalResponse == 0)
      return TCL_ERROR;

    int size = nodalResponse->Size();

    if (dof >= 0) {
      if (size < dof)
	return TCL_ERROR;

      double value = (*nodalResponse)(dof);
      
      // now we copy the value to the tcl string that is returned
      //      sprintf(interp->result,"%35.20f",value);
      char buffer [40];
      sprintf(buffer,"%35.20f", value);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    } else {

      char buffer[40];
      for (int i=0; i<size; i++) {
	sprintf(buffer,"%35.20f",(*nodalResponse)(i));
	Tcl_AppendResult(interp, buffer, NULL);
      }
    }
	
    return TCL_OK;
}

int 
setNodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 4) {
	opserr << "WARNING want - setNodeVel nodeTag? dof? value? <-commit>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;
    double value = 0.0;
    bool commit = false;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    

    Node *theNode = theDomain.getNode(tag);
    if (theNode == 0) {
      opserr << "WARNING setNodeVel -- node with tag " << tag << " not found" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read dof? \n";
      return TCL_ERROR;	        
    }        
    if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
      opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read value? \n";
      return TCL_ERROR;	        
    }        
    if (argc > 4 && strcmp(argv[4],"-commit") == 0)
      commit = true;

    dof--;

    int numDOF = theNode->getNumberDOF();

    if (dof >= 0 && dof < numDOF) {
      Vector vel(numDOF);
      vel = theNode->getVel();
      vel(dof) = value;
      theNode->setTrialVel(vel);
    }
    if (commit)
      theNode->commitState();
	
    return TCL_OK;
}

int 
setNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 4) {
	opserr << "WARNING want - setNodeDisp nodeTag? dof? value? <-commit>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;
    double value = 0.0;
    bool commit = false;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING setNodeDisp nodeTag? dof? value?- could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    

    Node *theNode = theDomain.getNode(tag);
    if (theNode == 0) {
      opserr << "WARNING setNodeDisp -- node with tag " << tag << " not found" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING setNodeDisp nodeTag? dof? value?- could not read dof? \n";
      return TCL_ERROR;	        
    }        
    if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
      opserr << "WARNING setNodeDisp nodeTag? dof? value?- could not read value? \n";
      return TCL_ERROR;	        
    }        
    if (argc > 4 && strcmp(argv[4],"-commit") == 0)
      commit = true;

    dof--;

    int numDOF = theNode->getNumberDOF();

    if (dof >= 0 && dof < numDOF) {
      Vector vel(numDOF);
      vel = theNode->getDisp();
      vel(dof) = value;
      theNode->setTrialDisp(vel);
    }
    if (commit)
      theNode->commitState();

    return TCL_OK;
}

int 
setNodeAccel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 4) {
	opserr << "WARNING want - setNodeAccel nodeTag? dof? value? <-commit>\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;
    double value = 0.0;
    bool commit = false;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING setNodeAccel nodeTag? dof? value?- could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    

    Node *theNode = theDomain.getNode(tag);
    if (theNode == 0) {
      opserr << "WARNING setNodeAccel -- node with tag " << tag << " not found" << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING setNodeDisp nodeTag? dof? value?- could not read dof? \n";
      return TCL_ERROR;	        
    }        
    if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
      opserr << "WARNING setNodeAccel nodeTag? dof? value?- could not read value? \n";
      return TCL_ERROR;	        
    }        
    if (argc > 4 && strcmp(argv[4],"-commit") == 0)
      commit = true;

    dof--;

    int numDOF = theNode->getNumberDOF();

    if (dof >= 0 && dof < numDOF) {
      Vector vel(numDOF);
      vel = theNode->getAccel();
      vel(dof) = value;
      theNode->setTrialAccel(vel);
    }
    if (commit)
      theNode->commitState();
	
    return TCL_OK;
}

int 
nodeAccel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 2) {
	opserr << "WARNING want - nodeAccel nodeTag? dof?\n";
	return TCL_ERROR;
   }    

    int tag;
    int dof = -1;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeAccel nodeTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    
    if (argc > 2) {
      if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeAccel nodeTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
      }   
    }     
    
    dof--;

    const Vector *nodalResponse = theDomain.getNodeResponse(tag, Accel);
    if (nodalResponse == 0)
      return TCL_ERROR;


    int size =  nodalResponse->Size();

    if (dof >= 0) {
      if (size < dof)
	return TCL_ERROR;

      double value = (*nodalResponse)(dof);
    
      // now we copy the value to the tcl string that is returned
      //sprintf(interp->result,"%35.20f",value);
      char buffer [40];
      sprintf(buffer,"%35.20f", value);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    } else {
      char buffer[40];
      for (int i=0; i<size; i++) {
	sprintf(buffer,"%35.20f",(*nodalResponse)(i));
	Tcl_AppendResult(interp, buffer, NULL);
      }
    }
	
    return TCL_OK;
}


int 
nodeResponse(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
    if (argc < 4) {
	opserr << "WARNING want - nodeResponse nodeTag? dof? responseID?\n";
	return TCL_ERROR;
   }    

    int tag, dof, responseID;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING nodeResponse nodeTag? dof? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING nodeResponse nodeTag? dof? - could not read dof? \n";
	return TCL_ERROR;	        
    }        
    if (Tcl_GetInt(interp, argv[3], &responseID) != TCL_OK) {
	opserr << "WARNING nodeResponse nodeTag? dof? responseID? - could not read responseID? \n";
	return TCL_ERROR;	        
    }        
    
    dof--;

    const Vector *nodalResponse = theDomain.getNodeResponse(tag, (NodeResponseType)responseID);
    if (nodalResponse == 0 || nodalResponse->Size() < dof || dof < 0)
      return TCL_ERROR;

    double value = (*nodalResponse)(dof);
    
    // now we copy the value to the tcl string that is returned
    //    sprintf(interp->result,"%35.20f",value);
    char buffer [40];
    sprintf(buffer,"%35.20f", value);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);	

    return TCL_OK;
}

int 
calculateNodalReactions(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    // make sure at least one other argument to contain type of system
  int incInertia = 0;

  if (argc == 2)  {
    if ((strcmp(argv[1],"-incInertia") == 0)
	|| (strcmp(argv[1],"-dynamical") == 0)
	|| (strcmp(argv[1],"-Dynamic") == 0)
	|| (strcmp(argv[1],"-dynamic") == 0))

      incInertia = 1;
    
    else if ((strcmp(argv[1],"-rayleigh") == 0))

      incInertia = 2;
  }

  theDomain.calculateNodalReactions(incInertia);

  return TCL_OK;
}


// AddingSensitivity:BEGIN ////////////////////////////////////
int 
sensNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

    // make sure at least one other argument to contain type of system
    if (argc < 4) {
      opserr << "WARNING want - sensNodeDisp nodeTag? dof? paramTag?\n";
      return TCL_ERROR;
    }    

     int tag, dof, paramTag;

     if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	 opserr << "WARNING nodeDisp nodeTag? dof? paramTag?- could not read nodeTag? ";
	 return TCL_ERROR;	        
     }    
     if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	 opserr << "WARNING nodeDisp nodeTag? dof? paramTag?- could not read dof? ";
	 return TCL_ERROR;	        
     }        
     if (Tcl_GetInt(interp, argv[3], &paramTag) != TCL_OK) {
	 opserr << "WARNING nodeDisp paramTag? dof? paramTag?- could not read paramTag? ";
	 return TCL_ERROR;	        
     }        

     Node *theNode = theDomain.getNode(tag);
     if (theNode == 0) {
       opserr << "sensNodeDisp: node " << tag << " not found" << endln;
       return TCL_ERROR;
     }

     Parameter *theParam = theDomain.getParameter(paramTag);
     if (theParam == 0) {
       opserr << "sensNodeDisp: parameter " << paramTag << " not found" << endln;
       return TCL_ERROR;
     }

     int gradIndex = theParam->getGradIndex();

     double value = theNode->getDispSensitivity(dof,gradIndex);

     char buffer[40];
     sprintf(buffer,"%35.20f",value);
     Tcl_SetResult(interp, buffer, TCL_VOLATILE);

     return TCL_OK;
 }

 int 
 sensNodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
 {

     // make sure at least one other argument to contain type of system
     if (argc < 4) {
       opserr << "WARNING want - sensNodeVel nodeTag? dof? paramTag?\n";
       return TCL_ERROR;
    }    

     int tag, dof, paramTag;

     if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	 opserr << "WARNING sensNodeVel nodeTag? dof? paramTag? - could not read nodeTag? \n";
	 return TCL_ERROR;	        
     }    
     if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	 opserr << "WARNING sensNodeVel nodeTag? dof? paramTag? - could not read dof? \n";
	 return TCL_ERROR;	        
     }        
     if (Tcl_GetInt(interp, argv[3], &paramTag) != TCL_OK) {
	 opserr << "WARNING sensNodeVel nodeTag? dof? paramTag? - could not read paramTag? \n";
	 return TCL_ERROR;	        
     }        

     Node *theNode = theDomain.getNode(tag);
     if (theNode == 0) {
       opserr << "sensNodeVel: node " << tag << " not found" << endln;
       return TCL_ERROR;
     }

     Parameter *theParam = theDomain.getParameter(paramTag);
     if (theParam == 0) {
       opserr << "sensNodeVel: parameter " << paramTag << " not found" << endln;
       return TCL_ERROR;
     }

     int gradIndex = theParam->getGradIndex();

     double value = theNode->getVelSensitivity(dof,gradIndex);

     char buffer[40];
     sprintf(buffer,"%35.20f",value);

     Tcl_SetResult(interp, buffer, TCL_VOLATILE);

     return TCL_OK;
 }

 int 
 sensNodeAccel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
 {

     // make sure at least one other argument to contain type of system
     if (argc < 4) {
       opserr <<  "WARNING want - sensNodeAccel nodeTag? dof? paramTag?\n";
       return TCL_ERROR;
   }    

    int tag, dof, paramTag;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING sensNodeAccel nodeTag? dof? paramTag? - could not read nodeTag? \n";
	return TCL_ERROR;	        
    }    
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
	opserr << "WARNING sensNodeAccel nodeTag? dof? paramTag? - could not read dof? \n";
	return TCL_ERROR;	        
    }        
    if (Tcl_GetInt(interp, argv[3], &paramTag) != TCL_OK) {
	opserr << "WARNING sendNodeAccel nodeTag? dof? paramTag? - could not read paramTag? \n";
	return TCL_ERROR;	        
    }        
    
    Node *theNode = theDomain.getNode(tag);
    if (theNode == 0) {
      opserr << "sensNodeAccel: node " << tag << " not found" << endln;
      return TCL_ERROR;
    }

    Parameter *theParam = theDomain.getParameter(paramTag);
    if (theParam == 0) {
      opserr << "sensNodeAccel: parameter " << paramTag << " not found" << endln;
      return TCL_ERROR;
    }

    int gradIndex = theParam->getGradIndex();
    
    double value = theNode->getAccSensitivity(dof,gradIndex);
    
    char buffer[40];
    sprintf(buffer,"%35.20f",value);

    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
}

int
sensNodePressure(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

    // make sure at least one other argument to contain type of system
    if (argc < 3) {
      opserr << "WARNING want - sensNodePressure nodeTag? paramTag?\n";
      return TCL_ERROR;
    }    

    int tag, paramTag;

    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
	opserr << "WARNING sensNodePressure nodeTag? paramTag?- could not read nodeTag? ";
	return TCL_ERROR;	        
    }    
    if (Tcl_GetInt(interp, argv[2], &paramTag) != TCL_OK) {
	opserr << "WARNING sensNodePressure paramTag? paramTag?- could not read paramTag? ";
	return TCL_ERROR;	        
    }        
    
    double dp = 0.0;
    Pressure_Constraint* thePC = theDomain.getPressure_Constraint(tag);
    if(thePC != 0) {
        // int ptag = thePC->getPressureNode();
        // Node* pNode = theDomain.getNode(ptag);
        Node* pNode = thePC->getPressureNode();
        if(pNode != 0) {

            Parameter *theParam = theDomain.getParameter(paramTag);
            if (theParam == 0) {
                opserr << "sensNodePressure: parameter " << paramTag << " not found" << endln;
                return TCL_ERROR;
            }

            int gradIndex = theParam->getGradIndex();
            dp = pNode->getVelSensitivity(1,gradIndex);
        }
    }
    
    char buffer[40];
    sprintf(buffer,"%35.20f",dp);

    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
}


int 
sensSectionForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
#ifdef _RELIABILITY
  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << "WARNING want - sensSectionForce eleTag? <secNum?> dof? paramTag?\n";
    return TCL_ERROR;
  }    
  
  //opserr << "sensSectionForce: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag, dof, paramTag;
  int secNum = 0;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sensSectionForce eleTag? secNum? dof? paramTag?- could not read eleTag? \n";
    return TCL_ERROR;	        
  }    

  // Make this work for zeroLengthSection too
  int currentArg = 2;
  if (argc > 4) {
    if (Tcl_GetInt(interp, argv[currentArg++], &secNum) != TCL_OK) {
      opserr << "WARNING sensSectionForce eleTag? secNum? dof? paramTag?- could not read secNum? \n";
      return TCL_ERROR;	        
    }
  }
  if (Tcl_GetInt(interp, argv[currentArg++], &dof) != TCL_OK) {
    opserr << "WARNING sensSectionForce eleTag? secNum? dof? paramTag?- could not read dof? \n";
    return TCL_ERROR;	        
  }        
  if (Tcl_GetInt(interp, argv[currentArg++], &paramTag) != TCL_OK) {
    opserr << "WARNING sensSectionForce eleTag? secNum? dof? paramTag?- could not read paramTag? \n";
    return TCL_ERROR;	        
  }

  ParameterIter &pIter = theDomain.getParameters();
  Parameter *theParam;
  while ((theParam = pIter()) != 0)
    theParam->activate(false);

  theParam = theDomain.getParameter(paramTag);
  int gradIndex = theParam->getGradIndex();
  theParam->activate(true);

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sensSectionForce element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "dsdh";
  const char *argvv[3];
  int argcc = 3;
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;
  if (argc < 5) { // For zeroLengthSection
    argcc = 2;
    argvv[1] = c;
  }

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    Tcl_SetResult(interp, "0.0", TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponseSensitivity(gradIndex);
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer,"%12.8g",theVec(dof-1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  theParam->activate(false);

  delete theResponse;
#endif
  return TCL_OK;
}

int 
sectionForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr <<  "WARNING want - sectionForce eleTag? <secNum?> dof? \n";
    return TCL_ERROR;
  }    
  
  //opserr << "sectionForce: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag, dof;
  int secNum = 0;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionForce eleTag? secNum? dof? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }

  // Make this work for zeroLengthSection too
  int currentArg = 2;
  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[currentArg++], &secNum) != TCL_OK) {
      opserr << "WARNING sectionForce eleTag? secNum? dof? - could not read secNum? \n";
      return TCL_ERROR;	        
    }
  }
  if (Tcl_GetInt(interp, argv[currentArg++], &dof) != TCL_OK) {
    opserr << "WARNING sectionForce eleTag? secNum? dof? - could not read dof? \n";
    return TCL_ERROR;	        
  }        

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionForce element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "force";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;
  if (argc < 4) { // For zeroLengthSection
    argcc = 2;
    argvv[1] = c;
  }

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer,"%12.8g",theVec(dof-1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  delete theResponse;

  return TCL_OK;
}

int 
sectionDeformation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr <<"WARNING want - sectionDeformation eleTag? secNum? dof? \n";
    return TCL_ERROR;
  }    
  
  //opserr << "sectionDeformation: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag, secNum, dof;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionDeformation eleTag? secNum? dof? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }    
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING sectionDeformation eleTag? secNum? dof? - could not read secNum? \n";
    return TCL_ERROR;	        
  }    
  if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK) {
    opserr << "WARNING sectionDeformation eleTag? secNum? dof? - could not read dof? \n";
    return TCL_ERROR;	        
  }        

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionDeformation element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "deformation";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer,"%12.8g",theVec(dof-1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  delete theResponse;

  return TCL_OK;
}


int 
sectionLocation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sectionLocation eleTag? secNum? \n";
    return TCL_ERROR;
  }    
  
  //opserr << "sectionDeformation: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionLocation eleTag? secNum? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }    
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING sectionLocation eleTag? secNum? - could not read secNum? \n";
    return TCL_ERROR;	        
  }    

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionLocation element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  int argcc = 1;
  char a[80] = "integrationPoints";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer [] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer,"%12.8g",theVec(secNum-1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  delete theResponse;

  return TCL_OK;
}

int 
sectionWeight(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sectionWeight eleTag? secNum? \n";
    return TCL_ERROR;
  }    
  
  //opserr << "sectionDeformation: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionWeight eleTag? secNum? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }    
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING sectionWeight eleTag? secNum? - could not read secNum? \n";
    return TCL_ERROR;	        
  }    

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionWeight element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  int argcc = 1;
  char a[80] = "integrationWeights";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  char buffer[40];
  sprintf(buffer,"%12.8g",theVec(secNum-1));

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  delete theResponse;

  return TCL_OK;
}

int 
sectionStiffness(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sectionStiffness eleTag? secNum? \n";
    return TCL_ERROR;
  }    
  
  //opserr << "sectionDeformation: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionStiffness eleTag? secNum? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }    
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING sectionStiffness eleTag? secNum? - could not read secNum? \n";
    return TCL_ERROR;	        
  }    

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionStiffness element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "stiffness";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer [] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Matrix &theMat = *(info.theMatrix);
  int nsdof = theMat.noCols();

  char buffer[200];
  for (int i = 0; i < nsdof; i++) {
    for (int j = 0; j < nsdof; j++) {
      sprintf(buffer,"%12.8g ",theMat(i,j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }
  
  delete theResponse;

  return TCL_OK;
}

int 
sectionFlexibility(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << "WARNING want - sectionFlexibility eleTag? secNum? \n";
    return TCL_ERROR;
  }    
  
  //opserr << "sectionDeformation: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING sectionFlexibility eleTag? secNum? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }    
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING sectionFlexibility eleTag? secNum? - could not read secNum? \n";
    return TCL_ERROR;	        
  }    

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING sectionFlexibility element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "flexibility";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Matrix &theMat = *(info.theMatrix);
  int nsdof = theMat.noCols();

  char buffer[200];
  for (int i = 0; i < nsdof; i++) {
    for (int j = 0; j < nsdof; j++) {
      sprintf(buffer,"%12.8g ",theMat(i,j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }
  
  delete theResponse;

  return TCL_OK;
}


int 
basicDeformation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - basicDeformation eleTag? \n";
    return TCL_ERROR;
  }    
  
  //opserr << "sectionDeformation: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING basicDeformation eleTag? dofNum? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }
  /*    
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING basicDeformation eleTag? dofNum? - could not read dofNum? \n";
    return TCL_ERROR;	        
  }    
  */

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING basicDeformation element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  int argcc = 1;
  char a[80] = "basicDeformation";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer [] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);
  int nbf = theVec.Size();

  char buffer[200];
  for (int i = 0; i < nbf; i++) {
    sprintf(buffer, "%12.8f ", theVec(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  delete theResponse;

  return TCL_OK;
}

int 
basicForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - basicForce eleTag? \n";
    return TCL_ERROR;
  }    
  
  //opserr << "sectionDeformation: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING basicForce eleTag? dofNum? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }
  /*    
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING basicDeformation eleTag? dofNum? - could not read dofNum? \n";
    return TCL_ERROR;	        
  }    
  */

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING basicDeformation element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  int argcc = 1;
  char a[80] = "basicForce";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);
  int nbf = theVec.Size();

  char buffer[200];
  for (int i = 0; i < nbf; i++) {
    sprintf(buffer, "%12.8f ", theVec(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  delete theResponse;

  return TCL_OK;
}

int 
basicStiffness(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING want - basicStiffness eleTag? \n";
    return TCL_ERROR;
  }    
  
  //opserr << "sectionDeformation: ";
  //for (int i = 0; i < argc; i++) 
  //  opserr << argv[i] << ' ' ;
  //opserr << endln;

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING basicStiffness eleTag? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }
  /*    
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << "WARNING basicDeformation eleTag? dofNum? - could not read dofNum? \n";
    return TCL_ERROR;	        
  }    
  */

  Element *theElement = theDomain.getElement(tag);
  if (theElement == 0) {
    opserr << "WARNING basicStiffness element with tag " << tag << " not found in domain \n";
    return TCL_ERROR; 
  }

  int argcc = 1;
  char a[80] = "basicStiffness";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == 0) {
    char buffer[] = "0.0";
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Matrix &theMatrix = *(info.theMatrix);
  int nbf = theMatrix.noCols();

  char buffer[200];
  for (int i = 0; i < nbf; i++) {
    for (int j = 0; j < nbf; j++) {
      sprintf(buffer, "%12.8f ", theMatrix(i,j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  delete theResponse;

  return TCL_OK;
}

// added by C.McGann, U.Washington
int
InitialStateAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "WARNING: Incorrect number of arguments for InitialStateAnalysis command" << endln;
    return TCL_ERROR;
  }
  
  if (strcmp(argv[1],"on") == 0) {
    opserr << "InitialStateAnalysis ON" << endln;

    // set global variable to true
    // FMK changes for parallel: 
    // ops_InitialStateAnalysis = true;

    Parameter *theP = new InitialStateParameter(true);
    theDomain.addParameter(theP);
    delete theP;

    return TCL_OK;
    
  } else if (strcmp(argv[1],"off") == 0) {
    opserr << "InitialStateAnalysis OFF" <<endln;
    
    // call revert to start to zero the displacements
    theDomain.revertToStart();
    
    // set global variable to false
    // FMK changes for parallel
    // ops_InitialStateAnalysis = false;
    Parameter *theP = new InitialStateParameter(false);
    theDomain.addParameter(theP);
    delete theP;

    return TCL_OK;
    
  } else {
    opserr << "WARNING: Incorrect arguments - want InitialStateAnalysis on, or InitialStateAnalysis off" << endln;
    
    return TCL_ERROR;		
  }
}

int 
computeGradients(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
#ifdef _RELIABILITY
  if (theSensitivityAlgorithm != 0)
    theSensitivityAlgorithm->computeSensitivities();
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



int 
modalDamping(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) { 
    opserr << "WARNING modalDamping ?factor - not enough arguments to command\n";
    return TCL_ERROR;
  }

  if (numEigen == 0 || theEigenSOE == 0) {
    opserr << "WARNING - modalDmping - eigen command needs to be called first - NO MODAL DAMPING APPLIED\n ";
  }

  int numModes = argc - 1;
  double factor;
  Vector modalDampingValues(numEigen);

  if (numModes != 1 && numModes != numEigen) {
    opserr << "WARNING modalDmping - same #damping factors as modes must be specified\n";
    opserr << "                    - same damping ratio will be applied to all\n";
  }

  // 
  // read in values and set factors
  //

  if (numModes == numEigen) {

    for (int i=0; i<numEigen; i++) {
      if (Tcl_GetDouble(interp, argv[1+i], &factor) != TCL_OK) {
	opserr << "WARNING modalDamping - could not read factor for model " << i+1 << endln;
	return TCL_ERROR;	        
      }        
      modalDampingValues[i] = factor;    
    } 

  } else {

    if (Tcl_GetDouble(interp, argv[1], &factor) != TCL_OK) {
      opserr << "WARNING modalDamping - could not read factor for all modes \n";
      return TCL_ERROR;	        
    }        

    for (int i=0; i<numEigen; i++)
      modalDampingValues[i] = factor;
  }

  // set factors in domain
  theDomain.setModalDampingFactors(&modalDampingValues, true);

  //opserr << "modalDamping Factors: " << modalDampingValues;

  return TCL_OK;
}

int 
modalDampingQ(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) { 
    opserr << "WARNING modalDamping ?factor - not enough arguments to command\n";
    return TCL_ERROR;
  }

  if (numEigen == 0 || theEigenSOE == 0) {
    opserr << "WARINING - modalDmping - eigen command needs to be called first - NO MODAL DAMPING APPLIED\n ";
  }


  int numModes = argc - 1;
  double factor = 0;
  Vector modalDampingValues(numEigen);

  if (numModes != 1 && numModes != numEigen) {
    opserr << "WARNING modalDmping - same #damping factors as modes must be specified\n";
    opserr << "                    - same damping ratio will be applied to all";
  }

  // 
  // read in values and set factors
  //

  if (numModes == numEigen) {

    // read in all factors one at a time
    for (int i=0; i<numEigen; i++) {
      if (Tcl_GetDouble(interp, argv[1+i], &factor) != TCL_OK) {
	opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaK? \n";
	return TCL_ERROR;	        
      }        
      modalDampingValues[i] = factor;    
    } 

  } else {

    //  read in one & set all factors to that value
    if (Tcl_GetDouble(interp, argv[1], &factor) != TCL_OK) {
      opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaK? \n";
      return TCL_ERROR;	        
    } 

    for (int i=0; i<numEigen; i++)
      modalDampingValues[i] = factor;
  } 

  // set factors in domain
  theDomain.setModalDampingFactors(&modalDampingValues, false);
  return TCL_OK;
}





int 
setElementRayleighDampingFactors(ClientData clientData, 
				 Tcl_Interp *interp, 
				 int argc, 
				 TCL_Char **argv)
{
  if (argc < 6) { 
    opserr << "WARNING setElementRayleighDampingFactors eleTag? alphaM? betaK? betaK0? betaKc? - not enough arguments to command\n";
    return TCL_ERROR;
  }
  int eleTag;
  double alphaM, betaK, betaK0, betaKc;

  if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read eleTag? \n";
    return TCL_ERROR;	        
  }    


  if (Tcl_GetDouble(interp, argv[2], &alphaM) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read alphaM? \n";
    return TCL_ERROR;	        
  }    
  if (Tcl_GetDouble(interp, argv[3], &betaK) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaK? \n";
    return TCL_ERROR;	        
  }        
  if (Tcl_GetDouble(interp, argv[4], &betaK0) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaK0? \n";
    return TCL_ERROR;	        
  }        
  if (Tcl_GetDouble(interp, argv[5], &betaKc) != TCL_OK) {
    opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read betaKc? \n";
    return TCL_ERROR;	        
  }        

  Element *theEle = theDomain.getElement(eleTag);
  theEle->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  return TCL_OK;
}


extern int
TclAddMeshRegion(ClientData clientData, Tcl_Interp *interp, int argc, 
	     TCL_Char **argv, Domain &theDomain);

int 
addRegion(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, &theDomain);
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
  bool echo = true;

  int cArg = 2;
  while (cArg < argc) {
    if (strcmp(argv[cArg],"-append") == 0) 
      mode = APPEND;
    if (strcmp(argv[cArg],"-noEcho") == 0) 
      echo = false;
    cArg++;
  }

  if (opserr.setFile(argv[1], mode, echo) < 0) 
    opserr << "WARNING logFile " << argv[1] << " failed to set the file\n";

  const char *pwd = getInterpPWD(interp);  
  simulationInfo.addOutputFile(argv[1], pwd);

  return TCL_OK;
}


int 
setPrecision(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

  if (argc < 2) { 
    opserr << "WARNING setPrecision precision? - no precision value supplied\n";
    return TCL_ERROR;
  }
  int precision;
  if (Tcl_GetInt(interp, argv[1], &precision) != TCL_OK) {
    opserr << "WARNING setPrecision precision? - error reading precision value supplied\n";
    return TCL_ERROR;
  }
  opserr.setPrecision(precision);

  return TCL_OK;
}


int 
exit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  Tcl_Finalize();
  return TCL_OK;
}



int 
getPID(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int pid = 0;
#ifdef _PARALLEL_INTERPRETERS
  if (theMachineBroker != 0)
    pid =  theMachineBroker->getPID();
#endif

#ifdef _PARALLEL_PROCESSING
  if (theMachineBroker != 0)
    pid =  theMachineBroker->getPID();
#endif

  // now we copy the value to the tcl string that is returned
  char buffer[30];
  sprintf(buffer,"%d",pid);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;  
}


int 
getNP(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int np = 1;
#ifdef _PARALLEL_INTERPRETERS
  if (theMachineBroker != 0)
    np =  theMachineBroker->getNP();
#endif

#ifdef _PARALLEL_PROCESSING
  if (theMachineBroker != 0)
    np =  theMachineBroker->getNP();
#endif

  // now we copy the value to the tcl string that is returned
  char buffer[30];
  sprintf(buffer,"%d",np);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;  
}

int
getEleTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  Element *theEle;
  ElementIter &eleIter = theDomain.getElements();

  char buffer[20];

  while ((theEle = eleIter()) != 0) {
    sprintf(buffer, "%d ", theEle->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
getNodeTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  Node *theEle;
  NodeIter &eleIter = theDomain.getNodes();

  char buffer[20];

  while ((theEle = eleIter()) != 0) {
    sprintf(buffer, "%d ", theEle->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
getParamTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  Parameter *theEle;
  ParameterIter &eleIter = theDomain.getParameters();

  char buffer[20];

  while ((theEle = eleIter()) != 0) {
    sprintf(buffer, "%d ", theEle->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
getParamValue(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 2) {
    opserr << "Insufficient arguments to getParamValue" << endln;
    return TCL_ERROR;
  }

  int paramTag;

  if (Tcl_GetInt(interp, argv[1], &paramTag) != TCL_OK) {
    opserr << "WARNING getParamValue -- could not read paramTag \n";
    return TCL_ERROR;	        
  }

  Parameter *theEle = theDomain.getParameter(paramTag);

  char buffer[40];

  sprintf(buffer, "%35.20f", theEle->getValue());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
sdfResponse(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 9) {
    opserr << "Insufficient arguments to sdfResponse" << endln;
    return TCL_ERROR;
  }

  double m, zeta, k, Fy, alpha, dtF, dt;
  if (Tcl_GetDouble(interp, argv[1], &m) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read mass \n";
    return TCL_ERROR;	        
  }
  if (Tcl_GetDouble(interp, argv[2], &zeta) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read zeta \n";
    return TCL_ERROR;	        
  }
  if (Tcl_GetDouble(interp, argv[3], &k) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read k \n";
    return TCL_ERROR;	        
  }
  if (Tcl_GetDouble(interp, argv[4], &Fy) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read Fy \n";
    return TCL_ERROR;	        
  }
  if (Tcl_GetDouble(interp, argv[5], &alpha) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read alpha \n";
    return TCL_ERROR;	        
  }
  if (Tcl_GetDouble(interp, argv[6], &dtF) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read dtF \n";
    return TCL_ERROR;	        
  }
  if (Tcl_GetDouble(interp, argv[8], &dt) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read dt \n";
    return TCL_ERROR;	        
  }
  double uresidual = 0.0;
  double umaxprev = 0.0;
  if (argc > 9) {
    if (Tcl_GetDouble(interp, argv[9], &uresidual) != TCL_OK) {
      opserr << "WARNING sdfResponse -- could not read uresidual \n";
      return TCL_ERROR;	        
    }
    if (Tcl_GetDouble(interp, argv[10], &umaxprev) != TCL_OK) {
      opserr << "WARNING sdfResponse -- could not read umaxprev \n";
      return TCL_ERROR;	        
    }    
  }

  double gamma = 0.5;
  double beta = 0.25;
  double tol = 1.0e-8;
  int maxIter = 10;

  std::ifstream infile(argv[7]);
 
  double c = zeta*2*sqrt(k*m);
  double Hkin = alpha/(1.0-alpha)*k;

  double p0 = 0.0;
  double u0 = uresidual;
  double v0 = 0.0;
  double fs0 = 0.0;
  double a0 = (p0-c*v0-fs0)/m;

  double a1 = m/(beta*dt*dt) + (gamma/(beta*dt))*c;
  double a2 = m/(beta*dt) + (gamma/beta-1.0)*c;
  double a3 = (0.5/beta-1.0)*m + dt*(0.5*gamma/beta-1.0)*c;

  double au = 1.0/(beta*dt*dt);
  double av = 1.0/(beta*dt);
  double aa = 0.5/beta-1.0;

  double vu = gamma/(beta*dt);
  double vv = 1.0-gamma/beta;
  double va = dt*(1-0.5*gamma/beta);
    
  double kT0 = k;

  double umax = fabs(umaxprev);
  double amax = 0.0; double tamax = 0.0;
  double up = uresidual; double up0 = up;
  int i = 0;
  double ft, u, du, v, a, fs, zs, ftrial, kT, kTeff, dg, phat, R, R0, accel;
  while (infile >> ft) {
    i++;
    
    u = u0;
      
    fs = fs0;
    kT = kT0;
    up = up0;
      
    phat = ft + a1*u0 + a2*v0 + a3*a0;
      
    R = phat - fs - a1*u;
    R0 = R;
    if (R0 == 0.0) {
      R0 = 1.0;
    }
    
    int iter = 0;

    while (iter < maxIter && fabs(R/R0) > tol) {
      iter++;

      kTeff = kT + a1;

      du = R/kTeff;

      u = u + du;

      fs = k*(u-up0);
      zs = fs-Hkin*up0;
      ftrial = fabs(zs)-Fy;
      if (ftrial > 0) {
	dg = ftrial/(k+Hkin);
	if (fs < 0) {
	  fs = fs + dg*k;
	  up = up0 - dg;
	} else {
	  fs = fs - dg*k;
	  up = up0 + dg;
	}
	kT = k*Hkin/(k+Hkin);
      } else {
	kT = k;
      }
      
      R = phat - fs - a1*u;
    }

    v = vu*(u-u0) + vv*v0 + va*a0;
    a = au*(u-u0) - av*v0 - aa*a0;

    u0 = u;
    v0 = v;
    a0 = a;
    fs0 = fs;
    kT0 = kT;
    up0 = up;

    if (fabs(u) > umax) {
      umax = fabs(u);
    }
    if (fabs(a) > amax) {
      amax = fabs(a);
      tamax = iter*dt;
    }
  }
  
  infile.close();
    
  
  char buffer[80];
  sprintf(buffer, "%f %f %f %f %f", umax, u, up, amax, tamax);

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int 
opsBarrier(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
#ifdef _PARALLEL_INTERPRETERS
  return MPI_Barrier(MPI_COMM_WORLD);
#endif

  return TCL_OK;  
}


int 
opsSend(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
#ifdef _PARALLEL_INTERPRETERS
  if (argc < 2)
    return TCL_OK;

  int otherPID = -1;
  int myPID =  theMachineBroker->getPID();
  int np =  theMachineBroker->getNP();
  const char *dataToSend = argv[argc-1];
  int msgLength = strlen(dataToSend)+1;
  
  const char *gMsg = dataToSend;
  //  strcpy(gMsg, dataToSend);
  
  if (strcmp(argv[1], "-pid") == 0 && argc > 3) {

    if (Tcl_GetInt(interp, argv[2], &otherPID) != TCL_OK) {
      opserr << "send -pid pid? data? - pid: " << argv[2] << " invalid\n";
      return TCL_ERROR;
    }

    if (otherPID > -1 && otherPID != myPID && otherPID < np) {
      
      MPI_Send((void *)(&msgLength), 1, MPI_INT, otherPID, 0, MPI_COMM_WORLD);
      MPI_Send((void *)gMsg, msgLength, MPI_CHAR, otherPID, 1, MPI_COMM_WORLD);

    } else {
      opserr << "send -pid pid? data? - pid: " << otherPID << " invalid\n";
      return TCL_ERROR;
    }

  } else {
    if (myPID == 0) { 
      MPI_Bcast((void *)(&msgLength), 1, MPI_INT,  0, MPI_COMM_WORLD);
      MPI_Bcast((void *)gMsg, msgLength, MPI_CHAR, 0, MPI_COMM_WORLD);
    } else {
      opserr << "send data - only process 0 can do a broadcast - you may need to kill the application";
      return TCL_ERROR;
    }
  }

#endif

  return TCL_OK;  
}


int 
opsRecv(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
#ifdef _PARALLEL_INTERPRETERS
  if (argc < 2)
    return TCL_OK;

  int otherPID = 0;
  int myPID =  theMachineBroker->getPID();
  int np =  theMachineBroker->getNP();
  TCL_Char *varToSet = argv[argc-1];

  int msgLength = 0;
  char *gMsg = 0;

  if (strcmp(argv[1], "-pid") == 0 && argc > 3) {

    bool fromAny = false;

    if ((strcmp(argv[2], "ANY") == 0) || (strcmp(argv[2],"ANY_SOURCE") == 0) ||
	(strcmp(argv[2], "MPI_ANY_SOURCE") == 0)) {
      fromAny = true;
    } else {
      if (Tcl_GetInt(interp, argv[2], &otherPID) != TCL_OK) {
	opserr << "recv -pid pid? data? - pid: " << argv[2] << " invalid\n";
	return TCL_ERROR;
      }
	}

    if (otherPID > -1 && otherPID < np) {
      MPI_Status status;
      
      if (fromAny == false)
		if (myPID != otherPID)
		MPI_Recv((void *)(&msgLength), 1, MPI_INT, otherPID, 0, MPI_COMM_WORLD, &status);
        else {
	  opserr << "recv -pid pid? data? - " << otherPID << " cant receive from self!\n";
	  return TCL_ERROR;
	}
      else {
	MPI_Recv((void *)(&msgLength), 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
	otherPID = status.MPI_SOURCE;
	  }

      if (msgLength > 0) {
	gMsg = new char [msgLength];
	
	  MPI_Recv((void *)gMsg, msgLength, MPI_CHAR, otherPID, 1, MPI_COMM_WORLD, &status);

	Tcl_SetVar(interp, varToSet, gMsg, TCL_LEAVE_ERR_MSG);
      }

    } else {
      opserr << "recv -pid pid? data? - " << otherPID << " invalid\n";
      return TCL_ERROR;
    }
  } else {

    if (myPID != 0) {
      MPI_Bcast((void *)(&msgLength), 1, MPI_INT, 0, MPI_COMM_WORLD);

      if (msgLength > 0) {
	gMsg = new char [msgLength];

	MPI_Bcast((void *)gMsg, msgLength, MPI_CHAR, 0, MPI_COMM_WORLD);

	Tcl_SetVar(interp, varToSet, gMsg, TCL_LEAVE_ERR_MSG);
      }

    } else {
      opserr << "recv data - only process 0 can do a broadcast - you may need to kill the application";
      return TCL_ERROR;
    }
  }

#endif

  return TCL_OK;  
}


int
defaultUnits(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
    if (argc < 7) {
        opserr << "defaultUnits - missing a unit type want: defaultUnits -Force type? -Length type? -Time type?\n";
        return -1;
    }

    const char *force = 0;
    const char *length = 0;
    const char *time = 0;
    const char *temperature = "N/A";

    int count = 1;
    while (count < argc) {
        if ((strcmp(argv[count], "-force") == 0) || (strcmp(argv[count], "-Force") == 0)
            || (strcmp(argv[count], "-FORCE") == 0)) {
            force = argv[count + 1];
        }
        else if ((strcmp(argv[count], "-length") == 0) || (strcmp(argv[count], "-Length") == 0)
            || (strcmp(argv[count], "-LENGTH") == 0)) {
            length = argv[count + 1];
        }
        else if ((strcmp(argv[count], "-time") == 0) || (strcmp(argv[count], "-Time") == 0)
            || (strcmp(argv[count], "-TIME") == 0)) {
            time = argv[count + 1];
        }
        else if ((strcmp(argv[count], "-temperature") == 0) || (strcmp(argv[count], "-Temperature") == 0)
            || (strcmp(argv[count], "-TEMPERATURE") == 0) || (strcmp(argv[count], "-temp") == 0)
            || (strcmp(argv[count], "-Temp") == 0) || (strcmp(argv[count], "-TEMP") == 0)) {
            temperature = argv[count + 1];
        }
        else {
            opserr << "defaultUnits - unrecognized unit: " << argv[count] << " want: defaultUnits -Force type? -Length type? -Time type?\n";
            return -1;
        }
        count += 2;
    }

    if (length == 0 || force == 0 || time == 0) {
        opserr << "defaultUnits - missing a unit type want: defaultUnits -Force type? -Length type? -Time type?\n";
        return -1;
    }

    double lb, kip, n, kn, mn, kgf, tonf;
    double in, ft, mm, cm, m;
    double sec, msec;

    if ((strcmp(force, "lb") == 0) || (strcmp(force, "lbs") == 0)) {
        lb = 1.0;
    }
    else if ((strcmp(force, "kip") == 0) || (strcmp(force, "kips") == 0)) {
        lb = 0.001;
    }
    else if ((strcmp(force, "N") == 0)) {
        lb = 4.4482216152605;
    }
    else if ((strcmp(force, "kN") == 0) || (strcmp(force, "KN") == 0) || (strcmp(force, "kn") == 0)) {
        lb = 0.0044482216152605;
    }
    else if ((strcmp(force, "mN") == 0) || (strcmp(force, "MN") == 0) || (strcmp(force, "mn") == 0)) {
        lb = 0.0000044482216152605;
    }
    else if ((strcmp(force, "kgf") == 0)) {
        lb = 4.4482216152605 / 9.80665;
    }
    else if ((strcmp(force, "tonf") == 0)) {
        lb = 4.4482216152605 / 9.80665 / 1000.0;
    }
    else {
        lb = 1.0;
        opserr << "defaultUnits - unknown force type, valid options: lb, kip, N, kN, MN, kgf, tonf\n";
        return TCL_ERROR;
    }

    if ((strcmp(length, "in") == 0) || (strcmp(length, "inch") == 0)) {
        in = 1.0;
    }
    else if ((strcmp(length, "ft") == 0) || (strcmp(length, "feet") == 0)) {
        in = 1.0 / 12.0;
    }
    else if ((strcmp(length, "mm") == 0)) {
        in = 25.4;
    }
    else if ((strcmp(length, "cm") == 0)) {
        in = 2.54;
    }
    else if ((strcmp(length, "m") == 0)) {
        in = 0.0254;
    }
    else {
        in = 1.0;
        opserr << "defaultUnits - unknown length type, valid options: in, ft, mm, cm, m\n";
        return TCL_ERROR;
    }

    if ((strcmp(time, "sec") == 0) || (strcmp(time, "Sec") == 0)) {
        sec = 1.0;
    }
    else if ((strcmp(time, "msec") == 0) || (strcmp(time, "mSec") == 0)) {
        sec = 1000.0;
    }
    else {
        sec = 1.0;
        opserr << "defaultUnits - unknown time type, valid options: sec, msec\n";
        return TCL_ERROR;
    }

    kip = lb / 0.001;
    n = lb / 4.4482216152605;
    kn = lb / 0.0044482216152605;
    mn = lb / 0.0000044482216152605;
    kgf = lb / (4.4482216152605 / 9.80665);
    tonf = lb / (4.4482216152605 / 9.80665 / 1000.0);

    ft = in * 12.0;
    mm = in / 25.4;
    cm = in / 2.54;
    m = in / 0.0254;

    msec = sec * 0.001;

    char string[50];

    sprintf(string, "set lb %.18e", lb);   Tcl_Eval(interp, string);
    sprintf(string, "set lbf %.18e", lb);   Tcl_Eval(interp, string);
    sprintf(string, "set kip %.18e", kip);   Tcl_Eval(interp, string);
    sprintf(string, "set N %.18e", n);   Tcl_Eval(interp, string);
    sprintf(string, "set kN %.18e", kn);   Tcl_Eval(interp, string);
    sprintf(string, "set Newton %.18e", n);   Tcl_Eval(interp, string);
    sprintf(string, "set kNewton %.18e", kn);   Tcl_Eval(interp, string);
    sprintf(string, "set MN %.18e", mn);   Tcl_Eval(interp, string);
    sprintf(string, "set kgf %.18e", kgf);   Tcl_Eval(interp, string);
    sprintf(string, "set tonf %.18e", tonf);   Tcl_Eval(interp, string);

    sprintf(string, "set in %.18e", in);   Tcl_Eval(interp, string);
    sprintf(string, "set inch %.18e", in);   Tcl_Eval(interp, string);
    sprintf(string, "set ft %.18e", ft);   Tcl_Eval(interp, string);
    sprintf(string, "set mm %.18e", mm);   Tcl_Eval(interp, string);
    sprintf(string, "set cm %.18e", cm);   Tcl_Eval(interp, string);
    sprintf(string, "set m  %.18e", m);   Tcl_Eval(interp, string);
    sprintf(string, "set meter  %.18e", m);   Tcl_Eval(interp, string);

    sprintf(string, "set sec %.18e", sec);   Tcl_Eval(interp, string);
    sprintf(string, "set msec %.18e", msec);   Tcl_Eval(interp, string);

    double g = 32.174049*ft / (sec*sec);
    sprintf(string, "set g %.18e", g);   Tcl_Eval(interp, string);
    sprintf(string, "set kg %.18e", n*sec*sec / m);   Tcl_Eval(interp, string);
    sprintf(string, "set Mg %.18e", 1e3*n*sec*sec / m);   Tcl_Eval(interp, string);
    sprintf(string, "set slug %.18e", lb*sec*sec / ft);   Tcl_Eval(interp, string);
    sprintf(string, "set Pa %.18e", n / (m*m));   Tcl_Eval(interp, string);
    sprintf(string, "set kPa %.18e", 1e3*n / (m*m));   Tcl_Eval(interp, string);
    sprintf(string, "set MPa %.18e", 1e6*n / (m*m));   Tcl_Eval(interp, string);
    sprintf(string, "set psi %.18e", lb / (in*in));   Tcl_Eval(interp, string);
    sprintf(string, "set ksi %.18e", kip / (in*in));   Tcl_Eval(interp, string);
    sprintf(string, "set psf %.18e", lb / (ft*ft));   Tcl_Eval(interp, string);
    sprintf(string, "set ksf %.18e", kip / (ft*ft));   Tcl_Eval(interp, string);
    sprintf(string, "set pcf %.18e", lb / (ft*ft*ft));   Tcl_Eval(interp, string);
    sprintf(string, "set in2 %.18e", in*in);   Tcl_Eval(interp, string);
    sprintf(string, "set ft2 %.18e", ft*ft);   Tcl_Eval(interp, string);
    sprintf(string, "set mm2 %.18e", mm*mm);   Tcl_Eval(interp, string);
    sprintf(string, "set cm2 %.18e", cm*cm);   Tcl_Eval(interp, string);
    sprintf(string, "set m2 %.18e", m*m);   Tcl_Eval(interp, string);
    sprintf(string, "set in4 %.18e", in*in*in*in);   Tcl_Eval(interp, string);
    sprintf(string, "set ft4 %.18e", ft*ft*ft*ft);   Tcl_Eval(interp, string);
    sprintf(string, "set mm4 %.18e", mm*mm*mm*mm);   Tcl_Eval(interp, string);
    sprintf(string, "set cm4 %.18e", cm*cm*cm*cm);   Tcl_Eval(interp, string);
    sprintf(string, "set m4 %.18e", m*m*m*m);   Tcl_Eval(interp, string);
    sprintf(string, "set pi %.18e", 2.0*asin(1.0));   Tcl_Eval(interp, string);
    sprintf(string, "set PI %.18e", 2.0*asin(1.0));   Tcl_Eval(interp, string);

    int res = simulationInfo.setForceUnit(force);
    res += simulationInfo.setLengthUnit(length);
    res += simulationInfo.setTimeUnit(time);
    res += simulationInfo.setTemperatureUnit(temperature);

    return res;
}



const char * getInterpPWD(Tcl_Interp *interp) {
  static char *pwd = 0;

  if (pwd != 0)
    delete [] pwd;

#ifdef _TCL84
  Tcl_Obj *cwd = Tcl_FSGetCwd(interp);
  if (cwd != NULL) {
    int length;
    const char *objPWD = Tcl_GetStringFromObj(cwd, &length);
    pwd = new char[length+1];
    strcpy(pwd, objPWD);
    Tcl_DecrRefCount(cwd);	
  }
#else

  Tcl_DString buf;
  const char *objPWD = Tcl_GetCwd(interp, &buf);

  pwd = new char[strlen(objPWD)+1];
  strcpy(pwd, objPWD);

  Tcl_DStringFree(&buf);

#endif
  return pwd;
}


int 
OpenSeesExit(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  theDomain.clearAll();

#ifdef _PARALLEL_PROCESSING
  //
  // mpi clean up
  //

  if (theMachineBroker != 0) {
    theMachineBroker->shutdown();
    fprintf(stderr, "Process Terminating\n");
    //    delete theMachineBroker;
    //    theMachineBroker = 0;
  }
  MPI_Finalize();
#endif

#ifdef _PARALLEL_INTERPRETERS
  //
  // mpi clean up
  //

  if (theMachineBroker != 0) {
    theMachineBroker->shutdown();
    fprintf(stderr, "Process Terminating\n");
    //    delete theMachineBroker;
    //    theMachineBroker = 0;
  }
  MPI_Finalize();
#endif

  if (simulationInfoOutputFilename != 0) {
    simulationInfo.end();
    XmlFileStream simulationInfoOutputFile;
    simulationInfoOutputFile.setFile(simulationInfoOutputFilename);
    simulationInfoOutputFile.open();
    simulationInfoOutputFile << simulationInfo;
    simulationInfoOutputFile.close();
    simulationInfoOutputFilename = 0;
  }

  int returnCode = 0;
  if (argc > 1) {
    if (Tcl_GetInt(interp, argv[1], &returnCode) != TCL_OK) {
      opserr << "WARNING: OpenSeesExit - failed to read return code\n";
    }
  }
  Tcl_Exit(returnCode);

  return 0;
}



int stripOpenSeesXML(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

  if (argc < 3) {
    opserr << "ERROR incorrect # args - stripXML input.xml output.dat <output.xml>\n";
    return -1;
  }

  const char *inputFile = argv[1];
  const char *outputDataFile = argv[2];
  const char *outputDescriptiveFile = 0;

  if (argc == 4)
    outputDescriptiveFile = argv[3];
  

  // open files
  ifstream theInputFile; 
  theInputFile.open(inputFile, ios::in);  
  if (theInputFile.bad()) {
    opserr << "stripXML - error opening input file: " << inputFile << endln;
    return -1;
  }

  ofstream theOutputDataFile; 
  theOutputDataFile.open(outputDataFile, ios::out);  
  if (theOutputDataFile.bad()) {
    opserr << "stripXML - error opening input file: " << outputDataFile << endln;
    return -1;
  }

  ofstream theOutputDescriptiveFile;
  if (outputDescriptiveFile != 0) {
    theOutputDescriptiveFile.open(outputDescriptiveFile, ios::out);  
    if (theOutputDescriptiveFile.bad()) {
      opserr << "stripXML - error opening input file: " << outputDescriptiveFile << endln;
      return -1;
    }
  }

  string line;
  bool spitData = false;
  while (! theInputFile.eof() ) {
    getline(theInputFile, line);
    const char *inputLine = line.c_str();

    if (spitData == true) {
      if (strstr(inputLine,"</Data>") != 0) 
	spitData = false;
      else 
	; //	theOutputDataFile << line << endln;
    } else {
      const char *inputLine = line.c_str();
      if (strstr(inputLine,"<Data>") != 0) 
	spitData = true;
      else if (outputDescriptiveFile != 0)
	; // theOutputDescriptiveFile << line << endln;
    }
  }      
  
  theInputFile.close();
  theOutputDataFile.close();

  if (outputDescriptiveFile != 0)
    theOutputDescriptiveFile.close();

  return 0;
}

extern int binaryToText(const char *inputFilename, const char *outputFilename);
extern int textToBinary(const char *inputFilename, const char *outputFilename);

int convertBinaryToText(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 3) {
    opserr << "ERROR incorrect # args - convertBinaryToText inputFile outputFile\n";
    return -1;
  }

  const char *inputFile = argv[1];
  const char *outputFile = argv[2];

  return binaryToText(inputFile, outputFile);
}


int convertTextToBinary(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  if (argc < 3) {
    opserr << "ERROR incorrect # args - convertTextToBinary inputFile outputFile\n";
    return -1;
  }

  const char *inputFile = argv[1];
  const char *outputFile = argv[2];

  return textToBinary(inputFile, outputFile);
}

int domainChange(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  theDomain.domainChange();
  return TCL_OK;
}


int record(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  theDomain.record(false);
  return TCL_OK;
}


extern 
int peerSearchNGA(const char *eq,
		  const char *soilType,
		  const char *fault,
		  const char *magLo,
		  const char *magHi,
		  const char *distLo,
		  const char *distHi,
		  const char *vsLo,
		  const char *vsHi,
		  const char *pgaLo,
		  const char *pgaHi,
		  const char *latSW,
		  const char *latNE,
		  const char *lngSW,
		  const char *lngNW,
		  StringContainer &recordNames);

int 
peerNGA(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  StringContainer ngaRecordNames;
  const char *eq =0;
  const char *soilType = 0;
  const char *fault =0;
  const char *magLo =0;
  const char *magHi =0;
  const char *distLo =0;
  const char *distHi =0;
  const char *vsLo =0;
  const char *vsHi =0;
  const char *pgaLo =0;
  const char *pgaHi =0;
  const char *latSW =0;
  const char *latNE =0;
  const char *lngSW =0;
  const char *lngNW =0;

  int currentArg = 1;
  while (currentArg+1 < argc) {
    if (strcmp(argv[currentArg],"-eq") == 0) {
      eq = argv[currentArg+1];
    } else if (strcmp(argv[currentArg],"-fault") == 0) {
      fault = argv[currentArg+1];
    } else if (strcmp(argv[currentArg],"-soil") == 0) {
      soilType = argv[currentArg+1];
    } else if (strcmp(argv[currentArg],"-magLo") == 0) {
      magLo = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-magHi") == 0) {
      magHi = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-distLo") == 0) {
      distLo  = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-distHi") == 0) {
      distHi = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-vsLo") == 0) {
      vsLo = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-vsHi") == 0) {
      vsHi = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-pgaLo") == 0) {
      pgaLo = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-pgaHi") == 0) {
      pgaHi = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-latSW") == 0) {
      latSW = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-latNE") == 0) {
      latNE = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-lngSW") == 0) {
      lngSW = argv[currentArg+1];      
    } else if (strcmp(argv[currentArg],"-lngNW") == 0) {
      lngNW = argv[currentArg+1];     
    } 
    // unrecognized
    currentArg+=2;
  }        

  peerSearchNGA(eq,
		soilType,
		fault,
		magLo,
		magHi,
		distLo,
		distHi,
		vsLo,
		vsHi,
		pgaLo,
		pgaHi,
		latSW,
		latNE,
		lngSW,
		lngNW,
		ngaRecordNames);

  int numStrings = ngaRecordNames.getNumStrings();
  for (int i=0; i<numStrings; i++) {
    Tcl_AppendResult(interp, ngaRecordNames.getString(i), NULL);
    Tcl_AppendResult(interp, " ", NULL);
  }

  return TCL_OK;
}

int
totalCPU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%f", theAlgorithm->getTotalTimeCPU());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
solveCPU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%f", theAlgorithm->getSolveTimeCPU());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
accelCPU(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%f", theAlgorithm->getAccelTimeCPU());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
numFact(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%d", theAlgorithm->getNumFactorizations());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
systemSize(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theSOE == 0) {
    sprintf(buffer, "NO SYSTEM SET");
    return TCL_OK;
  }


  sprintf(buffer, "%d", theSOE->getNumEqn());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
numIter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];

  if (theAlgorithm == 0)
    return TCL_ERROR;

  sprintf(buffer, "%d", theAlgorithm->getNumIterations());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
elementActivate(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
	int eleTag;
	int argLoc = 1;
	int Nelements = argc;
	ID activate_us(0, Nelements);

	while (argLoc < argc && Tcl_GetInt(interp, argv[argLoc], &eleTag) == TCL_OK) {
		activate_us.insert(eleTag);
		++argLoc;
	}

	theDomain.activateElements(activate_us);

  return TCL_OK;
}

int 
elementDeactivate(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

	int eleTag;
	int argLoc = 1;
	int Nelements = argc;
	ID deactivate_us(0, Nelements);

	while (argLoc < argc && Tcl_GetInt(interp, argv[argLoc], &eleTag) == TCL_OK) {
		deactivate_us.insert(eleTag);
		++argLoc;
	}

	theDomain.deactivateElements(deactivate_us);
  return TCL_OK;}

int
version(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  char buffer[20];


  sprintf(buffer, "%s", OPS_VERSION);
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

extern "C" int OpenSeesParseArgv(int argc, char **argv)
{
   if (argc > 1) {
	int currentArg = 1;
	while (currentArg < argc && argv[currentArg] != NULL) {

	  if ((strcmp(argv[currentArg], "-par") == 0) || (strcmp(argv[currentArg], "-Par") == 0)) {
	    
	    if (argc > (currentArg+2)) {
	      
	      char *parName = argv[currentArg+1];
	      char *parValue = argv[currentArg+2];
	      
	      // add a OpenSeesTcl_Parameter to end of list of parameters
	      OpenSeesTcl_Parameter *nextParam = new OpenSeesTcl_Parameter;
	      nextParam->name = new char [strlen(parName)+1];
	      strcpy(nextParam->name, parName);
	      nextParam->values = 0;
	      
	      if (theParameters == 0)
		theParameters = nextParam;
	      if (endParameters != 0)
		endParameters->next = nextParam;
	      nextParam->next = 0;
	      endParameters = nextParam;
	      
	      // now open par values files to create the values
	      char nextLine[1000];
	      FILE *valueFP = fopen(parValue,"r");
	      if (valueFP != 0) {
		OpenSeesTcl_ParameterValues *endValues = 0;
		
		while (fscanf(valueFP, "%s", nextLine) != EOF) {
		  
		  OpenSeesTcl_ParameterValues *nextValue = new OpenSeesTcl_ParameterValues;
		  nextValue->value = new char [strlen(nextLine)+1];
		  strcpy(nextValue->value, nextLine);
		  
		  if (nextParam->values == 0) {
		    nextParam->values = nextValue;
		  }
		if (endValues != 0)
		  endValues->next = nextValue;
		endValues = nextValue;
		nextValue->next = 0;	      
		}
		fclose(valueFP);
	      } else {
		
		OpenSeesTcl_ParameterValues *nextValue = new OpenSeesTcl_ParameterValues;		
		nextValue->value = new char [strlen(parValue)+1];
		
		strcpy(nextValue->value, parValue);
		
		nextParam->values = nextValue;
		nextValue->next = 0;
		
	      }
	      numParam++;
	    }
	    currentArg += 3;
	  } else if ((strcmp(argv[currentArg], "-info") == 0) || (strcmp(argv[currentArg], "-INFO") == 0)) {
	    if (argc > (currentArg+1)) {
	      simulationInfoOutputFilename = argv[currentArg+1];	    
	    }			   
	    currentArg+=2;
	  } else 
	    currentArg++;
	}
   }
   if (numParam != 0) {
     paramNames = new char *[numParam];
	 paramValues = new char *[numParam];
   }
	return numParam;
}


extern "C" int
EvalFileWithParameters(Tcl_Interp *interp, 
				char *tclStartupFileScript, 
				OpenSeesTcl_Parameter *theInputParameters,  
				int currentParam, 
				int rank, 
				int np)
{
  if (theInputParameters == 0)
	  theInputParameters = theParameters;

  if (currentParam < numParam) {
    OpenSeesTcl_Parameter *theCurrentParam = theInputParameters;
    OpenSeesTcl_Parameter *theNextParam = theParameters->next;
    char *paramName = theCurrentParam->name;
    paramNames[currentParam] = paramName;

    OpenSeesTcl_ParameterValues *theValue = theCurrentParam->values;
    int nextParam = currentParam+1;
    while (theValue != 0) {
      char *paramValue = theValue->value;
      paramValues[currentParam] = paramValue;
      EvalFileWithParameters(interp, 
			     tclStartupFileScript, 
			     theNextParam, 
			     nextParam, 
			     rank, 
			     np);

      theValue=theValue->next;
    } 
  } else {
    
    simulationInfo.start();
    static int count = 0;
    
    if ((count % np) == rank) {
      Tcl_Eval(interp, "wipe");
     
	  
      for (int i=0; i<numParam; i++) {
		  
	Tcl_SetVar(interp, paramNames[i], paramValues[i], TCL_GLOBAL_ONLY);	
    
	simulationInfo.addParameter(paramNames[i], paramValues[i]); 
     }

      count++;

      const char *pwd = getInterpPWD(interp);
      simulationInfo.addInputFile(tclStartupFileScript, pwd);

      int ok = Tcl_EvalFile(interp, tclStartupFileScript);

      simulationInfo.end();
      
      return ok;
    }
    else
      count++;
  }

  return 0;
}


int
setParameter(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int argLoc = 1;
  double newValue = 0.0;
  ID eleIDs(0, 32);
  int numEle = 0;
  int flag = 0;

  if (strstr(argv[argLoc],"-val") != 0) {
    if (Tcl_GetDouble(interp, argv[argLoc+1], &newValue) != TCL_OK) {
      opserr << "WARNING setParameter: invalid parameter value\n";
      return TCL_ERROR;
    } 
  } else {
    opserr << "WARNING setParameter:  -val not found " << endln;
    return TCL_ERROR;
  } 

  argLoc += 2;

  if (strstr(argv[argLoc],"-ele") != 0) {    

    if ((strcmp(argv[argLoc],"-ele") == 0) ||
	(strcmp(argv[argLoc],"-eles") == 0) ||
	(strcmp(argv[argLoc],"-element") == 0)) {
      
      //
      // read in a list of ele until end of command or other flag
      //
      
      argLoc++;
      int eleTag;
      
      while (argLoc < argc && Tcl_GetInt(interp, argv[argLoc], &eleTag) == TCL_OK) {
	eleIDs[numEle] = eleTag;
	numEle++;
	argLoc++;
      }

      if (numEle > 0)
	flag = 1;

    } else if (strcmp(argv[argLoc],"-eleRange") == 0) {
      
      flag = 2;

      // ensure no segmentation fault if user messes up
      if (argc < argLoc+3) {
	opserr << "WARNING recorder Element .. -eleRange start? end?  .. - no ele tags specified\n";
	return TCL_ERROR;
      }
      
      //
      // read in start and end tags of two elements & add set [start,end]
      //
      
      int start, end;
      if (Tcl_GetInt(interp, argv[argLoc+1], &start) != TCL_OK) {
	opserr << "WARNING recorder Element -eleRange start? end? - invalid start " << argv[argLoc+1] << endln;
	return TCL_ERROR;
      }      
      if (Tcl_GetInt(interp, argv[argLoc+2], &end) != TCL_OK) {
	opserr << "WARNING recorder Element -eleRange start? end? - invalid end " << argv[argLoc+2] << endln;
	return TCL_ERROR;
      }      
      if (start > end) {
	int swap = end;
	end = start;
	start = swap;
      }
      eleIDs[0] = start;
      eleIDs[1] = end;

      argLoc += 3;
    }

    ElementStateParameter theParameter(newValue, &argv[argLoc], argc-argLoc, flag, &eleIDs);

    theDomain.addParameter(&theParameter);
  }

  return TCL_OK;
}

int
maxOpenFiles(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  int maxOpenFiles;

  if (Tcl_GetInt(interp, argv[1], &maxOpenFiles) != TCL_OK) {
      return TCL_ERROR;
  } 

  #ifdef _WIN32
  int newMax = _setmaxstdio(maxOpenFiles);
  if (maxOpenFiles > 2048) {
		opserr << "setMaxOpenFiles: too many files specified (2048 max)\n";
  } else {
	 if (newMax != maxOpenFiles) {
		opserr << "setMaxOpenFiles FAILED: max allowed files: " << newMax;
		return TCL_ERROR;
	}
  }
  return TCL_OK;
  #endif

  opserr << "setMaxOpenFiles FAILED: - command not available on this machine\n";
  return TCL_OK;
}

// Talledo Start
int 
printModelGID(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
  // This function print's a file with node and elements in a format useful for GID
	int res = 0;
	bool hasLinear = 0;
	bool hasTri3  = 0;
	bool hasQuad4 = 0;
	bool hasQuad8 = 0;
	bool hasQuad9 = 0;
	bool hasBrick = 0;
	int startEle = 1;
	int endEle = 1;
	int eleRange = 0;
	int i = 2;

    FileStream outputFile;
    OPS_Stream *output = &opserr;

	if (argc < 2) {
		opserr << "WARNING printGID fileName? - no filename supplied\n";
		return TCL_ERROR;
	}
	openMode mode = OVERWRITE;
	if (argc >= 3)
	{
		if (strcmp(argv[i],"-append") == 0) 
		{
			mode = APPEND;
			i++;
		}
		if (strcmp(argv[i],"-eleRange") == 0) {
			//opserr<<"WARNING:commands: eleRange defined"<<endln;
			eleRange = 1;
			if (Tcl_GetInt(interp, argv[i+1], &startEle) != TCL_OK) {
				opserr << "WARNING print node failed to get integer: " << argv[i+1] << endln;
				return TCL_ERROR;
			}
			if (Tcl_GetInt(interp, argv[i+2], &endEle) != TCL_OK) {
				opserr << "WARNING print node failed to get integer: " << argv[i+2] << endln;
				return TCL_ERROR;
			}
			//opserr<<"startEle = "<<startEle<<" endEle = "<<endEle<<endln;
		}
	}
		
    if (outputFile.setFile(argv[1], mode) < 0) {
        opserr << "WARNING printGID " << argv[1] << " failed to set the file\n";
		return TCL_ERROR;
	}
    
	// Cycle over Elements to understand what type of elements are there
	ElementIter &theElements = theDomain.getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
		
		// Check type of Element with Number of Nodes
		// if 2 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 2) {
			hasLinear = 1;
		} else if (nNode == 4) {
			hasQuad4 = 1;
		} else if (nNode == 3) {
			hasTri3 = 1;
		} else if (nNode == 9) {
			hasQuad9 = 1;
		} else if (nNode == 8) {
			const char *name = theElement->getClassType();
			if (strcmp(name,"Brick") == 0) {
				hasBrick = 1;
			} else {
				hasQuad8 = 1;
			}
		}
	}
	// **** Linear Elements - 2 Nodes
	if (hasLinear == 1) {
		// Print HEADER
		outputFile << "MESH \"2NMESH\" dimension 3 ElemType Linear Nnode 2" << endln;
		outputFile << "#color 0 0 255" << endln << endln;
    
		// Print node coordinates
		outputFile << "Coordinates" << endln;
		NodeIter &theNodes = theDomain.getNodes();
		Node *theNode;
		while ((theNode = theNodes()) != 0) {
			int tag = theNode->getTag();
			const Vector &crds = theNode->getCrds();
			//outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
			int l_tmp = crds.Size();
			outputFile << tag << "\t\t";
			for (int ii = 0; ii<l_tmp; ii++) {
				outputFile << crds(ii) << "\t";
			};
			for (int ii = l_tmp; ii<3; ii++) {
				outputFile << 0.0 << "\t";
			};
			outputFile << endln;
		}
		outputFile << "End coordinates" << endln << endln;

		// Print elements connectivity
		outputFile << "Elements" << endln;
		ElementIter &theElements = theDomain.getElements();
		Element *theElement;
		while ((theElement = theElements()) != 0) {
			int tag = theElement->getTag();
			// Check if element tag is inside theRange
			if (((tag<=endEle) & (tag>=startEle)) || (eleRange == 0)) {
				// Check type of Element with Number of Nodes
				// if 2 Nodes print the Element
				int nNode = theElement->getNumExternalNodes();
				if (nNode == 2) {
					Node **NodePtrs;
					NodePtrs = theElement->getNodePtrs();		
					Vector tagNodes(nNode);
					for (int i = 0; i < nNode; i++) {
						tagNodes(i)=NodePtrs[i]->getTag();
					}
					outputFile << tag << "\t\t";
					for (int i = 0; i < nNode; i++) {
						outputFile << tagNodes(i) << "\t";
					}
					outputFile << endln;
				}
			}
		}
		outputFile << "End elements" << endln;
	}
	// **** Quadrilateral Elements - 4 Nodes
	if (hasQuad4 == 1) {
		// Print HEADER
		outputFile << "MESH \"4NMESH\" dimension 3 ElemType Quadrilateral Nnode 4" << endln;
		outputFile << "#color 0 255 0" << endln << endln;
    
		// Print node coordinates
		outputFile << "Coordinates" << endln;
		NodeIter &theNodes = theDomain.getNodes();
		Node *theNode;
		while ((theNode = theNodes()) != 0) {
			int tag = theNode->getTag();
			const Vector &crds = theNode->getCrds();
			//outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
			int l_tmp = crds.Size();
			outputFile << tag << "\t\t";
			for (int ii = 0; ii<l_tmp; ii++) {
				outputFile << crds(ii) << "\t";
			};
			for (int ii = l_tmp; ii<3; ii++) {
				outputFile << 0.0 << "\t";
			};
			outputFile << endln;
		}
		outputFile << "End coordinates" << endln << endln;

		// Print elements connectivity
		outputFile << "Elements" << endln;
		ElementIter &theElements = theDomain.getElements();
		Element *theElement;
		while ((theElement = theElements()) != 0) {
			int tag = theElement->getTag();
			// Check if element tag is inside theRange
			if (((tag<=endEle) & (tag>=startEle)) || (eleRange == 0)) {
		
				// Check type of Element with Number of Nodes
				// if 2 Nodes print the Element
				int nNode = theElement->getNumExternalNodes();
				if (nNode == 4) {
					Node **NodePtrs;
					NodePtrs = theElement->getNodePtrs();		
					Vector tagNodes(nNode);
					for (int i = 0; i < nNode; i++) {
						tagNodes(i)=NodePtrs[i]->getTag();
					}
					outputFile << tag << "\t\t";
					for (int i = 0; i < nNode; i++) {
						outputFile << tagNodes(i) << "\t";
					}
					outputFile << endln;
				}
			}
		}
		outputFile << "End elements" << endln;
	}
	// **** Triangular Elements - 3 Nodes
	if (hasTri3 == 1) {
		// Print HEADER
		outputFile << "MESH \"3NMESH\" dimension 3 ElemType Triangle Nnode 3" << endln;
		outputFile << "#color 0 255 0" << endln << endln;
    
		// Print node coordinates
		outputFile << "Coordinates" << endln;
		NodeIter &theNodes = theDomain.getNodes();
		Node *theNode;
		while ((theNode = theNodes()) != 0) {
			int tag = theNode->getTag();
			const Vector &crds = theNode->getCrds();
			//outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
			int l_tmp = crds.Size();
			outputFile << tag << "\t\t";
			for (int ii = 0; ii<l_tmp; ii++) {
				outputFile << crds(ii) << "\t";
			};
			for (int ii = l_tmp; ii<3; ii++) {
				outputFile << 0.0 << "\t";
			};
			outputFile << endln;
		}
		outputFile << "End coordinates" << endln << endln;

		// Print elements connectivity
		outputFile << "Elements" << endln;
		ElementIter &theElements = theDomain.getElements();
		Element *theElement;
		while ((theElement = theElements()) != 0) {
			int tag = theElement->getTag();
			// Check if element tag is inside theRange
			if (((tag<=endEle) & (tag>=startEle)) || (eleRange ==0)) {
		
				// Check type of Element with Number of Nodes
				// if 3 Nodes print the Element
				int nNode = theElement->getNumExternalNodes();
				if (nNode == 3) {
					Node **NodePtrs;
					NodePtrs = theElement->getNodePtrs();		
					Vector tagNodes(nNode);
					for (int i = 0; i < nNode; i++) {
						tagNodes(i)=NodePtrs[i]->getTag();
					}
					outputFile << tag << "\t\t";
					for (int i = 0; i < nNode; i++) {
						outputFile << tagNodes(i) << "\t";
					}
					outputFile << endln;
				}
			}
		}
		outputFile << "End elements" << endln;
	}
	// **** Quadrilateral Elements - 9 Nodes
	if (hasQuad9 == 1) {
		// Print HEADER
		outputFile << "MESH \"9NMESH\" dimension 3 ElemType Linear Nnode 9" << endln;
		outputFile << "#color 0 255 0" << endln << endln;
    
		// Print node coordinates
		outputFile << "Coordinates" << endln;
		NodeIter &theNodes = theDomain.getNodes();
		Node *theNode;
		while ((theNode = theNodes()) != 0) {
			int tag = theNode->getTag();
			const Vector &crds = theNode->getCrds();
			//outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
			int l_tmp = crds.Size();
			outputFile << tag << "\t\t";
			for (int ii = 0; ii<l_tmp; ii++) {
				outputFile << crds(ii) << "\t";
			};
			for (int ii = l_tmp; ii<3; ii++) {
				outputFile << 0.0 << "\t";
			};
			outputFile << endln;
		}
		outputFile << "End coordinates" << endln << endln;

		// Print elements connectivity
		outputFile << "Elements" << endln;
		ElementIter &theElements = theDomain.getElements();
		Element *theElement;
		while ((theElement = theElements()) != 0) {
			int tag = theElement->getTag();
			// Check if element tag is inside theRange
			if (((tag<=endEle) & (tag>=startEle)) || (eleRange ==0)) {
		
				// Check type of Element with Number of Nodes
				// if 2 Nodes print the Element
				int nNode = theElement->getNumExternalNodes();
				if (nNode == 9) {
					Node **NodePtrs;
					NodePtrs = theElement->getNodePtrs();		
					Vector tagNodes(nNode);
					for (int i = 0; i < nNode; i++) {
						tagNodes(i)=NodePtrs[i]->getTag();
					}
					outputFile << tag << "\t\t";
					for (int i = 0; i < nNode; i++) {
						outputFile << tagNodes(i) << "\t";
					}
					outputFile << endln;
				}
			}
		}
		outputFile << "End elements" << endln;
	}
	// **** Hexahedra Elements - 8 Nodes
	if (hasBrick == 1) {
		// Print HEADER
		outputFile << "MESH \"8NMESH\" dimension 3 ElemType Hexahedra Nnode 8" << endln;
		outputFile << "#color 255 0 0" << endln << endln;
    
		// Print node coordinates
		outputFile << "Coordinates" << endln;
		NodeIter &theNodes = theDomain.getNodes();
		MeshRegion *myRegion = theDomain.getRegion(0);
		Node *theNode;
		while ((theNode = theNodes()) != 0) {
			int tag = theNode->getTag();
			const Vector &crds = theNode->getCrds();
			//outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
			int l_tmp = crds.Size();
			outputFile << tag << "\t\t";
			for (int ii = 0; ii<l_tmp; ii++) {
				outputFile << crds(ii) << "\t";
			};
			for (int ii = l_tmp; ii<3; ii++) {
				outputFile << 0.0 << "\t";
			};
			outputFile << endln;
		}
		outputFile << "End coordinates" << endln << endln;

		// Print elements connectivity
		outputFile << "Elements" << endln;
		ElementIter &theElements = theDomain.getElements();
		Element *theElement;
		while ((theElement = theElements()) != 0) {
			int tag = theElement->getTag();
			// Check if element tag is inside theRange
			if (((tag<=endEle) & (tag>=startEle)) || (eleRange == 0)) {
		
				// Check type of Element with Number of Nodes
				// if 2 Nodes print the Element
				int nNode = theElement->getNumExternalNodes();
				if (nNode == 8) {
					Node **NodePtrs;
					NodePtrs = theElement->getNodePtrs();		
					Vector tagNodes(nNode);
					for (int i = 0; i < nNode; i++) {
						tagNodes(i)=NodePtrs[i]->getTag();
					}
					outputFile << tag << "\t\t";
					for (int i = 0; i < nNode; i++) {
						outputFile << tagNodes(i) << "\t";
					}
					outputFile << endln;
				}
			}
		}
		outputFile << "End elements" << endln;
	}

	outputFile.close();
	return res;
}
// Talledo End
