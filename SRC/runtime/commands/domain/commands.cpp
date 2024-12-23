//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified.
//
#include <tcl.h>
#include <G3_Logging.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <classTags.h>
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <assert.h>
#include <set>
#include <vector>
#include <algorithm>
//
#include "commands.h"
// Domain
#include <Domain.h>
#include <DOF_Group.h>
#include <Matrix.h>
#include <Node.h>
#include <Element.h>
#include <ElementIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
//
#include <Parameter.h>
#include <ParameterIter.h>
#include <InitialStateParameter.h>
#include <ElementStateParameter.h>
#include <Pressure_Constraint.h>
// Analysis
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>
#include <AnalysisModel.h>
#include <EquiSolnAlgo.h>
#include <Integrator.h>
#include <StaticIntegrator.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
// Other
#include <Recorder.h>
#include <DummyStream.h>
#include <Information.h>
#include <Response.h>
#include <packages.h>
//
// Global variables
//
class ModelBuilder;
ModelBuilder          *theBuilder         = nullptr;
VariableTimeStepDirectIntegrationAnalysis
                      *theVariableTimeStepTransientAnalysis = nullptr;
//
// Forward declarations
//
const char *getInterpPWD(Tcl_Interp *interp);
extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

Tcl_CmdProc TclCommand_record;
Tcl_CmdProc TclCommand_setLoadConst;
Tcl_CmdProc TclCommand_setCreep;


// TODO: reimplement defaultUnits and setParameter
// int defaultUnits(ClientData, Tcl_Interp *, int, TCL_Char ** const argv);
// int setParameter(ClientData, Tcl_Interp *, int, TCL_Char **);
int
G3_AddTclDomainCommands(Tcl_Interp *interp, Domain* the_domain)
{

  ClientData domain = (ClientData)the_domain;

  Tcl_CreateCommand(interp, "loadConst",           &TclCommand_setLoadConst,  domain, nullptr);
  Tcl_CreateCommand(interp, "recorder",            &TclAddRecorder,  domain,  nullptr);
  Tcl_CreateCommand(interp, "region",              &TclCommand_addMeshRegion, domain, nullptr);

  Tcl_CreateCommand(interp, "printGID",            &printModelGID, domain, nullptr);

  Tcl_CreateCommand(interp, "setTime",             &TclCommand_setTime,  domain, nullptr);
  Tcl_CreateCommand(interp, "getTime",             &TclCommand_getTime,  domain, nullptr);
  Tcl_CreateCommand(interp, "setCreep",            &TclCommand_setCreep, nullptr, nullptr);

  // DAMPING
  Tcl_CreateCommand(interp, "rayleigh",            &rayleighDamping, domain, nullptr);

  Tcl_CreateCommand(interp, "setElementRayleighDampingFactors", &TclCommand_addElementRayleigh, domain, nullptr);
  Tcl_CreateCommand(interp, "setElementRayleighFactors",        &TclCommand_addElementRayleigh, domain, nullptr);
  Tcl_CreateCommand(interp, "getLoadFactor",       &getLoadFactor, domain, nullptr);
  Tcl_CreateCommand(interp, "localForce",          &localForce,    domain, nullptr);
  Tcl_CreateCommand(interp, "eleType",             &eleType,       domain, nullptr);
  Tcl_CreateCommand(interp, "eleNodes",            &eleNodes,            domain, nullptr);
  Tcl_CreateCommand(interp, "getEleTags",          &TclCommand_getEleTags, domain, nullptr);
  Tcl_CreateCommand(interp, "basicDeformation",    &basicDeformation,    domain, nullptr);
  Tcl_CreateCommand(interp, "basicForce",          &basicForce,          domain, nullptr);
  Tcl_CreateCommand(interp, "basicStiffness",      &basicStiffness,      domain, nullptr);

  Tcl_CreateCommand(interp, "eleForce",            &eleForce,            domain, nullptr);
  Tcl_CreateCommand(interp, "eleResponse",         &eleResponse,         domain, nullptr);
  Tcl_CreateCommand(interp, "eleDynamicalForce",   &eleDynamicalForce,   domain, nullptr);

  Tcl_CreateCommand(interp, "nodeDOFs",            &nodeDOFs,            domain, nullptr);
  Tcl_CreateCommand(interp, "nodeCoord",           &nodeCoord,           domain, nullptr);
  Tcl_CreateCommand(interp, "nodeMass",            &nodeMass,            domain, nullptr);
  Tcl_CreateCommand(interp, "nodeVel",             &nodeVel,             domain, nullptr);
  Tcl_CreateCommand(interp, "nodeDisp",            &nodeDisp,            domain, nullptr);
  Tcl_CreateCommand(interp, "nodeAccel",           &nodeAccel,           domain, nullptr);
  Tcl_CreateCommand(interp, "nodeResponse",        &nodeResponse,        domain, nullptr);
  Tcl_CreateCommand(interp, "nodePressure",        &nodePressure,        domain, nullptr);
  Tcl_CreateCommand(interp, "nodeBounds",          &nodeBounds,          domain, nullptr);
  Tcl_CreateCommand(interp, "findNodeWithID",      &findID,              domain, nullptr);
  Tcl_CreateCommand(interp, "nodeUnbalance",       &nodeUnbalance,       domain, nullptr);
  Tcl_CreateCommand(interp, "nodeEigenvector",     &nodeEigenvector,     domain, nullptr);

  Tcl_CreateCommand(interp, "nodeReaction",        &nodeReaction,            domain, nullptr);
  Tcl_CreateCommand(interp, "reactions",           &calculateNodalReactions, domain, nullptr);

  Tcl_CreateCommand(interp, "setNodeVel",          &setNodeVel,              domain, nullptr);
  Tcl_CreateCommand(interp, "setNodeDisp",         &setNodeDisp,             domain, nullptr);
  Tcl_CreateCommand(interp, "setNodeAccel",        &setNodeAccel,            domain, nullptr);
  Tcl_CreateCommand(interp, "setNodeCoord",        &setNodeCoord,            domain, nullptr);

  Tcl_CreateCommand(interp, "getEleTags",          &getEleTags,              domain, nullptr);
  Tcl_CreateCommand(interp, "getNodeTags",         &getNodeTags,             domain, nullptr);

  Tcl_CreateCommand(interp, "getParamTags",        &getParamTags,            domain, nullptr);
  Tcl_CreateCommand(interp, "getParamValue",       &getParamValue,           domain, nullptr);
  Tcl_CreateCommand(interp, "parameter",           &TclCommand_parameter,    domain, nullptr);
  Tcl_CreateCommand(interp, "addToParameter",      &TclCommand_parameter,    domain, nullptr);
  Tcl_CreateCommand(interp, "updateParameter",     &TclCommand_parameter,    domain, nullptr);

  Tcl_CreateObjCommand(interp, "fixedNodes",          &fixedNodes,          domain, nullptr);
  Tcl_CreateObjCommand(interp, "fixedDOFs",           &fixedDOFs,           domain, nullptr);
  Tcl_CreateObjCommand(interp, "constrainedNodes",    &constrainedNodes,    domain, nullptr);
  Tcl_CreateObjCommand(interp, "constrainedDOFs",     &constrainedDOFs,     domain, nullptr);
  Tcl_CreateObjCommand(interp, "domainChange",        &domainChange,        domain, nullptr);
  Tcl_CreateObjCommand(interp, "remove",              &removeObject,        domain, nullptr);
  Tcl_CreateCommand(interp,    "retainedNodes",       &retainedNodes,       domain, nullptr);
  Tcl_CreateCommand(interp,    "retainedDOFs",        &retainedDOFs,        domain, nullptr);

  Tcl_CreateCommand(interp, "getNumElements",      &getNumElements,      domain, nullptr);
  Tcl_CreateCommand(interp, "getEleClassTags",     &getEleClassTags,     domain, nullptr);
  Tcl_CreateCommand(interp, "getEleLoadTags",      &getEleLoadTags,      domain, nullptr);
  Tcl_CreateCommand(interp, "getEleLoadData",      &getEleLoadData,      domain, nullptr);
  Tcl_CreateCommand(interp, "getEleLoadClassTags", &getEleLoadClassTags, domain, nullptr);


  Tcl_CreateCommand(interp, "sectionForce",        &sectionForce,        domain, nullptr);
  Tcl_CreateCommand(interp, "sectionTag",          &sectionTag,          domain, nullptr);
  Tcl_CreateCommand(interp, "sectionDisplacement", &sectionDisplacement, domain, nullptr);
  Tcl_CreateCommand(interp, "sectionDeformation",  &sectionDeformation,  domain, nullptr);
  Tcl_CreateCommand(interp, "sectionStiffness",    &sectionStiffness,    domain, nullptr);
  Tcl_CreateCommand(interp, "sectionFlexibility",  &sectionFlexibility,  domain, nullptr);
  Tcl_CreateCommand(interp, "sectionLocation",     &sectionLocation,     domain, nullptr);
  Tcl_CreateCommand(interp, "sectionWeight",       &sectionWeight,       domain, nullptr);

  Tcl_CreateCommand(interp, "recorderValue",       &OPS_recorderValue,   domain, nullptr);
  Tcl_CreateCommand(interp, "record",              &TclCommand_record,   domain, nullptr);

  Tcl_CreateCommand(interp, "updateElementDomain", &updateElementDomain, nullptr, nullptr);

  Tcl_CreateCommand(interp, "InitialStateAnalysis", &InitialStateAnalysis, nullptr, nullptr);


//   TODO: cmp, moved definition to packages/optimization; need to link in optionally
//   Tcl_CreateCommand(interp, "setParameter", &setParameter, nullptr, nullptr);

  // Tcl_CreateCommand(interp, "sdfResponse",      &sdfResponse, nullptr, nullptr);
  // Tcl_CreateCommand(interp, "database", &addDatabase, nullptr, nullptr);

  // wipeAnalysis(0, interp, 0, 0);
  return TCL_OK;
}



int
getLoadFactor(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData; 

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "no load pattern supplied -- getLoadFactor\n";
    return TCL_ERROR;
  }

  int pattern;
  if (Tcl_GetInt(interp, argv[1], &pattern) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "reading load pattern tag -- getLoadFactor\n";
    return TCL_ERROR;
  }

  LoadPattern *the_pattern = domain->getLoadPattern(pattern);
  if (the_pattern == nullptr) {
    opserr << G3_ERROR_PROMPT << "load pattern with tag " << pattern
           << " not found in domain -- getLoadFactor\n";
    return TCL_ERROR;
  }

  double factor = the_pattern->getLoadFactor();
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(factor));

  return TCL_OK;
}



// added by C.McGann, U.Washington
int
InitialStateAnalysis(ClientData clientData, Tcl_Interp *interp, int argc,
                     TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING: Incorrect number of arguments for InitialStateAnalysis "
              "command"
           << endln;
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "on") == 0) {
    opserr << "InitialStateAnalysis ON" << endln;

    // set global variable to true
    // FMK changes for parallel:
    // ops_InitialStateAnalysis = true;

    Parameter *theP = new InitialStateParameter(true);
    the_domain->addParameter(theP);
    delete theP;

    return TCL_OK;

  } else if (strcmp(argv[1], "off") == 0) {
    opserr << "InitialStateAnalysis OFF" << endln;

    // call revert to start to zero the displacements
    the_domain->revertToStart();

    // set global variable to false
    // FMK changes for parallel
    // ops_InitialStateAnalysis = false;
    Parameter *theP = new InitialStateParameter(false);
    the_domain->addParameter(theP);
    delete theP;

    return TCL_OK;

  } else {
    opserr << "WARNING: Incorrect arguments - want InitialStateAnalysis on, or "
              "InitialStateAnalysis off"
           << endln;

    return TCL_ERROR;
  }
}

int
rayleighDamping(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{

  if (argc < 3) {
    opserr << G3_ERROR_PROMPT
           << "rayleigh alphaM? betaK? betaK0? betaKc? - not enough "
              "arguments to command\n";
    return TCL_ERROR;
  }

  double alphaM, betaK, betaK0=0.0, betaKc=0.0;
  if (Tcl_GetDouble(interp, argv[1], &alphaM) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read alphaM? \n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[2], &betaK) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaK? \n";
    return TCL_ERROR;
  }

  if (argc > 3 && Tcl_GetDouble(interp, argv[3], &betaK0) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaK0? \n";
    return TCL_ERROR;
  }

  if (argc > 4 && Tcl_GetDouble(interp, argv[4], &betaKc) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaKc? \n";
    return TCL_ERROR;
  }

  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;
  the_domain->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  return TCL_OK;
}


int
getEleClassTags(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc == 1) {
    Element *theEle;
    ElementIter &eleIter = the_domain->getElements();

    char buffer[20];

    while ((theEle = eleIter()) != nullptr) {
      sprintf(buffer, "%d ", theEle->getClassTag());
      Tcl_AppendResult(interp, buffer, NULL);
    }
  } else if (argc == 2) {
    int eleTag;

    if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "getParamValue -- could not read paramTag \n";
      return TCL_ERROR;
    }

    Element *theEle = the_domain->getElement(eleTag);

    char buffer[20];
    sprintf(buffer, "%d ", theEle->getClassTag());
    Tcl_AppendResult(interp, buffer, NULL);

  } else {
    opserr << G3_ERROR_PROMPT << "want - getEleClassTags <eleTag?>\n" << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleLoadClassTags(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = the_domain->getLoadPatterns();

    char buffer[20];

    while ((thePattern = thePatterns()) != nullptr) {
      ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != nullptr) {
        sprintf(buffer, "%d ", theLoad->getClassTag());
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }

  } else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "getEleLoadClassTags -- could not read patternTag\n";
      return TCL_ERROR;
    }

    LoadPattern *thePattern = the_domain->getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << G3_ERROR_PROMPT << "load pattern with tag " << patternTag
             << " not found in domain -- getEleLoadClassTags\n";
      return TCL_ERROR;
    }

    ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
    ElementalLoad *theLoad;

    char buffer[20];

    while ((theLoad = theEleLoads()) != nullptr) {
      sprintf(buffer, "%d ", theLoad->getClassTag());
      Tcl_AppendResult(interp, buffer, NULL);
    }

  } else {
    opserr << G3_ERROR_PROMPT << "want - getEleLoadClassTags <patternTag?>\n" << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleLoadTags(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = the_domain->getLoadPatterns();

    char buffer[20];

    while ((thePattern = thePatterns()) != nullptr) {
      ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != nullptr) {
        sprintf(buffer, "%d ", theLoad->getElementTag());
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }

  } else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "getEleLoadTags -- could not read patternTag \n";
      return TCL_ERROR;
    }

    LoadPattern *thePattern = the_domain->getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << G3_ERROR_PROMPT << "load pattern with tag " << patternTag
             << " not found in domain -- getEleLoadTags\n";
      return TCL_ERROR;
    }

    ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
    ElementalLoad *theLoad;

    char buffer[20];

    while ((theLoad = theEleLoads()) != nullptr) {
      sprintf(buffer, "%d ", theLoad->getElementTag());
      Tcl_AppendResult(interp, buffer, NULL);
    }

  } else {
    opserr << G3_ERROR_PROMPT << "want - getEleLoadTags <patternTag?>\n" << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleLoadData(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = the_domain->getLoadPatterns();

    char buffer[40];
    int typeEL;

    while ((thePattern = thePatterns()) != nullptr) {
      ElementalLoadIter &theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != nullptr) {
        const Vector &eleLoadData = theLoad->getData(typeEL, 1.0);

        int eleLoadDataSize = eleLoadData.Size();
        opserr << "eleLoadDataSize: " << eleLoadDataSize << "\n";
        for (int i = 0; i < eleLoadDataSize; ++i) {
          sprintf(buffer, "%35.20f ", eleLoadData(i));
          Tcl_AppendResult(interp, buffer, NULL);
        }
      }
    }

  } else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "getEleLoadData -- could not read patternTag \n";
      return TCL_ERROR;
    }

    LoadPattern *thePattern = the_domain->getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << G3_ERROR_PROMPT << "load pattern with tag " << patternTag
             << " not found in domain -- getEleLoadData\n";
      return TCL_ERROR;
    }

    ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
    ElementalLoad *theLoad;

    int typeEL;
    char buffer[40];

    while ((theLoad = theEleLoads()) != nullptr) {
      const Vector &eleLoadData = theLoad->getData(typeEL, 1.0);

      int eleLoadDataSize = eleLoadData.Size();
      for (int i = 0; i < eleLoadDataSize; ++i) {
        sprintf(buffer, "%35.20f ", eleLoadData(i));
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }

  } else {
    opserr << G3_ERROR_PROMPT << "want - getEleLoadTags <patternTag?>\n" << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleTags(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  // NOTE: Maybe this can use a base class of ElementIter so we only need
  //       to work in terms of tagged object
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  Element *theEle;
  ElementIter &eleIter = the_domain->getElements();

  char buffer[20];

  while ((theEle = eleIter()) != nullptr) {
    sprintf(buffer, "%d ", theEle->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

