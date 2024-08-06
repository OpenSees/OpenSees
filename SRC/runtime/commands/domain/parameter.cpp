//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <OPS_Stream.h>
#include <G3_Logging.h>
#include <Domain.h>
#include <MovableObject.h>
#include <Element.h>
#include <Node.h>
#include <NodeData.h>

#include <BasicModelBuilder.h>

#include <Parameter.h>
#include <ParameterIter.h>
#include <ElementParameter.h>
#include <ElementStateParameter.h>

#include <NodeResponseParameter.h>
#include <LoadFactorParameter.h>
#include <LoadPattern.h>

#ifdef _RELIABILITY
#include <RandomVariable.h>
#include <RVParameter.h>
#include <ReliabilityDomain.h>

extern ReliabilityDomain *theReliabilityDomain;

#endif

int
TclCommand_parameter(ClientData clientData, Tcl_Interp *interp, int argc,
                     TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData; 


  // check at least two arguments so don't segemnt fault on strcmp
  if (argc < 2) {
    opserr << "WARNING need to specify a parameter tag\n";
    opserr << "Want: parameter tag <specific parameter args> .. see manual for "
              "valid parameter types and arguments\n";
    return TCL_ERROR;
  }

  // Figure out which parameter we are dealing with
  int paramTag;
  if (Tcl_GetInt(interp, argv[1], &paramTag) != TCL_OK) {
    return TCL_ERROR;
  }

  Parameter *theParameter = domain->getParameter(paramTag);
  int eleTag = -1;
  bool isele = false;

  // First, check special case of a blank parameter
  if (argc == 2 && strcmp(argv[0], "parameter") == 0) {
    Parameter *newParameter = new Parameter(paramTag, 0, 0, 0);

    domain->addParameter(newParameter);

    char buffer[40];
    sprintf(buffer, "%d", paramTag);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
  }

  // First, check special case of a blank parameter
  if (argc == 3 && strcmp(argv[0], "parameter") == 0) {
    Parameter *newParameter = new Parameter(paramTag, 0, 0, 0);

    double value;
    if (Tcl_GetDouble(interp, argv[2], &value) != TCL_OK)
      return TCL_ERROR;

    newParameter->setValue(value);

    domain->addParameter(newParameter);

    char buffer[40];
    sprintf(buffer, "%d", paramTag);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
  }

  if (argc >= 6 && strcmp(argv[0], "parameter") == 0 &&
      strcmp(argv[2], "node") == 0 && strcmp(argv[4], "disp") == 0) {

    int nodeTag;
    if (Tcl_GetInt(interp, argv[3], &nodeTag) != TCL_OK) {
      return TCL_ERROR;
    }
    Node *theNode = domain->getNode(nodeTag);

    int dof;
    if (Tcl_GetInt(interp, argv[5], &dof) != TCL_OK) {
      return TCL_ERROR;
    }

    Parameter *newParameter =
        new NodeResponseParameter(paramTag, theNode, NodeData::Disp, dof);

    domain->addParameter(newParameter);

    char buffer[40];
    sprintf(buffer, "%d", paramTag);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
  }

  if (argc >= 5 && strcmp(argv[0], "parameter") == 0 &&
      strcmp(argv[2], "pattern") == 0 && strcmp(argv[4], "lambda") == 0) {

    int patternTag;
    if (Tcl_GetInt(interp, argv[3], &patternTag) != TCL_OK) {
      return TCL_ERROR;
    }
    LoadPattern *thePattern = domain->getLoadPattern(patternTag);

    Parameter *newParameter = new LoadFactorParameter(paramTag, thePattern);

    domain->addParameter(newParameter);

    char buffer[40];
    sprintf(buffer, "%d", paramTag);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

    return TCL_OK;
  }


  // Now handle the parameter according to which command is invoked
  if (strcmp(argv[0], "parameter") == 0 ||
      strcmp(argv[0], "addToParameter") == 0) {
    // RandomVariable *theRV = 0;
    void *theRV = 0;

    MovableObject *theObject = nullptr;

    if (strstr(argv[2], "randomVariable") != 0) {
#ifdef _RELIABILITY
      int rvTag;
      if (Tcl_GetInt(interp, argv[3], &rvTag) != TCL_OK) {
        return TCL_ERROR;
      }

      if (theReliabilityDomain == 0) {
        opserr << "ERROR parameter " << paramTag
               << " -- reliability domain has not been created" << endln;
      }

      theRV = theReliabilityDomain->getRandomVariablePtr(rvTag);
      if (theRV == 0) {
        opserr << "ERROR parameter " << paramTag
               << " -- random variable with tag " << rvTag << " not defined"
               << endln;
        return TCL_ERROR;
      }
#endif
    }

    int argStart = (theRV) ? 4 : 2;

    if (argc > argStart && strstr(argv[argStart], "element") != 0) {

      if (argc < 4) {
        opserr << "WARNING parameter -- insufficient number of arguments for "
                  "parameter with tag "
               << paramTag << '\n';
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[argStart + 1], &eleTag) != TCL_OK) {
        opserr << "WARNING parameter -- invalid element tag\n";
        return TCL_ERROR;
      }
      isele = true;

      // Retrieve element from domain
      //  FMK theObject = (MovableObject *) theTclDomain->getElement(eleTag);
      theObject = static_cast<MovableObject *>(domain->getElement(eleTag));

      argStart = (theRV) ? 6 : 4;

    } else if (argc > argStart && strstr(argv[argStart], "node") != 0) {
      if (argc < 4) {
        opserr << "WARNING parameter -- insufficient number of arguments for "
                  "parameter with tag "
               << paramTag << '\n';
        return TCL_ERROR;
      }

      int nodeTag;
      if (Tcl_GetInt(interp, argv[argStart + 1], &nodeTag) != TCL_OK) {
        opserr << "WARNING parameter -- invalid node tag\n";
        return TCL_ERROR;
      }

      // Retrieve node from domain
      theObject = static_cast<MovableObject *>(domain->getNode(nodeTag));

      argStart = (theRV) ? 6 : 4;

    } else if (argc > argStart && strstr(argv[argStart], "loadPattern") != 0) {


      if (argc < 4) {
        opserr << "WARNING parameter -- insufficient number of arguments for "
                  "parameter with tag "
               << paramTag << '\n';
        return TCL_ERROR;
      }

      int loadTag;
      if (Tcl_GetInt(interp, argv[argStart + 1], &loadTag) != TCL_OK) {
        opserr << "WARNING parameter -- invalid load pattern tag\n";
        return TCL_ERROR;
      }

      // Retrieve element from domain
      theObject = static_cast<MovableObject *>(domain->getLoadPattern(loadTag));

      argStart = (theRV) ? 6 : 4;


    } else if (argc > argStart) {
      opserr << "WARNING - unable to assign parameter to object of type "
             << argv[2] << '\n';
      return TCL_ERROR;
    }

    ///////////////////////////////////

    // Create new parameter
    if (strcmp(argv[0], "parameter") == 0) {

      if (theParameter != 0) {
        opserr << "WARNING parameter -- parameter with tag " << paramTag
               << " already exists in domain\n";
        return TCL_ERROR;
      }

      Parameter *newParameter;
      if (argc > argStart) {
        if (isele == false) {
          newParameter =
              new Parameter(paramTag, theObject, (const char **)&argv[argStart],
                            argc - argStart);

        } else {
          newParameter = new ElementParameter(paramTag, eleTag,
                                              (const char **)&argv[argStart],
                                              argc - argStart);
        }
      } else
        newParameter = new Parameter(paramTag, 0, 0, 0);

      if (theRV != 0) {
#ifdef _RELIABILITY
        RVParameter *newRVParameter =
            new RVParameter(paramTag, theRV, newParameter);
        domain->addParameter(newRVParameter);
#else
        opserr << "ERROR: Reliability not compiled in\n";
#endif
      } else {
        domain->addParameter(newParameter);
      }

      char buffer[40];
      sprintf(buffer, "%d", paramTag);
      Tcl_SetResult(interp, buffer, TCL_VOLATILE);
    }

    // Add to an existing parameter
    else if (strcmp(argv[0], "addToParameter") == 0) {

      if (theParameter == 0) {
        opserr << "WARNING addToParameter -- parameter with tag " << paramTag
               << " not found in domain\n";
        return TCL_ERROR;

      } else {
        if (isele == false)
          theParameter->addComponent(theObject, (const char **)&argv[argStart],
                                     argc - argStart);
        else {
          theObject = static_cast<MovableObject *>(domain->getElement(eleTag));
          theParameter->addComponent(theObject, (const char **)&argv[argStart],
                                     argc - argStart);
          // Sorry, Frank, had to change this -- MHS
          // theParameter->addComponent(eleTag, (const char **)&argv[argStart],
          // argc-argStart);
        }
      }
    }

    return TCL_OK;
  }

  else if (strcmp(argv[0], "updateParameter") == 0) {

    // Cannot update a parameter that is not present
    if (theParameter == 0) {
      opserr << "WARNING updateParameter -- parameter with tag " << paramTag
             << " not found in domain\n";
      //  return TCL_ERROR;
    }

    double newValue;
    if (Tcl_GetDouble(interp, argv[2], &newValue) != TCL_OK) {
      opserr << "WARNING updateParameter -- invalid parameter value\n";
      return TCL_ERROR;
    }

    //    theParameter->update(newValue);
    domain->updateParameter(paramTag, newValue);
  }

  return TCL_OK;
}

int
getParamTags(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* the_domain = (Domain*)clientData; 

  Parameter *theParam;
  ParameterIter &paramIter = the_domain->getParameters();

  char buffer[20];

  while ((theParam = paramIter()) != nullptr) {
    sprintf(buffer, "%d ", theParam->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
getParamValue(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* the_domain = (Domain*)clientData; 

  if (argc < 2) {
    opserr << "Insufficient arguments to getParamValue" << endln;
    return TCL_ERROR;
  }

  int paramTag;

  if (Tcl_GetInt(interp, argv[1], &paramTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "getParamValue -- could not read paramTag \n";
    return TCL_ERROR;
  }

  Parameter *theEle = the_domain->getParameter(paramTag);

  char buffer[40];

  sprintf(buffer, "%35.20f", theEle->getValue());
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


int
setParameter(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  Domain *theDomain = (Domain*)clientData;

  int argLoc = 1;
  double newValue = 0.0;
  ID eleIDs(0, 32);
  int numEle = 0;
  int flag = 0;

  if (strstr(argv[argLoc], "-val") != 0) {
    if (Tcl_GetDouble(interp, argv[argLoc + 1], &newValue) != TCL_OK) {
      opserr << "WARNING setParameter: invalid parameter value\n";
      return TCL_ERROR;
    }
  } else {
    opserr << "WARNING setParameter:  -val not found " << endln;
    return TCL_ERROR;
  }

  argLoc += 2;

  if (strstr(argv[argLoc], "-ele") != 0) {

    if ((strcmp(argv[argLoc], "-ele") == 0) ||
        (strcmp(argv[argLoc], "-eles") == 0) ||
        (strcmp(argv[argLoc], "-element") == 0)) {

      //
      // read in a list of ele until end of command or other flag
      //

      argLoc++;
      int eleTag;

      while (argLoc < argc &&
             Tcl_GetInt(interp, argv[argLoc], &eleTag) == TCL_OK) {
        eleIDs[numEle] = eleTag;
        numEle++;
        argLoc++;
      }

      if (numEle > 0)
        flag = 1;

    } else if (strcmp(argv[argLoc], "-eleRange") == 0) {

      flag = 2;

      // ensure no segmentation fault if user messes up
      if (argc < argLoc + 3) {
        opserr << "WARNING recorder Element .. -eleRange start? end?  .. - no "
                  "ele tags specified\n";
        return TCL_ERROR;
      }

      //
      // read in start and end tags of two elements & add set [start,end]
      //

      int start, end;
      if (Tcl_GetInt(interp, argv[argLoc + 1], &start) != TCL_OK) {
        opserr
            << "WARNING recorder Element -eleRange start? end? - invalid start "
            << argv[argLoc + 1] << endln;
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[argLoc + 2], &end) != TCL_OK) {
        opserr
            << "WARNING recorder Element -eleRange start? end? - invalid end "
            << argv[argLoc + 2] << endln;
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

    ElementStateParameter theParameter(newValue, &argv[argLoc], argc - argLoc,
                                       flag, &eleIDs);

    theDomain->addParameter(&theParameter);
  }

  return TCL_OK;
}
