//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
#include <assert.h>
#include <tcl.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Vector.h>
#include <G3_Logging.h>

int
TclCommand_getEleTags(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  ElementIter &elemIter = the_domain->getElements();

  Element *elem;
  char buffer[128];
  while ((elem = elemIter()) != nullptr) {
    sprintf(buffer, "%d ", elem->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
getNumElements(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char ** const argv)
{
  assert(clientData != nullptr);

  Tcl_SetObjResult(interp, Tcl_NewIntObj(((Domain*)clientData)->getNumElements()));

  return TCL_OK;
}



int
TclCommand_addElementRayleigh(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  Domain* theTclDomain = (Domain*)clientData;

  // make sure corect number of arguments on command line
  if (argc < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: setElementRayleighFactors elementTag?  alphaM? $betaK? "
              "$betaKinit? $betaKcomm? \n";
    return TCL_ERROR;
  }

  int eleTag = 0;

  if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
    opserr << "WARNING: setElementRayleighFactors invalid eleTag: " << argv[1];
    opserr << " \n";
    return TCL_ERROR;
  }

  double alphaM, betaK, betaKinit, betaKcomm;

  if (Tcl_GetDouble(interp, argv[2], &alphaM) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "alphaM: " << argv[2] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &betaK) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "betaK: " << argv[3] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &betaKinit) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "betaKinit: " << argv[4] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &betaKcomm) != TCL_OK) {
    opserr << "WARNING : setElementRayleighFactors invalid ";
    opserr << "betaKcomm: " << argv[5] << "\n";
    return TCL_ERROR;
  }

  Element *elePtr = theTclDomain->getElement(eleTag);

  if (elePtr == nullptr)
    opserr << "WARNING : setElementRayleighFactors invalid eleTag: " << eleTag
           << " the element does not exist in the domain \n";

  if (elePtr->setRayleighDampingFactors(alphaM, betaK, betaKinit, betaKcomm) != 0) {
    opserr << "ERROR : setElementRayleighFactors: FAILED to add damping "
              "factors for element "
           << eleTag << "\n";
  }

  return TCL_OK;
}

#if 0
int
setElementRayleighDampingFactors(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 6) {
    opserr << G3_ERROR_PROMPT << "setElementRayleighDampingFactors eleTag? alphaM? betaK? "
              "betaK0? betaKc? - not enough arguments to command\n";
    return TCL_ERROR;
  }

  int eleTag;
  double alphaM, betaK, betaK0, betaKc;

  if (Tcl_GetInt(interp, argv[1], &eleTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read eleTag? \n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[2], &alphaM) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read alphaM? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &betaK) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaK? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &betaK0) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaK0? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &betaKc) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "rayleigh alphaM? betaK? betaK0? betaKc? - could not "
              "read betaKc? \n";
    return TCL_ERROR;
  }

  Element *theEle = the_domain->getElement(eleTag);
  theEle->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  return TCL_OK;
}
#endif

int
eleForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - eleForce eleTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "eleForce eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "eleForce eleTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

#if 0
  Element *theEle = the_domain->getElement(tag);
  if (theEle == 0)
    return TCL_ERROR;

  const Vector &force = theEle->getResistingForce();
#endif

  const char *myArgv[1];
  char myArgv0[8];
  strcpy(myArgv0, "forces");
  myArgv[0] = myArgv0;

  const Vector *force = domain->getElementResponse(tag, &myArgv[0], 1);
  if (force != 0) {
    int size = force->Size();

    if (dof >= 0) {

      if (size < dof)
        return TCL_ERROR;

      Tcl_SetObjResult(interp, Tcl_NewDoubleObj((*force)(dof)));

    } else {
      char buffer[128];
      for (int i = 0; i < size; ++i) {
        sprintf(buffer, "%35.20f", (*force)(i));
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }
  } else {
    opserr << G3_ERROR_PROMPT << "- failed to retrieve element force.\n";
    return TCL_ERROR;
  }
  return TCL_OK;
}

int
localForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - localForce eleTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "localForce eleTag? dof? - could not read eleTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "localForce eleTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const char *myArgv[1];
  char myArgv0[80];
  strcpy(myArgv0, "localForces");
  myArgv[0] = myArgv0;

  const Vector *force = theDomain->getElementResponse(tag, &myArgv[0], 1);

  if (force != 0) {
    int size = force->Size();

    if (dof >= 0) {

      if (size < dof)
        return TCL_ERROR;

      double value = (*force)(dof);
      Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));

    } else {
      char buffer[128];
      for (int i = 0; i < size; ++i) {
        sprintf(buffer, "%35.20f", (*force)(i));
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }
  }

  return TCL_OK;
}

int
eleDynamicalForce(ClientData clientData, Tcl_Interp *interp, int argc,
                  TCL_Char** const argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - eleForce eleTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "eleForce eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "eleForce eleTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;
  Element *theEle = theDomain->getElement(tag);
  if (theEle == nullptr)
    return TCL_ERROR;

  const Vector &force = theEle->getResistingForceIncInertia();
  int size = force.Size();

  if (dof >= 0) {

    if (size < dof)
      return TCL_ERROR;

    double value = force(dof);
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));

  } else {
    char buffer[40];
    for (int i = 0; i < size; ++i) {
      sprintf(buffer, "%35.20f", force(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
eleResponse(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char** const argv)
{
  Domain* the_domain = (Domain*)clientData; 

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - eleResponse eleTag? eleArgs...\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "eleResponse eleTag? args? - could not read eleTag? \n";
    return TCL_ERROR;
  }

  // TODO: Create element response function, remove from domain
  const Vector *data = the_domain->getElementResponse(tag, argv + 2, argc - 2);
  if (data != nullptr) {
    const int size = data->Size();
    Tcl_Obj* listPtr = Tcl_NewListObj(size, nullptr);
    for (int i = 0; i < size; ++i)
      Tcl_ListObjAppendElement(interp, listPtr, Tcl_NewDoubleObj((*data)(i)));

    Tcl_SetObjResult(interp,listPtr);	
  }
  return TCL_OK;
}

int
eleNodes(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - eleNodes eleTag?\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "eleNodes eleTag? \n";
    return TCL_ERROR;
  }

  char buffer[48];

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "eleNodes ele " << tag << " not found" << "\n";
    return TCL_ERROR;
  }
  int numTags = theElement->getNumExternalNodes();
  const ID &tags = theElement->getExternalNodes();
  for (int i = 0; i < numTags; ++i) {
    sprintf(buffer, "%d ", tags(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
eleType(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - eleType eleTag?\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "eleType eleTag? \n";
    return TCL_ERROR;
  }

  Element *theElement = the_domain->getElement(tag);

  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "eleType ele " << tag << " not found" << "\n";
    return TCL_ERROR;
  }
  const char *type = theElement->getClassType();

  Tcl_AppendResult(interp, type, NULL);

  return TCL_OK;
}

