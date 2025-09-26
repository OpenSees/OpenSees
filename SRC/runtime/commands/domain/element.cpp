//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#include <assert.h>
#include <tcl.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Vector.h>
#include <Logging.h>

namespace OpenSees {
namespace DomainCommands {

int
getEleTags(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char ** const argv)
{
  // NOTE: Maybe this can use a base class of ElementIter so we only need
  //       to work in terms of tagged object

  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  ElementIter &elemIter = the_domain->getElements();

  Tcl_Obj* result = Tcl_NewListObj(the_domain->getNumElements(), nullptr);

  Element *elem;
  while ((elem = elemIter()) != nullptr)
    Tcl_ListObjAppendElement(interp, result, Tcl_NewIntObj(elem->getTag()));

  Tcl_SetObjResult(interp, result);

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
addElementRayleigh(ClientData clientData, Tcl_Interp *interp,
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
    opserr << "WARNING : invalid ";
    opserr << "alphaM: " << argv[2] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &betaK) != TCL_OK) {
    opserr << "WARNING : invalid ";
    opserr << "betaK: " << argv[3] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &betaKinit) != TCL_OK) {
    opserr << "WARNING : invalid ";
    opserr << "betaKinit: " << argv[4] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &betaKcomm) != TCL_OK) {
    opserr << "WARNING : invalid ";
    opserr << "betaKcomm: " << argv[5] << "\n";
    return TCL_ERROR;
  }

  Element *elePtr = theTclDomain->getElement(eleTag);

  if (elePtr == nullptr)
    opserr << "WARNING : invalid eleTag: " << eleTag
           << " the element does not exist in the domain \n";

  if (elePtr->setRayleighDampingFactors(alphaM, betaK, betaKinit, betaKcomm) != 0) {
    opserr << "ERROR :: Failed to add damping "
              "factors for element "
           << eleTag << "\n";
  }

  return TCL_OK;
}


int
eleForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "want - eleForce eleTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "eleForce eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "eleForce eleTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

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
      Tcl_Obj* result = Tcl_NewListObj(size, nullptr);
      for (int i = 0; i < size; ++i)
        Tcl_ListObjAppendElement(interp, result, Tcl_NewDoubleObj((*force)(i)));

      Tcl_SetObjResult(interp, result);
    }

  } else {
    opserr << OpenSees::PromptValueError 
           << "- failed to retrieve element force."
           << OpenSees::SignalMessageEnd;
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
    opserr << OpenSees::PromptValueError 
           << "want - localForce eleTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "localForce eleTag? dof? - could not read eleTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "localForce eleTag? dof? - could not read dof? \n";
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
    opserr << OpenSees::PromptValueError << "want - eleForce eleTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "eleForce eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "eleForce eleTag? dof? - could not read dof? \n";
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
    opserr << OpenSees::PromptValueError << "want - eleResponse tag? eleArgs...\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "eleResponse tag? args? - could not read tag? \n";
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
    opserr << OpenSees::PromptValueError << "want - eleNodes eleTag?\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "eleNodes eleTag? \n";
    return TCL_ERROR;
  }

  char buffer[48];

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << OpenSees::PromptValueError << "eleNodes ele " << tag << " not found" << "\n";
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
    opserr << OpenSees::PromptValueError << "want - eleType eleTag?\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "eleType eleTag? \n";
    return TCL_ERROR;
  }

  Element *theElement = the_domain->getElement(tag);

  if (theElement == nullptr) {
    opserr << OpenSees::PromptValueError << "eleType ele " << tag << " not found" << "\n";
    return TCL_ERROR;
  }
  const char *type = theElement->getClassType();

  Tcl_AppendResult(interp, type, NULL);

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
      opserr << OpenSees::PromptValueError 
             << "getParamValue -- could not read paramTag \n";
      return TCL_ERROR;
    }

    Element *theEle = the_domain->getElement(eleTag);

    char buffer[20];
    sprintf(buffer, "%d ", theEle->getClassTag());
    Tcl_AppendResult(interp, buffer, NULL);

  } else {
    opserr << OpenSees::PromptValueError 
           << "want - getEleClassTags <eleTag?>\n" << "\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

} // namespace DomainCommands
} // namespace OpenSees
