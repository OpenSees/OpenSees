//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===// 
#include <tcl.h>
#include <assert.h>

#include <Node.h>
#include <Domain.h>
#include <Element.h>
#include <Logging.h>
#include <Parameter.h>
#include <Response.h>
#include <DummyStream.h>
#include <ParameterIter.h>
#include <BasicAnalysisBuilder.h>
#include <LoadPattern.h>
#include <Pressure_Constraint.h>
#include <InitialStateParameter.h>
#include <ElementStateParameter.h>


int
sensNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain *)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << OpenSees::PromptValueError << "want - sensNodeDisp nodeTag? dof? paramTag?\n";
    return TCL_ERROR;
  }

  int tag, dof, paramTag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr
        << OpenSees::PromptValueError << "nodeDisp nodeTag? dof? paramTag?- could not read nodeTag? ";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "nodeDisp nodeTag? dof? paramTag?- could not read dof? ";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &paramTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "nodeDisp paramTag? dof? paramTag?- could not read "
              "paramTag? ";
    return TCL_ERROR;
  }

  //
  //
  Node *theNode = theDomain->getNode(tag);
  if (theNode == nullptr) {
    opserr << "sensNodeDisp: node " << tag << " not found" << "\n";
    return TCL_ERROR;
  }

  Parameter *theParam = theDomain->getParameter(paramTag);
  if (theParam == nullptr) {
    opserr << "sensNodeDisp: parameter " << paramTag << " not found" << "\n";
    return TCL_ERROR;
  }

  int gradIndex = theParam->getGradIndex();

  double value = theNode->getDispSensitivity(dof, gradIndex);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));

  return TCL_OK;
}

int
sensNodeVel(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain *)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << OpenSees::PromptValueError << "want - sensNodeVel nodeTag? dof? paramTag?\n";
    return TCL_ERROR;
  }

  int tag, dof, paramTag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensNodeVel nodeTag? dof? paramTag? - could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensNodeVel nodeTag? dof? paramTag? - could not read "
              "dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &paramTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensNodeVel nodeTag? dof? paramTag? - could not read "
              "paramTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = theDomain->getNode(tag);
  if (theNode == nullptr) {
    opserr << "sensNodeVel: node " << tag << " not found" << "\n";
    return TCL_ERROR;
  }

  Parameter *theParam = theDomain->getParameter(paramTag);
  if (theParam == nullptr) {
    opserr << "sensNodeVel: parameter " << paramTag << " not found" << "\n";
    return TCL_ERROR;
  }

  int gradIndex = theParam->getGradIndex();

  double value = theNode->getVelSensitivity(dof, gradIndex);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));

  return TCL_OK;
}

int
sensNodeAccel(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain *)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr << OpenSees::PromptValueError << "want - sensNodeAccel nodeTag? dof? paramTag?\n";
    return TCL_ERROR;
  }

  int tag, dof, paramTag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensNodeAccel nodeTag? dof? paramTag? - could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensNodeAccel nodeTag? dof? paramTag? - could not read "
              "dof? \n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &paramTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sendNodeAccel nodeTag? dof? paramTag? - could not read "
              "paramTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = theDomain->getNode(tag);
  if (theNode == nullptr) {
    opserr << "sensNodeAccel: node " << tag << " not found" << "\n";
    return TCL_ERROR;
  }

  Parameter *theParam = theDomain->getParameter(paramTag);
  if (theParam == nullptr) {
    opserr << "sensNodeAccel: parameter " << paramTag << " not found" << "\n";
    return TCL_ERROR;
  }

  int gradIndex = theParam->getGradIndex();

  double value = theNode->getAccSensitivity(dof, gradIndex);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));

  return TCL_OK;
}

int
sensNodePressure(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain *)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 3) {
    opserr << OpenSees::PromptValueError << "want - sensNodePressure nodeTag? paramTag?\n";
    return TCL_ERROR;
  }

  int tag, paramTag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensNodePressure nodeTag? paramTag?- could not read "
              "nodeTag? ";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &paramTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensNodePressure paramTag? paramTag?- could not read "
              "paramTag? ";
    return TCL_ERROR;
  }

  double dp = 0.0;
  Pressure_Constraint *thePC = theDomain->getPressure_Constraint(tag);
  if (thePC != nullptr) {
    // int ptag = thePC->getPressureNode();
    // Node* pNode = theDomain->getNode(ptag);
    Node *pNode = thePC->getPressureNode();
    if (pNode != nullptr) {

      Parameter *theParam = theDomain->getParameter(paramTag);
      if (theParam == nullptr) {
        opserr << "sensNodePressure: parameter " << paramTag << " not found"
               << "\n";
        return TCL_ERROR;
      }

      int gradIndex = theParam->getGradIndex();
      dp = pNode->getVelSensitivity(1, gradIndex);
    }
  }

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(dp));

  return TCL_OK;
}

int
sensSectionForce(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv)
{
#if 1 // def _RELIABILITY
  assert(clientData != nullptr);
  Domain *theDomain = (Domain *)clientData;

  // make sure at least one other argument to contain type of system
  if (argc < 4) {
    opserr
        << OpenSees::PromptValueError << "want - sensSectionForce eleTag? <secNum?> dof? paramTag?\n";
    return TCL_ERROR;
  }

  int tag, dof, paramTag;
  int secNum = 0;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensSectionForce eleTag? secNum? dof? paramTag?- could "
              "not read eleTag? \n";
    return TCL_ERROR;
  }

  // Make this work for zeroLengthSection too
  int currentArg = 2;
  if (argc > 4) {
    if (Tcl_GetInt(interp, argv[currentArg++], &secNum) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "sensSectionForce eleTag? secNum? dof? paramTag?- "
                "could not read secNum? \n";
      return TCL_ERROR;
    }
  }
  if (Tcl_GetInt(interp, argv[currentArg++], &dof) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensSectionForce eleTag? secNum? dof? paramTag?- could "
              "not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[currentArg++], &paramTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "sensSectionForce eleTag? secNum? dof? paramTag?- could "
              "not read paramTag? \n";
    return TCL_ERROR;
  }

  ParameterIter &pIter = theDomain->getParameters();
  Parameter *theParam;
  while ((theParam = pIter()) != 0)
    theParam->activate(false);

  theParam = theDomain->getParameter(paramTag);
  int gradIndex = theParam->getGradIndex();
  theParam->activate(true);

  Element *theElement = theDomain->getElement(tag);
  if (theElement == 0) {
    opserr << OpenSees::PromptValueError << "sensSectionForce element with tag " << tag
           << " not found in domain \n";
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
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponseSensitivity(gradIndex);
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(theVec(dof - 1)));

  theParam->activate(false);

  delete theResponse;
#endif
  return TCL_OK;
}



