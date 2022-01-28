/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */

//
// Written by:
// Minjie Zhu (zhum@oregonstate.edu)
//

#include <Parameter.h>
#include <PythonEvaluator.h>
#include <RandomVariablePositionerIter.h>
#include <Vector.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>

PythonEvaluator::PythonEvaluator(
    ReliabilityDomain *passedReliabilityDomain,
    Domain *passedOpenSeesDomain, const char *passed_fileName)
    : FunctionEvaluator(),
      theReliabilityDomain(passedReliabilityDomain),
      theOpenSeesDomain(passedOpenSeesDomain) {
  theExpression = 0;
  int exprLen = strlen(passed_fileName);
  fileName = new char[exprLen + 1];
  strcpy(fileName, passed_fileName);
}

PythonEvaluator::PythonEvaluator(
    ReliabilityDomain *passedReliabilityDomain,
    Domain *passedOpenSeesDomain)
    : FunctionEvaluator(),
      theReliabilityDomain(passedReliabilityDomain),
      theOpenSeesDomain(passedOpenSeesDomain) {
  theExpression = 0;
  fileName = 0;
}

PythonEvaluator::~PythonEvaluator() {
  if (theExpression != 0) delete[] theExpression;
  if (fileName != 0) delete[] fileName;
}

int PythonEvaluator::setVariables() {
  // PyRun_SimpleString("print(dir())");

  // get module object
  PyObject *name = PyUnicode_FromString("opensees");
  PyObject *pymodule = PyImport_GetModule(name);

  if (pymodule == NULL) {
    opserr << "WARNING: module opensees is not imported\n";
    return -1;
  }

  // get module dict
  PyObject *moduleDict = PyModule_GetDict(pymodule);
  if (moduleDict == NULL) {
    opserr << "WARNING: module opensees dict is not available\n";
    return -1;
  }

  // get parameter variable
  PyObject *params = PyDict_GetItemString(moduleDict, "OpenSeesParameter");
  if (params == NULL) {
    opserr << "WARNING: variable OpenSeesParameter is not defined in "
              "module opensees\n ";
    return -1;
  }

  // clear parameter dict
  PyDict_Clear(params);

  // Set values of parameters in the Python interpreter
  int nparam = theOpenSeesDomain->getNumParameters();
  for (int i = 0; i < nparam; i++) {
    Parameter *theParam = theOpenSeesDomain->getParameterFromIndex(i);
    int paramTag = theParam->getTag();

    // now get parameter values directly
    double xval = theParam->getValue();

    // put existing parameters
    PyObject *key = PyLong_FromLong(paramTag);
    if (key == NULL) {
      opserr << "WARNING: failed to create parameter key\n";
      return -1;
    }
    PyObject *val = PyFloat_FromDouble(xval);
    if (val == NULL) {
      opserr << "WARNING: failed to create parameter value\n";
      return -1;
    }
    if (PyDict_SetItem(params, key, val) < 0) {
      opserr << "WARNING: failed to set parameter in Python\n";
      Py_DECREF(key);
      Py_DECREF(val);
      return -1;
    }
    Py_DECREF(key);
    Py_DECREF(val);
  }

  // clean up
  Py_DECREF(name);
  Py_DECREF(pymodule);

  return 0;
}

int PythonEvaluator::setExpression(const char *passedExpression) {
  if (theExpression != 0) delete[] theExpression;

  int exprLen = strlen(passedExpression);
  theExpression = new char[exprLen + 1];
  strcpy(theExpression, passedExpression);

  return 0;
}

int PythonEvaluator::addToExpression(const char *in) { return 0; }

double PythonEvaluator::evaluateExpression() {
  if (theExpression == 0) {
    opserr << "PythonEvaluator::evaluateExpression -- must set the "
              "expression before trying ";
    opserr << "to evaluate" << endln;
    return -1;
  }

  // if (Tcl_ExprDouble(theTclInterp, theExpression, &current_val) !=
  //     TCL_OK) {
  //   opserr << "PythonEvaluator::evaluateExpression -- expression \""
  //          << theExpression;
  //   opserr << "\" caused error:" << endln
  //          << Tcl_GetStringResult(theTclInterp) << endln;
  //   return -1;
  // }

  this->incrementEvaluations();
  return current_val;
}

int PythonEvaluator::runAnalysis() {
  // Let's just make a direct call since we have the pointer to OpenSees
  // domain This replaces above call to Tcl command; however, in the reset
  // command revertToStart() is also called on theTransientIntegrator --
  // MHS needs to check
  if (theOpenSeesDomain->revertToStart() != 0) {
    opserr << "ERROR PythonEvaluator -- error in resetting Domain"
           << endln;
    return -1;
  }

  // Source the code file that the user has provided
  if (fileName == 0) {
    // no source file provided, this is akin to the basic evaluator of days
    // gone by

  } else {
    // if (Tcl_Eval(theTclInterp, fileName) == TCL_ERROR) {
    //   opserr << "ERROR PythonEvaluator -- error in Tcl_Eval: "
    //          << Tcl_GetStringResult(theTclInterp) << endln;
    //   return -1;
    // }

    // make sure the parameter variables in the namespace update to reflect
    // the results of above analysis
    Parameter *theParam;

    // Set values of parameters in the Tcl interpreter
    int nparam = theOpenSeesDomain->getNumParameters();

    for (int i = 0; i < nparam; i++) {
      theParam = theOpenSeesDomain->getParameterFromIndex(i);
      if (theParam->isImplicit()) theParam->update(0.0);
    }
    this->setVariables();
  }

  return 0;
}

int PythonEvaluator::setResponseVariable(const char *label, int lsfTag,
                                         int rvTag, double value) {
  char theIndex[80];

  sprintf(theIndex, "%d,%d", lsfTag, rvTag);

  // if (Tcl_SetVar2Ex(theTclInterp, label, theIndex,
  // Tcl_NewDoubleObj(value),
  //                   TCL_LEAVE_ERR_MSG) == NULL) {
  //   opserr << "ERROR PythonEvaluator -- error in setResponseVariable for
  //   "
  //             "object with tag "
  //          << rvTag << endln;
  //   opserr << "of type " << Tcl_GetStringResult(theTclInterp) << endln;
  //   return -1;
  // }

  return 0;
}

int PythonEvaluator::setResponseVariable(const char *label, int lsfTag,
                                         double value) {
  char theIndex[80];

  sprintf(theIndex, "%d", lsfTag);

  // if (Tcl_SetVar2Ex(theTclInterp, label, theIndex,
  // Tcl_NewDoubleObj(value),
  //                   TCL_LEAVE_ERR_MSG) == NULL) {
  //   opserr << "ERROR PythonEvaluator -- error in setResponseVariable for
  //   "
  //             "object with tag "
  //          << lsfTag << endln;
  //   opserr << "of type " << Tcl_GetStringResult(theTclInterp) << endln;
  //   return -1;
  // }

  return 0;
}

double PythonEvaluator::getResponseVariable(const char *label, int lsfTag,
                                            int rvTag) {
  char theIndex[80];

  sprintf(theIndex, "%d,%d", lsfTag, rvTag);

  // Tcl_Obj *value =
  //     Tcl_GetVar2Ex(theTclInterp, label, theIndex, TCL_LEAVE_ERR_MSG);
  // if (value == NULL) {
  //   opserr << "ERROR PythonEvaluator -- error in getResponseVariable for
  //   "
  //             "object with tag "
  //          << rvTag << endln;
  //   opserr << "of type " << Tcl_GetStringResult(theTclInterp) << endln;
  //   return -1;
  // }

  double result;
  // Tcl_GetDoubleFromObj(theTclInterp, value, &result);

  return result;
}

double PythonEvaluator::getResponseVariable(const char *label,
                                            int lsfTag) {
  char theIndex[80];

  sprintf(theIndex, "%d", lsfTag);

  // Tcl_Obj *value =
  //     Tcl_GetVar2Ex(theTclInterp, label, theIndex, TCL_LEAVE_ERR_MSG);
  // if (value == NULL) {
  //   opserr << "ERROR PythonEvaluator -- error in getResponseVariable for
  //   "
  //             "object with tag "
  //          << lsfTag << endln;
  //   opserr << "of type " << Tcl_GetStringResult(theTclInterp) << endln;
  //   return -1;
  // }

  double result;
  // Tcl_GetDoubleFromObj(theTclInterp, value, &result);

  return result;
}
