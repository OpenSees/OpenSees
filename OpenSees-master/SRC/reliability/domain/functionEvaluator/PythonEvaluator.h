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

#ifndef PythonEvaluator_h
#define PythonEvaluator_h

#include <Domain.h>
#include <FunctionEvaluator.h>
#include <Python.h>
#include <ReliabilityDomain.h>
#include <vector>
#include <string>

class PythonEvaluator : public FunctionEvaluator {
 public:
  PythonEvaluator(ReliabilityDomain *passedReliabilityDomain,
                  Domain *passedOpenSeesDomain, const char *fileName);
  PythonEvaluator(ReliabilityDomain *passedReliabilityDomain,
                  Domain *passedOpenSeesDomain);
  ~PythonEvaluator();

  // pure virtual
  int setVariables(void);
  int setExpression(const char *expression);
  int addToExpression(const char *expression);
  double evaluateExpression(void);
  int runAnalysis(void);

  // MHS hack for reliability recorders ... set value in namespace
  int setResponseVariable(const char *label, int lsfTag, int rvTag,
                          double value);
  int setResponseVariable(const char *label, int lsfTag, double value);
  double getResponseVariable(const char *label, int lsfTag, int rvTag);
  double getResponseVariable(const char *label, int lsfTag);

 protected:

  // load module dict, return [pymodule, moduleDict]
  std::vector<PyObject*> loadModuleDict();

 private:
  ReliabilityDomain *theReliabilityDomain;
  Domain *theOpenSeesDomain;

  char *fileName;
  char *theExpression;

  double current_val;

  std::string moduleName;
};

#endif
