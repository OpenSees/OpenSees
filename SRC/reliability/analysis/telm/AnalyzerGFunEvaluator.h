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
                                                                        
// $Revision: 1.3 $
// $Date: 2010-09-13 21:39:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/AnalyzerGFunEvaluator.h,v $

//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef AnalyzerGFunEvaluator_h
#define AnalyzerGFunEvaluator_h

#include <Analyzer.h>
#include <FunctionEvaluator.h>
#include <ReliabilityDomain.h>
#include <Domain.h>
#include <TaggedObjectStorage.h>
#include <PerformanceFunctionCoeff.h>
#include <PerformanceFunctionCoefficientIter.h>
#include <Vector.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>
#include <tcl.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;


#include <fstream>
using std::ofstream;


class AnalyzerGFunEvaluator : public FunctionEvaluator
{

public:
  AnalyzerGFunEvaluator(Tcl_Interp *passedTclInterp,
			ReliabilityDomain *passedReliabilityDomain,
			Domain* passedDomain,
			Analyzer *passedAnalyzer);
  ~AnalyzerGFunEvaluator();
  
  int evaluateG(const Vector &x);
  double getG(void);
  int runGFunAnalysis(const Vector &x);
  int tokenizeSpecials(TCL_Char *theExpression, Tcl_Obj *paramList);
  
  void    setNsteps(int nsteps);
  double  getDt();
  int getNstep();
  //////////////////////////////////////////////////////////
  //// added by K Fujimura /////////////////////////////////
  //////////////////////////////////////////////////////////
  void activateSensitivty(void);
  void inactivateSensitivty(void);
  void setGFunEachStepEvaluator(GFunEachStepEvaluator *pGFunEachStepEvaluator);
  void inactivateGFunEachStepEvaluator();
  void setThreshold(double value){pfthreshold=value;}
  double getThreshold(){return pfthreshold;}
  void setPerformFuncCoeffs(TaggedObjectStorage* pobject)
  { thePerformFuncCoeffs=pobject;}
  void setPerformFuncCoeffIter(PerformanceFunctionCoefficientIter* pobject)
  { thePfCoeffIter=pobject;}
  Matrix* getEachStepResult();
  Matrix* getEachStepConvFlag();
  
  
 protected:
  
 private:
  int createRecorders();
  int removeRecorders();
  double PerformanceFunction();
  char *rec_node_occurrence(char tempchar[100], bool createRecorders, int &line, int &column);
  char *rec_element_occurrence(char tempchar[100], bool createRecorders, int &line, int &column);
  Analyzer *theAnalyzer;
  
  Domain* theDomain;
  TaggedObjectStorage* thePerformFuncCoeffs;
  PerformanceFunctionCoefficientIter* thePfCoeffIter;
  double pfthreshold;
  
  Tcl_Interp *theTclInterp;
  ReliabilityDomain *theReliabilityDomain;
  
  double g;
  int numberOfEvaluations;
};

#endif
