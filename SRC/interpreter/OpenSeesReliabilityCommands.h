/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of
California (Regents). All rights reserved.

Redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following
conditions are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided
with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and
documentation are those of the authors and should not be
interpreted as representing official policies, either expressed
or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".
REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

***************************************************************************
*/

// Written: Minjie

// Description: all opensees APIs are defined or declared here
//

#ifndef OpenSeesReliabilityCommands_h
#define OpenSeesReliabilityCommands_h

#include <Domain.h>
#include <FORMAnalysis.h>
#include <FOSMAnalysis.h>
#include <FindDesignPointAlgorithm.h>
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <ImportanceSamplingAnalysis.h>
#include <Integrator.h>
#include <MeritFunctionCheck.h>
#include <PolakHeSearchDirectionAndMeritFunction.h>
#include <ProbabilityTransformation.h>
#include <RandomNumberGenerator.h>
#include <ReliabilityConvergenceCheck.h>
#include <ReliabilityDomain.h>
#include <RootFinding.h>
#include <SQPsearchDirectionMeritFunctionAndHessian.h>
#include <SearchDirection.h>
#include <StepSizeRule.h>

class OpenSeesReliabilityCommands {
 public:
  explicit OpenSeesReliabilityCommands(Domain *structuralDomain);
  ~OpenSeesReliabilityCommands();

  ReliabilityDomain *getDomain();
  Domain *getStructuralDomain();

  void setProbabilityTransformation(
      ProbabilityTransformation *transform);
  ProbabilityTransformation *getProbabilityTransformation() {
    return theProbabilityTransformation;
  }
  void setRandomNumberGenerator(
      RandomNumberGenerator *generator);
  RandomNumberGenerator *getRandomNumberGenerator() {
    return theRandomNumberGenerator;
  }

  void setReliabilityConvergenceCheck(
      ReliabilityConvergenceCheck *check);
  ReliabilityConvergenceCheck *getReliabilityConvergenceCheck() {
    return theReliabilityConvergenceCheck;
  }

  void setSearchDirection(SearchDirection *search);
  SearchDirection *getSearchDirection() {
    return theSearchDirection;
  }

  void setMeritFunctionCheck(MeritFunctionCheck *merit);
  MeritFunctionCheck *getMeritFunctionCheck() {
    return theMeritFunctionCheck;
  }

  void setStepSizeRule(StepSizeRule *rule);
  StepSizeRule *getStepSizeRule() { return theStepSizeRule; }

  void setRootFinding(RootFinding *root);
  RootFinding *getRootFinding() { return theRootFinding; }

  void setFindDesignPointAlgorithm(
      FindDesignPointAlgorithm *algo);
  FindDesignPointAlgorithm *getFindDesignPointAlgorithm() {
    return theFindDesignPointAlgorithm;
  }

  void setFunctionEvaluator(FunctionEvaluator *eval);
  FunctionEvaluator *getFunctionEvaluator() {
    return theFunctionEvaluator;
  }

  void setGradientEvaluator(GradientEvaluator *eval);
  GradientEvaluator *getGradientEvaluator() {
    return theGradientEvaluator;
  }

  void setPolakHeDualPurpose(
      PolakHeSearchDirectionAndMeritFunction *fn) {
    thePolakHeDualPurpose = fn;
  }
  PolakHeSearchDirectionAndMeritFunction *
  getPolakHeDualPurpose() {
    return thePolakHeDualPurpose;
  }

  void setSQPtriplePurpose(
      SQPsearchDirectionMeritFunctionAndHessian *fn) {
    theSQPtriplePurpose = fn;
  }
  SQPsearchDirectionMeritFunctionAndHessian *
  getSQPtriplePurpose() {
    return theSQPtriplePurpose;
  }

  void setFOSMAnalysis(FOSMAnalysis *analysis);
  FOSMAnalysis *getFOSMAnalysis() { return theFOSMAnalysis; }

  void setFORMAnalysis(FORMAnalysis *analysis);
  FORMAnalysis *getFORMAnalysis() { return theFORMAnalysis; }

  void setImportanceSamplingAnalysis(
      ImportanceSamplingAnalysis *analysis);
  ImportanceSamplingAnalysis *getImportanceSamplingAnalysis() {
    return theImportanceSamplingAnalysis;
  }

  void setSensitivityAlgorithm(Integrator *inte) {
    theSensAlgo = inte;
  }
  Integrator *getSensitivityAlgorithm() { return theSensAlgo; }

  void wipe();

 private:
  ReliabilityDomain *theDomain;
  Domain *theStructuralDomain;

  ProbabilityTransformation *theProbabilityTransformation;
  RandomNumberGenerator *theRandomNumberGenerator;
  ReliabilityConvergenceCheck *theReliabilityConvergenceCheck;
  SearchDirection *theSearchDirection;
  MeritFunctionCheck *theMeritFunctionCheck;
  StepSizeRule *theStepSizeRule;
  RootFinding *theRootFinding;
  FindDesignPointAlgorithm *theFindDesignPointAlgorithm;

  FunctionEvaluator *theFunctionEvaluator;
  GradientEvaluator *theGradientEvaluator;

  // not independent pointers
  PolakHeSearchDirectionAndMeritFunction *thePolakHeDualPurpose;
  SQPsearchDirectionMeritFunctionAndHessian *theSQPtriplePurpose;

  // analysis
  FOSMAnalysis *theFOSMAnalysis;
  FORMAnalysis *theFORMAnalysis;
  ImportanceSamplingAnalysis *theImportanceSamplingAnalysis;

  // sensitivity algorithm
  Integrator *theSensAlgo;
};

ReliabilityDomain *OPS_GetReliabilityDomain();

#endif
