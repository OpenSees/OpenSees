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
                                                                        
// $Revision: 1.2 $
// $Date: 2007-07-11 23:51:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/direction/GradientProjectionSearchDirection.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef GradientProjectionSearchDirection_h
#define GradientProjectionSearchDirection_h

#include <SearchDirection.h>
#include <StepSizeRule.h>
#include <ProbabilityTransformation.h>
#include <FunctionEvaluator.h>
#include <RootFinding.h>
#include <Vector.h>

class GradientProjectionSearchDirection : public SearchDirection
{
 public:
  GradientProjectionSearchDirection(StepSizeRule *theStepSizeRule,
				    ProbabilityTransformation *theProbabilityTransformation,
				    FunctionEvaluator *theGFunEvaluator,
				    RootFinding *theRootFindingAlgorithm);
  ~GradientProjectionSearchDirection();
  
  int computeSearchDirection(int stepNumber, const Vector &u,
			     double g, const Vector &gradG);
  const Vector &getSearchDirection();
  
 protected:

 private:
  Vector searchDirection;
  StepSizeRule *theStepSizeRule;
  ProbabilityTransformation *theProbabilityTransformation;
  FunctionEvaluator *theFunctionEvaluator;
  RootFinding *theRootFindingAlgorithm;
};

#endif
