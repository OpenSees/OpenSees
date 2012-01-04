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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/direction/GradientProjectionSearchDirection.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <GradientProjectionSearchDirection.h>
#include <SearchDirection.h>
#include <StepSizeRule.h>
#include <ProbabilityTransformation.h>
#include <FunctionEvaluator.h>
#include <RootFinding.h>
#include <Vector.h>
#include <Matrix.h>


GradientProjectionSearchDirection::GradientProjectionSearchDirection(StepSizeRule *passedStepSizeRule,
								     ProbabilityTransformation *passedProbabilityTransformation,
								     FunctionEvaluator *passedGFunEvaluator,
								     RootFinding *passedRootFindingAlgorithm)
  :SearchDirection()
{
  theStepSizeRule = passedStepSizeRule;
  theProbabilityTransformation = passedProbabilityTransformation;
  theFunctionEvaluator = passedGFunEvaluator;
  theRootFindingAlgorithm = passedRootFindingAlgorithm;
}

GradientProjectionSearchDirection::~GradientProjectionSearchDirection()
{

}

const Vector&
GradientProjectionSearchDirection::getSearchDirection()
{
  return searchDirection;
}

int
GradientProjectionSearchDirection::computeSearchDirection(int stepNumber, 
							  const Vector &u, 
							  double passed_g, 
							  const Vector &gradG)
{
  // Initial declarations
  int i,j;
  Vector u_new;
  double initialStepSize;
  Vector Direction;
  
  
  // Problem size 
  int nrv = u.Size();
  
  
  // Unit matrix
  Matrix I(nrv,nrv);
  for (i=0; i<nrv; i++) {
    for (j=0; j<nrv; j++) {
      if (i==j) {
	I(i,j) = 1.0;
      }
      else {
	I(i,j) = 0.0;
      }
    }
  }
  
  
  // Matrix of "outer" product of gradient vector
  Matrix dGdG(nrv,nrv);
  for (i=0; i<nrv; i++) {
    for (j=0; j<nrv; j++) {
      dGdG(i,j) = gradG(i)*gradG(j);
    }
  }
  
  // Get initial step size from the step size algorithm
  initialStepSize = theStepSizeRule->getInitialStepSize();
  
  
  // As long as it is not the first step; do the usual thing
  // (shouldn't happen if the user restarts the search...)
  if (stepNumber != 1) {
    
    // Compute the initial search direction vector
    Vector direction = (-1)* (I - (1.0/(gradG^gradG))*dGdG ) * u;
    
    // Initial step 
    u_new = u + initialStepSize*direction;
    
    // Set the direction of the Newton search
    Direction = gradG;
  }
  // If it's the first step; do the Newton thing from the
  // start point in the direction of the iHLRF search direction. 
  else {
    
    u_new = u;
    
    // Compute the alpha-vector
    Vector alpha = gradG * ( (-1) / gradG.Norm() );
    
    // Compute the direction vector
    double alpha_times_u = alpha ^ u ;
    Vector direction = alpha * ( passed_g / gradG.Norm() + alpha_times_u ) - u;
    
    // Set the direction of the Newton search
    Direction = (-1)*direction;
  }
  
  
  // Do the search to bring the trial point 'u_new' onto the lsf surface
  double tangent = gradG.Norm();
  Vector u_newest = theRootFindingAlgorithm->findLimitStateSurface(2,passed_g, Direction, u_new);
  
  
  // Return the final search direction
  // (remember to scale it so that u_new = u_old + lambda*d really brings us to u_new
  searchDirection = (1.0/initialStepSize)*(u_newest-u);
  
  return 0;
}
