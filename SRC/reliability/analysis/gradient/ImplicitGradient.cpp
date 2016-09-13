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
                                                                        
// $Revision: 1.21 $
// $Date: 2010-09-13 21:38:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/ImplicitGradient.cpp,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <ImplicitGradient.h>
#include <Vector.h>
#include <GradientEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <SensitivityAlgorithm.h>
#include <string.h>
#include<Integrator.h>//Abbas


ImplicitGradient::ImplicitGradient(FunctionEvaluator *passedGFunEvaluator,
				   ReliabilityDomain *passedReliabilityDomain,
				   Domain *passedOpenSeesDomain,
				  // SensitivityAlgorithm *theAlgo)
				    Integrator *theAlgo)
  :GradientEvaluator(passedReliabilityDomain, passedGFunEvaluator), 
   theOpenSeesDomain(passedOpenSeesDomain)
{
	
  theSensAlgo = theAlgo;

  int nparam = theOpenSeesDomain->getNumParameters();
  grad_g = new Vector(nparam);
}


ImplicitGradient::~ImplicitGradient()
{
  if (grad_g != 0)
    delete grad_g;

}


const Vector&
ImplicitGradient::getGradient()
{
  return *grad_g;
}


int
ImplicitGradient::computeGradient(double g)
{
    // Compute gradients if this is a path-INdependent analysis
    // (This command only has effect if it IS path-independent.)
    //if (theSensAlgo != 0 && !(theSensAlgo->shouldComputeAtEachStep()) ) {
   if (theSensAlgo != 0){ 
        theSensAlgo->computeSensitivities();
   opserr<<" Implicit function SensAlgo is not zero"<<endln;
   }
    
    //}

    // Initialize gradient vector
    grad_g->Zero();

    // get limit-state function from reliability domain
    int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
    LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
    const char *lsfExpression = theLimitStateFunction->getExpression();
  
    // get parameters created in the domain
    int nparam = theOpenSeesDomain->getNumParameters();
    Vector partials(nparam);
  
  
    // first check for dg/dimplicit partials
    for (int i = 0; i < nparam; i++) {
        // get parameter tag
        Parameter *theParam = theOpenSeesDomain->getParameterFromIndex(i);
        int tag = theParam->getTag();
        if (theParam->isImplicit()) {

	  // check for analytic gradient first on dg/dimplicit
	  const char *gradExpression = theLimitStateFunction->getGradientExpression(tag);
	  if (gradExpression != 0) {
	    theFunctionEvaluator->setExpression(gradExpression);
            
	    if (theFunctionEvaluator->setVariables() < 0) {
	      opserr << "ERROR ImplicitGradient -- error setting variables in namespace" << endln;
	      return -1;
	    }
            
	    partials(i) = theFunctionEvaluator->evaluateExpression();
            
	    // Reset limit state function in evaluator -- subsequent calls could receive gradient expression
	    theFunctionEvaluator->setExpression(lsfExpression);
	  }
	  
	  // if no analytic gradient automatically do finite differences to get dg/dimplicit
	  else {
	    // use parameter defined perturbation after updating implicit parameter
	    theParam->update(0.0);
	    double h = theParam->getPerturbation();
	    double original = theParam->getValue();
	    theParam->setValue(original+h);
            
	    // set perturbed values in the variable namespace
	    if (theFunctionEvaluator->setVariables() < 0) {
	      opserr << "ERROR ImplicitGradient -- error setting variables in namespace" << endln;
	      return -1;
	    }
	    
	    // run analysis
	    //if (theFunctionEvaluator->runAnalysis() < 0) {
	    //    opserr << "ERROR ImplicitGradient -- error running analysis" << endln;
	    //    return -1;
	    //}
	    
	    // evaluate LSF and obtain result
	    theFunctionEvaluator->setExpression(lsfExpression);
            
	    // Add gradient contribution
	    double g_perturbed = theFunctionEvaluator->evaluateExpression();
	    partials(i) = (g_perturbed-g)/h;
            
	    // return values to previous state
	    theParam->update(0.0);
            
	    //opserr << "g_pert " << g_perturbed << ", g0 = " << g << endln;
	  }
        }
    }	
    
    //opserr << partials;

    // now loop through to create gradient vector
    // Mackie 7/31/2011: big consideration here is that you CANNOT have an explicit parameter appear in the 
    // same LSF as an implicit parameter.  For example, if there are two parameters: theta1 is modulus E and 
    // theta2 is nodal displacement, then g(theta) = theta1 + theta2
    // is not a viable LSF.  Note that the implicit computation would need to consider 
    // dg/dtheta1 + dg/du * du/dtheta1
    // If you insist on solving this, use FiniteDifferenceGradient
    for (int i = 0; i < nparam; i++) {
        // get parameter tag
        Parameter *theParam = theOpenSeesDomain->getParameterFromIndex(i);
        int tag = theParam->getTag();
        double result = 0;

        for (int j = 0; j < nparam; j++) {
            Parameter *theImplicit = theOpenSeesDomain->getParameterFromIndex(j);
            if (theImplicit->isImplicit()) {
                result = partials(j) * theImplicit->getSensitivity(i);
                (*grad_g)(i) += result;
            }
        }

    }
    
    return 0;
    
}
