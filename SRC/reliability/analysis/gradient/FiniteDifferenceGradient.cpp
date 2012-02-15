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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/FiniteDifferenceGradient.cpp,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <FiniteDifferenceGradient.h>
#include <Vector.h>
#include <GradientEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <string.h>


FiniteDifferenceGradient::FiniteDifferenceGradient(FunctionEvaluator *passedGFunEvaluator,
						   ReliabilityDomain *passedReliabilityDomain,
						   Domain *passedOpenSeesDomain)

:GradientEvaluator(passedReliabilityDomain, passedGFunEvaluator), 
theOpenSeesDomain(passedOpenSeesDomain)
{
	
	int nparam = theOpenSeesDomain->getNumParameters();
	grad_g = new Vector(nparam);
	
}


FiniteDifferenceGradient::~FiniteDifferenceGradient()
{
	if (grad_g != 0) 
		delete grad_g;
	
}


Vector
FiniteDifferenceGradient::getGradient()
{
	return *grad_g;
}


int
FiniteDifferenceGradient::computeGradient(double g)
{
	
	// Initialize gradient vector
	grad_g->Zero();
	
	// get limit-state function from reliability domain
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
	const char *lsfExpression = theLimitStateFunction->getExpression();
	
	// get parameters created in the domain
	int nparam = theOpenSeesDomain->getNumParameters();

	// initialize vector containing parameter values
	Vector param_vect(nparam);
	Vector param_pert(nparam);
	
	for (int i = 0; i < nparam; i++) {
		Parameter *theParam = theOpenSeesDomain->getParameterFromIndex(i);
		param_vect(i) = theParam->getValue();
	}
	
	// set perturbed values in the variable namespace
	if (theFunctionEvaluator->setVariables(param_vect) < 0) {
	  opserr << "ERROR FiniteDifferenceGradient -- error setting variables in namespace" << endln;
	  return -1;
	}
	
	// run analysis
	if (theFunctionEvaluator->runAnalysis(param_vect) < 0) {
	  opserr << "ERROR FiniteDifferenceGradient -- error running analysis" << endln;
	  return -1;
	}
	
	// evaluate LSF and obtain result
	theFunctionEvaluator->setExpression(lsfExpression);
	
	// Add gradient contribution
	g = theFunctionEvaluator->evaluateExpression();

	//opserr << "param vect " << param_vect;

	// now loop through to create gradient vector
	// note this is a for loop because there may be some conflict from a nested iterator already 
	// called at a higher level.
	for (int i = 0; i < nparam; i++) {

		// get parameter tag
		Parameter *theParam = theOpenSeesDomain->getParameterFromIndex(i);
		//opserr << theParam->getType() << endln;
		//if (strcmp(theParam->getType(),"FEModel") != 0)
		//continue;
		if (theParam->isImplicit())
		  continue;

		int tag = theParam->getTag();
		double result = 0;

		// check for analytic gradient first
		const char *gradExpression = theLimitStateFunction->getGradientExpression(tag);
		if (gradExpression != 0) {
		  theFunctionEvaluator->setExpression(gradExpression);

			if (theFunctionEvaluator->setVariables(param_vect) < 0) {
				opserr << "ERROR FiniteDifferenceGradient -- error setting variables in namespace" << endln;
				return -1;
			}
			
			result = theFunctionEvaluator->evaluateExpression();

			// Reset limit state function in evaluator -- subsequent calls could receive gradient expression
			theFunctionEvaluator->setExpression(lsfExpression);
		}
		
		// if no analytic gradient automatically do finite differences
		else {
			// use parameter defined perturbation
			double h = theParam->getPerturbation();
			//opserr << "h = " << h << endln;
			param_pert = param_vect;
			param_pert(i) += h;
			
			//opserr << "param_pert " << param_pert;

			theParam->setValue(param_pert(i));

			// set perturbed values in the variable namespace
			if (theFunctionEvaluator->setVariables(param_pert) < 0) {
				opserr << "ERROR FiniteDifferenceGradient -- error setting variables in namespace" << endln;
				return -1;
			}
			
			// run analysis
			if (theFunctionEvaluator->runAnalysis(param_pert) < 0) {
				opserr << "ERROR FiniteDifferenceGradient -- error running analysis" << endln;
				return -1;
			}
			
			theFunctionEvaluator->setVariables(param_pert);

			// evaluate LSF and obtain result
			theFunctionEvaluator->setExpression(lsfExpression);
			
			// Add gradient contribution
			double g_perturbed = theFunctionEvaluator->evaluateExpression();
			result = (g_perturbed-g)/h;
			theParam->setValue(param_vect(i));

			//opserr << "g_pert " << g_perturbed << ", g0 = " << g << endln;

			// reset original values in the variable namespace
			if (theFunctionEvaluator->setVariables(param_vect) < 0) {
			  opserr << "ERROR FiniteDifferenceGradient -- error setting variables in namespace" << endln;
			  return -1;
			}
		}
		
		(*grad_g)(i) = result;
		
	}
	

	return 0;
	
}

