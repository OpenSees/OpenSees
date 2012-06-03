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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/sensitivity/FiniteDifferenceHessian.cpp,v $

//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <FiniteDifferenceHessian.h>
#include <Matrix.h>
#include <HessianEvaluator.h>
#include <ReliabilityDomain.h>
#include <LimitStateFunction.h>
#include <string.h>


FiniteDifferenceHessian::FiniteDifferenceHessian(FunctionEvaluator *passedGFunEvaluator,
						   ReliabilityDomain *passedReliabilityDomain,
						   Domain *passedOpenSeesDomain)

:HessianEvaluator(passedReliabilityDomain, passedGFunEvaluator), 
theOpenSeesDomain(passedOpenSeesDomain)
{
	
	int nparam = theOpenSeesDomain->getNumParameters();
	grad_g = new Matrix(nparam,nparam);
	
}


FiniteDifferenceHessian::~FiniteDifferenceHessian()
{
	if (grad_g != 0) 
		delete grad_g;
	
}


Matrix
FiniteDifferenceHessian::getHessian()
{
	return *grad_g;
}


int
FiniteDifferenceHessian::computeHessian()
{
    // Initialize Hessian vector
	grad_g->Zero();
	
	// get limit-state function from reliability domain
    // NOTE: this needs to change in the future to allow treatment of other 
    // types of performanceFunctions (not just LSF)
	int lsf = theReliabilityDomain->getTagOfActiveLimitStateFunction();
	LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsf);
	const char *lsfExpression = theLimitStateFunction->getExpression();
	
	// get parameters created in the domain
	int nparam = theOpenSeesDomain->getNumParameters();
    Vector stored_g(nparam);
    
	// note FiniteDifferentHessian does not have g passed in so needs to establish base state
    if (theFunctionEvaluator->setVariables() < 0) {
        opserr << "ERROR FiniteDifferenceHessian -- error setting variables in namespace" << endln;
        return -1;
    }
    
    // run analysis
    if (theFunctionEvaluator->runAnalysis() < 0) {
        opserr << "ERROR FiniteDifferenceHessian -- error running analysis" << endln;
        return -1;
    }
    
    // evaluate LSF and obtain result
    theFunctionEvaluator->setExpression(lsfExpression);
    double g = theFunctionEvaluator->evaluateExpression();
    
	// now loop through to create hessian matrix
	// note this is a for loop because there may be some conflict from a nested iterator already 
	// called at a higher level.
	for (int i = 0; i < nparam; i++) {
		// get parameter tag
		Parameter *theParam_i = theOpenSeesDomain->getParameterFromIndex(i);
		double result = 0;

        // use parameter defined perturbation for forward difference --------------------------
        double hi = theParam_i->getPerturbation();
        double originali = theParam_i->getValue();
        theParam_i->update(originali+hi);

        // set perturbed values in the variable namespace
        if (theFunctionEvaluator->setVariables() < 0) {
            opserr << "ERROR FiniteDifferenceHessian -- error setting variables in namespace" << endln;
            return -1;
        }
        
        // run analysis
        if (theFunctionEvaluator->runAnalysis() < 0) {
            opserr << "ERROR FiniteDifferenceHessian -- error running analysis" << endln;
            return -1;
        }
        
        // evaluate LSF and obtain result
        theFunctionEvaluator->setExpression(lsfExpression);
        double g_perturbed_forward = theFunctionEvaluator->evaluateExpression();
        stored_g(i) = g_perturbed_forward;
        
        // use parameter defined perturbation for backward difference --------------------------
        theParam_i->update(originali-hi);
        
        // set perturbed values in the variable namespace
        if (theFunctionEvaluator->setVariables() < 0) {
            opserr << "ERROR FiniteDifferenceHessian -- error setting variables in namespace" << endln;
            return -1;
        }
        
        // run analysis
        if (theFunctionEvaluator->runAnalysis() < 0) {
            opserr << "ERROR FiniteDifferenceHessian -- error running analysis" << endln;
            return -1;
        }
        
        // evaluate LSF and obtain result
        theFunctionEvaluator->setExpression(lsfExpression);
        double g_perturbed_backward = theFunctionEvaluator->evaluateExpression();
        
        // diagonal Hessian from central difference approximation with 2*h perturbation
        result = (g_perturbed_forward - 2.0*g + g_perturbed_backward)/hi/hi;
        (*grad_g)(i,i) = result;
        
        // return values to previous state
        theParam_i->update(originali);
        
        // now get into off-diagonal perturbations
        for (int j = 0; j < i; j++) {
            // get parameter tag
            Parameter *theParam_j = theOpenSeesDomain->getParameterFromIndex(j);
            
            // use parameter defined perturbation for forward difference --------------------------
            double hj = theParam_j->getPerturbation();
            double originalj = theParam_j->getValue();
            theParam_j->update(originalj+hj);
            theParam_i->update(originali+hi);
            
            // set perturbed values in the variable namespace
            if (theFunctionEvaluator->setVariables() < 0) {
                opserr << "ERROR FiniteDifferenceHessian -- error setting variables in namespace" << endln;
                return -1;
            }
            
            // run analysis
            if (theFunctionEvaluator->runAnalysis() < 0) {
                opserr << "ERROR FiniteDifferenceHessian -- error running analysis" << endln;
                return -1;
            }
            
            // evaluate LSF and obtain result
            theFunctionEvaluator->setExpression(lsfExpression);
            double g_perturbed_forward_off = theFunctionEvaluator->evaluateExpression();
            
            // off-diagonal Hessian from forward difference approximation
            (*grad_g)(i,j) = ( g_perturbed_forward_off - stored_g(j) - stored_g(i) + g )/hi/hj;
            (*grad_g)(j,i) = (*grad_g)(i,j);
            
            // return values to previous state
            theParam_j->update(originalj);
            theParam_i->update(originali);
        }
    }

	return 0;
	
}

