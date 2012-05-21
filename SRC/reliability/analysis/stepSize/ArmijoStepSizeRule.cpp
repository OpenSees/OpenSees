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
                                                                        
// $Revision: 1.8 $
// $Date: 2010-09-13 21:35:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/stepSize/ArmijoStepSizeRule.cpp,v $


//
// Written by: 
// Kevin Mackie (kmackie@mail.ucf.edu)
// Michael Scott (mhscott@engr.orst.edu)
//

#include <ArmijoStepSizeRule.h>
#include <FunctionEvaluator.h>
#include <StepSizeRule.h>
#include <MeritFunctionCheck.h>
#include <RootFinding.h>
#include <math.h>
#include <Vector.h>

ArmijoStepSizeRule::ArmijoStepSizeRule(ReliabilityDomain *passedRelDomain, 
                                       FunctionEvaluator *passedGFunEvaluator,
				       ProbabilityTransformation *passedProbabilityTransformation,
				       MeritFunctionCheck *passedMeritFunctionCheck,
				       RootFinding *passedRootFindingAlgorithm, 
				       double Pbase,
				       int    PmaxNumReductions,
				       double Pb0,
				       int    PnumberOfShortSteps,
				       double Pradius,
				       double PsurfaceDistance,
				       double Pevolution,
				       int pprintFlag)
:StepSizeRule()
{
	theReliabilityDomain = passedRelDomain;
    theGFunEvaluator = passedGFunEvaluator;
	theProbabilityTransformation = passedProbabilityTransformation;
	theMeritFunctionCheck = passedMeritFunctionCheck;
	theRootFindingAlgorithm = passedRootFindingAlgorithm;

	base = Pbase; 
	maxNumReductions = PmaxNumReductions; 
	b0 = Pb0; 
	numberOfShortSteps = PnumberOfShortSteps; 
	radius = Pradius; 
	surfaceDistance = PsurfaceDistance; 
	evolution = Pevolution; 
	printFlag = pprintFlag;
    
    isCloseToSphere = false;
    isOutsideSphere = false;
    isSecondTime = false;
    FEconvergence = false;

}

ArmijoStepSizeRule::~ArmijoStepSizeRule()
{
    
}


int
ArmijoStepSizeRule::initialize()
{
    isCloseToSphere = false;
	numReduction = 0;
    FEconvergence = false;

    return 0;
}


double 
ArmijoStepSizeRule::getStepSize()
{
	return stepSize;
}



double 
ArmijoStepSizeRule::getInitialStepSize()
{
	return b0;
}


int
ArmijoStepSizeRule::computeStepSize(const Vector &u_old, 
				    const Vector &grad_G_old, 
				    double g_old, 
				    const Vector &dir_old,
				    int stepNumber,
					int reschk)
{

    // Check if this potentially could be the second time close to the surface
	if (isCloseToSphere)
		isSecondTime = true;
	else
		isSecondTime = false;
    
    
    // Set the first trial step size
    if (stepNumber <= numberOfShortSteps) {
        stepSize = b0;
        return 0;
    }
    
    // begin step size reduction
    numReduction++;
    stepSize = pow(base,numReduction);

    // Inform the user
    if (printFlag != 0) {
        opserr << "Armijo starting gFun evaluation at distance " << u_old.Norm() << "..." << endln
            << " .......: ";
    }

    if (u_old.Norm() > radius) {
        isOutsideSphere = true;

        // Set some dummy values to 'g' and 'gradG'; it doesn't matter
        // that the merit functions then will be 'wrong'; the step will
        // fail in any case because it is outside the sphere. 
        //g_new = g_old;
        //grad_G_new = grad_G_old;

        // Inform the user
        if (printFlag != 0) {
            opserr << "Armijo skipping gFun evaluation because of hyper sphere requirement..." << endln
                << " .......: ";
        }

    }
    else {
        isOutsideSphere = false;

        // Register where this u was close to the sphere surface
        if (u_old.Norm() > radius-surfaceDistance  && 
            u_old.Norm() < radius+surfaceDistance) {
            isCloseToSphere = true;
            if (isSecondTime) {
                radius = radius + evolution;
            }
        }
        else {
            isCloseToSphere = false;
        }
    }

    // KRM 4/9/2012
    // things that are missing here include previous call to theRootFindingAlgorithm and 
    // use of the meritFunctionCheck. Also need (potentially) to retain the ability to not
    // evaluate the next g function if deemed out of the hypersphere here in the stepSize.
    // Already possible, but designPointAlgorith needs some logic based on return value here.
    
    // check exit status
	if (numReduction > maxNumReductions) {
        return 0;
	}

    return 1;

}


/////S added by K Fujimura /////
int
ArmijoStepSizeRule::getNumReductions()
{
	return numReduction;

}
/////E added by K Fujimura /////
