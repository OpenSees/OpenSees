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
// $Date: 2003-10-27 23:45:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/stepSize/ArmijoStepSizeRule.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ArmijoStepSizeRule.h>
#include <GFunEvaluator.h>
#include <StepSizeRule.h>
#include <ProbabilityTransformation.h>
#include <MeritFunctionCheck.h>
#include <RootFinding.h>
#include <math.h>
#include <Vector.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

ArmijoStepSizeRule::ArmijoStepSizeRule(	GFunEvaluator *passedGFunEvaluator,
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
	theGFunEvaluator = passedGFunEvaluator;
	theProbabilityTransformation = passedProbabilityTransformation;
	theMeritFunctionCheck = passedMeritFunctionCheck;
	theRootFindingAlgorithm = passedRootFindingAlgorithm;
	gFunValue = 0;
	base = Pbase; 
	maxNumReductions = PmaxNumReductions; 
	b0 = Pb0; 
	numberOfShortSteps = PnumberOfShortSteps; 
	radius = Pradius; 
	surfaceDistance = PsurfaceDistance; 
	evolution = Pevolution; 
	isCloseToSphere = false;
	printFlag = pprintFlag;
}

ArmijoStepSizeRule::~ArmijoStepSizeRule()
{
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



double 
ArmijoStepSizeRule::getGFunValue()
{
	return gFunValue;
}


int
ArmijoStepSizeRule::computeStepSize(Vector u_old, 
									Vector grad_G_old, 
									double g_old, 
									Vector dir_old,
									int stepNumber)
{

	// Initial declarations
	bool isOutsideSphere;
	bool isSecondTime;
	bool FEconvergence;
	double result;
	Vector x_new;
	Matrix jacobian_x_u;


	// Inform user through log file
	static ofstream logfile( "ArmijoRuleLog.txt", ios::out );
	logfile << "Entering the Armijo rule to compute step size ..." << endln;
	logfile.flush();
	
	
	// Check that the norm of the gradient is not zero
	double gradNorm = grad_G_old.Norm();
	if (gradNorm == 0.0) {
		opserr << "ArmijoStepSizeRule::computeStepSize() - the norm " << endln
			<< " of the gradient is zero. " << endln;
		return -1;
	}


	// Check if this potentially could be the second time close to the surface
	if (isCloseToSphere) {
		isSecondTime = true;
	}
	else {
		isSecondTime = false;
	}


	// Set the first trial step size
	double lambda_new;
	if (stepNumber <= numberOfShortSteps) {
		lambda_new = b0;
	}
	else {
		lambda_new = 1.0;
	}
	

	// Take a trial step in standard normal space
	Vector u_new = u_old + dir_old * lambda_new;


	// Inform the user
	if (printFlag != 0) {
		opserr << "Armijo starting gFun evaluation at distance " << u_new.Norm() << "..." << endln
			<< " .......: ";
	}
	logfile << "Armijo starting gFun evaluation at distance " << u_new.Norm() << "..." << endln;
	logfile.flush();


	double g_new;
	Vector grad_G_new;
	if (u_new.Norm()>radius) {

		isOutsideSphere = true;

		// Set some dummy values to 'g' and 'gradG'; it doesn't matter
		// that the merit functions then will be 'wrong'; the step will
		// fail in any case because it is outside the sphere. 
		g_new = g_old;
		grad_G_new = grad_G_old;


		// Inform the user
		if (printFlag != 0) {
			opserr << "Armijo skipping gFun evaluation because of hyper sphere requirement..." << endln
				<< " .......: ";
		}
		logfile << "Armijo skipping gFun evaluation because of hyper sphere requirement..." << endln;
		logfile.flush();


	}
	else {

		isOutsideSphere = false;

		// Register where this u was close to the sphere surface
		if (u_new.Norm() > radius-surfaceDistance  && 
			u_new.Norm() < radius+surfaceDistance) {
			isCloseToSphere = true;
			if (isSecondTime) {
				radius = radius + evolution;
			}
		}
		else {
			isCloseToSphere = false;
		}


		// Transform the trial point into original space
		result = theProbabilityTransformation->set_u(u_new);
		if (result < 0) {
			opserr << "ArmijoStepSizeRule::computeStepSize() - could not set " << endln
				<< " vector u in the xu-transformation. " << endln;
			return -1;
		}
		result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) {
			opserr << "ArmijoStepSizeRule::computeStepSize() - could not  " << endln
				<< " transform u to x. " << endln;
			return -1;
		}
		x_new = theProbabilityTransformation->get_x();
		jacobian_x_u = theProbabilityTransformation->getJacobian_x_u();


		// Evaluate the limit-state function
		FEconvergence = true;
		result = theGFunEvaluator->runGFunAnalysis(x_new);
/*		if (result < 0) {
			// In this case the FE analysis did not converge
			// In this case; do not accept the new point!
			FEconvergence = false;
			opserr << "step size rejected due to non-converged FE analysis ..." << endln
				<< " .......: ";
		}
*/		result = theGFunEvaluator->evaluateG(x_new);
		if (result < 0) {
			opserr << "ArmijoStepSizeRule::computeStepSize() - could not  " << endln
				<< " tokenize the limit-state function. " << endln;
			return -1;
		}
		g_new = theGFunEvaluator->getG();

	}


	///////////////////////////////////////////
	// POSSIBLY START REDUCING THE STEP SIZE //
	///////////////////////////////////////////
	int i = 1;
	bool mustGoOn = false;

	if (theMeritFunctionCheck->check(u_old, g_old, grad_G_old, lambda_new, dir_old, g_new)<0) {
		mustGoOn = true;
	}
	if (!FEconvergence) {
		mustGoOn = true;
	}
	if (isOutsideSphere) {
		mustGoOn = true;
	}
	if (i>maxNumReductions) {
		mustGoOn = false;
	}


	while ( mustGoOn ) {

		
		// Notify user that step sizes are being reduced
		opserr << "Armijo trial point rejected; reducing step size..." << endln
			<< " .......: ";
		logfile << "Armijo trial point rejected; reducing step size..." << endln;
		logfile.flush();
		

		// Cut the step size in half (or whichever base the user has set) and try that instead
		if (stepNumber <= numberOfShortSteps) {
			lambda_new = b0 * pow(base,i);
		}
		else {
			lambda_new = pow(base,i);
		}
		u_new = u_old + dir_old * lambda_new;


		// Check if we are beyond the bounding sphere
		if (u_new.Norm()>radius) {

			isOutsideSphere = true;

			// Inform the user
			if (printFlag != 0) {
				opserr << "Armijo skipping gFun evaluation because of hyper sphere requirement..." << endln
					<< " .......: ";
			}
			logfile << "Armijo skipping gFun evaluation because of hyper sphere requirement..." << endln;
			logfile.flush();
		}
		else {

			isOutsideSphere = false;

			// Register where this u was close to the sphere surface
			if (u_new.Norm() > radius-surfaceDistance  && 
				u_new.Norm() < radius+surfaceDistance) {
				isCloseToSphere = true;
				if (isSecondTime) {
					radius = radius + evolution;
				}
			}
			else {
				isCloseToSphere = false;
			}
		
			if (printFlag != 0) {
				opserr << "Armijo starting gFun evaluation at distance " << u_new.Norm() << "..." << endln
					<< " .......: ";
			}
			logfile << "Armijo starting gFun evaluation at distance " << u_new.Norm() << "..." << endln;

			// Transform the trial point into original space
			double result = theProbabilityTransformation->set_u(u_new);
			if (result < 0) {
				opserr << "ArmijoStepSizeRule::computeStepSize() - could not set " << endln
					<< " vector u in the xu-transformation. " << endln;
				return -1;
			}
			result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
			if (result < 0) {
				opserr << "ArmijoStepSizeRule::computeStepSize() - could not  " << endln
					<< " transform u to x. " << endln;
				return -1;
			}
			x_new = theProbabilityTransformation->get_x();
			jacobian_x_u = theProbabilityTransformation->getJacobian_x_u();


			// Evaluate the limit-state function
			FEconvergence = true;
			result = theGFunEvaluator->runGFunAnalysis(x_new);
/*			if (result < 0) {
				// In this case the FE analysis did not converge
				// May still accept the new point!
				FEconvergence = false;
				opserr << "step size rejected due to non-converged FE analysis ..." << endln
					<< " .......: ";
				logfile << "step size rejected due to non-converged FE analysis ..." << endln;
			}
*/			result = theGFunEvaluator->evaluateG(x_new);
			if (result < 0) {
				opserr << "ArmijoStepSizeRule::computeStepSize() - could not  " << endln
					<< " tokenize the limit-state function. " << endln;
				return -1;
			}
			g_new = theGFunEvaluator->getG();


			// Possibly project the point onto the limit-state surface
			// (it is assumed that the user has provided a projection object 
			// if this is to be done)
			// Note that this has effect on the previously computed search direction
			// because the new values of the performance function and its gradient
			// is being used...
			// Q: Is the u_new point also being kept in the orchestrating algorithm?
			if (theRootFindingAlgorithm != 0) {
				theRootFindingAlgorithm->findLimitStateSurface(2,g_new, grad_G_old, u_new);
			}
		}

		// Increment counter
		i++;

		// Check if we need to go on
		mustGoOn = false;

		if (theMeritFunctionCheck->check(u_old, g_old, grad_G_old, lambda_new, dir_old, g_new)<0) {
			mustGoOn = true;
		}
		if (!FEconvergence) {
			mustGoOn = true;
		}
		if (isOutsideSphere) {
			mustGoOn = true;
		}
		if (i>maxNumReductions) {
			mustGoOn = false;
		}

 	
	}

	stepSize = lambda_new;
	gFunValue = g_new;

	return 1;

}
