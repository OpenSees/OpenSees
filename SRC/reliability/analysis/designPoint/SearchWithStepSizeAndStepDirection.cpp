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
                                                                        
// $Revision: 1.12 $
// $Date: 2007-07-13 19:54:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/designPoint/SearchWithStepSizeAndStepDirection.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <SearchWithStepSizeAndStepDirection.h>
#include <FindDesignPointAlgorithm.h>
#include <StepSizeRule.h>
#include <SearchDirection.h>
#include <ProbabilityTransformation.h>
#include <NatafProbabilityTransformation.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <RandomVariable.h>
#include <CorrelationCoefficient.h>
#include <MatrixOperations.h>
#include <HessianApproximation.h>
#include <ReliabilityConvergenceCheck.h>
#include <Matrix.h>
#include <Vector.h>
#include <GammaRV.h>

#include <fstream>
#include <iomanip>
#include <iostream>

using std::ifstream;
using std::ios;

using std::setw;
using std::setprecision;



SearchWithStepSizeAndStepDirection::SearchWithStepSizeAndStepDirection(
					int passedMaxNumberOfIterations, 
					GFunEvaluator *passedGFunEvaluator,
					GradGEvaluator *passedGradGEvaluator,
					StepSizeRule *passedStepSizeRule,
					SearchDirection *passedSearchDirection,
					ProbabilityTransformation *passedProbabilityTransformation,
					HessianApproximation *passedHessianApproximation,
					ReliabilityConvergenceCheck *passedReliabilityConvergenceCheck,
					int pprintFlag,
					char *pFileNamePrint,
					Vector *pStartPoint)
:FindDesignPointAlgorithm()
{
	maxNumberOfIterations			= passedMaxNumberOfIterations;
	theGFunEvaluator				= passedGFunEvaluator;
	theGradGEvaluator				= passedGradGEvaluator;
	theStepSizeRule					= passedStepSizeRule;
	theSearchDirection				= passedSearchDirection;
	theProbabilityTransformation	= passedProbabilityTransformation;
	theHessianApproximation			= passedHessianApproximation;
	theReliabilityConvergenceCheck  = passedReliabilityConvergenceCheck;
	startPoint						= pStartPoint;
	printFlag						= pprintFlag;
	numberOfEvaluations =0;
	if (printFlag != 0) {
		strcpy(fileNamePrint,pFileNamePrint);
	}
	else {
		strcpy(fileNamePrint,"searchpoints.out");
	}
}



SearchWithStepSizeAndStepDirection::~SearchWithStepSizeAndStepDirection()
{
  
}


int
SearchWithStepSizeAndStepDirection::findDesignPoint(ReliabilityDomain *passedReliabilityDomain)
{

	// Set the reliability domain (a data member of this class)
	theReliabilityDomain = passedReliabilityDomain;

	// Declaration of data used in the algorithm
	int numberOfRandomVariables = theReliabilityDomain->getNumberOfRandomVariables();
	int j;
	int zeroFlag;
	Vector dummy(numberOfRandomVariables);
	x = dummy;
	u = dummy;
	Vector u_old(numberOfRandomVariables);
	uSecondLast = dummy;
	Vector uNew(numberOfRandomVariables);
	alpha = dummy;
	gamma = dummy;
	alphaSecondLast = dummy;
	double gFunctionValue = 1.0;
	double gFunctionValue_old = 1.0;
	Vector gradientOfgFunction(numberOfRandomVariables);
	gradientInStandardNormalSpace = dummy;
	Vector gradientInStandardNormalSpace_old(numberOfRandomVariables);
	double normOfGradient =0;
	double stepSize;
	//Matrix jacobian_x_u(numberOfRandomVariables,numberOfRandomVariables);
	int evaluationInStepSize = 0;
	int result;
	theGFunEvaluator->initializeNumberOfEvaluations();

	
	// Prepare output file to store the search points
	ofstream outputFile2( fileNamePrint, ios::out );
	if (printFlag == 0) {
		outputFile2 << "The user has not specified to store any search points." << endln;
		outputFile2 << "This is just a dummy file. " << endln;
	}


	if (startPoint == 0) {

		// Here we want to start at the origin in the standard normal space. 
		// Hence, just leave 'u' as being initialized to zero. 
		// (Note; this is slightly different from the mean point, as done earlier)
	}
	else {

		// Get starting point
		x = (*startPoint);


		// Transform starting point into standard normal space
		result = theProbabilityTransformation->set_x(x);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not set x in the xu-transformation." << endln;
			return -1;
		}


		result = theProbabilityTransformation->transform_x_to_u();
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not transform from x to u." << endln;
			return -1;
		}
		u = theProbabilityTransformation->get_u();
	}

	// Loop to find design point
	i = 1;
	while ( i <= maxNumberOfIterations )
	{

		// Transform from u to x space
		result = theProbabilityTransformation->set_u(u);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not set u in the xu-transformation." << endln;
			return -1;
		}

		result = theProbabilityTransformation->transform_u_to_x_andComputeJacobian();
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not transform from u to x and compute Jacobian." << endln;
			return -1;
		}
		x = theProbabilityTransformation->get_x();
		const Matrix &jacobian_x_u = theProbabilityTransformation->getJacobian_x_u();


		// Possibly print the point to output file
		int iii;
		if (printFlag != 0) {

			if (printFlag == 1) {
				outputFile2.setf(ios::scientific, ios::floatfield);
				for (iii=0; iii<x.Size(); iii++) {
					outputFile2<<setprecision(5)<<setw(15)<<x(iii)<<endln;
				}
			}
			else if (printFlag == 2) {
				outputFile2.setf(ios::scientific, ios::floatfield);
				for (iii=0; iii<u.Size(); iii++) {
					outputFile2<<setprecision(5)<<setw(15)<<u(iii)<<endln;
				}
			}
		}


		// Evaluate limit-state function unless it has been done in 
		// a trial step by the "stepSizeAlgorithm"
		if (evaluationInStepSize == 0) {
			result = theGFunEvaluator->runGFunAnalysis(x);
			if (result < 0) {
				opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
					<< " could not run analysis to evaluate limit-state function. " << endln;
				return -1;
			}
			result = theGFunEvaluator->evaluateG(x);
			if (result < 0) {
				opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
					<< " could not tokenize limit-state function. " << endln;
				return -1;
			}
			gFunctionValue_old = gFunctionValue;
			gFunctionValue = theGFunEvaluator->getG();
		}


		// Set scale parameter
		if (i == 1)	{
			Gfirst = gFunctionValue;
			opserr << " Limit-state function value at start point, g=" << gFunctionValue << endln;
			opserr << " STEP #0: ";
			theReliabilityConvergenceCheck->setScaleValue(gFunctionValue);
		}


		// Gradient in original space
		result = theGradGEvaluator->computeGradG(gFunctionValue,x);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute gradients of the limit-state function. " << endln;
			return -1;
		}
		gradientOfgFunction = theGradGEvaluator->getGradG();


		// Check if all components of the vector is zero
		zeroFlag = 0;
		for (j=0; j<gradientOfgFunction.Size(); j++) {
			if (gradientOfgFunction[j] != 0.0) {
				zeroFlag = 1;
			}
		}
		if (zeroFlag == 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " all components of the gradient vector is zero. " << endln;
			return -1;
		}


		// Gradient in standard normal space
		gradientInStandardNormalSpace_old = gradientInStandardNormalSpace;
		//gradientInStandardNormalSpace = jacobian_x_u ^ gradientOfgFunction;
		gradientInStandardNormalSpace.addMatrixTransposeVector(0.0, jacobian_x_u, gradientOfgFunction, 1.0);


		// Compute the norm of the gradient in standard normal space
		normOfGradient = gradientInStandardNormalSpace.Norm();


		// Check that the norm is not zero
		if (normOfGradient == 0.0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " the norm of the gradient is zero. " << endln;
			return -1;
		}

		
		// Compute alpha-vector
		alpha = gradientInStandardNormalSpace *  ( (-1.0) / normOfGradient );


		// Check convergence
		result = theReliabilityConvergenceCheck->check(u,gFunctionValue,gradientInStandardNormalSpace);
		if (result > 0)  {
		

			// Inform the user of the happy news!
			opserr << "Design point found!" << endln;


			// Print the design point to file, if desired
			int iii;
			if (printFlag != 0) {
				if (printFlag == 3) {

					result = theProbabilityTransformation->set_u(u);
					if (result < 0) {
						opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
							<< " could not set u in the xu-transformation." << endln;
						return -1;
					}

					result = theProbabilityTransformation->transform_u_to_x();
					if (result < 0) {
						opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
							<< " could not transform from u to x." << endln;
						return -1;
					}
					x = theProbabilityTransformation->get_x();

					outputFile2.setf(ios::scientific, ios::floatfield);
					for (iii=0; iii<x.Size(); iii++) {
						outputFile2<<setprecision(5)<<setw(15)<<x(iii)<<endln;
					}
				}
				else if (printFlag == 4) {
					static ofstream outputFile2( fileNamePrint, ios::out );
					outputFile2.setf(ios::scientific, ios::floatfield);
					for (iii=0; iii<u.Size(); iii++) {
						outputFile2<<setprecision(5)<<setw(15)<<u(iii)<<endln;
					}
				}
			}


			// Compute the gamma vector
			const Matrix &jacobian_u_x = theProbabilityTransformation->getJacobian_u_x();

			//Vector tempProduct = jacobian_u_x ^ alpha;
			Vector tempProduct(numberOfRandomVariables);
			tempProduct.addMatrixTransposeVector(0.0, jacobian_u_x, alpha, 1.0);

			// Only diagonal elements of (J_xu*J_xu^T) are used
			for (j = 0; j < numberOfRandomVariables; j++) {
			  double sum = 0.0;
			  double jk;
			  for (int k = 0; k < numberOfRandomVariables; k++) {
			    jk = jacobian_x_u(j,k);
			    sum += jk*jk;
			  }
			  gamma(j) = sqrt(sum) * tempProduct(j);
			}
			
			Glast = gFunctionValue;

			numberOfEvaluations = theGFunEvaluator->getNumberOfEvaluations();

			return 1;
		}


		// Store 'u' and 'alpha' at the second last iteration point
		uSecondLast = u;
		alphaSecondLast = alpha;


		// Let user know that we have to take a new step
		opserr << " STEP #" << i <<": ";


		// Update Hessian approximation, if any
		if (  (theHessianApproximation!=0) && (i!=1)  ) {
			theHessianApproximation->updateHessianApproximation(u_old,
											    gFunctionValue_old,
											    gradientInStandardNormalSpace_old,
											    stepSize,
											    searchDirection,
											    gFunctionValue,
											    gradientInStandardNormalSpace);
		}


		// Determine search direction
		result = theSearchDirection->computeSearchDirection(i,
			u, gFunctionValue, gradientInStandardNormalSpace );
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute search direction. " << endln;
			return -1;
		}
		searchDirection = theSearchDirection->getSearchDirection();


		// Determine step size
		result = theStepSizeRule->computeStepSize(
			u, gradientInStandardNormalSpace, gFunctionValue, searchDirection, i);
		if (result < 0) {  // (something went wrong)
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not compute step size. " << endln;
			return -1;
		}
		else if (result == 0) {  // (nothing was evaluated in step size)
			evaluationInStepSize = 0;
		}
		else if (result == 1) {  // (the gfun was evaluated)
			evaluationInStepSize = 1;
			gFunctionValue_old = gFunctionValue;
			gFunctionValue = theStepSizeRule->getGFunValue();
		}
		stepSize = theStepSizeRule->getStepSize();


		// Determine new iteration point (take the step)
		u_old = u;
		//u = u_old + (searchDirection * stepSize);
		u = u_old;
		u.addVector(1.0, searchDirection, stepSize);


		// Increment the loop parameter
		i++;

	}


	// Print a message if max number of iterations was reached
	// (Note: in this case the last trial point was never transformed/printed)
	opserr << "Maximum number of iterations was reached before convergence." << endln;

	
	return -1;
}



const Vector&
SearchWithStepSizeAndStepDirection::get_x()
{
	return x;
}

const Vector &
SearchWithStepSizeAndStepDirection::get_u()
{
	return u;
}

const Vector&
SearchWithStepSizeAndStepDirection::get_alpha()
{
	return alpha;
}

const Vector&
SearchWithStepSizeAndStepDirection::get_gamma()
{
  if (gamma.Norm() > 1.0)
    gamma = gamma / gamma.Norm();

  return gamma;
}

int
SearchWithStepSizeAndStepDirection::getNumberOfSteps()
{
	return (i-1);
}

const Vector&
SearchWithStepSizeAndStepDirection::getSecondLast_u()
{
	return uSecondLast;
}

const Vector&
SearchWithStepSizeAndStepDirection::getSecondLast_alpha()
{
	return alphaSecondLast;
}

const Vector&
SearchWithStepSizeAndStepDirection::getLastSearchDirection()
{
	return searchDirection;
}

double
SearchWithStepSizeAndStepDirection::getFirstGFunValue()
{
	return Gfirst;
}

double
SearchWithStepSizeAndStepDirection::getLastGFunValue()
{
	return Glast;
}


const Vector&
SearchWithStepSizeAndStepDirection::getGradientInStandardNormalSpace()
{
	return gradientInStandardNormalSpace;
}



int
SearchWithStepSizeAndStepDirection::getNumberOfEvaluations()
{
	return numberOfEvaluations;
}






