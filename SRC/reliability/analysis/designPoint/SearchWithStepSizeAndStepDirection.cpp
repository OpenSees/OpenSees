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
                                                                        
// $Revision: 1.20 $
// $Date: 2010-09-13 21:34:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/designPoint/SearchWithStepSizeAndStepDirection.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <SearchWithStepSizeAndStepDirection.h>
#include <FindDesignPointAlgorithm.h>
#include <ReliabilityDomain.h>
#include <StepSizeRule.h>
#include <SearchDirection.h>
#include <ProbabilityTransformation.h>
#include <LimitStateFunction.h>
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>

#include <RandomVariable.h>
#include <Parameter.h>
#include <MatrixOperations.h>
#include <HessianEvaluator.h>
#include <ReliabilityConvergenceCheck.h>
#include <Matrix.h>
#include <Vector.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>

using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;


SearchWithStepSizeAndStepDirection::SearchWithStepSizeAndStepDirection(
					int passedMaxNumberOfIterations,
					ReliabilityDomain *passedReliabilityDomain,
                    Domain *passedOpenSeesDomain,
					FunctionEvaluator *passedFunctionEvaluator,
					GradientEvaluator *passedGradientEvaluator,
					StepSizeRule *passedStepSizeRule,
					SearchDirection *passedSearchDirection,
					ProbabilityTransformation *passedProbabilityTransformation,
					ReliabilityConvergenceCheck *passedReliabilityConvergenceCheck,
					int pprintFlag,
					char *pFileNamePrint)
  :FindDesignPointAlgorithm(), theReliabilityDomain(passedReliabilityDomain), 
        theOpenSeesDomain(passedOpenSeesDomain)
{
	maxNumberOfIterations			= passedMaxNumberOfIterations;
	theFunctionEvaluator			= passedFunctionEvaluator;
	theGradientEvaluator			= passedGradientEvaluator;
	theStepSizeRule					= passedStepSizeRule;
	theSearchDirection				= passedSearchDirection;
	theProbabilityTransformation	= passedProbabilityTransformation;
	theReliabilityConvergenceCheck  = passedReliabilityConvergenceCheck;
	printFlag						= pprintFlag;
	numberOfEvaluations = 0;
	if (printFlag != 0) {
		strcpy(fileNamePrint,pFileNamePrint);
	}
	
    // size objects
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();
	x = new Vector(nrv);
	u = new Vector(nrv);
	alpha = new Vector(nrv);
	gamma = new Vector(nrv);
	gradientInStandardNormalSpace = new Vector(nrv);
    gradientInOriginalSpace = new Vector(nrv);
	uSecondLast = new Vector(nrv);
	alphaSecondLast = new Vector(nrv);
	searchDirection = new Vector(nrv);
    Jux = new Matrix(nrv,nrv);
    Jxu = new Matrix(nrv,nrv);
}



SearchWithStepSizeAndStepDirection::~SearchWithStepSizeAndStepDirection()
{
	if (x != 0)
		delete x;
	if (u != 0)
		delete u;
	if (alpha != 0)
		delete alpha;
	if (gamma != 0)
		delete gamma;
	if (gradientInStandardNormalSpace != 0)
		delete gradientInStandardNormalSpace;
    if (gradientInOriginalSpace != 0)
        delete gradientInOriginalSpace;
	if (uSecondLast != 0)
		delete uSecondLast;
	if (alphaSecondLast != 0)
		delete alphaSecondLast;
    if (searchDirection != 0)
		delete searchDirection;
    if (Jux != 0)
        delete Jux;
    if (Jxu != 0)
        delete Jxu;
}


int
SearchWithStepSizeAndStepDirection::gradientStandardNormal(double gFunctionValue)
{
    // constants
    int numberOfRandomVariables = theReliabilityDomain->getNumberOfRandomVariables();
    int numberOfParameters = theOpenSeesDomain->getNumParameters();
    int result = 0;
    double normOfGradient = 0;
    
    // temporary storage
    Vector temp_grad(numberOfParameters);
    
    
    // Gradient in original space (make sure gFunctionValue is current)
    result = theGradientEvaluator->computeGradient(gFunctionValue);
    if (result < 0) {
        opserr << "SearchWithStepSizeAndStepDirection::gradientStandardNormal() - " << endln
               << " could not compute gradients of the limit-state function. " << endln;
        return -1;
    }
    temp_grad = theGradientEvaluator->getGradient();
    
    // map gradient from all parameters to just RVs
    for (int j = 0; j < numberOfRandomVariables; j++) {
        int param_indx = theReliabilityDomain->getParameterIndexFromRandomVariableIndex(j);
        (*gradientInOriginalSpace)(j) = temp_grad(param_indx);
    }
    
    // Get Jacobian x-space to u-space
    result = theProbabilityTransformation->getJacobian_x_to_u(*Jxu);
    //opserr << *Jxu;
    if (result < 0) {
        opserr << "SearchWithStepSizeAndStepDirection::gradientStandardNormal() - " << endln
               << " could not transform Jacobian from x to u." << endln;
        return -1;
    }
    
    // Gradient in standard normal space
    gradientInStandardNormalSpace->addMatrixTransposeVector(0.0, *Jxu, *gradientInOriginalSpace, 1.0);
    
    // Compute the norm of the gradient in standard normal space
    normOfGradient = gradientInStandardNormalSpace->Norm();
    
    // Check that the norm is not zero
    if (normOfGradient == 0.0) {
        opserr << "SearchWithStepSizeAndStepDirection::gradientStandardNormal() - " << endln
               << " the norm of the gradient is zero. " << endln;
        return -1;
    }
    
    // Compute alpha-vector
    alpha->addVector(0.0, *gradientInStandardNormalSpace, -1.0/normOfGradient );
    
    return 0;
}

    
int
SearchWithStepSizeAndStepDirection::findDesignPoint()
{
    // Declaration of data used in the algorithm
    int numberOfRandomVariables = theReliabilityDomain->getNumberOfRandomVariables();
    int result = 0;
    int evaluationInStepSize = 0;

    double gFunctionValue = 1.0;
    double gFunctionValue_old = 1.0;
    double stepSize = 1.0;

    Vector u_old(numberOfRandomVariables);
    Vector uNew(numberOfRandomVariables);	

    theFunctionEvaluator->initializeNumberOfEvaluations();
    theStepSizeRule->initialize();
    
    Jux->Zero();
    Jxu->Zero();
    
    // get starting x values from parameter directly
    for (int j = 0; j < numberOfRandomVariables; j++) {
        RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(j);
        int param_indx = theReliabilityDomain->getParameterIndexFromRandomVariableIndex(j);
        Parameter *theParam = theOpenSeesDomain->getParameterFromIndex(param_indx);
        
        double rvVal = theRV->getStartValue();
        (*x)(j) = rvVal;

        // now we should update the parameter value
        theParam->update(rvVal);
    }

    // Transform starting point into standard normal space to initialize u vector
    result = theProbabilityTransformation->transform_x_to_u(*u);
        if (result < 0) {
        opserr << "SearchWithStepSizeAndStepDirection::findDesignPoint() - " << endln
            << " could not transform from x to u." << endln;
        return -1;
    }
  
    // Loop to find design point
    steps = 1;
    while (steps <= maxNumberOfIterations) {

        // Evaluate limit-state function unless it has been done in 
        // a trial step by the "stepSizeAlgorithm"
        if (evaluationInStepSize == 0) {
            // set perturbed values in the variable namespace
            if (theFunctionEvaluator->setVariables() < 0) {
                opserr << "ERROR SearchWithStepSizeAndStepDirection -- error setting variables in namespace" << endln;
                return -1;
            }

            // run analysis
            if (theFunctionEvaluator->runAnalysis() < 0) {
                opserr << "ERROR SearchWithStepSizeAndStepDirection -- error running analysis" << endln;
                return -1;
            }

            // evaluate LSF and obtain result
            int lsfTag = theReliabilityDomain->getTagOfActiveLimitStateFunction();
            LimitStateFunction *theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(lsfTag);
            const char *lsfExpression = theLimitStateFunction->getExpression();
            theFunctionEvaluator->setExpression(lsfExpression);

            gFunctionValue_old = gFunctionValue;
            gFunctionValue = theFunctionEvaluator->evaluateExpression();
        }
    
        // Set scale parameter
        if (steps == 1)	{
            Gfirst = gFunctionValue;
            //opserr << "First g-fun value: " << Gfirst << endln;
            if (printFlag != 0) {
                opserr << " Limit-state function value at start point, g=" << gFunctionValue << endln;
                opserr << " STEP #0: ";
            }
            theReliabilityConvergenceCheck->setScaleValue(gFunctionValue);
        }


        // Compute gradient in standard normal variable space
        result = this->gradientStandardNormal(gFunctionValue);
        if (result < 0) {
            opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
                   << " could not compute gradient in standard normal space." << endln;
            return -1;
        }

        // Check convergence
        result = theReliabilityConvergenceCheck->check(*u,gFunctionValue,*gradientInStandardNormalSpace);
        
    
        if (result > 0 || steps == maxNumberOfIterations)  {
            // Inform the user of the happy news!
            opserr << "Design point found!" << endln;

            // Compute the gamma vector
            // Get Jacobian u-space to x-space
            result = theProbabilityTransformation->getJacobian_u_to_x(*u, *Jux);
            if (result < 0) {
                opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
                    << " could not transform from u to x." << endln;
                return -1;
            }

            Vector tempProduct(numberOfRandomVariables);
            tempProduct.addMatrixTransposeVector(0.0, *Jux, *alpha, 1.0);
            
            // Only diagonal elements of (J_xu*J_xu^T) are used
            for (int j = 0; j < numberOfRandomVariables; j++) {
                double sum = 0.0;
                double jk;
                for (int k = 0; k < numberOfRandomVariables; k++) {
                    jk = (*Jxu)(j,k);
                    sum += jk*jk;
                }
                (*gamma)(j) = sqrt(sum) * tempProduct(j);
                //beta += ((*alpha)(j))*((*u)(j));
            }
            gamma->Normalize();		      
      
            Glast = gFunctionValue;
            numberOfEvaluations = theFunctionEvaluator->getNumberOfEvaluations();
            

            // compute sensitivity pf beta wrt to LSF parameters (if they exist)
            // Note this will need to change because now parameters have multiple derived classes
            /*
            const Matrix &dGdPar = theGradientEvaluator->getDgDpar();
            int numPars = dGdPar.noRows();
            if (numPars > 0 && dGdPar.noCols() > 1) {
                opserr << "Sensitivity of beta wrt LSF parameters: ";
                Vector dBetadPar(numPars);

                for (int ib=0; ib<numPars; ib++)
                dBetadPar(ib) = dGdPar(ib,1) / gradientInStandardNormalSpace->Norm();

                opserr << dBetadPar;
            }
            */
            
            // compute the first principal curvature here (ADK and De Stefano 1992)
            // so we do not need to store or provide access to SecondLast variables
            double alphaProduct = *alphaSecondLast ^ *alpha;
            Vector uDiff = *u - *uSecondLast;
            double signumProduct = *alphaSecondLast ^ uDiff;
            
            firstCurvature = acos(alphaProduct) / uDiff.Norm();
            if (signumProduct < 0)
                firstCurvature *= -1.0;
            
            return 1;
        }

        
        // else we need to continue to perform iterations
        
        
        // Store 'u' and 'alpha' at the second last iteration point
        *uSecondLast = *u;
        *alphaSecondLast = *alpha;

        // Let user know that we have to take a new step
        if (printFlag != 0)
            opserr << " STEP #" << steps <<": ";  
        
        // Determine search direction
        result = theSearchDirection->computeSearchDirection(steps,
                        *u, gFunctionValue, *gradientInStandardNormalSpace );
        if (result < 0) {
            opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
                << " could not compute search direction. " << endln;
            return -1;
        }
        *searchDirection = theSearchDirection->getSearchDirection();


        // save the previous displacement before modifying
        u_old = *u;
        
        // Determine step size
        bool continueStepSize = true;
        while (continueStepSize) {
            result = theStepSizeRule->computeStepSize(*u, *gradientInStandardNormalSpace, 
                                                      gFunctionValue, *searchDirection, steps);
            if (result < 0) {  
                // something went wrong
                opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
                    << " could not compute step size. " << endln;
                return -1;
            }
            else if (result == 0) {  
                // step was successful, continue to next iteration
                evaluationInStepSize = 0;
                continueStepSize = false;
            }
            else if (result == 1) {  
                // stepSize requests step size reduction
                evaluationInStepSize = 0;
                continueStepSize = true;
            }
            stepSize = theStepSizeRule->getStepSize();
            
            u->addVector(1.0, *searchDirection, stepSize);
            
            // Transform to x-space
            result = theProbabilityTransformation->transform_u_to_x(*u, *x);
            if (result < 0) {
                opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
                    << " could not transform from u to x." << endln;
                return -1;
            }
            
            // update domain with new x values
            for (int j = 0; j < numberOfRandomVariables; j++) {
                RandomVariable *theRV = theReliabilityDomain->getRandomVariablePtrFromIndex(j);
                int param_indx = theReliabilityDomain->getParameterIndexFromRandomVariableIndex(j);
                Parameter *theParam = theOpenSeesDomain->getParameterFromIndex(param_indx);
                
                // now we should update the parameter value
                theParam->update( (*x)(j) );
            }
        }
        
        
        // Increment the loop parameter
        steps++;

    }


    // Print a message if max number of iterations was reached
    // (Note: in this case the last trial point was never transformed/printed)
    opserr << "Maximum number of iterations was reached before convergence." << endln;

    /*
    if (printFlag != 0)
    outputFile2.close();
    */

    return -1;
}


const Vector&
SearchWithStepSizeAndStepDirection::get_x()
{
	return *x;
}

const Vector &
SearchWithStepSizeAndStepDirection::get_u()
{
	return *u;
}

const Vector&
SearchWithStepSizeAndStepDirection::get_alpha()
{
	return *alpha;
}

const Vector&
SearchWithStepSizeAndStepDirection::get_gamma()
{
    return *gamma;
}

int
SearchWithStepSizeAndStepDirection::getNumberOfSteps()
{
	return (steps-1);
}

double
SearchWithStepSizeAndStepDirection::getFirstCurvature()
{
    return firstCurvature;
}

const Vector&
SearchWithStepSizeAndStepDirection::getLastSearchDirection()
{
	return *searchDirection;
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
	return *gradientInStandardNormalSpace;
}


const Vector&
SearchWithStepSizeAndStepDirection::getGradientInOriginalSpace()
{
	return *gradientInOriginalSpace;
}


int
SearchWithStepSizeAndStepDirection::getNumberOfEvaluations()
{
	return numberOfEvaluations;
}


