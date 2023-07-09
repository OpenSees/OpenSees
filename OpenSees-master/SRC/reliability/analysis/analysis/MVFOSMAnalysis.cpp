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
                                                                        
// $Revision: 1.3 $
// $Date: 2006-12-06 22:32:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/MVFOSMAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <MVFOSMAnalysis.h>
#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <RandomVariablePositioner.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <Matrix.h>
#include <Vector.h>
#include <tcl.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


MVFOSMAnalysis::MVFOSMAnalysis(ReliabilityDomain *passedReliabilityDomain,
							   GFunEvaluator *passedGFunEvaluator,
							   GradGEvaluator *passedGradGEvaluator,
							   Tcl_Interp *passedTclInterp,
							   const char *passedFileName)
:ReliabilityAnalysis()
{
	theReliabilityDomain	= passedReliabilityDomain;
	theGFunEvaluator		= passedGFunEvaluator;
	theGradGEvaluator = passedGradGEvaluator;
	theTclInterp			= passedTclInterp;
	strcpy(fileName,passedFileName);
}


MVFOSMAnalysis::~MVFOSMAnalysis()
{

}



int 
MVFOSMAnalysis::analyze(void)
{

	// Alert the user that the FORM analysis has started
	opserr << "MVFOSM Analysis is running ... " << endln;


	// Initial declarations
	int i,j,k;


	// Open output file
	ofstream outputFile( fileName, ios::out );


	// Get number of random variables 
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();

	
	// Get mean point
	RandomVariable *aRandomVariable;
	Vector meanVector(nrv);
	for (i=1; i<=nrv; i++)
	{
		aRandomVariable = theReliabilityDomain->getRandomVariablePtr(i);
		if (aRandomVariable == 0) {
			opserr << "MVFOSMAnalysis::analyze() -- Could not find" << endln
				<< " random variable with tag #" << i << "." << endln;
			return -1;
		}
		meanVector(i-1) = aRandomVariable->getMean();
	}


	// Establish vector of standard deviations
	Vector stdvVector(nrv);
	for (i=1; i<=nrv; i++)
	{
		aRandomVariable = theReliabilityDomain->getRandomVariablePtr(i);
		stdvVector(i-1) = aRandomVariable->getStdv();
	}


	// Evaluate limit-state function
	int result;
	result = theGFunEvaluator->runGFunAnalysis(meanVector);
	if (result < 0) {
		opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
			<< " could not run analysis to evaluate limit-state function. " << endln;
		return -1;
	}


	// Establish covariance matrix
	Matrix covMatrix(nrv,nrv);
	for (i=1; i<=nrv; i++) {
		covMatrix(i-1,i-1) = stdvVector(i-1)*stdvVector(i-1);
	}
	int ncorr = theReliabilityDomain->getNumberOfCorrelationCoefficients();
	CorrelationCoefficient *theCorrelationCoefficient;
	double covariance, correlation;
	int rv1, rv2;
	for (i=1 ; i<=ncorr ; i++) {
		theCorrelationCoefficient = theReliabilityDomain->getCorrelationCoefficientPtr(i);
		correlation = theCorrelationCoefficient->getCorrelation();
		rv1 = theCorrelationCoefficient->getRv1();
		rv2 = theCorrelationCoefficient->getRv2();
		covariance = correlation*stdvVector(rv1-1)*stdvVector(rv2-1);
		covMatrix(rv1-1,rv2-1) = covariance;
		covMatrix(rv2-1,rv1-1) = covariance;
	}


	// 'Before loop' declarations
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();
	Vector gradient(nrv);
	Matrix matrixOfGradientVectors(nrv,numLsf);
	Vector meanEstimates(numLsf);
	Vector responseStdv(numLsf);
	double responseVariance;
	

	// Loop over limit-state functions
	for (int lsf=1; lsf<=numLsf; lsf++ ) {

		// Inform the user which limit-state function is being evaluated
		opserr << "Limit-state function number: " << lsf << endln;


		// Set tag of active limit-state function
		theReliabilityDomain->setTagOfActiveLimitStateFunction(lsf);


		// Get limit-state function value (=estimation of the mean)
		result = theGFunEvaluator->evaluateG(meanVector);
		if (result < 0) {
			opserr << "SearchWithStepSizeAndStepDirection::doTheActualSearch() - " << endln
				<< " could not tokenize limit-state function. " << endln;
			return -1;
		}
		meanEstimates(lsf-1) = theGFunEvaluator->getG();


		// Evaluate (and store) gradient of limit-state function
		result = theGradGEvaluator->evaluateGradG(meanEstimates(lsf-1),meanVector);
		if (result < 0) {
			opserr << "MVFOSMAnalysis::analyze() -- could not" << endln
				<< " compute gradients of the limit-state function. " << endln;
			return -1;
		}
		gradient = theGradGEvaluator->getGradG();
		for (i=1 ; i<=nrv ; i++) {
			matrixOfGradientVectors(i-1,lsf-1) = gradient(i-1);
		}


		// Estimate of standard deviation
		responseVariance = (covMatrix^gradient)^gradient;
		if (responseVariance <= 0.0) {
			opserr << "ERROR: Response variance of limit-state function number "<< lsf << endln
				<< " is zero! " << endln;
		}
		else {
			responseStdv(lsf-1) = sqrt(responseVariance);
		}

	
		// Print MVFOSM results to the output file
		outputFile << "#######################################################################" << endln;
		outputFile << "#  MVFOSM ANALYSIS RESULTS, LIMIT-STATE FUNCTION NUMBER "
			<<setiosflags(ios::left)<<setprecision(1)<<setw(4)<<lsf <<"          #" << endln;
		outputFile << "#                                                                     #" << endln;
		
		outputFile << "#  Estimated mean: .................................... " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<meanEstimates(lsf-1) 
			<< "  #" << endln;
		outputFile << "#  Estimated standard deviation: ...................... " 
			<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<responseStdv(lsf-1) 
			<< "  #" << endln;
		outputFile << "#                                                                     #" << endln;
		outputFile << "#######################################################################" << endln << endln << endln;


		// Inform the user that we are done with this limit-state function
		opserr << "Done analyzing limit-state function " << lsf << endln;	
	}


	// Estimation of response covariance matrix
	Matrix responseCovMatrix(numLsf,numLsf);
	double responseCovariance;
	Vector gradientVector1(nrv), gradientVector2(nrv);
	for (i=1; i<=numLsf; i++) {
		for (k=0; k<nrv; k++) {
			gradientVector1(k) = matrixOfGradientVectors(k,i-1);
		}
		for (j=i+1; j<=numLsf; j++) {
			for (k=0; k<nrv; k++) {
				gradientVector2(k) = matrixOfGradientVectors(k,j-1);
			}
			responseCovariance = (covMatrix^gradientVector1)^gradientVector2;
			responseCovMatrix(i-1,j-1) = responseCovariance;
		}
	}
	for (i=1; i<=numLsf; i++) {
		for (j=1; j<i; j++) {
			responseCovMatrix(i-1,j-1) = responseCovMatrix(j-1,i-1);
		}
	}


	// Corresponding correlation matrix
	Matrix correlationMatrix(numLsf,numLsf);
	for (i=1; i<=numLsf; i++) {
		for (j=i+1; j<=numLsf; j++) {
			correlationMatrix(i-1,j-1) = responseCovMatrix(i-1,j-1)/(responseStdv(i-1)*responseStdv(j-1));
		}
	}
	for (i=1; i<=numLsf; i++) {
		for (j=1; j<i; j++) {
			correlationMatrix(i-1,j-1) = correlationMatrix(j-1,i-1);
		}
	}

	
	// Print correlation results
	outputFile << "#######################################################################" << endln;
	outputFile << "#  RESPONSE CORRELATION COEFFICIENTS                                  #" << endln;
	outputFile << "#                                                                     #" << endln;
	if (numLsf <=1) {
		outputFile << "#  Only one limit-state function!                                     #" << endln;
	}
	else {
		outputFile << "#   gFun   gFun     Correlation                                       #" << endln;
		outputFile.setf(ios::fixed, ios::floatfield);
		for (i=0; i<numLsf; i++) {
			for (j=i+1; j<numLsf; j++) {
//				outputFile.setf(ios::fixed, ios::floatfield);
				outputFile << "#    " <<setw(3)<<(i+1)<<"    "<<setw(3)<<(j+1)<<"     ";
				if (correlationMatrix(i,j)<0.0) { outputFile << "-"; }
				else { outputFile << " "; }
//				outputFile.setf(ios::scientific, ios::floatfield);
				outputFile <<setprecision(7)<<setw(11)<<fabs(correlationMatrix(i,j));
				outputFile << "                                      #" << endln;
			}
		}
	}
	outputFile << "#                                                                     #" << endln;
	outputFile << "#######################################################################" << endln << endln << endln;



	// Print summary of results to screen (more here!!!)
	opserr << "MVFOSMAnalysis completed." << endln;


	return 0;
}

