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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-04-28 20:51:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SystemAnalysis.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <SystemAnalysis.h>
#include <ReliabilityDomain.h>
#include <ReliabilityAnalysis.h>
#include <LimitStateFunction.h>
#include <MatrixOperations.h>
#include <NormalRV.h>
#include <Vector.h>
#include <Matrix.h>
#include <math.h>
#include <string.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

SystemAnalysis::SystemAnalysis(ReliabilityDomain *passedReliabilityDomain,
							   const char *passedFileName)
:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
	fileName = new char[256];
	strcpy(fileName,passedFileName);
}


SystemAnalysis::~SystemAnalysis()
{
	if (fileName != 0)
		delete [] fileName;
}


double 
SystemAnalysis::getLowerBound(void)
{
	return minLowerBound;
}



double 
SystemAnalysis::getUpperBound(void)
{
	return maxUpperBound;
}



int 
SystemAnalysis::analyze(void)
{

	// Alert the user that the system analysis has started
	opserr << "System Reliability Analysis is running ... " << endln;

	
	// Initial declarations
	double beta;
	double pf1;
	Vector alpha;
	LimitStateFunction *theLimitStateFunction;
	NormalRV aStdNormalRV(1, 0.0, 1.0, 0.0);
	int i, j, k, m, n;

	// Number of limit-state functions
	int numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();

	// Number of random variables
	int nrv = theReliabilityDomain->getNumberOfRandomVariables();

	// Allocate vectors to store ALL the betas and alphas
	Vector allBetas(numLsf);
	Vector allPf1s(numLsf);
	Matrix allAlphas(nrv,numLsf);

	// Loop over number of limit-state functions and collect results
	for (i=0; i<numLsf; i++ ) {

		// Get FORM results from the limit-state function
		theLimitStateFunction = theReliabilityDomain->getLimitStateFunctionPtr(i+1);
		beta = theLimitStateFunction->FORMReliabilityIndexBeta;
		pf1 = theLimitStateFunction->FORMProbabilityOfFailure_pf1;
		alpha = theLimitStateFunction->normalizedNegativeGradientVectorAlpha;

		// Put FORM results into vector of all betas and alphas
		// Note that the first index of 'allAlphas' here denote 'nrv'
		allBetas(i) = beta;
		allPf1s(i) = pf1;
		for (j=0; j<nrv; j++ ) {
			allAlphas(j,i) = alpha(j);
		}
	}


	// Compute vector of 'rhos', that is, dot products between the alphas
	// Note that only the upper diagonal is covered
	Matrix rhos(numLsf,numLsf);
	double dotProduct;
	for (i=0; i<numLsf; i++ ) {
		for (j=i+1; j<numLsf; j++ ) {
			dotProduct = 1.0;
			for (k=0; k<nrv; k++ ) {
				dotProduct=dotProduct*allAlphas(k,i)*allAlphas(k,j);
			}
			rhos(i,j) = dotProduct;
		}
	}


	// Compute the bi-variate normal distribution for all the pairs
	double beta1, beta2, integral;
	double a = 0.0;
	double h, fa, fb, sum_fx2j, sum_fx2j_1;
	int n_2;
	Matrix Pmn(numLsf,numLsf);
	for (m=0; m<numLsf; m++ ) {
		for (n=m+1; n<numLsf; n++ ) {
			beta1 = allBetas(m);
			beta2 = allBetas(n);
			// Composite Simpson's numerical integration
			n_2 = 400; // (half the number of intervals)
			h = rhos(m,n)/(2.0*n_2);
			fa = functionToIntegrate(a,beta1,beta2);
			fb = functionToIntegrate(rhos(m,n),beta1,beta2);
			sum_fx2j = 0.0;
			sum_fx2j_1 = 0.0;
			for (int j=1;  j<=n_2;  j++) {
				sum_fx2j   = sum_fx2j   + functionToIntegrate((a+(j*2.0)*h)    ,beta1,beta2);
				sum_fx2j_1 = sum_fx2j_1 + functionToIntegrate((a+(j*2.0-1.0)*h),beta1,beta2);
			}
			sum_fx2j = sum_fx2j - fb;
			integral = h/3.0*(fa + 2.0*sum_fx2j + 4.0*sum_fx2j_1 + fb);
			Pmn(m,n) = aStdNormalRV.getCDFvalue(-beta1)
				*aStdNormalRV.getCDFvalue(-beta2) + integral;
		}
	}

	// Arrange-envelope
	minLowerBound = 1.0;
	maxUpperBound = 0.0;
	Vector arrangement(numLsf);
	Matrix Pmn_arranged = Pmn;
	Vector allPf1s_arranged = allPf1s;

	for (i=0; i<(factorial(numLsf)); i++) {

		arrangement = arrange(numLsf);

		for (j=0; j<numLsf; j++) {
			allPf1s_arranged(j) = allPf1s(((int)arrangement(j)-1));
		}
		for (j=0; j<numLsf; j++) {
			for (k=0; k<numLsf; k++) {
				Pmn_arranged(j,k) = Pmn(((int)arrangement(j)-1),((int)arrangement(k)-1));
			}
		}

		// Compute lower probability bound
		double lowerBound = allPf1s_arranged(0);
		double temp1, temp2;
		for (m=2; m<=numLsf; m++) {
			temp1 = 0.0;
			for (n=1; n<=n-1; n++) {
				temp1 += Pmn_arranged(m-1,n-1);
			}
			temp2 = allPf1s_arranged(m-1) - temp1;
			if (temp2 < 0.0) {
				temp2 = 0.0;
			}
			lowerBound += temp2;
		}


		// Compute upper probability bound
		double upperBound = allPf1s_arranged(0);
		double maximus;
		for (m=2; m<=numLsf; m++) {
			maximus = 0.0;
			for (n=1; n<=m-1; n++) {
				if (Pmn_arranged(m-1,n-1) > maximus) {
					maximus = Pmn(m-1,n-1);
				}
			}
			upperBound += allPf1s_arranged(m-1) - maximus;
		}

		// Update bounds
		if (lowerBound < minLowerBound) {
			minLowerBound = lowerBound;
		}
		if (upperBound > maxUpperBound) {
			maxUpperBound = upperBound;
		}
	}


	// Print results  (should do this over all defined systems)
	ofstream outputFile( fileName, ios::out );
	
	outputFile << "#######################################################################" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile << "#  SYSTEM ANALYSIS RESULTS                                            #" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile.setf(ios::scientific, ios::floatfield);
	outputFile << "#  Lower probability bound: ........................... " 
		<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<minLowerBound
		<< "  #" << endln;
	outputFile << "#  Upper probability bound: ........................... " 
		<<setiosflags(ios::left)<<setprecision(5)<<setw(12)<<maxUpperBound
		<< "  #" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile << "#  Warning: All permutations of the limit-state functions should      #" << endln;
	outputFile << "#           be checked. Currently, this is not done automatically.    #" << endln;
	outputFile << "#           The user is therefore advised to do this manually.        #" << endln;
	outputFile << "#           Contact the developer for further details/plans.          #" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile << "#######################################################################" << endln << endln << endln;

	outputFile.close();


	// INform the user on screen
	opserr << "System reliability analysis completed." <<endln;

	return 0;
}



double
SystemAnalysis::functionToIntegrate(double rho, double beta1, double beta2)
{
	double pi = 3.14159265358979;
	return 1.0/(2.0*pi*sqrt(1.0-rho*rho)) 
		* exp(-(beta1*beta1+beta2*beta2-2*rho*beta1*beta2)
		/(2.0*(1.0-rho*rho)));
}

int
SystemAnalysis::factorial(int num)
{
	int i = num-1;
	int result = num;

	while (i>0) {
		result = result * i;
		i--;
	}

	return result;
}

Vector
SystemAnalysis::arrange(int num)
{
	Vector arrang(num);

	// This is not yet implemented

	for (int i=0; i<num; i++) {
		arrang(i) = i+1;
	}

	return arrang;
}
