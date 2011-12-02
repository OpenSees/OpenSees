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
                                                                        
// $Revision: 1.15 $
// $Date: 2008-04-10 18:10:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SystemAnalysis.cpp,v $


//
// Written by Kevin Mackie (kmackie@mail.ucf.edu)
//

#include <SystemAnalysis.h>
#include <ReliabilityDomain.h>
#include <ReliabilityAnalysis.h>
#include <LimitStateFunction.h>
#include <LimitStateFunctionIter.h>
#include <RandomNumberGenerator.h>
#include <CStdLibRandGenerator.h>
#include <NormalRV.h>

#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;

SystemAnalysis::SystemAnalysis(ReliabilityDomain *passedReliabilityDomain, TCL_Char *passedBeta, TCL_Char *passedRho)
	:ReliabilityAnalysis()
{
	theReliabilityDomain = passedReliabilityDomain;
	strcpy(betaFile,passedBeta);
	strcpy(rhoFile,passedRho);
	
	minLowerBound = 1.0;
	maxUpperBound = 0.0;
	
	int result = initialize();
	if (result < 0) {
		opserr << "SystemAnalysis::SystemAnalysis() ERROR - SystemAnalysis failed to initialize" << endln;
		exit(-1);
	}
}

SystemAnalysis::~SystemAnalysis()
{
	if (allBetas != 0)
		delete allBetas;
	if (rhos != 0 )
		delete rhos;
	if (allPf1s != 0)
		delete allPf1s;
	if (Pmn != 0 )
		delete Pmn;
}


int
SystemAnalysis::initialize()
{
	// Initial declarations
	double beta;
	double pf1;
	int i, j, k, m, n;
	
	if (strlen(betaFile) > 1 && strlen(rhoFile) > 1) {
		// read data for beta and rho from file, do not take from reliability domain
		ifstream inputBeta( betaFile, ios::in );
		if (!inputBeta) {
			opserr << "SystemAnalysis::initialize ERROR - could not open beta file " << betaFile << endln;
			return -1;
		}
		ifstream inputRho( rhoFile, ios::in );
		if (!inputRho) {
			opserr << "SystemAnalysis::initialize ERROR - could not open rho file " << rhoFile << endln;
			return -1;
		}
		
		numLsf = 0;
		while (inputBeta >> beta)
			numLsf++;
		inputBeta.clear();
		inputBeta.seekg(0);
		
		// allocate arrays
		allBetas = new Vector(numLsf);
		allPf1s = new Vector(numLsf);
		rhos = new Matrix(numLsf,numLsf);
		NormalRV uRV(1, 0.0, 1.0, 0.0);
		
		// read data from file into arrays
		for (i=0; i < numLsf; i++ ) {
			inputBeta >> (*allBetas)(i);
			(*allPf1s)(i) = 1.0 - uRV.getCDFvalue( (*allBetas)(i) );
			
			for (j=0; j < numLsf; j++ ) {
				if (!inputRho.eof() )
					inputRho >> (*rhos)(i,j);
				else {
					opserr << "SystemAnalysis::initialize ERROR - rho file does not contain the same data size as the beta file"
						   << endln;
					return -1;
				}
			}
		}
		
		inputBeta.close();
		inputRho.close();
		
	} else {
		// proceed with reliability domain data, requires previous FORM analysis objects

		// Number of limit-state functions
		numLsf = theReliabilityDomain->getNumberOfLimitStateFunctions();

		// Number of random variables
		int nrv = theReliabilityDomain->getNumberOfRandomVariables();

		// Allocate vectors to store ALL the betas and alphas
		allBetas = new Vector(numLsf);
		allPf1s = new Vector(numLsf);
		Matrix allAlphas(nrv,numLsf);

		// Loop over number of limit-state functions and collect results
		LimitStateFunction *theLimitStateFunction;
		LimitStateFunctionIter &lsfIter = theReliabilityDomain->getLimitStateFunctions();
		//for (i=0; i<numLsf; i++ ) {
		while ((theLimitStateFunction = lsfIter()) != 0) {

		  int tag = theLimitStateFunction->getTag();
		  //int i = theLimitStateFunction->getIndex();
		  int i = theReliabilityDomain->getLimitStateFunctionIndex(tag);

		  beta = theLimitStateFunction->getFORM_beta();
		  pf1 = theLimitStateFunction->getFORM_pf1();
		  const Vector &alpha = theLimitStateFunction->getFORM_alpha();

			// Put FORM results into vector of all betas and alphas
			// Note that the first index of 'allAlphas' here denote 'nrv'
			(*allBetas)(i) = beta;
			(*allPf1s)(i) = pf1;
			for (j=0; j<nrv; j++ ) {
				allAlphas(j,i) = alpha(j);
			}
		}
		
		// Compute vector of 'rhos', that is, dot products between the alphas
		rhos = new Matrix(numLsf,numLsf);
		
		double dotProduct;
		for (i=0; i<numLsf; i++ ) {
			for (j=i+1; j<numLsf; j++ ) {
				dotProduct = 0.0;
				for (k=0; k<nrv; k++ ) {
					dotProduct += allAlphas(k,i)*allAlphas(k,j);
				}
				(*rhos)(i,j) = dotProduct;
				(*rhos)(j,i) = dotProduct;
			}
		}
	}

	opserr << "B vector:" << *allBetas;
	opserr << "R matrix:" << *rhos;
	
	// Compute the bi-variate parallel probability for all the pairs
	// Note this is upper-diagonal only
	Pmn = new Matrix(numLsf,numLsf);
	for (m=0; m<numLsf; m++ ) {
		for (n=m+1; n<numLsf; n++ )
			(*Pmn)(m,n) = twoComponent(-(*allBetas)(m),-(*allBetas)(n),(*rhos)(m,n));
	}
	
	return 0;
}


int
SystemAnalysis::computeBounds(int aType)
{	
	double estimate, lowerBound, upperBound;
	int i, j, result;
			
	if (aType == 0) {
		// parallel system
		
		// uni-component lower bound from Boole parallel system
		estimate = 0;
		for (i=0; i<numLsf; i++)
			estimate += (*allPf1s)(i);
		estimate -= numLsf-1;
		minLowerBound = estimate;

		// uni-component upper bound from Boole parallel system
		upperBound = (*allPf1s)(0);
		for (i=1; i<numLsf; i++) {
			if ((*allPf1s)(i) < upperBound)
				upperBound = (*allPf1s)(i);
		}
		maxUpperBound = upperBound;

		if (minLowerBound < 0.0)
			minLowerBound = 0;
		if (maxUpperBound > 1.0)
			maxUpperBound = 1;
		
	} else if (aType == 1) {
		// series system

		// use simulation approach to limit factorial possibilities for the ordering
		RandomNumberGenerator *theRandomNumberGenerator = new CStdLibRandGenerator();
		result = theRandomNumberGenerator->generate_nIndependentUniformNumbers(numLsf,0,numLsf-1,time(NULL));
		if (result < 0)
			opserr << "SystemAnalysis::computeBounds() - could not generate random numbers for simulation." << endln;

		long int numPerms = factorial(numLsf+1);
		int multiplicationFactor = 500;
		if (numPerms > multiplicationFactor*numLsf)
			numPerms = multiplicationFactor*numLsf;
		ID permute(numLsf);
		
		for (long int arr = 1; arr <= numPerms; arr++ ) {
			arrange(numLsf,theRandomNumberGenerator,permute);
			
			// bi-component lower bound from KHD series system
			lowerBound = (*allPf1s)(permute(0));
			for (i=2; i <= numLsf; i++) {
				estimate = (*allPf1s)(permute(i-1));
				
				for (j=1; j <= i-1; j++)
					estimate -= (*Pmn)(permute(i-1),permute(j-1));
				
				if (estimate > 0)
					lowerBound += estimate;
			}

			// bi-component upper bound from KHD series system
			upperBound = (*allPf1s)(permute(0));
			for (i=2; i<= numLsf; i++) {
				estimate = 0;
				
				for (j=1; j < i; j++) {
					if ((*Pmn)(permute(i-1),permute(j-1)) > estimate)
						estimate = (*Pmn)(permute(i-1),permute(j-1));
				}
				
				upperBound += (*allPf1s)(permute(i-1)) - estimate;
			}

			// Update bounds
			if (lowerBound < minLowerBound && lowerBound >= 0.0)
				minLowerBound = lowerBound;
			if (upperBound > maxUpperBound && upperBound <= 1.0)
				maxUpperBound = upperBound;
		}
		
		delete theRandomNumberGenerator;
	
	} else {
	
	}
	
	return 0;
		
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
SystemAnalysis::getNumberLimitStateFunctions(void)
{
	return numLsf;
}

const Vector &
SystemAnalysis::getBeta(void)
{
	return *allBetas;
}

const Matrix &
SystemAnalysis::getRho(void)
{
	return *rhos;
}

double
SystemAnalysis::functionToIntegrate(double rho, double beta1, double beta2)
{
	static const double pi = acos(-1.0);
	if ( 1 - fabs(rho) < DBL_EPSILON ) {
		opserr << "SystemAnalysis::functionToIntegrate warning rho approx. 1 or -1" << endln;
		return 0.0;
	}
	else
		return 1.0/(2.0*pi*sqrt(1.0-rho*rho)) *
			exp(-(beta1*beta1+beta2*beta2-2.0*rho*beta1*beta2)/(2.0*(1.0-rho*rho)));
}

long int
SystemAnalysis::factorial(int num)
{
	int i = num-1;
	long int result = num;

	while (i > 0) {
		result = result * i;
		i--;
	}

	return result;
}

int
SystemAnalysis::arrange(int num, RandomNumberGenerator *randNum, ID &permutation)
{
	// randomly arranges the input pf and Pmn matrix entries
	int result;
	
	Vector index(num);
	
	randNum->generate_nIndependentUniformNumbers(num,0,num-1);
	const Vector &randomArray = randNum->getGeneratedNumbers();

	for (int i=0; i<num; i++) {
	  result = int(randomArray(i) + 0.5); // rounds to nearest int
		if (index(result) > 0) {
			// repeat found
			for (int j=0; j<num; j++) {
				if (index(j) == 0.0) {
					result = j;
					break;
				}
			}
		} 
		// unique assignment
		permutation(i) = result;
		index(result) = 1;
	}
	
	return 0;
}

double SystemAnalysis::twoComponent(double beta1, double beta2, double rho2)
{
	static NormalRV uRV(1, 0.0, 1.0, 0.0);
	
	double thresh = 0.99;
	double integral = 0;
	
	if ( rho2 > thresh ) {
		// create a refined domain near +1
		integral = Simpson(0,thresh,beta1,beta2,rho2);
		integral += Simpson(thresh,rho2,beta1,beta2,rho2);
	} else if ( rho2 < -thresh ) {
		// create a refined domain near -1
		integral = Simpson(0,-thresh,beta1,beta2,rho2);
		integral += Simpson(-thresh,rho2,beta1,beta2,rho2);
	} else {
		integral = Simpson(0,rho2,beta1,beta2,rho2);
	}
	
	return uRV.getCDFvalue(beta1)*uRV.getCDFvalue(beta2) + integral;
}

double
SystemAnalysis::Simpson(double a, double b, double beta1, double beta2, double rho)
{
	// Composite Simpson's numerical integration
	// n must be even, endpoints a and b
	int n = 1600;
	double h, fa, fb;
	double integral = 0;
	
	h = (b-a)/n;
	fa = functionToIntegrate(a,beta1,beta2);
	fb = functionToIntegrate(b,beta1,beta2);
	
	double sum_fx2j_1 = 0.0;
	double sum_fx2j_2 = 0.0;
	for (int j=2; j <= n/2;  j++)
		sum_fx2j_1 = sum_fx2j_1 + functionToIntegrate(a+h*(2*j-2),beta1,beta2);
	for (int j=1; j <= n/2;  j++)
		sum_fx2j_2 = sum_fx2j_2 + functionToIntegrate(a+h*(2*j-1),beta1,beta2);

	integral = h/3.0*(fa + 2.0*sum_fx2j_1 + 4.0*sum_fx2j_2 + fb);

	return integral;
}
