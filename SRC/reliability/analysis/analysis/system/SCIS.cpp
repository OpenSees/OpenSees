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
// $Date: 2007-11-08 20:12:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/system/SCIS.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <SCIS.h>
#include <SystemAnalysis.h>
#include <ReliabilityDomain.h>
#include <NormalRV.h>
#include <RandomNumberGenerator.h>
#include <CStdLibRandGenerator.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <float.h>
#include <limits.h>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

#ifndef isfinite
#define isfinite(x) ((x)-(x) == 0)
#endif

SCIS::SCIS(ReliabilityDomain *passedReliabilityDomain, TCL_Char *passedFileName, int aType, 
		 TCL_Char *passedBeta, TCL_Char *passedRho, long int passedN, double passedTol)
	:SystemAnalysis(passedReliabilityDomain, passedBeta, passedRho)
{
	strcpy(fileName,passedFileName);
	analysisType = aType;
	
	checkvals(passedN,passedTol);
}

SCIS::~SCIS()
{
}


void SCIS::checkvals(long int passedN, double passedTol) {
	// note that Nmax is limited to maximum size of Vector class (currently int)
	if (passedN < INT_MAX && passedN > 1)
		Nmax = (int)passedN;
	else if (passedN >= INT_MAX) {
		opserr << "SCIS::SCIS WARNING - SCIS Nmax input must be smaller than " << INT_MAX << endln;
		Nmax = INT_MAX;
	}
	else
	  //Nmax = 1e4;
	  Nmax = 10000;

	if (passedTol > DBL_EPSILON && passedTol < 1)
		errMax = passedTol;
	else if (passedTol <= DBL_EPSILON && passedTol > 0) {
		opserr << "SCIS::SCIS WARNING - SCIS tol input must be greater than " << DBL_EPSILON << endln;
		errMax = DBL_EPSILON*2.0;
	}
	else
		errMax = 1.0e-6;
}

int 
SCIS::analyze(void)
{

	// Alert the user that the system analysis has started
	opserr << "System Reliability Analysis (SCIS) is running ... " << endln;

	// Allocate beta and rho
	const Vector &allBetas = getBeta();
	const Matrix &rhos = getRho();
	
	// compute and get bounds
	int result = computeBounds(analysisType);
	if (result != 0)
		opserr << "SCIS::analyze WARNING - failed to compute system bounds" << endln;
	
	double uBound = getUpperBound();
	double lBound = getLowerBound();
	
	// perform analysis
	double pf = 0;
	if (analysisType == 0) {
		// parallel system
		pf = SCISfunc(allBetas,rhos,-1.0);
	} else {
		// series system
		pf = 1.0 - SCISfunc(allBetas,rhos,1.0);
	}
	
	// Print results  (should do this over all defined systems)
	ofstream outputFile( fileName, ios::out );
	
	outputFile << "#######################################################################" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile << "#  SYSTEM ANALYSIS RESULTS                                            #" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile.setf(ios::scientific, ios::floatfield);
	
	if (analysisType == 0) {
		// parallel system
		outputFile << "#  Parallel probability failure estimate (SCIS): ...... "
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf<< "  #" << endln;
		outputFile << "#  Lower parallel probability bound: .................. " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<lBound<< "  #" << endln;
		outputFile << "#  Upper parallel probability bound: .................. " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<uBound<< "  #" << endln;
	} else {
		// series system
		outputFile << "#  Series probability failure estimate (SCIS): ........ "
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf<< "  #" << endln;
		outputFile << "#  Lower series probability bound: .................... " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<lBound<< "  #" << endln;
		outputFile << "#  Upper series probability bound: .................... " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<uBound<< "  #" << endln;
	}
	outputFile << "#                                                                     #" << endln;
	outputFile << "#######################################################################" << endln << endln << endln;

	outputFile.close();

	// INform the user on screen
	opserr << "System reliability analysis completed." << endln;

	return 0;
}

double
SCIS::SCISfunc(const Vector &allbeta, const Matrix &rhoin, double modifier)
{
	int n = allbeta.Size();
	static NormalRV uRV(1, 0.0, 1.0, 0.0);
	Vector beta(n);
	Matrix rho(n,n);
	int i,ii,j,k,result,N;
	
	rho = rhoin;
	for (i=0; i < n; i++) {
		beta(i) = allbeta(i)*modifier;
		rho(i,i) = 1.0;
	}
	
	// SCIS allocation, note that Nmax is limited to maximum size of Vector class (currently int)
	srand ( time(NULL) );
	Matrix Dk(n,n);
	Matrix d(n,n);
	Vector y(Nmax);
	Vector x(n);
	double rand01;

	// Preprocessing
	// d : Matrix containing parts of inverse of covariance matrix 
	//     (Note) d(k,1:k) - k-th row of inverse matrix of C(1:k,1:k)
	for (ii = 1; ii <= n; ii++) {
		Matrix tempMat(ii,ii);
		for (i = 0; i < ii; i++) {
			for (j = 0; j < ii; j++)
				tempMat(i,j) = rho(i,j);
		}

		Matrix Cinv(ii,ii);
		tempMat.Invert(Cinv);

		for (i = 0; i < ii; i++) {
			for (j = 0; j < ii; j++)
				Dk(i,j) = Cinv(i,j);
		}
		
		for (j = 1; j <= ii; j++)
			d(ii-1,j-1) = Dk(ii-1,j-1); 
	}

	// ICIS algorithm (Ambartzumian and Der Kiureghian)
	// y    : Product of probabilities
	// sumy : Sum of products y
	// P    : Estimate of probability (sumy/N)
	// em   : Mean of conditional multinormal distribution
	// dd   : Standard deviation of conditional multinormal distribution
	// x    : Sampled value for each random variable
	double sumy = 0;
	double cov = DBL_MAX;
	double P = 0;
	double em, sm, dd, za, zb, zz;

	for (ii = 1; ii <= Nmax; ii++) {
		y(ii-1) = 1;

		// calculate probability of k-th variable
		for (k = 1; k <= n; k++) {
			// mean and std for k-th variable
			dd = 1/sqrt(d(k-1,k-1));
			em = 0;

			// mean of conditional distribution
			if (k != 1) {
				for (j = 1; j <= k-1; j++) {
					sm = d(k-1,j-1)*x(j-1)/d(k-1,k-1);
					em = em - sm;
				}
			}

			// conditional cumulative probability za, our interval goes from -inf to beta so za = 0
			//za = (-DBL_MAX-em)/dd;
			//za = uRV.getCDFvalue(za);
			za = 0;
			
			// conditional cumulative probability zb
			zb = (beta(k-1)-em)/dd;
			if ( isfinite(zb) )
				zb = uRV.getCDFvalue(zb);
			else
				zb = 0;
			y(ii-1) = y(ii-1)*(zb-za);

			// sample a random value x_k
			rand01 = (double)rand()/RAND_MAX;
			zz = uRV.getInverseCDFvalue( (zb-za)*(rand01)+za );
			x(k-1) = dd*zz + em;
		}

		// sum of y up to i-trials
		sumy = sumy + y(ii-1);
		N = ii;
		P = sumy / N;

		// calculation of c.o.v. of estimate P
		if (ii >= 5 && P > 0.0) {
			Vector *vectemp = new Vector(N);
			for (j = 0; j < N; j++)
				(*vectemp)(j) = y(j) - P;
				
			cov = vectemp->Norm()/N/P;
			delete vectemp;
		}
			
		// early exit
		if (cov < errMax)
			break;
	}

	return P;
}
