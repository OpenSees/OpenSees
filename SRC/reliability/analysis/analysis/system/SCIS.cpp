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
                                                                        
// $Revision: 1.9 $
// $Date: 2008-08-27 17:17:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/system/SCIS.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <SCIS.h>
#include <SystemAnalysis.h>
#include <Cutset.h>
#include <CutsetIter.h>
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

#include <math.h>
#include <time.h>

SCIS::SCIS(ReliabilityDomain *passedReliabilityDomain, FunctionEvaluator *passedEvaluator, 
           TCL_Char *passedFileName, int aType, 
           TCL_Char *passedBeta, TCL_Char *passedRho, long int passedN, double passedTol)
	:SystemAnalysis(passedReliabilityDomain, passedEvaluator, passedBeta, passedRho)
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
	int numCuts = 1;
	double uBound = 1;
	double lBound = 0;
	int result;
	
	if (analysisType == 0 || analysisType == 1) {
		// compute and get bounds
		result = computeBounds(analysisType);
		if (result != 0)
			opserr << "SCIS::analyze WARNING - failed to compute system bounds" << endln;
		
		uBound = getUpperBound();
		lBound = getLowerBound();
	} else if (analysisType == 2) {
		numCuts = theReliabilityDomain->getNumberOfCutsets();
		if (numCuts < 1)
			opserr << "SCIS::analyze WARNING - not enough cutsets available in domain" << endln;
		
		result = setCutsets();
		if (result != 0)
			opserr << "SCIS::analyze WARNING - failed to set cutset components" << endln;
	}
	
	// perform analysis
	double pf = 0;
	if (analysisType == 0) {
		// parallel system
		pf = SCISfunc(allBetas,rhos,-1.0);
	} else if (analysisType == 1) {
		// series system
		pf = 1.0 - SCISfunc(allBetas,rhos,1.0);
	} else if (analysisType == 2) {
		// general system (k cutset unions of parallel subsystems)
		int ck;
		
		// first-order terms in the inclusion-exclusion rule (always present)
		// upper bound to probability when cutsets are not disjoint
		Cutset *theCutset;
		CutsetIter &cutIter = theReliabilityDomain->getCutsets();
		while ((theCutset = cutIter()) != 0) {
			pf += SCISfunc(theCutset->getBetaCutset(),theCutset->getRhoCutset(),-1.0);
		}
		uBound = pf;
		
		// if necessary, carry out n-choose-k algorithm to find all remaining pairs of events
		// note, this is extremely dependent on ability to check n factorial pairs even with a 
		// greatly reduced decision tree
		for (ck = 2; ck <= numCuts; ck++) {
			int numPerms = getNumPermutations(ck,numCuts);
			result = setPermutations(ck,numCuts);
			if (result != 0) {
				opserr << "SCIS::analyze WARNING - failed to set permutations of cut sets" << endln;
				return -1;
			}
			
			double scaleF = pow(-1.0,ck-1);
			for (int ip = 1; ip <= numPerms; ip++) {
				result = setPermutedComponents(ck,ip-1);
				if (result != 0) {
					opserr << "SCIS::analyze WARNING - failed to set permuted beta and rho" << endln;
					return -1;
				}
				
				pf += scaleF * SCISfunc(getBetaPermutation(),getRhoPermutation(),-1.0);
			}
		}
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
	} else if (analysisType == 1) {
		// series system
		outputFile << "#  Series probability failure estimate (SCIS): ........ "
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf<< "  #" << endln;
		outputFile << "#  Lower series probability bound: .................... " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<lBound<< "  #" << endln;
		outputFile << "#  Upper series probability bound: .................... " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<uBound<< "  #" << endln;
	} else if (analysisType == 2) {
		// general system
		outputFile << "#  General probability failure estimate (SCIS): ....... "
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf<< "  #" << endln;
		outputFile << "#  General upper probability bound: ................... " 
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
	static NormalRV uRV(1, 0.0, 1.0);
	Vector beta(n);
	Matrix rho(n,n);
	int i,ii,j,k,N;
	
	rho = rhoin;
	for (i=0; i < n; i++) {
		beta(i) = allbeta(i)*modifier;
		rho(i,i) = 1.0;
	}
	
	if (n == 1)
		return uRV.getCDFvalue( beta(0) );
	else {
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

		// SCIS algorithm (Ambartzumian and Der Kiureghian)
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
	#ifdef _LINUX
				if ( isfinite(zb) )
					zb = uRV.getCDFvalue(zb);
				else
					zb = 0;
	#endif
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
}
