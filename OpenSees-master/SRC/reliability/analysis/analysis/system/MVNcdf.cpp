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
// $Date: 2008-08-27 17:17:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/system/MVNcdf.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <MVNcdf.h>
#include <SystemAnalysis.h>
#include <Cutset.h>
#include <CutsetIter.h>
#include <ReliabilityDomain.h>
#include <NormalRV.h>
#include <MatrixOperations.h>
#include <RandomNumberGenerator.h>
#include <CStdLibRandGenerator.h>
#include <stdlib.h>
#include <time.h>

#include <fstream>
#include <limits.h>
#include <float.h>
#include <iomanip>
#include <iostream>
#include <math.h>

using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;

MVNcdf::MVNcdf(ReliabilityDomain *passedReliabilityDomain, FunctionEvaluator *passedEvaluator, 
               TCL_Char *passedFileName, int aType, 
               TCL_Char *passedBeta, TCL_Char *passedRho, long int passedN, double passedTol)
	:SystemAnalysis(passedReliabilityDomain, passedEvaluator, passedBeta, passedRho)
{
	strcpy(fileName,passedFileName);
	analysisType = aType;
	
	checkvals(passedN,passedTol);
}


MVNcdf::~MVNcdf()
{
}


void MVNcdf::checkvals(long int passedN, double passedTol) {
	if (passedN < LONG_MAX && passedN > 1)
		Nmax = passedN;
	else if (passedN >= LONG_MAX) {
		opserr << "MVNcdf::MVNcdf WARNING - MVN Nmax input must be smaller than " << LONG_MAX << endln;
		Nmax = LONG_MAX;
	}
	else
	  //Nmax = 2e5;
	  Nmax = 200000;

	if (passedTol > DBL_EPSILON && passedTol < 1)
		errMax = passedTol;
	else if (passedTol <= DBL_EPSILON && passedTol > 0) {
		opserr << "MVNcdf::MVNcdf WARNING - MVN tol input must be greater than " << DBL_EPSILON << endln;
		errMax = DBL_EPSILON*2.0;
	}
	else
		errMax = 1.0e-6;
}

int 
MVNcdf::analyze(void)
{

	// Alert the user that the system analysis has started
	opserr << "System Reliability Analysis (MVNcdf) is running ... " << endln;

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
			opserr << "MVNcdf::analyze WARNING - failed to compute system bounds" << endln;
		
		uBound = getUpperBound();
		lBound = getLowerBound();
	} else if (analysisType == 2) {
		numCuts = theReliabilityDomain->getNumberOfCutsets();
		if (numCuts < 1)
			opserr << "MVNcdf::analyze WARNING - not enough cutsets available in domain" << endln;
		
		result = setCutsets();
		if (result != 0)
			opserr << "MVNcdf::analyze WARNING - failed to set cutset components" << endln;
	}
	
	// perform analysis
	double pf = 0;
	if (analysisType == 0) {
		// parallel system
		pf = MVNcdffunc(allBetas,rhos,-1.0);
	} else if (analysisType == 1) {
		// series system
		pf = 1.0 - MVNcdffunc(allBetas,rhos,1.0);
	} else if (analysisType == 2) {
		// general system (k cutset unions of parallel subsystems)
		int ck;
		
		// first-order terms in the inclusion-exclusion rule (always present)
		// upper bound to probability when cutsets are not disjoint
		Cutset *theCutset;
		CutsetIter &cutIter = theReliabilityDomain->getCutsets();
		while ((theCutset = cutIter()) != 0) {
			pf += MVNcdffunc(theCutset->getBetaCutset(),theCutset->getRhoCutset(),-1.0);
		}
		uBound = pf;
		opserr << "MVN upper bound = " << pf << endln;
		
		// if necessary, carry out n-choose-k algorithm to find all remaining pairs of events
		// note, this is extremely dependent on ability to check n factorial pairs even with a 
		// greatly reduced decision tree
		for (ck = 2; ck <= numCuts; ck++) {
			int numPerms = getNumPermutations(ck,numCuts);
			result = setPermutations(ck,numCuts);
			if (result != 0) {
				opserr << "MVNcdf::analyze WARNING - failed to set permutations of cut sets" << endln;
				return -1;
			}
			
			double scaleF = pow(-1.0,ck-1);
			for (int ip = 1; ip <= numPerms; ip++) {
				result = setPermutedComponents(ck,ip-1);
				if (result != 0) {
					opserr << "MVNcdf::analyze WARNING - failed to set permuted beta and rho" << endln;
					return -1;
				}
				
				pf += scaleF * MVNcdffunc(getBetaPermutation(),getRhoPermutation(),-1.0);
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
		outputFile << "#  Parallel probability failure estimate (MVN): ....... "
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf<< "  #" << endln;
		outputFile << "#  Lower parallel probability bound: .................. " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<lBound<< "  #" << endln;
		outputFile << "#  Upper parallel probability bound: .................. " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<uBound<< "  #" << endln;
	} else if (analysisType == 1) {
		// series system
		outputFile << "#  Series probability failure estimate (MVN): ......... "
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf<< "  #" << endln;
		outputFile << "#  Lower series probability bound: .................... " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<lBound<< "  #" << endln;
		outputFile << "#  Upper series probability bound: .................... " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<uBound<< "  #" << endln;
	} else if (analysisType == 2) {
		// general system
		outputFile << "#  General probability failure estimate (MVN): ........ "
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
MVNcdf::MVNcdffunc(const Vector &allbeta, const Matrix &rhoin, double modifier)
{
	int m = allbeta.Size();
	static NormalRV uRV(1, 0.0, 1.0);
	Vector beta(m);
	Matrix rho(m,m);
	int i,j,result;

	rho = rhoin;
	for (i=0; i < m; i++) {
		beta(i) = allbeta(i)*modifier;
		rho(i,i) = 1.0;
	}

	if (m == 1)
		return uRV.getCDFvalue( beta(0) );
	else {	
		Vector x = beta;
		double ci = 99.9;

		// Cholesky decomposition of correlation matrix
		MatrixOperations theMat(rho);
		result = theMat.computeCholeskyAndItsInverse();
		if (result < 0) {
			opserr << "SystemAnalysis::MVNcdf() - could not compute the Cholesky decomposition and " << endln
				<< "  its inverse for the correlation matrix." << endln;
		}
		Matrix C = theMat.getLowerCholesky();
		double alph = uRV.getInverseCDFvalue(ci/100.0);

		double intsum = 0, varsum = 0, q; 
		long int N = 0;
		
		// d is always zero for integration from -infinity to beta
		// so f = e in Genz
		Vector d(m); Vector e(m); Vector f(m);
		d(1-1) = 0.0;
		e(1-1) = uRV.getCDFvalue(x(1-1) / C(1-1,1-1));
		f(1-1) = e(1-1)-d(1-1);
		
		// initialize random number generator
		RandomNumberGenerator *theRandomNumberGenerator = new CStdLibRandGenerator();
		result = theRandomNumberGenerator->generate_nIndependentUniformNumbers(m-1,0,1,time(NULL));
		if (result < 0) {
			opserr << "SystemAnalysis::MVNcdf() - could not generate random numbers for simulation." << endln;
		}
		Vector w(m-1);

		// simulate
		Vector y(m);
		double err = 2.0 * errMax;
		while ( (err > errMax || N < 5) && (N < Nmax) ) {
			result = theRandomNumberGenerator->generate_nIndependentUniformNumbers(m-1,0,1);
			w = theRandomNumberGenerator->getGeneratedNumbers();

			for (i = 2; i <= m; i++) {
				y(i-1-1) = uRV.getInverseCDFvalue(d(i-1-1) + w(i-1-1) * (e(i-1-1)-d(i-1-1)) );
				q = 0.0;
				for (j = 1; j <= i-1; j++)
					q += C(i-1,j-1)*y(j-1);
					
				d(i-1) = 0.0;
				e(i-1) = uRV.getCDFvalue( (x(i-1) - q)/C(i-1,i-1) );
				f(i-1) = (e(i-1)-d(i-1)) * f(i-1-1);
			}

			N++;
			intsum += f(m-1);
			varsum += f(m-1)*f(m-1);
			err = alph*sqrt( (varsum/N - pow(intsum/N,2))/N );
			//opserr << "#" << N << " has err = " << err << " with alpha=" << alph << endln;
		}

		delete theRandomNumberGenerator;
		
		return intsum/N;
	}
}
