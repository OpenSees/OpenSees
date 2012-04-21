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
                                                                        
// $Revision: 1.6 $
// $Date: 2008-08-27 17:17:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/system/PCM.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <PCM.h>
#include <SystemAnalysis.h>
#include <Cutset.h>
#include <CutsetIter.h>
#include <CorrelatedStandardNormal.h>
#include <ReliabilityDomain.h>
#include <NormalRV.h>
#include <math.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


PCM::PCM(ReliabilityDomain *passedReliabilityDomain, FunctionEvaluator *passedEvaluator, 
         TCL_Char *passedFileName, int aType, 
		 TCL_Char *passedBeta, TCL_Char *passedRho)
	:SystemAnalysis(passedReliabilityDomain, passedEvaluator, passedBeta, passedRho)

{
	strcpy(fileName,passedFileName);
	analysisType = aType;
}

PCM::~PCM()
{
}


int 
PCM::analyze(void)
{

	// Alert the user that the system analysis has started
	opserr << "System Reliability Analysis (PCM) is running .. " << endln;

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
			opserr << "PCM::analyze WARNING - failed to compute system bounds" << endln;
		
		uBound = getUpperBound();
		lBound = getLowerBound();
	} else if (analysisType == 2) {
		numCuts = theReliabilityDomain->getNumberOfCutsets();
		if (numCuts < 1)
			opserr << "PCM::analyze WARNING - not enough cutsets available in domain" << endln;
		
		result = setCutsets();
		if (result != 0)
			opserr << "PCM::analyze WARNING - failed to set cutset components" << endln;
	}
	
	// perform analysis
	double pf = 0;
	if (analysisType == 0) {
		// parallel system (1 cut set with all n components)
		pf = PCMfunc(allBetas,rhos,-1.0);
	} else if (analysisType == 1) {
		// series system (n cut sets with single component)
		pf = 1.0 - PCMfunc(allBetas,rhos,1.0);
	} else if (analysisType == 2) {
		// general system (k cutset unions of parallel subsystems)
		int ck;
		
		// first-order terms in the inclusion-exclusion rule (always present)
		// upper bound to probability when cutsets are not disjoint
		Cutset *theCutset;
		CutsetIter &cutIter = theReliabilityDomain->getCutsets();
		while ((theCutset = cutIter()) != 0) {
			pf += PCMfunc(theCutset->getBetaCutset(),theCutset->getRhoCutset(),-1.0);
		}
		uBound = pf;
		
		// if necessary, carry out n-choose-k algorithm to find all remaining pairs of events
		// note, this is extremely dependent on ability to check n factorial pairs even with a 
		// greatly reduced decision tree
		for (ck = 2; ck <= numCuts; ck++) {
			int numPerms = getNumPermutations(ck,numCuts);
			result = setPermutations(ck,numCuts);
			if (result != 0) {
				opserr << "PCM::analyze WARNING - failed to set permutations of cut sets" << endln;
				return -1;
			}
			
			double scaleF = pow(-1.0,ck-1);
			for (int ip = 1; ip <= numPerms; ip++) {
				result = setPermutedComponents(ck,ip-1);
				if (result != 0) {
					opserr << "PCM::analyze WARNING - failed to set permuted beta and rho" << endln;
					return -1;
				}
				
				pf += scaleF * PCMfunc(getBetaPermutation(),getRhoPermutation(),-1.0);
			}
		}
	}
	
	// Print results 
	ofstream outputFile( fileName, ios::out );
	
	outputFile << "#######################################################################" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile << "#  SYSTEM ANALYSIS RESULTS                                            #" << endln;
	outputFile << "#                                                                     #" << endln;
	outputFile.setf(ios::scientific, ios::floatfield);
	
	if (analysisType == 0) {
		// parallel system
		outputFile << "#  Parallel probability failure estimate (PCM): ....... "
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf<< "  #" << endln;
		outputFile << "#  Lower parallel probability bound: .................. " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<lBound<< "  #" << endln;
		outputFile << "#  Upper parallel probability bound: .................. " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<uBound<< "  #" << endln;
	} else if (analysisType == 1) {
		// series system
		outputFile << "#  Series probability failure estimate (PCM): ......... "
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf<< "  #" << endln;
		outputFile << "#  Lower series probability bound: .................... " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<lBound<< "  #" << endln;
		outputFile << "#  Upper series probability bound: .................... " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<uBound<< "  #" << endln;
	} else if (analysisType == 2) {
		// general system
		outputFile << "#  General probability failure estimate (PCM): ........ "
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
PCM::PCMfunc(const Vector &allbeta, const Matrix &rhoin, double modifier)
{
	int n = allbeta.Size();
	static NormalRV uRV(1, 0.0, 1.0);
	Vector beta(n);
	Matrix rho(n,n);
	int i,ic,ir,j,k;
	
	rho = rhoin;
	for (i=0; i < n; i++) {
		beta(i) = allbeta(i)*modifier;
		rho(i,i) = beta(i);
	}
	
	if (n == 1)
		return uRV.getCDFvalue( beta(0) );
	else if (n == 2) {
		// check closed-form solution
		CorrelatedStandardNormal phi2(rhoin(1,0));
		double pcf = phi2.getCDF(beta(0),beta(1));
		//opserr << "pcf = " << pcf << " and PCM = " << exp(pf) << endln;
		return pcf;
	} else {
	
		// ÑÑÑÑ FIRST CYCLE ÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑ 
		double A1 = uRV.getPDFvalue(rho(1-1,1-1))/(uRV.getCDFvalue(rho(1-1,1-1)));
		double B1 = A1*(rho(1-1,1-1) + A1);
		for (k = 2; k <= n; k++) 
			rho(k-1,1-1) = (rho(k-1,k-1) + rho(1-1,k-1)*A1)/sqrt(1 - rho(1-1,k-1)*rho(1-1,k-1)*B1);

		for (ir = 2; ir <= n - 1; ir++) {
			for (ic = ir + 1; ic <= n; ic++)
				rho(ir-1,ic-1) = (rho(ir-1,ic-1) - rho(1-1,ir-1)*rho(1-1,ic-1)*B1)/sqrt((1-rho(1-1,ir-1)*rho(1-1,ir-1)*B1)*(1 - rho(1-1,ic-1)*rho(1-1,ic-1)*B1));
		}

		// ÑÑÑÑ- OTHER CYCLES ÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÐ 
		for (j = 2; j <= n - 1; j++) {
			A1 = uRV.getPDFvalue(rho(j-1,j-2))/(uRV.getCDFvalue(rho(j-1,j-2)));
			B1 = A1*(rho(j-1,j-2) + A1);
			for (k = j + 1; k <= n; k++) 
				rho(k-1,j-1) = (rho(k-1,j-2) + rho(j-1,k-1)*A1)/sqrt(1 - rho(j-1,k-1)*rho(j-1,k-1)*B1);	
			
			for (ir = j + 1; ir <= n - 1; ir++) {
				for (ic = ir + 1; ic <= n; ic++) 
					rho(ir-1,ic-1) = (rho(ir-1,ic-1) - rho(j-1,ir-1)*rho(j-1,ic-1)*B1)/sqrt((1 - rho(j-1,ir-1)*rho(j-1,ir-1)*B1)*(1 - rho(j-1,ic-1)*rho(j-1,ic-1)*B1));
			}
		}

		// ÑÑÑÐ Calculate the product of conditional marginals 
		double pf = log(uRV.getCDFvalue(rho(1-1,1-1))); 
		for (i = 2; i<=n; i++)
			pf = pf + log(uRV.getCDFvalue(rho(i-1,i-2)));
		
		return exp(pf);
	}
}
