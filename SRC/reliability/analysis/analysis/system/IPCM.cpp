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
                                                                        
// $Revision: 1.4 $
// $Date: 2007-11-08 20:12:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/system/IPCM.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <IPCM.h>
#include <SystemAnalysis.h>
#include <ReliabilityDomain.h>
#include <NormalRV.h>

#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;


#ifndef isnan
#define isnan(x) ((x)!=(x))
#endif

IPCM::IPCM(ReliabilityDomain *passedReliabilityDomain, TCL_Char *passedFileName, int aType, 
		 TCL_Char *passedBeta, TCL_Char *passedRho)
	:SystemAnalysis(passedReliabilityDomain, passedBeta, passedRho)
{
	strcpy(fileName,passedFileName);
	analysisType = aType;
}

IPCM::~IPCM()
{
}


int 
IPCM::analyze(void)
{

	// Alert the user that the system analysis has started
	opserr << "System Reliability Analysis (IPCM) is running ... " << endln;
	
	// Allocate beta and rho
	const Vector &allBetas = getBeta();
	const Matrix &rhos = getRho();
	
	// compute and get bounds
	int result = computeBounds(analysisType);
	if (result != 0)
		opserr << "IPCM::analyze WARNING - failed to compute system bounds" << endln;
	
	double uBound = getUpperBound();
	double lBound = getLowerBound();
	
	// perform analysis
	double pf = 0;
	if (analysisType == 0) {
		// parallel system
		pf = 1.0 - IPCMfunc(allBetas,rhos,-1.0);
	} else {
		// series system
		pf = IPCMfunc(allBetas,rhos,1.0);
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
		outputFile << "#  Parallel probability failure estimate (IPCM): ...... "
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<pf<< "  #" << endln;
		outputFile << "#  Lower parallel probability bound: .................. " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<lBound<< "  #" << endln;
		outputFile << "#  Upper parallel probability bound: .................. " 
				   << setiosflags(ios::left)<<setprecision(5)<<setw(12)<<uBound<< "  #" << endln;
	} else {
		// series system
		outputFile << "#  Series probability failure estimate (IPCM): ........ "
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
IPCM::IPCMfunc(const Vector &allbeta, const Matrix &rhoin, double modifier)
{
	int n = allbeta.Size();
	static NormalRV uRV(1, 0.0, 1.0, 0.0);
	Vector beta(n);
	Matrix rho(n,n);
	int i,ic,ir,j,k;
	
	rho = rhoin;
	for (i=0; i < n; i++) {
		beta(i) = allbeta(i)*modifier;
		rho(i,i) = beta(i);
	}
	
	double c1, c2, r, a, b, c21, jprob;
	double orig, newg, modval;
	
	// ÑÑÑÑ FIRST CYCLE ÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑ 
	double pdfc1 = uRV.getPDFvalue(rho(1-1,1-1)); 
	double cdfc1 = uRV.getCDFvalue(rho(1-1,1-1)); 
	double A1 = pdfc1/cdfc1;
	double B1 = A1*(rho(1-1,1-1) + A1);
	for (k = 2; k <= n; k++) {
		c1 = -rho(1-1,1-1); 
		c2 = -rho(k-1,k-1); 
		r = rho(1-1,k-1); 
		a = pdfc1/(1.0 - cdfc1); 
		b = a*(c1 + a); 
		c21 = (c2 + r*a)/sqrt(1.0 - r*r*b);
		
		// original
		orig = (1 - cdfc1)*uRV.getCDFvalue(c21);
		// bootstrap with binormal
		newg = twoComponent(c1,c2,r);
		
		jprob = uRV.getCDFvalue(c2) - newg;
		// might have a div by 0 problem here
		modval = 1.0 - jprob/cdfc1;
		if ( isnan(modval) ) {
			opserr << "SystemAnalysis::IPCM WARNING illegal norminv value input (nan)" << endln;
			opserr << "   cdfc1 = " << cdfc1 << " and c2 = " << c2 << endln;
		}
		
		rho(k-1,1-1) = uRV.getInverseCDFvalue(modval); 
	}

	for (ir = 2; ir <= n - 1; ir++) {
		for (ic = ir + 1; ic <= n ; ic++)
			rho(ir-1,ic-1) = (rho(ir-1,ic-1) - rho(1-1,ir-1)*rho(1-1,ic-1)*B1)/sqrt((1 - rho(1-1,ir-1)*rho(1-1,ir-1)*B1)*(1 - rho(1-1,ic-1)*rho(1-1,ic-1)*B1));
	}

	// ÑÑÑÑ- OTHER CYCLES ÑÑÑÑÑÑÑÑÑÑÑÑÑÑÑÐ 
	for (j = 2; j <= n - 1; j++) {
		pdfc1 = uRV.getPDFvalue(rho(j-1,j-2));
		cdfc1 = uRV.getCDFvalue(rho(j-1,j-2)); 
		A1 = pdfc1/cdfc1;
		B1 = A1*(rho(j-1,j-2) + A1);
		for (k = j + 1; k <= n; k++) {
			c1 = -rho(j-1,j-2); 
			c2 = -rho(k-1,j-2); 
			r = rho(j-1,k-1); 
			a = pdfc1/(1.0 - cdfc1); 
			b = a*(c1 + a); 
			c21 = (c2 + r*a)/sqrt(1.0 - r*r*b);
			
			// original
			orig = (1 - cdfc1)*uRV.getCDFvalue(c21);
			// bootstrap with binormal
			newg = twoComponent(c1,c2,r);
			
			jprob = uRV.getCDFvalue(c2) - newg;
			modval = 1.0 - jprob/cdfc1;
			if ( isnan(modval) ) {
				opserr << "SystemAnalysis::IPCM WARNING illegal norminv value input (nan)" << endln;
				opserr << "   cdfc1 = " << cdfc1 << " and c2 = " << c2 << endln;
			}
				
			rho(k-1,j-1) = uRV.getInverseCDFvalue(modval); 
		}
		
		for (ir = j + 1; ir <= n - 1; ir++) {
			for (ic = ir + 1; ic <= n; ic++)
				rho(ir-1,ic-1) = (rho(ir-1,ic-1) - rho(j-1,ir-1)*rho(j-1,ic-1)*B1)/sqrt((1 - rho(j-1,ir-1)*rho(j-1,ir-1)*B1)*(1 - rho(j-1,ic-1)*rho(j-1,ic-1)*B1)); 
		}
	}

	// ÑÑÑÐ Calculate the product of conditional marginals 
	double pf = log(uRV.getCDFvalue(rho(1-1,1-1))); 
	for (i = 2; i<=n; i++)
		pf = pf + log(uRV.getCDFvalue(rho(i-1,i-2)));
		
	// check closed-form solution
	double pcf = 1-twoComponent(beta(0),beta(1),rhoin(1,0));
	//opserr << "pcf = " << pcf << " and IPCM = " << 1-exp(pf) << endln;
	
	if (n == 2)
		return pcf;
	else
		return 1-exp(pf);
}
