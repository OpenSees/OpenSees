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
                                                                        
// $Revision: 1.2 $
// $Date: 2006-12-06 22:32:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/OptimizationAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef OptimizationAnalysis_h
#define OptimizationAnalysis_h

#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>

#include <fstream>
using std::ofstream;

class OptimizationAnalysis : public ReliabilityAnalysis
{

public:
	OptimizationAnalysis(ReliabilityDomain *passedReliabilityDomain, TCL_Char *fileName);
	virtual ~OptimizationAnalysis();

	int analyze(void);

protected:

private:
	double f0(Vector x, double assump, double betavar);
	Vector fjs(Vector x, double assump, double betavar, double betavarmax, double betavarmin);

	ReliabilityDomain *theReliabilityDomain;
	char fileName[256];
};

#endif
