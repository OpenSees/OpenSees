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
// $Date: 2003-10-27 23:45:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SystemAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef SystemAnalysis_h
#define SystemAnalysis_h

#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>

#include <fstream>
using std::ofstream;

class SystemAnalysis : public ReliabilityAnalysis
{

public:
	SystemAnalysis(ReliabilityDomain *passedReliabilityDomain,
				   TCL_Char *fileName);
	~SystemAnalysis();

	int		analyze(void);
	double	getLowerBound();
	double	getUpperBound();

protected:

private:
	double functionToIntegrate(double rho, double beta1, double beta2);
	int factorial(int num);
	Vector arrange(int num);

	double minLowerBound;
	double maxUpperBound;
	ReliabilityDomain *theReliabilityDomain;
	char *fileName;
};

#endif
