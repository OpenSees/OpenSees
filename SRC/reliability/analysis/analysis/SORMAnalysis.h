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
// $Date: 2006-12-06 22:32:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SORMAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef SORMAnalysis_h
#define SORMAnalysis_h

#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <FunctionEvaluator.h>
#include <FORMAnalysis.h>
#include <FindCurvatures.h>

#include <fstream>
using std::ofstream;

class SORMAnalysis : public ReliabilityAnalysis
{

public:
	SORMAnalysis(ReliabilityDomain *passedReliabilityDomain,
                 FunctionEvaluator *passedEvaluator,
                 FORMAnalysis *passedFORM,
                 FindCurvatures *passedCurvaturesAlgorithm,
                 TCL_Char *fileName);
	~SORMAnalysis();

	int analyze(void);

protected:

private:
	ReliabilityDomain *theReliabilityDomain;
    FunctionEvaluator *theFunctionEvaluator;
    FORMAnalysis *theFORManalysis;
	FindCurvatures *theCurvaturesAlgorithm;
	char fileName[256];
};

#endif
