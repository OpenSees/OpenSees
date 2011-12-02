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
// $Date: 2003-04-28 20:51:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/FORMAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef FORMAnalysis_h
#define FORMAnalysis_h

#include <ReliabilityAnalysis.h>
#include <FindDesignPointAlgorithm.h>
#include <ProbabilityTransformation.h>
#include <ReliabilityDomain.h>

#include <fstream>
using std::ofstream;

class FORMAnalysis : public ReliabilityAnalysis
{

public:
	FORMAnalysis(ReliabilityDomain *passedReliabilityDomain,
				 FindDesignPointAlgorithm *passedFindDesignPointAlgorithm,
				 ProbabilityTransformation *passedProbabilityTransformation,
				 const char *fileName,
				 int relSensTag);
	virtual ~FORMAnalysis();

	int analyze(void);

protected:

private:
	ReliabilityDomain *theReliabilityDomain;
	FindDesignPointAlgorithm *theFindDesignPointAlgorithm;
	ProbabilityTransformation *theProbabilityTransformation;
	char *fileName;
	int relSensTag;
};

#endif
