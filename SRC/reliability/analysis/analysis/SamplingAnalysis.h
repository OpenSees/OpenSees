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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-10-27 23:45:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/SamplingAnalysis.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef SamplingAnalysis_h
#define SamplingAnalysis_h

#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <ProbabilityTransformation.h>
#include <RandomNumberGenerator.h>
#include <GFunEvaluator.h>

#include <fstream>
using std::ofstream;

class SamplingAnalysis : public ReliabilityAnalysis
{

public:
	SamplingAnalysis(	ReliabilityDomain *passedReliabilityDomain,
						ProbabilityTransformation *passedProbabilityTransformation,
						GFunEvaluator *passedGFunEvaluator,
						RandomNumberGenerator *passedRandomNumberGenerator,
						int passedNumberOfSimulations,
						double passedTargetCOV,
						double samplingStdv,
						int printFlag,
						TCL_Char *fileName,
						Vector *startPoint,
						int analysisTypeTag);

	~SamplingAnalysis();

	int analyze(void);

protected:
	
private:
	ReliabilityDomain *theReliabilityDomain;
	ProbabilityTransformation *theProbabilityTransformation;
	GFunEvaluator *theGFunEvaluator;
	RandomNumberGenerator *theRandomNumberGenerator;
	int numberOfSimulations;
	double targetCOV;
	double samplingStdv;
	int printFlag;
	char *fileName;
	Vector *startPoint;
	int analysisTypeTag;
};

#endif
