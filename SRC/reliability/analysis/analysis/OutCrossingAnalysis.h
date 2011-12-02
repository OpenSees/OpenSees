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
// $Date: 2003-04-28 20:51:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/analysis/OutCrossingAnalysis.h,v $

//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef OutCrossingAnalysis_h
#define OutCrossingAnalysis_h

#include <ReliabilityDomain.h>
#include <GFunEvaluator.h>
#include <GradGEvaluator.h>
#include <FindDesignPointAlgorithm.h>
#include <ReliabilityAnalysis.h>
#include <tcl.h>

#include <fstream>
using std::ofstream;

class OutCrossingAnalysis : public ReliabilityAnalysis
{

public:
	OutCrossingAnalysis(
				ReliabilityDomain *theRelDom,
				GFunEvaluator *theGFunEval,
				GradGEvaluator *theSensEval,
				FindDesignPointAlgorithm *theFindDesPt,
				int analysisType,
				int p_sampleFreq,
				double p_littleDeltaT,
				const char *fileName);
	~OutCrossingAnalysis();

	int analyze(void);

protected:
	
private:
	ReliabilityDomain *theReliabilityDomain;
	GFunEvaluator *theGFunEvaluator;
	GradGEvaluator *theGradGEvaluator;
	FindDesignPointAlgorithm *theFindDesignPointAlgorithm;
	int analysisType;
	int sampleFreq;
	double littleDeltaT;
	char *fileName;

	// A private method:
	double functionToIntegrate(double rho, double beta1, double beta2);

};

#endif
