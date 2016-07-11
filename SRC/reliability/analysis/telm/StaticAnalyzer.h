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

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/StaticAnalyzer.h,v $
                                                                        
#ifndef StaticAnalyzer_h
#define StaticAnalyzer_h

#include <Analyzer.h>
#include <ReliabilityStaticAnalysis.h>

class StaticAnalyzer : public Analyzer
{
public:
	StaticAnalyzer(ReliabilityDomain*, 
					Domain* passedStructuralDomain,
					InitialStaticAnalysis* passedInitialStaticAnalysis,
					ReliabilityStaticAnalysis* passedStaticAnalysis,
				//	SensitivityAlgorithm* passedSensitivityAlgorithm,
				        Integrator *passedSensitivityAlgorithm,
				     //   SensitivityIntegrator *passedSensitivityIntegrator,//Abbas
				       Integrator *passedSensitivityIntegrator,//Abbas

					int passednstep,
					double passeddelta,
					int passednumLoadPatterns,
					int* passedStaticLoadPatterns,
					bool passedprint=false);
	~StaticAnalyzer();
	int analyze(Vector x);
	int analyzeMean();
	void activateSensitivty(void);
	void inactivateSensitivty(void);
	void setGFunEachStepEvaluator(GFunEachStepEvaluator *pGFunEachStepEvaluator);
	void inactivateGFunEachStepEvaluator();
	Matrix* getEachStepResult();
	Matrix* getEachStepConvFlag();


protected:
	ReliabilityStaticAnalysis* theStaticAnalysis;
private:
};
#endif
