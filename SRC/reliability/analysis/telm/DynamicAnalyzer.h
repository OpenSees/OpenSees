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
// $Date: 2008-02-29 19:43:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/DynamicAnalyzer.h,v $

#ifndef DynamicAnalyzer_h
#define DynamicAnalyzer_h

#include <Analyzer.h>
#include <ReliabilityDirectIntegrationAnalysis.h>
//#include<Integrator.h>

class DynamicAnalyzer : public Analyzer
{
public:
	DynamicAnalyzer(ReliabilityDomain*, 
					Domain* passedStructuralDomain,
					InitialStaticAnalysis* passedInitialStaticAnalysis,
					ReliabilityDirectIntegrationAnalysis* passedTransientAnalysis,
			                Integrator* passedSensitivityAlgorithm,
				//	SensitivityIntegrator* passedSensitivityIntegrator,
				        Integrator* passedSensitivityIntegrator,//Abbas

					int passednstep,
					double passeddelta,
					int passednumLoadPatterns,
					int* passedDynamicLoadPatterns,
					bool passedprint=false);
	~DynamicAnalyzer();
	int analyze(Vector x);
	int analyzeMean();
	void activateSensitivty(void);
	void inactivateSensitivty(void);
	void setGFunEachStepEvaluator(GFunEachStepEvaluator *pGFunEachStepEvaluator);
	void inactivateGFunEachStepEvaluator();
	virtual Matrix* getEachStepConvFlag();
	virtual Matrix* getEachStepResult();

protected:
	ReliabilityDirectIntegrationAnalysis* theTransientAnalysis;
	GFunEachStepEvaluator* theGFunEachStepEvaluator;
private:
};
#endif
