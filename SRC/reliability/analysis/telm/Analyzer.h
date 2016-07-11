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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/Analyzer.h,v $

#ifndef Analyzer_h
#define Analyzer_h

#include <Domain.h>
#include <ReliabilityDomain.h>
#include <TaggedObjectStorage.h>
#include <ArrayOfTaggedObjects.h>
#include <LoadPatternIter.h>
#include <LoadPattern.h>
#include <InitialStaticAnalysis.h>
#include <SensitivityAlgorithm.h>
#include <Integrator.h>
#include<SensitivityIntegrator.h>
#include <GFunEachStepEvaluator.h>
#include <Node.h>
#include <NodeIter.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;

class Analyzer
{
public:
	Analyzer(ReliabilityDomain* passedReliabilityDomain,
		     Domain* passedStructuralDomain,
		     InitialStaticAnalysis* passedInitialStaticAnalysis,
		     //SensitivityAlgorithm* passedSensitivityAlgorithm,
		     Integrator* passedSensitivityAlgorithm,//Abbas
		     //SensitivityIntegrator *passedSensitivityIntegrator,
		     Integrator *passedSensitivityIntegrator,//Abbas

		     int passedNumstep,
		     double passeddelta,
		     int passedNumLoadPatterns,
		     int* passedLoadPatterns,
		     bool passedprint);
	virtual ~Analyzer();
	virtual int analyze(Vector x)=0;
	virtual int analyzeMean()=0;
	virtual void activateSensitivty(void)=0;
	virtual void inactivateSensitivty(void)=0;
	void saveLoads();
	void modifyLoads();
	void recoverLoads();
	void printresults();
	void setNsteps(int);
	double getDt(void);
	int getNstep(void);
	virtual void setGFunEachStepEvaluator(GFunEachStepEvaluator *pGFunEachStepEvaluator)=0;
	virtual void inactivateGFunEachStepEvaluator()=0;
	virtual Matrix* getEachStepResult()=0;
	virtual Matrix* getEachStepConvFlag()=0;

protected:
	ReliabilityDomain*	theReliabilityDomain;
	Domain* theDomain;
	InitialStaticAnalysis* theInitialStaticAnalysis;
//	SensitivityAlgorithm* theSensitivityAlgorithm;
	Integrator* theSensitivityAlgorithm;

      // SensitivityIntegrator *theSensitivityIntegrator;
       Integrator *theSensitivityIntegrator;//Abbas

       bool activeSensitivity;
	bool print;
	ofstream output;

	int numOrgPatterns;
	LoadPattern** theOrgPatterns;

//	LoadPatternIter* theOrgPatternIter;
//	TaggedObjectStorage* theOrgPatterns;
	int Numstep;
	double delta;
	int NumLoadPatterns;
	int* LoadPatterns;


private:
};
#endif
