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
**   Quan Gu (qgu@ucsd.edu)                                           **
**   Joel P. Conte (jpconte@ucsd.edu)                                 **
** ****************************************************************** */
                                                                        
 
//
// Written by  Quan Gu UCSD
//

#if !defined ORTHOGONALPLANESAMPLINGANALYSIS_H__
#define      ORTHOGONALPLANESAMPLINGANALYSIS_H__

#include <SamplingAnalysis.h>
#include <ZeroFindingAlgorithm.h>

#include <ReliabilityAnalysis.h>
#include <ReliabilityDomain.h>
#include <ProbabilityTransformation.h>
#include <RandomNumberGenerator.h>
#include <FunctionEvaluator.h>

#include <fstream>
using std::ofstream;

class OrthogonalPlaneSamplingAnalysis : public SamplingAnalysis  
{
public:
	bool getContribution();
	double littleDt;
	double getG2FromG1(double,double);
//	int getNumOfSimulations();  // already in father class
	int getSeed();
	void setSeed(int);
	int getNumOfGFunEvaluations();
	double getSampledValue(double x);
	int analyze(void);
//	OrthogonalPlaneSamplingAnalysis();
	OrthogonalPlaneSamplingAnalysis(   Tcl_Interp *interp,
		                ReliabilityDomain *passedReliabilityDomain,
						ProbabilityTransformation *passedProbabilityTransformation,
						FunctionEvaluator *passedGFunEvaluator,
						RandomNumberGenerator *passedRandomNumberGenerator,
						int passedNumberOfSimulations,
						int passedMaxNumOfIterations,
						double passedTargetCOV,
						int pPrintFlag,
						TCL_Char *passedFileName,
						Vector * pDesignPoint,
						int pAnalysisTypeTag, 
						int zeroFindingType,
						double pFuncTol,
						double pVarTol,
						int pMaxIter,
						double pLittleDt);
 
	virtual ~OrthogonalPlaneSamplingAnalysis();


private:
	int analysisType;
	bool contribution;

	double beta;
	Vector * theDesignPointInUSpace;
	Vector *theDesignPoint; // NOTE: this design point is in Physical space !!
	Vector * uPrime;
	int seed;
	int numOfGFunEvaluations;
	ZeroFindingAlgorithm * theAlgorithm;

 	ReliabilityDomain *theReliabilityDomain;
	ProbabilityTransformation *theProbabilityTransformation;
	FunctionEvaluator *theGFunEvaluator;
	RandomNumberGenerator *theRandomNumberGenerator;
	int maxNumOfIterations; // in zerofinding algorithm
	double targetCOV;
	double samplingStdv;
	int printFlag;
	char *fileName;
	int analysisTypeTag;  // 1: timeInvariant. 2: upCrossing..... 
	Tcl_Interp *theInterp;
	double funcTol;
	double varTol;
	int maxIter;



};

#endif // 
