/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
** ****************************************************************** */

// $Revision: 1.2 $
// $Date: 2008-05-13 16:30:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/RandomVibrationAnalysis.h,v $

#ifndef RandomVibrationAnalysis_h
#define RandomVibrationAnalysis_h

#include <Domain.h>
#include <ReliabilityDomain.h>
#include <FindDesignPointAlgorithm.h>
#include <NewSearchWithStepSizeAndStepDirection.h>
#include <ReliabilityConvergenceCheck.h>
#include <InitialPointBuilder.h>
#include <CrossingRateAnalyzer.h>
#include <FirstPassageAnalyzer.h>
 
#include <RandomProcess.h>
#include <TimePoints.h>
#include <ReliabilityConvergenceCheck.h>
#include <OutCrossingResults.h>
#include <MatrixOperations.h>

#include <NormalRV.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;

class RandomVibrationAnalysis
{
  public:
	RandomVibrationAnalysis(ReliabilityDomain* passedReliabilityDomain,
							FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
							//NewSearchWithStepSizeAndStepDirection* passedFindDesignPointAlgorithm,
							Domain* passedDomain, 
							InitialPointBuilder* passedInitialPointBuilder,
							CrossingRateAnalyzer* passedCrossingRateAnalyzer,
							FirstPassageAnalyzer* passedFirstPassageAnalyzer,
							FunctionEvaluator* passedGFunEvaluator,
							GradientEvaluator* passedGradGEvaluator,
							ReliabilityConvergenceCheck* passedReliabilityConvergenceCheck,
							double passedStartTime,
							double passedEndTime,
							double passedTimeInterval,
							double passedStartAnalysis,
							double passedFragMin,
							double passedFragInt,
							int passednFrag,
							int passeddesignpoint,
							bool stationary,
							bool passedmirrorimage,
							bool passedinitialpoint,
							bool passedfirstpassage,
						    TCL_Char *passedFileName,
							char* passedfileBinary,
							Tcl_Interp* passedTclInterp,
	 						bool passedprint);
	~RandomVibrationAnalysis();

	void analyze(void);

  protected: 
    
  private:

    // Functions //
//	void nonStatDesignPoints(void);
//	void StatDesignPoints(void);
	void DesignPoints(void);
	void RdDesignPoints(void);
	double crossingRate(void);
	void saveDesResults(int loc, int presult, int plsf, double pthre);
	void analyzeComponentFirstPasssage(int,int,int);
	void setScale(double scale);
	double setLimitState(LimitStateFunction* theLSF);
	void printPerformanceFunction(int lsf);
	void PrintTitle(void);
	void findTELS(void);

	// classes //
	Domain* theDomain; 
	ReliabilityDomain* theReliabilityDomain;
	FindDesignPointAlgorithm* theFindDesignPointAlgorithm;
	//NewSearchWithStepSizeAndStepDirection* theFindDesignPointAlgorithm;
    FunctionEvaluator* theGFunEvaluator;
    GradientEvaluator* theGradGEvaluator;
	ReliabilityConvergenceCheck* theReliabilityConvergenceCheck;

	InitialPointBuilder* theInitialPointBuilder;
	CrossingRateAnalyzer* theCrossingRateAnalyzer;
	FirstPassageAnalyzer* theFirstPassageAnalyzer;

	RandomProcess* theRandomProcess;
	TimePoints* theTimePoints;
	LimitStateFunction *theLSF;
	OutCrossingResults* theOutCrossingResults;
    TaggedObjectStorage* thePerformFuncCoeffs;
	PerformanceFunctionCoefficientIter* thePfCoeffIter;

	// variables //
	double delta;
	double delta_Pulse;
	int designpoint;
	bool stationary;
	bool mirrorimage;
	bool initialpoint;
	bool firstpassage;
	bool print;
	char* fileName;
	char* fileBinary;

	int numRV;
	int numRVPos;
	int numLsf;
	int NumTotalPulse;
	int NumTimePoints;
	int StartSteps;

	int numFragility;
	Vector* Fragility;

 	double Scorg;
	bool Scfixorg;

	Vector* xDesTmp;
	Vector* uDesTmp;
	Vector* alphaTmp;
	Vector* hfunc;
	Vector* ExcTmp;
	double betaTmp;
	double pfTmp;
	double RateTmp;
//	double Pmn;
	ofstream output;
	ofstream outputTELs;
	ofstream outputFile;
	Tcl_Interp *theTclInterp;
	double Performthreshold;

	int mirStep;
	Vector* mirExcitation;
	double mirTime;

	int baseStep;
	int baseStep0;
	Vector* baseExcitation;
	Vector* baseExcitation0;
	double crossRate;

	int AnalysisStep;
	double AnalysisTime;

	Vector* xinitial;

	int numAnaLineSearch;
	int numAnaIncSens;
	int numEvalinFORM;
	int numIterinFORM;

	int iresult;
	int numDesPoints;
	Matrix** FirstPassageMatrix;

};

#endif
