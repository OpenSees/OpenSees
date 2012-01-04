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
                                                                        
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/FirstPassageAnalyzer.h,v $

#ifndef FirstPassageAnalyzer_h
#define FirstPassageAnalyzer_h

#include <ReliabilityDomain.h>
#include <FindDesignPointAlgorithm.h>
#include <FunctionEvaluator.h>
#include <RandomProcess.h>
#include <FOSeriesSimulation.h>
#include <OutCrossingResults.h>
#include <TimePoints.h>

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

class FirstPassageAnalyzer
{
  public:
  FirstPassageAnalyzer(ReliabilityDomain* passedReliabilityDomain,
		       FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
		       FunctionEvaluator* passedGFunEvaluator,
		       FOSeriesSimulation* passedFOSeriesSimulation,
		       int passedAnalysisType,
		       bool passedtwoside);
  virtual ~FirstPassageAnalyzer();
  
  void setRandomProcess(RandomProcess* passedRandomProcess);
  virtual Vector componentFisrtPassage
    (Vector* udes, Vector* xdes, Vector* alpha, Vector* hfunc,
     double beta,double pfTmp, double AnalysisTime, int AnalysisStep,
     TimePoints* passedTimePoints,
     int passedlsf, int passedidfragility,
     ofstream& outputFile)=0; 
  
  virtual void systemFisrtPassage(void)=0;
  virtual bool stationary(void)=0;
  
 protected: 
  ReliabilityDomain* theReliabilityDomain;
  RandomProcess* theRandomProcess;
  FindDesignPointAlgorithm* theFindDesignPointAlgorithm;
  FunctionEvaluator* theGFunEvaluator;
  FOSeriesSimulation* theFOSeriesSimulation;
  //	OutCrossingResults* theOutCrossingResult;
  //	OutCrossingResults** thecomponentResults;
  TimePoints* theTimePoints;
  
  int analysisType;
  int numRV;
  int numRVPos;
  int numTotalPulses;
  int numOtherRVs;
  int numTimePoints;
  int numLsf;
  double delta_Pulse;
  double delta;
  bool twoside;
  Vector* timepoints;
  Vector* udesres;
  Vector* xdesres;
  Vector* alphadesres;
  Vector* hfuncdesres;
  double betadesres;
  double pfdesres;
  double analysistime;
  int analysisstep;
  void setComponentResult(Vector* udes, Vector* xdes, Vector* alpha, 
			  Vector* hfunc, double beta, double pf,
			  double analysistime, int analysisstep,
			  TimePoints* passedTimePoints,
			  int passedlsf, int passedidfragility,
			  bool stationary);
  //	void setSystemResults(OutCrossingResults* passedOutCrossingResults,
  //						  TimePoints* passedTimePoints,
  //						  int passedlsf, int passedidfragility,
  //						  bool stationary);
  
 private:
  
};

#endif
