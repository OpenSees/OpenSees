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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NonStatFirstPassageAnalyzer.h,v $

#ifndef NonStatFirstPassageAnalyzer_h
#define NonStatFirstPassageAnalyzer_h

#include <FirstPassageAnalyzer.h>
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

class NonStatFirstPassageAnalyzer : public FirstPassageAnalyzer 
{
 public:
  NonStatFirstPassageAnalyzer(ReliabilityDomain* passedReliabilityDomain,
    			      FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
			      FunctionEvaluator* passedGFunEvaluator,
			      FOSeriesSimulation* passedFOSeriesSimulation,
			      int passedanalysisType,
			      int interpolationType,
			      bool passedtwoside,
			      bool passedprint);
  
  ~NonStatFirstPassageAnalyzer();
  Vector componentFisrtPassage
    (Vector* udes, Vector* xdes, Vector* alpha, Vector* hfunc, 
     double beta, double pf, double analysistime, int analysisstep,
     TimePoints* passedTimePoints,
     int passedlsf, int passedidfragility,
     ofstream& outputFile); 
  
  void systemFisrtPassage(void);
  bool stationary(void){return false;}
  
  
 protected: 
  
 private:
  bool print;
  int interpolationType;
  void componentFisrtPassage1();
  void componentFisrtPassage2();
  Vector* FPprob;
  Vector* covres;
  Vector* betares;
  int* numSim;
  ofstream output;
};

#endif
