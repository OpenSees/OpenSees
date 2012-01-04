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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/CrossingRateAnalyzer.h,v $


#ifndef CrossingRateAnalyzer_h
#define CrossingRateAnalyzer_h

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ReliabilityDomain.h>
#include <FindDesignPointAlgorithm.h>
#include <Domain.h>
#include <FunctionEvaluator.h>
#include <GradientEvaluator.h>
#include <NormalRV.h>
#include <RandomProcess.h>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;

class CrossingRateAnalyzer
{
  public:
	CrossingRateAnalyzer(ReliabilityDomain* passedReliabilityDomain,
						 FindDesignPointAlgorithm* passedFindDesignPointAlgorithm,
						 FunctionEvaluator* passedGFunEvaluator,
						 GradientEvaluator* passedGradGEvaluator,
						 int passedanalysisType=2,
						 double passedlittleDt=0.1,
  						 bool passedprint=false);
	~CrossingRateAnalyzer();

	void setAnalysis(int passedNstep, Vector& xdes, 
					 Vector& alphades, double betades);
	void setRandomProcess(RandomProcess* passedRandomProcess);
	void setdelta(double passeddelta){ delta=passeddelta;}
//	void setdelta_pulse(double passeddelta){ delta_pulse=passeddelta;}
	double computeRate(void);
	void clear();

  protected: 
    
  private:
	ReliabilityDomain* theReliabilityDomain;
	FindDesignPointAlgorithm* theFindDesignPointAlgorithm;
	FunctionEvaluator*   theGFunEvaluator;
	GradientEvaluator* theGradGEvaluator;
	RandomProcess* theRandomProcess;
	double computeRate1(void);
	double computeRate2(void);
	double functionToIntegrate(double rho, double beta1, double beta2);
	int analysisType;
	double littleDt;
	double delta;
	double delta_pulse;
	bool print;

	int numSteps;
	int numRV;
	Vector* uDesign0;
	Vector* uShifted;
	double beta0;
	double betaShifted;
	double beta1;

	ofstream output;

};

#endif
