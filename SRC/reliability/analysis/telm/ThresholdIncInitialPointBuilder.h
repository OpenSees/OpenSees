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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/ThresholdIncInitialPointBuilder.h,v $

#ifndef ThresholdIncInitialPointBuilder_h
#define ThresholdIncInitialPointBuilder_h

#include <fstream>
#include <iomanip>
#include <iostream>
#include <InitialPointBuilder.h>
#include <FindDesignPointAlgorithm.h>
#include <NewSearchWithStepSizeAndStepDirection.h>
#include <OutCrossingResults.h>
#include <NormalRV.h>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;

class ThresholdIncInitialPointBuilder : public InitialPointBuilder
{
  public:
	ThresholdIncInitialPointBuilder(ReliabilityDomain *theReliabilityDomain,
						FunctionEvaluator* theGFunEvaluator,
						FindDesignPointAlgorithm* theFindDesignPointAlgorithm,
						//NewSearchWithStepSizeAndStepDirection* theFindDesignPointAlgorithm,
						int maxDivide=20,
						double eps=0.05,
						bool start_mirror=true,
						int MaxLineSearch=100,
						double thresholdRatio=0.01, 
						bool print=false);
	~ThresholdIncInitialPointBuilder();

	Vector buildInitialPoint(int);
	int getnumAna(void){ return numAna; };
	int getType(void){ return 1; }
	void clear();
	void setPrevResults(int, double, double, const Vector&);
	void setOutCrossingResults
		(ofstream& outputFile, OutCrossingResults* theOutCrossingResults);
	void setMirrorImageExcitation(int, Vector);

  protected: 
    
  private:
	//NewSearchWithStepSizeAndStepDirection* theFindDesignPointAlgorithm;
	FindDesignPointAlgorithm* theFindDesignPointAlgorithm;
	int numAna;
	bool print;
	ofstream output;
	ofstream outputFile;
	Vector* prevDesign;
	double eps;
	double completedSteps;
	double completedTime;
	double completedThreshold;
	Vector* xtemp;
	Vector* xcompleted;
	int maxDivide;
	OutCrossingResults* theOutCrossingResults;
};

#endif
