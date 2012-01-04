
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/InitialPointBuilder.h,v $


#ifndef InitialPointBuilder_h
#define InitialPointBuilder_h

#include <fstream>
#include <iomanip>
#include <iostream>
#include <RandomProcess.h>
#include <OutCrossingResults.h>
#include <ReliabilityDomain.h>
#include <FunctionEvaluator.h>
#include <Matrix.h>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;


class InitialPointBuilder
{
 public:
  InitialPointBuilder(ReliabilityDomain *theReliabilityDomain,
		      FunctionEvaluator* theGFunEvaluator);
 virtual  ~InitialPointBuilder();
	
 void setRandomProcess(RandomProcess* passedRandomProcess);
 void setDt(double passedDt){ dt=passedDt; }
 void setthreshold(double value){ threshold = value; }
 bool startMirror(){ return start_mirror; }
 void setLsf(int passedLsf){ lsf=passedLsf; }
 virtual Vector buildInitialPoint(int)=0;
 virtual void setPrevResults(int, double, double, const Vector&) =0;
 virtual int getnumAna(void)=0;
 virtual void clear()=0;
 virtual	void setMirrorImageExcitation(int, Vector)=0;
 virtual int getType(void)=0;
 virtual void setOutCrossingResults
   (ofstream& outputFile, OutCrossingResults* theOutCrossingResults) =0;
 
 protected:
 
 void clearInitialPointBuilder();
 RandomProcess* theRandomProcess;
 FunctionEvaluator* theGFunEvaluator;
 ReliabilityDomain* theReliabilityDomain;
 int NumTotalPulse;
 double dt;
 double dpulse;
 int numRV;
 Vector* xmean;
 double threshold;
 RandomVariable* theRV;
 Vector* xinitial;
 bool start_mirror;
 int lsf;
 
    
 private:

};

#endif
