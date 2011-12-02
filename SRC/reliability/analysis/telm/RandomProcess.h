// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/RandomProcess.h,v $

#ifndef RandomProcess_h
#define RandomProcess_h

#include<LoadPatternIter.h>
#include<LoadPattern.h>
#include<UniformExcitation.h>
#include<GroundMotion.h>
#include<TimeSeries.h>
#include<Domain.h>
#include<ReliabilityDomain.h>
#include<NewDiscretizedRandomProcessSeries.h>
#include<fstream>
#include<iomanip>
#include<iostream>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;

class RandomProcess
{
  public:
   	RandomProcess(Domain *theDomain,
				  ReliabilityDomain *theReliabilityDomain,
				  bool print);
	~RandomProcess();

	int getNumTotalPulse(void){ return numTotalPulse; }
	int getNumOfActivePulses(double time);
	double getPulseTime(int ipulse);
	double getFactorSensitivity(double ctime, double ktime);
	int getRVID(int ipulse){ return iRVPulse[ipulse]; }
	double getDeltaPulse(void){ return delta_Pulse; }
	Vector getExcitation(int Nstep, double delta);
	int getRVseqid(int ipulse);
  protected: 
    
  private:
	void AnalyzeProcess(void);
	void FindRandomProcess(Domain*);
	int getDirection(void){ return DirExcitation; }

	NewDiscretizedRandomProcessSeries* getProcessptr(void)
	{ return theDiscreSeriese; }

	ReliabilityDomain* theReliabilityDomain;
	NewDiscretizedRandomProcessSeries* theDiscreSeriese;
	  int numTotalPulse;
	  int DirExcitation;
	  int* iRVPulse; 
	  Vector* timePulse;
	  double delta_Pulse;
	  bool print;
	  ofstream output;

};

#endif
