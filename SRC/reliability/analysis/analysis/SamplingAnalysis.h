// SamplingAnalysis.h: interface for the SamplingAnalysis class.
//
//////////////////////////////////////////////////////////////////////

#if !defined   SAMPLINGANALYSIS_H_100 
#define        SAMPLINGANALYSIS_H_100


#include "ReliabilityAnalysis.h"

class SamplingAnalysis : public ReliabilityAnalysis  
{
public:
	virtual bool getContribution();
	virtual double getProbability();
	virtual double getCov();
	virtual int getNumOfSimulations();

	virtual int analyze(void) =0;

	virtual double getSampledValue(double)=0;
	SamplingAnalysis();
	virtual ~SamplingAnalysis();

protected:
	double maxNumOfIterations;
	int numOfSimulations;
	double cov;
	double probability;

private:



};

#endif //  
