// MonteCarloResponseAnalysis.h: interface for the MonteCarloResponseAnalysis class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MONTECARLORESPONSEANALYSIS
#define MONTECARLORESPONSEANALYSIS

#include <ReliabilityDomain.h>
#include <ProbabilityTransformation.h>
#include <RandomNumberGenerator.h>


class MonteCarloResponseAnalysis  
{
public:
	MonteCarloResponseAnalysis(ReliabilityDomain *passedReliabilityDomain,
						Tcl_Interp *passedTclInterp,
						ProbabilityTransformation *passedProbabilityTransformation,
						RandomNumberGenerator *passedRandomNumberGenerator,
						int passedNumberOfSimulations,
						int printFlag,
						TCL_Char *outputFileName,
						TCL_Char *tclFileToRunFileName,
						int seed
						);


	virtual ~MonteCarloResponseAnalysis();
	int analyze();

private:
	ReliabilityDomain *theReliabilityDomain;
	Tcl_Interp *theTclInterp;
	ProbabilityTransformation *theProbabilityTransformation;
	RandomNumberGenerator *theRandomNumberGenerator;
	int numberOfSimulations;
	int printFlag;
	char fileName[25];
	char * tclFileToRun;
	int seed;



};

#endif 
