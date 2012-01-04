
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/InitialPointBuilder.cpp,v $


#include <InitialPointBuilder.h>
InitialPointBuilder::InitialPointBuilder
					(ReliabilityDomain* passedReliabilityDomain,
	                 FunctionEvaluator* passedGFunEvaluator)
{
	theReliabilityDomain=passedReliabilityDomain;
	numRV = theReliabilityDomain->getNumberOfRandomVariables();
	xmean=0;
	xmean = new Vector(numRV);
	if (xmean ==0 ){
		opserr << " Insufficient memory /n";
		opserr << " InitialPointBuilder xmean/n";
		exit(-1);
	}
	for(int i=0; i<numRV; i++){
		theRV=theReliabilityDomain->getRandomVariablePtr(i+1);
		(*xmean)(i)=theRV->getMean();
	}
	xinitial=0;
	theGFunEvaluator=passedGFunEvaluator;
} 
InitialPointBuilder::~InitialPointBuilder()
{
	if( xinitial !=0 ){ delete xinitial; xinitial=0;}
	if( xmean !=0 ){ delete xmean; xmean=0;}
} 
void InitialPointBuilder::clearInitialPointBuilder()
{
	if( xinitial !=0 ){ delete xinitial; xinitial=0;}
	if( xmean !=0 ){ delete xmean; xmean=0;}
} 
void InitialPointBuilder::setRandomProcess(RandomProcess* passedRandomProcess)
{
	theRandomProcess=passedRandomProcess;
	NumTotalPulse=theRandomProcess->getNumTotalPulse();
	dpulse=theRandomProcess->getDeltaPulse();
}
