// SamplingAnalysis.cpp: implementation of the SamplingAnalysis class.
//
//////////////////////////////////////////////////////////////////////

#include "SamplingAnalysis.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SamplingAnalysis::SamplingAnalysis()
{
	cov = 0;
	numOfSimulations = 0;
	probability = 0;
}

SamplingAnalysis::~SamplingAnalysis()
{

}



int SamplingAnalysis::getNumOfSimulations()
{
	return numOfSimulations;
}

double SamplingAnalysis::getCov()
{
	return cov;
}

double SamplingAnalysis::getProbability()
{
	return probability;
}

bool SamplingAnalysis::getContribution()
{
	return true; 
}
