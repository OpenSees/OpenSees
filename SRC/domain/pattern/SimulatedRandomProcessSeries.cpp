/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.3 $
// $Date: 2010-02-04 18:32:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/SimulatedRandomProcessSeries.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu), February 2002
// Revised: 
//

#include <SimulatedRandomProcessSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <Spectrum.h>
#include <RandomNumberGenerator.h>
#include <classTags.h>
//#include <fstream>
#include <math.h>

#include <elementAPI.h>
#include <Spectrum.h>
#include <OpenSeesReliabilityCommands.h>

CStdLibRandGenerator SimulatedRandomProcessSeries::randGenerator;

void *
OPS_SimulatedRandomProcessSeries(void)
{
  // Pointer that will be returned
  TimeSeries *theSeries = 0;

  int tag;
  double mean;
  int spectrumTag, NfreqIntervals;

  int numData = 1;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs < 4) {
    opserr << "ERROR: Insufficient arguments for SimulatedRandomProcess series" << endln;
    return 0;
  }

  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid tag for SimulatedRandomProcess" << endln;
    return 0;
  }
  if (OPS_GetInt(&numData, &spectrumTag) != 0) {
    opserr << "WARNING invalid spectrum tag for SimulatedRandomProcess with tag: " << tag << endln;
    return 0;
  }  
  if (OPS_GetDouble(&numData, &mean) != 0) {
    opserr << "WARNING invalid mean for SimulatedRandomProcess with tag: " << tag << endln;
    return 0;
  }
  if (OPS_GetInt(&numData, &NfreqIntervals) != 0) {
    opserr << "WARNING invalid num frequency internvals for SimulatedRandomProcess with tag: " << tag << endln;
    return 0;
  }  

  ReliabilityDomain *theReliabilityDomain = OPS_GetReliabilityDomain();
  if (theReliabilityDomain == 0) {
    opserr << "ERROR SimulatedRandomProcess -- reliability domain not defined" << endln;
    return 0;
  }

  Spectrum *theSpectrum = theReliabilityDomain->getSpectrum(spectrumTag);
  if (theSpectrum == 0) {
    opserr << "ERROR SimulatedRandomProcess -- spectrum with tag " << spectrumTag
	   << " not found" << endln;
    return 0;
  }

  /*
  RandomNumberGenerator *theRandNumGenerator = cmds->getRandomNumberGenerator();
  if (theRandNumGenerator == 0) {
    opserr << "ERROR SimulatedRandomProcess -- random number generator not defined" << endln;
    return 0;
  }
  */
  RandomNumberGenerator *theRandNumGenerator = 0;
  
  theSeries = new SimulatedRandomProcessSeries(tag, theRandNumGenerator,
					       theSpectrum, NfreqIntervals, mean);
  
  if (theSeries == 0) {
    opserr << "WARNING ran out of memory creating DiscretizedRandomProcess series with tag: " << tag << endln;
    return 0;
  }

  return theSeries;
}

SimulatedRandomProcessSeries::SimulatedRandomProcessSeries(int tag,
							   RandomNumberGenerator *theRandNumGenerator,
							   Spectrum *theSpectr,
							   int numFreqInt,
							   double pmean)
  :TimeSeries(tag, TSERIES_TAG_SimulatedRandomProcessSeries),
   theSpectrum(theSpectr), numFreqIntervals(numFreqInt), mean(pmean),
   deltaW(0.0), theta(numFreqInt), A(numFreqInt)
{
  //theRandomNumberGenerator = theRandNumGenerator;
	theRandomNumberGenerator = &randGenerator;
	
	// Generate random numbers, uniformly distributed between 0 and 2pi
	double pi = 3.14159265358979;
	theRandomNumberGenerator->generate_nIndependentUniformNumbers(numFreqIntervals,0.0,(2*pi));
	theta = theRandomNumberGenerator->getGeneratedNumbers();

	// Generate standard normal random numbers
	theRandomNumberGenerator->generate_nIndependentStdNormalNumbers(numFreqIntervals);
	A = theRandomNumberGenerator->getGeneratedNumbers();

	// Length of each interval
	deltaW = (theSpectrum->getMaxFrequency()-theSpectrum->getMinFrequency())/numFreqIntervals;
}



TimeSeries *
SimulatedRandomProcessSeries::getCopy(void) {

  SimulatedRandomProcessSeries *theCopy =
    new SimulatedRandomProcessSeries(this->getTag(),
				     theRandomNumberGenerator,
				     theSpectrum, numFreqIntervals, mean);
  theCopy->theta = theta;
  theCopy->A = A;

  return theCopy;
}

SimulatedRandomProcessSeries::~SimulatedRandomProcessSeries()
{

}


double
SimulatedRandomProcessSeries::getFactor(double time)
{
//static ofstream outputFile( "simulated_process.out" , ios::out );

	if (time == 0.0) {
		return 0.0;
	}
	else {

		// Add up over all frequency intervals
		double factor = 0.0;
		double W, S;
		for (int i=0; i<numFreqIntervals; i++) {
			W = (i+0.5)*deltaW+theSpectrum->getMinFrequency();
			S = theSpectrum->getAmplitude(W);
			factor += sqrt(2.0*S*deltaW) * A(i) * cos(W*time+theta(i));
		}

//outputFile << (mean+factor) << endl;

		return (mean + factor);
	}
}




int
SimulatedRandomProcessSeries::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}


int 
SimulatedRandomProcessSeries::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	return 0;    
}


void
SimulatedRandomProcessSeries::Print(OPS_Stream &s, int flag)
{
}
