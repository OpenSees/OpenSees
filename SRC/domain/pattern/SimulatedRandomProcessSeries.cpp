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

SimulatedRandomProcessSeries::SimulatedRandomProcessSeries(int tag,
							   RandomNumberGenerator *theRandNumGenerator,
							   Spectrum *theSpectr,
							   int numFreqInt,
							   double pmean)
 :TimeSeries(tag, TSERIES_TAG_SimulatedRandomProcessSeries)
{
	theRandomNumberGenerator = theRandNumGenerator;
	theSpectrum = theSpectr;
	numFreqIntervals = numFreqInt;
	mean = pmean;

	
	// Generate random numbers, uniformly distributed between 0 and 2pi
	double pi = 3.14159265358979;
	theRandomNumberGenerator->generate_nIndependentUniformNumbers(numFreqIntervals,0.0,(2*pi));
	Vector theta1 = theRandomNumberGenerator->getGeneratedNumbers();
	theta = new Vector(theta1);
	

	// Generate standard normal random numbers
	theRandomNumberGenerator->generate_nIndependentStdNormalNumbers(numFreqIntervals);
	Vector A1 = theRandomNumberGenerator->getGeneratedNumbers();
	A = new Vector(A1);


	// Length of each interval
	deltaW = (theSpectrum->getMaxFrequency()-theSpectrum->getMinFrequency())/numFreqIntervals;


}



TimeSeries *
SimulatedRandomProcessSeries::getCopy(void) {
  opserr << "SimulatedRandomProcessSeries::getCopy(void) - not yet implemented\n";
  return 0;
}

SimulatedRandomProcessSeries::~SimulatedRandomProcessSeries()
{
	if (theta != 0)
		delete theta;
	if (A != 0) 
		delete A;
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
			factor += sqrt(2.0*S*deltaW) * (*A)(i) * cos(W*time+(*theta)(i));
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
