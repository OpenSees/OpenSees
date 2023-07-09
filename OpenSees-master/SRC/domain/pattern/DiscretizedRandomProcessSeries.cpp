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
                                                                        
// $Revision: 1.7 $
// $Date: 2010-02-04 18:32:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/DiscretizedRandomProcessSeries.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu), February 2002
// Revised: 
//

#include <DiscretizedRandomProcessSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <ModulatingFunction.h>
#include <Filter.h>
#include <classTags.h>
#include <Parameter.h>
#include <math.h>
#include <stdlib.h>

DiscretizedRandomProcessSeries::DiscretizedRandomProcessSeries(int tag,
							       int num, 
							       ModulatingFunction **theModFuncs,
							       double p_mean,
							       double p_maxStdv)
 :TimeSeries(tag, TSERIES_TAG_DiscretizedRandomProcessSeries)
{
  randomVariables = 0;
  kickInTimes = 0;
  theModulatingFunctions = theModFuncs;
  numModFuncs = num;
  mean = p_mean;
  maxStdv = p_maxStdv;
  
  c = 0.0;
}

TimeSeries *
DiscretizedRandomProcessSeries::getCopy(void) 
{
  opserr << "DiscretizedRandomProcessSeries::getCopy() - not yet implemented\n";
  return 0;
}


DiscretizedRandomProcessSeries::~DiscretizedRandomProcessSeries()
{
	if (randomVariables != 0) 
		delete randomVariables;

	if (kickInTimes != 0) 
		delete kickInTimes;
}


double
DiscretizedRandomProcessSeries::getFactor(double time)
{
	if (time == 0.0) {
		return 0.0;
	}
	else if (randomVariables == 0 || kickInTimes == 0) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactor(): " << endln
			<< " random variables or kick-in times vector(s) do not exist. " << endln;
		return 0.0;
	}
	else if (kickInTimes->Size() != randomVariables->Size() ) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactor(): " << endln
			<< " number of random variables is not the same as kick-in times. " << endln;
		return 0.0;
	}
	else {
		double sum1;
		double sum2;
		int nrv = 0;
		double modFuncAmplitude, filterAmplitude;
		Filter *theFilter;

		// Loop over all modulating functions
		sum1 = 0.0;
		for (int k=0; k<numModFuncs; k++) {

			// Get value of modulating function number k at time t
			modFuncAmplitude = theModulatingFunctions[k]->getAmplitude(time);
			theFilter = theModulatingFunctions[k]->getFilter();

			// Number of discretizing random variables
			nrv = randomVariables->Size();

			// Loop over all active rv's 
			sum2 = 0.0;
			for (int i=0; i<nrv; i++) {

				// Get value of filter for argument (t-ti)
				filterAmplitude = theFilter->getAmplitude(time-(*kickInTimes)(i));
				
				// Add contribution 'ui * hi'
				sum2 += (*randomVariables)(i) * filterAmplitude;

				// Break when we get to inactive rv's
				if (time-(*kickInTimes)(i) < 0.0) {
					break;
				}
			}

			sum1 += sum2*modFuncAmplitude;
		}

		double result = mean + c*sum1;
		return result;
	}
}


double
DiscretizedRandomProcessSeries::getFactorSensitivity(double time)
{
	// The parameterID has been set to the number of 
	// the random variable in question

	// So, do the same thing as above, just set x(i-1) equal to 1.0
	// for i==parameterID

	if (time == 0.0) {
		return 0.0;
	}
	else if (randomVariables == 0 || kickInTimes == 0) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " random variables or kick-in times vector(s) do not exist. " << endln;
		return 0.0;
	}
	else if (kickInTimes->Size() != randomVariables->Size() ) {
		opserr << "ERROR in DiscretizedRandomProcessSeries::getFactorSensitivity(): " << endln
			<< " number of random variables is not the same as kick-in times. " << endln;
		return 0.0;
	}
	else {

		double sum1;
		double sum2;
		int nrv = 0;
		double modFuncAmplitude;
		Filter *theFilter;

		// Loop over all modulating functions
		sum1 = 0.0;
		for (int k=0; k<numModFuncs; k++) {

			// Get value of modulating function number k at time t
			modFuncAmplitude = theModulatingFunctions[k]->getAmplitude(time);
			theFilter = theModulatingFunctions[k]->getFilter();

			// Number of discretizing random variables
			nrv = randomVariables->Size();

			// Loop over all rv's (even though some may be zero at this time)
			sum2 = theFilter->getAmplitude(time-(*kickInTimes)(parameterID-1));
			sum1 += sum2*modFuncAmplitude;
		}

		double result = mean + c*sum1;
		return result;
	}
}


int
DiscretizedRandomProcessSeries::sendSelf(int commitTag, Channel &theChannel)
{
	return 0;
}


int 
DiscretizedRandomProcessSeries::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	return 0;    
}


void
DiscretizedRandomProcessSeries::Print(OPS_Stream &s, int flag)
{
}

int
DiscretizedRandomProcessSeries::setParameter(const char **argv, int argc,
					     Parameter &param)
{
  if (argc < 1)
    return -1;

  // **** MHS needs to fix this!!
  //int rvNumber = info.theInt;
  int rvNumber = 1;  // to get it to compile for now
  // **********************

	// The second argument tells when the random variable "kicks in".
	// Store this in a table...
	// In case the vector doesn't exist
	if (kickInTimes == 0) {
		kickInTimes = new Vector(rvNumber);
		(*kickInTimes)(rvNumber-1) = (double)atof(argv[0]);

		// Assume more than one random variable, so don't 
		// update factor 'c' here.
	}
	// In case the vector isn't big enough
	else if (kickInTimes->Size() < rvNumber) {

		// Store old values in a temporary vector
		Vector temp(*kickInTimes);

		// Create a large enough vector
		delete kickInTimes;
		kickInTimes = new Vector(rvNumber);

		// Put in old values
		for (int i=0; i<temp.Size(); i++) {
			(*kickInTimes)(i) = temp(i);
		}

		// Put in new value
		(*kickInTimes)(rvNumber-1) = (double)atof(argv[0]);


		/////// Update factor 'c' /////////

		// Number of discretizing random variables
		int nrv = kickInTimes->Size();

		double new_c;
		double denominator;

		// Loop over all time instances
		for (int t=1; t<=nrv; t++) {

			denominator = 0.0;

			// Loop over modulating functions
			for (int k=0; k<numModFuncs; k++) {

				// Get value of modulating function number k at time t
				Filter *theFilter_k = theModulatingFunctions[k]->getFilter();
				double modFuncAmplitude_k = theModulatingFunctions[k]->getAmplitude((*kickInTimes)(t-1));
		
				
				// Loop over modulating functions
				for (int l=0; l<numModFuncs; l++) {


					// Get value of modulating function number l at time t
					Filter *theFilter_l = theModulatingFunctions[l]->getFilter();
					double modFuncAmplitude_l = theModulatingFunctions[l]->getAmplitude((*kickInTimes)(t-1));


					// Loop over all rv's (even though some may be zero at this time)
					for (int i=0; i<nrv; i++) {

						// Get value of filters for argument (t-ti)
						double filterAmplitude_k = theFilter_k->getAmplitude(((*kickInTimes)(t-1))-(*kickInTimes)(i));
						double filterAmplitude_l = theFilter_l->getAmplitude(((*kickInTimes)(t-1))-(*kickInTimes)(i));
						
						// Add contribution 'ui * hi'
						denominator += filterAmplitude_k*filterAmplitude_l
									 * modFuncAmplitude_k*modFuncAmplitude_l;
					}
				}
			}
		}

		new_c = sqrt(maxStdv*maxStdv/denominator);

		if (c==0.0) {
			c = new_c;
		}
		else if (new_c < c) {
			c = new_c;
		}

c = maxStdv;
opserr << "c: " << c << endln;

		//////////////////////////////////////
	}
	else {
		(*kickInTimes)(rvNumber-1) = (double)atof(argv[0]);
	}

	// The random variable number is returned as a parameter ID
	return param.addObject(rvNumber, this);
}

int
DiscretizedRandomProcessSeries::updateParameter(int parameterID, Information &info)
{
	// In case the vector doesn't exist
	if (randomVariables == 0) {
		randomVariables = new Vector(parameterID);
		(*randomVariables)(parameterID-1) = info.theDouble;
	}
	// In case the vector isn't big enough
	else if (randomVariables->Size() < parameterID) {

		// Store old values in a temporary vector
		Vector temp(*randomVariables);

		// Create a large enough vector
		delete randomVariables;
		randomVariables = new Vector(parameterID);

		// Put in old values
		for (int i=0; i<temp.Size(); i++) {
			(*randomVariables)(i) = temp(i);
		}

		// Put in new value
		(*randomVariables)(parameterID-1) = info.theDouble;
	}
	else {
		(*randomVariables)(parameterID-1) = info.theDouble;
	}

	return 0;
}

int
DiscretizedRandomProcessSeries::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}
