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
                                                                        
// $Revision: 1.1 $
// $Date: 2000-12-12 07:40:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TrapezoidalTimeSeriesIntegrator.cpp,v $
                                                                        
                                                                        
// File: ~/domain/pattern/TrapezoidalTimeSeriesIntegrator.cpp
// 
// Written: MHS
// Created: 10/99
// Revision: A
//
// Description: This file contains the class definition for 
// a TrapezoidalTimeSeriesIntegrator, which integrates a
// ground motion TimeSeries using the trapezoidal rule.
//
// What: "@(#) TrapezoidalTimeSeriesIntegrator.cpp, revA"

#include <iostream.h>
#include <TrapezoidalTimeSeriesIntegrator.h>
#include <Vector.h>
#include <Channel.h>
#include <PathSeries.h>

TrapezoidalTimeSeriesIntegrator::TrapezoidalTimeSeriesIntegrator() 
  :TimeSeriesIntegrator(TIMESERIES_INTEGRATOR_TAG_Trapezoidal)
{

}

TrapezoidalTimeSeriesIntegrator::~TrapezoidalTimeSeriesIntegrator()
{

}

TimeSeries*
TrapezoidalTimeSeriesIntegrator::integrate(TimeSeries *theSeries, double delta)
{	
  // Check for zero time step, before dividing to get number of steps
  if (delta <= 0.0) {
    g3ErrorHandler->warning("TrapezoidalTimeSeriesIntegrator::integrate() %s %d %s\n",
			   "Attempting to integrate time step", delta, "<= 0");
    return 0;
   }

  // check a TimeSeries object was passed
  if (theSeries == 0) {
    g3ErrorHandler->warning("TrapezoidalTimeSeriesIntegrator::integrate() %s\n",
			     "- no TimeSeries passed");
    return 0;
  }

  // Add one to get ceiling out of type cast
  int numSteps = (int)(theSeries->getDuration()/delta + 1.0);

  Vector *theIntegratedValues = new Vector (numSteps);

  // Check that the Vector was allocated properly
  if (theIntegratedValues == 0 || theIntegratedValues->Size() == 0) {
    g3ErrorHandler->warning("TrapezoidalTimeSeriesIntegrator::integrate() %s %d\n",
			    "Ran out of memory allocating Vector of size ", numSteps);

    if (theIntegratedValues != 0)
      delete theIntegratedValues;

    return 0;
  }
  else {
    int i;                // Counter for indexing
    double dummyTime;     // Dummy variable for integrating
    double previousValue; // Temporary storage to avoid accessing same value twice
	                        // through identical method calls
    double currentValue;
      
    // Set the first point
    // Assuming initial condition is zero, i.e. F(0) = 0
    (*theIntegratedValues)[0] = theSeries->getFactor(0.0) * delta * 0.5;

    previousValue = (*theIntegratedValues)[0];
    
    dummyTime = delta;
    
    for (i = 1; i < numSteps; i++, dummyTime += delta) {
      currentValue = theSeries->getFactor(dummyTime);
      
      // Apply the trapezoidal rule to update the integrated value
      (*theIntegratedValues)[i] = (*theIntegratedValues)[i-1] +
	delta*0.5 * (currentValue + previousValue);
      
      previousValue = currentValue;
    }
    
    // Set the last point
    (*theIntegratedValues)[i] = (*theIntegratedValues)[i-1] +
      delta*0.5 * (theSeries->getFactor(dummyTime));
  }

  // Set the method return value
  PathSeries *returnSeries = new PathSeries (*theIntegratedValues, delta);

  if (returnSeries == 0) {
    g3ErrorHandler->warning("TrapezoidalTimeSeriesIntegrator::integrate() %s\n",
			    "Ran out of memory creating PathSeries");

    return 0;
   }

  return returnSeries;
}

int
TrapezoidalTimeSeriesIntegrator::sendSelf(int commitTag, Channel &theChannel)
{
   // Need to implement, return dummy value for now
   return 1;
}

int
TrapezoidalTimeSeriesIntegrator::recvSelf(int commitTag, Channel &theChannel, 
		           FEM_ObjectBroker &theBroker)
{
   // Need to implement, return dummy value for now
   return 1;
}

void
TrapezoidalTimeSeriesIntegrator::Print(ostream &s, int flag)
{
   // Need to implement, return for now
   return;
}
