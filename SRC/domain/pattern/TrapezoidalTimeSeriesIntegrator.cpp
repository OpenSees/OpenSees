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
                                                                        
// $Revision: 1.5 $
// $Date: 2010-02-04 00:36:46 $
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

#include <TrapezoidalTimeSeriesIntegrator.h>
#include <Vector.h>
#include <Channel.h>
#include <PathSeries.h>

void* OPS_TrapezoidalTimeSeriesIntegrator()
{
    return new TrapezoidalTimeSeriesIntegrator();
}

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
    opserr << "TrapezoidalTimeSeriesIntegrator::integrate() Attempting to integrate usiing time step" <<
      delta << "<= 0\n";
    return 0;
   }

  // check a TimeSeries object was passed
  if (theSeries == 0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::integrate() - - no TimeSeries passed\n";
    return 0;
  }

  // Add one to get ceiling out of type cast
  long long numSteps = (long long)theSeries->getDuration() / delta + 1.0;

  Vector *theInt = new Vector (numSteps);

  // Check that the Vector was allocated properly
  if (theInt == 0 || theInt->Size() == 0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::integrate() Ran out of memory allocating Vector " << endln;


    if (theInt != 0)
      delete theInt;

    return 0;
  }

  double dummyTime;     // Dummy variable for integrating
  double fi, fj;        //function values
  double F;             //intergral value

  dummyTime = theSeries->getStartTime();

  F = 0.0;

  fi = 0.0;

  for (long long i = 0; i < numSteps; i++, dummyTime += delta) {
    fj = theSeries->getFactor(dummyTime);

    // Apply the trapezoidal rule to update the integral
    F = F + 0.5 * delta * (fi + fj);

    (*theInt)[i] = F;

    fi = fj;
  }

  // Set the method return value
  PathSeries *returnSeries = new PathSeries (0, *theInt, delta, true);
  delete theInt;

  if (returnSeries == 0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::integrate() Ran out of memory creating PathSeries\n";

    return 0;
   }

  return returnSeries;
}

TimeSeries*
TrapezoidalTimeSeriesIntegrator::differentiate(TimeSeries *theSeries, double delta)
{	
  // Check for zero time step, before dividing to get number of steps
  if (delta <= 0.0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::differentiate() Attempting to differentiate using time step" <<
      delta << "<= 0\n";
    return 0;
   }

  // check a TimeSeries object was passed
  if (theSeries == 0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::differentiate() - - no TimeSeries passed\n";
    return 0;
  }

  // Add one to get ceiling out of type cast
  long long numSteps = (long long)theSeries->getDuration() / delta + 1.0;

  Vector *theDif = new Vector (numSteps);

  // Check that the Vector was allocated properly
  if (theDif == 0 || theDif->Size() == 0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::differentiate() Ran out of memory allocating Vector " << endln;

    if (theDif != 0)
      delete theDif;

    return 0;
  }

  double dummyTime;     // Dummy variable for integrating
  double Fi, Fj;        // function values
  double f;             // derivative value
      
  dummyTime = theSeries->getStartTime();

  f = 0.0;

  Fi = 0.0;

  for (long long i = 0; i < numSteps; i++, dummyTime += delta) {
    Fj = theSeries->getFactor(dummyTime);

    // Apply the trapezoidal rule to update the derivative
    f = 2.0 * (Fj - Fi) / delta - f;

    (*theDif)[i] = f;

    Fi = Fj;
  }

  // Set the method return value
  PathSeries *returnSeries = new PathSeries (0, *theDif, delta, true);
  delete theDif;

  if (returnSeries == 0) {
    opserr << "TrapezoidalTimeSeriesIntegrator::differentiate() Ran out of memory creating PathSeries\n";

    return 0;
   }

  return returnSeries;
}

int
TrapezoidalTimeSeriesIntegrator::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int
TrapezoidalTimeSeriesIntegrator::recvSelf(int commitTag, Channel &theChannel, 
		           FEM_ObjectBroker &theBroker)
{
  return 0;
}

void
TrapezoidalTimeSeriesIntegrator::Print(OPS_Stream &s, int flag)
{
   // Need to implement, return for now
   return;
}
