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

// $Revision$
// $Date$
// $URL$

// Written: Andreas Schellenberg
// Created: 05/15
// Revision: A
//
// Description: This file contains the class definition for 
// a SimpsonTimeSeriesIntegrator, which integrates a
// ground motion TimeSeries using the trapezoidal rule.


#include <SimpsonTimeSeriesIntegrator.h>
#include <Vector.h>
#include <Channel.h>
#include <PathSeries.h>

void* OPS_SimpsonTimeSeriesIntegrator()
{
  return new SimpsonTimeSeriesIntegrator();
}

SimpsonTimeSeriesIntegrator::SimpsonTimeSeriesIntegrator()
  : TimeSeriesIntegrator(TIMESERIES_INTEGRATOR_TAG_Simpson)
{
  
}


SimpsonTimeSeriesIntegrator::~SimpsonTimeSeriesIntegrator()
{
  
}


TimeSeries* SimpsonTimeSeriesIntegrator::integrate(TimeSeries *theSeries, double delta)
{
  // check for zero time step, before dividing to get number of steps
  if (delta <= 0.0)  {
    opserr << "SimpsonTimeSeriesIntegrator::integrate() - attempting to integrate using time step "
      << delta << "<= 0.0.\n";
    return 0;
  }
  
  // check a TimeSeries object was passed
  if (theSeries == 0)  {
    opserr << "SimpsonTimeSeriesIntegrator::integrate() - no TimeSeries passed.\n";
    return 0;
  }
  
  // add one to get ceiling out of type cast
  long long numSteps = (long long)(theSeries->getDuration()/delta + 1.0);
  
  // create new vector for integrated values
  Vector *theInt = new Vector(numSteps);
  
  // check that the Vector was allocated properly
  if (theInt == 0 || theInt->Size() == 0)  {
    opserr << "SimpsonTimeSeriesIntegrator::integrate() - ran out of memory allocating Vector " << endln;

    if (theInt != 0)
      delete theInt;
    
    return 0;
  }
  
  double dummyTime;
  double fi, fj, fk;  //function values
  double Fi, Fj, Fk;  //intergral values

  dummyTime = theSeries->getStartTime();

  Fi = 0.0;
  Fj = 0.0;

  fi = 0.0;
  fj = 0.0;

  for (long long i = 0; i < numSteps; i++, dummyTime += delta)  {

    fk = theSeries->getFactor(dummyTime);
    // Apply the Simpson's rule to update the integral
    Fk = Fi + delta / 3.0 * (fi + 4.0 * fj + fk);

    (*theInt)[i] = Fk;

    fi = fj;
    fj = fk;

    Fi = Fj;
    Fj = Fk;
  }

  // set the method return value
  PathSeries *returnSeries = new PathSeries(0, *theInt, delta, true);
  delete theInt;
  
  if (returnSeries == 0)  {
    opserr << "SimpsonTimeSeriesIntegrator::integrate() - ran out of memory creating PathSeries.\n";
    return 0;
  }
  
  return returnSeries;
}


TimeSeries* SimpsonTimeSeriesIntegrator::differentiate(TimeSeries *theSeries, double delta)
{
  // check for zero time step, before dividing to get number of steps
  if (delta <= 0.0)  {
    opserr << "SimpsonTimeSeriesIntegrator::differentiate() - attempting to differentiate using time step "
      << delta << "<= 0.0.\n";
    return 0;
  }
  
  // check a TimeSeries object was passed
  if (theSeries == 0)  {
    opserr << "SimpsonTimeSeriesIntegrator::differentiate() - no TimeSeries passed.\n";
    return 0;
  }
  
  // add one to get ceiling out of type cast
  long long numSteps = (long long)(theSeries->getDuration()/delta + 1.0);
  
  // create new vector for integrated values
  Vector *theDif = new Vector(numSteps);
  
  // check that the Vector was allocated properly
  if (theDif == 0 || theDif->Size() == 0)  {
    opserr << "SimpsonTimeSeriesIntegrator::differentiate() - ran out of memory allocating Vector " << endln;

    if (theDif != 0)
      delete theDif;
    
    return 0;
  }
  
  double dummyTime;
  double Fi, Fj, Fk;  //function values
  double fi, fj, fk;  //derivative values

  dummyTime = theSeries->getStartTime();

  Fi = 0.0;
  Fj = 0.0;

  fi = 0.0;
  fj = 0.0;

  for (long long i = 0; i < numSteps; i++, dummyTime += delta)  {

    Fk = theSeries->getFactor(dummyTime);
    // Apply the Simpson's rule to update the derivative
    fk = 3.0 * (Fk - Fi) / delta - fi - 4.0 * fj;

    (*theDif)[i] = fk;

    fi = fj;
    fj = fk;

    Fi = Fj;
    Fj = Fk;
  }

  // set the method return value
  PathSeries *returnSeries = new PathSeries(0, *theDif, delta, true);
  delete theDif;
  
  if (returnSeries == 0)  {
    opserr << "SimpsonTimeSeriesIntegrator::differentiate() - ran out of memory creating PathSeries.\n";
    return 0;
  }
  
  return returnSeries;
}


int SimpsonTimeSeriesIntegrator::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}


int SimpsonTimeSeriesIntegrator::recvSelf(int commitTag, Channel &theChannel, 
  FEM_ObjectBroker &theBroker)
{
  return 0;
}


void SimpsonTimeSeriesIntegrator::Print(OPS_Stream &s, int flag)
{
  return;
}
