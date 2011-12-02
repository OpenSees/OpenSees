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
                                                                        
// $Revision: 1.2 $
// $Date: 2000-12-12 07:30:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/GroundMotion.cpp,v $
                                                                        
                                                                        
// File: ~/earthquake/GroundMotion.C
//
// Written: fmk 
// Created: 05/98
// Revision: A

#include <iostream.h>

#include <GroundMotion.h>
#include <TimeSeriesIntegrator.h>
#include <G3Globals.h>


GroundMotion::GroundMotion(TimeSeries *dispSeries, 
			   TimeSeries *velSeries, 
			   TimeSeries *accelSeries,
			   TimeSeriesIntegrator *theIntegratr,
			   double dTintegration)
:MovableObject(GROUND_MOTION_TAG_GroundMotion), 
 theAccelSeries(accelSeries), theVelSeries(velSeries), 
 theDispSeries(dispSeries), theIntegrator(theIntegratr),
 data(3), delta(dTintegration)
{

}



GroundMotion::GroundMotion(int theClassTag)
:MovableObject(theClassTag), 
 theAccelSeries(0), theVelSeries(0), theDispSeries(0), theIntegrator(0),
 data(3), delta(0.0)
{

}

GroundMotion::~GroundMotion()
{
  if (theAccelSeries != 0)
    delete theAccelSeries;
  if (theVelSeries != 0)
    delete theVelSeries;
  if (theDispSeries != 0)
    delete theDispSeries;
  if (theIntegrator != 0)
    delete theIntegrator;
}


void
GroundMotion::setIntegrator(TimeSeriesIntegrator *integrator)
{
  if (theIntegrator != 0)
    delete theIntegrator;

  theIntegrator = integrator;
}

TimeSeries*
GroundMotion::integrate(TimeSeries *theSeries, double delta)
{
  // check that an integrator & accel series exist
  if(theIntegrator == 0) {
    return 0;
  }

  if(theSeries == 0) {
    g3ErrorHandler->warning("GroundMotion::integrate %s\n",
			    "- no TimeSeries specified");
    return 0;
  }

  // integrate the series, if no vel series exists set it to new one
  TimeSeries *theNewSeries = theIntegrator->integrate(theSeries, delta);


  if(theNewSeries == 0) {
    g3ErrorHandler->warning("GroundMotion::integrate %s\n",
			    "- no TimeSeriesIntegrator failed to integrate");
    return 0;
  }

  return theNewSeries;
}

double 
GroundMotion::getDuration(void)
{
  if (theAccelSeries != 0)
    return theAccelSeries->getDuration();
  else
    return 0.0;
}

double 
GroundMotion::getPeakAccel(void)
{
  if (theAccelSeries != 0)
    return theAccelSeries->getPeakFactor();
  else
    return 0.0;

}

double 
GroundMotion::getPeakVel(void)
{
  if (theVelSeries != 0)
    return theVelSeries->getPeakFactor();

  // if theAccel is not 0, integrate accel series to get a vel series
  else if (theAccelSeries != 0) {
    theVelSeries = this->integrate(theAccelSeries, delta);
    if (theVelSeries != 0)
      return theVelSeries->getPeakFactor();      
    else
      return 0.0;
  }
  return 0.0;
}

double 
GroundMotion::getPeakDisp(void)
{
  if (theDispSeries != 0)
    return theDispSeries->getPeakFactor();

  // if theVel is not 0, integrate vel series to get disp series
  else if (theVelSeries != 0) {
    theDispSeries = this->integrate(theVelSeries, delta);
    if (theDispSeries != 0)
      return theDispSeries->getPeakFactor();      
    else
      return 0.0;
  }

  // if theAccel is not 0, integrate vel series to get disp series
  else if (theAccelSeries != 0) {
    theVelSeries = this->integrate(theAccelSeries, delta);
    if (theVelSeries != 0) {
      theDispSeries = this->integrate(theVelSeries, delta);
      if (theDispSeries != 0)
	return theDispSeries->getPeakFactor();      
      else
	return 0.0;
    } else
      return 0.0;
  }

  return 0.0;
}

double 
GroundMotion::getAccel(double time)
{
  if (time < 0.0)
    return 0.0;
  
  if (theAccelSeries != 0)
    return theAccelSeries->getFactor(time);
  else
    return 0.0;
}     

double 
GroundMotion::getVel(double time)
{
  if (time < 0.0)
    return 0.0;

  if (theVelSeries != 0) 
    return theVelSeries->getFactor(time);      
  
  // if theAccel is not 0, integrate accel series to get a vel series
  else if (theAccelSeries != 0) {
    theVelSeries = this->integrate(theAccelSeries, delta);
    if (theVelSeries != 0)
      return theVelSeries->getFactor(time);      
    else
      return 0.0;
  }

  return 0.0;
}

double 
GroundMotion::getDisp(double time)
{
  if (time < 0.0)
    return 0.0;

  if (theDispSeries != 0)
    return theDispSeries->getFactor(time);


  // if theVel is not 0, integrate vel series to get disp series
  else if (theVelSeries != 0) {
    theDispSeries = this->integrate(theVelSeries, delta);
    if (theDispSeries != 0)
      return theDispSeries->getFactor(time);      
    else
      return 0.0;
  }

  // if theAccel is not 0, integrate vel series to get disp series
  else if (theAccelSeries != 0) {
    theVelSeries = this->integrate(theAccelSeries, delta);
    if (theVelSeries != 0) {
      theDispSeries = this->integrate(theVelSeries, delta);
      if (theDispSeries != 0)
	return theDispSeries->getFactor(time);      
      else
	return 0.0;
    } else
      return 0.0;
  }

  return 0.0;
}

const Vector &
GroundMotion::getDispVelAccel(double time)
{
  if (time < 0.0) {
    data(0) = 0.0;
    data(1) = 0.0;
    data(2) = 0.0;
    return data;
  }

  if (theAccelSeries != 0 && theVelSeries != 0 && theDispSeries != 0) {
    data(0) = theDispSeries->getFactor(time);
    data(1) = theVelSeries->getFactor(time);
    data(2) = theAccelSeries->getFactor(time);
  } else {
    data(2) = this->getAccel(time);
    data(1) = this->getVel(time);
    data(0) = this->getDisp(time);
  }

  return data;
}


int 
GroundMotion::sendSelf(int commitTag, Channel &theChannel)
{
  cerr << "GroundMotion::sendSelf() -- not yet implemented" << endl;
  return -1;
}


int 
GroundMotion::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  cerr << "GroundMotion::recvSelf() -- not yet impelemented" << endl;
  return -1;
}





