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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/GroundMotionRecord.cpp,v $
                                                                        
                                                                        
// File: ~/earthquake/GroundMotionRecord.C
// 
// Written: fmk 
// Created: 05/98
// Revision: A
//
// Description: This file contains the class definition for 
// GroundMotionRecord. The units are (g).
//
// What: "@(#) GroundMotionRecord.C, revA"

#include <GroundMotionRecord.h>
#include <PathSeries.h>
#include <PathTimeSeries.h>
#include <stdlib.h>
#include <math.h>
#include <iostream.h>
#include <classTags.h>
#include <Vector.h>
#include <Channel.h>

GroundMotionRecord::GroundMotionRecord()
  :GroundMotion(GROUND_MOTION_TAG_GroundMotionRecord),
   theAccelTimeSeries(0), theVelTimeSeries(0), theDispTimeSeries(0),
   data(3), delta(0.0)
{

}

GroundMotionRecord::GroundMotionRecord(char *fileNameAccel,
				       double timeStep,
				       double theFactor,
				       double dT)
  :GroundMotion(GROUND_MOTION_TAG_GroundMotionRecord),
   theAccelTimeSeries(0), theVelTimeSeries(0), theDispTimeSeries(0),
   data(3), delta(dT)
{

  theAccelTimeSeries = new PathSeries(fileNameAccel, timeStep, theFactor);

  if (theAccelTimeSeries == 0) {
    g3ErrorHandler->warning("GroundMotionRecord::GroundMotionRecord() - %s\n",
			    "unable to create PathSeries");
  }  

}

GroundMotionRecord::GroundMotionRecord(char *fileNameAccel,
				       char *fileNameTime,
				       double theFactor,
				       double dT)
  :GroundMotion(GROUND_MOTION_TAG_GroundMotionRecord),
   theAccelTimeSeries(0), theVelTimeSeries(0), theDispTimeSeries(0),
   data(3), delta(dT)
{

  theAccelTimeSeries = new PathTimeSeries(fileNameAccel, fileNameTime, theFactor);

  if (theAccelTimeSeries == 0) {
    g3ErrorHandler->warning("GroundMotionRecord::GroundMotionRecord() - %s\n",
			    "unable to create PathSeries");
  }
}

GroundMotionRecord::~GroundMotionRecord()
{
  if (theAccelTimeSeries != 0)
    delete theAccelTimeSeries;
  if (theVelTimeSeries != 0)
    delete theVelTimeSeries;
  if (theDispTimeSeries != 0)
    delete theDispTimeSeries;
}

double 
GroundMotionRecord::getDuration(void)
{
  if (theAccelTimeSeries != 0)
    return theAccelTimeSeries->getDuration();
  else
    return 0.0;
}

double 
GroundMotionRecord::getPeakAccel(void)
{
  if (theAccelTimeSeries != 0)
    return theAccelTimeSeries->getPeakFactor();
  else
    return 0.0;

}

double 
GroundMotionRecord::getPeakVel(void)
{
  if (theVelTimeSeries != 0)
    return theVelTimeSeries->getPeakFactor();

  // if theAccel is not 0, integrate accel series to get a vel series
  else if (theAccelTimeSeries != 0) {
    theVelTimeSeries = this->integrate(theAccelTimeSeries, delta);
    if (theVelTimeSeries != 0)
      return theVelTimeSeries->getPeakFactor();      
    else
      return 0.0;
  }
  return 0.0;
}

double 
GroundMotionRecord::getPeakDisp(void)
{
  if (theDispTimeSeries != 0)
    return theDispTimeSeries->getPeakFactor();

  // if theVel is not 0, integrate vel series to get disp series
  else if (theVelTimeSeries != 0) {
    theDispTimeSeries = this->integrate(theVelTimeSeries, delta);
    if (theDispTimeSeries != 0)
      return theDispTimeSeries->getPeakFactor();      
    else
      return 0.0;
  }

  // if theAccel is not 0, integrate vel series to get disp series
  else if (theAccelTimeSeries != 0) {
    theVelTimeSeries = this->integrate(theAccelTimeSeries, delta);
    if (theVelTimeSeries != 0) {
      theDispTimeSeries = this->integrate(theVelTimeSeries, delta);
      if (theDispTimeSeries != 0)
	return theDispTimeSeries->getPeakFactor();      
      else
	return 0.0;
    } else
      return 0.0;
  }

  return 0.0;
}

double 
GroundMotionRecord::getAccel(double time)
{
  if (time < 0.0)
    return 0.0;
  
  if (theAccelTimeSeries != 0)
    return theAccelTimeSeries->getFactor(time);
  else
    return 0.0;
}     

double 
GroundMotionRecord::getVel(double time)
{
  if (time < 0.0)
    return 0.0;

  if (theVelTimeSeries != 0) 
    return theVelTimeSeries->getFactor(time);      
  
  // if theAccel is not 0, integrate accel series to get a vel series
  else if (theAccelTimeSeries != 0) {
    theVelTimeSeries = this->integrate(theAccelTimeSeries, delta);
    if (theVelTimeSeries != 0)
      return theVelTimeSeries->getFactor(time);      
    else
      return 0.0;
  }

  return 0.0;
}

double 
GroundMotionRecord::getDisp(double time)
{
  if (time < 0.0)
    return 0.0;

  if (theDispTimeSeries != 0)
    return theDispTimeSeries->getFactor(time);

  // if theVel is not 0, integrate vel series to get disp series
  else if (theVelTimeSeries != 0) {
    theDispTimeSeries = this->integrate(theVelTimeSeries, delta);
    if (theDispTimeSeries != 0)
      return theDispTimeSeries->getFactor(time);      
    else
      return 0.0;
  }

  // if theAccel is not 0, integrate vel series to get disp series
  else if (theAccelTimeSeries != 0) {
    theVelTimeSeries = this->integrate(theAccelTimeSeries, delta);
    if (theVelTimeSeries != 0) {
      theDispTimeSeries = this->integrate(theVelTimeSeries, delta);
      if (theDispTimeSeries != 0)
	return theDispTimeSeries->getFactor(time);      
      else
	return 0.0;
    } else
      return 0.0;
  }

  return 0.0;
}

const Vector &
GroundMotionRecord::getDispVelAccel(double time)
{
  if (time < 0.0) {
    data(0) = 0.0;
    data(1) = 0.0;
    data(2) = 0.0;
    return data;
  }

  if (theAccelTimeSeries != 0 && theVelTimeSeries != 0 && theDispTimeSeries != 0) {
    data(0) = theDispTimeSeries->getFactor(time);
    data(1) = theVelTimeSeries->getFactor(time);
    data(2) = theAccelTimeSeries->getFactor(time);
  } else {
    data(2) = this->getAccel(time);
    data(1) = this->getVel(time);
    data(0) = this->getDisp(time);
  }

  return data;
}


int 
GroundMotionRecord::sendSelf(int commitTag, Channel &theChannel)
{
  // FIX -- need to send the TimeSeries's, etc.	
  cerr << "GroundMotionRecord::sendSelf() -- not yet implemented" << endl;
  return -1;
  // return theChannel.sendVector(this->getDbTag(), commitTag, data);
}


int 
GroundMotionRecord::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // FIX -- need to receive the TimeSeries's, etc.
  cerr << "GroundMotionRecord::recvSelf() -- not yet impelemented" << endl;
  return -1;
  // theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  // return 0;
}
