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
                                                                        
// $Revision: 1.6 $
// $Date: 2010-02-04 00:37:42 $
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
#include <classTags.h>
#include <Vector.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

GroundMotionRecord::GroundMotionRecord()
  :GroundMotion(GROUND_MOTION_TAG_GroundMotionRecord),
   theAccelTimeSeries(0), theVelTimeSeries(0), theDispTimeSeries(0),
   data(3), delta(0.0)
{

}

GroundMotionRecord::GroundMotionRecord(const char *fileNameAccel,
				       double timeStep,
				       double theFactor,
				       double dT)
  :GroundMotion(GROUND_MOTION_TAG_GroundMotionRecord),
   theAccelTimeSeries(0), theVelTimeSeries(0), theDispTimeSeries(0),
   data(3), delta(dT)
{

  theAccelTimeSeries = new PathSeries(0, fileNameAccel, timeStep, theFactor);

  if (theAccelTimeSeries == 0) {
    opserr << "GroundMotionRecord::GroundMotionRecord() - unable to create PathSeries\n";
  }  

}

GroundMotionRecord::GroundMotionRecord(const char *fileNameAccel,
				       const char *fileNameTime,
				       double theFactor,
				       double dT)
  :GroundMotion(GROUND_MOTION_TAG_GroundMotionRecord),
   theAccelTimeSeries(0), theVelTimeSeries(0), theDispTimeSeries(0),
   data(3), delta(dT)
{

  theAccelTimeSeries = new PathTimeSeries(0, fileNameAccel, fileNameTime, theFactor);

  if (theAccelTimeSeries == 0) {
    opserr << "GroundMotionRecord::GroundMotionRecord() - unable to create PathSeries\n";
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
  int dbTag = this->getDbTag();

  static ID idData(6);
  
  if (theAccelTimeSeries != 0) {
    idData(0) = theAccelTimeSeries->getClassTag();
    int seriesDbTag = theAccelTimeSeries->getDbTag();
    if (seriesDbTag == 0) {
      seriesDbTag = theChannel.getDbTag();
      theAccelTimeSeries->setDbTag(seriesDbTag);
    }
    idData(1) = seriesDbTag;
  } else
    idData(0) = -1;

  if (theVelTimeSeries != 0) {
    idData(2) = theVelTimeSeries->getClassTag();
    int seriesDbTag = theVelTimeSeries->getDbTag();
    if (seriesDbTag == 0) {
      seriesDbTag = theChannel.getDbTag();
      theVelTimeSeries->setDbTag(seriesDbTag);
    }
    idData(3) = seriesDbTag;
  } else
    idData(2) = -1;

  if (theDispTimeSeries != 0) {
    idData(4) = theDispTimeSeries->getClassTag();
    int seriesDbTag = theDispTimeSeries->getDbTag();
    if (seriesDbTag == 0) {
      seriesDbTag = theChannel.getDbTag();
      theDispTimeSeries->setDbTag(seriesDbTag);
    }
    idData(5) = seriesDbTag;
  } else
    idData(4) = -1;

  int res = theChannel.sendID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "GroundMotionRecord::sendSelf() - channel failed to send data\n";
    return res;
  }

  // now send the series
  if (theAccelTimeSeries != 0) {
    res = theAccelTimeSeries->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "GroundMotionRecord::sendSelf - failed to send accel series\n";
      return res;
    }
  }

  if (theVelTimeSeries != 0) {
    res = theVelTimeSeries->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "GroundMotionRecord::sendSelf - failed to send velocity series\n";
      return res;
    }
  }

  if (theDispTimeSeries != 0) {
    res = theDispTimeSeries->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "GroundMotionRecord::sendSelf - failed to send disp series\n";
      return res;
    }
  }

  return 0;
}


int 
GroundMotionRecord::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  static ID idData(6);
  int res = theChannel.recvID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "GroundMotionRecord::recvSelf() - channel failed to send data\n";
    return res;
  }

  int seriesClassTag = idData(0);
  if (seriesClassTag != -1) {
    int seriesDbTag = idData(1);
    if (theAccelTimeSeries == 0 || theAccelTimeSeries->getClassTag() != seriesClassTag) {
      if (theAccelTimeSeries != 0)
	delete theAccelTimeSeries;
      theAccelTimeSeries = theBroker.getNewTimeSeries(seriesClassTag);
      if (theAccelTimeSeries == 0) {
	opserr << "GroundMotionRecord::recvSelf - could not create a TimeSeries object\n";
	return -2;
      }
    }
    theAccelTimeSeries->setDbTag(seriesDbTag);
    res = theAccelTimeSeries->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "GroundMotionRecord::recvSelf() - accel series failed to send data\n";
      return res;
    }
  }

  seriesClassTag = idData(2);
  if (seriesClassTag != -1) {
    int seriesDbTag = idData(3);
    if (theVelTimeSeries == 0 || theVelTimeSeries->getClassTag() != seriesClassTag) {
      if (theVelTimeSeries != 0)
	delete theVelTimeSeries;
      theVelTimeSeries = theBroker.getNewTimeSeries(seriesClassTag);
      if (theVelTimeSeries == 0) {
	opserr << "GroundMotionRecord::recvSelf - could not create a TimeSeries object\n";
	return -2;
      }
    }
    theVelTimeSeries->setDbTag(seriesDbTag);
    res = theVelTimeSeries->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "GroundMotionRecord::recvSelf() - accel series failed to send data\n";
      return res;
    }
  }

  seriesClassTag = idData(4);
  if (seriesClassTag != -1) {
    int seriesDbTag = idData(5);
    if (theDispTimeSeries == 0 || theDispTimeSeries->getClassTag() != seriesClassTag) {
      if (theDispTimeSeries != 0)
	delete theDispTimeSeries;
      theDispTimeSeries = theBroker.getNewTimeSeries(seriesClassTag);
      if (theDispTimeSeries == 0) {
	opserr << "GroundMotionRecord::recvSelf - could not create a TimeSeries object\n";
	return -2;
      }
    }
    theDispTimeSeries->setDbTag(seriesDbTag);
    res = theDispTimeSeries->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "GroundMotionRecord::recvSelf() - accel series failed to send data\n";
      return res;
    }
  }

  return 0;
}

