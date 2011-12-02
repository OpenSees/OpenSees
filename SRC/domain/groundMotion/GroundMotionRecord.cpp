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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:18 $
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

GroundMotionRecord::GroundMotionRecord(double theFactor)
:GroundMotion(GROUND_MOTION_TAG_GroundMotionRecord)
{
  //cerr << "GroundMotionRecord::GroundMotionRecord() -";
  //cerr << " does not determine Vel or Disp - just returns 0.0\n";

  // Probably don't need this constructor since it makes little sense
  // not to construct a GMR without some data
  // Will create default series for now though
  //theAccelSeries = new PathSeries ();
  //theVelSeries = new PathSeries ();
  //theDispSeries = new PathSeries ();

  theAccelSeries = 0;
  theVelSeries = 0;
  theDispSeries = 0;
}

GroundMotionRecord::GroundMotionRecord (char *fileNameAccel,
										  double timeStep,
										  double theFactor)
:GroundMotion(GROUND_MOTION_TAG_GroundMotionRecord)
{

  //cerr << "GroundMotionRecord::GroundMotionRecord() -";
  //cerr << " does not determine Vel or Disp - just returns 0.0\n";
  
  theAccelSeries = new PathSeries (fileNameAccel, timeStep, theFactor);
  
  // Integrate to get velocities and displacements
  // Call default constructors for now, will change later
  // theVelSeries = new PathSeries ();
  // theDispSeries = new PathSeries ();
  theVelSeries = 0;
  theDispSeries = 0;
  
  if (theAccelSeries == 0 /* || theVelSeries == 0 || theDispSeries == 0 */)
  {
	cerr << "GroundMotionRecord::GroundMotionRecord() - "
		 << "unable to create PathSeries" << endl;
  }
}

GroundMotionRecord::GroundMotionRecord (char *fileNameAccel,
										  char *fileNameTime,
										  double theFactor)
:GroundMotion(GROUND_MOTION_TAG_GroundMotionRecord)
{
  //cerr << "GroundMotionRecord::GroundMotionRecord() -";
  //cerr << " does not determine Vel or Disp - just returns 0.0\n";

  theAccelSeries = new PathTimeSeries (fileNameAccel, fileNameTime, theFactor);

  // Integrate to get velocities and displacements
  // Call default constructors for now, will change later
  // theVelSeries = new PathTimeSeries ();
  // theDispSeries = new PathTimeSeries ();
  theVelSeries = 0;
  theDispSeries = 0;

  if (theAccelSeries == 0 /* || theVelSeries == 0 || theDispSeries == 0 */)
  {
	cerr << "GroundMotionRecord::GroundMotionRecord() - "
		 << "unable to create PathTimeSeries" << endl;
  }
}

GroundMotionRecord::~GroundMotionRecord()
{
  if (theAccelSeries != 0)
    delete theAccelSeries;
  if (theVelSeries != 0)
    delete theVelSeries;
  if (theDispSeries != 0)
    delete theDispSeries;

}

double 
GroundMotionRecord::getDuration(void)
{
  return theAccelSeries->getDuration();
}

double 
GroundMotionRecord::getPeakAccel(void)
{
  return theAccelSeries->getPeakFactor();
}

double 
GroundMotionRecord::getPeakVel(void)
{
  // return theVelSeries->getPeakFactor();
  return 0.0;
}

double 
GroundMotionRecord::getPeakDisp(void)
{
  // return theDispSeries->getPeakFactor();
  return 0.0;
}

double 
GroundMotionRecord::getAccel(double time)
{
  if (time < 0.0)
    return 0.0;
  else
	return theAccelSeries->getFactor(time);
}     

double 
GroundMotionRecord::getVel(double time)
{
  //if (time < 0.0)
  //  return 0.0;
  //else
  //  return theVelSeries->getFactor(time);
  return 0.0;
}

double 
GroundMotionRecord::getDisp(double time)
{
  //if (time < 0.0)
  //  return 0.0;
  //else
  //  return theDispSeries->getFactor(time);
  return 0.0;
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
