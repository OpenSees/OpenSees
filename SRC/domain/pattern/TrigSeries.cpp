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
                                                                        
// $Revision: 1.4 $
// $Date: 2010-04-06 20:16:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/TrigSeries.cpp,v $

// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class implementation of TrigSeries.
//
// What: "@(#) TrigSeries.C, revA"


#include <TrigSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>

#include <math.h>

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_NewTrigSeries(void)
{
  // Pointer to a uniaxial material that will be returned
  TimeSeries *theSeries = 0;
  
  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  
  if (numRemainingArgs < 3) {
    opserr << "WARNING: invalid num args TrigSeries <tag?> $tStart $tFinish $period <-shift shift> <-factor cFactor?>\n";
    return 0;
  }

  int tag = 0;     // default tag = 0
  double dData[5];
  dData[3] = 0.0; // default shift = 0.0
  dData[4] = 1.0; // default cFactor = 1.0
  int numData = 0;

  // get tag if provided
  if (numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
      opserr << "WARNING invalid series tag in Trig tag?" << endln;
      return 0;
    }
    numRemainingArgs -= 1;
  }
  
  numData = 3;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid double data for TrigSeries with tag: " << tag << endln;
    return 0;
  }    
  numRemainingArgs -= 3;


  while (numRemainingArgs > 1) {
    char argvS[10];
    if (OPS_GetString(argvS, 10) != 0) {
      opserr << "WARNING invalid string in Trig <tag?> <-factor cFactor?>" << endln;
      return 0;
    } 
  
    if (strcmp(argvS,"-shift") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &dData[3]) != 0) {
	opserr << "WARNING invalid shift in Trig Series with tag?" << tag << endln;
	return 0;
      }
    } else if (strcmp(argvS,"-factor") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &dData[4]) != 0) {
	opserr << "WARNING invalid shift in Trig Series with tag?" << tag << endln;
	return 0;
      }
    } else {
      opserr << "WARNING unknown option: " << argvS << "  in Rectangular Series with tag?" << tag << endln;      
      return 0;
    }      
    numRemainingArgs -= 2;
  }
  
  theSeries = new TrigSeries(tag, dData[0], dData[1], dData[2], dData[3], dData[4]);

  if (theSeries == 0) {
    opserr << "WARNING ran out of memory creating TrigSeries with tag: " << tag << "\n";
    return 0;
  }

  return theSeries;
}


TrigSeries::TrigSeries(int tag,
		       double startTime, 
		       double finishTime,
		       double T, 
		       double phi, 
		       double theFactor)
  :TimeSeries(tag, TSERIES_TAG_TrigSeries),
   tStart(startTime),tFinish(finishTime),
   period(T),shift(phi),cFactor(theFactor)
{
  if (period == 0.0) {
    opserr << "TrigSeries::TrigSeries -- input period is zero, setting period to PI\n";
    period = 2*asin(1.0);
  }
}


TrigSeries::TrigSeries()
  :TimeSeries(TSERIES_TAG_TrigSeries),
   tStart(0.0),tFinish(0.0),period(1.0),shift(0.0),cFactor(1.0)
{
  // does nothing
}


TrigSeries::~TrigSeries()
{
  // does nothing
}

TimeSeries *
TrigSeries::getCopy(void) {
  return new TrigSeries(this->getTag(), tStart, tFinish, period, shift, cFactor);
}

  
double
TrigSeries::getFactor(double pseudoTime)
{
  static double twopi = 4*asin(1.0);
  
  if (pseudoTime >= tStart && pseudoTime <= tFinish)
    return cFactor*sin(twopi*(pseudoTime-tStart)/period + shift);
  else
    return 0.0;
}

int
TrigSeries::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  Vector data(5);
  data(0) = cFactor;
  data(1) = tStart;	
  data(2) = tFinish;
  data(3) = period;
  data(4) = shift;
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "TrigSeries::sendSelf() - channel failed to send data\n";
    return result;
  }
  return 0;
}


int 
TrigSeries::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(5);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "TrigSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    tStart= 0.0;

    tFinish = 0.0;
	period = 1.0;
	shift = 0.0;
    return result;
  }
  cFactor = data(0);
  tStart = data(1);
  tFinish = data(2);
  period = data(3);
  shift = data(4);

  return 0;    
}

void
TrigSeries::Print(OPS_Stream &s, int flag)
{
    s << "Trig Series" << endln;
	s << "\tFactor: " << cFactor << endln;
	s << "\ttStart: " << tStart << endln;
	s << "\ttFinish: " << tFinish << endln;
	s << "\tPeriod: " << period << endln;
	s << "\tPhase Shift: " << shift << endln;
}
