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
                                                                        
// $Revision: 1.3 $
// $Date: 2010-02-04 00:34:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/RectangularSeries.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class implementation of RectangularSeries.
//
// What: "@(#) RectangularSeries.C, revA"


#include <RectangularSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <Parameter.h>

#include <elementAPI.h>

void *
OPS_RectangularSeries(void)
{
  // Pointer to a uniaxial material that will be returned
  TimeSeries *theSeries = 0;
  
  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  
  if (numRemainingArgs < 2) {
    opserr << " Rectangular <tag?> tStart tFinish <-factor cFactor>\n";
    return 0;
  }

  int tag = 0;     // default tag = 0
  double dData[3];
  dData[2] = 1.0; // default cFactor = 1.0
  int numData = 0;

  // get tag if provided
  if (numRemainingArgs == 3 || numRemainingArgs == 5) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
      opserr << "WARNING invalid series tag in Rectangular tag? tStart tFinish <-factor cFactor>\n";
      return 0;
    }
    numRemainingArgs -= 1;
  }
  
  numData = 2;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid double data for RectangularSeries with tag: " << tag << endln;
    return 0;
  }    
  numRemainingArgs -= 2;

  while (numRemainingArgs > 1) {
    const char *argvS = OPS_GetString();

    if (strcmp(argvS,"-factor") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &dData[2]) != 0) {
	opserr << "WARNING invalid shift in Trig Series with tag?" << tag << endln;
	return 0;
      }
    } else {
      opserr << "WARNING unknown option: " << argvS << "  in Rectangular Series with tag?" << tag << endln;      
      return 0;
    }      
    numRemainingArgs -= 2;
  }

  theSeries = new RectangularSeries(tag, dData[0], dData[1], dData[2]);

  if (theSeries == 0) {
    opserr << "WARNING ran out of memory creating RectangularSeries with tag: " << tag << "\n";
    return 0;
  }

  return theSeries;
}


RectangularSeries::RectangularSeries(int tag, 
				     double startTime, 
				     double finishTime,
				     double theFactor)
  :TimeSeries(tag, TSERIES_TAG_RectangularSeries),
   tStart(startTime),tFinish(finishTime),cFactor(theFactor)
{
  // does nothing
}


RectangularSeries::RectangularSeries()
  :TimeSeries(TSERIES_TAG_RectangularSeries),
   tStart(0),tFinish(0),cFactor(0)
{
  // does nothing
}


RectangularSeries::~RectangularSeries()
{
  // does nothing
}

TimeSeries *
RectangularSeries::getCopy(void) {
  return new RectangularSeries(this->getTag(), tStart, tFinish, cFactor);
}

double
RectangularSeries::getFactor(double pseudoTime)
{	
  if (pseudoTime >= tStart && pseudoTime <= tFinish)
    return cFactor;
  else
    return 0;
}

int
RectangularSeries::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  Vector data(3);
  data(0) = cFactor;
  data(1) = tStart;	
  data(2) = tFinish;
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "RectangularSeries::sendSelf() - channel failed to send data\n";
    return result;
  }
  return 0;
}


int 
RectangularSeries::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(3);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "RectangularSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    tStart= 0.0;
    tFinish = 0.0;
    return result;
  }
  cFactor = data(0);
  tStart = data(1);
  tFinish = data(2);

  return 0;    
}

void
RectangularSeries::Print(OPS_Stream &s, int flag)
{
    s << "Linear Series: constant factor: " << cFactor;
    s << "  tStart: " << tStart << "  tFinish: " << tFinish << endln;

}

double
RectangularSeries::getFactorSensitivity(double pseudoTime)
{
  if (pseudoTime >= tStart && pseudoTime <= tFinish) {
    if (parameterID == 1)
      return 1.0;
  }

  return 0.0;
}

int 
RectangularSeries::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strncmp(argv[0],"factor",80) == 0) {
    param.setValue(cFactor);
    return param.addObject(1, this);
  }

  return -1;
}
   
int 
RectangularSeries::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    cFactor = info.theDouble;
    return 0;
  }

  return -1;
}

int
RectangularSeries::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}
