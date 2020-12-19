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
// $Date: 2010-02-04 00:34:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/LinearSeries.cpp,v $


// Written: fmk
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for LinearSeries.
// LinearSeries is a concrete class. A LinearSeries object provides
// a linear time series. the factor is given by the pseudoTime and
// a constant factor provided in the constructor.
//
// What: "@(#) LinearSeries.C, revA"


#include <LinearSeries.h>
#include <Vector.h>
#include <Channel.h>

#include <elementAPI.h>

void *
OPS_LinearSeries(void)
{
  // Pointer to a uniaxial material that will be returned
  TimeSeries *theSeries = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  int tag = 0;
  double cFactor = 1.0;
  double tStart = 0.0;
  int numData = 0;

  if (numRemainingArgs != 0) {

    if (numRemainingArgs == 1 || numRemainingArgs == 3 || numRemainingArgs == 5) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &tag) != 0) {
	opserr << "WARNING invalid series tag in LinearSeries tag? <-factor factor?> <-tStart tStart?" << endln;
	return 0;
      }
      numRemainingArgs--;
    }

    while (numRemainingArgs > 1) {
      const char *argvS = OPS_GetString();
      if (argvS == 0) {
	opserr << "WARNING string error in LinearSeries with tag: " << tag << endln;
	return 0;
      }
      if (strcmp(argvS,"-factor") == 0) {
	numData = 1;
	if (OPS_GetDouble(&numData, &cFactor) != 0) {
	  opserr << "WARNING invalid factor in LinearSeries with tag: " << tag << endln;
	  return 0;
	}
	numRemainingArgs -= 2;
      }
      if (strcmp(argvS,"-tStart") == 0) {
	numData = 1;
	if (OPS_GetDouble(&numData, &tStart) != 0) {
	  opserr << "WARNING invalid tStart in LinearSeries with tag: " << tag << endln;
	  return 0;
	}
	numRemainingArgs -= 2;
      }      
    }
  }

  theSeries = new LinearSeries(tag, cFactor, tStart);

  if (theSeries == 0) {
    opserr << "WARNING ran out of memory creating ConstantTimeSeries with tag: " << tag << "\n";
    return 0;
  }

  return theSeries;
}


LinearSeries::LinearSeries(int tag, double theFactor, double startTime)
  :TimeSeries(tag, TSERIES_TAG_LinearSeries),
   cFactor(theFactor), tStart(startTime)
{
  // does nothing
}


LinearSeries::~LinearSeries()
{
  // does nothing
}

TimeSeries *
LinearSeries::getCopy(void) {
  return new LinearSeries(this->getTag(), cFactor, tStart);
}

double
LinearSeries::getFactor(double pseudoTime)
{
  if (pseudoTime < tStart)
    return 0.0;
  else
    return cFactor*(pseudoTime-tStart);
}


int
LinearSeries::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  Vector data(2);
  data(0) = cFactor;
  data(1) = tStart;  
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "LinearSeries::sendSelf() - channel failed to send data\n";
    return result;
  }
  return 0;
}


int
LinearSeries::recvSelf(int commitTag, Channel &theChannel,
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(2);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "LinearSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    tStart = 0.0;
    return result;
  }
  cFactor = data(0);
  tStart = data(1);

  return 0;
}


void
LinearSeries::Print(OPS_Stream &s, int flag)
{
  s << "Linear Series: constant factor: " << cFactor << ", start time: " << tStart << "\n";

}
