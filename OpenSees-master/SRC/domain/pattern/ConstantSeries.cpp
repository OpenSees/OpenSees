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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/ConstantSeries.cpp,v $

// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for ConstantSeries.
// ConstantSeries is a concrete class. A ConstantSeries object provides
// a linear time series. the factor is given by the pseudoTime and 
// a constant factor provided in the constructor.
//
// What: "@(#) ConstantSeries.C, revA"


#include <ConstantSeries.h>
#include <Vector.h>
#include <Channel.h>
#include <classTags.h>
#include <Parameter.h>

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_ConstantSeries(void)
{
  // Pointer to a uniaxial material that will be returned
  TimeSeries *theSeries = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  int tag = 0;
  double cFactor = 1.0;
  int numData = 0;

  if (numRemainingArgs != 0) {
  
    if (numRemainingArgs == 1 || numRemainingArgs == 3) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &tag) != 0) {
	opserr << "WARNING invalid series tag in ConstantSeries tag? <-factor factor?>" << endln;
	return 0;
      }
      numRemainingArgs --;
    }

    if (numRemainingArgs > 1) {
      const char *argvS = OPS_GetString();
	  if (argvS == 0) {
		  opserr << "WARNING string error in  ConstantSeries with tag: " << tag << endln;
		return 0;
	  }
      numData = 1;
      if (OPS_GetDouble(&numData, &cFactor) != 0) {
	  opserr << "WARNING invalid factor in  ConstantSeries with tag: " << tag << endln;
	  return 0;
      }
    }
  }

  theSeries = new ConstantSeries(tag, cFactor);

  if (theSeries == 0) {
    opserr << "WARNING ran out of memory creating ConstantTimeSeries with tag: " << tag << "\n";
    return 0;
  }

  return theSeries;
}


ConstantSeries::ConstantSeries(int tag, double theFactor)
  :TimeSeries(tag,  TSERIES_TAG_ConstantSeries),
   cFactor(theFactor), parameterID(0)
{
  // does nothing
}


TimeSeries *
ConstantSeries::getCopy(void) {
  return new ConstantSeries(this->getTag(), cFactor);
}

ConstantSeries::~ConstantSeries()
{
  // does nothing
}



int
ConstantSeries::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  Vector data(1);
  data(0) = cFactor;
  int result = theChannel.sendVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "ConstantSeries::sendSelf() - channel failed to send data\n";
    return result;
  }
  return 0;
}


int 
ConstantSeries::recvSelf(int commitTag, Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  Vector data(1);
  int result = theChannel.recvVector(dbTag,commitTag, data);
  if (result < 0) {
    opserr << "ConstantSeries::sendSelf() - channel failed to receive data\n";
    cFactor = 1.0;
    return result;
  }
  cFactor = data(0);

  return 0;    
}


void
ConstantSeries::Print(OPS_Stream &s, int flag)
{
    s << "Constant Series: factor: " << cFactor << "\n";
}

    // AddingSensitivity:BEGIN //////////////////////////////////////////
double
ConstantSeries::getFactorSensitivity(double pseudoTime)
{
  if (parameterID == 1)
    return 1.0;
  else
    return 0.0;
}

int 
ConstantSeries::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strstr(argv[0],"factor") != 0) {
    param.setValue(cFactor);
    return param.addObject(1, this);
  }

  return -1;
}
   
int 
ConstantSeries::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    cFactor = info.theDouble;
    return 0;
  }

  return -1;
}

int
ConstantSeries::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}
