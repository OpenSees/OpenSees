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
// $Source$

// Written: MHS
// Created: August 2000
//
// Description: This file contains the implementation of 
// ConstantStrengthDegradation, which models hysteretic strength
// degradation with the simple relation f(i) = alpha*f(i-1) + beta.

#include <ConstantStrengthDegradation.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_ConstantStrengthDegradation(void)
{
  StrengthDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: strengthDegradation Constant tag? alpha? beta?" << endln;
    return 0;
  }

  int iData[1];
  double dData[2];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for strengthDegradation Constant" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for strengthDegradation Constant" << endln;
    return 0;
  }

  theDegradation = new ConstantStrengthDegradation(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create ConstantStrengthDegradation\n";
    return 0;
  }

  return theDegradation;
}

ConstantStrengthDegradation::ConstantStrengthDegradation(int tag, double a, double b):
  StrengthDegradation(tag,DEG_TAG_STRENGTH_Constant), alpha(a), beta(b)
{
  this->revertToStart();
  this->revertToLastCommit();
}

ConstantStrengthDegradation::ConstantStrengthDegradation():
  StrengthDegradation(0,DEG_TAG_STRENGTH_Constant), alpha(0.0), beta(0.0), Cfactor(0.0)
{
  
}

ConstantStrengthDegradation::~ConstantStrengthDegradation()
{

}

const char*
ConstantStrengthDegradation::getMeasure(void)
{
  return "nothing";
}

int
ConstantStrengthDegradation::setTrialMeasure(double measure)
{
  Tfactor = alpha*Cfactor;
  
  return 0;
}

double
ConstantStrengthDegradation::getValue(void)
{
  return Tfactor + beta;
}

int
ConstantStrengthDegradation::commitState(void)
{
  Cfactor = Tfactor;

  return 0;
}
 
int
ConstantStrengthDegradation::revertToLastCommit(void)
{
  Tfactor = Cfactor;
  
  return 0;
}

int
ConstantStrengthDegradation::revertToStart(void)
{
  Cfactor = 1.0;
  
  return 0;
}

StrengthDegradation*
ConstantStrengthDegradation::getCopy(void)
{
  ConstantStrengthDegradation *theCopy =
    new ConstantStrengthDegradation (this->getTag(), alpha, beta);
  
  theCopy->Cfactor = Cfactor;
  
  return theCopy;
}

int
ConstantStrengthDegradation::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(4);
  
  data(0) = this->getTag();
  data(1) = alpha;
  data(2) = beta;
  data(3) = Cfactor;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "ConstantStrengthDegradation::sendSelf() - failed to send data\n";
  
  return res;
}

int
ConstantStrengthDegradation::recvSelf(int commitTag, Channel &theChannel, 
				      FEM_ObjectBroker &theBroker)
{
  static Vector data(4);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "ConstantStrengthDegradation::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    alpha = data(1);
    beta = data(2);
    Cfactor = data(3);
  }
  
  return res;
}

void
ConstantStrengthDegradation::Print(OPS_Stream &s, int flag)
{
  s << "ConstantStrengthDegradation, tag: " << this->getTag() << endln;
  s << "\talpha: " << alpha << endln;
  s << "\tbeta: " << beta << endln;
}
