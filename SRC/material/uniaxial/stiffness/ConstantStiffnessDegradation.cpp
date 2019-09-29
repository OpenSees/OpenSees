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
// ConstantStiffnessDegradation, which models hysteretic stiffness
// degradation with the simple relation f(i) = alpha*f(i-1) + beta.

#include <ConstantStiffnessDegradation.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_ConstantStiffnessDegradation(void)
{
  StiffnessDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: stiffnessDegradation Constant tag? alpha? beta?" << endln;
    return 0;
  }

  int iData[1];
  double dData[2];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for stiffnessDegradation Constant" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for stiffnessDegradation Constant" << endln;
    return 0;
  }

  theDegradation = new ConstantStiffnessDegradation(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create ConstantStiffnessDegradation\n";
    return 0;
  }

  return theDegradation;
}

ConstantStiffnessDegradation::ConstantStiffnessDegradation(int tag, double a, double b):
  StiffnessDegradation(tag,DEG_TAG_STIFF_Constant), alpha(a), beta(b)
{
  this->revertToStart();
  this->revertToLastCommit();
}

ConstantStiffnessDegradation::ConstantStiffnessDegradation():
  StiffnessDegradation(0,DEG_TAG_STIFF_Constant), alpha(0.0), beta(0.0), Cfactor(0.0)
{
  
}

ConstantStiffnessDegradation::~ConstantStiffnessDegradation()
{
  
}

const char*
ConstantStiffnessDegradation::getMeasure(void)
{
  return "nothing";
}

int
ConstantStiffnessDegradation::setTrialMeasure(double measure)
{
  return 0;
}

double
ConstantStiffnessDegradation::getValue(void)
{
  Tfactor = alpha*Cfactor;
  
  return Tfactor + beta;
}

int
ConstantStiffnessDegradation::commitState(void)
{
  Cfactor = Tfactor;
  
  return 0;
}
 
int
ConstantStiffnessDegradation::revertToLastCommit(void)
{
  Tfactor = Cfactor;
  
  return 0;
}

int
ConstantStiffnessDegradation::revertToStart(void)
{
  Cfactor = 1.0;

  return 0;
}

StiffnessDegradation*
ConstantStiffnessDegradation::getCopy(void)
{
  ConstantStiffnessDegradation *theCopy =
    new ConstantStiffnessDegradation (this->getTag(), alpha, beta);
  
  theCopy->Cfactor = Cfactor;
  
  return theCopy;
}

int
ConstantStiffnessDegradation::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(4);
  
  data(0) = this->getTag();
  data(1) = alpha;
  data(2) = beta;
  data(3) = Cfactor;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "ConstantStiffnessDegradation::sendSelf() - failed to send data\n";
  
  return res;
}

int
ConstantStiffnessDegradation::recvSelf(int commitTag, Channel &theChannel, 
				       FEM_ObjectBroker &theBroker)
{
  static Vector data(4);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "ConstantStiffnessDegradation::recvSelf() - failed to receive data\n";
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
ConstantStiffnessDegradation::Print(OPS_Stream &s, int flag)
{
  s << "ConstantStiffnessDegradation, tag: " << this->getTag() << endln;
  s << "\talpha: " << alpha << endln;
  s << "\tbeta: " << beta << endln;
}
