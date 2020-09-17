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
// ConstantUnloadingRule, which models hysteretic unloading stiffness
// degradation with the simple relation f(i) = alpha*f(i-1) + beta.

#include <ConstantUnloadingRule.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_ConstantUnloadingRule(void)
{
  UnloadingRule *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: unloadingRule Constant tag? alpha? beta?" << endln;
    return 0;
  }

  int iData[1];
  double dData[2];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for unloadingRule Constant" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for unloadingRule Constant" << endln;
    return 0;
  }

  theDegradation = new ConstantUnloadingRule(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create ConstantUnloadingRule\n";
    return 0;
  }

  return theDegradation;
}

ConstantUnloadingRule::ConstantUnloadingRule(int tag, double a, double b):
  UnloadingRule(tag,DEG_TAG_UNLOAD_Constant), alpha(a), beta(b)
{
  this->revertToStart();
  this->revertToLastCommit();
}

ConstantUnloadingRule::ConstantUnloadingRule():
  UnloadingRule(0,DEG_TAG_UNLOAD_Constant), 
  alpha(0.0), beta(0.0), Cfactor(0.0)
{

}

ConstantUnloadingRule::~ConstantUnloadingRule()
{

}

const char*
ConstantUnloadingRule::getMeasure(void)
{
  return "nothing";
}

int
ConstantUnloadingRule::setTrialMeasure(double measure)
{
  return 0;
}

double
ConstantUnloadingRule::getValue(void)
{
  Tfactor = alpha*Cfactor + beta;
  
  return Tfactor;
}

int
ConstantUnloadingRule::commitState(void)
{
  Cfactor = Tfactor;
  
  return 0;
}
 
int
ConstantUnloadingRule::revertToLastCommit(void)
{
  Tfactor = Cfactor;
  
  return 0;
}

int
ConstantUnloadingRule::revertToStart(void)
{
  Cfactor = 1.0;
  
  return 0;
}

UnloadingRule*
ConstantUnloadingRule::getCopy(void)
{
  ConstantUnloadingRule *theCopy =
    new ConstantUnloadingRule (this->getTag(), alpha, beta);
  
  theCopy->Cfactor = Cfactor;
  
  return theCopy;
}

int
ConstantUnloadingRule::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(4);
  
  data(0) = this->getTag();
  data(1) = alpha;
  data(2) = beta;
  data(3) = Cfactor;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "ConstantUnloadingRule::sendSelf() - failed to send data\n";
  
  return res;
}

int
ConstantUnloadingRule::recvSelf(int commitTag, Channel &theChannel, 
				FEM_ObjectBroker &theBroker)
{
  static Vector data(4);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "ConstantUnloadingRule::recvSelf() - failed to receive data\n";
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
ConstantUnloadingRule::Print(OPS_Stream &s, int flag)
{
  s << "ConstantUnloadingRule, tag: " << this->getTag() << endln;
  s << "\talpha: " << alpha << endln;
  s << "\tbeta: " << beta << endln;
}
