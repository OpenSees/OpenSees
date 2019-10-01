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
// TakedaUnloadingRule, which models hysteretic unloading stiffness
// degradation with the Takeda type relation C*u^beta, based on
// deformation ductility.

#include <TakedaUnloadingRule.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_TakedaUnloadingRule(void)
{
  UnloadingRule *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: unloadingRule Takeda(Ductility) tag? alpha? beta?" << endln;
    return 0;
  }

  int iData[1];
  double dData[2];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for unloadingRule Takeda(Ductility)" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for unloadingRule Takeda(Ductility)" << endln;
    return 0;
  }

  theDegradation = new TakedaUnloadingRule(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create TakedaUnloadingRule\n";
    return 0;
  }

  return theDegradation;
}

TakedaUnloadingRule::TakedaUnloadingRule
(int tag, double c, double b):
  UnloadingRule(tag,DEG_TAG_UNLOAD_Takeda),
  isNegative(false), C(c), beta(b)
{
  this->revertToStart();
  this->revertToLastCommit();
}
 
TakedaUnloadingRule::TakedaUnloadingRule():
  UnloadingRule(0,DEG_TAG_UNLOAD_Takeda), 
  isNegative(false), C(0.0), beta(0.0), Cductility(0.0)
{

}

TakedaUnloadingRule::~TakedaUnloadingRule()
{

}

const char*
TakedaUnloadingRule::getMeasure(void)
{
  if (isNegative)
    return "negDuctility";
  else
    return "posDuctility";
}

int
TakedaUnloadingRule::setTrialMeasure(double measure)
{
  Tductility = measure;
  
  if (Tductility < Cductility)
    Tductility = Cductility;
  
  return 0;
}

double
TakedaUnloadingRule::getValue(void)
{
  if (Tductility < 1.0)
    return 1.0;
  else
    return C*pow(1.0/Tductility,beta);
}

void
TakedaUnloadingRule::setNegative(bool flag)
{
  isNegative = flag;
}

int
TakedaUnloadingRule::commitState(void)
{
  Cductility = Tductility;
  
  return 0;
}
 
int
TakedaUnloadingRule::revertToLastCommit(void)
{
  Tductility = Cductility;
  
  return 0;
}

int
TakedaUnloadingRule::revertToStart(void)
{
  Cductility = 0.0;

  return 0;
}

UnloadingRule*
TakedaUnloadingRule::getCopy(void)
{
  TakedaUnloadingRule *theCopy =
    new TakedaUnloadingRule (this->getTag(), C, beta);
  
  theCopy->Cductility = Cductility;
  
  return theCopy;
}

int
TakedaUnloadingRule::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(5);
  
  data(0) = this->getTag();
  data(1) = C;
  data(2) = beta;
  data(3) = Cductility;
  data(4) = (isNegative) ? -1.0 : 1.0;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "TakedaUnloadingRule::sendSelf() - failed to send data\n";
  
  return res;
}

int
TakedaUnloadingRule::recvSelf(int commitTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  static Vector data(5);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "TakedaUnloadingRule::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    C = data(1);
    beta = data(2);
    Cductility = data(3);
    isNegative = (data(4) < 0.0) ? true : false;
  }
  
  return res;
}

void
TakedaUnloadingRule::Print(OPS_Stream &s, int flag)
{
  s << "TakedaUnloadingRule, tag: " << this->getTag() << endln;
  s << "C: " << C << endln;
  s << "beta: " << beta << endln;
}
