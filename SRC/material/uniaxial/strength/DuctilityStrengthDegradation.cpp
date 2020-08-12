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
// DuctilityStrengthDegradation, which models hysteretic Strength
// degradation with the relation alpha*(u-1.0), based on
// deformation ductility.

#include <DuctilityStrengthDegradation.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_DuctilityStrengthDegradation(void)
{
  StrengthDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: strengthDegradation Ductility tag? alpha? beta?" << endln;
    return 0;
  }

  int iData[1];
  double dData[2];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for strengthDegradation Ductility" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for strengthDegradation Ductility" << endln;
    return 0;
  }

  theDegradation = new DuctilityStrengthDegradation(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create DuctilityStrengthDegradation\n";
    return 0;
  }

  return theDegradation;
}

DuctilityStrengthDegradation::DuctilityStrengthDegradation
(int tag, double a, double b):
  StrengthDegradation(tag,DEG_TAG_STRENGTH_Ductility),
  isNegative(false), alpha(a), beta(b)
{
  this->revertToStart();
  this->revertToLastCommit();
}
 
DuctilityStrengthDegradation::DuctilityStrengthDegradation():
  StrengthDegradation(0,DEG_TAG_STRENGTH_Ductility),
  isNegative(false), alpha(0.0), beta(0.0), Cductility(0.0)
{

}

DuctilityStrengthDegradation::~DuctilityStrengthDegradation()
{

}

const char*
DuctilityStrengthDegradation::getMeasure(void)
{
  if (!isNegative)
    return "negDuctility";
  else
    return "posDuctility";
}

int
DuctilityStrengthDegradation::setTrialMeasure(double measure)
{
  Tductility = measure;
  
  //if (Tductility < Cductility)
  //	Tductility = Cductility;
  
  return 0;
}

double
DuctilityStrengthDegradation::getValue(void)
{
  if (Tductility < Cductility) {
    Tductility = Cductility;
    return 1.0;
  }
  else if (Tductility > beta)
    return 1.0 - alpha*(Tductility-beta);
  else
    return 1.0;
}

void
DuctilityStrengthDegradation::setNegative(bool flag)
{
  isNegative = flag;
}

int
DuctilityStrengthDegradation::commitState(void)
{
  Cductility = Tductility;
  
  return 0;
}
 
int
DuctilityStrengthDegradation::revertToLastCommit(void)
{
  Tductility = Cductility;
  
  return 0;
}

int
DuctilityStrengthDegradation::revertToStart(void)
{
  Cductility = 0.0;
  
  return 0;
}

StrengthDegradation*
DuctilityStrengthDegradation::getCopy(void)
{
  DuctilityStrengthDegradation *theCopy =
    new DuctilityStrengthDegradation (this->getTag(), alpha, beta);
  
  theCopy->Cductility = Cductility;
  
  return theCopy;
}

int
DuctilityStrengthDegradation::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(5);
  
  data(0) = this->getTag();
  data(1) = alpha;
  data(2) = beta;
  data(3) = Cductility;
  data(4) = (isNegative) ? -1.0 : 1.0;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "DuctilityStrengthDegradation::sendSelf() - failed to send data\n";
  
  return res;
}

int
DuctilityStrengthDegradation::recvSelf(int commitTag, Channel &theChannel, 
				       FEM_ObjectBroker &theBroker)
{
  static Vector data(5);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "DuctilityStrengthDegradation::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    alpha = data(1);
    beta = data(2);
    Cductility = data(3);
    isNegative = (data(4) < 0.0) ? true : false;
  }
  
  return res;
}

void
DuctilityStrengthDegradation::Print(OPS_Stream &s, int flag)
{
  s << "DuctilityStrengthDegradation, tag: " << this->getTag() << endln;
  s << "alpha: " << alpha << endln;
  s << "beta: " << beta << endln;
}
