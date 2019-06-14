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
// DuctilityStiffnessDegradation, which models hysteretic stiffness
// degradation with the relation alpha*(u-1.0), based on
// deformation ductility.

#include <DuctilityStiffnessDegradation.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_DuctilityStiffnessDegradation(void)
{
  StiffnessDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: stiffnessDegradation Ductility tag? alpha? beta?" << endln;
    return 0;
  }

  int iData[1];
  double dData[2];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for stiffnessDegradation Ductility" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for stiffnessDegradation Ductility" << endln;
    return 0;
  }

  theDegradation = new DuctilityStiffnessDegradation(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create DuctilityStiffnessDegradation\n";
    return 0;
  }

  return theDegradation;
}

DuctilityStiffnessDegradation::DuctilityStiffnessDegradation
(int tag, double a, double b):
  StiffnessDegradation(tag,DEG_TAG_STIFF_Ductility),
  isNegative(false), alpha(a), beta(b)
{
  this->revertToStart();
  this->revertToLastCommit();
}

DuctilityStiffnessDegradation::DuctilityStiffnessDegradation():
  StiffnessDegradation(0,DEG_TAG_STIFF_Ductility),
  isNegative(false), alpha(0.0), beta(0.0), Cductility(0.0)
{
  
}

DuctilityStiffnessDegradation::~DuctilityStiffnessDegradation()
{
  
}

const char*
DuctilityStiffnessDegradation::getMeasure(void)
{
  if (!isNegative)
    return "negDuctility";
  else
    return "posDuctility";
}

int
DuctilityStiffnessDegradation::setTrialMeasure(double measure)
{
  Tductility = measure;
  
  //if (Tductility < Cductility)
  //	Tductility = Cductility;
  
  return 0;
}

double
DuctilityStiffnessDegradation::getValue(void)
{
  if (Tductility < Cductility) {
    Tductility = Cductility;
    return 1.0;
  }
  if (Tductility > beta)
    return 1.0 + alpha*(Tductility-beta);
  else
    return 1.0;
}

void
DuctilityStiffnessDegradation::setNegative(bool flag)
{
  isNegative = flag;
}

int
DuctilityStiffnessDegradation::commitState(void)
{
  Cductility = Tductility;
  
  return 0;
}
 
int
DuctilityStiffnessDegradation::revertToLastCommit(void)
{
  Tductility = Cductility;
  
  return 0;
}

int
DuctilityStiffnessDegradation::revertToStart(void)
{
  Cductility = 0.0;
  
  return 0;
}

StiffnessDegradation*
DuctilityStiffnessDegradation::getCopy(void)
{
  DuctilityStiffnessDegradation *theCopy =
    new DuctilityStiffnessDegradation (this->getTag(), alpha, beta);
  
  theCopy->Cductility = Cductility;
  
  return theCopy;
}

int
DuctilityStiffnessDegradation::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(5);
  
  data(0) = this->getTag();
  data(1) = alpha;
  data(2) = beta;
  data(3) = Cductility;
  data(4) = (isNegative) ? -1.0 : 1.0;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "DuctilityStiffnessDegradation::sendSelf() - failed to send data\n";
  
  return res;
}

int
DuctilityStiffnessDegradation::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(5);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "DuctilityStiffnessDegradation::recvSelf() - failed to receive data\n";
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
DuctilityStiffnessDegradation::Print(OPS_Stream &s, int flag)
{
  s << "DuctilityStiffnessDegradation, tag: " << this->getTag() << endln;
  s << "alpha: " << alpha << endln;
  s << "beta: " << beta << endln;
}
