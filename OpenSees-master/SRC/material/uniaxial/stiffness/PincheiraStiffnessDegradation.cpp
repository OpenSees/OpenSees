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
// Created: March 2001
//
// Description: This file contains the implementation of 
// PincheiraStiffnessDegradation, which models hysteretic stiffness
// degradation with the relation alpha*(u-1.0), based on
// deformation ductility.

#include <PincheiraStiffnessDegradation.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_PincheiraStiffnessDegradation(void)
{
  StiffnessDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: stiffnessDegradation Pincheira tag? alpha? beta? eta? nu?" << endln;
    return 0;
  }

  int iData[1];
  double dData[4];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for stiffnessDegradation Pincheira" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for stiffnessDegradation Pincheira" << endln;
    return 0;
  }

  theDegradation = new PincheiraStiffnessDegradation(iData[0], dData[0], dData[1], dData[2], dData[3]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create PincheiraStiffnessDegradation\n";
    return 0;
  }

  return theDegradation;
}

PincheiraStiffnessDegradation::PincheiraStiffnessDegradation
(int tag, double a, double b, double e, double n):
  StiffnessDegradation(tag,DEG_TAG_STIFF_Pincheira),
  isNegative(false), alpha(a), beta(b), eta(e), nu(n)
{
  this->revertToStart();
  this->revertToLastCommit();
}
 
PincheiraStiffnessDegradation::PincheiraStiffnessDegradation():
  StiffnessDegradation(0,DEG_TAG_STIFF_Pincheira),
  isNegative(false), alpha(0.0), beta(0.0), eta(0.0), nu(0.0)
{

}

PincheiraStiffnessDegradation::~PincheiraStiffnessDegradation()
{

}

const char*
PincheiraStiffnessDegradation::getMeasure(void)
{
  if (!isNegative)
    return "negDuctility";
  else
    return "posDuctility";
}

int
PincheiraStiffnessDegradation::setTrialMeasure(double measure)
{
  Tductility = measure;
  
  return 0;
}

double
PincheiraStiffnessDegradation::getValue(void)
{
  TmaxDuctility = CmaxDuctility;
  TnumCycles = CnumCycles;
  TcycleFlag = CcycleFlag;

  if (Tductility > CmaxDuctility) {
    TmaxDuctility = Tductility;
    TcycleFlag = false;
    TnumCycles = 0;
    return 1.0 + beta*(Tductility-alpha);
  }
  else if (Tductility > alpha) {
    TcycleFlag = true;
    TnumCycles++;
    return 1.0 + eta*pow(nu,TnumCycles-1)*(Tductility-alpha);
    //TnumCycles = 1;
    //return 1.0;
  }
  else
    return 1.0;
}

void
PincheiraStiffnessDegradation::setNegative(bool flag)
{
  isNegative = flag;
}

int
PincheiraStiffnessDegradation::commitState(void)
{
  CmaxDuctility = TmaxDuctility;
  CnumCycles = TnumCycles;
  CcycleFlag = TcycleFlag;

  //opserr << Tductility << ' ' << CmaxDuctility << ' ' << CnumCycles << endln;

  return 0;
}
 
int
PincheiraStiffnessDegradation::revertToLastCommit(void)
{
  TmaxDuctility = CmaxDuctility;
  TnumCycles = CnumCycles;
  TcycleFlag = CcycleFlag;

  return 0;
}

int
PincheiraStiffnessDegradation::revertToStart(void)
{
  CmaxDuctility = alpha;
  CnumCycles = 0;
  CcycleFlag = false;
  
  return 0;
}

StiffnessDegradation*
PincheiraStiffnessDegradation::getCopy(void)
{
  PincheiraStiffnessDegradation *theCopy =
    new PincheiraStiffnessDegradation (this->getTag(), alpha, beta, eta, nu);

  theCopy->CmaxDuctility = CmaxDuctility;
  theCopy->CnumCycles = CnumCycles;
  theCopy->CcycleFlag = CcycleFlag;
  
  return theCopy;
}

int
PincheiraStiffnessDegradation::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(7);
  
  data(0) = this->getTag();
  data(1) = alpha;
  data(2) = beta;
  data(3) = eta;
  data(4) = nu;
  data(5) = CmaxDuctility;
  data(6) = (isNegative) ? -1.0 : 1.0;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "PincheiraStiffnessDegradation::sendSelf() - failed to send data\n";
  
  return res;
}

int
PincheiraStiffnessDegradation::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(7);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "DuctilityStiffnessDegradation::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    alpha = data(1);
    beta = data(2);
    eta = data(3);
    nu = data(4);
    CmaxDuctility = data(5);
    isNegative = (data(6) < 0.0) ? true : false;
  }
  
  return res;
}

void
PincheiraStiffnessDegradation::Print(OPS_Stream &s, int flag)
{
  s << "PincheiraStiffnessDegradation, tag: " << this->getTag() << endln;
  s << "\talpha: " << alpha << endln;
  s << "\tbeta: " << beta << endln;
  s << "\teta: " << eta << endln;
  s << "\tnu: " << nu << endln;
}
