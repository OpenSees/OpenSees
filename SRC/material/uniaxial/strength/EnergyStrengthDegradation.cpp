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
// EnergyStrengthDegradation, which models hysteretic strength
// degradation with the deterioration parameter developed by
// Rahnama and Krawinkler (1993).

#include <EnergyStrengthDegradation.h>
#include <UniaxialMaterial.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_EnergyStrengthDegradation(void)
{
  StrengthDegradation *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: strengthDegradation Energy tag? Et? c?" << endln;
    return 0;
  }

  int iData[1];
  double dData[2];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for strengthDegradation Energy" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for strengthDegradation Energy" << endln;
    return 0;
  }

  theDegradation = new EnergyStrengthDegradation(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create EnergyStrengthDegradation\n";
    return 0;
  }

  return theDegradation;
}

EnergyStrengthDegradation::EnergyStrengthDegradation(int tag, double et, double C):
  StrengthDegradation(tag,DEG_TAG_STRENGTH_Energy),
  Et(et), c(C), energyExcursion(0.0)
{
  this->revertToStart();
  this->revertToLastCommit();
}

EnergyStrengthDegradation::EnergyStrengthDegradation():
  StrengthDegradation(0,DEG_TAG_STRENGTH_Energy),
  Et(0.0), c(0.0), energyExcursion(0.0), Cfactor(0.0)
{

}

EnergyStrengthDegradation::~EnergyStrengthDegradation()
{
  
}

const char*
EnergyStrengthDegradation::getMeasure(void)
{
  return "energyExcursion";
}

int
EnergyStrengthDegradation::setTrialMeasure(double measure)
{
  energyExcursion = measure;
  TenergySum = CenergySum + energyExcursion;
  
  return 0;
}

double
EnergyStrengthDegradation::getValue(void)
{
  if (TenergySum >= Et)
    return Cfactor;
  else {
    // Beta = (E_i/(E_t-sum(E_i))^c
    double beta = pow(energyExcursion/(Et-TenergySum),c);
    if (beta > 1.0) {
      opserr << "Beta: " << beta << endln;
      beta = 1.0;
    }
    Tfactor = (1.0-beta)*Cfactor;
    return Tfactor;
  }
}

int
EnergyStrengthDegradation::commitState(void)
{
  Cfactor = Tfactor;
  CenergySum = TenergySum;
  
  return 0;
}
 
int
EnergyStrengthDegradation::revertToLastCommit(void)
{
  Tfactor = Cfactor;
  TenergySum = CenergySum;
  
  return 0;
}

int
EnergyStrengthDegradation::revertToStart(void)
{
  Cfactor = 1.0;
  CenergySum = 0.0;
  
  return 0;
}

StrengthDegradation*
EnergyStrengthDegradation::getCopy(void)
{
  EnergyStrengthDegradation *theCopy =
    new EnergyStrengthDegradation (this->getTag(), Et, c);
  
  theCopy->Cfactor = Cfactor;
  theCopy->CenergySum = CenergySum;
  
  return theCopy;
}

int
EnergyStrengthDegradation::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(5);
  
  data(0) = this->getTag();
  data(1) = Et;
  data(2) = c;
  data(3) = Cfactor;
  data(4) = CenergySum;
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) 
    opserr << "EnergyStrengthDegradation::sendSelf() - failed to send data\n";
  
  return res;
}

int
EnergyStrengthDegradation::recvSelf(int commitTag, Channel &theChannel, 
				    FEM_ObjectBroker &theBroker)
{
  static Vector data(5);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
  if (res < 0) {
    opserr << "EnergyStrengthDegradation::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    Et = data(1);
    c = data(2);
    Cfactor = data(3);
    CenergySum = data(4);
  }
  
  return res;
}

void
EnergyStrengthDegradation::Print(OPS_Stream &s, int flag)
{
  s << "EnergyStrengthDegradation, tag: " << this->getTag() << endln;
  s << "\tEt: " << Et << endln;
  s << "\tc: " << c << endln;
}
