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
                                                                        
// Written: MHS
// Created: March 2001
//
// Description: This file contains the implementation of 
// KarsanUnloadingRule, which models concrete unloading stiffness
// degradation.

#include <KarsanUnloadingRule.h>
#include <G3Globals.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_KarsanUnloadingRule(void)
{
  UnloadingRule *theDegradation = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid number of args, want: unloadingRule Karsan tag? epsc? epsu?" << endln;
    return 0;
  }

  int iData[1];
  double dData[2];
  
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for unloadingRule Karsan" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for unloadingRule Karsan" << endln;
    return 0;
  }

  theDegradation = new KarsanUnloadingRule(iData[0], dData[0], dData[1]);
  if (theDegradation == 0) {
    opserr << "WARNING could not create KarsanUnloadingRule\n";
    return 0;
  }

  return theDegradation;
}

KarsanUnloadingRule::KarsanUnloadingRule
(int tag, double e, double eu):
  UnloadingRule(tag,DEG_TAG_UNLOAD_Karsan),
  epsc(e), epscu(eu)
{
  this->revertToStart();
  this->revertToLastCommit();
}
 
KarsanUnloadingRule::KarsanUnloadingRule():
  UnloadingRule(0,DEG_TAG_UNLOAD_Karsan), 
  epsc(0.0), epscu(0.0), CminStrain(0.0)
{

}

KarsanUnloadingRule::~KarsanUnloadingRule()
{

}

const char*
KarsanUnloadingRule::getMeasure(void)
{
  return "minStrain";
}

int
KarsanUnloadingRule::setTrialMeasure(double measure)
{
  TminStrain = measure;
  
  if (TminStrain > CminStrain)
    TminStrain = CminStrain;
  
  if (TminStrain < epscu)
    TminStrain = epscu;
  
  return 0;
}

double
KarsanUnloadingRule::getValue(void)
{
  double eta = TminStrain/epsc;
  
  double ratio = 0.707*(eta-2.0) + 0.834;
  
  if (eta < 2.0)
    ratio = 0.145*eta*eta + 0.13*eta;
  
  return ratio*epsc;
}

int
KarsanUnloadingRule::commitState(void)
{
  CminStrain = TminStrain;
  
  return 0;
}
 
int
KarsanUnloadingRule::revertToLastCommit(void)
{
  TminStrain = CminStrain;
  
  return 0;
}

int
KarsanUnloadingRule::revertToStart(void)
{
  CminStrain = 0.0;

  return 0;
}

UnloadingRule*
KarsanUnloadingRule::getCopy(void)
{
  KarsanUnloadingRule *theCopy =
    new KarsanUnloadingRule (this->getTag(), epsc, epscu);
  
  theCopy->CminStrain = CminStrain;
  
  return theCopy;
}

int
KarsanUnloadingRule::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
KarsanUnloadingRule::recvSelf(int commitTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
KarsanUnloadingRule::Print(OPS_Stream &s, int flag)
{
  s << "KarsanUnloadingRule, tag: " << this->getTag() << endln;
}
