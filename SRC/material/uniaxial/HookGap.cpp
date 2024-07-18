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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HookGap.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// HookGap. 
//
// What: "@(#) HookGap.C, revA"

#include <HookGap.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <OPS_Globals.h>

#include <elementAPI.h>

void *
OPS_HookGap(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid #args,  want: uniaxialMaterial Elastic tag? E? gap? ... " << endln;
    return 0;
  }
  
  int iData[1];
  double dData[3];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial HookGapMaterial" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData >= 3) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxialMaterial HookGap " << iData[0] << endln;
      return 0;	
    }
  } else {
    numData = 2;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid data for uniaxialMaterial HookGap " << iData[0] << endln;
      return 0;	
    }
    dData[2] = dData[1];
    dData[1] = -dData[1];
  }

  // Parsing was successful, allocate the material
  theMaterial = new HookGap(iData[0], dData[0], dData[1], dData[2]);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type HookGap\n";
    return 0;
  }

  return theMaterial;
}


HookGap::HookGap(int tag, double e, double gap)
:UniaxialMaterial(tag,MAT_TAG_HookGap),
 trialStrain(0.0),  
 E(e), gapN(-gap), gapP(gap)
{

}

HookGap::HookGap(int tag, double e, double gapNeg, double gapPos)
:UniaxialMaterial(tag,MAT_TAG_HookGap),
 trialStrain(0.0),  
 E(e), gapN(gapNeg), gapP(gapPos)
{
  // Make sure gapN is negative
  if (gapN > 0.0)
    gapN = -gapN;
}

HookGap::HookGap()
:UniaxialMaterial(0,MAT_TAG_HookGap),
 trialStrain(0.0),  
 E(0.0), gapN(0.0), gapP(0)
{

}

HookGap::~HookGap()
{
  // does nothing
}

int 
HookGap::setTrialStrain(double strain, double strainRate)
{
    trialStrain     = strain;
    return 0;
}


int 
HookGap::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
    trialStrain     = strain;
    return 0;
}

double 
HookGap::getStress(void)
{
  if (gapN <= trialStrain && trialStrain <= gapP)
    return 0;
  else if (trialStrain > gapP) {
    return E*(trialStrain-gapP);
  } else {
    return E*(trialStrain-gapN);
  }
}

double 
HookGap::getTangent(void)
{
  if (gapN < trialStrain && trialStrain < gapP)
    return 0;
  else 
    return E;
}

double 
HookGap::getInitialTangent(void)
{
  return 0;
}


int 
HookGap::commitState(void)
{
    return 0;
}


int 
HookGap::revertToLastCommit(void)
{
    return 0;
}


int 
HookGap::revertToStart(void)
{
    trialStrain      = 0.0;
    return 0;
}

UniaxialMaterial *
HookGap::getCopy(void)
{
  HookGap *theCopy = new HookGap(this->getTag(), E, gapN, gapP);
    theCopy->trialStrain     = trialStrain;
    return theCopy;
}

int 
HookGap::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(4);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = gapN;
  data(3) = gapP;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "HookGap::sendSelf() - failed to send data\n";

  return res;
}

int 
HookGap::recvSelf(int cTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(4);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "HookGap::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag(int(data(0)));
    E   = data(1);
    gapN = data(2);
    gapP = data(3);
  }
    
  return res;
}

void 
HookGap::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "HookGap tag: " << this->getTag() << endln;
        s << "  E: " << E << " gapN: " << gapN << " gapP: " << gapP << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"HookGap\", ";
        s << "\"E\": " << E << ", ";
        s << "\"gapN\": " << gapN << ", ";
        s << "\"gapP\": " << gapP << "}";
    }
}

int
HookGap::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);
  
  else if (strcmp(argv[0],"eta") == 0)
    return param.addObject(2, this);

  return -1;
}

int 
HookGap::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    gapP = info.theDouble;
    gapN = -gapP;
    return 0;
  default:
    return -1;
  }
}

