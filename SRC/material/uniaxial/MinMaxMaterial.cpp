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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-04-02 22:02:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/MinMaxMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition for 
// MinMaxMaterial.  MinMaxMaterial wraps a UniaxialMaterial
// and imposes min and max strain limits.

#include <stdlib.h>

#include <MinMaxMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

MinMaxMaterial::MinMaxMaterial(int tag, UniaxialMaterial &material,
			       double min, double max)
  :UniaxialMaterial(tag,MAT_TAG_MinMax), theMaterial(0),
   minStrain(min), maxStrain(max), Tfailed(false), Cfailed(false)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "MinMaxMaterial::MinMaxMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

MinMaxMaterial::MinMaxMaterial()
  :UniaxialMaterial(0,MAT_TAG_MinMax), theMaterial(0),
   minStrain(0.0), maxStrain(0.0), Tfailed(false), Cfailed(false)
{

}

MinMaxMaterial::~MinMaxMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
MinMaxMaterial::setTrialStrain(double strain, double strainRate)
{
  if (Cfailed)
    return 0;
  
  if (strain >= maxStrain || strain <= minStrain) {
    Tfailed = true;
    return 0;
  }
  else {
    Tfailed = false;
    return theMaterial->setTrialStrain(strain, strainRate);
  }
}

double 
MinMaxMaterial::getStress(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getStress();
}

double 
MinMaxMaterial::getTangent(void)
{
  if (Tfailed)
    //return 0.0;
    return 1.0e-8*theMaterial->getInitialTangent();
  else
    return theMaterial->getTangent();
}

double 
MinMaxMaterial::getDampTangent(void)
{
  if (Tfailed)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}



double 
MinMaxMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
MinMaxMaterial::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
MinMaxMaterial::commitState(void)
{	
  Cfailed = Tfailed;

  // Check if failed at current step
  if (Tfailed)
    return 0;
  else
    return theMaterial->commitState();
}

int 
MinMaxMaterial::revertToLastCommit(void)
{
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

int 
MinMaxMaterial::revertToStart(void)
{
  Cfailed = false;
  Tfailed = false;
  
  return theMaterial->revertToStart();
}

UniaxialMaterial *
MinMaxMaterial::getCopy(void)
{
  MinMaxMaterial *theCopy = 
    new MinMaxMaterial(this->getTag(), *theMaterial, minStrain, maxStrain);
        
  theCopy->Cfailed = Cfailed;
  theCopy->Tfailed = Tfailed;
  
  return theCopy;
}

int 
MinMaxMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if ( matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "MinMaxMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(3);
  dataVec(0) = minStrain;
  dataVec(1) = maxStrain;
  if (Cfailed == true)
    dataVec(2) = 1.0;
  else
    dataVec(2) = 0.0;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "MinMaxMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "MinMaxMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
MinMaxMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "MinMaxMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "MinMaxMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(3);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "MinMaxMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  minStrain = dataVec(0);
  maxStrain = dataVec(1);
  
  if (dataVec(2) == 1.0)
    Cfailed = true;
  else
    Cfailed = false;

  Tfailed = Cfailed;

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "MinMaxMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
MinMaxMaterial::Print(OPS_Stream &s, int flag)
{
  s << "MinMaxMaterial tag: " << this->getTag() << endln;
  s << "\tMaterial: " << theMaterial->getTag() << endln;
  s << "\tMin strain: " << minStrain << endln;
  s << "\tMax strain: " << maxStrain << endln;
}
