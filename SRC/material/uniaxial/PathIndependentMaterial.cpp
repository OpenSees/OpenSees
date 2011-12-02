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
                                                                        
// $Revision: 1.7 $
// $Date: 2007-02-12 20:37:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PathIndependentMaterial.cpp,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the class definition for 
// PathIndependentMaterial.  PathIndependentMaterial uses a UniaxialMaterial
// object to represent a path-independent uniaxial material.  Since
// it is path-independent, no state information is stored by
// PathIndependentMaterial.

#include <stdlib.h>
#include <PathIndependentMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

PathIndependentMaterial::PathIndependentMaterial(int tag, UniaxialMaterial &material)
:UniaxialMaterial(tag,MAT_TAG_PathIndependent), theMaterial(0)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "PathIndependentMaterial::PathIndependentMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

PathIndependentMaterial::PathIndependentMaterial()
:UniaxialMaterial(0,MAT_TAG_PathIndependent), theMaterial(0)
{

}

PathIndependentMaterial::~PathIndependentMaterial()
{
	if (theMaterial)
		delete theMaterial;
}

int 
PathIndependentMaterial::setTrialStrain(double strain, double strainRate)
{
    return theMaterial->setTrialStrain(strain, strainRate);
}

double 
PathIndependentMaterial::getStress(void)
{
    return theMaterial->getStress();
}


double 
PathIndependentMaterial::getTangent(void)
{
    return theMaterial->getTangent();
}

double 
PathIndependentMaterial::getDampTangent(void)
{
    return theMaterial->getDampTangent();
}

double 
PathIndependentMaterial::getInitialTangent(void)
{
    return theMaterial->getInitialTangent();
}

double 
PathIndependentMaterial::getStrain(void)
{
    return theMaterial->getStrain();
}

double 
PathIndependentMaterial::getStrainRate(void)
{
    return theMaterial->getStrainRate();
}

int 
PathIndependentMaterial::commitState(void)
{
  return 0; // commit nothing, path independent
}

int 
PathIndependentMaterial::revertToLastCommit(void)
{
    return 0;
}

int 
PathIndependentMaterial::revertToStart(void)
{
    return theMaterial->revertToStart();
}

UniaxialMaterial *
PathIndependentMaterial::getCopy(void)
{
    PathIndependentMaterial *theCopy = 
		new PathIndependentMaterial(this->getTag(), *theMaterial);
        
	return theCopy;
}


int 
PathIndependentMaterial::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;

	static ID classTags(3);

	int clTag = theMaterial->getClassTag();
	int dbTag = theMaterial->getDbTag();
	
	classTags(0) = clTag;

	if (dbTag == 0) {
		dbTag = theChannel.getDbTag();
		if (dbTag != 0)
			theMaterial->setDbTag(dbTag);
	}

	classTags(1) = dbTag;
	classTags(2) = this->getTag();

	res = theChannel.sendID(dbTag, cTag, classTags);
	if (res < 0) {
	  opserr << "PathIndependentMaterial::sendSelf -- could not send ID\n";
	  return res;
	}
    
	res = theMaterial->sendSelf(cTag, theChannel);
	if (res < 0) {
	  opserr << "PathIndependentMaterial::sendSelf -- could not send UniaxialMaterial\n";
	  return res;
	}

	return res;
}

int 
PathIndependentMaterial::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID classTags(3);

  int dbTag = this->getDbTag();

  res = theChannel.recvID(dbTag, cTag, classTags);
  if (res < 0) {
    opserr << "PathIndependentMaterial::recvSelf -- could not receive ID\n";
    return res;
  }

  this->setTag(int(classTags(2)));

  // Check if the material is null; if so, get a new one
  if (theMaterial == 0) {
    theMaterial = theBroker.getNewUniaxialMaterial(classTags(0));
    if (theMaterial == 0) {
      opserr << " PathIndependentMaterial::recvSelf -- could not get a UniaxialMaterial\n";
      return -1;
    }
  }
  // Check that the material is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theMaterial->getClassTag() != classTags(0)) {
    delete theMaterial;
    theMaterial = theBroker.getNewUniaxialMaterial(classTags(0));
    if (theMaterial == 0) {
      opserr << "PathIndependentMaterial::recvSelf -- could not get a UniaxialMaterial\n";
      exit(-1);
    }
  }
  
  // Now, receive the material
  theMaterial->setDbTag(classTags(1));
  res += theMaterial->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "PathIndependentMaterial::recvSelf -- could not receive UniaxialMaterial\n";
    return res;
  }
  
  return res;
}

void 
PathIndependentMaterial::Print(OPS_Stream &s, int flag)
{
    s << "PathIndependentMaterial tag: " << this->getTag() << endln;
    s << "\tmaterial: " << theMaterial->getTag() << endln;
}

double
PathIndependentMaterial::getStressSensitivity(int gradNumber, bool conditional)
{
  return theMaterial->getStressSensitivity(gradNumber, conditional);
}

double
PathIndependentMaterial::getStrainSensitivity(int gradNumber)
{
  return theMaterial->getStrainSensitivity(gradNumber);
}

double
PathIndependentMaterial::getInitialTangentSensitivity(int gradNumber)
{
  return theMaterial->getInitialTangentSensitivity(gradNumber);
}

double
PathIndependentMaterial::getDampTangentSensitivity(int gradNumber)
{
  return theMaterial->getDampTangentSensitivity(gradNumber);
}

double
PathIndependentMaterial::getRhoSensitivity(int gradNumber)
{
  return theMaterial->getRhoSensitivity(gradNumber);
}

int   
PathIndependentMaterial::commitSensitivity(double strainGradient,
					   int gradNumber, int numGrads)
{
  return 0; // commit nothing, path independent
}
