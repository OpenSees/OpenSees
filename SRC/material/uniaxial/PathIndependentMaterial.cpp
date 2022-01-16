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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-08-26 16:33:12 $
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
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_PathIndependentMaterial)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 2) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial PathIndependent tag? matTag?" << endln;
	return 0;
    }

    int tag[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata,tag) < 0) {
	return 0;
    }

    UniaxialMaterial* theMat = OPS_getUniaxialMaterial(tag[1]);
    if (theMat == 0) {
	opserr << "WARNING material does not exist\n";
	opserr << "material: " << tag[1]; 
	opserr << "\nuniaxialMaterial PathIndependent: " << tag[0] << endln;
	return 0;
    }

    UniaxialMaterial* mat = new PathIndependentMaterial(tag[0],*theMat);
    if (mat == 0) {
	opserr << "WARNING: failed to create PathIndependentmaterial material\n";
	return 0;
    }

    return mat;
}

PathIndependentMaterial::PathIndependentMaterial(int tag, UniaxialMaterial &material)
:UniaxialMaterial(tag,MAT_TAG_PathIndependent), theMaterial(0)
{
  theMaterial = material.getCopy();

  if (theMaterial == 0) {
    opserr <<  "PathIndependentMaterial::PathIndependentMaterial -- failed to get copy of material\n";
    //exit(-1);
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
  if (theMaterial)
    return theMaterial->setTrialStrain(strain, strainRate);
  else
    return -1;
}

double 
PathIndependentMaterial::getStress(void)
{
  if (theMaterial)
    return theMaterial->getStress();
  else
    return 0.0;
}


double 
PathIndependentMaterial::getTangent(void)
{
  if (theMaterial)
    return theMaterial->getTangent();
  else
    return 0.0;
}

double 
PathIndependentMaterial::getDampTangent(void)
{
  if (theMaterial)
    return theMaterial->getDampTangent();
  else
    return 0.0;
}

double 
PathIndependentMaterial::getInitialTangent(void)
{
  if (theMaterial)
    return theMaterial->getInitialTangent();
  else
    return 0.0;
}

double 
PathIndependentMaterial::getStrain(void)
{
  if (theMaterial)
    return theMaterial->getStrain();
  else
    return 0.0;
}

double 
PathIndependentMaterial::getStrainRate(void)
{
  if (theMaterial)
    return theMaterial->getStrainRate();
  else
    return 0.0;
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
  if (theMaterial)
    return theMaterial->revertToStart();
  else
    return -1;
}

UniaxialMaterial *
PathIndependentMaterial::getCopy(void)
{
  PathIndependentMaterial *theCopy = 0;
  if (theMaterial)
    theCopy = new PathIndependentMaterial(this->getTag(), *theMaterial);
        
  return theCopy;
}


int 
PathIndependentMaterial::sendSelf(int cTag, Channel &theChannel)
{
  if (theMaterial == 0) {
    opserr << "PathIndependentMaterial::sendSelf() - theMaterial is null, nothing to send\n";
    return -1;
  }
  
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
      //exit(-1);
      return -1;
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
    if (theMaterial)
      s << "\tMaterial: " << theMaterial->getTag() << endln;
    else
      s << "\tMaterial is NULL" << endln;
}

int
PathIndependentMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (theMaterial)
    return theMaterial->setParameter(argv, argc, param);
  else
    return -1;
}

int
PathIndependentMaterial::updateParameter(int parameterID, Information &info)
{
  return 0;
}

double
PathIndependentMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  if (theMaterial)
    return theMaterial->getStressSensitivity(gradIndex, conditional);
  else
    return 0.0;
}

double
PathIndependentMaterial::getStrainSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getStrainSensitivity(gradIndex);
  else
    return 0.0;
}

double
PathIndependentMaterial::getInitialTangentSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getInitialTangentSensitivity(gradIndex);
  else
    return 0.0;
}

double
PathIndependentMaterial::getDampTangentSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getDampTangentSensitivity(gradIndex);
  else
    return 0.0;
}

double
PathIndependentMaterial::getRhoSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getRhoSensitivity(gradIndex);
  else
    return 0.0;
}

int   
PathIndependentMaterial::commitSensitivity(double strainGradient,
					   int gradIndex, int numGrads)
{
  return 0; // commit nothing, path independent
}
