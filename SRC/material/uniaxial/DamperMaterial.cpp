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
                                                                        
// Written: fmk 
// Created: 12/17
// Revision: A
//
// Description: This file contains the class definition for 
// DamperModel. DamperModel is an aggregation
// of UniaxialMaterial objects all considered acting in parallel.
//
// What: "@(#) DamperModel.C, revA"

#include <DamperMaterial.h>
#include <ID.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <stdlib.h>
#include <string.h>
#include <MaterialResponse.h>

#include <elementAPI.h>

void *
OPS_DamperMaterial(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int argc = OPS_GetNumRemainingInputArgs();

  if (argc < 2) {
    opserr << "Invalid #args,  want: uniaxialMaterial Damper $tag $tag1 $tag2 ... <-factors $fact1 $fact2 ...>" << endln;
    return 0;
  }

  int numData = 2;
  int *iData = new int[2];
  
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid data for uniaxialMaterial Damper" << endln;
    return 0;
  }

  UniaxialMaterial *theMat = OPS_getUniaxialMaterial(iData[1]);
  if (theMat == 0) {
    opserr << "WARNING no existing material with tag " << iData[1] 
	   << " for uniaxialMaterial Damper" << iData[0] << endln;
    delete [] iData;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new DamperMaterial(iData[0], theMat);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type Damper\n";
    return 0;
  }
  
  delete [] iData;

  return theMaterial;
}



DamperMaterial::DamperMaterial(int tag, 
			       UniaxialMaterial *theMaterialModel)
:UniaxialMaterial(tag, MAT_TAG_DamperMaterial),
 trialStrain(0.0), trialStrainRate(0.0), theMaterial(0)
{
  theMaterial = theMaterialModel->getCopy();
}



// this constructor is used for a DamperMaterailModel object that
// needs to be constructed in a remote actor process. recvSelf() needs
// to be called on this object
DamperMaterial::DamperMaterial()
:UniaxialMaterial(0,MAT_TAG_DamperMaterial),
 trialStrain(0.0), trialStrainRate(0.0), theMaterial(0)
{

}


DamperMaterial::~DamperMaterial()
{
  delete theMaterial;
}



int 
DamperMaterial::setTrialStrain(double strain, double strainRate)
{
    // set the trialStrain and the trialStrain in each of the
    // local MaterialModel objects 
    trialStrain = strain;
    trialStrainRate = strainRate;

    theMaterial->setTrialStrain(strainRate, 0);

    return 0;
}


double 
DamperMaterial::getStrain(void)
{
  return trialStrain;
}

double 
DamperMaterial::getStrainRate(void)
{
    return trialStrainRate;
}

double 
DamperMaterial::getStress(void)
{
  return theMaterial->getStress();
}



double 
DamperMaterial::getTangent(void)
{
  return 0;
}

double 
DamperMaterial::getInitialTangent(void)
{
  return 0;
}

double 
DamperMaterial::getDampTangent(void)
{
  return theMaterial->getTangent();
}

int 
DamperMaterial::commitState(void)
{
  return theMaterial->commitState();
}

int 
DamperMaterial::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}


int 
DamperMaterial::revertToStart(void)
{
    trialStrain = 0.0;
    trialStrainRate = 0.0;
    theMaterial->revertToStart();
    
    return 0;    
}



UniaxialMaterial *
DamperMaterial::getCopy(void)
{
  DamperMaterial *theCopy = new 
    DamperMaterial(this->getTag(), theMaterial);
  
  theCopy->trialStrain = trialStrain;
  theCopy->trialStrainRate = trialStrainRate;
  
  return theCopy;
}


int 
DamperMaterial::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;

    static ID data(3);

    // send ID of size 3 so no possible conflict with classTags ID
    int dbTag = this->getDbTag();
    data(0) = this->getTag();
    data(1) = theMaterial->getClassTag();
    int matDbTag = theMaterial->getDbTag();
    if (matDbTag == 0) {
      matDbTag  = theChannel.getDbTag();
      if (matDbTag != 0)
	    theMaterial->setDbTag(matDbTag);
    }

    data(2) = matDbTag;

    res = theChannel.sendID(dbTag, cTag, data);
    if (res < 0) {
      opserr << "DamperMaterial::sendSelf() - failed to send data\n";
      return res;
    }

    return theMaterial->sendSelf(cTag, theChannel);
}

int 
DamperMaterial::recvSelf(int cTag, Channel &theChannel, 
				FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static ID data(3);
    int dbTag = this->getDbTag();

    res = theChannel.recvID(dbTag, cTag, data);
    if (res < 0) {
      opserr << "DamperMaterial::recvSelf() - failed to receive data\n";
      return res;
    }

    this->setTag(int(data(0)));
    if (theMaterial != 0)
      delete theMaterial;

    int matClassTag = data[1];
    dbTag = data[2];
    UniaxialMaterial *theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial  == 0) {
      opserr << "FATAL DamperMaterial::recvSelf() ";
      opserr << " could not get a UniaxialMaterial \n";
      exit(-1);
    }    	    
    theMaterial->setDbTag(dbTag);
    theMaterial->recvSelf(cTag, theChannel, theBroker);

    return 0;
}

void 
DamperMaterial::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
      s << "DamperMaterial tag: " << this->getTag() << endln;
      s << " ";
      theMaterial->Print(s, flag);
    }
}


