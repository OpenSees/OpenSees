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
// Created: Dec 2012

#include <stdlib.h>
#include <math.h>

#include <InitStressNDMaterial.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>

#include <OPS_Globals.h>

#include <elementAPI.h>
#define OPS_Export 

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_InitStressNDMaterial)
{
  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;
  NDMaterial *theOtherMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3) {
    opserr << "Want: nDMaterial InitStress tag? otherTag? sig0? <nDim?>" << endln;
  }

  int    iData[2];
  double dData[1];
  int    dim[1];
  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid nDMaterial InitStressNDMaterial $tag $otherTag $nDim" << endln;
    return 0;
  }

  theOtherMaterial = OPS_getNDMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "Could not find material with tag: " << iData[1] << "nDMaterial InitStress $tag $otherTag $nDim $sig0" << endln;
    return 0;	
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: nDMaterial InitStress $tag $otherTag $nDim $sig0" << endln;
    return 0;	
  }

  if (numArgs == 4) {
    if (OPS_GetIntInput(&numData, dim) != 0) {
        return 0;
    }
  } else {
    dim[0] = 3;
  }

  Vector sig0(3*dim[0]-3);
  if (dim[0] == 3) {
    sig0(0) = dData[0];
    sig0(1) = dData[0];
    sig0(2) = dData[0];
  } else if (dim[0] == 2) {
    sig0(0) = dData[0];
    sig0(1) = dData[0];
  } else {
    opserr << "nDMaterial InitStress - Invalid number of dimensions: want 2 or 3" << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  if (numArgs == 4) {
    theMaterial = new InitStressNDMaterial(iData[0], *theOtherMaterial, sig0, dim[0]);
  } else {
    theMaterial = new InitStressNDMaterial(iData[0], *theOtherMaterial, sig0);
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type InitStressNDMaterial\n";
    return 0;
  }

  return theMaterial;
}


InitStressNDMaterial::InitStressNDMaterial(int tag, NDMaterial &material, const Vector &sigini, int ndim)
  :NDMaterial(tag,ND_TAG_InitStressNDMaterial), theMaterial(0),
   epsInit(3*ndim-3), sigInit(sigini)
{

  numDim = ndim;
  // get copy of the main material
  if (numDim == 2) {
    theMaterial = material.getCopy("PlaneStrain");
  } else if (numDim == 3) {
    theMaterial = material.getCopy("ThreeDimensional");
  } else {
    opserr << "nDMaterial InitStress - Invalid number of dimensions: want 2 or 3" << endln;
  }

  if (theMaterial == 0) {
    opserr <<  "InitStressNDMaterial::InitStressNDMaterial -- failed to get copy of material\n";
    exit(-1);
  }

  // determine the initial strain
  int mDim = 3*numDim-3;
  double tol=1e-12;
  Vector dSig(sigInit);
  Vector dStrain(mDim);
  Vector tStrain(mDim);
  Vector tStress(mDim);
  Matrix K(mDim,mDim);
  int count = 0;

  do {
    count++;
    K = theMaterial->getTangent();
    K.Solve(dSig,dStrain);
    tStrain += dStrain;
    theMaterial->setTrialStrain(tStrain);
    tStress = theMaterial->getStress();
    dSig = sigInit-tStress;
    dStrain = tStress-sigInit;
  } while (dStrain.Norm() > tol && (count <= 100));

  epsInit = tStrain;

  if (dStrain.Norm() < tol) 
    theMaterial->setTrialStrain(epsInit);
  else {
    opserr << "WARNING: InitStressNDMaterial - could not find initStrain to within tol for material: " << tag;
    opserr << " wanted sigInit: " << sigInit << " using tStress: " << theMaterial->getStress() << endln;
  }

  theMaterial->commitState();
}

InitStressNDMaterial::InitStressNDMaterial()
  :NDMaterial(0,ND_TAG_InitStressNDMaterial), theMaterial(0),
   epsInit(6), sigInit(6)
{

}

InitStressNDMaterial::~InitStressNDMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
InitStressNDMaterial::setTrialStrain(const Vector &strain) 
{
  return theMaterial->setTrialStrain(strain+epsInit);
}

int 
InitStressNDMaterial::setTrialStrain(const Vector &strain, 
				     const Vector &strainRate)
{
  return theMaterial->setTrialStrain(strain+epsInit, strainRate);
}

int 
InitStressNDMaterial::setTrialStrainIncr(const Vector &strain) 
{
  return theMaterial->setTrialStrainIncr(strain);
}

int 
InitStressNDMaterial::setTrialStrainIncr(const Vector &strain, 
					 const Vector &strainRate)
{
  return theMaterial->setTrialStrainIncr(strain, strainRate);
}

const Vector &
InitStressNDMaterial::getStress(void)
{
  return theMaterial->getStress();
}

const Matrix &
InitStressNDMaterial::getTangent(void)
{
  return theMaterial->getTangent();  
}

const Matrix &
InitStressNDMaterial::getInitialTangent(void)
{
  return theMaterial->getInitialTangent();  
}

const Vector & 
InitStressNDMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

int 
InitStressNDMaterial::commitState(void)
{	
  return theMaterial->commitState();
}

int 
InitStressNDMaterial::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}

int 
InitStressNDMaterial::revertToStart(void)
{
  int res = 0;
  res = theMaterial->revertToStart();
  res += theMaterial->setTrialStrain(epsInit);
  res += theMaterial->commitState();
  return res;
}

double 
InitStressNDMaterial::getRho(void) 
{
	return theMaterial->getRho();
}

NDMaterial *
InitStressNDMaterial::getCopy(void)
{
  InitStressNDMaterial *theCopy = 
    new InitStressNDMaterial(this->getTag(), *theMaterial, sigInit, numDim);
        
  return theCopy;
}

NDMaterial *
InitStressNDMaterial::getCopy(const char *type)
{
  /*if (strcmp(type,"ThreeDimensional") == 0) {
    InitStressNDMaterial *theCopy = 
      new InitStressNDMaterial(this->getTag(), *theMaterial, sigInit);

    return theCopy;
  }

  return NDMaterial::getCopy(type);*/
  
  return this->getCopy();
}

const char*
InitStressNDMaterial::getType(void) const
{
  return theMaterial->getType();
}

int 
InitStressNDMaterial::sendSelf(int cTag, Channel &theChannel)
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
    opserr << "InitStressNDMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(1);
  //dataVec(0) = epsInit;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStressNDMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "InitStressNDMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
InitStressNDMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "InitStressNDMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "InitStressNDMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(1);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStressNDMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  //epsInit = dataVec(0);
  
  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "InitStressNDMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
InitStressNDMaterial::Print(OPS_Stream &s, int flag)
{
  s << "InitStressNDMaterial tag: " << this->getTag() << endln;
  s << "\tMaterial: " << theMaterial->getTag() << endln;
  s << "\tinitital strain: " << epsInit << endln;
}

int 
InitStressNDMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

const Vector &
InitStressNDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  return theMaterial->getStressSensitivity(gradIndex, conditional);
}

int
InitStressNDMaterial::commitSensitivity(const Vector &depsdh, 
					int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(depsdh, gradIndex, numGrads);
}
