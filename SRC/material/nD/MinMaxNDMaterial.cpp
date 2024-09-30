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

// MHS

#include <stdlib.h>
#include <math.h>

#include <MinMaxNDMaterial.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <Parameter.h>
#include <OPS_Globals.h>

#include <elementAPI.h>

void*
OPS_MinMaxNDMaterial(void)
{
  // Pointer to an ND material that will be returned
  NDMaterial* theMaterial = 0;
  NDMaterial* theOtherMaterial = 0;
  double minStrain = -1.0e16;
  double maxStrain = 1.0e16;
  int    iData[2];
  
  int argc = OPS_GetNumRemainingInputArgs();
  if (argc < 2) {
    opserr << "WARNING invalid nDMaterial MinMax $tag $otherTag <-min $minStrain> <-max $maxStrain>" << endln;
    return 0;
  }

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid nDMaterial MinMax $tag $otherTag" << endln;
    return 0;
  }

  theOtherMaterial = OPS_getNDMaterial(iData[1]);
  if (theOtherMaterial == 0) {
    opserr << "WARNING invalid otherTag nDMaterial MinMax tag: " << iData[0] << endln;
    return 0;	
  }

  argc = OPS_GetNumRemainingInputArgs();  
  while (argc > 1) {
    const char *argvLoc = OPS_GetString();
    numData = 1;

    if ((strcmp(argvLoc, "-min") == 0) || (strcmp(argvLoc, "-Min") == 0) || (strcmp(argvLoc, "-MIN") == 0)) {
      if (OPS_GetDouble(&numData, &minStrain) != 0) {      
	opserr << "WARNING invalid min value  nDMaterial MinMax tag: " << iData[0] << endln;	
	return 0;
      }
    } else if ((strcmp(argvLoc, "-max") == 0) || (strcmp(argvLoc, "-Max") == 0) || (strcmp(argvLoc, "-MAX") == 0)) {
      if (OPS_GetDouble(&numData, &maxStrain) != 0) {      
	opserr << "WARNING invalid min value  nDMaterial MinMax tag: " << iData[0] << endln;  
	return 0;
      }
    } else {
      opserr << "WARNING invalid option:" << argvLoc << " uniaxialMaterial MinMax tag: " << iData[0] << endln;  
      return 0;
    }
    
    argc = OPS_GetNumRemainingInputArgs();
  }

  // Parsing was successful, allocate the material
  theMaterial = new MinMaxNDMaterial(iData[0], *theOtherMaterial, minStrain, maxStrain);

  if (theMaterial == 0) {
    opserr << "WARNING could not create nDMaterial of type MinMax\n";
    return 0;
  }
    
  return theMaterial;
}

MinMaxNDMaterial::MinMaxNDMaterial(int tag, NDMaterial& material, double epsmin, double epsmax)
  : NDMaterial(tag, ND_TAG_MinMaxNDMaterial),
    theMaterial(0), minStrain(epsmin), maxStrain(epsmax),
    Tfailed(false), Cfailed(false), myType(ThreeDimensional)
{
  // If passed a base class material, get a 3D copy
  if (strncmp(material.getType(),"BaseClass",80) == 0) {
    opserr << "MinMaxNDMaterial - base class material passed, assuming 3D" << endln;
    theMaterial = material.getCopy("ThreeDimensional");
  }
  else {
    // get copy of the main material
    theMaterial = material.getCopy();
  }

  if (theMaterial == 0) {
    opserr << "MinMaxNDMaterial::MinMaxNDMaterial -- failed to get copy of material (a 3D material is required)" << endln;
    exit(-1);
  }

  const char *type = theMaterial->getType();
  if (strncmp(type,"ThreeDimensional",80) == 0)
    myType = ThreeDimensional;
  if (strncmp(type,"PlateFiber",80) == 0)
    myType = PlateFiber;
  if (strncmp(type,"PlaneStress",80) == 0 ||
      strncmp(type,"PlaneStress2D",80) == 0)
    myType = PlaneStress2d;
  if (strncmp(type,"BeamFiber",80) == 0 ||
      strncmp(type,"TimoshenkoFiber",80) == 0)
    myType = BeamFiber;
  if (strncmp(type,"BeamFiber2d",80) == 0 ||
      strncmp(type,"TimoshenkoFiber2d",80) == 0)
    myType = BeamFiber2d;    
  if (strncmp(type,"BeamFiber2dPS",80) == 0 ||
      strncmp(type,"TimoshenkoFiber2dPS",80) == 0)
    myType = BeamFiber2dPS;    
}

MinMaxNDMaterial::MinMaxNDMaterial()
  : NDMaterial(0, ND_TAG_MinMaxNDMaterial),
    theMaterial(0), minStrain(0.0), maxStrain(0.0),
    Tfailed(false), Cfailed(false), myType(ThreeDimensional)
{

}

MinMaxNDMaterial::~MinMaxNDMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int
MinMaxNDMaterial::setTrialStrain(const Vector& strain)
{
  if (Cfailed)
    return 0;

  double eps11 = strain(0);
  double eps22 = 0.0;
  double eps12 = 0.0;
  if (myType == PlaneStress2d) {
    eps22 = strain(1);
    eps12 = strain(2);
  }
  if (myType == BeamFiber || myType == BeamFiber2d || myType == BeamFiber2dPS) {
    eps12 = strain(1);
  }
  if (myType == PlateFiber) {
    eps22 = strain(1);
    eps12 = strain(2);
  }
  if (myType == ThreeDimensional) {
    eps22 = strain(1);
    eps12 = strain(3);
  }    

  double a = 0.5*(eps11-eps22);
  double b = 0.5*eps12;
  double ab = pow(a*a+b*b,0.5);
  double eps1 = 0.5*(eps11+eps22) + ab;
  double eps2 = 0.5*(eps11+eps22) - ab;

  if (eps1 >= maxStrain || eps2 <= minStrain) {
    Tfailed = true;
    return 0;
  }
  else {
    Tfailed = false;
    return theMaterial->setTrialStrain(strain);
  }
}

int
MinMaxNDMaterial::setTrialStrain(const Vector& strain, const Vector& /*strainRate*/)
{
  return this->setTrialStrain(strain);
}

int
MinMaxNDMaterial::setTrialStrainIncr(const Vector& strain)
{
  if (Cfailed)
    return 0;

  double eps11 = strain(0);
  double eps22 = strain(1);
  double eps12 = strain(3);
  double a = 0.5*(eps11-eps22);
  double b = 0.5*eps12;
  double ab = pow(a*a+b*b,0.5);
  double eps1 = 0.5*(eps11+eps22) + ab;
  double eps2 = 0.5*(eps11+eps22) - ab;
  
  if (eps1 >= maxStrain || eps2 <= minStrain) {
    Tfailed = true;
    return 0;
  }
  else {
    Tfailed = false;
    static Vector strain_from_ele(6);
    strain_from_ele = theMaterial->getStrain();
    strain_from_ele.addVector(1.0, strain, 1.0);
    return setTrialStrain(strain_from_ele);
  }
}

int
MinMaxNDMaterial::setTrialStrainIncr(const Vector& strain, const Vector& /*strainRate*/)
{
  return setTrialStrainIncr(strain);
}

const Vector&
MinMaxNDMaterial::getStress(void)
{
  if (Tfailed) {
    static Vector zeroStress(6);
    return zeroStress;
  }
  else
    return theMaterial->getStress();
}

const Matrix&
MinMaxNDMaterial::getTangent(void)
{
  if (Tfailed) {
    static Matrix zeroTangent(6,6);
    zeroTangent = theMaterial->getInitialTangent();
    zeroTangent *= 1e-8;
    return zeroTangent;
  }
  else
    return theMaterial->getTangent();
}

const Matrix&
MinMaxNDMaterial::getInitialTangent(void)
{
  return theMaterial->getInitialTangent();
}

const Vector&
MinMaxNDMaterial::getStrain(void)
{
  return theMaterial->getStrain();
}

int
MinMaxNDMaterial::commitState(void)
{
  Cfailed = Tfailed;

  if (Tfailed)
    return 0;
  else
    return theMaterial->commitState();
}

int
MinMaxNDMaterial::revertToLastCommit(void)
{
  Tfailed = Cfailed;
  
  // Check if failed at last step
  if (Cfailed)
    return 0;
  else
    return theMaterial->revertToLastCommit();
}

int
MinMaxNDMaterial::revertToStart(void)
{
  Cfailed = false;
  Tfailed = false;
  
  return theMaterial->revertToStart();
}

double
MinMaxNDMaterial::getRho(void)
{
  return theMaterial->getRho();
}

NDMaterial*
MinMaxNDMaterial::getCopy(void)
{
  MinMaxNDMaterial* theCopy = new MinMaxNDMaterial(this->getTag(), *theMaterial, minStrain, maxStrain);

  theCopy->Cfailed = Cfailed;
  theCopy->Tfailed = Tfailed;
  
  return theCopy;
}

NDMaterial *
MinMaxNDMaterial::getCopy(const char *type)
{
  NDMaterial *copy = theMaterial->getCopy(type);
  MinMaxNDMaterial *theCopy = new MinMaxNDMaterial(this->getTag(), *copy, minStrain, maxStrain);
  delete copy;

  if (strncmp(type,"PlaneStress",80) == 0 ||
      strncmp(type,"PlaneStress2D",80) == 0)
    theCopy->myType = PlaneStress2d;
  if (strncmp(type,"BeamFiber",80) == 0 ||
      strncmp(type,"TimoshenkoFiber",80) == 0)
    theCopy->myType = BeamFiber;
  if (strncmp(type,"BeamFiber2d",80) == 0 ||
      strncmp(type,"TimoshenkoFiber2d",80) == 0)
    theCopy->myType = BeamFiber2d;    
  if (strncmp(type,"BeamFiber2dPS",80) == 0 ||
      strncmp(type,"TimoshenkoFiber2dPS",80) == 0)
    theCopy->myType = BeamFiber2dPS;
  if (strncmp(type,"PlateFiber",80) == 0)
    theCopy->myType = PlateFiber;
  if (strncmp(type,"ThreeDimensonal",80) == 0)
    theCopy->myType = ThreeDimensional;      
  
  theCopy->Cfailed = Cfailed;
  theCopy->Tfailed = Tfailed;
  
  return theCopy;  
}

const char*
MinMaxNDMaterial::getType(void) const
{
  return theMaterial->getType();
}

int MinMaxNDMaterial::getOrder(void) const
{
  return theMaterial->getOrder();
}

int
MinMaxNDMaterial::sendSelf(int cTag, Channel& theChannel)
{
  if (theMaterial == 0) {
    opserr << "MinMaxNDMaterial::sendSelf() - theMaterial is null, nothing to send" << endln;
    return -1;
  }

  int dbTag = this->getDbTag();
  
  static ID dataID(4);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "MinMaxNDMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }
  dataID(3) = myType;
  
  static Vector dataVec(3);
  dataVec(0) = minStrain;
  dataVec(1) = maxStrain;
  if (Cfailed == true)
    dataVec(2) = 1.0;
  else
    dataVec(2) = 0.0;
  
  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "MinMaxNDMaterial::sendSelf() - failed to send epsInit\n";
    return -2;
  }
  
  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "MinMaxNDMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }
  
  return 0;
}

int
MinMaxNDMaterial::recvSelf(int cTag, Channel& theChannel,
			   FEM_ObjectBroker& theBroker)
{
  int dbTag = this->getDbTag();
  
  static ID dataID(4);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "MinMaxNDMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(dataID(0));
  
  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = dataID(1);
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "MinMaxNDMaterial::recvSelf() - failed to create Material with classTag "
	     << matClassTag << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));
  myType = (type)dataID(3);
  
  static Vector dataVec(3);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "MinMaxNDMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  minStrain = dataVec(0);
  maxStrain = dataVec(1);
  
  if (dataVec(2) == 1.0)
    Cfailed = true;
  else
    Cfailed = false;

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "MinMaxNDMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }

  this->revertToLastCommit();
  
  return 0;
}

void
MinMaxNDMaterial::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"MinMaxNDMaterial\", ";    
    s << "\"material\": \"" << theMaterial->getTag() << "\", ";
    s << "\"minStrain\": " << minStrain << ", ";
    s << "\"maxStrain\": " << maxStrain;
    s << "\t\t\t}";
  }
  else {
    s << "MinMaxNDMaterial tag: " << this->getTag() << endln;
    s << "\tMaterial: " << theMaterial->getTag() << endln;
    s << "  min strain: " << minStrain << endln;
    s << "  max strain: " << maxStrain << endln;
  }
}

int
MinMaxNDMaterial::setParameter(const char** argv, int argc, Parameter& param)
{
  //
  // I suppose epsMin and epsMax could be parameters, but not for now -- MHS
  //
  return theMaterial->setParameter(argv, argc, param);  
}

int MinMaxNDMaterial::updateParameter(int parameterID, Information& info)
{
  return 0;
}

Response* MinMaxNDMaterial::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  return theMaterial->setResponse(argv, argc, output);
}

const Vector&
MinMaxNDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  return theMaterial->getStressSensitivity(gradIndex, conditional);
}

int
MinMaxNDMaterial::commitSensitivity(const Vector& depsdh,
				    int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(depsdh, gradIndex, numGrads);
}
