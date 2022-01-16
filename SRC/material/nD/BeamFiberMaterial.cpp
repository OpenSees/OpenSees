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

// $Revision: 1.4 $
// $Date: 2003-02-14 23:01:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BeamFiberMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class definition of BeamFiberMaterial.
// The BeamFiberMaterial class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11, 12, and 13
// stress components which can then be integrated over an area to model a
// shear flexible 3D beam.


#include <BeamFiberMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <elementAPI.h>

Vector BeamFiberMaterial::stress(3);
Matrix BeamFiberMaterial::tangent(3,3);

//      0  1  2  3  4  5
// ND: 11 22 33 12 23 31
// BF: 11 12 31 22 33 23


void * OPS_ADD_RUNTIME_VPV(OPS_BeamFiberMaterial)
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial BeamFiber tag? matTag?" << endln;
	return 0;
    }

    int tags[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, tags) < 0) {
	opserr << "WARNING invalid nDMaterial BeamFiber tag or matTag" << endln;
	return 0;
    }

    int tag = tags[0];
    int matTag = tags[1];

    NDMaterial *threeDMaterial = OPS_getNDMaterial(matTag);
    if (threeDMaterial == 0) {
	opserr << "WARNING nD material does not exist\n";
	opserr << "nD material: " << matTag;
	opserr << "\nBeamFiber nDMaterial: " << tag << endln;
	return 0;
    }

    return new BeamFiberMaterial(tag, *threeDMaterial);
}

BeamFiberMaterial::BeamFiberMaterial(void)
: NDMaterial(0, ND_TAG_BeamFiberMaterial),
Tstrain22(0.0), Tstrain33(0.0), Tgamma23(0.0),
Cstrain22(0.0), Cstrain33(0.0), Cgamma23(0.0),
theMaterial(0), strain(3)
{
	// Nothing to do
}

BeamFiberMaterial::BeamFiberMaterial(int tag, NDMaterial &theMat)
: NDMaterial(tag, ND_TAG_BeamFiberMaterial),
Tstrain22(0.0), Tstrain33(0.0), Tgamma23(0.0),
Cstrain22(0.0), Cstrain33(0.0), Cgamma23(0.0),
theMaterial(0), strain(3)

{
  // Get a copy of the material
  theMaterial = theMat.getCopy("ThreeDimensional");
  
  if (theMaterial == 0) {
    opserr << "BeamFiberMaterial::BeamFiberMaterial -- failed to get copy of material\n";
    exit(-1);
  }
}

BeamFiberMaterial::~BeamFiberMaterial(void) 
{ 
  if (theMaterial != 0)
    delete theMaterial;
} 

NDMaterial*
BeamFiberMaterial::getCopy(void) 
{
  BeamFiberMaterial *theCopy =
    new BeamFiberMaterial(this->getTag(), *theMaterial);

  theCopy->Tstrain22 = this->Tstrain22;
  theCopy->Tstrain33 = this->Tstrain33;
  theCopy->Tgamma23  = this->Tgamma23;
  theCopy->Cstrain22 = this->Cstrain22;
  theCopy->Cstrain33 = this->Cstrain33;
  theCopy->Cgamma23  = this->Cgamma23;
  
  return theCopy;
}

NDMaterial* 
BeamFiberMaterial::getCopy(const char *type)
{
  if (strcmp(type, "BeamFiber") == 0)
    return this->getCopy();
  else
    return 0;
}

int 
BeamFiberMaterial::getOrder(void) const
{
  return 3;
}

const char*
BeamFiberMaterial::getType(void) const 
{
  return "BeamFiber";
}

int 
BeamFiberMaterial::commitState(void)
{
  Cstrain22 = Tstrain22;
  Cstrain33 = Tstrain33;
  Cgamma23 = Tgamma23;

  return theMaterial->commitState();
}

int 
BeamFiberMaterial::revertToLastCommit(void)
{
  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma23 = Cgamma23;
  
  return theMaterial->revertToLastCommit();
}

int
BeamFiberMaterial::revertToStart()
{
  this->Tstrain22 = 0.0;
  this->Tstrain33 = 0.0;
  this->Tgamma23  = 0.0;
  this->Cstrain22 = 0.0;
  this->Cstrain33 = 0.0;
  this->Cgamma23  = 0.0;

  return theMaterial->revertToStart();
}

double
BeamFiberMaterial::getRho(void)
{
  return theMaterial->getRho();
}


//receive the strain
int 
BeamFiberMaterial::setTrialStrain(const Vector &strainFromElement)
{
  static const double tolerance = 1.0e-08;

  strain(0) = strainFromElement(0);
  strain(1) = strainFromElement(1);
  strain(2) = strainFromElement(2);

  //newton loop to solve for out-of-plane strains

  double norm;
  static Vector condensedStress(3);
  static Vector strainIncrement(3);
  static Vector threeDstrain(6);
  static Matrix dd22(3,3);

  int count = 0;
  const int maxCount = 20;
  double norm0;

  do {

    //set three dimensional strain
    threeDstrain(0) = this->strain(0);
    threeDstrain(1) = this->Tstrain22;
    threeDstrain(2) = this->Tstrain33;
    threeDstrain(3) = this->strain(1); 
    threeDstrain(4) = this->Tgamma23;
    threeDstrain(5) = this->strain(2);

    if (theMaterial->setTrialStrain(threeDstrain) < 0) {
      opserr << "BeamFiberMaterial::setTrialStrain - setStrain failed in material with strain " << threeDstrain;
      return -1;   
    }

    //three dimensional stress
    const Vector &threeDstress = theMaterial->getStress();

    //three dimensional tangent 
    const Matrix &threeDtangent = theMaterial->getTangent();

    //NDmaterial strain order        = 11, 22, 33, 12, 23, 31  
    //BeamFiberMaterial strain order = 11, 12, 31, 22, 33, 23

    condensedStress(0) = threeDstress(1);
    condensedStress(1) = threeDstress(2);
    condensedStress(2) = threeDstress(4);

    dd22(0,0) = threeDtangent(1,1);
    dd22(1,0) = threeDtangent(2,1);
    dd22(2,0) = threeDtangent(4,1);

    dd22(0,1) = threeDtangent(1,2);
    dd22(1,1) = threeDtangent(2,2);
    dd22(2,1) = threeDtangent(4,2);

    dd22(0,2) = threeDtangent(1,4);
    dd22(1,2) = threeDtangent(2,4);
    dd22(2,2) = threeDtangent(4,4);

    //set norm
    norm = condensedStress.Norm();
    if (count == 0)
      norm0 = norm;

    //condensation 
    dd22.Solve(condensedStress, strainIncrement);

    //update out of plane strains
    this->Tstrain22 -= strainIncrement(0);
    this->Tstrain33 -= strainIncrement(1);
    this->Tgamma23  -= strainIncrement(2);

  } while (count++ < maxCount && norm > tolerance);

  return 0;
}

const Vector& 
BeamFiberMaterial::getStrain(void)
{
  return strain;
}

const Vector&  
BeamFiberMaterial::getStress()
{
  const Vector &threeDstress = theMaterial->getStress();

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(3);
  stress(2) = threeDstress(5);

  return stress;
}

const Vector& 
BeamFiberMaterial::getStressSensitivity(int gradIndex,
					bool conditional)
{
  const Vector &threeDstress = theMaterial->getStressSensitivity(gradIndex, conditional);

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(3);
  stress(2) = threeDstress(5);

  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix dd12(3,3);
  dd12(0,0) = threeDtangent(0,1);
  dd12(1,0) = threeDtangent(3,1);
  dd12(2,0) = threeDtangent(5,1);

  dd12(0,1) = threeDtangent(0,2);
  dd12(1,1) = threeDtangent(3,2);
  dd12(2,1) = threeDtangent(5,2);

  dd12(0,2) = threeDtangent(0,4);
  dd12(1,2) = threeDtangent(3,4);
  dd12(2,2) = threeDtangent(5,4);


  static Matrix dd22(3,3);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);
  
  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);
  
  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);

  
  static Vector sigma2(3);
  sigma2(0) = threeDstress(1);
  sigma2(1) = threeDstress(2);
  sigma2(2) = threeDstress(4);

  static Vector dd22sigma2(3);
  dd22.Solve(sigma2,dd22sigma2);

  stress.addMatrixVector(1.0, dd12, dd22sigma2, -1.0);

  return stress;
}

const Matrix&  
BeamFiberMaterial::getTangent()
{
  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix dd11(3,3);
  dd11(0,0) = threeDtangent(0,0);
  dd11(1,0) = threeDtangent(3,0);
  dd11(2,0) = threeDtangent(5,0);

  dd11(0,1) = threeDtangent(0,3);
  dd11(1,1) = threeDtangent(3,3);
  dd11(2,1) = threeDtangent(5,3);

  dd11(0,2) = threeDtangent(0,5);
  dd11(1,2) = threeDtangent(3,5);
  dd11(2,2) = threeDtangent(5,5);


  static Matrix dd12(3,3);
  dd12(0,0) = threeDtangent(0,1);
  dd12(1,0) = threeDtangent(3,1);
  dd12(2,0) = threeDtangent(5,1);

  dd12(0,1) = threeDtangent(0,2);
  dd12(1,1) = threeDtangent(3,2);
  dd12(2,1) = threeDtangent(5,2);

  dd12(0,2) = threeDtangent(0,4);
  dd12(1,2) = threeDtangent(3,4);
  dd12(2,2) = threeDtangent(5,4);

  static Matrix dd21(3,3);
  dd21(0,0) = threeDtangent(1,0);
  dd21(1,0) = threeDtangent(2,0);
  dd21(2,0) = threeDtangent(4,0);

  dd21(0,1) = threeDtangent(1,3);
  dd21(1,1) = threeDtangent(2,3);
  dd21(2,1) = threeDtangent(4,3);

  dd21(0,2) = threeDtangent(1,5);
  dd21(1,2) = threeDtangent(2,5);
  dd21(2,2) = threeDtangent(4,5);

  static Matrix dd22(3,3);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);

  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);

  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);


  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(3,3);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11; 

  return tangent;
}

const Matrix&  
BeamFiberMaterial::getInitialTangent()
{
  const Matrix &threeDtangent = theMaterial->getInitialTangent();

  static Matrix dd11(3,3);
  dd11(0,0) = threeDtangent(0,0);
  dd11(1,0) = threeDtangent(3,0);
  dd11(2,0) = threeDtangent(5,0);

  dd11(0,1) = threeDtangent(0,3);
  dd11(1,1) = threeDtangent(3,3);
  dd11(2,1) = threeDtangent(5,3);

  dd11(0,2) = threeDtangent(0,5);
  dd11(1,2) = threeDtangent(3,5);
  dd11(2,2) = threeDtangent(5,5);


  static Matrix dd12(3,3);
  dd12(0,0) = threeDtangent(0,1);
  dd12(1,0) = threeDtangent(3,1);
  dd12(2,0) = threeDtangent(5,1);

  dd12(0,1) = threeDtangent(0,2);
  dd12(1,1) = threeDtangent(3,2);
  dd12(2,1) = threeDtangent(5,2);

  dd12(0,2) = threeDtangent(0,4);
  dd12(1,2) = threeDtangent(3,4);
  dd12(2,2) = threeDtangent(5,4);

  static Matrix dd21(3,3);
  dd21(0,0) = threeDtangent(1,0);
  dd21(1,0) = threeDtangent(2,0);
  dd21(2,0) = threeDtangent(4,0);

  dd21(0,1) = threeDtangent(1,3);
  dd21(1,1) = threeDtangent(2,3);
  dd21(2,1) = threeDtangent(4,3);

  dd21(0,2) = threeDtangent(1,5);
  dd21(1,2) = threeDtangent(2,5);
  dd21(2,2) = threeDtangent(4,5);

  static Matrix dd22(3,3);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);

  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);

  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);


  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(3,3);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11; 

  return tangent;
}

void  
BeamFiberMaterial::Print(OPS_Stream &s, int flag)
{
  s << "BeamFiberMaterial, tag: " << this->getTag() << endln;
  s << "\tWrapped material: "<< theMaterial->getTag() << endln;

  theMaterial->Print(s, flag);
}

int 
BeamFiberMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  // put tag and associated materials class and database tags into an id and send it
  static ID idData(3);
  idData(0) = this->getTag();
  idData(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;

  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "BeamFiberMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(3);
  vecData(0) = Cstrain22;
  vecData(1) = Cstrain33;
  vecData(2) = Cgamma23;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "BeamFiberMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "BeamFiberMaterial::sendSelf() - failed to send vector material\n";

  return res;
}

int 
BeamFiberMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "BeamFiberMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  
  // if the associated material has not yet been created or is of the wrong type
  // create a new material for recvSelf later
  if (theMaterial == 0 || theMaterial->getClassTag() != matClassTag) {
    if (theMaterial != 0)
      delete theMaterial;
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "BeamFiberMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(3);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "BeamFiberMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);
  Cstrain33 = vecData(1);
  Cgamma23  = vecData(2);

  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma23  = Cgamma23;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "BeamFiberMaterial::sendSelf() - failed to send vector material\n";
  
  return res;
}

int
BeamFiberMaterial::setParameter(const char **argv, int argc,
				Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}
