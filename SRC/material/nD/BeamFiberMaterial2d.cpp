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
// Created: Aug 2001
//
// Description: This file contains the class definition of BeamFiberMaterial.
// The BeamFiberMaterial class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11 and 12
// stress components which can then be integrated over an area to model a
// shear flexible 2D beam.


#include <BeamFiberMaterial2d.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <elementAPI.h>

Vector BeamFiberMaterial2d::stress(2);
Matrix BeamFiberMaterial2d::tangent(2,2);

void* OPS_BeamFiberMaterial2d()
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial BeamFiber2d tag? matTag?" << endln;
	return 0;
    }

    int tags[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, tags) < 0) {
	opserr << "WARNING invalid nDMaterial BeamFiber2d tag or matTag" << endln;
	return 0;
    }

    int tag = tags[0];
    int matTag = tags[1];

    NDMaterial *threeDMaterial = OPS_getNDMaterial(matTag);
    if (threeDMaterial == 0) {
	opserr << "WARNING nD material does not exist\n";
	opserr << "nD material: " << matTag;
	opserr << "\nBeamFiber2d nDMaterial: " << tag << endln;
	return 0;
    }

    return new BeamFiberMaterial2d(tag, *threeDMaterial);
}

//      0  1  2  3  4  5
// ND: 11 22 33 12 23 31
// BF: 11 12 22 33 23 31

BeamFiberMaterial2d::BeamFiberMaterial2d(void)
  :NDMaterial(0, ND_TAG_BeamFiberMaterial2d),
   Tstrain22(0.0), Tstrain33(0.0), Tgamma31(0.0), Tgamma23(0.0),
   Cstrain22(0.0), Cstrain33(0.0), Cgamma31(0.0), Cgamma23(0.0),
   theMaterial(0), strain(2)
{
	// Nothing to do
}

BeamFiberMaterial2d::BeamFiberMaterial2d(int tag, NDMaterial &theMat)
  :NDMaterial(tag, ND_TAG_BeamFiberMaterial2d),
   Tstrain22(0.0), Tstrain33(0.0), Tgamma31(0.0), Tgamma23(0.0),
   Cstrain22(0.0), Cstrain33(0.0), Cgamma31(0.0), Cgamma23(0.0),
   theMaterial(0), strain(2)

{
  // Get a copy of the material
  theMaterial = theMat.getCopy("ThreeDimensional");
  
  if (theMaterial == 0) {
    opserr << "BeamFiberMaterial2d::BeamFiberMaterial2d -- failed to get copy of material\n";
    exit(-1);
  }
}

BeamFiberMaterial2d::~BeamFiberMaterial2d(void) 
{ 
  if (theMaterial != 0)
    delete theMaterial;
} 

NDMaterial*
BeamFiberMaterial2d::getCopy(void) 
{
  BeamFiberMaterial2d *theCopy =
    new BeamFiberMaterial2d(this->getTag(), *theMaterial);
  
  theCopy->Tstrain22 = this->Tstrain22;
  theCopy->Tstrain33 = this->Tstrain33;
  theCopy->Tgamma31  = this->Tgamma31;
  theCopy->Tgamma23  = this->Tgamma23;
  theCopy->Cstrain22 = this->Cstrain22;
  theCopy->Cstrain33 = this->Cstrain33;
  theCopy->Cgamma31  = this->Cgamma31;
  theCopy->Cgamma23  = this->Cgamma23;
  
  return theCopy;
}

NDMaterial* 
BeamFiberMaterial2d::getCopy(const char *type)
{
  if (strcmp(type, "BeamFiber2d") == 0)
    return this->getCopy();
  else
    return 0;
}

int 
BeamFiberMaterial2d::getOrder(void) const
{
  return 2;
}

const char*
BeamFiberMaterial2d::getType(void) const 
{
  return "BeamFiber2d";
}

int 
BeamFiberMaterial2d::commitState(void)
{
  Cstrain22 = Tstrain22;
  Cstrain33 = Tstrain33;
  Cgamma31 = Tgamma31;
  Cgamma23 = Tgamma23;

  return theMaterial->commitState();
}

int 
BeamFiberMaterial2d::revertToLastCommit(void)
{
  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma31 = Cgamma31;
  Tgamma23 = Cgamma23;
  
  return theMaterial->revertToLastCommit();
}

int
BeamFiberMaterial2d::revertToStart()
{
  this->Tstrain22 = 0.0;
  this->Tstrain33 = 0.0;
  this->Tgamma31  = 0.0;
  this->Tgamma23  = 0.0;
  this->Cstrain22 = 0.0;
  this->Cstrain33 = 0.0;
  this->Cgamma31  = 0.0;
  this->Cgamma23  = 0.0;

  strain.Zero();

  return theMaterial->revertToStart();
}

double
BeamFiberMaterial2d::getRho(void)
{
  return theMaterial->getRho();
}


//receive the strain
int 
BeamFiberMaterial2d::setTrialStrain(const Vector &strainFromElement)
{
  static const double tolerance = 1.0e-12;

  strain(0) = strainFromElement(0);
  strain(1) = strainFromElement(1);

  //newton loop to solve for out-of-plane strains

  double norm;
  static Vector condensedStress(4);
  static Vector strainIncrement(4);
  static Vector threeDstrain(6);
  static Matrix dd22(4,4);

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
    threeDstrain(5) = this->Tgamma31;

    if (theMaterial->setTrialStrain(threeDstrain) < 0) {
      opserr << "BeamFiberMaterial2d::setTrialStrain - setStrain failed in material with strain " << threeDstrain;
      return -1;   
    }

    //three dimensional stress
    const Vector &threeDstress = theMaterial->getStress();

    //three dimensional tangent 
    const Matrix &threeDtangent = theMaterial->getTangent();

    condensedStress(0) = threeDstress(1);
    condensedStress(1) = threeDstress(2);
    condensedStress(2) = threeDstress(4);
    condensedStress(3) = threeDstress(5);

    dd22(0,0) = threeDtangent(1,1);
    dd22(1,0) = threeDtangent(2,1);
    dd22(2,0) = threeDtangent(4,1);
    dd22(3,0) = threeDtangent(5,1);

    dd22(0,1) = threeDtangent(1,2);
    dd22(1,1) = threeDtangent(2,2);
    dd22(2,1) = threeDtangent(4,2);
    dd22(3,1) = threeDtangent(5,2);

    dd22(0,2) = threeDtangent(1,4);
    dd22(1,2) = threeDtangent(2,4);
    dd22(2,2) = threeDtangent(4,4);
    dd22(3,2) = threeDtangent(5,4);

    dd22(0,3) = threeDtangent(1,5);
    dd22(1,3) = threeDtangent(2,5);
    dd22(2,3) = threeDtangent(4,5);
    dd22(3,3) = threeDtangent(5,5);

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
    this->Tgamma31  -= strainIncrement(3);

  } while (count++ < maxCount && norm > tolerance*norm0);

  return 0;
}

const Vector& 
BeamFiberMaterial2d::getStrain(void)
{
  return strain;
}

const Vector&  
BeamFiberMaterial2d::getStress()
{
  const Vector &threeDstress = theMaterial->getStress();

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(3);

  return stress;
}

const Vector& 
BeamFiberMaterial2d::getStressSensitivity(int gradIndex,
					  bool conditional)
{
  const Vector &threeDstress = theMaterial->getStressSensitivity(gradIndex, conditional);

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(3);

  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix dd12(2,4);
  dd12(0,0) = threeDtangent(0,1);
  dd12(1,0) = threeDtangent(3,1);

  dd12(0,1) = threeDtangent(0,2);
  dd12(1,1) = threeDtangent(3,2);

  dd12(0,2) = threeDtangent(0,4);
  dd12(1,2) = threeDtangent(3,4);

  dd12(0,3) = threeDtangent(0,5);
  dd12(1,3) = threeDtangent(3,5);


  static Matrix dd22(4,4);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);
  dd22(3,0) = threeDtangent(5,1);
  
  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);
  dd22(3,1) = threeDtangent(5,2);
  
  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);
  dd22(3,2) = threeDtangent(5,4);
  
  dd22(0,3) = threeDtangent(1,5);
  dd22(1,3) = threeDtangent(2,5);
  dd22(2,3) = threeDtangent(4,5);
  dd22(3,3) = threeDtangent(5,5);
  
  static Vector sigma2(4);
  sigma2(0) = threeDstress(1);
  sigma2(1) = threeDstress(2);
  sigma2(2) = threeDstress(4);
  sigma2(3) = threeDstress(5);

  static Vector dd22sigma2(4);
  dd22.Solve(sigma2,dd22sigma2);

  stress.addMatrixVector(1.0, dd12, dd22sigma2, -1.0);

  return stress;
}

int
BeamFiberMaterial2d::commitSensitivity(const Vector &depsdh, int gradIndex,
				       int numGrads)
{
  static Vector dstraindh(6);

  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix dd22(4,4);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);
  dd22(3,0) = threeDtangent(5,1);
  
  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);
  dd22(3,1) = threeDtangent(5,2);
  
  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);
  dd22(3,2) = threeDtangent(5,4);
  
  dd22(0,3) = threeDtangent(1,5);
  dd22(1,3) = threeDtangent(2,5);
  dd22(2,3) = threeDtangent(4,5);
  dd22(3,3) = threeDtangent(5,5);

  static Matrix dd21(4,2);
  dd21(0,0) = threeDtangent(1,0);
  dd21(1,0) = threeDtangent(2,0);
  dd21(2,0) = threeDtangent(4,0);
  dd21(3,0) = threeDtangent(5,0);
  
  dd21(0,1) = threeDtangent(1,3);
  dd21(1,1) = threeDtangent(2,3);
  dd21(2,1) = threeDtangent(4,3);
  dd21(3,1) = threeDtangent(5,3);
  
  static Vector sigma2(4);
  sigma2.addMatrixVector(0.0, dd21, depsdh, -1.0);

  const Vector &threeDstress = theMaterial->getStressSensitivity(gradIndex, true);
  //opserr << threeDstress;
  sigma2(0) -= threeDstress(1);
  sigma2(1) -= threeDstress(2);
  sigma2(2) -= threeDstress(4);
  sigma2(3) -= threeDstress(5);

  //const Vector &threeDstress2 = theMaterial->getStressSensitivity(gradIndex, false);
  //opserr << threeDstress2;
  //sigma2(0) += threeDstress2(1);
  //sigma2(1) += threeDstress2(2);
  //sigma2(2) += threeDstress2(4);
  //sigma2(3) += threeDstress2(5);


  static Vector strain2(4);
  dd22.Solve(sigma2,strain2);


  dstraindh(0) = depsdh(0);
  dstraindh(1) = strain2(0);
  dstraindh(2) = strain2(1);
  dstraindh(3) = depsdh(1);
  dstraindh(4) = strain2(2);
  dstraindh(5) = strain2(3);

  return theMaterial->commitSensitivity(dstraindh, gradIndex, numGrads);
}

const Matrix&  
BeamFiberMaterial2d::getTangent()
{
  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix dd11(2,2);
  dd11(0,0) = threeDtangent(0,0);
  dd11(1,0) = threeDtangent(3,0);

  dd11(0,1) = threeDtangent(0,3);
  dd11(1,1) = threeDtangent(3,3);


  static Matrix dd12(2,4);
  dd12(0,0) = threeDtangent(0,1);
  dd12(1,0) = threeDtangent(3,1);

  dd12(0,1) = threeDtangent(0,2);
  dd12(1,1) = threeDtangent(3,2);

  dd12(0,2) = threeDtangent(0,4);
  dd12(1,2) = threeDtangent(3,4);

  dd12(0,3) = threeDtangent(0,5);
  dd12(1,3) = threeDtangent(3,5);


  static Matrix dd21(4,2);
  dd21(0,0) = threeDtangent(1,0);
  dd21(1,0) = threeDtangent(2,0);
  dd21(2,0) = threeDtangent(4,0);
  dd21(3,0) = threeDtangent(5,0);

  dd21(0,1) = threeDtangent(1,3);
  dd21(1,1) = threeDtangent(2,3);
  dd21(2,1) = threeDtangent(4,3);
  dd21(3,1) = threeDtangent(5,3);


  static Matrix dd22(4,4);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);
  dd22(3,0) = threeDtangent(5,1);
  
  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);
  dd22(3,1) = threeDtangent(5,2);
  
  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);
  dd22(3,2) = threeDtangent(5,4);
  
  dd22(0,3) = threeDtangent(1,5);
  dd22(1,3) = threeDtangent(2,5);
  dd22(2,3) = threeDtangent(4,5);
  dd22(3,3) = threeDtangent(5,5);


  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(4,2);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11;

  return tangent;
}

const Matrix&  
BeamFiberMaterial2d::getInitialTangent()
{
  const Matrix &threeDtangent = theMaterial->getInitialTangent();

  static Matrix dd11(2,2);
  dd11(0,0) = threeDtangent(0,0);
  dd11(1,0) = threeDtangent(3,0);

  dd11(0,1) = threeDtangent(0,3);
  dd11(1,1) = threeDtangent(3,3);


  static Matrix dd12(2,4);
  dd12(0,0) = threeDtangent(0,1);
  dd12(1,0) = threeDtangent(3,1);

  dd12(0,1) = threeDtangent(0,2);
  dd12(1,1) = threeDtangent(3,2);

  dd12(0,2) = threeDtangent(0,4);
  dd12(1,2) = threeDtangent(3,4);

  dd12(0,3) = threeDtangent(0,5);
  dd12(1,3) = threeDtangent(3,5);


  static Matrix dd21(4,2);
  dd21(0,0) = threeDtangent(1,0);
  dd21(1,0) = threeDtangent(2,0);
  dd21(2,0) = threeDtangent(4,0);
  dd21(3,0) = threeDtangent(5,0);

  dd21(0,1) = threeDtangent(1,3);
  dd21(1,1) = threeDtangent(2,3);
  dd21(2,1) = threeDtangent(4,3);
  dd21(3,1) = threeDtangent(5,3);


  static Matrix dd22(4,4);
  dd22(0,0) = threeDtangent(1,1);
  dd22(1,0) = threeDtangent(2,1);
  dd22(2,0) = threeDtangent(4,1);
  dd22(3,0) = threeDtangent(5,1);
  
  dd22(0,1) = threeDtangent(1,2);
  dd22(1,1) = threeDtangent(2,2);
  dd22(2,1) = threeDtangent(4,2);
  dd22(3,1) = threeDtangent(5,2);
  
  dd22(0,2) = threeDtangent(1,4);
  dd22(1,2) = threeDtangent(2,4);
  dd22(2,2) = threeDtangent(4,4);
  dd22(3,2) = threeDtangent(5,4);
  
  dd22(0,3) = threeDtangent(1,5);
  dd22(1,3) = threeDtangent(2,5);
  dd22(2,3) = threeDtangent(4,5);
  dd22(3,3) = threeDtangent(5,5);

  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(4,2);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11;

  return tangent;
}

void  
BeamFiberMaterial2d::Print(OPS_Stream &s, int flag)
{
  s << "BeamFiberMaterial2d, tag: " << this->getTag() << endln;
  s << "\tWrapped material: "<< theMaterial->getTag() << endln;

  theMaterial->Print(s, flag);
}

int 
BeamFiberMaterial2d::sendSelf(int commitTag, Channel &theChannel) 
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
    opserr << "BeamFiberMaterial2d::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(4);
  vecData(0) = Cstrain22;
  vecData(1) = Cstrain33;
  vecData(2) = Cgamma31;
  vecData(3) = Cgamma23;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "BeamFiberMaterial2d::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "BeamFiberMaterial2d::sendSelf() - failed to send vector material\n";

  return res;
}

int 
BeamFiberMaterial2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "BeamFiberMaterial2d::sendSelf() - failed to send id data\n";
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
      opserr << "BeamFiberMaterial2d::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(4);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "BeamFiberMaterial2d::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);
  Cstrain33 = vecData(1);
  Cgamma31  = vecData(2);
  Cgamma23  = vecData(3);

  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma31  = Cgamma31;
  Tgamma23  = Cgamma23;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "BeamFiberMaterial2d::sendSelf() - failed to send vector material\n";
  
  return res;
}

int
BeamFiberMaterial2d::setParameter(const char **argv, int argc,
				  Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}
