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
// Created: June 2019
//
// Description: This file contains the class definition of BeamFiberMaterial2dPS.
// The BeamFiberMaterial2dPS class is a wrapper class that performs static
// condensation on a plane stress material model to give the 11 and 12
// stress components which can then be integrated over an area to model a
// shear flexible 2D beam.


#include <BeamFiberMaterial2dPS.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <elementAPI.h>

Vector BeamFiberMaterial2dPS::stress(2);
Matrix BeamFiberMaterial2dPS::tangent(2,2);

void* OPS_BeamFiberMaterial2dPS()
{
    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial BeamFiber2dPS tag? matTag?" << endln;
	return 0;
    }

    int tags[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, tags) < 0) {
	opserr << "WARNING invalid nDMaterial BeamFiber2dPS tag or matTag" << endln;
	return 0;
    }

    int tag = tags[0];
    int matTag = tags[1];

    NDMaterial *psMaterial = OPS_getNDMaterial(matTag);
    if (psMaterial == 0) {
	opserr << "WARNING nD material does not exist\n";
	opserr << "nD material: " << matTag;
	opserr << "\nBeamFiber2d nDMaterial: " << tag << endln;
	return 0;
    }

    return new BeamFiberMaterial2dPS(tag, *psMaterial);
}

//      0  1  2  3  4  5
// ND: 11 22 33 12 23 31
// BF: 11 12 22 33 23 31

BeamFiberMaterial2dPS::BeamFiberMaterial2dPS(void)
  :NDMaterial(0, ND_TAG_BeamFiberMaterial2dPS),
   Tstrain22(0.0), Cstrain22(0.0), theMaterial(0), strain(2)
{
	// Nothing to do
}

BeamFiberMaterial2dPS::BeamFiberMaterial2dPS(int tag, NDMaterial &theMat)
  :NDMaterial(tag, ND_TAG_BeamFiberMaterial2dPS),
   Tstrain22(0.0), Cstrain22(0.0), theMaterial(0), strain(2)

{
  // Get a copy of the material
  theMaterial = theMat.getCopy("PlaneStress");
  
  if (theMaterial == 0) {
    opserr << "BeamFiberMaterial2dPS::BeamFiberMaterial2dPS -- failed to get copy of material\n";
    exit(-1);
  }
}

BeamFiberMaterial2dPS::~BeamFiberMaterial2dPS(void) 
{ 
  if (theMaterial != 0)
    delete theMaterial;
} 

NDMaterial*
BeamFiberMaterial2dPS::getCopy(void) 
{
  BeamFiberMaterial2dPS *theCopy =
    new BeamFiberMaterial2dPS(this->getTag(), *theMaterial);
  
  theCopy->Tstrain22 = this->Tstrain22;
  theCopy->Cstrain22 = this->Cstrain22;
  
  return theCopy;
}

NDMaterial* 
BeamFiberMaterial2dPS::getCopy(const char *type)
{
  if (strcmp(type, "BeamFiber2d") == 0)
    return this->getCopy();
  else
    return 0;
}

int 
BeamFiberMaterial2dPS::getOrder(void) const
{
  return 2;
}

const char*
BeamFiberMaterial2dPS::getType(void) const 
{
  return "BeamFiber2d";
}

int 
BeamFiberMaterial2dPS::commitState(void)
{
  Cstrain22 = Tstrain22;

  return theMaterial->commitState();
}

int 
BeamFiberMaterial2dPS::revertToLastCommit(void)
{
  Tstrain22 = Cstrain22;
  
  return theMaterial->revertToLastCommit();
}

int
BeamFiberMaterial2dPS::revertToStart()
{
  this->Tstrain22 = 0.0;
  this->Cstrain22 = 0.0;

  strain.Zero();

  return theMaterial->revertToStart();
}

double
BeamFiberMaterial2dPS::getRho(void)
{
  return theMaterial->getRho();
}


//receive the strain
int 
BeamFiberMaterial2dPS::setTrialStrain(const Vector &strainFromElement)
{
  static const double tolerance = 1.0e-12;

  strain(0) = strainFromElement(0);
  strain(1) = strainFromElement(1);

  //newton loop to solve for out-of-plane strains

  double norm;
  static Vector condensedStress(1);
  static Vector strainIncrement(1);
  static Vector PSstrain(3);
  static Matrix dd22(1,1);

  int count = 0;
  const int maxCount = 20;
  double norm0;

  do {

    //set three dimensional strain
    PSstrain(0) = this->strain(0); // 11
    PSstrain(1) = this->Tstrain22; // 22
    PSstrain(2) = this->strain(1); // 12

    if (theMaterial->setTrialStrain(PSstrain) < 0) {
      opserr << "BeamFiberMaterial2dPS::setTrialStrain - setStrain failed in material with strain " << PSstrain;
      return -1;   
    }

    //three dimensional stress
    const Vector &PSstress = theMaterial->getStress();

    //three dimensional tangent 
    const Matrix &PStangent = theMaterial->getTangent();

    condensedStress(0) = PSstress(1);

    dd22(0,0) = PStangent(1,1);

    //set norm
    norm = condensedStress.Norm();
    if (count == 0)
      norm0 = norm;

    //condensation 
    dd22.Solve(condensedStress, strainIncrement);

    //update out of plane strains
    this->Tstrain22 -= strainIncrement(0);

  } while (count++ < maxCount && norm > tolerance*norm0);

  return 0;
}

const Vector& 
BeamFiberMaterial2dPS::getStrain(void)
{
  return strain;
}

const Vector&  
BeamFiberMaterial2dPS::getStress()
{
  const Vector &PSstress = theMaterial->getStress();

  stress(0) = PSstress(0);
  stress(1) = PSstress(2);

  return stress;
}

const Vector& 
BeamFiberMaterial2dPS::getStressSensitivity(int gradIndex,
					  bool conditional)
{
  const Vector &PSstress = theMaterial->getStressSensitivity(gradIndex, conditional);

  stress(0) = PSstress(0);
  stress(1) = PSstress(2);

  const Matrix &PStangent = theMaterial->getTangent();

  static Matrix dd12(2,1);
  dd12(0,0) = PStangent(0,1);
  dd12(1,0) = PStangent(2,1);

  static Matrix dd22(1,1);
  dd22(0,0) = PStangent(1,1);
  
  static Vector sigma2(1);
  sigma2(0) = PSstress(1);

  static Vector dd22sigma2(1);
  dd22.Solve(sigma2,dd22sigma2);

  stress.addMatrixVector(1.0, dd12, dd22sigma2, -1.0);

  return stress;
}

int
BeamFiberMaterial2dPS::commitSensitivity(const Vector &depsdh, int gradIndex,
				       int numGrads)
{
  static Vector dstraindh(6);

  const Matrix &PStangent = theMaterial->getTangent();

  static Matrix dd22(1,1);
  dd22(0,0) = PStangent(1,1);

  static Matrix dd21(1,2);
  dd21(0,0) = PStangent(1,0);
  dd21(0,1) = PStangent(1,2);
  
  static Vector sigma2(1);
  sigma2.addMatrixVector(0.0, dd21, depsdh, -1.0);

  const Vector &PSstress = theMaterial->getStressSensitivity(gradIndex, true);
  //opserr << PSstress;
  sigma2(0) -= PSstress(1);

  //const Vector &PSstress2 = theMaterial->getStressSensitivity(gradIndex, false);
  //opserr << PSstress2;
  //sigma2(0) += PSstress2(1);

  static Vector strain2(1);
  dd22.Solve(sigma2,strain2);

  dstraindh(0) = depsdh(0);
  dstraindh(1) = strain2(0);
  dstraindh(2) = depsdh(1);

  return theMaterial->commitSensitivity(dstraindh, gradIndex, numGrads);
}

const Matrix&  
BeamFiberMaterial2dPS::getTangent()
{
  const Matrix &PStangent = theMaterial->getTangent();

  static Matrix dd11(2,2);
  dd11(0,0) = PStangent(0,0);
  dd11(1,0) = PStangent(2,0);

  dd11(0,1) = PStangent(0,2);
  dd11(1,1) = PStangent(2,2);

  static Matrix dd12(2,1);
  dd12(0,0) = PStangent(0,1);
  dd12(1,0) = PStangent(2,1);

  static Matrix dd21(1,2);
  dd21(0,0) = PStangent(1,0);
  dd21(0,1) = PStangent(1,2);

  static Matrix dd22(1,1);
  dd22(0,0) = PStangent(1,1);

  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(1,2);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11;

  return tangent;
}

const Matrix&  
BeamFiberMaterial2dPS::getInitialTangent()
{
  const Matrix &PStangent = theMaterial->getInitialTangent();

  static Matrix dd11(2,2);
  dd11(0,0) = PStangent(0,0);
  dd11(1,0) = PStangent(2,0);

  dd11(0,1) = PStangent(0,2);
  dd11(1,1) = PStangent(2,2);

  static Matrix dd12(2,1);
  dd12(0,0) = PStangent(0,1);
  dd12(1,0) = PStangent(2,1);

  static Matrix dd21(1,2);
  dd21(0,0) = PStangent(1,0);
  dd21(0,1) = PStangent(1,2);

  static Matrix dd22(1,1);
  dd22(0,0) = PStangent(1,1);

  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(1,2);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11; 
  //this->tangent  -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11;

  return tangent;
}

void  
BeamFiberMaterial2dPS::Print(OPS_Stream &s, int flag)
{
  s << "BeamFiberMaterial2dPS, tag: " << this->getTag() << endln;
  s << "\tWrapped material: "<< theMaterial->getTag() << endln;

  theMaterial->Print(s, flag);
}

int 
BeamFiberMaterial2dPS::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  // put tag and assocaited materials class and database tags into an id and send it
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
    opserr << "BeamFiberMaterial2dPS::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(1);
  vecData(0) = Cstrain22;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "BeamFiberMaterial2dPS::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "BeamFiberMaterial2dPS::sendSelf() - failed to send vector material\n";

  return res;
}

int 
BeamFiberMaterial2dPS::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containg the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "BeamFiberMaterial2dPS::sendSelf() - failed to send id data\n";
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
      opserr << "BeamFiberMaterial2dPS::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(1);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "BeamFiberMaterial2dPS::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);
  Tstrain22 = Cstrain22;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "BeamFiberMaterial2dPS::sendSelf() - failed to send vector material\n";
  
  return res;
}

int
BeamFiberMaterial2dPS::setParameter(const char **argv, int argc,
				  Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}
