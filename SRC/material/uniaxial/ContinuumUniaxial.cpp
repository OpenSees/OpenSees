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

// $Revision: 1.1 $
// $Date: 2007-10-26 04:29:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ContinuumUniaxial.cpp,v $

// Written: MHS
// Created: June 2002
//
// Description: This file contains the class definition of ContinuumUniaxial.
// The ContinuumUniaxial class is a wrapper class that performs static
// condensation on a three-dimensional material model to give the 11
// uniaxial stress component.

#include <ContinuumUniaxial.h>
#include <NDMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

ContinuumUniaxial::ContinuumUniaxial(void):
  UniaxialMaterial(0, MAT_TAG_ContinuumUniaxial), strain11(0.0),
  Tstrain22(0.0),Tstrain33(0.0),Tgamma12(0.0),Tgamma23(0.0),Tgamma31(0.0),
  Cstrain22(0.0),Cstrain33(0.0),Cgamma12(0.0),Cgamma23(0.0),Cgamma31(0.0),
  initialTangent(0.0), theMaterial(0)
{
  // Nothing to do
}

ContinuumUniaxial::ContinuumUniaxial(int tag, NDMaterial &theMat):
  UniaxialMaterial(tag, MAT_TAG_ContinuumUniaxial), strain11(0.0),
  Tstrain22(0.0),Tstrain33(0.0),Tgamma12(0.0),Tgamma23(0.0),Tgamma31(0.0),
  Cstrain22(0.0),Cstrain33(0.0),Cgamma12(0.0),Cgamma23(0.0),Cgamma31(0.0),
  initialTangent(0.0), theMaterial(0)

{
  // Get a copy of the material
  theMaterial = theMat.getCopy("ThreeDimensional");
  
  if (theMaterial == 0)
    opserr << "ContinuumUniaxial::ContinuumUniaxial -- failed to get copy of material" << endln;

  initialTangent = this->getTangent();
}

ContinuumUniaxial::~ContinuumUniaxial(void) 
{ 
  if (theMaterial != 0)
    delete theMaterial;
} 

UniaxialMaterial*
ContinuumUniaxial::getCopy(void) 
{
  ContinuumUniaxial *theCopy =
    new ContinuumUniaxial(this->getTag(), *theMaterial);
  
  theCopy->Tstrain22 = Tstrain22;
  theCopy->Tstrain33 = Tstrain33;
  theCopy->Tgamma12  = Tgamma12;
  theCopy->Tgamma23  = Tgamma23;
  theCopy->Tgamma31  = Tgamma31;

  theCopy->Cstrain22 = Cstrain22;
  theCopy->Cstrain33 = Cstrain33;
  theCopy->Cgamma12  = Cgamma12;
  theCopy->Cgamma23  = Cgamma23;
  theCopy->Cgamma31  = Cgamma31;
  
  return theCopy;
}

int 
ContinuumUniaxial::commitState(void)
{
  Cstrain22 = Tstrain22;
  Cstrain33 = Tstrain33;
  Cgamma12 = Tgamma12;
  Cgamma23 = Tgamma23;
  Cgamma31 = Tgamma31;

  return theMaterial->commitState();
}

int 
ContinuumUniaxial::revertToLastCommit(void)
{
  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma12 = Cgamma12;
  Tgamma23 = Cgamma23;
  Tgamma31 = Cgamma31;
  
  return theMaterial->revertToLastCommit();
}

int
ContinuumUniaxial::revertToStart(void)
{
  Tstrain22 = 0.0;
  Tstrain33 = 0.0;
  Tgamma12  = 0.0;
  Tgamma23  = 0.0;
  Tgamma31  = 0.0;

  Cstrain22 = 0.0;
  Cstrain33 = 0.0;
  Cgamma12  = 0.0;
  Cgamma23  = 0.0;
  Cgamma31  = 0.0;

  return theMaterial->revertToStart();
}

int 
ContinuumUniaxial::setTrialStrain(double strain, double strainRate)
{
  static const double tolerance = 1.0e-08;

  strain11 = strain;

  double norm;
  static Vector condensedStress(5);
  static Vector strainIncrement(5);
  static Vector threeDstrain(6);
  static Matrix dd22(5,5);

  //newton loop to solve for out-of-plane strains
  do {
    //set three dimensional strain
    threeDstrain(0) = strain11;
    threeDstrain(1) = Tstrain22;
    threeDstrain(2) = Tstrain33;
    threeDstrain(3) = Tgamma12;
    threeDstrain(4) = Tgamma23;
    threeDstrain(5) = Tgamma31;

    if (theMaterial->setTrialStrain(threeDstrain) < 0) {
      opserr << "ContinuumUniaxial::setTrialStrain -- setTrialStrain() failed on NDMaterial" << endln;
      return -1;   
    }

    //three dimensional stress
    const Vector &threeDstress = theMaterial->getStress();

    //three dimensional tangent 
    const Matrix &threeDtangent = theMaterial->getTangent();

    //out of plane stress and tangents
    for (int i=0; i<5; i++) {

      condensedStress(i) = threeDstress(i+1);

      for (int j=0; j<5; j++) 
	dd22(i,j) = threeDtangent(i+1,j+1);

    }

    //set norm
    norm = condensedStress.Norm();

    //condensation 
    dd22.Solve(condensedStress, strainIncrement);

    //update out of plane strains
    Tstrain22 -= strainIncrement(0);
    Tstrain33 -= strainIncrement(1);
    Tgamma12  -= strainIncrement(2);
    Tgamma23  -= strainIncrement(3);
    Tgamma31  -= strainIncrement(4);

  } while (norm > tolerance);

  return 0;
}

double
ContinuumUniaxial::getStrain(void)
{
  return strain11;
}

double
ContinuumUniaxial::getStress()
{
  const Vector &threeDstress = theMaterial->getStress();

  return threeDstress(0);
}

double
ContinuumUniaxial::getTangent()
{
  static Matrix dd11(1,1);
  static Matrix dd12(1,5);
  static Matrix dd21(5,1);
  static Matrix dd22(5,5);
  static Matrix dd22invdd21(5,1);

  const Matrix &threeDtangent = theMaterial->getTangent();

  dd11(0,0) = threeDtangent(0,0);

  for (int i=0; i<5; i++) {
    dd12(0,i) = threeDtangent(0,i+1);
    dd21(i,0) = threeDtangent(i+1,0);
    for (int j=0; j<5; j++) {
      dd22(i,j) = threeDtangent(i+1,j+1);
    }
  }

  //condensation 
  dd22.Solve(dd21, dd22invdd21);
  //dd11 -= (dd12*dd22invdd21);
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);

  return dd11(0,0);
}

double
ContinuumUniaxial::getInitialTangent(void)
{
  return initialTangent;
}

void  
ContinuumUniaxial::Print(OPS_Stream &s, int flag)
{
  s << "ContinuumUniaxial, tag: " << this->getTag() << endln;
  s << "\tWrapped material: "<< theMaterial->getTag() << endln;
  
  theMaterial->Print(s, flag);
}

int 
ContinuumUniaxial::sendSelf(int commitTag, Channel &theChannel) 
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
    opserr << "ContinuumUniaxial::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(5);
  vecData(0) = Cstrain22;
  vecData(1) = Cstrain33;
  vecData(2) = Cgamma12;
  vecData(3) = Cgamma23;
  vecData(4) = Cgamma31;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "ContinuumUniaxial::sendSelf() - failed to send vector data" << endln;
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "ContinuumUniaxial::sendSelf() - failed to send vector material" << endln;

  return res;
}

int 
ContinuumUniaxial::recvSelf(int commitTag, Channel &theChannel,
			     FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containg the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "ContinuumUniaxial::sendSelf() - failed to send id data" << endln;
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
      opserr << "ContinuumUniaxial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(5);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "ContinuumUniaxial::sendSelf() - failed to send vector data" << endln;
    return res;
  }

  Cstrain22 = vecData(0);
  Cstrain33 = vecData(1);
  Cgamma12  = vecData(2);
  Cgamma23  = vecData(3);
  Cgamma31  = vecData(4);

  Tstrain22 = Cstrain22;
  Tstrain33 = Cstrain33;
  Tgamma12  = Cgamma12;
  Tgamma23  = Cgamma23;
  Tgamma31  = Cgamma31;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "ContinuumUniaxial::sendSelf() - failed to send vector material" << endln;
  
  return res;
}
