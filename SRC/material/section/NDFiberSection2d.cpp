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
// Created: 2012
//
// Description: This file contains the class implementation of FiberSection2d.

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <NDFiberSection2d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <NDMaterial.h>
#include <SectionIntegration.h>
#include <Parameter.h>
#include <elementAPI.h>

ID NDFiberSection2d::code(3);
Matrix NDFiberSection2d::fs(3,3);

void* OPS_NDFiberSection2d()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	opserr<<"insufficient arguments for NDFiberSection2d\n";
	return 0;
    }

    numData = 1;
    int tag;
    if (OPS_GetIntInput(&numData,&tag) < 0) return 0;

    bool computeCentroid = true;
    if (OPS_GetNumRemainingInputArgs() > 0) {
      const char* opt = OPS_GetString();
      if (strcmp(opt, "-noCentroid") == 0)
	computeCentroid = false;
    }
    
    int num = 30;
    return new NDFiberSection2d(tag, num, computeCentroid);
}

// constructors:
NDFiberSection2d::NDFiberSection2d(int tag, int num, Fiber **fibers, double a, bool compCentroid): 
  SectionForceDeformation(tag, SEC_TAG_NDFiberSection2d),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), Abar(0.0), yBar(0.0), computeCentroid(compCentroid),
  alpha(a), sectionIntegr(0), e(3), s(0), ks(0), 
  parameterID(0), dedh(3)
{
  if (numFibers != 0) {
    theMaterials = new NDMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "NDFiberSection2d::NDFiberSection2d -- failed to allocate Material pointers";
      exit(-1);
    }

    matData = new double [numFibers*2];

    if (matData == 0) {
      opserr << "NDFiberSection2d::NDFiberSection2d -- failed to allocate double array for material data\n";
      exit(-1);
    }


    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();
      Abar  += Area;
      QzBar += yLoc*Area;
      matData[i*2] = yLoc;
      matData[i*2+1] = Area;
      NDMaterial *theMat = theFiber->getNDMaterial();
      theMaterials[i] = theMat->getCopy("BeamFiber2d");

      if (theMaterials[i] == 0) {
	opserr << "NDFiberSection2d::NDFiberSection2d -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    

    if (computeCentroid)
      yBar = QzBar/Abar;  
  }

  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;
  kData[4] = 0.0;
  kData[5] = 0.0;
  kData[6] = 0.0;
  kData[7] = 0.0;
  kData[8] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_VY;
}

NDFiberSection2d::NDFiberSection2d(int tag, int num, double a, bool compCentroid): 
    SectionForceDeformation(tag, SEC_TAG_NDFiberSection2d),
    numFibers(0), sizeFibers(num), theMaterials(0), matData(0),
    QzBar(0.0), Abar(0.0), yBar(0.0), computeCentroid(compCentroid),
    alpha(a), sectionIntegr(0), e(3), s(0), ks(0), 
    parameterID(0), dedh(3)
{
    if (sizeFibers != 0) {
	theMaterials = new NDMaterial *[sizeFibers];

	if (theMaterials == 0) {
	    opserr << "NDFiberSection2d::NDFiberSection2d -- failed to allocate Material pointers";
	    exit(-1);
	}

	matData = new double [sizeFibers*2];

	if (matData == 0) {
	    opserr << "NDFiberSection2d::NDFiberSection2d -- failed to allocate double array for material data\n";
	    exit(-1);
	}


	for (int i = 0; i < sizeFibers; i++) {
	    matData[i*2] = 0.0;
	    matData[i*2+1] = 0.0;
	    theMaterials[i] = 0;
	}    
    }

    s = new Vector(sData, 3);
    ks = new Matrix(kData, 3, 3);

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;

    kData[0] = 0.0;
    kData[1] = 0.0;
    kData[2] = 0.0;
    kData[3] = 0.0;
    kData[4] = 0.0;
    kData[5] = 0.0;
    kData[6] = 0.0;
    kData[7] = 0.0;
    kData[8] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_VY;
}

NDFiberSection2d::NDFiberSection2d(int tag, int num, NDMaterial **mats,
				   SectionIntegration &si, double a, bool compCentroid):
  SectionForceDeformation(tag, SEC_TAG_NDFiberSection2d),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), Abar(0.0), yBar(0.0), computeCentroid(compCentroid),
  alpha(a), sectionIntegr(0), e(3), s(0), ks(0), 
  parameterID(0), dedh(3)
{
  if (numFibers != 0) {
    theMaterials = new NDMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "NDFiberSection2d::NDFiberSection2d -- failed to allocate Material pointers";
      exit(-1);
    }
    matData = new double [numFibers*2];

    if (matData == 0) {
      opserr << "NDFiberSection2d::NDFiberSection2d -- failed to allocate double array for material data\n";
      exit(-1);
    }
  }

  sectionIntegr = si.getCopy();
  if (sectionIntegr == 0) {
    opserr << "Error: NDFiberSection2d::NDFiberSection2d: could not create copy of section integration object" << endln;
    exit(-1);
  }

  static double fiberLocs[10000];
  sectionIntegr->getFiberLocations(numFibers, fiberLocs);
  
  static double fiberArea[10000];
  sectionIntegr->getFiberWeights(numFibers, fiberArea);

  for (int i = 0; i < numFibers; i++) {

    Abar  += fiberArea[i];
    QzBar += fiberLocs[i]*fiberArea[i];

    theMaterials[i] = mats[i]->getCopy("BeamFiber2d");
    
    if (theMaterials[i] == 0) {
      opserr << "NDFiberSection2d::NDFiberSection2d -- failed to get copy of a Material\n";
      exit(-1);
    }
  }    

  if (computeCentroid)
    yBar = QzBar/Abar;  

  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);
  
  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  
  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;
  kData[4] = 0.0;
  kData[5] = 0.0;
  kData[6] = 0.0;
  kData[7] = 0.0;
  kData[8] = 0.0;
  
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_VY;
}

// constructor for blank object that recvSelf needs to be invoked upon
NDFiberSection2d::NDFiberSection2d():
  SectionForceDeformation(0, SEC_TAG_NDFiberSection2d),
  numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
  QzBar(0.0), Abar(0.0), yBar(0.0), computeCentroid(true),
  alpha(1.0), sectionIntegr(0), 
  e(3), s(0), ks(0), parameterID(0), dedh(3)
{
  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;
  kData[4] = 0.0;
  kData[5] = 0.0;
  kData[6] = 0.0;
  kData[7] = 0.0;
  kData[8] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_VY;
}

int
NDFiberSection2d::addFiber(Fiber &newFiber)
{
  // need to create larger arrays
  if(numFibers == sizeFibers) {
      int newSize = 2*sizeFibers;
      NDMaterial **newArray = new NDMaterial *[newSize]; 
      double *newMatData = new double [2 * newSize];
      if (newArray == 0 || newMatData == 0) {
	  opserr <<"NDFiberSection2d::addFiber -- failed to allocate Fiber pointers\n";
	  return -1;
      }
      
      // copy the old pointers and data
      for (int i = 0; i < numFibers; i++) {
	  newArray[i] = theMaterials[i];
	  newMatData[2*i] = matData[2*i];
	  newMatData[2*i+1] = matData[2*i+1];
      }

      // initial memory
      for (int i = numFibers; i < newSize; i++) {
	  newArray[i] = 0;
	  newMatData[2*i] = 0.0;
	  newMatData[2*i+1] = 0.0;
      }
      sizeFibers = newSize;

      // set new memory
      if (theMaterials != 0) {
	  delete [] theMaterials;
	  delete [] matData;
      }

      theMaterials = newArray;
      matData = newMatData;
  }

  // set the new pointers and data
  double yLoc, zLoc, Area;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  matData[numFibers*2] = yLoc;
  matData[numFibers*2+1] = Area;
  NDMaterial *theMat = newFiber.getNDMaterial();
  theMaterials[numFibers] = theMat->getCopy("BeamFiber2d");

  if (theMaterials[numFibers] == 0) {
    opserr <<"NDFiberSection2d::addFiber -- failed to get copy of a Material\n";
    return -1;
  }

  numFibers++;

  // Recompute centroid
  if (computeCentroid) {
    Abar  += Area;
    QzBar += yLoc*Area;
    yBar = QzBar/Abar;
  }

  return 0;
}


// destructor:
NDFiberSection2d::~NDFiberSection2d()
{
  if (theMaterials != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];
      
    delete [] theMaterials;
  }

  if (matData != 0)
    delete [] matData;

  if (s != 0)
    delete s;

  if (ks != 0)
    delete ks;

  if (sectionIntegr != 0)
    delete sectionIntegr;
}

int
NDFiberSection2d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;

  e = deforms;

  kData[0] = 0.0; 
  kData[1] = 0.0; 
  kData[2] = 0.0; 
  kData[3] = 0.0;
  kData[4] = 0.0;
  kData[5] = 0.0;
  kData[6] = 0.0;
  kData[7] = 0.0;
  kData[8] = 0.0;
  
  sData[0] = 0.0; 
  sData[1] = 0.0;
  sData[2] = 0.0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);

  static double fiberLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      fiberLocs[i] = matData[2*i];
      fiberArea[i] = matData[2*i+1];
    }
  }

  static Vector eps(2);

  double rootAlpha = 1.0;
  eps(1) = d2;
  if (alpha != 1.0) {
    rootAlpha = sqrt(alpha);
    eps(1) *= rootAlpha;
  }

  for (int i = 0; i < numFibers; i++) {
    NDMaterial *theMat = theMaterials[i];
    double y = fiberLocs[i] - yBar;
    double A = fiberArea[i];

    // determine material strain and set it
    eps(0) = d0 - y*d1;

    res += theMat->setTrialStrain(eps);

    const Vector &stress = theMat->getStress();
    const Matrix &tangent = theMat->getTangent();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;

    double ks1 = d00*-y;
    kData[0] += d00;
    kData[1] += ks1;
    kData[4] += ks1*-y;

    kData[2] += d10;
    kData[6] += d01;

    kData[5] += d10*-y;
    kData[7] += d01*-y;

    kData[8] += d11;

    double fs0 = stress(0)*A;
    double fs1 = stress(1)*A;
    sData[0] += fs0;
    sData[1] += fs0*-y;
    sData[2] += fs1;
  }

  kData[3] = kData[1];

  if (alpha != 1.0) {
    sData[2] *= rootAlpha;
    
    kData[2] *= rootAlpha;
    kData[6] *= rootAlpha;
    
    kData[5] *= rootAlpha;
    kData[7] *= rootAlpha;
    
    kData[8] *= alpha;
  }

  return res;
}

const Vector&
NDFiberSection2d::getSectionDeformation(void)
{
  return e;
}

const Matrix&
NDFiberSection2d::getInitialTangent(void)
{
  static double kInitial[9];
  static Matrix kInitialMatrix(kInitial, 3, 3);
  kInitial[0] = 0.0; 
  kInitial[1] = 0.0; 
  kInitial[2] = 0.0; 
  kInitial[3] = 0.0;
  kInitial[4] = 0.0;
  kInitial[5] = 0.0;
  kInitial[6] = 0.0;
  kInitial[7] = 0.0;
  kInitial[8] = 0.0;

  static double fiberLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      fiberLocs[i] = matData[2*i];
      fiberArea[i] = matData[2*i+1];
    }
  }

  for (int i = 0; i < numFibers; i++) {
    NDMaterial *theMat = theMaterials[i];
    double y = fiberLocs[i] - yBar;
    double A = fiberArea[i];

    const Matrix &tangent = theMat->getInitialTangent();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;

    double ks1 = d00 * -y;
    kInitial[0] += d00;
    kInitial[1] += ks1;
    kInitial[4] += ks1*-y;

    kInitial[2] += d10;
    kInitial[6] += d01;

    kInitial[5] += d10*-y;
    kInitial[7] += d01*-y;

    kInitial[8] += d11;
  }

  kInitial[3] = kInitial[1];

  if (alpha != 1.0) {
    double rootAlpha = sqrt(alpha);

    kInitial[2] *= rootAlpha;
    kInitial[6] *= rootAlpha;
    
    kInitial[5] *= rootAlpha;
    kInitial[7] *= rootAlpha;
    
    kInitial[8] *= alpha;
  }

  return kInitialMatrix;
}

const Matrix&
NDFiberSection2d::getSectionTangent(void)
{
  //opserr << *ks << endln;
  return *ks;
}

const Vector&
NDFiberSection2d::getStressResultant(void)
{
  return *s;
}

SectionForceDeformation*
NDFiberSection2d::getCopy(void)
{
  NDFiberSection2d *theCopy = new NDFiberSection2d ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;
  theCopy->sizeFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new NDMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr <<"NDFiberSection2d::getCopy -- failed to allocate Material pointers\n";
      exit(-1);
    }
  
    theCopy->matData = new double [numFibers*2];

    if (theCopy->matData == 0) {
      opserr << "NDFiberSection2d::getCopy -- failed to allocate double array for material data\n";
      exit(-1);
    }
			    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*2] = matData[i*2];
      theCopy->matData[i*2+1] = matData[i*2+1];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy("BeamFiber2d");

      if (theCopy->theMaterials[i] == 0) {
	opserr <<"NDFiberSection2d::getCopy -- failed to get copy of a Material";
	exit(-1);
      }
    }  
  }

  theCopy->e = e;
  theCopy->QzBar = QzBar;
  theCopy->Abar = Abar;
  theCopy->yBar = yBar;
  theCopy->computeCentroid = computeCentroid;
  theCopy->alpha = alpha;
  theCopy->parameterID = parameterID;

  theCopy->kData[0] = kData[0];
  theCopy->kData[1] = kData[1];
  theCopy->kData[2] = kData[2];
  theCopy->kData[3] = kData[3];
  theCopy->kData[4] = kData[4];
  theCopy->kData[5] = kData[5];
  theCopy->kData[6] = kData[6];
  theCopy->kData[7] = kData[7];
  theCopy->kData[8] = kData[8];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];

  if (sectionIntegr != 0)
    theCopy->sectionIntegr = sectionIntegr->getCopy();
  else
    theCopy->sectionIntegr = 0;

  return theCopy;
}

const ID&
NDFiberSection2d::getType ()
{
  return code;
}

int
NDFiberSection2d::getOrder () const
{
  return 3;
}

int
NDFiberSection2d::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  return err;
}

int
NDFiberSection2d::revertToLastCommit(void)
{
  int err = 0;

  kData[0] = 0.0; 
  kData[1] = 0.0; 
  kData[2] = 0.0; 
  kData[3] = 0.0;
  kData[4] = 0.0;
  kData[5] = 0.0;
  kData[6] = 0.0;
  kData[7] = 0.0;
  kData[8] = 0.0;

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  
  static double fiberLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      fiberLocs[i] = matData[2*i];
      fiberArea[i] = matData[2*i+1];
    }
  }

  for (int i = 0; i < numFibers; i++) {
    NDMaterial *theMat = theMaterials[i];
    double y = fiberLocs[i] - yBar;
    double A = fiberArea[i];

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    // get material stress & tangent for this strain and determine ks and fs
    const Matrix &tangent = theMat->getTangent();
    const Vector &stress = theMat->getStress();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;

    double ks1 = d00*-y;
    kData[0] += d00;
    kData[1] += ks1;
    kData[4] += ks1*-y;

    kData[2] += d10;
    kData[6] += d01;

    kData[5] += d10*-y;
    kData[7] += d01*-y;

    kData[8] += d11;

    double fs0 = stress(0)*A;
    double fs1 = stress(1)*A;
    sData[0] += fs0;
    sData[1] += fs0*-y;
    sData[2] += fs1;
  }

  kData[3] = kData[1];

  if (alpha != 1.0) {
    double rootAlpha = sqrt(alpha);

    sData[2] *= rootAlpha;
    
    kData[2] *= rootAlpha;
    kData[6] *= rootAlpha;
    
    kData[5] *= rootAlpha;
    kData[7] *= rootAlpha;
    
    kData[8] *= alpha;
  }

  return err;
}

int
NDFiberSection2d::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;

  e.Zero();

  kData[0] = 0.0; 
  kData[1] = 0.0; 
  kData[2] = 0.0; 
  kData[3] = 0.0;
  kData[4] = 0.0;
  kData[5] = 0.0;
  kData[6] = 0.0;
  kData[7] = 0.0;
  kData[8] = 0.0;
  
  sData[0] = 0.0; 
  sData[1] = 0.0;
  sData[2] = 0.0;
  
  static double fiberLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      fiberLocs[i] = matData[2*i];
      fiberArea[i] = matData[2*i+1];
    }
  }

  for (int i = 0; i < numFibers; i++) {
    NDMaterial *theMat = theMaterials[i];
    double y = fiberLocs[i] - yBar;
    double A = fiberArea[i];

    // invoke revertToLast on the material
    err += theMat->revertToStart();

    // get material stress & tangent for this strain and determine ks and fs
    const Matrix &tangent = theMat->getTangent();
    const Vector &stress = theMat->getStress();

    double d00 = tangent(0,0)*A;
    double d01 = tangent(0,1)*A;
    double d10 = tangent(1,0)*A;
    double d11 = tangent(1,1)*A;

    double ks1 = d00*-y;
    kData[0] += d00;
    kData[1] += ks1;
    kData[4] += ks1*-y;

    kData[2] += d10;
    kData[6] += d01;

    kData[5] += d10*-y;
    kData[7] += d01*-y;

    kData[8] += d11;

    double fs0 = stress(0)*A;
    double fs1 = stress(1)*A;
    sData[0] += fs0;
    sData[1] += fs0*-y;
    sData[2] += fs1;
  }

  kData[3] = kData[1];

  if (alpha != 1.0) {
    double rootAlpha = sqrt(alpha);

    sData[2] *= rootAlpha;
    
    kData[2] *= rootAlpha;
    kData[6] *= rootAlpha;
    
    kData[5] *= rootAlpha;
    kData[7] *= rootAlpha;
    
    kData[8] *= alpha;
  }

  return err;
}

int
NDFiberSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 3 so no conflict with matData below if just 1 fiber
  static ID data(3);
  data(0) = this->getTag();
  data(1) = numFibers;
  data(2) = computeCentroid ? 1 : 0; // Now the ID data is really 3  
  int dbTag = this->getDbTag();
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "NDFiberSection2d::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (numFibers != 0) {
    
    // create an id containingg classTag and dbTag for each material & send it
    ID materialData(2*numFibers);
    for (int i=0; i<numFibers; i++) {
      NDMaterial *theMat = theMaterials[i];
      materialData(2*i) = theMat->getClassTag();
      int matDbTag = theMat->getDbTag();
      if (matDbTag == 0) {
	matDbTag = theChannel.getDbTag();
	if (matDbTag != 0)
	  theMat->setDbTag(matDbTag);
      }
      materialData(2*i+1) = matDbTag;
    }    
    
    res += theChannel.sendID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "NDFiberSection2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 2*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "NDFiberSection2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);

  }

  return res;
}

int
NDFiberSection2d::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "NDFiberSection2d::recvSelf - failed to recv ID data\n";
    return res;
  }    
  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "NDFiberSection2d::recvSelf - failed to recv material data\n";
      return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != data(1)) {
      // delete old stuff if outa date
      if (theMaterials != 0) {
	for (int i=0; i<numFibers; i++)
	  delete theMaterials[i];
	delete [] theMaterials;
	if (matData != 0)
	  delete [] matData;
	matData = 0;
	theMaterials = 0;
      }

      // create memory to hold material pointers and fiber data
      numFibers = data(1);
      sizeFibers = data(1);
      if (numFibers != 0) {
	theMaterials = new NDMaterial *[numFibers];
	
	if (theMaterials == 0) {
	  opserr <<"NDFiberSection2d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}
	
	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;

	matData = new double [numFibers*2];

	if (matData == 0) {
	  opserr <<"NDFiberSection2d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, 2*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "NDFiberSection2d::recvSelf - failed to recv material data\n";
      return res;
    }    

    int i;
    for (i=0; i<numFibers; i++) {
      int classTag = materialData(2*i);
      int dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of corrcet type, 
      // release old and create a new one
      if (theMaterials[i] == 0)
	theMaterials[i] = theBroker.getNewNDMaterial(classTag);
      else if (theMaterials[i]->getClassTag() != classTag) {
	delete theMaterials[i];
	theMaterials[i] = theBroker.getNewNDMaterial(classTag);      
      }

      if (theMaterials[i] == 0) {
	opserr <<"NDFiberSection2d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    QzBar = 0.0;
    Abar  = 0.0;
    double yLoc, Area;

    computeCentroid = data(2) ? true : false;
    
    // Recompute centroid
    for (i = 0; computeCentroid && i < numFibers; i++) {
      yLoc = matData[2*i];
      Area = matData[2*i+1];
      Abar  += Area;
      QzBar += yLoc*Area;
    }

    if (computeCentroid)
      yBar = QzBar/Abar;
    else
      yBar = 0.0;
  }    

  return res;
}

void
NDFiberSection2d::Print(OPS_Stream &s, int flag)
{
  s << "\nNDFiberSection2d, tag: " << this->getTag() << endln;
  s << "\tSection code: " << code;
  s << "\tNumber of Fibers: " << numFibers << endln;
  s << "\tCentroid: " << yBar << endln;
  s << "\tShape factor, alpha = " << alpha << endln;

  if (flag == 1) {
    for (int i = 0; i < numFibers; i++) {
      s << "\nLocation (y) = (" << matData[2*i] << ")";
      s << "\nArea = " << matData[2*i+1] << endln;
      theMaterials[i]->Print(s, flag);
    }
  }
}

Response*
NDFiberSection2d::setResponse(const char **argv, int argc,
			      OPS_Stream &output)
{
  Response *theResponse =0;

  if (argc > 2 || strcmp(argv[0],"fiber") == 0) {
    
    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 3) {		  // fiber number was input directly
      
      key = atoi(argv[1]);
      
    } else if (argc > 4) {  // find fiber closest to coord. with mat tag
      
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      
      double closestDist = 0;
      double ySearch, dy;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < numFibers; j++) {
	if (matTag == theMaterials[j]->getTag()) {
	  ySearch = matData[2*j];
	  dy = ySearch-yCoord;
	  closestDist = fabs(dy);
	  key = j;
	  break;
	}
      }
      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
	if (matTag == theMaterials[j]->getTag()) {
	  ySearch = matData[2*j];
	  dy = ySearch-yCoord;
	  distance = fabs(dy);
	  if (distance < closestDist) {
	    closestDist = distance;
	    key = j;
	  }
	}
      }
      passarg = 4;
    }
    
    else {                  // fiber near-to coordinate specified
      
      double yCoord = atof(argv[1]);
      double closestDist;
      double ySearch, dy;
      double distance;
      
      ySearch = matData[0];
      dy = ySearch-yCoord;
      closestDist = fabs(dy);
      key = 0;
      for (int j = 1; j < numFibers; j++) {
	ySearch = matData[2*j];
	dy = ySearch-yCoord;
	
	distance = fabs(dy);
	if (distance < closestDist) {
	  closestDist = distance;
	  key = j;
	}
      }
      passarg = 3;
    }
    
    if (key < numFibers && key >= 0) {
      output.tag("FiberOutput");
      output.attr("yLoc",matData[2*key]);
      output.attr("zLoc",0.0);
      output.attr("area",matData[2*key+1]);
      
      theResponse =  theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }

    return theResponse;
  }

  // If not a fiber response, call the base class method
  return SectionForceDeformation::setResponse(argv, argc, output);
}


int 
NDFiberSection2d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}



// AddingSensitivity:BEGIN ////////////////////////////////////
int
NDFiberSection2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  if (strstr(argv[0],"alpha") != 0)
    return param.addObject(1, this);

  // Check if the parameter belongs to the material (only option for now)
  if (strstr(argv[0],"material") != 0) {
    
    if (argc < 3)
      return 0;

    // Get the tag of the material
    int materialTag = atoi(argv[1]);
    
    // Loop over fibers to find the right material
    for (int i = 0; i < numFibers; i++)
      if (materialTag == theMaterials[i]->getTag()) {
	int ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
	if (ok != -1)
	  result = ok;
      }
    return result;
  }

  // Check if it belongs to the section integration
  else if (strstr(argv[0],"integration") != 0) {
    if (sectionIntegr != 0)
      return sectionIntegr->setParameter(&argv[1], argc-1, param);
    else
      return -1;
  }

  int ok = 0;
  
  for (int i = 0; i < numFibers; i++) {
    ok = theMaterials[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  if (sectionIntegr != 0) {
    ok = sectionIntegr->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  return result;
}

int
NDFiberSection2d::updateParameter(int paramID, Information &info)
{
  switch(paramID) {
  case 1:
    alpha = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
NDFiberSection2d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector &
NDFiberSection2d::getSectionDeformationSensitivity(int gradIndex)
{
  return dedh;
}

const Vector &
NDFiberSection2d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(3);
  
  ds.Zero();
  
  double y, A;
  static Vector stress(2);
  static Vector dsigdh(2);
  static Vector sig_dAdh(2);
  static Matrix tangent(2,2);

  static double fiberLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      fiberLocs[i] = matData[2*i];
      fiberArea[i] = matData[2*i+1];
    }
  }

  static double locsDeriv[10000];
  static double areaDeriv[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getLocationsDeriv(numFibers, locsDeriv);  
    sectionIntegr->getWeightsDeriv(numFibers, areaDeriv);
  }
  else {
    for (int i = 0; i < numFibers; i++) {
      locsDeriv[i] = 0.0;
      areaDeriv[i] = 0.0;
    }
  }
  
  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  double drootAlphadh = 0.0;
  if (parameterID == 1)
    drootAlphadh = 0.5/rootAlpha;

  for (int i = 0; i < numFibers; i++) {
    y = fiberLocs[i] - yBar;
    A = fiberArea[i];
    
    dsigdh = theMaterials[i]->getStressSensitivity(gradIndex,conditional);
    ds(0) += dsigdh(0)*A;
    ds(1) += -y*dsigdh(0)*A;
    ds(2) += rootAlpha*dsigdh(1)*A;

    if (areaDeriv[i] != 0.0 || locsDeriv[i] != 0.0 || parameterID == 1)
      stress = theMaterials[i]->getStress();

    if (locsDeriv[i] != 0.0 || parameterID == 1)
      tangent = theMaterials[i]->getTangent();

    if (areaDeriv[i] != 0.0) {
      sig_dAdh(0) = stress(0)*areaDeriv[i];
      sig_dAdh(1) = stress(1)*areaDeriv[i];
      
      ds(0) += sig_dAdh(0);
      ds(1) += -y*sig_dAdh(0);
      ds(2) += rootAlpha*sig_dAdh(1);
    }

    if (locsDeriv[i] != 0.0) {
      //ds(0) += 0.0;
      ds(1) += -locsDeriv[i] * (stress(0)*A);
      //ds(2) += 0.0;

      ds(0) += (-locsDeriv[i]*tangent(0,0)*e(1))*A;
      ds(1) += -y*(-locsDeriv[i]*tangent(0,0)*e(1))*A;
      ds(2) += rootAlpha*(-locsDeriv[i]*tangent(1,0)*e(1))*A;
    }

    if (parameterID == 1) {
      //ds(0) += 0.0;
      //ds(1) += 0.0;
      ds(2) += drootAlphadh * (stress(1)*A);

      ds(0) += (drootAlphadh*tangent(0,1)*e(2))*A;
      ds(1) += -y*(drootAlphadh*tangent(0,1)*e(2))*A;
      //ds(2) += rootAlpha*(drootAlphadh*tangent(1,1)*e(2))*A;
      ds(2) += 0.5*tangent(1,1)*e(2)*A;
    }
  }

  return ds;
}

const Matrix &
NDFiberSection2d::getInitialTangentSensitivity(int gradIndex)
{
  static Matrix dksdh(3,3);
  
  dksdh.Zero();
  /*
  double y, A, dydh, dAdh, tangent, dtangentdh;

  static double fiberLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      fiberLocs[i] = matData[2*i];
      fiberArea[i] = matData[2*i+1];
    }
  }

  static double locsDeriv[10000];
  static double areaDeriv[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getLocationsDeriv(numFibers, locsDeriv);  
    sectionIntegr->getWeightsDeriv(numFibers, areaDeriv);
  }
  else {
    for (int i = 0; i < numFibers; i++) {
      locsDeriv[i] = 0.0;
      areaDeriv[i] = 0.0;
    }
  }
  
  for (int i = 0; i < numFibers; i++) {
    y = fiberLocs[i] - yBar;
    A = fiberArea[i];
    dydh = locsDeriv[i];
    dAdh = areaDeriv[i];
    
    tangent = theMaterials[i]->getInitialTangent();
    dtangentdh = theMaterials[i]->getInitialTangentSensitivity(gradIndex);

    dksdh(0,0) += dtangentdh*A + tangent*dAdh;

    dksdh(0,1) += -y*(dtangentdh*A+tangent*dAdh) - dydh*(tangent*A);

    dksdh(1,1) += 2*(y*dydh*tangent*A) + y*y*(dtangentdh*A+tangent*dAdh);
  }

  dksdh(1,0) = dksdh(0,1);
  */
  return dksdh;
}

int
NDFiberSection2d::commitSensitivity(const Vector& defSens,
				    int gradIndex, int numGrads)
{
  double d0 = defSens(0);
  double d1 = defSens(1);
  double d2 = defSens(2);

  dedh = defSens;

  static double fiberLocs[10000];

  if (sectionIntegr != 0)
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
  else {
    for (int i = 0; i < numFibers; i++)
      fiberLocs[i] = matData[2*i];
  }

  static double locsDeriv[10000];
  static double areaDeriv[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getLocationsDeriv(numFibers, locsDeriv);  
    sectionIntegr->getWeightsDeriv(numFibers, areaDeriv);
  }
  else {
    for (int i = 0; i < numFibers; i++) {
      locsDeriv[i] = 0.0;
      areaDeriv[i] = 0.0;
    }
  }

  double y;
  double kappa = e(1);
  double gamma = e(2);

  static Vector depsdh(2);

  double rootAlpha = 1.0;
  if (alpha != 1.0)
    rootAlpha = sqrt(alpha);

  double drootAlphadh = 0.0;
  if (parameterID == 1)
    drootAlphadh = 0.5/rootAlpha;

  depsdh(1) = rootAlpha*d2 + drootAlphadh*gamma;

  for (int i = 0; i < numFibers; i++) {
    NDMaterial *theMat = theMaterials[i];
    y = fiberLocs[i] - yBar;

    // determine material strain and set it
    depsdh(0) = d0 - y*d1 - locsDeriv[i]*kappa;
    theMat->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////
