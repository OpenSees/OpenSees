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
                                                                        
// $Revision: 1.35 $
// $Date: 2008-11-04 21:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSection2d.cpp,v $
                                                                        
// Written: fmk
// Created: 04/04
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
#include <FiberSection2d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <SectionIntegration.h>
#include <elementAPI.h>

ID FiberSection2d::code(2);

void* OPS_FiberSection2d()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	opserr<<"insufficient arguments for FiberSection2d\n";
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
    return new FiberSection2d(tag, num, computeCentroid);
}

// constructors:
FiberSection2d::FiberSection2d(int tag, int num, Fiber **fibers, bool compCentroid): 
  SectionForceDeformation(tag, SEC_TAG_FiberSection2d),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(compCentroid),
  sectionIntegr(0), e(2), s(0), ks(0), dedh(2)

{
  if (numFibers > 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSection2d::FiberSection2d -- failed to allocate Material pointers";
      exit(-1);
    }

    matData = new double [numFibers*2];

    if (matData == 0) {
      opserr << "FiberSection2d::FiberSection2d -- failed to allocate double array for material data\n";
      exit(-1);
    }

    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();
      ABar  += Area;
      QzBar += yLoc*Area;
      matData[i*2] = yLoc;
      matData[i*2+1] = Area;
      UniaxialMaterial *theMat = theFiber->getMaterial();
      theMaterials[i] = theMat->getCopy();

      if (theMaterials[i] == 0) {
	opserr << "FiberSection2d::FiberSection2d -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    

    if (computeCentroid)
      yBar = QzBar/ABar;
  }

  s = new Vector(sData, 2);
  ks = new Matrix(kData, 2, 2);

  sData[0] = 0.0;
  sData[1] = 0.0;

  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
}

// allocate memory for fibers
FiberSection2d::FiberSection2d(int tag, int num, bool compCentroid): 
  SectionForceDeformation(tag, SEC_TAG_FiberSection2d),
  numFibers(0), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(compCentroid),
  sectionIntegr(0), e(2), s(0), ks(0), dedh(2)
{
    if(sizeFibers > 0) {
	theMaterials = new UniaxialMaterial *[sizeFibers];

	if(theMaterials == 0) {
	    opserr << "FiberSection2d::FiberSection2d -- failed to allocate Material pointers";
	    exit(-1);
	}

	matData = new double [sizeFibers*2];

	if(matData == 0) {
	    opserr << "FiberSection2d::FiberSection2d -- failed to allocate double array for material data\n";
	    exit(-1);
	}

	for(int i = 0; i < sizeFibers; i++) {
	    matData[i*2] = 0.0;
	    matData[i*2+1] = 0.0;
	    theMaterials[i] = 0;
	}    
    }

    s = new Vector(sData, 2);
    ks = new Matrix(kData, 2, 2);

    sData[0] = 0.0;
    sData[1] = 0.0;

    kData[0] = 0.0;
    kData[1] = 0.0;
    kData[2] = 0.0;
    kData[3] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
}

FiberSection2d::FiberSection2d(int tag, int num, UniaxialMaterial **mats,
			       SectionIntegration &si, bool compCentroid):
  SectionForceDeformation(tag, SEC_TAG_FiberSection2d),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(compCentroid),
  sectionIntegr(0), e(2), s(0), ks(0), dedh(2)
{
  if (numFibers != 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSection2d::FiberSection2d -- failed to allocate Material pointers";
      exit(-1);
    }
    matData = new double [numFibers*2];

    if (matData == 0) {
      opserr << "FiberSection2d::FiberSection2d -- failed to allocate double array for material data\n";
      exit(-1);
    }
  }

  sectionIntegr = si.getCopy();
  if (sectionIntegr == 0) {
    opserr << "Error: FiberSection2d::FiberSection2d: could not create copy of section integration object" << endln;
    exit(-1);
  }

  static double fiberLocs[10000];
  sectionIntegr->getFiberLocations(numFibers, fiberLocs);
  
  static double fiberArea[10000];
  sectionIntegr->getFiberWeights(numFibers, fiberArea);

  for (int i = 0; i < numFibers; i++) {

    ABar  += fiberArea[i];
    QzBar += fiberLocs[i]*fiberArea[i];

    theMaterials[i] = mats[i]->getCopy();
    
    if (theMaterials[i] == 0) {
      opserr << "FiberSection2d::FiberSection2d -- failed to get copy of a Material\n";
      exit(-1);
    }
  }    

  if (computeCentroid)
    yBar = QzBar/ABar;
  
  s = new Vector(sData, 2);
  ks = new Matrix(kData, 2, 2);
  
  sData[0] = 0.0;
  sData[1] = 0.0;
  
  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;
  
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSection2d::FiberSection2d():
  SectionForceDeformation(0, SEC_TAG_FiberSection2d),
  numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(true),
  sectionIntegr(0), e(2), s(0), ks(0), dedh(2)
{
  s = new Vector(sData, 2);
  ks = new Matrix(kData, 2, 2);

  sData[0] = 0.0;
  sData[1] = 0.0;

  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
}

int
FiberSection2d::addFiber(Fiber &newFiber)
{
  // need to create larger arrays
  if(numFibers == sizeFibers) {
      int newsize = 2*sizeFibers;
      if(newsize == 0) newsize = 30;
      UniaxialMaterial **newArray = new UniaxialMaterial *[newsize]; 
      double *newMatData = new double [2 * newsize];
      if (newArray == 0 || newMatData == 0) {
	  opserr <<"FiberSection2d::addFiber -- failed to allocate Fiber pointers\n";
	  return -1;
      }

      // copy the old pointers and data
      int i;
      for (i = 0; i < sizeFibers; i++) {
	  newArray[i] = theMaterials[i];
	  newMatData[2*i] = matData[2*i];
	  newMatData[2*i+1] = matData[2*i+1];
      }

      // initialize new memory
      for(i = sizeFibers; i<newsize; i++) {
	  newArray[i] = 0;
	  newMatData[2*i] = 0.0;
	  newMatData[2*i+1] = 0.0;
      }

      sizeFibers = newsize;

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
  UniaxialMaterial *theMat = newFiber.getMaterial();
  theMaterials[numFibers] = theMat->getCopy();

  if(theMaterials[numFibers] == 0) {
    opserr <<"FiberSection2d::addFiber -- failed to get copy of a Material\n";
    return -1;
  }

  numFibers++;

  // Recompute centroid
  if (computeCentroid) {
    ABar += Area;
    QzBar += yLoc*Area;
    yBar = QzBar/ABar;
  }
  
  return 0;
}


// destructor:
FiberSection2d::~FiberSection2d()
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
FiberSection2d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;

  e = deforms;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;

  double d0 = deforms(0);
  double d1 = deforms(1);

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
    UniaxialMaterial *theMat = theMaterials[i];
    double y = fiberLocs[i] - yBar;
    double A = fiberArea[i];

    // determine material strain and set it
    double strain = d0 - y*d1;
    double tangent, stress;
    res += theMat->setTrial(strain, stress, tangent);

    double ks0 = tangent * A;
    double ks1 = ks0 * -y;
    kData[0] += ks0;
    kData[1] += ks1;
    kData[3] += ks1 * -y;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * -y;
  }

  kData[2] = kData[1];

  return res;
}

const Vector&
FiberSection2d::getSectionDeformation(void)
{
  return e;
}

const Matrix&
FiberSection2d::getInitialTangent(void)
{
  static double kInitial[4];
  static Matrix kInitialMatrix(kInitial, 2, 2);
  kInitial[0] = 0.0; kInitial[1] = 0.0; kInitial[2] = 0.0; kInitial[3] = 0.0;

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
    UniaxialMaterial *theMat = theMaterials[i];
    double y = fiberLocs[i] - yBar;
    double A = fiberArea[i];

    double tangent = theMat->getInitialTangent();

    double ks0 = tangent * A;
    double ks1 = ks0 * -y;
    kInitial[0] += ks0;
    kInitial[1] += ks1;
    kInitial[3] += ks1 * -y;
  }

  kInitial[2] = kInitial[1];

  return kInitialMatrix;
}

const Matrix&
FiberSection2d::getSectionTangent(void)
{
  return *ks;
}

const Vector&
FiberSection2d::getStressResultant(void)
{
  return *s;
}

SectionForceDeformation*
FiberSection2d::getCopy(void)
{
  FiberSection2d *theCopy = new FiberSection2d ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;
  theCopy->sizeFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr <<"FiberSection2d::getCopy -- failed to allocate Material pointers\n";
      exit(-1);
    }
  
    theCopy->matData = new double [numFibers*2];

    if (theCopy->matData == 0) {
      opserr << "FiberSection2d::getCopy -- failed to allocate double array for material data\n";
      exit(-1);
    }
			    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*2] = matData[i*2];
      theCopy->matData[i*2+1] = matData[i*2+1];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr <<"FiberSection2d::getCopy -- failed to get copy of a Material";
	exit(-1);
      }
    }  
  }

  theCopy->e = e;
  theCopy->QzBar = QzBar;
  theCopy->ABar = ABar;
  theCopy->yBar = yBar;

  theCopy->kData[0] = kData[0];
  theCopy->kData[1] = kData[1];
  theCopy->kData[2] = kData[2];
  theCopy->kData[3] = kData[3];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];

  theCopy->computeCentroid = computeCentroid;
  
  if (sectionIntegr != 0)
    theCopy->sectionIntegr = sectionIntegr->getCopy();
  else
    theCopy->sectionIntegr = 0;

  return theCopy;
}

const ID&
FiberSection2d::getType ()
{
  return code;
}

int
FiberSection2d::getOrder () const
{
  return 2;
}

int
FiberSection2d::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  return err;
}

int
FiberSection2d::revertToLastCommit(void)
{
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;
  
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
    UniaxialMaterial *theMat = theMaterials[i];
    double y = fiberLocs[i] - yBar;
    double A = fiberArea[i];

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    // get material stress & tangent for this strain and determine ks and fs
    double tangent = theMat->getTangent();
    double stress = theMat->getStress();
    double ks0 = tangent * A;
    double ks1 = ks0 * -y;
    kData[0] += ks0;
    kData[1] += ks1;
    kData[3] += ks1 * -y;

    double fs0 = stress * A;
    sData[0] = fs0;
    sData[1] = fs0 * -y;
  }

  kData[2] = kData[1];

  return err;
}

int
FiberSection2d::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;
  
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
    UniaxialMaterial *theMat = theMaterials[i];
    double y = fiberLocs[i] - yBar;
    double A = fiberArea[i];

    // invoke revertToLast on the material
    err += theMat->revertToStart();

    // get material stress & tangent for this strain and determine ks and fs
    double tangent = theMat->getTangent();
    double stress = theMat->getStress();
    double ks0 = tangent * A;
    double ks1 = ks0 * -y;
    kData[0] += ks0;
    kData[1] += ks1;
    kData[3] += ks1 * -y;

    double fs0 = stress * A;
    sData[0] = fs0;
    sData[1] = fs0 * -y;
  }

  kData[2] = kData[1];

  return err;
}

int
FiberSection2d::sendSelf(int commitTag, Channel &theChannel)
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
    opserr <<  "FiberSection2d::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (numFibers != 0) {
    
    // create an id containingg classTag and dbTag for each material & send it
    ID materialData(2*numFibers);
    for (int i=0; i<numFibers; i++) {
      UniaxialMaterial *theMat = theMaterials[i];
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
      opserr <<  "FiberSection2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 2*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "FiberSection2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);

  }

  return res;
}

int
FiberSection2d::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "FiberSection2d::recvSelf - failed to recv ID data\n";
    return res;
  }    
  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "FiberSection2d::recvSelf - failed to recv material data\n";
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
	theMaterials = new UniaxialMaterial *[numFibers];
	
	if (theMaterials == 0) {
	  opserr <<"FiberSection2d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}
	
	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;

	matData = new double [numFibers*2];

	if (matData == 0) {
	  opserr <<"FiberSection2d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, 2*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "FiberSection2d::recvSelf - failed to recv material data\n";
      return res;
    }    

    int i;
    for (i=0; i<numFibers; i++) {
      int classTag = materialData(2*i);
      int dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of corrcet type, 
      // release old and create a new one
      if (theMaterials[i] == 0)
	theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);
      else if (theMaterials[i]->getClassTag() != classTag) {
	delete theMaterials[i];
	theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);      
      }

      if (theMaterials[i] == 0) {
	opserr <<"FiberSection2d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    QzBar = 0.0;
    ABar  = 0.0;
    double yLoc, Area;

    computeCentroid = data(2) ? true : false;
    
    // Recompute centroid
    for (i = 0; computeCentroid && i < numFibers; i++) {
      yLoc = matData[2*i];
      Area = matData[2*i+1];
      ABar  += Area;
      QzBar += yLoc*Area;
    }

    if (computeCentroid)
      yBar = QzBar/ABar;
    else
      yBar = 0.0;
  }    

  return res;
}

void
FiberSection2d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "\nFiberSection2d, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: " << yBar << endln;
    
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
      for (int i = 0; i < numFibers; i++) {
	s << "\nLocation (y) = (" << matData[2*i] << ")";
	s << "\nArea = " << matData[2*i+1] << endln;
	theMaterials[i]->Print(s, flag);
      }
    }
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
	s << "\"name\": \"" << this->getTag() << "\", ";
	s << "\"type\": \"FiberSection2d\", ";
    s << "\"fibers\": [\n";
    for (int i = 0; i < numFibers; i++) {
      s << "\t\t\t\t{\"coord\": [" << matData[2*i] << ", 0.0], ";
      s << "\"area\": " << matData[2*i+1] << ", ";
      s << "\"material\": \"" << theMaterials[i]->getTag() << "\"";
      if (i < numFibers-1)
	s << "},\n";
      else
	s << "}\n";	
    }
    s << "\t\t\t]}";
  }
}

Response*
FiberSection2d::setResponse(const char **argv, int argc,
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

  } else if (strcmp(argv[0],"fiberData") == 0) {
    int numData = numFibers*5;
    for (int j = 0; j < numFibers; j++) {
      output.tag("FiberOutput");
      output.attr("yLoc", matData[2*j]);
      output.attr("zLoc", 0.0);
      output.attr("area", matData[2*j+1]);    
      output.tag("ResponseType","yCoord");
      output.tag("ResponseType","zCoord");
      output.tag("ResponseType","area");
      output.tag("ResponseType","stress");
      output.tag("ResponseType","strain");
      output.endTag();
    }
    Vector theResponseData(numData);
    return theResponse = new MaterialResponse(this, 5, theResponseData);

  
  } else if ((strcmp(argv[0],"numFailedFiber") == 0) || (strcmp(argv[0],"numFiberFailed") == 0)) {
    int count = 0;
    return theResponse = new MaterialResponse(this, 6, count);

  } else if ((strcmp(argv[0],"sectionFailed") == 0) || 
	     (strcmp(argv[0],"hasSectionFailed") == 0) ||
	     (strcmp(argv[0],"hasFailed") == 0)) {
    int count = 0;
    return theResponse = new MaterialResponse(this, 7, count);
  }
  //by SAJalali
  else if ((strcmp(argv[0], "energy") == 0) || (strcmp(argv[0], "Energy") == 0)) {
	  return theResponse = new MaterialResponse(this, 8, getEnergy());
  }

// If not a fiber response, call the base class method
return SectionForceDeformation::setResponse(argv, argc, output);
}


int 
FiberSection2d::getResponse(int responseID, Information &sectInfo)
{
  if (responseID == 5) {
    int numData = 5*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      double yLoc, zLoc, A, stress, strain;
      yLoc = matData[2*j];
      zLoc = 0.0;
      A = matData[2*j+1];
      stress = theMaterials[j]->getStress();
      strain = theMaterials[j]->getStrain();
      data(count) = yLoc; data(count+1) = zLoc; data(count+2) = A;
      data(count+3) = stress; data(count+4) = strain;
      count += 5;
    }
    return sectInfo.setVector(data);	
  } else  if (responseID == 6) {
    int count = 0;
    for (int j = 0; j < numFibers; j++) {    
      if (theMaterials[j]->hasFailed() == true)
	count++;
    }
    return sectInfo.setInt(count);

  } else  if (responseID == 7) {
    int count = 0;
    for (int j = 0; j < numFibers; j++) {    
      if (theMaterials[j]->hasFailed() == true) {
	count+=1;
      }
    }
    if (count == numFibers)
      count = 1;
    else
      count = 0;

    return sectInfo.setInt(count);
  } 
  //by SAJalali
  else if (responseID == 8) {
	  return sectInfo.setDouble(getEnergy());
  }

  return SectionForceDeformation::getResponse(responseID, sectInfo);
}



// AddingSensitivity:BEGIN ////////////////////////////////////
int
FiberSection2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  // Check if the parameter belongs to the material
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

  // Check if the parameter belongs to a fiber
  // unlike setResponse, only allowing 'fiber y z matTag ...' because
  // the setResponse logic breaks down with the trailing arguments
  if (strstr(argv[0],"fiber") != 0) {
    
    int key = numFibers;
    int passarg = 2;
    
    if (argc < 5)
      return 0;

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
      passarg = 4;
    }
    
    // Finally, call setParameter
    if (key >= 0 && key < numFibers)
      return theMaterials[key]->setParameter(&argv[passarg], argc-passarg, param);
  }

  // Check if it belongs to the section integration
  if (strstr(argv[0],"integration") != 0) {
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

const Vector &
FiberSection2d::getSectionDeformationSensitivity(int gradIndex)
{
  static Vector dummy(2);

  return dummy;
}

const Vector &
FiberSection2d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(2);
  
  ds.Zero();
  
  double y, A;
  double stressGradient = 0.0;
  double stress = 0.0;
  double tangent = 0.0;
  double sig_dAdh = 0.0;

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
    
    stressGradient = theMaterials[i]->getStressSensitivity(gradIndex,true);
    stressGradient = stressGradient * A;

    ds(0) += stressGradient;
    ds(1) += stressGradient * -y;

    if (areaDeriv[i] != 0.0 || locsDeriv[i] != 0.0)
      stress = theMaterials[i]->getStress();

    if (areaDeriv[i] != 0.0) {
      sig_dAdh = stress*areaDeriv[i];
      
      ds(0) += sig_dAdh;
      ds(1) += sig_dAdh * -y;
    }

    if (locsDeriv[i] != 0.0) {
      //ds(0) += 0.0;
      ds(1) += (stress*A) * -locsDeriv[i];
      
      tangent = theMaterials[i]->getTangent();
      tangent = tangent * A * e(1);
      
      ds(0) += -locsDeriv[i]*tangent;
      ds(1) += fiberLocs[i]*locsDeriv[i]*tangent;
    }

    //opserr << locsDeriv[i] << ' ' << areaDeriv[i] << endln;
  }
  
  return ds;
}

const Matrix &
FiberSection2d::getInitialTangentSensitivity(int gradIndex)
{
  static Matrix dksdh(2,2);
  
  dksdh.Zero();

  double y, A, dydh, dAdh;
  double tangent = 0.0;
  double dtangentdh = 0.0;

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

  return dksdh;
}

int
FiberSection2d::commitSensitivity(const Vector& defSens,
				  int gradIndex, int numGrads)
{
  double d0 = defSens(0);
  double d1 = defSens(1);

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

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    y = fiberLocs[i] - yBar;

    // determine material strain and set it
    double strainSens = d0 - y*d1 - locsDeriv[i]*kappa;
    theMat->commitSensitivity(strainSens,gradIndex,numGrads);
  }

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////

//by SAJalali
double FiberSection2d::getEnergy() const
{
	static double fiberArea[10000];

	if (sectionIntegr != 0) {
		sectionIntegr->getFiberWeights(numFibers, fiberArea);
	}
	else {
		for (int i = 0; i < numFibers; i++) {
			fiberArea[i] = matData[2 * i + 1];
		}
	}
	double energy = 0;
	for (int i = 0; i < numFibers; i++)
	{
		double A = fiberArea[i];
		energy += A * theMaterials[i]->getEnergy();
	}
	return energy;
}
