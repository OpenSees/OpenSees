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
                                                                        
// $Revision: 1.32 $
// $Date: 2010-08-16 05:05:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSection3d.cpp,v $
                                                                        
// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSection3d.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSection3d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include <SectionIntegration.h>
#include <elementAPI.h>
#include <string.h>

ID FiberSection3d::code(4);

void* OPS_FiberSection3d()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	    opserr<<"insufficient arguments for FiberSection3d\n";
	    return 0;
    }
    
    numData = 1;
    int tag;
    if (OPS_GetIntInput(&numData, &tag) < 0) return 0;

    if (OPS_GetNumRemainingInputArgs() < 2) {
      opserr << "WARNING torsion not specified for FiberSection\n";
      opserr << "Use either -GJ $GJ or -torsion $matTag\n";
      opserr << "\nFiberSection3d section: " << tag << endln;
      return 0;
    }
    
    UniaxialMaterial *torsion = 0;
    bool deleteTorsion = false;
    bool computeCentroid = true;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      const char* opt = OPS_GetString();
      if (strcmp(opt,"-noCentroid") == 0) {
	computeCentroid = false;
      }
      if (strcmp(opt, "-GJ") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
	numData = 1;
	double GJ;
	if (OPS_GetDoubleInput(&numData, &GJ) < 0) {
	  opserr << "WARNING: failed to read GJ\n";
	  return 0;
	}
	torsion = new ElasticMaterial(0,GJ);
	deleteTorsion = true;
      }
      if (strcmp(opt, "-torsion") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
	numData = 1;
	int torsionTag;
	if (OPS_GetIntInput(&numData, &torsionTag) < 0) {
	  opserr << "WARNING: failed to read torsion\n";
	  return 0;
	}
	torsion = OPS_getUniaxialMaterial(torsionTag);
      }
    }

    if (torsion == 0) {
      opserr << "WARNING torsion not specified for FiberSection\n";
      opserr << "\nFiberSection3d section: " << tag << endln;
      return 0;
    }
    
    int num = 30;
    SectionForceDeformation *section = new FiberSection3d(tag, num, *torsion, computeCentroid);
    if (deleteTorsion)
      delete torsion;
    return section;
}

// constructors:
FiberSection3d::FiberSection3d(int tag, int num, Fiber **fibers,
			       UniaxialMaterial &torsion, bool compCentroid): 
  SectionForceDeformation(tag, SEC_TAG_FiberSection3d),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
  sectionIntegr(0), e(4), s(0), ks(0), theTorsion(0)
{
  if (numFibers != 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate Material pointers\n";
      exit(-1);
    }

    matData = new double [numFibers*3];

    if (matData == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate double array for material data\n";
      exit(-1);
    }

    double yLoc, zLoc, Area;
    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();

      QzBar += yLoc*Area;
      QyBar += zLoc*Area;
      Abar  += Area;

      matData[i*3] = yLoc;
      matData[i*3+1] = zLoc;
      matData[i*3+2] = Area;
      UniaxialMaterial *theMat = theFiber->getMaterial();
      theMaterials[i] = theMat->getCopy();

      if (theMaterials[i] == 0) {
	opserr << "FiberSection3d::FiberSection3d -- failed to get copy of a Material\n";
	exit(-1);
      }
    }

    if (computeCentroid) {
      yBar = QzBar/Abar;
      zBar = QyBar/Abar;
    }
  }

  theTorsion = torsion.getCopy();
  if (theTorsion == 0)
    opserr << "FiberSection3d::FiberSection3d -- failed to get copy of torsion material\n";

  s = new Vector(sData, 4);
  ks = new Matrix(kData, 4, 4);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;

  for (int i=0; i<16; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
}

FiberSection3d::FiberSection3d(int tag, int num, UniaxialMaterial &torsion, bool compCentroid): 
    SectionForceDeformation(tag, SEC_TAG_FiberSection3d),
    numFibers(0), sizeFibers(num), theMaterials(0), matData(0),
    QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
    sectionIntegr(0), e(4), s(0), ks(0), theTorsion(0)
{
    if(sizeFibers != 0) {
	theMaterials = new UniaxialMaterial *[sizeFibers];

	if (theMaterials == 0) {
	    opserr << "FiberSection3d::FiberSection3d -- failed to allocate Material pointers\n";
	    exit(-1);
	}

	matData = new double [sizeFibers*3];

	if (matData == 0) {
	    opserr << "FiberSection3d::FiberSection3d -- failed to allocate double array for material data\n";
	    exit(-1);
	}

	for (int i = 0; i < sizeFibers; i++) {
	    matData[i*3] = 0.0;
	    matData[i*3+1] = 0.0;
	    matData[i*3+2] = .0;
	    theMaterials[i] = 0;
	}
    }

    theTorsion = torsion.getCopy();
    if (theTorsion == 0) 
      opserr << "FiberSection3d::FiberSection3d -- failed to get copy of torsion material\n";

    s = new Vector(sData, 4);
    ks = new Matrix(kData, 4, 4);

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;
    sData[3] = 0.0;

    for (int i=0; i<16; i++)
	kData[i] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_T;
}

FiberSection3d::FiberSection3d(int tag, int num, UniaxialMaterial **mats,
			       SectionIntegration &si, UniaxialMaterial &torsion,
			       bool compCentroid):
  SectionForceDeformation(tag, SEC_TAG_FiberSection3d),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
  sectionIntegr(0), e(4), s(0), ks(0), theTorsion(0)
{
  if (numFibers != 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate Material pointers";
      exit(-1);
    }
    matData = new double [numFibers*3];

    if (matData == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate double array for material data\n";
      exit(-1);
    }
  }

  sectionIntegr = si.getCopy();
  if (sectionIntegr == 0) {
    opserr << "Error: FiberSection3d::FiberSection3d: could not create copy of section integration object" << endln;
    exit(-1);
  }

  static double yLocs[10000];
  static double zLocs[10000];
  sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
  
  static double fiberArea[10000];
  sectionIntegr->getFiberWeights(numFibers, fiberArea);
  
  for (int i = 0; i < numFibers; i++) {

    Abar  += fiberArea[i];
    QzBar += yLocs[i]*fiberArea[i];
    QyBar += zLocs[i]*fiberArea[i];

    theMaterials[i] = mats[i]->getCopy();
    
    if (theMaterials[i] == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to get copy of a Material\n";
      exit(-1);
    }
  }    

  if (computeCentroid) {
    yBar = QzBar/Abar;  
    zBar = QyBar/Abar;  
  }
  
  theTorsion = torsion.getCopy();
  if (theTorsion == 0)
    opserr << "FiberSection3d::FiberSection3d -- failed to get copy of torsion material\n";

  s = new Vector(sData, 4);
  ks = new Matrix(kData, 4, 4);
  
  for (int i = 0; i < 4; i++)
    sData[i] = 0.0;

  for (int i = 0; i < 16; i++)
    kData[i] = 0.0;
  
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSection3d::FiberSection3d():
  SectionForceDeformation(0, SEC_TAG_FiberSection3d),
  numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), computeCentroid(true),
  sectionIntegr(0), e(4), s(0), ks(0), theTorsion(0)
{
  s = new Vector(sData, 4);
  ks = new Matrix(kData, 4, 4);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;

  for (int i=0; i<16; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
}

int
FiberSection3d::addFiber(Fiber &newFiber)
{
  // need to create a larger array
  if(numFibers == sizeFibers) {
      int newSize = 2*sizeFibers;
      UniaxialMaterial **newArray = new UniaxialMaterial *[newSize]; 
      double *newMatData = new double [3 * newSize];
      
      if (newArray == 0 || newMatData == 0) {
	  opserr << "FiberSection3d::addFiber -- failed to allocate Fiber pointers\n";
	  exit(-1);
      }

      // copy the old pointers
      for (int i = 0; i < numFibers; i++) {
	  newArray[i] = theMaterials[i];
	  newMatData[3*i] = matData[3*i];
	  newMatData[3*i+1] = matData[3*i+1];
	  newMatData[3*i+2] = matData[3*i+2];
      }

      // initialize new memomry
      for (int i = numFibers; i < newSize; i++) {
	  newArray[i] = 0;
	  newMatData[3*i] = 0.0;
	  newMatData[3*i+1] = 0.0;
	  newMatData[3*i+2] = 0.0;
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
	    
  // set the new pointers
  double yLoc, zLoc, Area;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  matData[numFibers*3] = yLoc;
  matData[numFibers*3+1] = zLoc;
  matData[numFibers*3+2] = Area;
  UniaxialMaterial *theMat = newFiber.getMaterial();
  theMaterials[numFibers] = theMat->getCopy();

  if (theMaterials[numFibers] == 0) {
    opserr << "FiberSection3d::addFiber -- failed to get copy of a Material\n";
    return -1;
  }

  numFibers++;

  // Recompute centroid
  if (computeCentroid) {
    Abar  += Area;
    QzBar += yLoc*Area;
    QyBar += zLoc*Area;
    
    yBar = QzBar/Abar;
    zBar = QyBar/Abar;
  }
  
  return 0;
}



// destructor:
FiberSection3d::~FiberSection3d()
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

  if (theTorsion != 0)
    delete theTorsion;
}

int
FiberSection3d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;
 
  for (int i = 0; i < 4; i++)
    sData[i] = 0.0;
  for (int i = 0; i < 16; i++)
    kData[i] = 0.0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);
  double d3 = deforms(3);

  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];
 
  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
	
    for (int i = 0; i < numFibers; i++) {
		
      yLocs[i] = matData[3*i];
      zLocs[i] = matData[3*i+1];
      fiberArea[i] = matData[3*i+2];
    }
  }
 
  double tangent, stress;
  for (int i = 0; i < numFibers; i++) {
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];

    // determine material strain and set it
    double strain = d0 - y*d1 + z*d2;
    res += theMaterials[i]->setTrial(strain, stress, tangent);

    double value = tangent * A;
    double vas1 = -y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;
    
    kData[5] += vas1 * -y;
    kData[6] += vas1as2;
    
    kData[10] += vas2 * z; 

    double fs0 = stress * A;

    sData[0] += fs0;
    sData[1] += fs0 * -y;
    sData[2] += fs0 * z;
  }

  kData[4] = kData[1];
  kData[8] = kData[2];
  kData[9] = kData[6];
 
  if (theTorsion != 0) {
    res += theTorsion->setTrial(d3, stress, tangent);
    sData[3] = stress;
    kData[15] = tangent;
  }

  return res;
}

const Matrix&
FiberSection3d::getInitialTangent(void)
{
  static double kInitialData[16];
  static Matrix kInitial(kInitialData, 4, 4);
  
  kInitial.Zero();

  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[3*i];
      zLocs[i] = matData[3*i+1];
      fiberArea[i] = matData[3*i+2];
    }
  }

  for (int i = 0; i < numFibers; i++) {
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];

    double tangent = theMaterials[i]->getInitialTangent();

    double value = tangent * A;
    double vas1 = -y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kInitialData[0] += value;
    kInitialData[1] += vas1;
    kInitialData[2] += vas2;
    
    kInitialData[5] += vas1 * -y;
    kInitialData[6] += vas1as2;
    
    kInitialData[10] += vas2 * z; 
  }

  kInitialData[4] = kInitialData[1];
  kInitialData[8] = kInitialData[2];
  kInitialData[9] = kInitialData[6];

  if (theTorsion != 0)
    kInitialData[15] = theTorsion->getInitialTangent();

  return kInitial;
}

const Vector&
FiberSection3d::getSectionDeformation(void)
{
  return e;
}

const Matrix&
FiberSection3d::getSectionTangent(void)
{
  return *ks;
}

const Vector&
FiberSection3d::getStressResultant(void)
{
  return *s;
}

SectionForceDeformation*
FiberSection3d::getCopy(void)
{
  FiberSection3d *theCopy = new FiberSection3d ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;
  theCopy->sizeFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate Material pointers\n";
      exit(-1);			    
    }

    theCopy->matData = new double [numFibers*3];

    if (theCopy->matData == 0) {
      opserr << "FiberSection3d::FiberSection3d -- failed to allocate double array for material data\n";
      exit(-1);
    }
			    
    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*3] = matData[i*3];
      theCopy->matData[i*3+1] = matData[i*3+1];
      theCopy->matData[i*3+2] = matData[i*3+2];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr << "FiberSection3d::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    
  }

  theCopy->e = e;
  theCopy->QzBar = QzBar;
  theCopy->QyBar = QyBar;
  theCopy->Abar = Abar;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;
  theCopy->computeCentroid = computeCentroid;
  
  for (int i=0; i<16; i++)
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];
  theCopy->sData[3] = sData[3];

  if (theTorsion != 0)
    theCopy->theTorsion = theTorsion->getCopy();
  else
    theCopy->theTorsion = 0;

  if (sectionIntegr != 0)
    theCopy->sectionIntegr = sectionIntegr->getCopy();
  else
    theCopy->sectionIntegr = 0;

  return theCopy;
}

const ID&
FiberSection3d::getType ()
{
  return code;
}

int
FiberSection3d::getOrder () const
{
  return 4;
}

int
FiberSection3d::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  if (theTorsion != 0)
    err += theTorsion->commitState();

  return err;
}

int
FiberSection3d::revertToLastCommit(void)
{
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0; 
  kData[15] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; sData[3] = 0.0;

  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[3*i];
      zLocs[i] = matData[3*i+1];
      fiberArea[i] = matData[3*i+2];
    }
  }

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = -y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;
    
    kData[5] += vas1 * -y;
    kData[6] += vas1as2;
    
    kData[10] += vas2 * z; 

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * -y;
    sData[2] += fs0 * z;
  }

  kData[4] = kData[1];
  kData[8] = kData[2];
  kData[9] = kData[6];

  if (theTorsion != 0) {
    err += theTorsion->revertToLastCommit();
    kData[15] = theTorsion->getTangent();
  } else
    kData[15] = 0.0;

  return err;
}

int
FiberSection3d::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0;
  kData[15] = 0.0; 
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; sData[3] = 0.0;

  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[3*i];
      zLocs[i] = matData[3*i+1];
      fiberArea[i] = matData[3*i+2];
    }
  }

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];

    // invoke revertToStart on the material
    err += theMat->revertToStart();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = -y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;
    
    kData[5] += vas1 * -y;
    kData[6] += vas1as2;
    
    kData[10] += vas2 * z; 

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * -y;
    sData[2] += fs0 * z;
  }

  kData[4] = kData[1];
  kData[8] = kData[2];
  kData[9] = kData[6];

  if (theTorsion != 0) {
    err += theTorsion->revertToStart();
    kData[15] = theTorsion->getTangent();
    sData[3] = theTorsion->getStress();
  } else {
    kData[15] = 0.0;
    sData[3] = 0.0;
  }

  return err;
}

int
FiberSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 5 so no conflict with matData below if just 2 fibers
  static ID data(5);
  data(0) = this->getTag();
  data(1) = numFibers;
  data(2) = (theTorsion != 0) ? 1 : 0;
  int dbTag = this->getDbTag();
  if (theTorsion != 0) {
    theTorsion->setDbTag(dbTag);
    data(3) = theTorsion->getClassTag();
  }
  data(4) = computeCentroid ? 1 : 0; // Now the ID data is really 5
  
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FiberSection3d::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (theTorsion != 0)
    theTorsion->sendSelf(commitTag, theChannel);

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
     opserr << "FiberSection3d::sendSelf - failed to send material data\n";
     return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 3*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSection3d::sendSelf - failed to send fiber data\n";
     return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);
  }

  return res;
}

int
FiberSection3d::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(5);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);

  if (res < 0) {
   opserr << "FiberSection3d::recvSelf - failed to recv ID data\n";
   return res;
  } 
   
  this->setTag(data(0));

  if (data(2) == 1 && theTorsion == 0) {	
    int cTag = data(3);
    theTorsion = theBroker.getNewUniaxialMaterial(cTag);
    if (theTorsion == 0) {
      opserr << "FiberSection3d::recvSelf - failed to get torsion material \n";
      return -1;
    }
    theTorsion->setDbTag(dbTag);
  }

  if (theTorsion->recvSelf(commitTag, theChannel, theBroker) < 0) {
	   opserr << "FiberSection3d::recvSelf - torsion failed to recvSelf \n";
       return -2;
  }
  
  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "FiberSection3d::recvSelf - failed to recv material data\n";
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
	  opserr << "FiberSection3d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}

	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;
	
	matData = new double [numFibers*3];

	if (matData == 0) {
	  opserr << "FiberSection3d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, 3*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSection3d::recvSelf - failed to recv fiber data\n";
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
	opserr << "FiberSection3d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    QzBar = 0.0;
    QyBar = 0.0;
    Abar  = 0.0;
    double yLoc, zLoc, Area;

    computeCentroid = data(4) ? true : false;
    
    // Recompute centroid
    for (i = 0; computeCentroid && i < numFibers; i++) {
      yLoc = matData[3*i];
      zLoc = matData[3*i+1];
      Area = matData[3*i+2];
      Abar  += Area;
      QzBar += yLoc*Area;
      QyBar += zLoc*Area;
    }

    if (computeCentroid) {
      yBar = QzBar/Abar;
      zBar = QyBar/Abar;
    } else {
      yBar = 0.0;
      zBar = 0.0;      
    }
  }    

  return res;
}

void
FiberSection3d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "\nFiberSection3d, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: (" << yBar << ", " << zBar << ')' << endln;
    if (theTorsion != 0)
        theTorsion->Print(s, flag);    

    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
      for (int i = 0; i < numFibers; i++) {
	s << "\nLocation (y, z) = (" << matData[3*i] << ", " << matData[3*i+1] << ")";
	s << "\nArea = " << matData[3*i+2] << endln;
	theMaterials[i]->Print(s, flag);
	
      }
    }
  }
  if (flag == 3) {
    for (int i = 0; i < numFibers; i++) {
      s << theMaterials[i]->getTag() << " " << matData[3*i] << " "  << matData[3*i+1] << " "  << matData[3*i+2] << " " ;
      s << theMaterials[i]->getStress() << " "  << theMaterials[i]->getStrain() << endln;
    } 
  }
    
  if (flag == 4) {
    for (int i = 0; i < numFibers; i++) {
      s << "add fiber # " << i+1 << " using material # " << theMaterials[i]->getTag() << " to section # 1\n";
      s << "fiber_cross_section = " << matData[3*i+2] << "*m^2\n";
      s << "fiber_location = (" << matData[3*i] << "*m, " << matData[3*i+1] << "*m);\n\n";
    }
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
	  s << "\t\t\t{";
	  s << "\"name\": \"" << this->getTag() << "\", ";
	  s << "\"type\": \"FiberSection3d\", ";
	  if (theTorsion != 0)
	    s << "\"torsion\": " << theTorsion->getInitialTangent() << ", ";
	  s << "\"fibers\": [\n";
	  for (int i = 0; i < numFibers; i++) {
		  s << "\t\t\t\t{\"coord\": [" << matData[3*i] << ", " << matData[3*i+1] << "], ";
		  s << "\"area\": " << matData[3*i+2] << ", ";
		  s << "\"material\": \"" << theMaterials[i]->getTag() << "\"";
		  if (i < numFibers - 1)
			  s << "},\n";
		  else
			  s << "}\n";
	  }
	  s << "\t\t\t]}";
  }
}

Response*
FiberSection3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  const ID &type = this->getType();
  int typeSize = this->getOrder();

  Response *theResponse =0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 1, this->getSectionDeformation());
  
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new MaterialResponse(this, 2, this->getStressResultant());
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","kappaZ");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","eps");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","gammaY");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","kappaY");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","gammaZ");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","theta");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }
    for (int j=0; j<typeSize; j++) {
      int code = type(j);
      switch (code){
      case SECTION_RESPONSE_MZ:
	output.tag("ResponseType","Mz");
	break;
      case SECTION_RESPONSE_P:
	output.tag("ResponseType","P");
	break;
      case SECTION_RESPONSE_VY:
	output.tag("ResponseType","Vy");
	break;
      case SECTION_RESPONSE_MY:
	output.tag("ResponseType","My");
	break;
      case SECTION_RESPONSE_VZ:
	output.tag("ResponseType","Vz");
	break;
      case SECTION_RESPONSE_T:
	output.tag("ResponseType","T");
	break;
      default:
	output.tag("ResponseType","Unknown");
      }
    }

    theResponse =  new MaterialResponse(this, 4, Vector(2*this->getOrder()));
  
  } else if (strcmp(argv[0],"fiberData") == 0) {
    int numData = numFibers*5;
    for (int j = 0; j < numFibers; j++) {
      output.tag("FiberOutput");
      output.attr("yLoc", matData[3*j]);
      output.attr("zLoc", matData[3*j+1]);
      output.attr("area", matData[3*j+2]);    
      output.tag("ResponseType","yCoord");
      output.tag("ResponseType","zCoord");
      output.tag("ResponseType","area");
      output.tag("ResponseType","stress");
      output.tag("ResponseType","strain");
      output.endTag();
    }
    Vector theResponseData(numData);
    theResponse = new MaterialResponse(this, 5, theResponseData);

  } else if ((strcmp(argv[0],"numFailedFiber") == 0) || 
	     (strcmp(argv[0],"numFiberFailed") == 0)) {
    int count = 0;
    theResponse = new MaterialResponse(this, 6, count);

  } else if ((strcmp(argv[0],"sectionFailed") == 0) ||
	     (strcmp(argv[0],"hasSectionFailed") == 0) ||
	     (strcmp(argv[0],"hasFailed") == 0)) {

    int count = 0;
    return theResponse = new MaterialResponse(this, 7, count);
  }
  //by SAJalali
  else if ((strcmp(argv[0], "energy") == 0) || (strcmp(argv[0], "Energy") == 0)) {
	  return theResponse = new MaterialResponse(this, 10, getEnergy());
  }


  else {
    if (argc > 2 || strcmp(argv[0],"fiber") == 0) {
      
      int key = numFibers;
      int passarg = 2;
      
      if (argc <= 3)	{  // fiber number was input directly
	
	key = atoi(argv[1]);
	
      } else if (argc > 4) {         // find fiber closest to coord. with mat tag
	int matTag = atoi(argv[3]);
	double yCoord = atof(argv[1]);
	double zCoord = atof(argv[2]);
	double closestDist = 0.0;
	double ySearch, zSearch, dy, dz;
	double distance;
	int j;
	
	// Find first fiber with specified material tag
	for (j = 0; j < numFibers; j++) {
	  if (matTag == theMaterials[j]->getTag()) {
	    ySearch = matData[3*j];
	    zSearch = matData[3*j+1];
	    dy = ySearch-yCoord;
	    dz = zSearch-zCoord;
	    closestDist = sqrt(dy*dy + dz*dz);
	    key = j;
	    break;
	  }
	}
	
	// Search the remaining fibers
	for ( ; j < numFibers; j++) {
	  if (matTag == theMaterials[j]->getTag()) {
	    ySearch = matData[3*j];
	    zSearch = matData[3*j+1];
	    dy = ySearch-yCoord;
	    dz = zSearch-zCoord;
	    distance = sqrt(dy*dy + dz*dz);
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
	double zCoord = atof(argv[2]);
	double closestDist;
	double ySearch, zSearch, dy, dz;
	double distance;
	ySearch = matData[0];
	zSearch = matData[1];
	dy = ySearch-yCoord;
	dz = zSearch-zCoord;
	closestDist = sqrt(dy*dy + dz*dz);
	key = 0;
	for (int j = 1; j < numFibers; j++) {
	  ySearch = matData[3*j];
	  zSearch = matData[3*j+1];
	  dy = ySearch-yCoord;
	  dz = zSearch-zCoord;
	  distance = sqrt(dy*dy + dz*dz);
	  if (distance < closestDist) {
	    closestDist = distance;
	    key = j;
	  }
	}
	passarg = 3;
      }
      
      if (key < numFibers && key >= 0) {
	output.tag("FiberOutput");
	output.attr("yLoc",matData[3*key]);
	output.attr("zLoc",matData[3*key+1]);
	output.attr("area",matData[3*key+2]);
	
	theResponse =  theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
	
	output.endTag();
      }
    }
  }

  output.endTag();
  return theResponse;
}


int 
FiberSection3d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  if (responseID == 5) {
    int numData = 5*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      double yLoc, zLoc, A, stress, strain;
      yLoc = matData[3*j];
      zLoc = matData[3*j+1];
      A = matData[3*j+2];
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
  else  if (responseID == 10) {

	  return sectInfo.setDouble(getEnergy());
  }

  return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int
FiberSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  // A material parameter
  if (strstr(argv[0],"material") != 0) {

    // Get the tag of the material
    int paramMatTag = atoi(argv[1]);

    // Loop over fibers to find the right material(s)
    int ok = 0;
    for (int i = 0; i < numFibers; i++)
      if (paramMatTag == theMaterials[i]->getTag()) {
	ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
	if (ok != -1)
	  result = ok;
      }
    
    if (paramMatTag == theTorsion->getTag()) {
	ok = theTorsion->setParameter(&argv[2], argc-2, param);
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
  
  // loop over every material
  for (int i = 0; i < numFibers; i++) {
    ok = theMaterials[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  // Don't really need to do this in "default" mode
  //ok = theTorsion->setParameter(argv, argc, param);
  //if (ok != -1)
  //  result = ok;

  if (sectionIntegr != 0) {
    ok = sectionIntegr->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  return result;
}

const Vector &
FiberSection3d::getSectionDeformationSensitivity(int gradIndex)
{
  static Vector dummy(4);
  
  dummy.Zero();
  
  return dummy;
}

   
const Vector &
FiberSection3d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(4);
  
  ds.Zero();
  
  double y, z, A;
  double stress = 0;
  double dsigdh = 0;
  double sig_dAdh = 0;
  double tangent = 0;

  static double yLocs[10000];
  static double zLocs[10000];
  static double fiberArea[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }  
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[3*i];
      zLocs[i] = matData[3*i+1];
      fiberArea[i] = matData[3*i+2];
    }
  }

  static double dydh[10000];
  static double dzdh[10000];
  static double areaDeriv[10000];

  if (sectionIntegr != 0) {
    sectionIntegr->getLocationsDeriv(numFibers, dydh, dzdh);  
    sectionIntegr->getWeightsDeriv(numFibers, areaDeriv);
  }
  else {
    for (int i = 0; i < numFibers; i++) {
      dydh[i] = 0.0;
      dzdh[i] = 0.0;
      areaDeriv[i] = 0.0;
    }
  }
  
  for (int i = 0; i < numFibers; i++) {
    y = yLocs[i] - yBar;
    z = zLocs[i] - zBar;
    A = fiberArea[i];
    
    dsigdh = theMaterials[i]->getStressSensitivity(gradIndex, conditional);

    ds(0) += dsigdh*A;
    ds(1) += -y*dsigdh*A;
    ds(2) +=  z*dsigdh*A;

    if (areaDeriv[i] != 0.0 || dydh[i] != 0.0 ||  dzdh[i] != 0.0)
      stress = theMaterials[i]->getStress();

    if (dydh[i] != 0.0 || dzdh[i] != 0.0)
      tangent = theMaterials[i]->getTangent();

    if (areaDeriv[i] != 0.0) {
      sig_dAdh = stress*areaDeriv[i];
      
      ds(0) += sig_dAdh;
      ds(1) += -y*sig_dAdh;
      ds(2) +=  z*sig_dAdh;
    }

    if (dydh[i] != 0.0)
      ds(1) += -dydh[i] * (stress*A);

    if (dzdh[i] != 0.0)
      ds(2) +=  dzdh[i] * (stress*A);

    static Matrix as(1,3);
    as(0,0) = 1;
    as(0,1) = -y;
    as(0,2) = z;
    
    static Matrix dasdh(1,3);
    dasdh(0,1) = -dydh[i];
    dasdh(0,2) = dzdh[i];
    
    static Matrix tmpMatrix(3,3);
    tmpMatrix.addMatrixTransposeProduct(0.0, as, dasdh, tangent);
    
    //ds.addMatrixVector(1.0, tmpMatrix, e, A);
    ds(0) += (tmpMatrix(0,0)*e(0) + tmpMatrix(0,1)*e(1) + tmpMatrix(0,2)*e(2))*A;
    ds(1) += (tmpMatrix(1,0)*e(0) + tmpMatrix(1,1)*e(1) + tmpMatrix(1,2)*e(2))*A;
    ds(2) += (tmpMatrix(2,0)*e(0) + tmpMatrix(2,1)*e(1) + tmpMatrix(2,2)*e(2))*A;
  }

  ds(3) = theTorsion->getStressSensitivity(gradIndex, conditional);

  return ds;
}

const Matrix &
FiberSection3d::getSectionTangentSensitivity(int gradIndex)
{
  static Matrix something(4,4);
  
  something.Zero();

  something(3,3) = theTorsion->getTangentSensitivity(gradIndex);
  
  return something;
}

int
FiberSection3d::commitSensitivity(const Vector& defSens, int gradIndex, int numGrads)
{

  double d0 = defSens(0);
  double d1 = defSens(1);
  double d2 = defSens(2);
  double d3 = defSens(3);

  //dedh = defSens;

  static double yLocs[10000];
  static double zLocs[10000];

  if (sectionIntegr != 0)
    sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
  else {
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[3*i];
      zLocs[i] = matData[3*i+1];
    }
  }

  static double dydh[10000];
  static double dzdh[10000];

  if (sectionIntegr != 0)
    sectionIntegr->getLocationsDeriv(numFibers, dydh, dzdh);  
  else {
    for (int i = 0; i < numFibers; i++) {
      dydh[i] = 0.0;
      dzdh[i] = 0.0;
    }
  }

  double y, z;

  double depsdh = 0;

  for (int i = 0; i < numFibers; i++) {
    y = yLocs[i] - yBar;
    z = zLocs[i] - zBar;

    // determine material strain and set it
    depsdh = d0 - y*d1 + z*d2 - dydh[i]*e(1) + dzdh[i]*e(2);

    theMaterials[i]->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  theTorsion->commitSensitivity(d3, gradIndex, numGrads);

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////


//by SAJalali
double FiberSection3d::getEnergy() const
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
