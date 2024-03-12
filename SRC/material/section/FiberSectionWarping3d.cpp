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
                                                                        
// $Revision: 1.31 $
// $Date: 2009/09/28 22:48:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionWarping3d.cpp,v $
                                                                        
// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSectionWarping3d.
// Modified by Xi Zhang from University of Sydney, Australia (include warping degrees of freedom). Refer to 
// Formulation and Implementation of Three-dimensional Doubly Symmetric Beam-Column Analyses with Warping Effects in OpenSees
// Research Report R917, School of Civil Engineering, University of Sydney.

#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSectionWarping3d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <SectionIntegration.h>
#include <elementAPI.h>
using std::string;
using namespace std;


ID FiberSectionWarping3d::code(6);

void* OPS_FiberSectionWarping3d()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	    opserr<<"insufficient arguments for FiberSectionWarping3d\n";
	    return 0;
    }
    
    numData = 1;
    int tag;
    if (OPS_GetIntInput(&numData, &tag) < 0) {
      opserr << "FiberSectionWarping3d - unable to read tag" << endln;
      return 0;
    }
    
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
    SectionForceDeformation *section = new FiberSectionWarping3d(tag, num, *torsion);
    if (deleteTorsion)
      delete torsion;
    return section;
}

// constructors:
FiberSectionWarping3d::FiberSectionWarping3d(int tag, int num, Fiber **fibers,
					     UniaxialMaterial &torsion): 
  SectionForceDeformation(tag, SEC_TAG_FiberSectionWarping3d),
  numFibers(num), sizeFibers(num),  theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(8), eCommit(8), s(0), ks(0), theTorsion(0)
{
  if (numFibers != 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSectionWarping3d::FiberSectionWarping3d -- failed to allocate Material pointers\n";
      exit(-1);
    }

    matData = new double [numFibers*4];

    if (matData == 0) {
      opserr << "FiberSectionWarping3d::FiberSectionWarping3d -- failed to allocate double array for material data\n";
      exit(-1);
    }

    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;
    double Heightt;
    
    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();
	  Heightt = theFiber->getd();

      Qz += yLoc*Area;
      Qy += zLoc*Area;
      A  += Area;

      matData[i*4] = yLoc;
      matData[i*4+1] = zLoc;
      matData[i*4+2] = Area;
	  matData[i*4+3] = Heightt;
      UniaxialMaterial *theMat = theFiber->getMaterial();
      theMaterials[i] = theMat->getCopy();

      if (theMaterials[i] == 0) {
	opserr << "FiberSectionWarping3d::FiberSectionWarping3d -- failed to get copy of a Material\n";
	exit(-1);
      }
    }

    yBar = -Qz/A;
    zBar = Qy/A;
  }

  theTorsion = torsion.getCopy();
  if (theTorsion == 0)
    opserr << "FiberSectionWarping3d::FiberSectionWarping3d -- failed to get copy of torsion material\n";
  
  s = new Vector(sData, 6);
  ks = new Matrix(kData, 6, 6);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0;
  sData[5] = 0.0;  

  for (int i=0; i<36; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_W;
  code(4) = SECTION_RESPONSE_B;  
  code(5) = SECTION_RESPONSE_T;  

 // AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
  SHVs=0;
  // AddingSensitivity:END //////////////////////////////////////
}

FiberSectionWarping3d::FiberSectionWarping3d(int tag, int num, UniaxialMaterial &torsion): 
    SectionForceDeformation(tag, SEC_TAG_FiberSectionWarping3d),
    numFibers(0), sizeFibers(num), theMaterials(0), matData(0),
    yBar(0.0), zBar(0.0),
    e(8), eCommit(8), s(0), ks(0), theTorsion(0)
{
    if(sizeFibers != 0) {
	theMaterials = new UniaxialMaterial *[sizeFibers];

	if (theMaterials == 0) {
	    opserr << "FiberSection3d::FiberSection3d -- failed to allocate Material pointers\n";
	    exit(-1);
	}

	matData = new double [sizeFibers*4];

	if (matData == 0) {
	    opserr << "FiberSection3d::FiberSection3d -- failed to allocate double array for material data\n";
	    exit(-1);
	}

	for (int i = 0; i < sizeFibers; i++) {
	    matData[i*4] = 0.0;
	    matData[i*4+1] = 0.0;
	    matData[i*4+2] = 0.0;
	    matData[i*4+3] = 0.0;	    
	    theMaterials[i] = 0;
	}
    }

    theTorsion = torsion.getCopy();
    if (theTorsion == 0) 
      opserr << "FiberSection3d::FiberSection3d -- failed to get copy of torsion material\n";

    s = new Vector(sData, 6);
    ks = new Matrix(kData, 6, 6);

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;
    sData[3] = 0.0;
    sData[4] = 0.0;
    sData[5] = 0.0;        

    for (int i=0; i<36; i++)
	kData[i] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_W;
    code(4) = SECTION_RESPONSE_B;  
    code(5) = SECTION_RESPONSE_T;    
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSectionWarping3d::FiberSectionWarping3d():
  SectionForceDeformation(0, SEC_TAG_FiberSectionWarping3d),
  numFibers(0), theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(8), eCommit(8), s(0), ks(0), theTorsion(0)
{
  s = new Vector(sData, 6);
  ks = new Matrix(kData, 6, 6);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0;
  sData[5] = 0.0;  

  for (int i=0; i<36; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_W;
  code(4) = SECTION_RESPONSE_B;  
  code(5) = SECTION_RESPONSE_T;  

 // AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
  SHVs=0;
  // AddingSensitivity:END //////////////////////////////////////
}

int
FiberSectionWarping3d::addFiber(Fiber &newFiber)
{
  // need to create a larger array
  int newSize = numFibers+1;

  UniaxialMaterial **newArray = new UniaxialMaterial *[newSize]; 
  double *newMatData = new double [4 * newSize];
  
  if (newArray == 0 || newMatData == 0) {
    opserr << "FiberSectionWarping3d::addFiber -- failed to allocate Fiber pointers\n";
    exit(-1);
  }

  // copy the old pointers
  int i;
  for (i = 0; i < numFibers; i++) {
    newArray[i] = theMaterials[i];
    newMatData[4*i] = matData[4*i];
    newMatData[4*i+1] = matData[4*i+1];
    newMatData[4*i+2] = matData[4*i+2];
    newMatData[4*i+3] = matData[4*i+3];
  }
  // set the new pointers
  double yLoc, zLoc, Area, Height;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  Height = newFiber.getd();
  newMatData[numFibers*4] = -yLoc;
  newMatData[numFibers*4+1] = zLoc;
  newMatData[numFibers*4+2] = Area;
  newMatData[numFibers*4+3] = Height;
  UniaxialMaterial *theMat = newFiber.getMaterial();
  newArray[numFibers] = theMat->getCopy();

  if (newArray[numFibers] == 0) {
    opserr << "FiberSectionWarping3d::addFiber -- failed to get copy of a Material\n";
    exit(-1);

    delete [] newArray;
    delete [] newMatData;
    return -1;
  }

  numFibers++;
  
  if (theMaterials != 0) {
    delete [] theMaterials;
    delete [] matData;
  }

  theMaterials = newArray;
  matData = newMatData;

  double Qz = 0.0;
  double Qy = 0.0;
  double A  = 0.0;

  // Recompute centroid
  for (i = 0; i < numFibers; i++) {
    yLoc = -matData[4*i];
    zLoc = matData[4*i+1];
    Area = matData[4*i+2];
    Height = matData[4*i+3];
    A  += Area;
    Qz += yLoc*Area;
    Qy += zLoc*Area;
  }

  yBar = -Qz/A;
  zBar = Qy/A;

  return 0;
}



// destructor:
FiberSectionWarping3d::~FiberSectionWarping3d()
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

  if (theTorsion != 0)
    delete theTorsion;  
}

int
FiberSectionWarping3d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;

  for (int i = 0; i < 36; i++)
    kData[i] = 0.0;

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0;
  sData[5] = 0.0;  

  int loc = 0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);
  double d3 = deforms(3);
  double d4 = deforms(4);
  double d5 = deforms(5);
  double d6 = deforms(6);
  double d7 = deforms(7);
  double d8 = 0.0; // Torsion?? -- MHS
  
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
    double Height=matData[loc++];
	
    // calculate sectorial area
    double omig=0.0;
    if (y>0.0)
      omig = -z*(y-Height);
    else
      omig = -z* (y+Height);
    

    // determine material strain and set it, include second order terms
    double strain = d0 - y*d1 - z*d2 - omig*d3 + 0.5*d5*d5 + 0.5*d6*d6 + 0.5*(y*y+z*z)*d4*d4 - y*d7*d2 + z*d7*d1;
    double tangent, stress;
    res += theMat->setTrial(strain, stress, tangent);

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    // section stiffness matrix k, refer to Alemdar
    kData[0] += value;
    kData[3] += (y*y+z*z)*value;
    kData[6] += vas1 * y;
    kData[12] += vas2 * z; 
    kData[15] += (y*y+z*z)*value;
    kData[18] += (y*y+z*z)*(y*y+z*z)*value;
    kData[24] += omig*omig*value;
    
    // section force vector D, refer to Alemdar
    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += -1.0 * fs0 * y;
    sData[2] += -1.0 * fs0 * z;
    sData[3] += fs0 * (y*y+z*z);
    sData[4] += -fs0 * omig;
  }

  if (theTorsion != 0) {
    double stress, tangent;
    res += theTorsion->setTrial(d8, stress, tangent);
    sData[5] = stress;
    kData[35] = tangent;
  }
  
  return res;
}

const Matrix&
FiberSectionWarping3d::getInitialTangent(void)
{
  static double kInitialData[36];
  static Matrix kInitial(kInitialData, 6, 6);
  for (int i=0; i<36; i++)
    kInitialData[i]=0.0;

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
    double Height = matData[loc++];
    // calculate sectorial area
    double omig;
    if (y>0.0)
      omig = -z*(y-Height);
    else
      omig = -z* (y+Height);
    
    double tangent = theMat->getInitialTangent();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    // section stiffness matrix k, refer to Alemdar
    kInitialData[0] += value;
    kInitialData[3] += (y*y+z*z)*value;
    kInitialData[6] += vas1 * y;
    kInitialData[12] += vas2 * z; 
    kInitialData[15] += (y*y+z*z)*value;
    kInitialData[18] += (y*y+z*z)*(y*y+z*z)*value;
    kInitialData[24] += omig*omig*value;
  }

  if (theTorsion != 0)
    kInitialData[35] = theTorsion->getInitialTangent();

  return kInitial;
}

const Vector&
FiberSectionWarping3d::getSectionDeformation(void)
{
  return e;
}

const Matrix&
FiberSectionWarping3d::getSectionTangent(void)
{
  return *ks;
}

const Vector&
FiberSectionWarping3d::getStressResultant(void)
{
  return *s;
}

SectionForceDeformation*
FiberSectionWarping3d::getCopy(void)
{
  FiberSectionWarping3d *theCopy = new FiberSectionWarping3d ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "FiberSectionWarping3d::FiberSectionWarping3d -- failed to allocate Material pointers\n";
      exit(-1);			    
    }

    theCopy->matData = new double [numFibers*4];

    if (theCopy->matData == 0) {
      opserr << "FiberSectionWarping3d::FiberSectionWarping3d -- failed to allocate double array for material data\n";
      exit(-1);
    }
			    
    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*4] = matData[i*4];
      theCopy->matData[i*4+1] = matData[i*4+1];
      theCopy->matData[i*4+2] = matData[i*4+2];
      theCopy->matData[i*4+3] = matData[i*4+3];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr << "FiberSectionWarping3d::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    
  }

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;

  for (int i=0; i<36; i++)
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];
  theCopy->sData[3] = sData[3];
  theCopy->sData[4] = sData[4];
  theCopy->sData[5] = sData[5];  

  if (theTorsion != 0)
    theCopy->theTorsion = theTorsion->getCopy();
  else
    theCopy->theTorsion = 0;
  
  return theCopy;
}

const ID&
FiberSectionWarping3d::getType ()
{
  return code;
}

int
FiberSectionWarping3d::getOrder () const
{
  return 5;
}

int
FiberSectionWarping3d::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  if (theTorsion != 0)
    err += theTorsion->commitState();
  
  eCommit = e;

  return err;
}

int
FiberSectionWarping3d::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;

  for (int i = 0; i < 36; i++)
    kData[i] = 0.0;

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0;
  sData[5] = 0.0;

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
    double Height = matData[loc++];

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;
    double omig;
    if (y>0.0)
      omig = -z*(y-Height);
    else
      omig = -z* (y+Height);
    
    kData[0] += value;
    kData[3] += (y*y+z*z)*value;
    kData[6] += vas1 * y;
    kData[12] += vas2 * z; 
    kData[15] += (y*y+z*z)*value;
    kData[18] += (y*y+z*z)*(y*y+z*z)*value;
    kData[24] += omig*omig*value;
    
    double fs0 = stress * A;
    
    sData[0] += fs0;
    sData[1] += -1.0 * fs0 * y;
    sData[2] += -1.0 * fs0 * z;
    sData[3] += fs0 * (y*y+z*z);
    sData[4] += -fs0 * omig;
  }

  if (theTorsion != 0) {
    err += theTorsion->revertToLastCommit();
    sData[5] = theTorsion->getStress();
    kData[35] = theTorsion->getTangent();
  } else {
    sData[5] = 0.0;
    kData[35] = 0.0;
  }
  
  return err;
}

int
FiberSectionWarping3d::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;

  for (int i = 0; i < 36; i++)
    kData[i] = 0.0;

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0;
  sData[5] = 0.0;  
  
  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
    double Height = matData[loc++];
    double omig;
    if (y>0.0)
      omig = -z*(y-Height);
    else
      omig = -z* (y+Height);
    
    // invoke revertToStart on the material
    err += theMat->revertToStart();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[3] += (y*y+z*z)*value;
    kData[6] += vas1 * y;
    kData[12] += vas2 * z; 
    kData[15] += (y*y+z*z)*value;
    kData[18] += (y*y+z*z)*(y*y+z*z)*value;
    kData[24] += omig*omig*value;
    
    double fs0 = stress * A;

    sData[0] += fs0;
    sData[1] += -1.0 * fs0 * y;
    sData[2] += -1.0 * fs0 * z;
    sData[3] += fs0 * (y*y+z*z);
    sData[4] += -fs0 * omig;
  }

  if (theTorsion != 0) {
    err += theTorsion->revertToStart();
    kData[35] = theTorsion->getTangent();
    sData[5] = theTorsion->getStress();
  } else {
    kData[35] = 0.0;
    sData[5] = 0.0;
  }  

  return err;
}

int
FiberSectionWarping3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 5 so no conflict with matData below if just 1 fiber
  static ID data(5);
  data(0) = this->getTag();
  data(1) = numFibers;
  data(2) = (theTorsion != 0) ? 1 : 0;  
  int dbTag = this->getDbTag();
  if (theTorsion != 0) {
    theTorsion->setDbTag(dbTag);
    data(3) = theTorsion->getClassTag();
  }  
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FiberSectionWarping3d::sendSelf - failed to send ID data\n";
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
     opserr << "FiberSectionWarping3d::sendSelf - failed to send material data\n";
     return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 4*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSectionWarping3d::sendSelf - failed to send material data\n";
     return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);
  }

  return res;
}

int
FiberSectionWarping3d::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(5);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);

  if (res < 0) {
   opserr << "FiberSectionWarping3d::sendSelf - failed to recv ID data\n";
   return res;
  } 
   
  this->setTag(data(0));

  if (data(2) == 1 && theTorsion == 0) {	
    int cTag = data(3);
    theTorsion = theBroker.getNewUniaxialMaterial(cTag);
    if (theTorsion == 0) {
      opserr << "FiberSectionWarping3d::recvSelf - failed to get torsion material \n";
      return -1;
    }
    theTorsion->setDbTag(dbTag);
  }

  if (theTorsion->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "FiberSectionWarping3d::recvSelf - torsion failed to recvSelf \n";
    return -2;
  }
  
  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "FiberSectionWarping3d::sendSelf - failed to send material data\n";
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
	  opserr << "FiberSectionWarping3d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}

	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;
	
	matData = new double [numFibers*4];

	if (matData == 0) {
	  opserr << "FiberSectionWarping3d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, 4*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSectionWarping3d::sendSelf - failed to send material data\n";
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
	opserr << "FiberSectionWarping3d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;
    double yLoc, zLoc, Area;//, Height;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = -matData[4*i];
      zLoc = matData[4*i+1];
      Area = matData[4*i+2];
      //Height=matData[4*i+3];
      A  += Area;
      Qz += yLoc*Area;
      Qy += zLoc*Area;
    }
    
    yBar = -Qz/A;
    zBar = Qy/A;
  }    

  return res;
}

void
FiberSectionWarping3d::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {
    for (int i = 0; i < numFibers; i++) {
      s << -matData[4*i] << " "  << matData[4*i+1] << " "  << matData[4*i+2] << " " ;
      s << theMaterials[i]->getStress() << " "  << theMaterials[i]->getStrain() << endln;
    } 
  } else {
    s << "\nFiberSectionWarping3d, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: (" << -yBar << ", " << zBar << ')' << endln;
    
    if (flag == 1) {
      for (int i = 0; i < numFibers; i++) {
	s << "\nLocation (y, z) = (" << -matData[4*i] << ", " << matData[4*i+1] << ")";
	s << "\nArea = " << matData[4*i+2] << endln;
      theMaterials[i]->Print(s, flag);
      }
    }
  }
}

Response*
FiberSectionWarping3d::setResponse(const char **argv, int argc, OPS_Stream &output)
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
	double closestDist;
	double ySearch, zSearch, dy, dz;
	double distance;
	int j;
	
	// Find first fiber with specified material tag
	for (j = 0; j < numFibers; j++) {
	  if (matTag == theMaterials[j]->getTag()) {
	    ySearch = -matData[4*j];
	    zSearch =  matData[4*j+1];
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
	    ySearch = -matData[4*j];
	    zSearch =  matData[4*j+1];
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
	ySearch = -matData[0];
	zSearch =  matData[1];
	dy = ySearch-yCoord;
	dz = zSearch-zCoord;
	closestDist = sqrt(dy*dy + dz*dz);
	key = 0;
	for (int j = 1; j < numFibers; j++) {
	  ySearch = -matData[4*j];
	  zSearch =  matData[4*j+1];
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
	output.attr("yLoc",-matData[4*key]);
	output.attr("zLoc",matData[4*key+1]);
	output.attr("area",matData[4*key+2]);
	
	theResponse =  theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
	
	output.endTag();
      }
    }
  }

  output.endTag();
  return theResponse;
}


int 
FiberSectionWarping3d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int
FiberSectionWarping3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 3)
    return -1;


  int result = 0;

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
    
    return result;
  }    

  int ok = 0;
  
  // loop over every material
  for (int i = 0; i < numFibers; i++) {
    ok = theMaterials[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  return result;
}

const Vector &
FiberSectionWarping3d::getSectionDeformationSensitivity(int gradIndex)
{
  static Vector dummy(3);
  dummy.Zero();
  if (SHVs !=0) {
    dummy(0) = (*SHVs)(0,gradIndex);
    dummy(1) = (*SHVs)(1,gradIndex);
    dummy(2) = (*SHVs)(2,gradIndex);
  }
  return dummy;
}

   
const Vector &
FiberSectionWarping3d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  
  static Vector ds(3);
  
  ds.Zero();
  
  double  stressGradient;
  int loc = 0;
  
  
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
    stressGradient = theMaterials[i]->getStressSensitivity(gradIndex,conditional);
    stressGradient *=  A;
    ds(0) += stressGradient;
    ds(1) += stressGradient * y;
    ds(2) += stressGradient * z;
    
  } 
  
  return ds;
}

const Matrix &
FiberSectionWarping3d::getSectionTangentSensitivity(int gradIndex)
{
  static Matrix something(2,2);
  
  something.Zero();
  
  return something;
}

int
FiberSectionWarping3d::commitSensitivity(const Vector& defSens, int gradIndex, int numGrads)
{

  // here add SHVs to store the strain sensitivity.

  if (SHVs == 0) {
    SHVs = new Matrix(3,numGrads);
  }
  
  (*SHVs)(0,gradIndex) = defSens(0);
  (*SHVs)(1,gradIndex) = defSens(1);
  (*SHVs)(2,gradIndex) = defSens(2);

  int loc = 0;

  double d0 = defSens(0);
  double d1 = defSens(1);
  double d2 = defSens(2);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    loc++;   // skip A data.
    
    double strainSens = d0 + y*d1 + z*d2;
    
    theMat->commitSensitivity(strainSens,gradIndex,numGrads);
  }

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////


