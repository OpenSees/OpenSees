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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionAsym3d.cpp,v $
                                                                        
// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSection2d.

// Modified by: Xinlong Du and Jerome F. Hajjar, Northeastern University, USA; Year 2019
// Description: Modified FiberSection3d.cpp (from version 3.0.0 on 11/5/2018) 
//              to include shear center coordinates and high-order longitudinal strain terms.
// References:
// Du, X., & Hajjar, J. (2021). Three-dimensional nonlinear displacement-based beam element
// for members with angle and tee sections. Engineering Structures, 239, 112239.
// Du, X., & Hajjar, J. F. (2021). Three-dimensional nonlinear mixed 6-DOF beam element 
// for thin-walled members. Thin-Walled Structures, 164, 107817.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSectionAsym3d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include <SectionIntegration.h>
#include <elementAPI.h>
#include <string.h>

ID FiberSectionAsym3d::code(5);

void* OPS_FiberSectionAsym3d()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	    opserr<<"insufficient arguments for FiberSectionAsym3d\n";
	    return 0;
    }
    
    numData = 1;
    int tag;
    if (OPS_GetIntInput(&numData, &tag) < 0) {
      opserr << "FiberSectionAsym3d - unable to read tag" << endln;
      return 0;
    }
    
    numData = 2;
    double dData[2];
    if (OPS_GetDoubleInput(&numData, dData) < 0) {
      opserr << "FiberSectionAsym3d - unable to read shear center data" << endln;
      return 0;
    }

    double GJ = 0.0;
    UniaxialMaterial *torsion = 0;
    bool deleteTorsion = false;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();
        if (strcmp(opt, "-GJ") == 0  && OPS_GetNumRemainingInputArgs() > 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &GJ) < 0) {
	      opserr << "WARNING: failed to read GJ" << endln;
	      return 0;
	    }
            torsion = new ElasticMaterial(0, GJ);
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
    
    int num = 30;
    SectionForceDeformation* section = new FiberSectionAsym3d(tag, num, torsion, dData[0], dData[1]);
    if (deleteTorsion)
        delete torsion;
    return section;
}

// constructors:
FiberSectionAsym3d::FiberSectionAsym3d(int tag, int num, Fiber **fibers, UniaxialMaterial *torsion, double yss, double zss):
  SectionForceDeformation(tag, SEC_TAG_FiberSectionAsym3d),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), sectionIntegr(0), e(5), s(0), ks(0), theTorsion(0), ys(yss), zs(zss)   //Xinlong
{
  if (numFibers != 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to allocate Material pointers\n";
      exit(-1);
    }

    matData = new double [numFibers*3];

    if (matData == 0) {
      opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to allocate double array for material data\n";
      exit(-1);
    }

    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
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
	opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to get copy of a Material\n";
	exit(-1);
      }
    }

    yBar = QzBar/Abar;
    zBar = QyBar/Abar;
  }

  if (torsion != 0) {
    theTorsion = torsion->getCopy();
    if (theTorsion == 0) {
      opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to get copy of torsion material\n";
    }
  }

  s = new Vector(sData, 5);              //Xinlong
  ks = new Matrix(kData, 5, 5);          //Xinlong

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0; //Xinlong

  for (int i=0; i<25; i++)      //Xinlong
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
  code(4) = SECTION_RESPONSE_W;
}

FiberSectionAsym3d::FiberSectionAsym3d(int tag, int num, UniaxialMaterial *torsion, double yss, double zss):    //Xinlong 
    SectionForceDeformation(tag, SEC_TAG_FiberSectionAsym3d),
    numFibers(0), sizeFibers(num), theMaterials(0), matData(0),
    QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), sectionIntegr(0), e(5), s(0), ks(0), theTorsion(0), ys(yss), zs(zss)  //Xinlong
{
    if(sizeFibers != 0) {
	theMaterials = new UniaxialMaterial *[sizeFibers];

	if (theMaterials == 0) {
	    opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to allocate Material pointers\n";
	    exit(-1);
	}

	matData = new double [sizeFibers*3];

	if (matData == 0) {
	    opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to allocate double array for material data\n";
	    exit(-1);
	}

	for (int i = 0; i < sizeFibers; i++) {
	    matData[i*3] = 0.0;
	    matData[i*3+1] = 0.0;
	    matData[i*3+2] = .0;
	    theMaterials[i] = 0;
	}
    }

    if (torsion != 0) {
      theTorsion = torsion->getCopy();
      if (theTorsion == 0) {
	opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to get copy of torsion material\n";
      }
    }

    s = new Vector(sData, 5);                //Xinlong
    ks = new Matrix(kData, 5, 5);            //Xinlong

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;
    sData[3] = 0.0;
	sData[4] = 0.0; //Xinlong

    for (int i=0; i<25; i++)  //Xinlong
	kData[i] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_T;
    code(4) = SECTION_RESPONSE_W;
}

FiberSectionAsym3d::FiberSectionAsym3d(int tag, int num, UniaxialMaterial **mats,
			       SectionIntegration &si, UniaxialMaterial *torsion, double yss, double zss):                   //Xinlong
  SectionForceDeformation(tag, SEC_TAG_FiberSectionAsym3d),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), sectionIntegr(0), e(5), s(0), ks(0), theTorsion(0), ys(yss), zs(zss)    //Xinlong
{
  if (numFibers != 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to allocate Material pointers";
      exit(-1);
    }
    matData = new double [numFibers*3];

    if (matData == 0) {
      opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to allocate double array for material data\n";
      exit(-1);
    }
  }

  sectionIntegr = si.getCopy();
  if (sectionIntegr == 0) {
    opserr << "Error: FiberSectionAsym3d::FiberSectionAsym3d: could not create copy of section integration object" << endln;
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
      opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to get copy of a Material\n";
      exit(-1);
    }
  }    
  
  yBar = QzBar/Abar;  
  zBar = QyBar/Abar;  

  if (torsion != 0) {
    theTorsion = torsion->getCopy();
    if (theTorsion == 0) {
      opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to get copy of torsion material\n";
    }
  }

  s = new Vector(sData, 5);      //Xinlong
  ks = new Matrix(kData, 5, 5);  //Xinlong
  
  for (int i = 0; i < 5; i++)    //Xinlong
    sData[i] = 0.0;

  for (int i = 0; i < 25; i++)   //Xinlong
    kData[i] = 0.0;
  
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
  code(4) = SECTION_RESPONSE_W;
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSectionAsym3d::FiberSectionAsym3d():
  SectionForceDeformation(0, SEC_TAG_FiberSectionAsym3d),
  numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), Abar(0.0), yBar(0.0), zBar(0.0), sectionIntegr(0), e(5), s(0), ks(0), theTorsion(0), ys(0.0), zs(0.0)     //Xinlong
{
  s = new Vector(sData, 5);     //Xinlong
  ks = new Matrix(kData, 5, 5); //Xinlong

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;
  sData[3] = 0.0;
  sData[4] = 0.0; //Xinlong

  for (int i=0; i<25; i++)  //Xinlong
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
  code(4) = SECTION_RESPONSE_W;
}

int
FiberSectionAsym3d::addFiber(Fiber &newFiber)
{
  // need to create a larger array
  if(numFibers == sizeFibers) {
      int newSize = 2*sizeFibers;
      UniaxialMaterial **newArray = new UniaxialMaterial *[newSize]; 
      double *newMatData = new double [3 * newSize];
      
      if (newArray == 0 || newMatData == 0) {
	  opserr << "FiberSectionAsym3d::addFiber -- failed to allocate Fiber pointers\n";
	  exit(-1);
      }

      // copy the old pointers
      for (int i = 0; i < numFibers; i++) {
	  newArray[i] = theMaterials[i];
	  newMatData[3*i] = matData[3*i];
	  newMatData[3*i+1] = matData[3*i+1];
	  newMatData[3*i+2] = matData[3*i+2];
      }

      // initialize new memory
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
    opserr << "FiberSectionAsym3d::addFiber -- failed to get copy of a Material\n";
    return -1;
  }

  numFibers++;

  // Recompute centroid
  Abar  += Area;
  QzBar += yLoc*Area;
  QyBar += zLoc*Area;

  yBar = QzBar/Abar;
  zBar = QyBar/Abar;

  return 0;
}



// destructor:
FiberSectionAsym3d::~FiberSectionAsym3d()
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
FiberSectionAsym3d::setTrialSectionDeformation (const Vector &deforms)
{
  //double zs = 0.6385; //z coord of shear center w.r.t. centroid
  //double ys = -0.6741; //y coord of shear center w.r.t. centroid

  int res = 0;
  e = deforms;
 
  for (int i = 0; i < 5; i++)  //Xinlong
    sData[i] = 0.0;
  for (int i = 0; i < 25; i++) //Xinlong
    kData[i] = 0.0;
  /*
  double d0 = deforms(0);  //u'   Xinlong
  double d1 = deforms(1);  //v"   Xinlong
  double d2 = deforms(2);  //w"   Xinlong
  double d3 = deforms(3);  //phi' Xinlong
  double d4 = deforms(4);  //v'   Xinlong
  double d5 = deforms(5);  //w'   Xinlong
  double d6 = deforms(6);  //phi  Xinlong
  double d7 = deforms(7);  //theta_Iz
  double d8 = deforms(8);  //theta_Iy
  double d9 = deforms(9);  //theta_Jz
  double d10 = deforms(10); //theta_Jy
  */
  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);
  double d3 = deforms(3);
  double d4 = deforms(4); //Phi'

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
    UniaxialMaterial *theMat = theMaterials[i];
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];

    // determine material strain and set it
	double pSquare = (y - ys)*(y - ys) + (z - zs)*(z - zs);
	//double strain = d0 - y * d1 - z * d2 + (4.0*d7*d7+4.0*d8*d8+4.0*d9*d9+4.0*d10*d10-2.0*d7*d9-2.0*d8*d10)/60.0 + 0.5*pSquare*d3*d3 + (zs*d4 - ys * d5)*d3 + (z*d1 - y * d2)*d6;
	double strain = d0 - y * d1 + z * d2 + pSquare * d3;
	res += theMat->setTrial(strain, stress, tangent);

    double value = tangent * A; //EA
    double vas1 = -y*value;     //-yEA
    double vas2 = z*value;      //zEA
    double vas1as2 = vas1*z;    //-yzEA

    kData[0] += value;           //EA
    kData[1] += vas1;            //-yEA
    kData[2] += vas2;            //zEA
	kData[3] += pSquare * value; //p2EA
    
    kData[6] += vas1 * -y;       //y2EA
    kData[7] += vas1as2;         //-yzEA
	kData[8] += pSquare * vas1;  //-p2yEA

	kData[12] += vas2 * z;       //z2EA
	kData[13] += pSquare * vas2; //p2zEA

	kData[18] += pSquare * pSquare*value; //p4EA
    
    double fs0 = stress * A; //sigma*A

    sData[0] += fs0;           //N
    sData[1] += fs0 * -y;      //Mz
    sData[2] += fs0 * z;       //My
	sData[3] += fs0 * pSquare; //W
  }

  kData[5] = kData[1];
  kData[10] = kData[2];
  kData[15] = kData[3];
  kData[11] = kData[7];
  kData[16] = kData[8];
  kData[17] = kData[13];

  if (theTorsion != 0) {
      res += theTorsion->setTrial(d4, stress, tangent);
      sData[4] = stress;   //T
      kData[24] = tangent; //GJ
  }

  return res;
}

const Matrix&
FiberSectionAsym3d::getInitialTangent(void)
{
  //double zs = 0.6385; //z coord of shear center w.r.t. centroid
  //double ys = -0.6741; //y coord of shear center w.r.t. centroid
  static double kInitialData[25];
  static Matrix kInitial(kInitialData, 5, 5);
  
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
    UniaxialMaterial *theMat = theMaterials[i];
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;
    double A = fiberArea[i];
	double pSquare = (y - ys)*(y - ys) + (z - zs)*(z - zs);

    double tangent = theMat->getInitialTangent();

    double value = tangent * A;
    double vas1 = -y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

	kInitialData[0] += value;           //EA
	kInitialData[1] += vas1;            //-yEA
	kInitialData[2] += vas2;            //zEA
	kInitialData[3] += pSquare * value; //p2EA

	kInitialData[6] += vas1 * -y;       //y2EA
	kInitialData[7] += vas1as2;         //-yzEA
	kInitialData[8] += pSquare * vas1;  //-p2yEA

	kInitialData[12] += vas2 * z;       //z2EA
	kInitialData[13] += pSquare * vas2; //p2zEA

	kInitialData[18] += pSquare * pSquare*value; //p4EA
  }

  kInitialData[5] = kInitialData[1];
  kInitialData[10] = kInitialData[2];
  kInitialData[15] = kInitialData[3];
  kInitialData[11] = kInitialData[7];
  kInitialData[16] = kInitialData[8];
  kInitialData[17] = kInitialData[13];

  if (theTorsion != 0)
      kInitialData[24] = theTorsion->getInitialTangent();

  return kInitial;
}

const Vector&
FiberSectionAsym3d::getSectionDeformation(void)
{
  return e;
}

const Matrix&
FiberSectionAsym3d::getSectionTangent(void)
{
  return *ks;
}

const Vector&
FiberSectionAsym3d::getStressResultant(void)
{
  return *s;
}

SectionForceDeformation*
FiberSectionAsym3d::getCopy(void)
{
  FiberSectionAsym3d *theCopy = new FiberSectionAsym3d ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;
  theCopy->sizeFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to allocate Material pointers\n";
      exit(-1);			    
    }

    theCopy->matData = new double [numFibers*3];

    if (theCopy->matData == 0) {
      opserr << "FiberSectionAsym3d::FiberSectionAsym3d -- failed to allocate double array for material data\n";
      exit(-1);
    }
			    
    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*3] = matData[i*3];
      theCopy->matData[i*3+1] = matData[i*3+1];
      theCopy->matData[i*3+2] = matData[i*3+2];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr << "FiberSectionAsym3d::getCopy -- failed to get copy of a Material\n";
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
  theCopy->ys = ys; //Xinlong 11/20/2019
  theCopy->zs = zs; //Xinlong 11/20/2019

  for (int i=0; i<25; i++)       //Xinlong
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];
  theCopy->sData[3] = sData[3];
  theCopy->sData[4] = sData[4]; //Xinlong

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
FiberSectionAsym3d::getType ()
{
  return code;
}

int
FiberSectionAsym3d::getOrder () const
{
  return 5;
}

int
FiberSectionAsym3d::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  if (theTorsion != 0)
      err += theTorsion->commitState();

  return err;
}

int
FiberSectionAsym3d::revertToLastCommit(void)
{
  //double zs = 0.6385; //z coord of shear center w.r.t. centroid
  //double ys = -0.6741; //y coord of shear center w.r.t. centroid

  int err = 0;

  for (int i = 0; i < 5; i++)  //Xinlong
	  sData[i] = 0.0;
  for (int i = 0; i < 25; i++) //Xinlong
	  kData[i] = 0.0;

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

	double pSquare = (y - ys)*(y - ys) + (z - zs)*(z - zs);

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = -y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;


	kData[0] += value;           //EA
	kData[1] += vas1;            //-yEA
	kData[2] += vas2;            //zEA
	kData[3] += pSquare * value; //p2EA

	kData[6] += vas1 * -y;       //y2EA
	kData[7] += vas1as2;         //-yzEA
	kData[8] += pSquare * vas1;  //-p2yEA

	kData[12] += vas2 * z;       //z2EA
	kData[13] += pSquare * vas2; //p2zEA

	kData[18] += pSquare * pSquare*value; //p4EA

    double fs0 = stress * A;

	sData[0] += fs0;           //N
	sData[1] += fs0 * -y;      //Mz
	sData[2] += fs0 * z;       //My
	sData[3] += fs0 * pSquare; //W
  }

  kData[5] = kData[1];
  kData[10] = kData[2];
  kData[15] = kData[3];
  kData[11] = kData[7];
  kData[16] = kData[8];
  kData[17] = kData[13];

  if (theTorsion != 0) {
      err += theTorsion->revertToLastCommit();
      kData[24] = theTorsion->getTangent();
  }
  else
      kData[24] = 0.0;
  //why do not have sData[4] here?
  return err;
}

int
FiberSectionAsym3d::revertToStart(void)
{
  //double zs = 0.6385; //z coord of shear center w.r.t. centroid
  //double ys = -0.6741; //y coord of shear center w.r.t. centroid

  // revert the fibers to start    
  int err = 0;

  for (int i = 0; i < 5; i++)  //Xinlong
	  sData[i] = 0.0;
  for (int i = 0; i < 25; i++) //Xinlong
	  kData[i] = 0.0;

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

	double pSquare = (y - ys)*(y - ys) + (z - zs)*(z - zs);

    // invoke revertToStart on the material
    err += theMat->revertToStart();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = -y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

	kData[0] += value;           //EA
	kData[1] += vas1;            //-yEA
	kData[2] += vas2;            //zEA
	kData[3] += pSquare * value; //p2EA

	kData[6] += vas1 * -y;       //y2EA
	kData[7] += vas1as2;         //-yzEA
	kData[8] += pSquare * vas1;  //-p2yEA

	kData[12] += vas2 * z;       //z2EA
	kData[13] += pSquare * vas2; //p2zEA

	kData[18] += pSquare * pSquare*value; //p4EA

    double fs0 = stress * A;

	sData[0] += fs0;           //N
	sData[1] += fs0 * -y;      //Mz
	sData[2] += fs0 * z;       //My
	sData[3] += fs0 * pSquare; //W
  }

  kData[5] = kData[1];
  kData[10] = kData[2];
  kData[15] = kData[3];
  kData[11] = kData[7];
  kData[16] = kData[8];
  kData[17] = kData[13];

  if (theTorsion != 0) {
      err += theTorsion->revertToStart();
      kData[24] = theTorsion->getTangent();
      sData[4] = theTorsion->getStress();
  }
  else {
      kData[24] = 0.0;
      sData[4] = 0.0;
  }

  return err;
}

int
FiberSectionAsym3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 6 so no conflict with matData below if just 2 fibers
  static Vector data(7);
  data(0) = this->getTag();
  data(1) = numFibers;
  data(2) = (theTorsion != 0) ? 1 : 0;
  int dbTag = this->getDbTag();
  if (theTorsion != 0) {
    data(3) = theTorsion->getClassTag();
    int torsionDbTag = theTorsion->getDbTag();
    if (torsionDbTag == 0) {
      torsionDbTag = theChannel.getDbTag();
      if (torsionDbTag != 0)
	theTorsion->setDbTag(torsionDbTag);
    }
    data(4) = torsionDbTag;    
  }
  data(5) = ys;
  data(6) = zs;

  res += theChannel.sendVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FiberSectionAsym3d::sendSelf - failed to send Vector data\n";
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
     opserr << "FiberSectionAsym3d::sendSelf - failed to send material data\n";
     return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 3*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSectionAsym3d::sendSelf - failed to send fiber data\n";
     return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);
  }

  return res;
}

int
FiberSectionAsym3d::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(7);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvVector(dbTag, commitTag, data);
  ys = data(5);
  zs = data(6);

  if (res < 0) {
   opserr << "FiberSectionAsym3d::recvSelf - failed to recv Vector data\n";
   return res;
  } 
   
  this->setTag((int)data(0));

  if ((int)data(2) == 1 && theTorsion == 0) {
    int torsionClassTag = data(3);
    int torsionDbTag = data(4);
    theTorsion = theBroker.getNewUniaxialMaterial(torsionClassTag);
    if (theTorsion == 0) {
      opserr << "FiberSection3d::recvSelf - failed to get torsion material \n";
      return -1;
    }
    theTorsion->setDbTag(torsionDbTag);
  }
  if ((int)data(2) == 1) {
    if (theTorsion->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "FiberSection3d::recvSelf - torsion failed to recvSelf \n";
      return -2;
    }
  }
  else
    theTorsion = 0;
  
  // recv data about materials objects, classTag and dbTag
  if ((int)data(1) != 0) {
    ID materialData(2*(int)data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "FiberSectionAsym3d::recvSelf - failed to send material data\n";
     return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != (int)data(1)) {
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
      numFibers = (int)data(1);
      sizeFibers = (int)data(1);
      if (numFibers != 0) {

	theMaterials = new UniaxialMaterial *[numFibers];
	
	if (theMaterials == 0) {
	  opserr << "FiberSectionAsym3d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}

	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;
	
	matData = new double [numFibers*3];

	if (matData == 0) {
	  opserr << "FiberSectionAsym3d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, 3*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSectionAsym3d::recvSelf - failed to recv fiber data\n";
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
	opserr << "FiberSectionAsym3d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    QzBar = 0.0;
    QyBar = 0.0;
    Abar  = 0.0;
    double yLoc, zLoc, Area;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = matData[3*i];
      zLoc = matData[3*i+1];
      Area = matData[3*i+2];
      Abar  += Area;
      QzBar += yLoc*Area;
      QyBar += zLoc*Area;
    }
    
    yBar = QzBar/Abar;
    zBar = QyBar/Abar;
  }    

  return res;
}

void
FiberSectionAsym3d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "\nFiberSectionAsym3d, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: (" << -yBar << ", " << zBar << ')' << endln;
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
	  s << "\"type\": \"FiberSectionAsym3d\", ";
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
FiberSectionAsym3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  
  if (argc > 2 && strcmp(argv[0],"fiber") == 0) {

    static double yLocs[10000];
    static double zLocs[10000];
    
    if (sectionIntegr != 0) {
      sectionIntegr->getFiberLocations(numFibers, yLocs, zLocs);
    }  
    else {
      for (int i = 0; i < numFibers; i++) {
	yLocs[i] = matData[3*i];
	zLocs[i] = matData[3*i+1];
      }
    }
    
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
	  //ySearch = matData[3*j];
	  //zSearch = matData[3*j+1];
	  ySearch = yLocs[j];
	  zSearch = zLocs[j];	    	  
	  dy = ySearch-yCoord;
	  dz = zSearch-zCoord;
	  closestDist = dy*dy + dz*dz;
	  key = j;
	  break;
	}
      }
      
      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
	if (matTag == theMaterials[j]->getTag()) {
	  //ySearch = matData[3*j];
	  //zSearch = matData[3*j+1];
	  ySearch = yLocs[j];
	  zSearch = zLocs[j];	    	    	  
	  dy = ySearch-yCoord;
	  dz = zSearch-zCoord;
	  distance = dy*dy + dz*dz;
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
      //ySearch = matData[0];
      //zSearch = matData[1];
      ySearch = yLocs[0];
      zSearch = zLocs[0];      
      dy = ySearch-yCoord;
      dz = zSearch-zCoord;
      closestDist = dy*dy + dz*dz;
      key = 0;
      for (int j = 1; j < numFibers; j++) {
	//ySearch = matData[3*j];
	//zSearch = matData[3*j+1];
	ySearch = yLocs[j];
	zSearch = zLocs[j];	    	    	  	
	dy = ySearch-yCoord;
	dz = zSearch-zCoord;
	distance = dy*dy + dz*dz;
	if (distance < closestDist) {
	  closestDist = distance;
	  key = j;
	}
      }
      passarg = 3;
    }
    
    if (key < numFibers && key >= 0) {
      output.tag("FiberOutput");
      output.attr("yLoc",-matData[3*key]);
      output.attr("zLoc",matData[3*key+1]);
      output.attr("area",matData[3*key+2]);
      
      theResponse = theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }
  
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
    theResponse = new MaterialResponse(this, 7, count);
  }
  else if (strcmp(argv[0],"centroid") == 0) 
    theResponse = new MaterialResponse(this, 20, Vector(2));
  
  if (theResponse == 0)
    return SectionForceDeformation::setResponse(argv, argc, output);

  return theResponse;
}


int 
FiberSectionAsym3d::getResponse(int responseID, Information &sectInfo)
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
  else if (responseID == 20) {
    static Vector centroid(2);
    centroid(0) = yBar;
    centroid(1) = zBar;
    return sectInfo.setVector(centroid);
  }
  
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int
FiberSectionAsym3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
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
FiberSectionAsym3d::getSectionDeformationSensitivity(int gradIndex)
{
  static Vector dummy(4);
  
  dummy.Zero();
  
  return dummy;
}

   
const Vector &
FiberSectionAsym3d::getStressResultantSensitivity(int gradIndex, bool conditional)
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
    ds(0) += (tmpMatrix(0,0)*e(0) + tmpMatrix(0,1)*e(1) + tmpMatrix(0,2)*e(2))*A;//Xinlong: may need to be modified because e is different now.
    ds(1) += (tmpMatrix(1,0)*e(0) + tmpMatrix(1,1)*e(1) + tmpMatrix(1,2)*e(2))*A;
    ds(2) += (tmpMatrix(2,0)*e(0) + tmpMatrix(2,1)*e(1) + tmpMatrix(2,2)*e(2))*A;
  }

  ds(3) = theTorsion->getStressSensitivity(gradIndex, conditional);

  return ds;
}

const Matrix &
FiberSectionAsym3d::getSectionTangentSensitivity(int gradIndex)
{
  static Matrix something(4,4);
  
  something.Zero();

  something(3,3) = theTorsion->getTangentSensitivity(gradIndex);
  
  return something;
}

int
FiberSectionAsym3d::commitSensitivity(const Vector& defSens, int gradIndex, int numGrads)
{

  double d0 = defSens(0);  //Xinlong: seems this should be replaced by d0~d7
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
    UniaxialMaterial *theMat = theMaterials[i];
    y = yLocs[i] - yBar;
    z = zLocs[i] - zBar;

    // determine material strain and set it
    depsdh = d0 - y*d1 + z*d2 - dydh[i]*e(1) + dzdh[i]*e(2); //Xinlong: seems this should be replaced by d0~d7

    theMat->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  theTorsion->commitSensitivity(d3, gradIndex, numGrads);

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////


