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

// $Revision: 1.8 $
// $Date: 2007-11-30 23:34:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionGJThermal.cpp,v $

// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSectionGJThermal .

#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSectionGJThermal.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <math.h>

ID FiberSectionGJThermal::code(4);

// constructors:
FiberSectionGJThermal::FiberSectionGJThermal(int tag, int num, Fiber **fibers, double gj):
  SectionForceDeformation(tag, SEC_TAG_FiberSectionGJThermal),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),dataMixed(25),
  QzBar(0.0), QyBar(0.0), ABar(0.0),
  yBar(0.0), zBar(0.0), e(4), eCommit(4), GJ(gj),sT(3), Fiber_ElongP(0), AverageThermalElong(3)
{
  if (numFibers > 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSectionGJThermal::FiberSectionGJThermal -- failed to allocate Material pointers\n";
      exit(-1);
    }

    matData = new double [numFibers*3];

    if (matData == 0) {
      opserr << "FiberSectionGJThermal::FiberSectionGJThermal -- failed to allocate double array for material data\n";
      exit(-1);
    }

    Fiber_ElongP = new double [numFibers];

    if (Fiber_ElongP == 0) {
      opserr << "FiberSectionGJThermal::FiberSectionGJThermal -- failed to allocate double array for fiber data\n";
      exit(-1);
    }
    
    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;

    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();

      QzBar += yLoc*Area;
      QyBar += zLoc*Area;
      ABar  += Area;

      matData[i*3] = -yLoc;
      matData[i*3+1] = zLoc;
      matData[i*3+2] = Area;
      UniaxialMaterial *theMat = theFiber->getMaterial();
      theMaterials[i] = theMat->getCopy();

      if (theMaterials[i] == 0) {
	opserr << "FiberSectionGJThermal::FiberSectionGJThermal -- failed to get copy of a Material\n";
	exit(-1);
      }

      Fiber_ElongP[i] = 0;      
    }

    yBar = -QzBar/ABar;
    zBar = QyBar/ABar;
  }

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

FiberSectionGJThermal::FiberSectionGJThermal(int tag, int num, double gj):
  SectionForceDeformation(tag, SEC_TAG_FiberSectionGJThermal),
  numFibers(0), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), ABar(0.0),
  yBar(0.0), zBar(0.0), e(4), eCommit(4), GJ(gj), dataMixed(25), sT(3), Fiber_ElongP(0), AverageThermalElong(3)
{
  if(sizeFibers > 0) {
    theMaterials = new UniaxialMaterial *[sizeFibers];
    
    if (theMaterials == 0) {
      opserr << "FiberSectionGJThermal::FiberSectionGJThermal -- failed to allocate Material pointers\n";
      exit(-1);
    }
    
    matData = new double [sizeFibers*3];
    
    if (matData == 0) {
      opserr << "FiberSectionGJThermal::FiberSectionGJThermal -- failed to allocate double array for material data\n";
      exit(-1);
    }

    Fiber_ElongP = new double [numFibers];
    if (Fiber_ElongP == 0) {
      opserr << "FiberSectionGJThermal::FiberSectionGJThermal -- failed to allocate double array for fiber data\n";
      exit(-1);
    }
    
    for (int i = 0; i < sizeFibers; i++) {
      matData[i*3] = 0.0;
      matData[i*3+1] = 0.0;
      matData[i*3+2] = 0.0;
      theMaterials[i] = 0;
      Fiber_ElongP[i] = 0.0;
    }
  }
  
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

// constructor for blank object that recvSelf needs to be invoked upon
FiberSectionGJThermal::FiberSectionGJThermal():
  SectionForceDeformation(0, SEC_TAG_FiberSectionGJThermal),
  numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), ABar(0.0),
  yBar(0.0), zBar(0.0), e(4), eCommit(4), GJ(1.0),dataMixed(25),sT(3), Fiber_ElongP(0), AverageThermalElong(3)
{
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
FiberSectionGJThermal::addFiber(Fiber &newFiber)
{
      // need to create a larger array
  if(numFibers == sizeFibers) {
      int newSize = 2*sizeFibers;
      UniaxialMaterial **newArray = new UniaxialMaterial *[newSize];
      double *newMatData = new double [3 * newSize];
      double *newFiberElongP = new double [newSize];
      
      if (newArray == 0 || newMatData == 0 || newFiberElongP == 0) {
	  opserr << "FiberSection3d::addFiber -- failed to allocate Fiber pointers\n";
	  exit(-1);
      }

      // copy the old pointers
      for (int i = 0; i < numFibers; i++) {
	  newArray[i] = theMaterials[i];
	  newMatData[3*i] = matData[3*i];
	  newMatData[3*i+1] = matData[3*i+1];
	  newMatData[3*i+2] = matData[3*i+2];
	  newFiberElongP[i] = Fiber_ElongP[i];
      }

      // initialize new memory
      for (int i = numFibers; i < newSize; i++) {
	  newArray[i] = 0;
	  newMatData[3*i] = 0.0;
	  newMatData[3*i+1] = 0.0;
	  newMatData[3*i+2] = 0.0;
	  newFiberElongP[i] = 0.0;
      }
      sizeFibers = newSize;

      // set new memory
      if (theMaterials != 0)
	  delete [] theMaterials;
      if (matData != 0)
	delete [] matData;
      if (Fiber_ElongP != 0)
	delete [] Fiber_ElongP;

      theMaterials = newArray;
      matData = newMatData;
      Fiber_ElongP = newFiberElongP;
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
  ABar  += Area;
  QzBar += yLoc*Area;
  QyBar += zLoc*Area;

  yBar = QzBar/ABar;
  zBar = QyBar/ABar;

  return 0;
}



// destructor:
FiberSectionGJThermal::~FiberSectionGJThermal()
{
  if (theMaterials != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];

    delete [] theMaterials;
  }

  if (matData != 0)
    delete [] matData;

  if(Fiber_ElongP != 0 )
    delete [] Fiber_ElongP;
}

int
FiberSectionGJThermal::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;

  for (int i = 0; i < 4; i++)
    sData[i] = 0.0;
  for (int i = 0; i < 16; i++)
    kData[i] = 0.0;

  int loc = 0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);
  double d3 = deforms(3);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double yi = matData[i*3]-yBar;
	double zi = matData[i*3+1]-zBar;
    double A = matData[i*3+2];
	double FiberTemperature = 0 ; 
	double FiberTempMax = 0;

	FiberTemperature = this->determineFiberTemperature( dataMixed, -yi-yBar,zi+zBar);

	//---Calculating the Fiber Temperature---end

	double strain = d0 + yi*d1 + zi*d2;  //axial strain d0, rotational degree d1,d2;
    double tangent =0.0;
	double stress = 0.0;
	double ThermalElongation = 0.0;
	 static Vector tData(4);
    static Information iData(tData);
    tData(0) = FiberTemperature;
	tData(1) = tangent;
	tData(2) = ThermalElongation;
    tData(3) = FiberTempMax;
    iData.setVector(tData);
    theMat->getVariable("ElongTangent", iData);
    tData = iData.getData();
    tangent = tData(1);
    ThermalElongation = tData(2);

    // determine material strain and set it
    strain = d0 + yi*d1 + zi*d2 - ThermalElongation;
    res += theMat->setTrial(strain, FiberTemperature, stress, tangent, ThermalElongation);

    double value = tangent * A;
    double vas1 = yi*value;
    double vas2 = zi*value;
    double vas1as2 = vas1*zi;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;

    kData[5] += vas1 * yi;
    kData[6] += vas1as2;

    kData[10] += vas2 * zi;
	//if (FiberTemperature> 450)
		//opserr << "Trial strain: " << strain << "   tangent: " << tangent << "   Tstress: " << stress << "Thelong: " << ThermalElongation <<"added s"<< sData[0]<< endln;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * yi;
    sData[2] += fs0 * zi;
  }

  kData[4] = kData[1];
  kData[8] = kData[2];
  kData[9] = kData[6];

  sData[3] = GJ*d3;
  kData[15] = GJ;
  
  return res;
}

const Matrix&
FiberSectionGJThermal::getInitialTangent(void)
{
  static double kInitialData[16];
  static Matrix kInitial(kInitialData, 4, 4);
  
  kInitial.Zero();  

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];

    double tangent = theMat->getInitialTangent();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kInitialData[0] += value;
    kInitialData[1] += vas1;
    kInitialData[2] += vas2;

    kInitialData[5] += vas1 * y;
    kInitialData[6] += vas1as2;

    kInitialData[10] += vas2 * z;
  }

  kInitialData[4] = kInitialData[1];
  kInitialData[8] = kInitialData[2];
  kInitialData[9] = kInitialData[6];

  kInitial(3,3) = GJ;

  return kInitial;
}

const Vector&
FiberSectionGJThermal::getSectionDeformation(void)
{
  return e;
}

const Matrix&
FiberSectionGJThermal::getSectionTangent(void)
{
  return *ks;
}

const Vector&
FiberSectionGJThermal::getStressResultant(void)
{
  return *s;
}
//JJadd--12.2010---to get section force due to thermal load----start-----
//---Liming modified the following block---
const Vector&
FiberSectionGJThermal::getTemperatureStress(const Vector& DataMixed)
{
   AverageThermalElong.Zero();
  dataMixed = DataMixed;

  double ThermalTangent[1000];
  double DeltaThermalElong[1000];

  for (int i = 0; i < numFibers; i++) {
       ThermalTangent[i]=0;
      DeltaThermalElong[i]=0;
  }

  for (int i = 0; i < numFibers; i++) {

    UniaxialMaterial *theMat = theMaterials[i];

	//double seefiberlocs1,seefiberlocs2;
    //seefiberlocs1 = fiberLocsZ[i];
	int jy;
	int jz;
    jy = i*3; //retrieve temp along y
	jz = i*3+1; //retrieve temp along z
	double yi = matData[jy];
	double zi = matData[jz];
	double FiberTemperature = 0 ; 
	double FiberTempMax = 0;

	FiberTemperature= this->determineFiberTemperature( dataMixed, -yi, zi);
    // determine material strain and set it
	double tangent =0.0;
	double ThermalElongation =0.0;
    static Vector tData(4);
    static Information iData(tData);
    tData(0) = FiberTemperature;
	tData(1) = tangent;
	tData(2) = ThermalElongation;
    tData(3) = FiberTempMax;
    iData.setVector(tData);
    theMat->getVariable("ElongTangent", iData);
    tData = iData.getData();
	FiberTemperature = tData(0);
    tangent = tData(1);
    ThermalElongation = tData(2);
	FiberTempMax = tData(3);

	DeltaThermalElong[i]= ThermalElongation-Fiber_ElongP[i];
	Fiber_ElongP[i]= ThermalElongation;

    //  double strain = -ThermalElongation;
    //  theMat->setTrialTemperature(strain, FiberTemperature, stress, tangent, ThermalElongation);
    ThermalTangent[i] = tangent;
  }

 // calculate section resisting force due to thermal load
 double FiberForce;
  double SectionArea=0;
  double ThermalForce=0;
  double ThermalMoment1 =0;
  double SectionMomofArea1 =0;
  double ThermalMoment2 =0;
  double SectionMomofArea2 =0;
 //Liming add this for calculating average Thermal Elongation
  sT.Zero();
  for (int i = 0; i < numFibers; i++) {
	  FiberForce = ThermalTangent[i]*matData[3*i+2]*DeltaThermalElong[i];
	  SectionArea +=matData[3*i+2];

	  SectionMomofArea1 += (matData[3*i+2]*(matData[3*i] - yBar));
	  SectionMomofArea2 += (matData[3*i+2]*(matData[3*i+1] - zBar));
	  ThermalForce += Fiber_ElongP[i]*matData[3*i+2];
	  ThermalMoment1 += Fiber_ElongP[i]*matData[3*i+2]*(matData[3*i] - yBar);
	  ThermalMoment2 += Fiber_ElongP[i]*matData[3*i+2]*(matData[3*i+1] - zBar);

	  sT(0) += FiberForce;
	  sT(1) += FiberForce*(matData[3*i] - yBar);
	  sT(2) += FiberForce*(matData[3*i+1] - zBar);
  }
  AverageThermalElong(0) = ThermalForce/SectionArea;
 AverageThermalElong(1) = ThermalMoment1/SectionMomofArea1;
  AverageThermalElong(2) = ThermalMoment2/SectionMomofArea2;
  //double ThermalMoment;
  //ThermalMoment = abs(sTData[1]);
  //sTData[1] = ThermalMoment;
  sT(1) = abs(sT(1));

  return sT;
}

const Vector&
FiberSectionGJThermal::getThermalElong(void)
{
  return AverageThermalElong;
}

//------Liming-Modified for beamThermalAction3d-----
//JJadd--12.2010---to get section force due to thermal load----end-----
SectionForceDeformation*
FiberSectionGJThermal::getCopy(void)
{
  FiberSectionGJThermal *theCopy = new FiberSectionGJThermal();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;

  if (numFibers > 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "FiberSectionGJThermal::getCopy -- failed to allocate Material pointers\n";
      exit(-1);
    }

    theCopy->matData = new double [numFibers*3];

    if (theCopy->matData == 0) {
      opserr << "FiberSectionGJThermal::getCopy -- failed to allocate double array for material data\n";
      exit(-1);
    }
    theCopy->Fiber_ElongP = new double [numFibers];
    if (theCopy->Fiber_ElongP == 0) {
      opserr << "FiberSectionGJThermal::getCopy -- failed to allocate double array for fiber data\n";
      exit(-1);
    }
    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*3] = matData[i*3];
      theCopy->matData[i*3+1] = matData[i*3+1];
      theCopy->matData[i*3+2] = matData[i*3+2];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr << "FiberSectionGJThermal::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }
      theCopy->Fiber_ElongP[i] = Fiber_ElongP[i];
    }
  }

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;

  for (int i=0; i<16; i++)
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];
  theCopy->sData[3] = sData[3];

  theCopy->GJ = GJ;
  theCopy->sT = sT;
  
  return theCopy;
}

const ID&
FiberSectionGJThermal::getType ()
{
  return code;
}

int
FiberSectionGJThermal::getOrder () const
{
  return 4;
}

int
FiberSectionGJThermal::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  eCommit = e;

  return err;
}

int
FiberSectionGJThermal::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0; 
  kData[15] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; sData[3] = 0.0;

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;

    kData[5] += vas1 * y;
    kData[6] += vas1as2;

    kData[10] += vas2 * z;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  kData[4] = kData[1];
  kData[8] = kData[2];
  kData[9] = kData[6];

  kData[15] = GJ;
  sData[3] = GJ*e(3);
  
  return err;
}

int
FiberSectionGJThermal::revertToStart(void)
{
  // revert the fibers to start
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0;
  kData[15] = 0.0; 
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; sData[3] = 0.0;
  
  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];

    // invoke revertToStart on the material
    err += theMat->revertToStart();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;

    kData[5] += vas1 * y;
    kData[6] += vas1as2;

    kData[10] += vas2 * z;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  kData[4] = kData[1];
  kData[8] = kData[2];
  kData[9] = kData[6];

  kData[15] = GJ;
  sData[3] = GJ*e(3);
  
  return err;
}

int
FiberSectionGJThermal::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers,
  //     size 4 so no conflict with fiberData below if just 1 fiber
  static Vector data(4);
  data(0) = this->getTag();
  data(1) = numFibers;
  data(2) = GJ;
  int dbTag = this->getDbTag();
  res += theChannel.sendVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FiberSection2d::sendSelf - failed to send ID data\n";
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
      opserr << "FiberSection2d::sendSelf- failed to send material data\n";
      return res;
    }

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 3*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr << "FiberSection2d::sendSelf - failed to send material data\n";
      return res;
    }

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);
  }

  return res;
}

int
FiberSectionGJThermal::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(4);

  int dbTag = this->getDbTag();
  res += theChannel.recvVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FiberSection2d::recvSelf - failed to recv ID data\n";
    return res;
  }
  this->setTag((int)data(0));
  GJ = data(2);

  // recv data about materials objects, classTag and dbTag
  numFibers = (int)data(1);
  if (numFibers != 0) {
    ID materialData(2*numFibers);
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr << "FiberSection2d::recvSelf - failed to send material data\n";
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
      if (numFibers != 0) {

	theMaterials = new UniaxialMaterial *[numFibers];

	if (theMaterials == 0) {
	  opserr << "FiberSection2d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}


	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;

	matData = new double [numFibers*3];

	if (matData == 0) {
	  opserr << "FiberSection2d::recvSelf -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }

    Vector fiberData(matData, 3*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr << "FiberSection2d::recvSelf - failed to send material data\n";

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
	opserr << "FiberSection2d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    double Qz = 0.0;
    double Qy = 0.0;
    double A  = 0.0;
    double yLoc, zLoc, Area;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = -matData[2*i];
      zLoc = matData[2*i+1];
      Area = matData[2*i+2];
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
FiberSectionGJThermal::Print(OPS_Stream &s, int flag)
{
  s << "\nFiberSectionGJThermal, tag: " << this->getTag() << endln;
  s << "\tSection code: " << code;
  s << "\tNumber of Fibers: " << numFibers << endln;
  s << "\tCentroid: (" << -yBar << ", " << zBar << ')' << endln;
  s << "\tTorsional Stiffness: " << GJ << endln;

  if (flag == 1) {
    int loc = 0;
    for (int i = 0; i < numFibers; i++) {
      s << "\nLocation (y, z) = (" << -matData[loc++] << ", " << matData[loc++] << ")";
      s << "\nArea = " << matData[loc++] << endln;
      theMaterials[i]->Print(s, flag);
    }
  }
}

Response*
FiberSectionGJThermal::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  // See if the response is one of the defaults
  Response *theResponse = SectionForceDeformation::setResponse(argv, argc, output);
  if (theResponse != 0)
    return theResponse;


  if (argc <=2 || strcmp(argv[0],"fiber") != 0)
    return 0;

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
	ySearch = -matData[3*j];
	zSearch =  matData[3*j+1];
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
	ySearch = -matData[3*j];
	zSearch =  matData[3*j+1];
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
      ySearch = -matData[3*j];
      zSearch =  matData[3*j+1];
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
    output.attr("yLoc",-matData[2*key]);
    output.attr("zLoc",matData[2*key+1]);
    output.attr("area",matData[2*key+2]);

    theResponse =  theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);

    output.endTag();
  }

  return theResponse;
}


int
FiberSectionGJThermal::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int
FiberSectionGJThermal::setParameter (const char **argv, int argc, Parameter &param)
{
  if (argc < 3)
    return 0;

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
  }

  return result;
}

double  
FiberSectionGJThermal::determineFiberTemperature(const Vector& DataMixed, double fiberLocy, double fiberLocz) 
{
	double FiberTemperature = 0;
	if(DataMixed.Size()==18){
	//--------------if temperature Data has 18 elements--------------------
		if ( fabs(DataMixed(1)) <= 1e-10 && fabs(DataMixed(17)) <= 1e-10 ) //no tempe load
		{
			return 0 ;
		}
		
		double dataTempe[18]; //PK changed 18 to 27 to pass max temps
		for (int i = 0; i < 18; i++) { 
			dataTempe[i] = DataMixed(i);
		}

		if (  fiberLocy <= dataTempe[1])
		{
			opserr <<"FiberSectionGJThermal "<<this->getTag()<<":: fiber locy "<< fiberLocy <<" is out of the section below "<< dataTempe[1]<<endln;
		}
		else if (fiberLocy <= dataTempe[3])
		{
			FiberTemperature = dataTempe[0] - (dataTempe[1] - fiberLocy) * (dataTempe[0] - dataTempe[2])/(dataTempe[1] - dataTempe[3]);
		}
		else if (   fiberLocy <= dataTempe[5] )
		{
			FiberTemperature = dataTempe[2] - (dataTempe[3] - fiberLocy) * (dataTempe[2] - dataTempe[4])/(dataTempe[3] - dataTempe[5]);
		}
		else if ( fiberLocy <= dataTempe[7] )
		{
			FiberTemperature = dataTempe[4] - (dataTempe[5] - fiberLocy) * (dataTempe[4] - dataTempe[6])/(dataTempe[5] - dataTempe[7]);
		}
		else if ( fiberLocy <= dataTempe[9] )
		{
			FiberTemperature = dataTempe[6] - (dataTempe[7] - fiberLocy) * (dataTempe[6] - dataTempe[8])/(dataTempe[7] - dataTempe[9]);
		}
		else if (fiberLocy <= dataTempe[11] )
		{
			FiberTemperature = dataTempe[8] - (dataTempe[9] - fiberLocy) * (dataTempe[8] - dataTempe[10])/(dataTempe[9] - dataTempe[11]);
		}
		else if (fiberLocy <= dataTempe[13] )
		{
			FiberTemperature = dataTempe[10] - (dataTempe[11] - fiberLocy) * (dataTempe[10] - dataTempe[12])/(dataTempe[11] - dataTempe[13]);
		}
		else if (fiberLocy <= dataTempe[15] )
		{
			FiberTemperature = dataTempe[12] - (dataTempe[13] - fiberLocy) * (dataTempe[12] - dataTempe[14])/(dataTempe[13] - dataTempe[15]);
		}
		else if ( fiberLocy <= dataTempe[17] )
		{
			FiberTemperature = dataTempe[14] - (dataTempe[15] - fiberLocy) * (dataTempe[14] - dataTempe[16])/(dataTempe[15] - dataTempe[17]);
		}
		else 
		{
			opserr <<"FiberSectionGJThermal " << this->getTag() << " :: fiber loc " <<fiberLocy<<" is out of the section over" << dataTempe[17] << endln;
		}
	}
	else if(DataMixed.Size()==25){
	//---------------if temperature Data has 25 elements--------------------
	
		double dataTempe[25]; //
		for (int i = 0; i < 25; i++) { //
			dataTempe[i] = DataMixed(i);
		}

		if ( fabs(dataTempe[0]) <= 1e-10 && fabs(dataTempe[10]) <= 1e-10 &&fabs(dataTempe[11]) <= 1e-10) //no tempe load
		{
			return 0;
		}
	
	//calculate the fiber tempe, T=T1-(Y-Y1)*(T1-T2)/(Y1-Y2)
	//first for bottom flange if existing
		if (  fiberLocy <= dataTempe[1]) 
		{
			if (fiberLocz <= dataTempe[12]){
			opserr<<"WARNING: FiberSectionGJThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<endln;
			}
			else if (fiberLocz<= dataTempe[15]){
			FiberTemperature = dataTempe[10] - (dataTempe[10] - dataTempe[13])*(dataTempe[12] - fiberLocz) /(dataTempe[12] - dataTempe[15]);
			}
			else if (fiberLocz<= dataTempe[18]){
			FiberTemperature = dataTempe[13] - (dataTempe[13] - dataTempe[16])*(dataTempe[15] - fiberLocz) /(dataTempe[15] - dataTempe[18]);
			}
			else if (fiberLocz<= dataTempe[21]){
			FiberTemperature = dataTempe[16] - (dataTempe[16] - dataTempe[19])*(dataTempe[18] - fiberLocz) /(dataTempe[18] - dataTempe[21]);
			}
			else if (fiberLocz<= dataTempe[24]){
			FiberTemperature = dataTempe[19] - (dataTempe[19] - dataTempe[22])*(dataTempe[21] - fiberLocz) /(dataTempe[21] - dataTempe[24]);
			}
			else {
			opserr<<"WARNING: FiberSectionGJThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<endln;
			}
		}
		else if (fiberLocy <= dataTempe[3])
		{
			FiberTemperature = dataTempe[0] - (dataTempe[1] - fiberLocy) * (dataTempe[0] - dataTempe[2])/(dataTempe[1] - dataTempe[3]);
		}
		else if (   fiberLocy <= dataTempe[5] )
		{
			FiberTemperature = dataTempe[2] - (dataTempe[3] - fiberLocy) * (dataTempe[2] - dataTempe[4])/(dataTempe[3] - dataTempe[5]);
		}
		else if ( fiberLocy <= dataTempe[7] )
		{
			FiberTemperature = dataTempe[4] - (dataTempe[5] - fiberLocy) * (dataTempe[4] - dataTempe[6])/(dataTempe[5] - dataTempe[7]);
		}
		else if ( fiberLocy <= dataTempe[9] )
		{
			FiberTemperature = dataTempe[6] - (dataTempe[7] - fiberLocy) * (dataTempe[6] - dataTempe[8])/(dataTempe[7] - dataTempe[9]);
		}
		else {
			if (fiberLocz <= dataTempe[12]){
			opserr<<"WARNING: FiberSectionGJThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<endln;
			}
			else if (fiberLocz<= dataTempe[15]){
			FiberTemperature = dataTempe[11] - (dataTempe[11] - dataTempe[14])*(dataTempe[12] - fiberLocz) /(dataTempe[12] - dataTempe[15]);
			}
			else if (fiberLocz<= dataTempe[18]){
			FiberTemperature = dataTempe[14] - (dataTempe[14] - dataTempe[17])*(dataTempe[15] - fiberLocz) /(dataTempe[15] - dataTempe[18]);
			}
			else if (fiberLocz<= dataTempe[21]){
			FiberTemperature = dataTempe[17] - (dataTempe[17] - dataTempe[20])*(dataTempe[18] - fiberLocz) /(dataTempe[18] - dataTempe[21]);
			}
			else if (fiberLocz<= dataTempe[24]){
			FiberTemperature = dataTempe[20] - (dataTempe[20] - dataTempe[23])*(dataTempe[21] - fiberLocz) /(dataTempe[21] - dataTempe[24]);
			}
			else {
			opserr<<"WARNING: FiberSectionGJThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<endln;
			}
		}
	}
	return FiberTemperature;
}
