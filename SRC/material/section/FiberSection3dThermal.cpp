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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSection3dThermal.cpp,v $

// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSection2d.
// Modified for SIF modelling by Jian Jiang,Liming Jiang [http://openseesforfire.github.io]


#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSection3dThermal.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include <math.h>
#include <elementAPI.h>

ID FiberSection3dThermal::code(3);

void* OPS_FiberSection3dThermal()
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
    SectionForceDeformation *section = new FiberSection3dThermal(tag, num, computeCentroid);
    if (deleteTorsion)
      delete torsion;
    return section;
}

// constructors:
FiberSection3dThermal::FiberSection3dThermal(int tag, int num, Fiber **fibers, bool compCentroid):
  SectionForceDeformation(tag, SEC_TAG_FiberSection3dThermal),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), ABar(0.0), yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
  e(3), eCommit(3), s(0), ks(0), sT(3), Fiber_T(0), Fiber_TMax(0),
  parameterID(0), SHVs(0)
{
  if (numFibers > 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to allocate Material pointers\n";
      exit(-1);
    }

    matData = new double [numFibers*3];

    if (matData == 0) {
      opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to allocate double array for material data\n";
      exit(-1);
    }

    Fiber_T = new double [numFibers];
    if (Fiber_T == 0) {
      opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to allocate double array for fiber data\n";
      exit(-1);
    }

    Fiber_TMax = new double [numFibers];
    if (Fiber_TMax == 0) {
      opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to allocate double array for fiber data\n";
      exit(-1);
    }    
    
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
	opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to get copy of a Material\n";
	exit(-1);
      }

      Fiber_T[i] = 0.0;
      Fiber_TMax[i] = 0.0;
    }

    if (computeCentroid) {
      yBar = QzBar/ABar;
      zBar = QyBar/ABar;
    }
  }

  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  for (int i=0; i<9; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;

 // AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
  SHVs=0;
  // AddingSensitivity:END //////////////////////////////////////
}

FiberSection3dThermal::FiberSection3dThermal(int tag, int num, bool compCentroid):
  SectionForceDeformation(tag, SEC_TAG_FiberSection3dThermal),
  numFibers(0), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), ABar(0.0), yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
  e(3), eCommit(3), s(0), ks(0),
  sT(3), Fiber_T(0), Fiber_TMax(0),
  parameterID(0), SHVs(0)
{
  if(sizeFibers > 0) {
    theMaterials = new UniaxialMaterial *[sizeFibers];
    
    if (theMaterials == 0) {
      opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to allocate Material pointers\n";
      exit(-1);
    }
    
    matData = new double [sizeFibers*3];
    
    if (matData == 0) {
      opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to allocate double array for material data\n";
      exit(-1);
    }

    Fiber_T = new double [numFibers];
    if (Fiber_T == 0) {
      opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to allocate double array for fiber data\n";
      exit(-1);
    }

    Fiber_TMax = new double [numFibers];
    if (Fiber_TMax == 0) {
      opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to allocate double array for fiber data\n";
      exit(-1);
    }
    
    for (int i = 0; i < sizeFibers; i++) {
      matData[i*3] = 0.0;
      matData[i*3+1] = 0.0;
      matData[i*3+2] = 0.0;
      theMaterials[i] = 0;
      Fiber_T[i] = 0.0;
      Fiber_TMax[i] = 0.0;
    }
  }
  
  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  for (int i=0; i<9; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;

 // AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
  SHVs=0;
  // AddingSensitivity:END //////////////////////////////////////
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSection3dThermal::FiberSection3dThermal():
  SectionForceDeformation(0, SEC_TAG_FiberSection3dThermal),
  numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
  QzBar(0.0), QyBar(0.0), ABar(0.0), yBar(0.0), zBar(0.0), computeCentroid(true),
  e(3), eCommit(3), s(0), ks(0),
  sT(3), Fiber_T(0), Fiber_TMax(0),
  parameterID(0), SHVs(0)
{
  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  for (int i=0; i<9; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;

 // AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
  SHVs=0;
  // AddingSensitivity:END //////////////////////////////////////
}

int
FiberSection3dThermal::addFiber(Fiber &newFiber)
{
    // need to create a larger array
  if(numFibers == sizeFibers) {
      int newSize = 2*sizeFibers;
      UniaxialMaterial **newArray = new UniaxialMaterial *[newSize];
      double *newMatData = new double [3 * newSize];
      double *newFiberT = new double [newSize];
      double *newFiberTMax = new double [newSize];
      if (newArray == 0 || newMatData == 0 || newFiberT == 0 || newFiberTMax == 0) {
	  opserr << "FiberSection3dThermal::addFiber -- failed to allocate Fiber pointers\n";
	  exit(-1);
      }

      // copy the old pointers
      for (int i = 0; i < numFibers; i++) {
	  newArray[i] = theMaterials[i];
	  newMatData[3*i] = matData[3*i];
	  newMatData[3*i+1] = matData[3*i+1];
	  newMatData[3*i+2] = matData[3*i+2];
	  newFiberT[i] = Fiber_T[i];
	  newFiberTMax[i] = Fiber_TMax[i];
      }

      // initialize new memory
      for (int i = numFibers; i < newSize; i++) {
	  newArray[i] = 0;
	  newMatData[3*i] = 0.0;
	  newMatData[3*i+1] = 0.0;
	  newMatData[3*i+2] = 0.0;
	  newFiberT[i] = 0.0;
	  newFiberTMax[i] = 0.0;
      }
      sizeFibers = newSize;

      // set new memory
      if (theMaterials != 0)
	  delete [] theMaterials;
      if (matData != 0)
	  delete [] matData;
      if (Fiber_T != 0)
	delete [] Fiber_T;
      if (Fiber_TMax != 0)
	delete [] Fiber_TMax;      

      theMaterials = newArray;
      matData = newMatData;
      Fiber_T = newFiberT;
      Fiber_TMax = newFiberTMax;
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
    opserr << "FiberSection3dThermal::addFiber -- failed to get copy of a Material\n";
    return -1;
  }

  numFibers++;

  // Recompute centroid
  if (computeCentroid) {
    ABar  += Area;
    QzBar += yLoc*Area;
    QyBar += zLoc*Area;
    
    yBar = QzBar/ABar;
    zBar = QyBar/ABar;
  }
  
  return 0;
}



// destructor:
FiberSection3dThermal::~FiberSection3dThermal()
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

  //if (TemperatureTangent != 0)
    //delete [] TemperatureTangent;

  if (Fiber_T != 0)
    delete [] Fiber_T;
  if (Fiber_TMax != 0)
    delete [] Fiber_TMax;  
}

int
FiberSection3dThermal::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0;

  int loc = 0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);


  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];

	double FiberTemperature = Fiber_T[i]; //Added by Liming to obtain fiber T;
    double FiberTempMax= Fiber_TMax[i]; //Maximum Temp;


    int jy;
	int jz;
    jy = i*3; //retrieve temp along y
	jz = i*3+1; //retrieve temp along z
	double yi;
	double zi;
	yi = matData[jy];
    zi = matData[jz];


	//---Calculating the Fiber Temperature---end

	double strain = d0 + y*d1 + z*d2;  //axial strain d0, rotational degree d1,d2;
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
    strain = d0 + y*d1 + z*d2 - ThermalElongation;
    res += theMat->setTrial(strain, FiberTemperature, stress, tangent, ThermalElongation);

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;

    kData[4] += vas1 * y;
    kData[5] += vas1as2;

    kData[8] += vas2 * z;

    double fs0 = stress * A;

    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  kData[3] = kData[1];
  kData[6] = kData[2];
  kData[7] = kData[5];

  return res;
}

const Matrix&
FiberSection3dThermal::getInitialTangent(void)
{
  static double kInitialData[9];
  static Matrix kInitial(kInitialData, 3, 3);

  kInitialData[0] = 0.0; kInitialData[1] = 0.0;
  kInitialData[2] = 0.0; kInitialData[3] = 0.0;
  kInitialData[4] = 0.0; kInitialData[5] = 0.0;
  kInitialData[6] = 0.0; kInitialData[7] = 0.0;
  kInitialData[8] = 0.0;

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

    kInitialData[4] += vas1 * y;
    kInitialData[5] += vas1as2;

    kInitialData[8] += vas2 * z;
  }

  kInitialData[3] = kInitialData[1];
  kInitialData[6] = kInitialData[2];
  kInitialData[7] = kInitialData[5];

  return kInitial;
}

const Vector&
FiberSection3dThermal::getSectionDeformation(void)
{
  return e;
}

const Matrix&
FiberSection3dThermal::getSectionTangent(void)
{
  return *ks;
}

const Vector&
FiberSection3dThermal::getStressResultant(void)
{
  return *s;
}


//JJadd--12.2010---to get section force due to thermal load----start-----
const Vector&
FiberSection3dThermal::getTemperatureStress(const Vector& dataMixed)
{
  sT.Zero();

  //JJadd, 12/2010, updata yBar = Ai*Ei*yi/(Ai*E*)  start
  double ThermalTangent[1000];
  double ThermalElong[1000];
  for (int i = 0; i < numFibers; i++) {
          ThermalTangent[i]=0;
          ThermalElong[i]=0;
  }

  for (int i = 0; i < numFibers; i++) {

    UniaxialMaterial *theMat = theMaterials[i];

	//double seefiberlocs1,seefiberlocs2;
    //seefiberlocs1 = fiberLocsZ[i];
	int jy;
	int jz;
    jy = i*3; //retrieve temp along y
	jz = i*3+1; //retrieve temp along z
	double yi;
	double zi;
	yi = matData[jy];
    zi = matData[jz];

	double FiberTemperature = 0 ; //JZ
	double FiberTempMax=0; //PK add for max temp

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

    //  double strain = -ThermalElongation;
    //  theMat->setTrialTemperature(strain, FiberTemperature, stress, tangent, ThermalElongation);
    Fiber_T[i] = FiberTemperature;
	Fiber_TMax[i] = FiberTempMax;
    ThermalTangent[i] = tangent;
	ThermalElong[i] = ThermalElongation;

  }

 // calculate section resisting force due to thermal load

  double FiberForce;
  for (int i = 0; i < numFibers; i++) {
	  FiberForce = ThermalTangent[i]*matData[3*i+2]*ThermalElong[i];
	  sT(0) += FiberForce;
	  sT(1) += FiberForce*(matData[3*i] - yBar);
	  sT(2) += FiberForce*(matData[3*i+1] - zBar);
  }
  //double ThermalMoment;
  //ThermalMoment = abs(sTData[1]);
 // sTData[1] = ThermalMoment;

  return sT;
}
//JJadd--12.2010---to get section force due to thermal load----end-----



SectionForceDeformation*
FiberSection3dThermal::getCopy(void)
{
  FiberSection3dThermal *theCopy = new FiberSection3dThermal ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;

  if (numFibers > 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to allocate Material pointers\n";
      exit(-1);
    }

    theCopy->matData = new double [numFibers*3];

    if (theCopy->matData == 0) {
      opserr << "FiberSection3dThermal::getCopy -- failed to allocate double array for material data\n";
      exit(-1);
    }

    theCopy->Fiber_T = new double [numFibers];
    theCopy->Fiber_TMax = new double [numFibers];
    if (theCopy->Fiber_TMax == 0 || theCopy->Fiber_TMax) {
      opserr << "FiberSection3dThermal::getCopy -- failed to allocate double array for fiber data\n";
      exit(-1);
    }        

    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*3] = matData[i*3];
      theCopy->matData[i*3+1] = matData[i*3+1];
      theCopy->matData[i*3+2] = matData[i*3+2];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr << "FiberSection3dThermal::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }

      theCopy->Fiber_T[i] = Fiber_T[i];
      theCopy->Fiber_TMax[i] = Fiber_TMax[i];
    }
  }

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->QzBar = QzBar;
  theCopy->QyBar = QyBar;
  theCopy->ABar = ABar;  
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;
  theCopy->computeCentroid = computeCentroid;
  
  for (int i=0; i<9; i++)
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];
  theCopy->sT = sT;
  
  return theCopy;
}

const ID&
FiberSection3dThermal::getType ()
{
  return code;
}

int
FiberSection3dThermal::getOrder () const
{
  return 3;
}

int
FiberSection3dThermal::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  eCommit = e;

  return err;
}

int
FiberSection3dThermal::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;


  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0;

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

    kData[4] += vas1 * y;
    kData[5] += vas1as2;

    kData[8] += vas2 * z;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  kData[3] = kData[1];
  kData[6] = kData[2];
  kData[7] = kData[5];

  return err;
}

int
FiberSection3dThermal::revertToStart(void)
{
  // revert the fibers to start
  int err = 0;


  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0;

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

    kData[4] += vas1 * y;
    kData[5] += vas1as2;

    kData[8] += vas2 * z;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  kData[3] = kData[1];
  kData[6] = kData[2];
  kData[7] = kData[5];

  return err;
}

int
FiberSection3dThermal::sendSelf(int commitTag, Channel &theChannel)
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
     opserr << "FiberSection2d::sendSelf - failed to send material data\n";
     return res;
    }

    // send the fiber data, i.e. area and loc, T, and Tmax
    Vector fiberData(5*numFibers);
    for (int i = 0; i < numFibers; i++) {
      fiberData(            i) = matData[i];
      fiberData(  numFibers+i) = matData[numFibers+i];
      fiberData(2*numFibers+i) = matData[2*numFibers+i];      
      fiberData(3*numFibers+i) = Fiber_T[i];
      fiberData(4*numFibers+i) = Fiber_TMax[i];
    }    
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
FiberSection3dThermal::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);

  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);

  if (res < 0) {
   opserr << "FiberSection2d::sendSelf - failed to recv ID data\n";
   return res;
  }

  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "FiberSection2d::sendSelf - failed to send material data\n";
     return res;
    }

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != data(1)) {
      // delete old stuff if outa date
      if (theMaterials != 0) {
	for (int i=0; i<numFibers; i++)
	  delete theMaterials[i];
	delete [] theMaterials;
	theMaterials = 0;	
      }
      if (matData != 0)
	delete [] matData;
      if (Fiber_T != 0)
	delete [] Fiber_T;
      if (Fiber_TMax != 0)
	delete [] Fiber_TMax;	
      Fiber_T = 0;
      Fiber_TMax = 0;	
      matData = 0;

      // create memory to hold material pointers and fiber data
      numFibers = data(1);
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
	  opserr << "FiberSection2d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}

	Fiber_T = new double [numFibers];
	if (Fiber_T == 0) {
	  opserr <<"FiberSection3dThermal::recvSelf  -- failed to allocate double array for fiber T\n";
	  exit(-1);
	}
	Fiber_TMax = new double [numFibers];
	if (Fiber_TMax == 0) {
	  opserr <<"FiberSection3dThermal::recvSelf  -- failed to allocate double array for fiber TMax\n";
	  exit(-1);
	}			
      }
    }

    Vector fiberData(5*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSection2d::sendSelf - failed to send material data\n";
     return res;
    }
    for (int i = 0; i < numFibers; i++) {
      matData[i] =             fiberData(          i);
      matData[numFibers+i] =   fiberData(numFibers+i);
      matData[2*numFibers+i] = fiberData(2*numFibers+i);
      Fiber_T[i]    = fiberData(3*numFibers+i);
      Fiber_TMax[i] = fiberData(4*numFibers+i);
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

    QzBar = 0.0;
    QyBar = 0.0;
    ABar = 0.0;

    computeCentroid = data(2) ? true : false;
    
    // Recompute centroid
    double yLoc, zLoc, Area;    
    for (i = 0; computeCentroid && i < numFibers; i++) {
      yLoc = matData[3*i];
      zLoc = matData[3*i+1];
      Area = matData[3*i+2];
      ABar  += Area;
      QzBar += yLoc*Area;
      QyBar += zLoc*Area;
    }

    if (computeCentroid) {
      yBar = QzBar/ABar;
      zBar = QyBar/ABar;
    } else {
      yBar = 0.0;
      zBar = 0.0;
    }
  }

  return res;
}

void
FiberSection3dThermal::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {
    for (int i = 0; i < numFibers; i++) {
      s << -matData[3*i] << " "  << matData[3*i+1] << " "  << matData[3*i+2] << " " ;
      s << theMaterials[i]->getStress() << " "  << theMaterials[i]->getStrain() << endln;
    }
  } else {
    s << "\nFiberSection3dThermal, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: (" << yBar << ", " << zBar << ')' << endln;

    if (flag == 1) {
      for (int i = 0; i < numFibers; i++) {
	s << "\nLocation (y, z) = (" << -matData[3*i] << ", " << matData[3*i+1] << ")";
	s << "\nArea = " << matData[3*i+2] << endln;
      theMaterials[i]->Print(s, flag);
      }
    }
  }
}

Response*
FiberSection3dThermal::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  
  if (argc > 2 && strcmp(argv[0],"fiber") == 0) {

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
      output.attr("yLoc",matData[3*key]);
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
  }

  if (theResponse == 0)
    return SectionForceDeformation::setResponse(argv, argc, output);

  return theResponse;
}


int
FiberSection3dThermal::getResponse(int responseID, Information &sectInfo)
{
  if (responseID == 5) {
    int numData = 5*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      double yLoc, zLoc, A, stress, strain;
      yLoc = -matData[3*j];
      zLoc =  matData[3*j+1];
      A = matData[3*j+2];
      stress = theMaterials[j]->getStress();
      strain = theMaterials[j]->getStrain();
      data(count) = yLoc; data(count+1) = zLoc; data(count+2) = A;
      data(count+3) = stress; data(count+4) = strain;
      count += 5;
    }
    return sectInfo.setVector(data);
  } else
    return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int
FiberSection3dThermal::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 3)
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
FiberSection3dThermal::getSectionDeformationSensitivity(int gradIndex)
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
FiberSection3dThermal::getStressResultantSensitivity(int gradIndex, bool conditional)
{

  static Vector ds(3);

  ds.Zero();

  double  stressGradient;
  int loc = 0;


  for (int i = 0; i < numFibers; i++) {
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];
    stressGradient = theMaterials[i]->getStressSensitivity(gradIndex,conditional);
    stressGradient *=  A;
    ds(0) += stressGradient;
    ds(1) += stressGradient * y;
    ds(2) += stressGradient * z;

  }  //for

  return ds;
}

const Matrix &
FiberSection3dThermal::getSectionTangentSensitivity(int gradIndex)
{
  static Matrix something(3,3);

  something.Zero();

  return something;
}

int
FiberSection3dThermal::commitSensitivity(const Vector& defSens, int gradIndex, int numGrads)
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
   	double y = matData[loc++] - yBar;
	double z = matData[loc++] - zBar;
	loc++;   // skip A data.

	double strainSens = d0 + y*d1 + z*d2;

	theMaterials[i]->commitSensitivity(strainSens,gradIndex,numGrads);
  }

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////


double
FiberSection3dThermal::determineFiberTemperature(const Vector& DataMixed, double fiberLocy, double fiberLocz)
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
			opserr <<"FiberSection2dThermal::setTrialSectionDeformationTemperature -- fiber loc is out of the section";
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
			opserr <<"FiberSection3dThermal::setTrialSectionDeformation -- fiber loc " <<fiberLocy<<" is out of the section"<<endln;
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
			opserr<<"WARNING: FiberSection3dThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<endln;
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
			opserr<<"WARNING: FiberSection3dThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<endln;
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
			opserr<<"WARNING: FiberSection3dThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<endln;
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
			opserr<<"WARNING: FiberSection3dThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<endln;
			}
		}
	}
	return FiberTemperature;
}
