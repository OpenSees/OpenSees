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

//Modified by Jian Zhang, [University of Edinburgh]
//Modified by Panagiotis Kotsovinos, [University of Edinburgh]
//Modified by Liming Jiang, [University of Edinburgh]

// Description: This file contains the class implementation of FiberSection2dThermal.

#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSection2dThermal.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <SectionIntegration.h>
#include <math.h> //JZ
#include <elementAPI.h>

//#include "ThermalField.h"
//#include "ThermalField2d.h"

ID FiberSection2dThermal::code(2);

void* OPS_FiberSection2dThermal()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	opserr<<"insufficient arguments for FiberSection2dThermal\n";
	return 0;
    }

    numData = 1;
    int tag;
    if(OPS_GetIntInput(&numData,&tag) < 0) return 0;

    int num = 30;
    return new FiberSection2dThermal(tag,num);
}

// constructors:
FiberSection2dThermal::FiberSection2dThermal(int tag, int num, Fiber **fibers, bool compCentroid):
  SectionForceDeformation(tag, SEC_TAG_FiberSection2dThermal),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(compCentroid),
  sectionIntegr(0), e(2), eCommit(2), s(0), ks(0),
  DataMixed(27), sT(2), Fiber_Tangent(0), Fiber_ElongP(0), AverageThermalElong(2),
  dedh(2)//,theTemperatures(temperatures),theTemperatureFactor(0)
{
  if (numFibers > 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to allocate Material pointers";
      exit(-1);
    }

    matData = new double [numFibers*2];

    if (matData == 0) {
      opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to allocate double array for material data\n";
      exit(-1);
    }

    Fiber_Tangent = new double [numFibers];
    Fiber_ElongP = new double [numFibers];
    if (Fiber_Tangent == 0 || Fiber_ElongP == 0) {
      opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to allocate double array for fiber data data\n";
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
	opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to get copy of a Material\n";
	exit(-1);
      }

      Fiber_Tangent[i] = 0;
      Fiber_ElongP[i] = 0;      
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
FiberSection2dThermal::FiberSection2dThermal(int tag, int num, bool compCentroid):
  SectionForceDeformation(tag, SEC_TAG_FiberSection2dThermal),
  numFibers(0), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(compCentroid),
  sectionIntegr(0), e(2), eCommit(2), s(0), ks(0),
  DataMixed(27), sT(2), Fiber_Tangent(0), Fiber_ElongP(0), AverageThermalElong(2),
  dedh(2)
{
    if(sizeFibers > 0) {
	theMaterials = new UniaxialMaterial *[sizeFibers];

	if(theMaterials == 0) {
	    opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to allocate Material pointers";
	    exit(-1);
	}

	matData = new double [sizeFibers*2];

	if(matData == 0) {
	    opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to allocate double array for material data\n";
	    exit(-1);
	}


	Fiber_Tangent = new double [sizeFibers];
	Fiber_ElongP = new double [sizeFibers];
	if (Fiber_Tangent == 0 || Fiber_ElongP == 0) {
	    opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to allocate double array for fiber data data\n";
	    exit(-1);
	}

	for(int i = 0; i < sizeFibers; i++) {
	    matData[i*2] = 0.0;
	    matData[i*2+1] = 0.0;
	    theMaterials[i] = 0;
	    Fiber_Tangent[i] = 0.0;
	    Fiber_ElongP[i] = 0.0;
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


FiberSection2dThermal::FiberSection2dThermal(int tag, int num, UniaxialMaterial **mats,
					     SectionIntegration &si, bool compCentroid):
  SectionForceDeformation(tag, SEC_TAG_FiberSection2dThermal),
  numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(compCentroid),
  sectionIntegr(0), e(2), eCommit(2), s(0), ks(0),
  DataMixed(27), sT(2), Fiber_Tangent(0), Fiber_ElongP(0), AverageThermalElong(2),
  dedh(2)//,theTemperature(0)
{
  if (numFibers > 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to allocate Material pointers";
      exit(-1);
    }
    matData = new double [numFibers*2];

    if (matData == 0) {
      opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to allocate double array for material data\n";
      exit(-1);
    }

    Fiber_Tangent = new double [numFibers];
    Fiber_ElongP = new double [numFibers];
    if (Fiber_Tangent == 0 || Fiber_ElongP == 0) {
      opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to allocate double array for fiber data data\n";
      exit(-1);
    }    
  }

  sectionIntegr = si.getCopy();
  if (sectionIntegr == 0) {
    opserr << "Error: FiberSection2dThermal::FiberSection2dThermal: could not create copy of section integration object" << endln;
    exit(-1);
  }

  double fiberLocs[10000];
  sectionIntegr->getFiberLocations(numFibers, fiberLocs);

  double fiberArea[10000];
  sectionIntegr->getFiberWeights(numFibers, fiberArea);

  for (int i = 0; i < numFibers; i++) {

    ABar  += fiberArea[i];
    QzBar += fiberLocs[i]*fiberArea[i];

    theMaterials[i] = mats[i]->getCopy();

    if (theMaterials[i] == 0) {
      opserr << "FiberSection2dThermal::FiberSection2dThermal -- failed to get copy of a Material\n";
      exit(-1);
    }

    Fiber_Tangent[i] = 0.0;
    Fiber_ElongP[i] = 0.0;
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
FiberSection2dThermal::FiberSection2dThermal():
  SectionForceDeformation(0, SEC_TAG_FiberSection2dThermal),
  numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(true),
  sectionIntegr(0), e(2), eCommit(2), s(0), ks(0),
  DataMixed(27), Fiber_Tangent(0), Fiber_ElongP(0), AverageThermalElong(2),
  dedh(2)//, theTemperatures(0),theTemperatureFactor(0)
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
FiberSection2dThermal::addFiber(Fiber &newFiber)
{
    // need to create larger arrays
    if(numFibers == sizeFibers) {
	int newsize = 2*sizeFibers;
	if(newsize == 0) newsize = 30;
	UniaxialMaterial **newArray = new UniaxialMaterial *[newsize];
	double *newMatData = new double [2 * newsize];
	double *newFiberElongP = new double [newsize];
	double *newFiberTangent = new double [newsize];
	if (newArray == 0 || newMatData == 0 || newFiberElongP == 0 || newFiberTangent == 0) {
	    opserr <<"FiberSection2dThermal::addFiber -- failed to allocate Fiber pointers\n";
	    return -1;
	}

	// copy the old pointers and data
	int i;
	for (i = 0; i < sizeFibers; i++) {
	    newArray[i] = theMaterials[i];
	    newMatData[2*i] = matData[2*i];
	    newMatData[2*i+1] = matData[2*i+1];
	    newFiberElongP[i] = Fiber_ElongP[i];
	    newFiberTangent[i] = Fiber_Tangent[i];	    
	}

	// initialize new memory
	for(i = sizeFibers; i<newsize; i++) {
	    newArray[i] = 0;
	    newMatData[2*i] = 0.0;
	    newMatData[2*i+1] = 0.0;
	    newFiberElongP[i] = 0.0;
	    newFiberTangent[i] = 0.0;
	}

	sizeFibers = newsize;

	// set new memory
	if (theMaterials != 0)
	    delete [] theMaterials;
	if (matData != 0)
	    delete [] matData;
	if (Fiber_ElongP != 0)
	  delete [] Fiber_ElongP;
	if (Fiber_Tangent != 0)
	  delete [] Fiber_Tangent;	

	theMaterials = newArray;
	matData = newMatData;
	Fiber_ElongP = newFiberElongP;
	Fiber_Tangent = newFiberTangent;
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
	opserr <<"FiberSection2dThermal::addFiber -- failed to get copy of a Material\n";
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
FiberSection2dThermal::~FiberSection2dThermal()
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

  if (Fiber_Tangent != 0)
    delete [] Fiber_Tangent;

  if (Fiber_ElongP != 0)
    delete [] Fiber_ElongP;
}



int
FiberSection2dThermal::setTrialSectionDeformation(const Vector& deforms)
{

  int res = 0;
  e = deforms;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;

  double d0 = deforms(0);
  double d1 = deforms(1);

  //d0 = d0 - 0.0087; //JZ

  double fiberLocs[10000];
  double fiberArea[10000];

  if (sectionIntegr != 0)
  {
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
    sectionIntegr->getFiberWeights(numFibers, fiberArea);
  }
  else
  {
    for (int i = 0; i < numFibers; i++) {
      fiberLocs[i] = matData[2*i];
      fiberArea[i] = matData[2*i+1];
    }
  }

	for (int i = 0; i < numFibers; i++) {

		// initializing material strain and set it
		UniaxialMaterial *theMat = theMaterials[i];
		double tangent =0.0;
		double ThermalElongation = 0.0;
		double FiberTemperature = 0;
		double FiberTempMax = 0;
		if ( fabs(DataMixed(1)) <= 1e-10 && fabs(DataMixed(17)) <= 1e-10 ) //no tempe load
		{
			FiberTemperature = 0;
			FiberTempMax=0;
		}
		else
		{
			//calculate the fiber tempe, T=T1-(Y-Y1)*(T1-T2)/(Y1-Y2)
			Vector TempV= this->determineFiberTemperature( DataMixed, fiberLocs[i]);
			FiberTemperature = TempV(0);
			FiberTempMax= TempV(1);
		}
		// get the data from thermal material
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
		Fiber_Tangent[i]=tangent;

		double y = fiberLocs[i] - yBar;
		double A = fiberArea[i];
		double stress = 0.0;
		double strain = d0 - y*d1;   //axial strain d0, rotational degree d1;

		strain = strain - ThermalElongation;  //Mechanical strain calculated by subtracting total strain with thermal strain

		//opserr<<"Total strain "<<strain+ThermalElongation<<"Updated mechanical strain "<<strain <<endln;

		res += theMat->setTrial(strain, FiberTemperature, stress, tangent, ThermalElongation);//***JZ

		Fiber_Tangent[i]=tangent;
		double ks0 = tangent * A;
		double ks1 = ks0 * -y;
		kData[0] += ks0;
		kData[1] += ks1;
		kData[3] += ks1 * -y;

		//force and temperature load
		double fs0 = stress * A ;//Stress resultant
		sData[0] += fs0;
		sData[1] += fs0 * -y;
  }

  kData[2] = kData[1];

  return res;
}


const Vector&
FiberSection2dThermal::getSectionDeformation(void)
{
  return e;
}

const Matrix&
FiberSection2dThermal::getInitialTangent(void)
{
  static double kInitial[4];
  static Matrix kInitialMatrix(kInitial, 2, 2);
  kInitial[0] = 0.0; kInitial[1] = 0.0; kInitial[2] = 0.0; kInitial[3] = 0.0;

  double fiberLocs[10000];
  double fiberArea[10000];

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
FiberSection2dThermal::getSectionTangent(void)
{
  return *ks;
}

const Vector&
FiberSection2dThermal::getStressResultant(void)
{
  return *s;
}


//by UoE, this member function is used when applying thermal action loading
// Here we update the fiber temperature but use the last committed fiber material properties to generate thermal forces
const Vector&
FiberSection2dThermal::getTemperatureStress(const Vector &dataMixed)
{

  AverageThermalElong.Zero();
  DataMixed = dataMixed;

  double fiberLocs[10000];
  double fiberArea[10000];
  sT.Zero();
  
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


  //------updata fiber initial Modulus corresponding to its temperature----------------

  double DeltaThermalElong[1000];
  for( int i=0; i< numFibers; i++) {
	DeltaThermalElong[i]=0;
	UniaxialMaterial *theMat = theMaterials[i];
    //Updating the fibre temperature  ---UoE Group
    double FiberTemperature = 0 ;
    double FiberTempMax=0; //PK add for max temp
    //if locY1 and locY9 are not less than zoro

    if ( fabs(dataMixed(1)) <= 1e-10 && fabs(dataMixed(17)) <= 1e-10 ) //no tempe load
    {
		FiberTemperature = 0;
		FiberTempMax=0;
    }
    else
	{
		//calculate the fiber tempe, T=T1-(Y-Y1)*(T1-T2)/(Y1-Y2)
		Vector TempV= this->determineFiberTemperature( dataMixed, fiberLocs[i]);
		FiberTemperature = TempV(0);FiberTempMax= TempV(1);
	}

    // obtaining new thermal Elongation
    double tangent =0.0;
	double ThermalElongation =0.0;
    static Vector tData(4);
    static Information iData(tData);
    tData(0) = FiberTemperature;
	tData(1) = tangent;
	tData(2) = ThermalElongation;
    tData(3) = FiberTempMax;
    iData.setVector(tData);
    theMat->getVariable("ElongTangent", iData);   //Actually here update initial tangent and thermalElongation corresponding  to current temperature
    tData = iData.getData();
    tangent = tData(1);
    ThermalElongation = tData(2);

	DeltaThermalElong[i]= ThermalElongation-Fiber_ElongP[i];
	Fiber_ElongP[i]= ThermalElongation;
  }

 // calculate section resisting force due to thermal load
  double FiberForce;
  double SectionArea=0;
  double ThermalForce=0;
  double ThermalMoment =0;
  double SectionMomofArea =0;
  for (int i = 0; i < numFibers; i++) {
	  FiberForce = Fiber_Tangent[i]*fiberArea[i]*DeltaThermalElong[i];
	  SectionArea +=fiberArea[i];
	  SectionMomofArea += (fiberArea[i]*(fiberLocs[i] - yBar)*(fiberLocs[i] - yBar));
	  ThermalForce += Fiber_ElongP[i]*fiberArea[i];
	  ThermalMoment += Fiber_ElongP[i]*fiberArea[i]*(fiberLocs[i] - yBar);
	  sT(0) += FiberForce;
	  sT(1) -= FiberForce*(fiberLocs[i] - yBar);
  }
  AverageThermalElong(0) = ThermalForce/SectionArea;
  AverageThermalElong(1) = ThermalMoment/SectionMomofArea;
  return sT;

}
//UoE group///Calculating Thermal stresses at each /////////////////////////////////////////////////////end
const Vector&
FiberSection2dThermal::getThermalElong(void)
{
  return AverageThermalElong;
}
//Retuning ThermalElongation



SectionForceDeformation*
FiberSection2dThermal::getCopy(void)
{
  FiberSection2dThermal *theCopy = new FiberSection2dThermal ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;

  if (numFibers > 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr <<"FiberSection2dThermal::getCopy -- failed to allocate Material pointers\n";
      exit(-1);
    }

    theCopy->matData = new double [numFibers*2];

    if (theCopy->matData == 0) {
      opserr << "FiberSection2dThermal::getCopy -- failed to allocate double array for material data\n";
      exit(-1);
    }

    theCopy->Fiber_Tangent = new double [numFibers];
    theCopy->Fiber_ElongP = new double [numFibers];
    if (theCopy->Fiber_Tangent == 0 || theCopy->Fiber_ElongP == 0) {
      opserr << "FiberSection2dThermal::getCopy -- failed to allocate double array for fiber data\n";
      exit(-1);
    }
    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*2] = matData[i*2];
      theCopy->matData[i*2+1] = matData[i*2+1];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr <<"FiberSection2dThermal::getCopy -- failed to get copy of a Material";
	exit(-1);
      }
      theCopy->Fiber_Tangent[i] = Fiber_Tangent[i];
      theCopy->Fiber_ElongP[i] = Fiber_ElongP[i];
    }
    //retrieve temperatures
	//theCopy->theTemperatures = theTemperatures->getCopy();

//	if (theCopy->theTemperatures == 0) {
//		opserr <<"FiberSection2dThermal::getCopy -- failed to get copy of the temperatures";
//		exit(-1);
//	}

  }


  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->sT = sT;

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
FiberSection2dThermal::getType ()
{
  return code;
}

int
FiberSection2dThermal::getOrder () const
{
  return 2;
}

int
FiberSection2dThermal::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  eCommit = e;

  return err;
}

int
FiberSection2dThermal::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;


  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;

  double fiberLocs[10000];
  double fiberArea[10000];

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
FiberSection2dThermal::revertToStart(void)
{
  // revert the fibers to start
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;

  double fiberLocs[10000];
  double fiberArea[10000];

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
FiberSection2dThermal::sendSelf(int commitTag, Channel &theChannel)
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
    opserr <<  "FiberSection2dThermal::sendSelf - failed to send ID data\n";
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
      opserr <<  "FiberSection2dThermal::sendSelf - failed to send material data\n";
      return res;
    }

    // send the fiber data, i.e. area and loc, tangent, and elongP
    Vector fiberData(4*numFibers);
    for (int i = 0; i < numFibers; i++) {
      fiberData(            i) = matData[i];
      fiberData(  numFibers+i) = matData[numFibers+i];
      fiberData(2*numFibers+i) = Fiber_Tangent[i];
      fiberData(3*numFibers+i) = Fiber_ElongP[i];
    }
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "FiberSection2dThermal::sendSelf - failed to send material data\n";
      return res;
    }

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);

  }

  return res;
}

int
FiberSection2dThermal::recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(3);

  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "FiberSection2dThermal::recvSelf - failed to recv ID data\n";
    return res;
  }
  this->setTag(data(0));

  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "FiberSection2dThermal::recvSelf - failed to recv material data\n";
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
      if (Fiber_Tangent != 0)
	delete [] Fiber_Tangent;
      if (Fiber_ElongP != 0)
	delete [] Fiber_ElongP;	
      Fiber_Tangent = 0;
      Fiber_ElongP = 0;	
      matData = 0;

      // create memory to hold material pointers and fiber data
      numFibers = data(1);
      if (numFibers != 0) {
	theMaterials = new UniaxialMaterial *[numFibers];

	if (theMaterials == 0) {
	  opserr <<"FiberSection2dThermal::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}

	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;

	matData = new double [numFibers*2];

	if (matData == 0) {
	  opserr <<"FiberSection2dThermal::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}

	Fiber_Tangent = new double [numFibers];
	if (Fiber_Tangent == 0) {
	  opserr <<"FiberSection2dThermal::recvSelf  -- failed to allocate double array for fiber tangent\n";
	  exit(-1);
	}
	Fiber_ElongP = new double [numFibers];
	if (Fiber_ElongP == 0) {
	  opserr <<"FiberSection2dThermal::recvSelf  -- failed to allocate double array for fiber elongation\n";
	  exit(-1);
	}		
      }
    }

    Vector fiberData(4*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "FiberSection2dThermal::recvSelf - failed to recv material data\n";
      return res;
    }
    for (int i = 0; i < numFibers; i++) {
      matData[i] =           fiberData(          i);
      matData[numFibers+i] = fiberData(numFibers+i);
      Fiber_Tangent[i] = fiberData(2*numFibers+i);
      Fiber_ElongP[i]  = fiberData(3*numFibers+i);
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
	opserr <<"FiberSection2dThermal::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    double Qz = 0.0;
    double A  = 0.0;
    double yLoc, Area;

    computeCentroid = data(2) ? true : false;
    
    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = matData[2*i];
      Area = matData[2*i+1];
      A  += Area;
      Qz += yLoc*Area;
    }

    if (computeCentroid)
      yBar = Qz/A;
    else
      yBar = 0.0;
  }

  return res;
}

void
FiberSection2dThermal::Print(OPS_Stream &s, int flag)
{
  s << "\nFiberSection2dThermal, tag: " << this->getTag() << endln;
  s << "\tSection code: " << code;
  s << "\tNumber of Fibers: " << numFibers << endln;
  s << "\tCentroid: " << yBar << endln;

  if (flag == 1) {
    for (int i = 0; i < numFibers; i++) {
      s << "\nLocation (y) = (" << matData[2*i] << ")";
      s << "\nArea = " << matData[2*i+1] << endln;
      theMaterials[i]->Print(s, flag);
    }
  }
}

Response*
FiberSection2dThermal::setResponse(const char **argv, int argc,
			    OPS_Stream &output)
{
  Response *theResponse = 0;
  
  if (argc > 2 && strcmp(argv[0],"fiber") == 0) {
    
    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 3) {		  // fiber number was input directly

      key = atoi(argv[1]);

    } else if (argc > 4) {  // find fiber closest to coord. with mat tag

      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      
      double closestDist;
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
      
      theResponse = theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }
  }

  if (theResponse == 0)
    return SectionForceDeformation::setResponse(argv, argc, output);

  return theResponse;
}


int
FiberSection2dThermal::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}



// AddingSensitivity:BEGIN ////////////////////////////////////
int
FiberSection2dThermal::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

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
FiberSection2dThermal::getSectionDeformationSensitivity(int gradIndex)
{
  static Vector dummy(2);

  return dummy;
}

const Vector &
FiberSection2dThermal::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(2);

  ds.Zero();

  double y, A, stressGradient, stress, tangent, sig_dAdh;

  double fiberLocs[10000];
  double fiberArea[10000];

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

  double locsDeriv[10000];
  double areaDeriv[10000];

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
FiberSection2dThermal::getInitialTangentSensitivity(int gradIndex)
{
  static Matrix dksdh(2,2);

  dksdh.Zero();

  double y, A, dydh, dAdh, tangent, dtangentdh;

  double fiberLocs[10000];
  double fiberArea[10000];

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

  double locsDeriv[10000];
  double areaDeriv[10000];

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
FiberSection2dThermal::commitSensitivity(const Vector& defSens,
				  int gradIndex, int numGrads)
{
  double d0 = defSens(0);
  double d1 = defSens(1);

  dedh = defSens;

  double fiberLocs[10000];

  if (sectionIntegr != 0)
    sectionIntegr->getFiberLocations(numFibers, fiberLocs);
  else {
    for (int i = 0; i < numFibers; i++)
      fiberLocs[i] = matData[2*i];
  }

  double locsDeriv[10000];
  double areaDeriv[10000];

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


const Vector&
FiberSection2dThermal::determineFiberTemperature(const Vector& DataMixed, double fiberLoc)
{
		double FiberTemperature = 0;
		double FiberTempMax = 0;

		double dataTempe[27]; //PK changed 18 to 27 to pass max temps
		for (int i = 0; i < 27; i++) {
			dataTempe[i] = DataMixed(i);
		}

		if (  fiberLoc <= dataTempe[1])
		{
			opserr <<"FiberSection2dThermal::setTrialSectionDeformationTemperature -- fiber loc is out of the section";
		}
		else if (fiberLoc <= dataTempe[3])
		{
			FiberTemperature = dataTempe[0] - (dataTempe[1] - fiberLoc) * (dataTempe[0] - dataTempe[2])/(dataTempe[1] - dataTempe[3]);
			//FiberTempMax = dataTempe[18] - (dataTempe[1] - fiberLoc) * (dataTempe[18] - dataTempe[19])/(dataTempe[1] - dataTempe[3]);
		}
		else if (   fiberLoc <= dataTempe[5] )
		{
			FiberTemperature = dataTempe[2] - (dataTempe[3] - fiberLoc) * (dataTempe[2] - dataTempe[4])/(dataTempe[3] - dataTempe[5]);
			//FiberTempMax = dataTempe[19] - (dataTempe[3] - fiberLoc) * (dataTempe[19] - dataTempe[20])/(dataTempe[3] - dataTempe[5]);
		}
		else if ( fiberLoc <= dataTempe[7] )
		{
			FiberTemperature = dataTempe[4] - (dataTempe[5] - fiberLoc) * (dataTempe[4] - dataTempe[6])/(dataTempe[5] - dataTempe[7]);
			//FiberTempMax = dataTempe[20] - (dataTempe[5] - fiberLoc) * (dataTempe[20] - dataTempe[21])/(dataTempe[5] - dataTempe[7]);
		}
		else if ( fiberLoc <= dataTempe[9] )
		{
			FiberTemperature = dataTempe[6] - (dataTempe[7] - fiberLoc) * (dataTempe[6] - dataTempe[8])/(dataTempe[7] - dataTempe[9]);
			//FiberTempMax = dataTempe[21] - (dataTempe[7] - fiberLoc) * (dataTempe[21] - dataTempe[22])/(dataTempe[7] - dataTempe[9]);
		}
		else if (fiberLoc <= dataTempe[11] )
		{
			FiberTemperature = dataTempe[8] - (dataTempe[9] - fiberLoc) * (dataTempe[8] - dataTempe[10])/(dataTempe[9] - dataTempe[11]);
			//FiberTempMax = dataTempe[22] - (dataTempe[9] - fiberLoc) * (dataTempe[22] - dataTempe[23])/(dataTempe[9] - dataTempe[11]);
		}
		else if (fiberLoc <= dataTempe[13] )
		{
			FiberTemperature = dataTempe[10] - (dataTempe[11] - fiberLoc) * (dataTempe[10] - dataTempe[12])/(dataTempe[11] - dataTempe[13]);
			//FiberTempMax = dataTempe[23] - (dataTempe[11] - fiberLoc) * (dataTempe[23] - dataTempe[24])/(dataTempe[11] - dataTempe[13]);
		}
		else if (fiberLoc <= dataTempe[15] )
		{
			FiberTemperature = dataTempe[12] - (dataTempe[13] - fiberLoc) * (dataTempe[12] - dataTempe[14])/(dataTempe[13] - dataTempe[15]);
			//FiberTempMax = dataTempe[24] - (dataTempe[13] - fiberLoc) * (dataTempe[24] - dataTempe[25])/(dataTempe[13] - dataTempe[15]);
		}
		else if ( fiberLoc <= dataTempe[17] )
		{
			FiberTemperature = dataTempe[14] - (dataTempe[15] - fiberLoc) * (dataTempe[14] - dataTempe[16])/(dataTempe[15] - dataTempe[17]);
			//FiberTempMax = dataTempe[25] - (dataTempe[15] - fiberLoc) * (dataTempe[25] - dataTempe[26])/(dataTempe[15] - dataTempe[17]);
		}
		else
		{
			opserr <<"FiberSection2dThermal::setTrialSectionDeformation -- fiber loc is out of the section";
		}

		static Vector returnedTemperature(2);
		returnedTemperature(0)=FiberTemperature;
		returnedTemperature(1)=FiberTempMax;
		return returnedTemperature;
}
