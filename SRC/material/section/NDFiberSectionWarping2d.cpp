
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
#include <NDFiberSectionWarping2d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <NDMaterial.h>
#include <SectionIntegration.h>
#include <Parameter.h>
#include <elementAPI.h>

ID NDFiberSectionWarping2d::code(5);

void* OPS_NDFiberSectionWarping2d()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
	opserr<<"insufficient arguments for NDFiberSectionWarping2d\n";
	return 0;
    }

    numData = 1;
    int tag;
    if (OPS_GetIntInput(&numData,&tag) < 0) return 0;

    double alpha = 1.0;
    bool computeCentroid = true;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      const char* opt = OPS_GetString();
      if (strcmp(opt, "-noCentroid") == 0)
	computeCentroid = false;
      if (strcmp(opt, "-alpha") == 0 || strcmp(opt, "-shape") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 1)
	  break;
	numData = 1;
	if (OPS_GetDoubleInput(&numData,&alpha) < 0)
	  return 0;
      }
    }
    
    int num = 30;
    return new NDFiberSectionWarping2d(tag, num, alpha);
}

// constructors:
NDFiberSectionWarping2d::NDFiberSectionWarping2d(int tag, int num, Fiber **fibers, double a): 
SectionForceDeformation(tag, SEC_TAG_NDFiberSectionWarping2d),
numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
yBar(0.0), alpha(a), yBarZero(0.0), DeltaYbar(0.0), sectionIntegr(0),
e(5), eCommit(5), s(0), ks(0), parameterID(0), dedh(5)
{ 
	if (numFibers != 0) {
		theMaterials = new NDMaterial *[numFibers]; 

		if (theMaterials == 0) {
			opserr << "NDFiberSectionWarping2d::NDFiberSectionWarping2d -- failed to allocate Material pointers";
			exit(-1);
		}

		matData = new double [numFibers*2];

		if (matData == 0) {
			opserr << "NDFiberSectionWarping2d::NDFiberSectionWarping2d -- failed to allocate double array for material data\n";
			exit(-1);
		}


		double Qz = 0.0;
		double A  = 0.0;

		for (int i = 0; i < numFibers; i++) {
			Fiber *theFiber = fibers[i];
			double yLoc, zLoc, Area;
			theFiber->getFiberLocation(yLoc, zLoc);
			Area = theFiber->getArea();

			NDMaterial *theMat = theFiber->getNDMaterial();

			A  += Area;
			Qz += yLoc * Area;
			matData[i*2] = yLoc;
			matData[i*2+1] = Area;

			theMaterials[i] = theMat->getCopy("BeamFiber2d"); 

			if (theMaterials[i] == 0) {
				opserr << "NDFiberSectionWarping2d::NDFiberSectionWarping2d -- failed to get copy of a Material\n";
				exit(-1);
			}
		}    

		yBar = Qz/A; 
		yBarZero = Qz/A;
	}

	s = new Vector(sData, 5);
	ks = new Matrix(kData, 5, 5);

	sData[0] = 0.0;
	sData[1] = 0.0;
	sData[2] = 0.0;
	sData[3] = 0.0;
	sData[4] = 0.0;

	kData[0] = 0.0;
	kData[1] = 0.0;
	kData[2] = 0.0;
	kData[3] = 0.0;
	kData[4] = 0.0;
	kData[5] = 0.0;
	kData[6] = 0.0;
	kData[7] = 0.0;
	kData[8] = 0.0;
	kData[9] = 0.0;
	kData[10] = 0.0;
	kData[11] = 0.0;
	kData[12] = 0.0;
	kData[13] = 0.0;
	kData[14] = 0.0;
	kData[15] = 0.0;
	kData[16] = 0.0;
	kData[17] = 0.0;
	kData[18] = 0.0;
	kData[19] = 0.0;
	kData[20] = 0.0;
	kData[21] = 0.0;
	kData[22] = 0.0;
	kData[23] = 0.0;
	kData[24] = 0.0;

	code(0) = SECTION_RESPONSE_P;
	code(1) = SECTION_RESPONSE_MZ;
	code(2) = SECTION_RESPONSE_VY;
	code(3) = SECTION_RESPONSE_R;
	code(4) = SECTION_RESPONSE_Q;
}

NDFiberSectionWarping2d::NDFiberSectionWarping2d(int tag, int num, double a): 
    SectionForceDeformation(tag, SEC_TAG_NDFiberSectionWarping2d),
    numFibers(0), sizeFibers(num), theMaterials(0), matData(0),
    yBar(0.0), alpha(a), yBarZero(0.0), DeltaYbar(0.0), sectionIntegr(0),
    e(5), eCommit(5), s(0), ks(0), parameterID(0), dedh(5)
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

    int order = 5;
    s = new Vector(sData, order);
    ks = new Matrix(kData, order, order);

    for (int i = 0; i < order; i++) {
      sData[i] = 0.0;
      for (int j = 0; j < order; j++)
	kData[5*i + j] = 0.0;
    }
    
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_VY;
    code(3) = SECTION_RESPONSE_R;
    code(4) = SECTION_RESPONSE_Q;    
}

NDFiberSectionWarping2d::NDFiberSectionWarping2d(int tag, int num, NDMaterial **mats,
	SectionIntegration &si, double a):
SectionForceDeformation(tag, SEC_TAG_NDFiberSectionWarping2d),
numFibers(num), sizeFibers(num), theMaterials(0), matData(0),
yBar(0.0), alpha(a), yBarZero(0.0), DeltaYbar(0.0), sectionIntegr(0),
e(5), eCommit(5), s(0), ks(0), parameterID(0), dedh(5)
{
	if (numFibers != 0) {
		theMaterials = new NDMaterial *[numFibers];

		if (theMaterials == 0) {
			opserr << "NDFiberSectionWarping2d::NDFiberSectionWarping2d -- failed to allocate Material pointers";
			exit(-1);
		}
		matData = new double [numFibers*2];

		if (matData == 0) {
			opserr << "NDFiberSectionWarping2d::NDFiberSectionWarping2d -- failed to allocate double array for material data\n";
			exit(-1);
		}
	}

	sectionIntegr = si.getCopy();
	if (sectionIntegr == 0) {
		opserr << "Error: NDFiberSectionWarping2d::NDFiberSectionWarping2d: could not create copy of section integration object" << endln;
		exit(-1);
	}

	static double fiberLocs[10000];
	sectionIntegr->getFiberLocations(numFibers, fiberLocs);

	static double fiberArea[10000];
	sectionIntegr->getFiberWeights(numFibers, fiberArea);

	double Qz = 0.0;
	double A  = 0.0;

	for (int i = 0; i < numFibers; i++) {

		A  +=  fiberArea[i];
		Qz +=  fiberLocs[i]*fiberArea[i];

		theMaterials[i] = mats[i]->getCopy("BeamFiber2d");

		if (theMaterials[i] == 0) {
			opserr << "NDFiberSectionWarping2d::NDFiberSectionWarping2d -- failed to get copy of a Material\n";
			exit(-1);
		}
	}    

	yBar = Qz/A; 
	yBarZero = Qz/A;

	s = new Vector(sData, 5);
	ks = new Matrix(kData, 5, 5); 

	sData[0] = 0.0;
	sData[1] = 0.0;
	sData[2] = 0.0;
	sData[3] = 0.0;
	sData[4] = 0.0;

	kData[0] = 0.0;
	kData[1] = 0.0;
	kData[2] = 0.0;
	kData[3] = 0.0;
	kData[4] = 0.0;
	kData[5] = 0.0;
	kData[6] = 0.0;
	kData[7] = 0.0;
	kData[8] = 0.0;
	kData[9] = 0.0;
	kData[10] = 0.0;
	kData[11] = 0.0;
	kData[12] = 0.0;
	kData[13] = 0.0;
	kData[14] = 0.0;
	kData[15] = 0.0;
	kData[16] = 0.0;
	kData[17] = 0.0;
	kData[18] = 0.0;
	kData[19] = 0.0;
	kData[20] = 0.0;
	kData[21] = 0.0;
	kData[22] = 0.0;
	kData[23] = 0.0;
	kData[24] = 0.0;

	code(0) = SECTION_RESPONSE_P;
	code(1) = SECTION_RESPONSE_MZ;
	code(2) = SECTION_RESPONSE_VY;
	code(3) = SECTION_RESPONSE_R;
	code(4) = SECTION_RESPONSE_Q;
}

// constructor for blank object that recvSelf needs to be invoked upon
NDFiberSectionWarping2d::NDFiberSectionWarping2d():
SectionForceDeformation(0, SEC_TAG_NDFiberSectionWarping2d),
numFibers(0), sizeFibers(0), theMaterials(0), matData(0),
yBar(0.0), alpha(0.0), yBarZero(0.0), DeltaYbar(0.0), sectionIntegr(0),
e(5), eCommit(5), s(0), ks(0), parameterID(0), dedh(5)
{
	s = new Vector(sData, 5);
	ks = new Matrix(kData, 5, 5);

	sData[0] = 0.0;
	sData[1] = 0.0;
	sData[2] = 0.0;
	sData[3] = 0.0;
	sData[4] = 0.0;

	kData[0] = 0.0;
	kData[1] = 0.0;
	kData[2] = 0.0;
	kData[3] = 0.0;
	kData[4] = 0.0;
	kData[5] = 0.0;
	kData[6] = 0.0;
	kData[7] = 0.0;
	kData[8] = 0.0;
	kData[9] = 0.0;
	kData[10] = 0.0;
	kData[11] = 0.0;
	kData[12] = 0.0;
	kData[13] = 0.0;
	kData[14] = 0.0;
	kData[15] = 0.0;
	kData[16] = 0.0;
	kData[17] = 0.0;
	kData[18] = 0.0;
	kData[19] = 0.0;
	kData[20] = 0.0;
	kData[21] = 0.0;
	kData[22] = 0.0;
	kData[23] = 0.0;
	kData[24] = 0.0;

	code(0) = SECTION_RESPONSE_P;
	code(1) = SECTION_RESPONSE_MZ;
	code(2) = SECTION_RESPONSE_VY;
	code(3) = SECTION_RESPONSE_R;
	code(4) = SECTION_RESPONSE_Q;
}

int
NDFiberSectionWarping2d::addFiber(Fiber &newFiber)
{
	// need to create larger arrays
  if (numFibers == sizeFibers) {
	int newSize = numFibers+1;
	NDMaterial **newArray = new NDMaterial *[newSize]; 
	double *newMatData = new double [2 * newSize];
	if (newArray == 0 || newMatData == 0) {
		opserr <<"NDFiberSectionWarping2d::addFiber -- failed to allocate Fiber pointers\n";
		return -1;
	}

	// copy the old pointers and data
	int i;
	for (i = 0; i < numFibers; i++) {
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
    opserr <<"NDFiberSectionWarping2d::addFiber -- failed to get copy of a Material\n";
    return -1;
  }

  numFibers++;

	double Qz = 0.0;
	double A  = 0.0;

	// Recompute centroid
	for (int i = 0; i < numFibers; i++) {
		yLoc = -matData[2*i];
		Area = matData[2*i+1];
		A  += Area;
		Qz += yLoc*Area;
	}

	yBar = Qz/A;
	yBarZero = Qz/A;

	return 0;
}


// destructor:
NDFiberSectionWarping2d::~NDFiberSectionWarping2d()
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
	NDFiberSectionWarping2d::setTrialSectionDeformation (const Vector &deforms)
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
	kData[9] = 0.0;
	kData[10] = 0.0;
	kData[11] = 0.0;
	kData[12] = 0.0;
	kData[13] = 0.0;
	kData[14] = 0.0;
	kData[15] = 0.0;
	kData[16] = 0.0;
	kData[17] = 0.0;
	kData[18] = 0.0;
	kData[19] = 0.0;
	kData[20] = 0.0;
	kData[21] = 0.0;
	kData[22] = 0.0;
	kData[23] = 0.0;
	kData[24] = 0.0;

	sData[0] = 0.0;
	sData[1] = 0.0;
	sData[2] = 0.0;
	sData[3] = 0.0;
	sData[4] = 0.0;

	double d0 = deforms(0);
	double d1 = deforms(1);
	double d2 = deforms(2);
	double d3 = deforms(3);
	double d4 = deforms(4);

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

	// h ~ parameter equals Height/2 for symmetric cases
	double maxLoc(fiberLocs[1] - yBarZero), minLoc(fiberLocs[1] - yBarZero); 
	for (int i = 0; i < numFibers; i++) {
		if (fiberLocs[i] - yBarZero > maxLoc) { maxLoc = fiberLocs[i] - yBarZero;}
		if (fiberLocs[i] - yBarZero < minLoc) { minLoc = fiberLocs[i] - yBarZero;}
	}
	double h (maxLoc);  //opserr<<"   h:   "<<h<<endln;

	for (int i = 0; i < numFibers; i++) {

		NDMaterial *theMat = theMaterials[i];
		double y = fiberLocs[i] - yBar; 
		double A = fiberArea[i];

		// determine material strain and set it
		eps(0) = d0 - y*d1 + (y*y*y/(h*h*h) - 0.6*y/h)*d4;

		double rootAlpha = 1.0;
		eps(1) = d2 + (3*y*y/(h*h*h) - 0.6/h)*d3; 

		if (alpha != 1.0) {
			rootAlpha = sqrt(alpha);
			eps(1) *= rootAlpha;
		}

		res += theMat->setTrialStrain(eps);
		const Vector &stress = theMat->getStress();
		const Matrix &tangent = theMat->getTangent();

		double d00 = tangent(0,0)*A;
		double d01 = tangent(0,1)*A;
		double d10 = tangent(1,0)*A;
		double d11 = tangent(1,1)*A; 

		double omega      = y*y*y/(h*h*h) - 0.6*y/h; 
		double omegaprime = 3*y*y/(h*h*h) - 0.6/h;       

		kData[0] += d00;
		kData[1] += -y * d00;
		kData[2] += d01;
		kData[3] += omegaprime * d01;
		kData[4] += omega * d00;

		kData[5] += -y * d00;
		kData[6] += y * y * d00;
		kData[7] += -y * d01;
		kData[8] += -y * omegaprime * d01;
		kData[9] += -y * omega * d00;

		kData[10] += d10;
		kData[11] += -y * d10;
		kData[12] += d11;
		kData[13] += omegaprime * d11;
		kData[14] += omega * d10;

		kData[15] += omegaprime * d10;
		kData[16] += -y * omegaprime * d10;
		kData[17] += omegaprime * d11;
		kData[18] += omegaprime * omegaprime * d11;
		kData[19] += omega * omegaprime * d10;

		kData[20] += omega * d00;
		kData[21] += -y * omega * d00;
		kData[22] += omega * d01;
		kData[23] += omegaprime * omega * d01;
		kData[24] += omega * omega * d00;

		double fs0 = stress(0)*A;
		double fs1 = stress(1)*A;

		sData[0] += fs0;
		sData[1] += fs0*-y;
		sData[2] += fs1;
		sData[3] += fs1*omegaprime;
		sData[4] += fs0*omega;
	}

	double rootAlpha = 1.0;

	if (alpha != 1.0) {
		rootAlpha = sqrt(alpha);
		eps(1) *= rootAlpha;
	}

	if (alpha != 1.0) {
		sData[2] *= rootAlpha;
		sData[3] *= rootAlpha;

		kData[2] *= rootAlpha;
		kData[3] *= rootAlpha;

		kData[7] *= rootAlpha;
		kData[8] *= rootAlpha;

		kData[10] *= rootAlpha;
		kData[11] *= rootAlpha;
		kData[14] *= rootAlpha;

		kData[15] *= rootAlpha;
		kData[16] *= rootAlpha;
		kData[19] *= rootAlpha;

		kData[22] *= rootAlpha;
		kData[23] *= rootAlpha;

		kData[12] *= alpha;
		kData[13] *= alpha;

		kData[17] *= alpha;
		kData[18] *= alpha;

	}

	return res;
}

const Vector&
	NDFiberSectionWarping2d::getSectionDeformation(void)
{
	return e;
}

double
	NDFiberSectionWarping2d::getRho(void)
{
	return 0.0;
}

const Matrix&
	NDFiberSectionWarping2d::getInitialTangent(void)
{
	static double kInitial[25];
	static Matrix kInitialMatrix(kInitial, 5, 5);
	kInitial[0] = 0.0; 
	kInitial[1] = 0.0;
	kInitial[2] = 0.0;
	kInitial[3] = 0.0;
	kInitial[4] = 0.0;
	kInitial[5] = 0.0;
	kInitial[6] = 0.0;
	kInitial[7] = 0.0;
	kInitial[8] = 0.0;
	kInitial[9] = 0.0;
	kInitial[10] = 0.0;
	kInitial[11] = 0.0;
	kInitial[12] = 0.0;
	kInitial[13] = 0.0;
	kInitial[14] = 0.0;
	kInitial[15] = 0.0;
	kInitial[16] = 0.0;
	kInitial[17] = 0.0;
	kInitial[18] = 0.0;
	kInitial[19] = 0.0;
	kInitial[20] = 0.0;
	kInitial[21] = 0.0;
	kInitial[22] = 0.0;
	kInitial[23] = 0.0;
	kInitial[24] = 0.0;

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

		// h ~ section height / 2 in Linear case
		double maxLoc(fiberLocs[1] - yBarZero), minLoc(fiberLocs[1] - yBarZero);
		for (int i = 0; i < numFibers; i++) {
			if (fiberLocs[i] - yBarZero > maxLoc) { maxLoc = fiberLocs[i] - yBarZero;}
			if (fiberLocs[i] - yBarZero < minLoc) { minLoc = fiberLocs[i] - yBarZero;}
		}
		double h (maxLoc);
		double omega      = y*y*y/(h*h*h)-0.6*y/h;
		double omegaprime = 3*y*y/(h*h*h)-0.6/h;

		kInitial[0] += d00;
		kInitial[1] += -y * d00;
		kInitial[2] += d01;
		kInitial[3] += omegaprime * d01;
		kInitial[4] += omega * d00;

		kInitial[5] += -y * d00;
		kInitial[6] += y * y * d00;
		kInitial[7] += -y * d01;
		kInitial[8] += -y * omegaprime * d01;
		kInitial[9] += -y * omega * d00;

		kInitial[10] += d10;
		kInitial[11] += -y * d10;
		kInitial[12] += d11;
		kInitial[13] += omegaprime * d11;
		kInitial[14] += omega * d10;

		kInitial[15] += omegaprime * d10;
		kInitial[16] += -y * omegaprime * d10;
		kInitial[17] += omegaprime * d11;
		kInitial[18] += omegaprime * omegaprime * d11;
		kInitial[19] += omega * omegaprime * d10;

		kInitial[20] += omega * d00;
		kInitial[21] += -y * omega * d00;
		kInitial[22] += omega * d01;
		kInitial[23] += omegaprime * omega * d01;
		kInitial[24] += omega * omega * d00;
	}

	if (alpha != 1.0) {
		double rootAlpha = sqrt(alpha);

		kInitial[2] *= rootAlpha;
		kInitial[3] *= rootAlpha;

		kInitial[7] *= rootAlpha;
		kInitial[8] *= rootAlpha;

		kInitial[10] *= rootAlpha;
		kInitial[11] *= rootAlpha;
		kInitial[14] *= rootAlpha;

		kInitial[15] *= rootAlpha;
		kInitial[16] *= rootAlpha;
		kInitial[19] *= rootAlpha;

		kInitial[22] *= rootAlpha;
		kInitial[23] *= rootAlpha;

		kInitial[12] *= alpha;
		kInitial[13] *= alpha;

		kInitial[17] *= alpha;
		kInitial[18] *= alpha;
	}

	return kInitialMatrix;
}

const Matrix&
	NDFiberSectionWarping2d::getSectionTangent(void)
{
	return *ks;
}

const Vector&
	NDFiberSectionWarping2d::getStressResultant(void)
{
	return *s;
}

SectionForceDeformation*
	NDFiberSectionWarping2d::getCopy(void)
{
	NDFiberSectionWarping2d *theCopy = new NDFiberSectionWarping2d ();
	theCopy->setTag(this->getTag());

	theCopy->numFibers = numFibers;
	theCopy->sizeFibers = numFibers;	

	if (numFibers != 0) {
		theCopy->theMaterials = new NDMaterial *[numFibers];

		if (theCopy->theMaterials == 0) {
			opserr <<"NDFiberSectionWarping2d::getCopy -- failed to allocate Material pointers\n";
			exit(-1);
		}

		theCopy->matData = new double [numFibers*2];

		if (theCopy->matData == 0) {
			opserr << "NDFiberSectionWarping2d::getCopy -- failed to allocate double array for material data\n";
			exit(-1);
		}

		for (int i = 0; i < numFibers; i++) {
			theCopy->matData[i*2] = matData[i*2];
			theCopy->matData[i*2+1] = matData[i*2+1];
			theCopy->theMaterials[i] = theMaterials[i]->getCopy("BeamFiber2d");

			if (theCopy->theMaterials[i] == 0) {
				opserr <<"NDFiberSectionWarping2d::getCopy -- failed to get copy of a Material";
				exit(-1);
			}
		}  
	}

	theCopy->eCommit = eCommit;
	theCopy->e = e;
	theCopy->yBar = yBar;
	theCopy->yBarZero = yBarZero;
	theCopy->DeltaYbar = DeltaYbar;	

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
	theCopy->kData[9] = kData[9];
	theCopy->kData[10] = kData[10];
	theCopy->kData[11] = kData[11];
	theCopy->kData[12] = kData[12];
	theCopy->kData[13] = kData[13];
	theCopy->kData[14] = kData[14];
	theCopy->kData[15] = kData[15];
	theCopy->kData[16] = kData[16];
	theCopy->kData[17] = kData[17];
	theCopy->kData[18] = kData[18];
	theCopy->kData[19] = kData[19];
	theCopy->kData[20] = kData[20];
	theCopy->kData[21] = kData[21];
	theCopy->kData[22] = kData[22];
	theCopy->kData[23] = kData[23];
	theCopy->kData[24] = kData[24];


	theCopy->sData[0] = sData[0];
	theCopy->sData[1] = sData[1];
	theCopy->sData[2] = sData[2];
	theCopy->sData[3] = sData[3];
	theCopy->sData[4] = sData[4];

	if (sectionIntegr != 0)
		theCopy->sectionIntegr = sectionIntegr->getCopy();
	else
		theCopy->sectionIntegr = 0;

	return theCopy;
}

const ID&
	NDFiberSectionWarping2d::getType ()
{
	return code;
}

int
	NDFiberSectionWarping2d::getOrder () const
{
	return 5;
}

int
	NDFiberSectionWarping2d::commitState(void)
{
	int err = 0;

	for (int i = 0; i < numFibers; i++)
		err += theMaterials[i]->commitState();

	eCommit = e;

	return err;
}

int
	NDFiberSectionWarping2d::revertToLastCommit(void)
{
	int err = 0;

	// Last committed section deformations
	e = eCommit;

	kData[0] = 0.0;
	kData[1] = 0.0;
	kData[2] = 0.0;
	kData[3] = 0.0;
	kData[4] = 0.0;
	kData[5] = 0.0;
	kData[6] = 0.0;
	kData[7] = 0.0;
	kData[8] = 0.0;
	kData[9] = 0.0;
	kData[10] = 0.0;
	kData[11] = 0.0;
	kData[12] = 0.0;
	kData[13] = 0.0;
	kData[14] = 0.0;
	kData[15] = 0.0;
	kData[16] = 0.0;
	kData[17] = 0.0;
	kData[18] = 0.0;
	kData[19] = 0.0;
	kData[20] = 0.0;
	kData[21] = 0.0;
	kData[22] = 0.0;
	kData[23] = 0.0;
	kData[24] = 0.0;

	sData[0] = 0.0;
	sData[1] = 0.0;
	sData[2] = 0.0;
	sData[3] = 0.0;
	sData[4] = 0.0;

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

	// h ~ section height / 2 in Linear case
	double maxLoc(fiberLocs[1] - yBarZero), minLoc(fiberLocs[1] - yBarZero);
	for (int i = 0; i < numFibers; i++) {
		if (fiberLocs[i] - yBarZero > maxLoc) { maxLoc = fiberLocs[i] - yBarZero;}
		if (fiberLocs[i] - yBarZero < minLoc) { minLoc = fiberLocs[i] - yBarZero;}
	}
	double h (maxLoc);

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

		double omega      = y*y*y/(h*h*h) - 0.6*y/h;
		double omegaprime = 3*y*y/(h*h*h) - 0.6/h;

		kData[0] += d00;
		kData[1] += -y * d00;
		kData[2] += d01;
		kData[3] += omegaprime * d01;
		kData[4] += omega * d00;

		kData[5] += -y * d00;
		kData[6] += y * y * d00;
		kData[7] += -y * d01;
		kData[8] += -y * omegaprime * d01;
		kData[9] += -y * omega * d00;

		kData[10] += d10;
		kData[11] += -y * d10;
		kData[12] += d11;
		kData[13] += omegaprime * d11;
		kData[14] += omega * d10;

		kData[15] += omegaprime * d10;
		kData[16] += -y * omegaprime * d10;
		kData[17] += omegaprime * d11;
		kData[18] += omegaprime * omegaprime * d11;
		kData[19] += omega * omegaprime * d10;

		kData[20] += omega * d00;
		kData[21] += -y * omega * d00;
		kData[22] += omega * d01;
		kData[23] += omegaprime * omega * d01;
		kData[24] += omega * omega * d00;

		double fs0 = stress(0)*A;
		double fs1 = stress(1)*A;

		sData[0] += fs0;
		sData[1] += fs0*-y;
		sData[2] += fs1;
		sData[3] += fs1*omegaprime;
		sData[4] += fs0*omega;

	}


	double rootAlpha = 1.0;

	if (alpha != 1.0) {
		rootAlpha = sqrt(alpha);

		sData[2] *= rootAlpha;
		sData[3] *= rootAlpha;

		kData[2] *= rootAlpha;
		kData[3] *= rootAlpha;

		kData[7] *= rootAlpha;
		kData[8] *= rootAlpha;

		kData[10] *= rootAlpha;
		kData[11] *= rootAlpha;
		kData[14] *= rootAlpha;

		kData[15] *= rootAlpha;
		kData[16] *= rootAlpha;
		kData[19] *= rootAlpha;

		kData[22] *= rootAlpha;
		kData[23] *= rootAlpha;

		kData[12] *= alpha;
		kData[13] *= alpha;

		kData[17] *= alpha;
		kData[18] *= alpha;
	}

	return err;
}

int
	NDFiberSectionWarping2d::revertToStart(void)
{
	// revert the fibers to start    
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
	kData[9] = 0.0;
	kData[10] = 0.0;
	kData[11] = 0.0;
	kData[12] = 0.0;
	kData[13] = 0.0;
	kData[14] = 0.0;
	kData[15] = 0.0;
	kData[16] = 0.0;
	kData[17] = 0.0;
	kData[18] = 0.0;
	kData[19] = 0.0;
	kData[20] = 0.0;
	kData[21] = 0.0;
	kData[22] = 0.0;
	kData[23] = 0.0;
	kData[24] = 0.0;

	sData[0] = 0.0;
	sData[1] = 0.0;
	sData[2] = 0.0;
	sData[3] = 0.0;
	sData[4] = 0.0;

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

		// h ~ section height / 2 in Linear case
		double maxLoc(fiberLocs[1] - yBarZero), minLoc(fiberLocs[1] - yBarZero);
		for (int i = 0; i < numFibers; i++) {
			if (fiberLocs[i] - yBarZero > maxLoc) { maxLoc = fiberLocs[i] - yBarZero;}
			if (fiberLocs[i] - yBarZero < minLoc) { minLoc = fiberLocs[i] - yBarZero;}
		}
		double h (maxLoc);
		double omega      = y*y*y/(h*h*h) - 0.6*y/h;
		double omegaprime = 3*y*y/(h*h*h) - 0.6/h;

		kData[0] += d00;
		kData[1] += -y * d00;
		kData[2] += d01;
		kData[3] += omegaprime * d01;
		kData[4] += omega * d00;

		kData[5] += -y * d00;
		kData[6] += y * y * d00;
		kData[7] += -y * d01;
		kData[8] += -y * omegaprime * d01;
		kData[9] += -y * omega * d00;

		kData[10] += d10;
		kData[11] += -y * d10;
		kData[12] += d11;
		kData[13] += omegaprime * d11;
		kData[14] += omega * d10;

		kData[15] += omegaprime * d10;
		kData[16] += -y * omegaprime * d10;
		kData[17] += omegaprime * d11;
		kData[18] += omegaprime * omegaprime * d11;
		kData[19] += omega * omegaprime * d10;

		kData[20] += omega * d00;
		kData[21] += -y * omega * d00;
		kData[22] += omega * d01;
		kData[23] += omegaprime * omega * d01;
		kData[24] += omega * omega * d00;

		double fs0 = stress(0)*A;
		double fs1 = stress(1)*A;

		sData[0] += fs0;
		sData[1] += fs0*-y;
		sData[2] += fs1;
		sData[3] += fs1*omegaprime;
		sData[4] += fs0*omega;
	}

	if (alpha != 1.0) {
		double rootAlpha = sqrt(alpha);

		sData[2] *= rootAlpha;
		sData[3] *= rootAlpha;

		kData[2] *= rootAlpha;
		kData[3] *= rootAlpha;

		kData[7] *= rootAlpha;
		kData[8] *= rootAlpha;

		kData[10] *= rootAlpha;
		kData[11] *= rootAlpha;
		kData[14] *= rootAlpha;

		kData[15] *= rootAlpha;
		kData[16] *= rootAlpha;
		kData[19] *= rootAlpha;

		kData[22] *= rootAlpha;
		kData[23] *= rootAlpha;

		kData[12] *= alpha;
		kData[13] *= alpha;

		kData[17] *= alpha;
		kData[18] *= alpha;
	}

	return err;
}

int
NDFiberSectionWarping2d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  //     size 3 so no conflict with matData below if just 1 fiber
  static ID data(3);
  data(0) = this->getTag();
  data(1) = numFibers;
  int dbTag = this->getDbTag();
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "NDFiberSectionWarping2d::sendSelf - failed to send ID data\n";
    return res;
  }    

  static Vector ddata(4+1); // +1 so no conflict when 2 fibers
  ddata(0) = yBar;
  ddata(1) = alpha;
  ddata(2) = yBarZero;
  ddata(3) = DeltaYbar;
  res += theChannel.sendVector(dbTag, commitTag, ddata);
  if (res < 0) {
    opserr <<  "NDFiberSectionWarping2d::sendSelf - failed to send Vector data\n";
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
      opserr <<  "NDFiberSectionWarping2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 2*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "NDFiberSectionWarping2d::sendSelf - failed to send material data\n";
      return res;
    }    

    // now invoke send(0 on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);
    
  }
  
  return res;
}

int
NDFiberSectionWarping2d::recvSelf(int commitTag, Channel &theChannel,
				  FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static ID data(3);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);
  if (res < 0) {
    opserr <<  "NDFiberSectionWarping2d::recvSelf - failed to recv ID data\n";
    return res;
  }    
  this->setTag(data(0));

  static Vector ddata(5);
  res += theChannel.recvVector(dbTag, commitTag, ddata);
  if (res < 0) {
    opserr <<  "NDFiberSectionWarping2d::recvSelf - failed to recv Vector data\n";
    return res;
  }

  yBar = ddata(0);
  alpha = ddata(1);
  yBarZero = ddata(2);
  DeltaYbar = ddata(3);
  
  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
      opserr <<  "NDFiberSectionWarping2d::recvSelf - failed to recv material data\n";
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
	  opserr <<"NDFiberSectionWarping2d::recvSelf -- failed to allocate Material pointers\n";
	  exit(-1);
	}
	
	for (int j=0; j<numFibers; j++)
	  theMaterials[j] = 0;
	
	matData = new double [numFibers*2];
	
	if (matData == 0) {
	  opserr <<"NDFiberSectionWarping2d::recvSelf  -- failed to allocate double array for material data\n";
	  exit(-1);
	}
      }
    }
    
    Vector fiberData(matData, 2*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
      opserr <<  "NDFiberSectionWarping2d::recvSelf - failed to recv material data\n";
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
	opserr <<"NDFiberSectionWarping2d::recvSelf -- failed to allocate double array for material data\n";
	exit(-1);
      }
      
      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }
    
    double Qz = 0.0;
    double A  = 0.0;
    double yLoc, Area;
    
    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = matData[2*i];
      Area = matData[2*i+1];
      A  += Area;
      Qz += yLoc*Area;
    }
    
    yBar = Qz/A;
  }    
  
  return res;
}

void
	NDFiberSectionWarping2d::Print(OPS_Stream &s, int flag)
{
	s << "\nNDFiberSectionWarping2d, tag: " << this->getTag() << endln;
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
	NDFiberSectionWarping2d::setResponse(const char **argv, int argc,
	OPS_Stream &output)
{
	Response *theResponse =0;

	if (argc > 2 && strcmp(argv[0],"fiber") == 0) {

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
					closestDist = dy*dy;
					key = j;
					break;
				}
			}
			// Search the remaining fibers
			for ( ; j < numFibers; j++) {
				if (matTag == theMaterials[j]->getTag()) {
					ySearch = matData[2*j];
					dy = ySearch-yCoord;
					distance = dy*dy;
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
			closestDist = dy*dy;
			key = 0;
			for (int j = 1; j < numFibers; j++) {
				ySearch = matData[2*j];
				dy = ySearch-yCoord;

				distance = dy*dy;
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

	}

	if (theResponse == 0)
	  return SectionForceDeformation::setResponse(argv, argc, output);

	return theResponse;
}


int 
NDFiberSectionWarping2d::getResponse(int responseID, Information &sectInfo)
{
	// Just call the base class method ... don't need to define
	// this function, but keeping it here just for clarity
	return SectionForceDeformation::getResponse(responseID, sectInfo);
}



// AddingSensitivity:BEGIN ////////////////////////////////////
int
	NDFiberSectionWarping2d::setParameter(const char **argv, int argc, Parameter &param)
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
	NDFiberSectionWarping2d::updateParameter(int paramID, Information &info)
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
	NDFiberSectionWarping2d::activateParameter(int paramID)
{
	parameterID = paramID;

	return 0;
}

const Vector &
	NDFiberSectionWarping2d::getSectionDeformationSensitivity(int gradIndex)
{
	return dedh;
}

const Vector &
	NDFiberSectionWarping2d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
	static Vector ds(5);

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

		dsigdh = theMaterials[i]->getStressSensitivity(gradIndex,true);

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
	NDFiberSectionWarping2d::getInitialTangentSensitivity(int gradIndex)
{
	static Matrix dksdh(5,5);

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
	NDFiberSectionWarping2d::commitSensitivity(const Vector& defSens,
	int gradIndex, int numGrads)
{
	double d0 = defSens(0);
	double d1 = defSens(1);
	double d2 = defSens(2);
	double d3 = defSens(3);
	double d4 = defSens(4);

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
	double kappa    = e(1);
	double gamma    = e(2);
	double phi      = e(3);
	double phiprime = e(4);

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
