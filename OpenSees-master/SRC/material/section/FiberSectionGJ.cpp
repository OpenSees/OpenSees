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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSectionGJ.cpp,v $
                                                                        
// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSection2d.

#include <stdlib.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <Fiber.h>
#include <classTags.h>
#include <FiberSectionGJ.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>

#include <string.h>

ID FiberSectionGJ::code(4);
Vector FiberSectionGJ::s(4);
Matrix FiberSectionGJ::ks(4,4);

// constructors:
FiberSectionGJ::FiberSectionGJ(int tag, int num, Fiber **fibers, double gj): 
  SectionForceDeformation(tag, SEC_TAG_FiberSectionGJ),
  numFibers(num), theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(4), GJ(gj)
{
  if (numFibers != 0) {
    theMaterials = new UniaxialMaterial *[numFibers];

    if (theMaterials == 0) {
      opserr << "FiberSectionGJ::FiberSectionGJ -- failed to allocate Material pointers\n";
      exit(-1);
    }

    matData = new double [numFibers*3];

    if (matData == 0) {
      opserr << "FiberSectionGJ::FiberSectionGJ -- failed to allocate double array for material data\n";
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

      Qz += yLoc*Area;
      Qy += zLoc*Area;
      A  += Area;

      matData[i*3] = -yLoc;
      matData[i*3+1] = zLoc;
      matData[i*3+2] = Area;
      UniaxialMaterial *theMat = theFiber->getMaterial();
      theMaterials[i] = theMat->getCopy();

      if (theMaterials[i] == 0) {
	opserr << "FiberSectionGJ::FiberSectionGJ -- failed to get copy of a Material\n";
	exit(-1);
      }
    }

    yBar = -Qz/A;
    zBar = Qy/A;
  }

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  for (int i=0; i<6; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSectionGJ::FiberSectionGJ():
  SectionForceDeformation(0, SEC_TAG_FiberSectionGJ),
  numFibers(0), theMaterials(0), matData(0),
  yBar(0.0), zBar(0.0), e(4), GJ(1.0)
{
  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  for (int i=0; i<6; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
}

int
FiberSectionGJ::addFiber(Fiber &newFiber)
{
  // need to create a larger array
  int newSize = numFibers+1;

  UniaxialMaterial **newArray = new UniaxialMaterial *[newSize]; 
  double *newMatData = new double [3 * newSize];
  
  if (newArray == 0 || newMatData == 0) {
    opserr << "FiberSectionGJ::addFiber -- failed to allocate Fiber pointers\n";
    return -1;
  }

  // copy the old pointers
  int i;
  for (i = 0; i < numFibers; i++) {
    newArray[i] = theMaterials[i];
    newMatData[3*i] = matData[3*i];
    newMatData[3*i+1] = matData[3*i+1];
    newMatData[3*i+2] = matData[3*i+2];
  }
  // set the new pointers
  double yLoc, zLoc, Area;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  newMatData[numFibers*3] = -yLoc;
  newMatData[numFibers*3+1] = zLoc;
  newMatData[numFibers*3+2] = Area;
  UniaxialMaterial *theMat = newFiber.getMaterial();
  newArray[numFibers] = theMat->getCopy();

  if (newArray[numFibers] == 0) {
    opserr << "FiberSectionGJ::addFiber -- failed to get copy of a Material\n";
			  

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
    yLoc = -matData[3*i];
    zLoc = matData[3*i+1];
    Area = matData[3*i+2];
    A  += Area;
    Qz += yLoc*Area;
    Qy += zLoc*Area;
  }

  yBar = -Qz/A;
  zBar = Qy/A;

  return 0;
}



// destructor:
FiberSectionGJ::~FiberSectionGJ()
{
  if (theMaterials != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theMaterials[i] != 0)
	delete theMaterials[i];
      
    delete [] theMaterials;
  }

  if (matData != 0)
    delete [] matData;
}

int
FiberSectionGJ::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0; 

  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0; 

  int loc = 0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double z = matData[loc++] - zBar;
    double A = matData[loc++];

    // determine material strain and set it
    double strain = d0 + y*d1 + z*d2;
    double tangent, stress;
    res = theMat->setTrial(strain, stress, tangent);

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;
    
    kData[3] += vas1 * y;
    kData[4] += vas1as2;
    
    kData[5] += vas2 * z; 

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  return res;
}

const Matrix&
FiberSectionGJ::getInitialTangent(void)
{
  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;

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

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;
    
    kData[3] += vas1 * y;
    kData[4] += vas1as2;
    
    kData[5] += vas2 * z; 
  }

  ks(0,0) = kData[0];
  ks(0,1) = ks(1,0) = kData[1];
  ks(0,2) = ks(2,0) = kData[2];
  ks(1,1) = kData[3];
  ks(1,2) = ks(2,1) = kData[4];
  ks(2,2) = kData[5];

  ks(3,3) = GJ;

  return ks;
}

const Vector&
FiberSectionGJ::getSectionDeformation(void)
{
  return e;
}

const Matrix&
FiberSectionGJ::getSectionTangent(void)
{
  ks(0,0) = kData[0];
  ks(0,1) = ks(1,0) = kData[1];
  ks(0,2) = ks(2,0) = kData[2];
  ks(1,1) = kData[3];
  ks(1,2) = ks(2,1) = kData[4];
  ks(2,2) = kData[5];

  ks(3,3) = GJ;

  return ks;
}

const Vector&
FiberSectionGJ::getStressResultant(void)
{
  s(0) = sData[0];
  s(1) = sData[1];
  s(2) = sData[2];

  s(3) = GJ*e(3);

  return s;
}

SectionForceDeformation*
FiberSectionGJ::getCopy(void)
{
  FiberSectionGJ *theCopy = new FiberSectionGJ();
  theCopy->setTag(this->getTag());

  theCopy->numFibers = numFibers;

  if (numFibers != 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    if (theCopy->theMaterials == 0) {
      opserr << "FiberSectionGJ::FiberSectionGJ -- failed to allocate Material pointers\n";
      exit(-1);
    }

    theCopy->matData = new double [numFibers*3];

    if (theCopy->matData == 0) {
      opserr << "FiberSectionGJ::FiberSectionGJ -- failed to allocate double array for material data\n";
      exit(-1);
    }    
    for (int i = 0; i < numFibers; i++) {
      theCopy->matData[i*3] = matData[i*3];
      theCopy->matData[i*3+1] = matData[i*3+1];
      theCopy->matData[i*3+2] = matData[i*3+2];
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == 0) {
	opserr << "FiberSectionGJ::getCopy -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    
  }

  theCopy->e = e;
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;

  for (int i=0; i<6; i++)
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];

  theCopy->GJ = GJ;

  return theCopy;
}

const ID&
FiberSectionGJ::getType ()
{
  return code;
}

int
FiberSectionGJ::getOrder () const
{
  return 4;
}

int
FiberSectionGJ::commitState(void)
{
  int err = 0;

  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  return err;
}

int
FiberSectionGJ::revertToLastCommit(void)
{
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;

  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0; 

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
    
    kData[3] += vas1 * y;
    kData[4] += vas1as2;
    
    kData[5] += vas2 * z; 

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  return err;
}

int
FiberSectionGJ::revertToStart(void)
{
  // revert the fibers to start    
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0;
  kData[3] = 0.0; kData[4] = 0.0; kData[5] = 0.0;

  sData[0] = 0.0; sData[1] = 0.0; sData[2] = 0.0; 

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
    
    kData[3] += vas1 * y;
    kData[4] += vas1as2;
    
    kData[5] += vas2 * z; 

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  return err;
}

int
FiberSectionGJ::sendSelf(int commitTag, Channel &theChannel)
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
FiberSectionGJ::recvSelf(int commitTag, Channel &theChannel,
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
FiberSectionGJ::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "\nFiberSectionGJ, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: (" << -yBar << ", " << zBar << ')' << endln;
    s << "\tTorsional Stiffness: " << GJ << endln;
    
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
      int loc = 0;
      for (int i = 0; i < numFibers; i++) {
	s << "\nLocation (y, z) = (" << -matData[loc] << ", " << matData[loc+1] << ")";
	s << "\nArea = " << matData[loc+2] << endln;
	theMaterials[i]->Print(s, flag);
	loc+= 3;
      }
    }
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << "\t\t\t{";
      s << "\"name\": \"" << this->getTag() << "\", ";
      s << "\"type\": \"FiberSectionGJ\", ";
      s << "\"GJ\": " << GJ << ", ";
      s << "\"fibers\": [\n";
      for (int i = 0; i < numFibers; i++) {
          s << "\t\t\t\t{\"coord\": [" << matData[3 * i] << ", " << matData[3 * i + 1] << "], ";
          s << "\"area\": " << matData[3 * i + 2] << ", ";
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
FiberSectionGJ::setResponse(const char **argv, int argc, OPS_Stream &output)
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
    double closestDist = 0.0;
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
FiberSectionGJ::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}

int
FiberSectionGJ::setParameter (const char **argv, int argc, Parameter &param)
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
