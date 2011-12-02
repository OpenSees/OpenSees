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
                                                                        
// $Revision: 1.21 $
// $Date: 2004-08-25 22:17:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/FiberSection2d.cpp,v $
                                                                        
// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSection2d.

#include <stdlib.h>

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

ID FiberSection2d::code(2);

// constructors:
FiberSection2d::FiberSection2d(int tag, int num, Fiber **fibers): 
  SectionForceDeformation(tag, SEC_TAG_FiberSection2d),
  numFibers(num), theMaterials(0), matData(0),
  yBar(0.0), e(2), eCommit(2), s(0), ks(0)
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


    double Qz = 0.0;
    double A  = 0.0;
    
    for (int i = 0; i < numFibers; i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();
      A  += Area;
      Qz += yLoc*Area;
      matData[i*2] = -yLoc;
      matData[i*2+1] = Area;
      UniaxialMaterial *theMat = theFiber->getMaterial();
      theMaterials[i] = theMat->getCopy();

      if (theMaterials[i] == 0) {
	opserr << "FiberSection2d::FiberSection2d -- failed to get copy of a Material\n";
	exit(-1);
      }
    }    

    yBar = -Qz/A;  
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

// AddingSensitivity:BEGIN ////////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////

}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSection2d::FiberSection2d():
  SectionForceDeformation(0, SEC_TAG_FiberSection2d),
  numFibers(0), theMaterials(0), matData(0),
  yBar(0.0), e(2), eCommit(2), s(0), ks(0)
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

// AddingSensitivity:BEGIN ////////////////////////////////////
	parameterID = 0;
// AddingSensitivity:END //////////////////////////////////////
}

int
FiberSection2d::addFiber(Fiber &newFiber)
{
  // need to create larger arrays
  int newSize = numFibers+1;
  UniaxialMaterial **newArray = new UniaxialMaterial *[newSize]; 
  double *newMatData = new double [2 * newSize];
  if (newArray == 0 || newMatData == 0) {
    opserr <<"FiberSection2d::addFiber -- failed to allocate Fiber pointers\n";
    return -1;
  }

  // copy the old pointers and data
  int i;
  for (i = 0; i < numFibers; i++) {
    newArray[i] = theMaterials[i];
    newMatData[2*i] = matData[2*i];
    newMatData[2*i+1] = matData[2*i+1];
  }

  // set the new pointers and data
  double yLoc, zLoc, Area;
  newFiber.getFiberLocation(yLoc, zLoc);
  Area = newFiber.getArea();
  newMatData[numFibers*2] = -yLoc;
  newMatData[numFibers*2+1] = Area;
  UniaxialMaterial *theMat = newFiber.getMaterial();
  newArray[numFibers] = theMat->getCopy();

  if (newArray[numFibers] == 0) {
    opserr <<"FiberSection2d::addFiber -- failed to get copy of a Material\n";
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
  double A  = 0.0;

  // Recompute centroid
  for (i = 0; i < numFibers; i++) {
    yLoc = -matData[2*i];
    Area = matData[2*i+1];
    A  += Area;
    Qz += yLoc*Area;
  }

  yBar = -Qz/A;

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
}

int FiberSection2d::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;

  e = deforms;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;

  int loc = 0;

  double d0 = deforms(0);
  double d1 = deforms(1);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double A = matData[loc++];

    // determine material strain and set it
    double strain = d0 + y*d1;
    double tangent, stress;
    res = theMat->setTrial(strain, stress, tangent);

    double ks0 = tangent * A;
    double ks1 = ks0 * y;
    kData[0] += ks0;
    kData[1] += ks1;
    kData[3] += ks1 * y;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
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
  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;

  int loc = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double A = matData[loc++];

    double tangent = theMat->getInitialTangent();

    double ks0 = tangent * A;
    double ks1 = ks0 * y;
    kData[0] += ks0;
    kData[1] += ks1;
    kData[3] += ks1 * y;
  }

  kData[2] = kData[1];

  return *ks;
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

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->yBar = yBar;

  theCopy->kData[0] = kData[0];
  theCopy->kData[1] = kData[1];
  theCopy->kData[2] = kData[2];
  theCopy->kData[3] = kData[3];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];

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

  eCommit = e;

  return err;
}

int
FiberSection2d::revertToLastCommit(void)
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;


  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;
  
  int loc = 0;
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double A = matData[loc++];

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    // get material stress & tangent for this strain and determine ks and fs
    double tangent = theMat->getTangent();
    double stress = theMat->getStress();
    double ks0 = tangent * A;
    double ks1 = ks0 * y;
    kData[0] += ks0;
    kData[1] += ks1;
    kData[3] += ks1 * y;

    double fs0 = stress * A;
    sData[0] = fs0;
    sData[1] = fs0 * y;
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
  
  int loc = 0;
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++] - yBar;
    double A = matData[loc++];

    // invoke revertToLast on the material
    err += theMat->revertToStart();

    // get material stress & tangent for this strain and determine ks and fs
    double tangent = theMat->getTangent();
    double stress = theMat->getStress();
    double ks0 = tangent * A;
    double ks1 = ks0 * y;
    kData[0] += ks0;
    kData[1] += ks1;
    kData[3] += ks1 * y;

    double fs0 = stress * A;
    sData[0] = fs0;
    sData[1] = fs0 * y;
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

    double Qz = 0.0;
    double A  = 0.0;
    double yLoc, Area;

    // Recompute centroid
    for (i = 0; i < numFibers; i++) {
      yLoc = -matData[2*i];
      Area = matData[2*i+1];
      A  += Area;
      Qz += yLoc*Area;
    }
    
    yBar = -Qz/A;
  }    

  return res;
}

void
FiberSection2d::Print(OPS_Stream &s, int flag)
{
  s << "\nFiberSection2d, tag: " << this->getTag() << endln;
  s << "\tSection code: " << code;
  s << "\tNumber of Fibers: " << numFibers << endln;
  s << "\tCentroid: " << -yBar << endln;

  if (flag == 1) {
    for (int i = 0; i < numFibers; i++) {
      s << "\nLocation (y) = (" << -matData[2*i] << ")";
      s << "\nArea = " << matData[2*i+1] << endln;
      theMaterials[i]->Print(s, flag);
    }
  }
}

Response*
FiberSection2d::setResponse(const char **argv, int argc, Information &sectInfo)
{
  // See if the response is one of the defaults
  Response *res = SectionForceDeformation::setResponse(argv, argc, sectInfo);
  if (res != 0)
    return res;
  
  // Check if fiber response is requested
  else if (strcmp(argv[0],"fiber") == 0) {
    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 2)          // not enough data input
      return 0;

    if (argc <= 3) {		  // fiber number was input directly
      key = atoi(argv[1]);
      if (key < numFibers && key >= 0)
	return theMaterials[key]->setResponse(&argv[passarg], argc-passarg, sectInfo);
      else 
	return 0;
    }

    if (argc > 4) {  // find fiber closest to coord. with mat tag
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double closestDist;
      double ySearch, dy;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < numFibers; j++) {
	if (matTag == theMaterials[j]->getTag()) {
	  ySearch = -matData[2*j];
	  dy = ySearch-yCoord;
	  closestDist = fabs(dy);
	  key = j;
	  break;
	}
      }
      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
	if (matTag == theMaterials[j]->getTag()) {
	  ySearch = -matData[2*j];
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
      ySearch = -matData[0];
      dy = ySearch-yCoord;
      closestDist = fabs(dy);
      key = 0;
      for (int j = 1; j < numFibers; j++) {
	ySearch = -matData[2*j];
	dy = ySearch-yCoord;
	distance = fabs(dy);
	if (distance < closestDist) {
	  closestDist = distance;
	  key = j;
	}
      }
      passarg = 3;
    }
    
    if (key < numFibers)
      return theMaterials[key]->setResponse(&argv[passarg],argc-passarg,sectInfo);
    else
      return 0;
  }
  
  // otherwise response quantity is unknown for the FiberSection class
  else
    return 0;
}


int 
FiberSection2d::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}



// AddingSensitivity:BEGIN ////////////////////////////////////
int
FiberSection2d::setParameter (const char **argv, int argc, Information &info)
{
	// Initial declarations
	int parameterID = -1;

	// Check if the parameter belongs to the material (only option for now)
	if (strcmp(argv[0],"-material") == 0 || strcmp(argv[0],"material") == 0) {

		// Get the tag of the material
		int materialTag = atoi(argv[1]);

		// Loop over fibers to find the right material
		for (int i=0; i<numFibers; i++) {
			if (materialTag == theMaterials[i]->getTag()) {
				parameterID = theMaterials[i]->setParameter(&argv[2], argc-2, info);
			}
		}

		// Check that the parameterID is valid
		if (parameterID < 0) {
			opserr << "FiberSection2d::setParameter() - could not set parameter. " << endln;
			return -1;
		}
		
		else {
			return (parameterID + 100000*this->getTag() + 1000*materialTag);
		}
	} 

	else {
		opserr << "FiberSection2d::setParameter() - could not set parameter. " << endln;
		return -1;
	}
}

int
FiberSection2d::updateParameter (int parameterID, Information &info)
{

	// Check if it is a material parameter (only option for now)
	if (parameterID > 1000) {

		// Get section number and material number
		int sectionTag = this->getTag();
		parameterID = parameterID - sectionTag*100000;
		int materialTag = (int)( floor((double)parameterID) / (1000) );
		parameterID = parameterID - materialTag*1000;

		int ok = -1;
		for (int i=0; i<numFibers; i++) {
			if (materialTag == theMaterials[i]->getTag()) {
				ok = theMaterials[i]->updateParameter(parameterID, info);
			}
		}

		if (ok < 0) {
			opserr << "FiberSection2d::updateParameter() - could not update parameter. " << endln;
			return ok;
		}

		else {
			return ok;
		}
	}
	else {
		opserr << "FiberSection2d::updateParameter() - could not update parameter. " << endln;
		return -1;
	}
}

int
FiberSection2d::activateParameter(int passedParameterID)
{
	// Note that the parameteID that is stored here at the 
	// section level contains all information about section
	// and material tag number:
	parameterID = passedParameterID;

	if (passedParameterID == 0 ) {

		// "Zero out" the identifier in all materials
		for (int i=0; i<numFibers; i++) {
			theMaterials[i]->activateParameter(passedParameterID);
		}
	}

	else {

		// Extract section and material tags
		int activeSectionTag = (int)( floor((double)passedParameterID) / (100000) );
		passedParameterID -= activeSectionTag*100000;
		int activeMaterialTag = (int)( floor((double)passedParameterID) / (1000) );
		passedParameterID -= activeMaterialTag*1000;

		// Go down to the sections and set appropriate flags
		for (int i=0; i<numFibers; i++) {
			if (activeMaterialTag == theMaterials[i]->getTag()) {
				theMaterials[i]->activateParameter(passedParameterID);
			}
		}
	}

	return 0;
}

const Vector &
FiberSection2d::getSectionDeformationSensitivity(int gradNumber)
{
	static Vector dummy(2);
	return dummy;
}

const Vector &
FiberSection2d::getStressResultantSensitivity(int gradNumber, bool conditional)
{
	static Vector ds(2);

	ds.Zero();

	double y, A, stressGradient;
	int loc = 0;

	for (int i = 0; i < numFibers; i++) {
		y = matData[loc++];
		A = matData[loc++];

		stressGradient = theMaterials[i]->getStressSensitivity(gradNumber,true);
		stressGradient = stressGradient * A;
		ds(0) += stressGradient;
		ds(1) += stressGradient * y;
	}

	return ds;
}

const Matrix &
FiberSection2d::getSectionTangentSensitivity(int gradNumber)
{
	static Matrix something(2,2);

	something.Zero();

	return something;
}

int
FiberSection2d::commitSensitivity(const Vector& defSens, int gradNumber, int numGrads)
{

  int loc = 0;

  double d0 = defSens(0);
  double d1 = defSens(1);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = matData[loc++];
    loc++;

    // determine material strain and set it
    double strainSens = d0 + y*d1;
	theMat->commitSensitivity(strainSens,gradNumber,numGrads);
  }

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////
