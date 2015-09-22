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
// Revision: 
//
// Description: This file contains the implementation for the
// NDFiber3d class. NDFiber3d provides the abstraction of a
// uniaXial Fiber whose position  is defined with only one coordinate.
// The NDFiber3d is subjected to a stress state with 
// only one nonzero axial stress and corresponding axial strain.
//
// What: "@(#) NDFiber3d.h, revA"

#include <stdlib.h>

#include <NDMaterial.h>
#include <NDFiber3d.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
#include <SectionForceDeformation.h>
#include <Information.h>
#include <Parameter.h>
#include <FiberResponse.h>
#include <elementAPI.h>

Matrix NDFiber3d::ks(6,6); 
Vector NDFiber3d::fs(6); 
ID NDFiber3d::code(6);

static int numNDFiber3d = 0;

void* OPS_NDFiber3d()
{
    if(OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"insufficient arguments for NDFiber3d\n";
	return 0;
    }

    // get data
    int numData = 3;
    double data[3];
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) return 0;

    // get mat tag
    int tag;
    numData = 1;
    if(OPS_GetIntInput(&numData,&tag) < 0) return 0;

    // get material
    NDMaterial* theMat = OPS_getNDMaterial(tag);
    if(theMat == 0) {
	opserr<<"invalid NDMaterial tag\n";
	return 0;
    }

    return new NDFiber3d(numNDFiber3d++,*theMat,data[2],data[0],data[1]);
}


// constructor:
NDFiber3d::NDFiber3d(int tag, NDMaterial &theMat,
		     double Area, double yy, double zz):
  Fiber(tag, FIBER_TAG_ND3d),
  theMaterial(0), area(Area), y(yy), z(zz)
{
  theMaterial = theMat.getCopy("BeamFiber");
  
  if (theMaterial == 0) {
    opserr << "NDFiber3d::NDFiber3d -- failed to get copy of NDMaterial\n";
    exit(-1);
  }
  
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_VY;
    code(4) = SECTION_RESPONSE_VZ;
    code(5) = SECTION_RESPONSE_T;
  }
}

// constructor for blank object that recvSelf needs to be invoked upon
NDFiber3d::NDFiber3d(): 
  Fiber(0, FIBER_TAG_ND3d),
  theMaterial(0), area(0), y(0.0), z(0.0)
{
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_VY;
    code(4) = SECTION_RESPONSE_VZ;
    code(5) = SECTION_RESPONSE_T;
  }
}


// Destructor: 
NDFiber3d::~NDFiber3d ()
{
  if (theMaterial != 0)
    delete theMaterial;
}


int   
NDFiber3d::setTrialFiberStrain(const Vector &vs)
{
  static Vector strain(3);
  strain(0) = 0;
  strain(1) = 0;
  strain(2) = 0;
  
  opserr << "NDFiber3d::setTrialFiberStrain() -- not implemented" << endln;

  return theMaterial->setTrialStrain(strain);
}



// get fiber stress resultants 
Vector &
NDFiber3d::getFiberStressResultants (void) 
{
  fs.Zero();
  
  opserr << "NDFiber3d::getFiberStressResultants() -- not implemented" << endln;

  return fs;
}



// get contribution of fiber to section tangent stiffness
Matrix &
NDFiber3d::getFiberTangentStiffContr(void) 
{
  ks.Zero();

  opserr << "NDFiber3d::getFiberTangentStiffContr() -- not implemented" << endln;

  return ks;
}

Fiber*
NDFiber3d::getCopy (void)
{
   // make a copy of the fiber 
  NDFiber3d *theCopy = new NDFiber3d (this->getTag(), 
				      *theMaterial, area, y, z);

  return theCopy;
}  

int
NDFiber3d::getOrder(void)
{
  return 6;
}

const ID&
NDFiber3d::getType(void)
{
  return code;
}

int   
NDFiber3d::commitState(void)
{
  return theMaterial->commitState();
}


int   
NDFiber3d::revertToLastCommit(void)
{
  return theMaterial->revertToLastCommit();
}


int   
NDFiber3d::revertToStart(void)
{
  return theMaterial->revertToStart();
}


int   
NDFiber3d::sendSelf(int commitTag, Channel &theChannel)
{
  // 
  // store tag and material info in an ID and send it
  //
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID idData(3);
  
  idData(0) = this->getTag();
  idData(1) = theMaterial->getClassTag();
  
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  
  idData(2) = matDbTag;
  
  res += theChannel.sendID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "NDFiber3d::sendSelf - failed to send ID data\n";
    return res;
  }    
  
  // 
  // store area and position data in a vector and send it
  //
  static Vector dData(3);
  
  dData(0) = area;
  dData(1) = y;
  dData(2) = z;
  
  res += theChannel.sendVector(dbTag, commitTag, dData);
  if (res < 0) {
    opserr << "NDFiber3d::sendSelf - failed to send Vector data\n";
    return res;
  }    

  // now invoke sendSelf on the material
  res += theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "NDFiber3d::sendSelf - failed to send UniaxialMaterial\n";
      return res;
  }
    
  return res;
}


int   
NDFiber3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // 
  // get tag and material info from an ID
  //
  
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID idData(3);
    
  res += theChannel.recvID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "NDFiber3d::recvSelf - failed to receive ID data\n";
    return res;
  }    
  
  this->setTag(idData(0));

  // 
  // get area from a vector received from channel
  //
  
  static Vector dData(3);
  
  res += theChannel.recvVector(dbTag, commitTag, dData);
  if (res < 0) {
      opserr << "NDFiber3d::recvSelf - failed to receive Vector data\n";
      return res;
  }
  
  area = dData(0);
  y = dData(1);
  z = dData(2);

  //
  // now we do the material stuff
  //
  
  int matClassTag = idData(1);    
  
    // if we have a material, check it is of correct type
  if (theMaterial != 0) {
    if (matClassTag != theMaterial->getClassTag()) {
      delete theMaterial;
      theMaterial = 0;
    }
    }
  
  // if no material we need to get one,
  // NOTE: not an else if in case deleted in if above
  if (theMaterial == 0) {
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "NDFiber3d::recvSelf() - " <<
	  "failed to get a NDMaterial of type " << matClassTag << endln;
      return -1;
    }
  }
  
    // set the materials dbTag and invoke recvSelf on the material
  theMaterial->setDbTag(idData(2));
  
  // now invoke recvSelf on the material
  res += theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "NDFiber3d::recvSelf() - the material failed in recvSelf()\n";
    return res;
  }    	
  
  return res;
}


void NDFiber3d::Print(OPS_Stream &s, int flag)
{
  s << "\nNDFiber3d, tag: " << this->getTag() << endln;
  s << "\tArea: " << area << endln; 
  s << "\tLocation (y,z): " << y << " " << z << endln; 
  s << "\tMaterial, tag: " << theMaterial->getTag() << endln;
}

Response*
NDFiber3d::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  if (argc == 0)
    return 0;
  
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
    return new FiberResponse(this, 1, Vector(2));
  
  else
    return theMaterial->setResponse(argv, argc, s);
}

int
NDFiber3d::getResponse(int responseID, Information &fibInfo)
{
  switch(responseID) {
  case 1:
    return fibInfo.setVector(this->getFiberStressResultants());
    
  default:
    return -1;
  }
}

void 
NDFiber3d::getFiberLocation(double &yLoc, double &zLoc)
{
  yLoc = y;
  zLoc = z;
}

int
NDFiber3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"A") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"y") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"z") == 0)
    return param.addObject(3, this);

  else
    return theMaterial->setParameter(argv, argc, param);
}

int
NDFiber3d::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    area = info.theDouble;
    return 0;
  case 2:
    y = info.theDouble;
    return 0;
  case 3:
    z = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
NDFiber3d::activateParameter(int parameterID)
{
  return -1;
}

const Vector&
NDFiber3d::getFiberSensitivity(int gradNumber, bool cond)
{
  return Fiber::getFiberSensitivity(gradNumber, cond);
}

int 
NDFiber3d::commitSensitivity(const Vector &dedh, int gradNumber,
			     int numGrads)
{
  return -1;
}

