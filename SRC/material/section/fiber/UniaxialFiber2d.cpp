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
                                                                        
// $Revision: 1.11 $
// $Date: 2007-02-02 01:18:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/UniaxialFiber2d.cpp,v $
                                                                        
                                                                        
// Written: Remo Magalhaes de Souza
// Created: 10/98
// Revision: 
//
// Description: This file contains the implementation for the
// UniaxialFiber2d class. UniaxialFiber2d provides the abstraction of a
// uniaXial Fiber whose position  is defined with only one coordinate.
// The UniaxialFiber2d is subjected to a stress state with 
// only one nonzero axial stress and corresponding axial strain.
//
// What: "@(#) UniaxialFiber2d.h, revA"

#include <stdlib.h>

#include <UniaxialMaterial.h>
#include <UniaxialFiber2d.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
#include <SectionForceDeformation.h>
#include <Information.h>
#include <Parameter.h>
#include <FiberResponse.h>
#include <elementAPI.h>

Matrix UniaxialFiber2d::ks(2,2); 
Vector UniaxialFiber2d::fs(2); 
ID UniaxialFiber2d::code(2);

static int numUniaxialFiber2d = 0;

void* OPS_UniaxialFiber2d()
{
    if(OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"insufficient arguments for UniaxialFiber2d\n";
	return 0;
    }

    // get data
    int numData = 3;
    double data[3];
    if(OPS_GetDoubleInput(&numData,&data[0]) < 0) {
	opserr<<"WARNING failed to read double\n";
	return 0;
    }

    // get mat tag
    int tag;
    numData = 1;
    if(OPS_GetIntInput(&numData,&tag) < 0) {
	opserr<<"WARNING failed to read int\n";
	return 0;
    }

    // get material
    UniaxialMaterial* theMat = OPS_getUniaxialMaterial(tag);
    if(theMat == 0) {
	opserr<<"invalid UniaxialMaterial tag\n";
	return 0;
    }

    return new UniaxialFiber2d(numUniaxialFiber2d++,*theMat,data[2],data[0]);
}


// constructor:
UniaxialFiber2d::UniaxialFiber2d(int tag, 
                                 UniaxialMaterial &theMat,
                                 double Area, double position):
                                 Fiber(tag, FIBER_TAG_Uniaxial2d),
                                 theMaterial(0), area(Area), y(-position)
{
  theMaterial = theMat.getCopy();  // get a copy of the MaterialModel
  
  if (theMaterial == 0) {
    opserr <<"UniaxialFiber2d::UniaxialFiber2d  -- failed to get copy of UniaxialMaterial\n";
    exit(-1);
  }
  
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
  }
}

// constructor for blank object that recvSelf needs to be invoked upon
UniaxialFiber2d::UniaxialFiber2d(): Fiber(0, FIBER_TAG_Uniaxial2d),
                                    theMaterial(0), area(0), y(0.0)
{
  if (code(0) != SECTION_RESPONSE_P) {
    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
  }
}


// Destructor: 
UniaxialFiber2d::~UniaxialFiber2d ()
{
  if (theMaterial != 0)
    delete theMaterial;
}


int   
UniaxialFiber2d::setTrialFiberStrain(const Vector &vs)
{
  // Use the section kinematic matrix to get the fiber strain
  // eps = as * vs;
  double strain = vs(0) + y*vs(1); // fiber strain
  
  return theMaterial->setTrialStrain(strain);
}



// get fiber stress resultants 
Vector &
UniaxialFiber2d::getFiberStressResultants (void) 
{
  // Use the section kinematic matrix to get the fiber 
  // stress resultant vector
  // fs = as^ * area * sigma;
  double df = theMaterial->getStress() * area;
  
  fs(0) = df;
  fs(1) = y * df;
  
  return fs;
}



// get contribution of fiber to section tangent stiffness
Matrix &
UniaxialFiber2d::getFiberTangentStiffContr(void) 
{
  // Use the section kinematic matrix to get the fiber 
  // tangent stiffness matrix
  // ks = (as^as) * area * Et;
  double value = theMaterial->getTangent() * area;
  double value_as1 = value*y;
  
    ks(0,0) = value;
    ks(0,1) = value_as1;
    ks(1,0) = value_as1;
    ks(1,1) = value_as1 * y;
    
    return ks;
}

Fiber*
UniaxialFiber2d::getCopy (void)
{
   // make a copy of the fiber 
   UniaxialFiber2d *theCopy = new UniaxialFiber2d (this->getTag(), 
                                                   *theMaterial, area, -y);
   return theCopy;
}  

int
UniaxialFiber2d::getOrder(void)
{
	return 2;
}

const ID&
UniaxialFiber2d::getType(void)
{
	return code;
}

int   
UniaxialFiber2d::commitState(void)
{
   return theMaterial->commitState();
}


int   
UniaxialFiber2d::revertToLastCommit(void)
{
   return theMaterial->revertToLastCommit();
}


int   
UniaxialFiber2d::revertToStart(void)
{
   return theMaterial->revertToStart();
}


int   
UniaxialFiber2d::sendSelf(int commitTag, Channel &theChannel)
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
    opserr << "UniaxialFiber2d::sendSelf - failed to send ID data\n";
    return res;
  }    
  
  // 
  // store area and position data in a vector and send it
  //
  static Vector dData(2);
  
  dData(0) = area;
  dData(1) = y;
  
  res += theChannel.sendVector(dbTag, commitTag, dData);
  if (res < 0) {
    opserr << "UniaxialFiber2d::sendSelf - failed to send Vector data\n";
    return res;
  }    

  // now invoke sendSelf on the material
  res += theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "UniaxialFiber2d::sendSelf - failed to send UniaxialMaterial\n";
      return res;
  }
    
  return res;
}


int   
UniaxialFiber2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  // 
  // get tag and material info from an ID
  //
  
  int res = 0;
  
  int dbTag = this->getDbTag();
  
  static ID idData(3);
    
  res += theChannel.recvID(dbTag, commitTag, idData);
  if (res < 0) {
    opserr << "UniaxialFiber2d::rcvSelf - failed to receive ID data\n";
    return res;
  }    
  
  this->setTag(idData(0));

  // 
  // get area from a vector received from channel
  //
  
  static Vector dData(2);
  
  res += theChannel.recvVector(dbTag, commitTag, dData);
  if (res < 0) {
      opserr << "UniaxialFiber2d::recvSelf - failed to receive Vector data\n";
      return res;
  }
  
  area = dData(0);
  y = dData(1);

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
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "UniaxialFiber2d::recvSelf() - " <<
	  "failed to get a UniaxialMaterial of type " << matClassTag << endln;
      return -1;
    }
  }
  
    // set the materials dbTag and invoke recvSelf on the material
  theMaterial->setDbTag(idData(2));
  
  // now invoke recvSelf on the material
  res += theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "UniaxialFiber2d::recvSelf() - the material failed in recvSelf()\n";
    return res;
  }    	
  
  return res;
}


void UniaxialFiber2d::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
  s << "\nUniaxialFiber2d, tag: " << this->getTag() << endln;
  s << "\tArea: " << area << endln; 
  s << "\tMatrix as: " << 1.0 << " " << y << endln; 
  s << "\tMaterial, tag: " << theMaterial->getTag() << endln;
  }
if (flag == OPS_PRINT_PRINTMODEL_JSON) {
	 s << "\t\t\t\t{\"type\": \"fiber\", \"material\": "<<theMaterial->getTag()<<", ";
	 s << "\"area\": "<<area<<", \"coord\": [" << y << ", 0.0]}";
}
}

Response*
UniaxialFiber2d::setResponse(const char **argv, int argc, OPS_Stream &s)
{
  if (argc == 0)
    return 0;
  
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
    return new FiberResponse(this, 1, Vector(2));
  
  else
    return theMaterial->setResponse(argv, argc, s);
}

int
UniaxialFiber2d::getResponse(int responseID, Information &fibInfo)
{
  switch(responseID) {
  case 1:
    return fibInfo.setVector(this->getFiberStressResultants());
    
  default:
    return -1;
  }
}

void 
UniaxialFiber2d::getFiberLocation(double &yLoc, double &zLoc)
{
  yLoc = -y;
  zLoc = 0.0;
}

int
UniaxialFiber2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"A") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"y") == 0)
    return param.addObject(2, this);

  else
    return theMaterial->setParameter(argv, argc, param);
}

int
UniaxialFiber2d::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    area = info.theDouble;
    return 0;
  case 2:
    y = -info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
UniaxialFiber2d::activateParameter(int parameterID)
{
  return -1;
}

const Vector&
UniaxialFiber2d::getFiberSensitivity(int gradNumber, bool cond)
{
  return Fiber::getFiberSensitivity(gradNumber, cond);
}

int 
UniaxialFiber2d::commitSensitivity(const Vector &dedh, int gradNumber,
				   int numGrads)
{
  return -1;
}

