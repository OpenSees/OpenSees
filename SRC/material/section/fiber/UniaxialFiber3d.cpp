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
                                                                        
// $Revision: 1.2 $
// $Date: 2000-12-18 10:44:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/UniaxialFiber3d.cpp,v $
                                                                        
                                                                        
// File: ~/section/UniaxialFiber3d.C
//
// Written: Remo Magalhaes de Souza
// Created: 10/98
// Revision: 
//
// Description: This file contains the implementation for the
// UniaxialFiber3d class. UniaxialFiber3d provides the abstraction of a
// uniaxial fiber that forms a fiber section for 3d frame elements.
// The UniaxialFiber3d is subjected to a stress state with
// only one nonzero axial stress and corresponding axial strain.
//
// What: "@(#) UniaxialFiber3d.C, revA"

#include <stdlib.h>
#include <stdio.h>

#include <UniaxialMaterial.h>
#include <UniaxialFiber3d.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ID.h>
#include <SectionForceDeformation.h>
#include <Information.h>
#include <FiberResponse.h>

Matrix UniaxialFiber3d::ks(3,3); 
Vector UniaxialFiber3d::fs(3); 
ID UniaxialFiber3d::code(3); 

// constructor:
UniaxialFiber3d::UniaxialFiber3d()
:Fiber(0, FIBER_TAG_Uniaxial3d),
 theMaterial(0), area(0.0)
{
	if (code(0) != SECTION_RESPONSE_P) {
		code(0) = SECTION_RESPONSE_P;
		code(1) = SECTION_RESPONSE_MZ;
		code(2) = SECTION_RESPONSE_MY;
	}

   as[0] = 0.0;
   as[1] = 0.0;
}

UniaxialFiber3d::UniaxialFiber3d(int tag, 
                                 UniaxialMaterial &theMat,
                                 double Area, const Vector &position)
:Fiber(tag, FIBER_TAG_Uniaxial3d),
 theMaterial(0), area(Area)
{
	theMaterial = theMat.getCopy();  // get a copy of the MaterialModel

	if (theMaterial == 0)
		g3ErrorHandler->fatal("%s -- failed to get copy of UniaxialMaterial",
			"UniaxialFiber3d::UniaxialFiber2d");

   	if (code(0) != SECTION_RESPONSE_P) {
		code(0) = SECTION_RESPONSE_P;
		code(1) = SECTION_RESPONSE_MZ;
		code(2) = SECTION_RESPONSE_MY;
	}

	as[0] = -position(0);
	as[1] =  position(1);
}

// destructor:
UniaxialFiber3d::~UniaxialFiber3d ()
{
   if (theMaterial != 0)
      delete theMaterial;
}


int   
UniaxialFiber3d::setTrialFiberStrain(const Vector &vs)
{
  double strain = vs(0) + as[0]*vs(1) + as[1]*vs(2);

  if (theMaterial != 0)
      return theMaterial->setTrialStrain(strain);
  else {
      g3ErrorHandler->fatal("UniaxialFiber3d::setTrialFiberStrain() - %s\n",
			    "recvSelf() has not been invoked");
      return -1; // in case fatal does not exit
  }
}



// get fiber stress resultants 
Vector &
UniaxialFiber3d::getFiberStressResultants (void)
{
    double df = theMaterial->getStress() * area;

    // fs = as^ df;
    fs(0) = df;
    fs(1) = as[0]*df;
    fs(2) = as[1]*df;

    return fs;
}



// get contribution of fiber to section tangent stiffness
Matrix &
UniaxialFiber3d::getFiberTangentStiffContr(void) 
{
    // ks = (as^as) * area * Et;
    double value = theMaterial->getTangent() * area;

    double as1 = as[0];
    double as2 = as[1];
    double vas1 = as1*value;
    double vas2 = as2*value;
    double vas1as2 = vas1*as2;

    ks(0,0) = value;
    ks(0,1) = vas1;
    ks(0,2) = vas2;
    
    ks(1,0) = vas1;
    ks(1,1) = vas1*as1;
    ks(1,2) = vas1as2;
    
    ks(2,0) = vas2;
    ks(2,1) = vas1as2;
    ks(2,2) = vas2*as2;

    return ks;
}


Matrix &
UniaxialFiber3d::getFiberSecantStiffContr (void) 
{
    // ks = (as^as) * area * S;
    double value = theMaterial->getSecant() * area;

    double as1 = as[0];
    double as2 = as[1];
    double vas1 = as1*value;
    double vas2 = as2*value;

    ks(0,0) = value;
    ks(0,1) = vas1;
    ks(0,2) = vas2;
    
    ks(1,0) = vas1;
    ks(1,1) = vas1*as1;
    ks(1,2) = vas1*as2;
    
    ks(2,0) = vas2;
    ks(2,1) = vas2*as1;
    ks(2,2) = vas2*as2;
 
    return ks;
}


Fiber*
UniaxialFiber3d::getCopy (void)
{
   // make a copy of the fiber 
   static Vector position(2);

   position(0) = -as[0];
   position(1) =  as[1];

   UniaxialFiber3d *theCopy = new UniaxialFiber3d (this->getTag(), 
                                                   *theMaterial, area, 
                                                   position);
   return theCopy;
}  

int
UniaxialFiber3d::getOrder(void)
{
	return 3;
}

const ID&
UniaxialFiber3d::getType(void)
{
	return code;
}

int   
UniaxialFiber3d::commitState(void)
{
   return theMaterial->commitState();
}


int   
UniaxialFiber3d::revertToLastCommit(void)
{
   return theMaterial->revertToLastCommit();
}

int   
UniaxialFiber3d::revertToStart(void)
{
   return theMaterial->revertToStart();
}


int   
UniaxialFiber3d::sendSelf(int commitTag, Channel &theChannel)
{
    // 
    // store tag and material info in an ID and send it
    //

    static ID idData(3);
    int dbTag = this->getDbTag();
    idData(0) = this->getTag();
    idData(1) = theMaterial->getClassTag();
    int matDbTag = theMaterial->getDbTag();
    if (matDbTag == 0) {
	matDbTag = theChannel.getDbTag();
	if (matDbTag != 0)
	    theMaterial->setDbTag(matDbTag);
    }
    idData(2) = matDbTag;
    
    if (theChannel.sendID(dbTag, commitTag, idData) < 0)  {
	g3ErrorHandler->warning("UniaxialFiber3d::sendSelf() - %s\n",
			      "failed to send ID data");
	return -1;
    }    
    
    // 
    // store area and position data in a vector and send it
    //
    
    static Vector dData(3);
    dData(0) = area;
    dData(1) = as[0];
    dData(2) = as[1];
    if (theChannel.sendVector(dbTag, commitTag, dData) < 0)  {
	g3ErrorHandler->warning("UniaxialFiber3d::sendSelf() - %s\n",
				"failed to send Vector data");
	return -2;
    }    

    // now invoke sendSelf on the material
    if (theMaterial->sendSelf(commitTag, theChannel) < 0) {
	g3ErrorHandler->warning("UniaxialFiber3d::sendSelf() - %s\n",
				"the material failed in sendSelf()");
	return -3;
    }    	
    
    return 0;
}


int   
UniaxialFiber3d::recvSelf(int commitTag, Channel &theChannel, 
			  FEM_ObjectBroker &theBroker)
{
    // 
    // get tag and material info from an ID
    //

    static ID idData(3);
    int dbTag = this->getDbTag();
    
    if (theChannel.recvID(dbTag, commitTag, idData) < 0)  {
	g3ErrorHandler->warning("UniaxialFiber3d::recvSelf() - %s\n",
			      "failed to recv ID data");
	return -1;
    }    

    this->setTag(idData(0));

    // 
    // get area and position datafrom a vector
    //
    
    static Vector dData(3);
    if (theChannel.recvVector(dbTag, commitTag, dData) < 0)  {
	g3ErrorHandler->warning("UniaxialFiber3d::recvSelf() - %s\n",
				"failed to recv Vector data");
	return -2;
    }        
    area = dData(0);
    as[0] = dData(1);
    as[1] = dData(2);

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
	    g3ErrorHandler->warning("UniaxialFiber3d::recvSelf() - %s %d\n",
				    "failed to get a UniaxialMaterial of type",
				    matClassTag);
	    return -3;
	}
    }

    // set the materials dbTag and invoke recvSelf on the material
    theMaterial->setDbTag(idData(2));

    // now invoke recvSelf on the material
    if (theMaterial->recvSelf(commitTag, theChannel, theBroker) < 0) {
	g3ErrorHandler->warning("UniaxialFiber3d::recvSelf() - %s\n",
				"the material failed in recvSelf()");
	return -4;
    }    	

    return 0;
}


void UniaxialFiber3d::Print(ostream &s, int flag)
{
    s << "\nUniaxialFiber3d, tag: " << this->getTag() << endl;
    s << "\tArea: " << area << endl; 
    s << "\tMatrix as: " << 1.0 << ' ' << as[0] << ' ' << as[1] << endl; 
    s << "\tMaterial, tag: " << theMaterial->getTag() << endl;
}

Response*
UniaxialFiber3d::setResponse(char **argv, int argc, Information &info)
{
	if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)		
		return new FiberResponse(this, 1, Vector(3));

	else
		return theMaterial->setResponse(argv, argc, info);
}

int
UniaxialFiber3d::getResponse(int responseID, Information &fibInfo)
{
	switch(responseID) {
		case 1:
			return fibInfo.setVector(this->getFiberStressResultants());

		default:
			return -1;
	}
}

void 
UniaxialFiber3d::getFiberLocation(double &yLoc, double &zLoc)
{
	yLoc = -as[0];
	zLoc = as[1];
}
