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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/UniaxialFiber2d.cpp,v $
                                                                        
                                                                        
// File: ~/section/UniaxialFiber2d.C
//
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

Matrix UniaxialFiber2d::ks(2,2); 
Vector UniaxialFiber2d::fs(2); 
ID UniaxialFiber2d::code(2);

// constructor:
UniaxialFiber2d::UniaxialFiber2d(int tag, 
                                 UniaxialMaterial &theMat,
                                 double Area, double position):
                                 Fiber(tag, FIBER_TAG_Uniaxial2d),
                                 theMaterial(0), area(Area), y(-position)
{
   theMaterial = theMat.getCopy();  // get a copy of the MaterialModel

   if (theMaterial == 0)
	   g3ErrorHandler->fatal("%s -- failed to get copy of UniaxialMaterial",
		   "UniaxialFiber2d::UniaxialFiber2d");

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

Matrix &
UniaxialFiber2d::getFiberSecantStiffContr (void) 
{
    // Use the section kinematic matrix to get the fiber 
    // secant stiffness matrix
    // ks = (as^as) * area * S;
    double value = theMaterial->getSecant() * area;
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
		g3ErrorHandler->warning("%s - failed to send ID data",
			"UniaxialFiber2d::sendSelf");
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
		g3ErrorHandler->warning("%s - failed to send Vector data",
			"UniaxialFiber2d::sendSelf");
		return res;
    }    

    // now invoke sendSelf on the material
    res += theMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
		g3ErrorHandler->warning("%s - failed to send UniaxialMaterial",
			"UniaxialFiber2d::sendSelf");
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
		g3ErrorHandler->warning("%s - failed to receive ID data",
		      "UniaxialFiber2d::recvSelf");
		return res;
    }    

    this->setTag(idData(0));

    // 
    // get area from a vector received from channel
    //
    
    static Vector dData(2);

    res += theChannel.recvVector(dbTag, commitTag, dData);
	if (res < 0) {
		g3ErrorHandler->warning("%s - failed to receive Vector data",
		      "UniaxialFiber2d::recvSelf");
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
			g3ErrorHandler->warning("UniaxialFiber2d::recvSelf() - %s %d\n",
			    "failed to get a UniaxialMaterial of type", matClassTag);
		    return -1;
		}
    }

    // set the materials dbTag and invoke recvSelf on the material
    theMaterial->setDbTag(idData(2));

    // now invoke recvSelf on the material
    res += theMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("UniaxialFiber2d::recvSelf() - %s\n",
			"the material failed in recvSelf()");
		return res;
    }    	

	return res;
}


void UniaxialFiber2d::Print(ostream &s, int flag)
{
    s << "\nUniaxialFiber2d, tag: " << this->getTag() << endl;
    s << "\tArea: " << area << endl; 
    s << "\tMatrix as: " << 1.0 << ' ' << y << endl; 
    s << "\tMaterial, tag: " << theMaterial->getTag() << endl;
}
