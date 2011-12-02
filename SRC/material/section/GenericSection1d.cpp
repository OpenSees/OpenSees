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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/GenericSection1d.cpp,v $
                                                                        
                                                                        
// File: ~/material/GenericSection1d.C
//
// Written: MHS 
// Created: Apr 2000
// Revision: A
//
// Description: This file contains the class implementation for GenericSection1d.
//
// What: "@(#) GenericSection1d.C, revA"

#include <GenericSection1d.h>
#include <UniaxialMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <G3Globals.h>
#include <classTags.h>

#include <string.h>

GenericSection1d::GenericSection1d(int tag, UniaxialMaterial &m, int c)
:SectionForceDeformation(tag,SEC_TAG_Generic1d), code(c)
{
    theModel = m.getCopy();

    if (!theModel) {
		g3ErrorHandler->fatal("%s -- failed to get copy of material model",
			"GenericSection1d::GenericSection1d");
    }
}

GenericSection1d::GenericSection1d()
:SectionForceDeformation(0,SEC_TAG_Generic1d), theModel(0), code(0)
{
    
}

GenericSection1d::~GenericSection1d()
{
    if (theModel)
		delete theModel;
}

int
GenericSection1d::setTrialSectionDeformation (const Vector& def)
{
	return theModel->setTrialStrain(def(0));
}

const Vector&
GenericSection1d::getSectionDeformation ()
{
	static Vector e(1);	// static for class-wide returns

	e(0) = theModel->getStrain();

    return e;
}

const Vector&
GenericSection1d::getStressResultant ()
{
	static Vector s(1);	// static for class-wide returns

    s(0) = theModel->getStress();

    return s;
}

const Matrix&
GenericSection1d::getSectionTangent ()
{
	static Matrix k(1,1);	// static for class-wide returns

    k(0,0) = theModel->getTangent();

    return k;
}

const Matrix&
GenericSection1d::getSectionFlexibility ()
{
	static Matrix f(1,1);	// static for class-wide returns

    double tangent = theModel->getTangent();

    if (tangent != 0.0)
		f(0,0) = 1.0/tangent;
	else
		f(0,0) = 1.0e12;

    return f;
}

int
GenericSection1d::commitState ()
{
    return theModel->commitState();
}

int
GenericSection1d::revertToLastCommit ()
{
    return theModel->revertToLastCommit();
}

int
GenericSection1d::revertToStart ()
{
    return theModel->revertToStart();
}

const ID&
GenericSection1d::getType () const
{
	static ID c(1);	// static for class-wide returns

	c(0) = code;

    return c;
}

int
GenericSection1d::getOrder () const
{
    return 1;
}

SectionForceDeformation*
GenericSection1d::getCopy ()
{
    GenericSection1d *theCopy = new GenericSection1d (this->getTag(),
						      *theModel, code);

    return theCopy;
}

int
GenericSection1d::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;

	static ID data(4);

	data(0) = this->getTag();
	data(1) = code;
	data(2) = theModel->getClassTag();

	int dbTag = theModel->getDbTag();

	if (dbTag == 0) {
		dbTag = theChannel.getDbTag();
		if (dbTag != 0)
			theModel->setDbTag(dbTag);
	}

	data(3) = dbTag;

	// Send the ID vector
	res += theChannel.sendID(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send ID",
			"GenericSection1d::sendSelf");
		return res;
	}
    
	// Ask the UniaxialMaterial to send itself
	res += theModel->sendSelf(cTag, theChannel);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send UniaxialMaterial",
			"GenericSection1d::sendSelf");
		return res;
	}

    return res;
}

int
GenericSection1d::recvSelf(int cTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
	int res = 0;

    static ID data(4);

    res += theChannel.recvID(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive ID",
			"GenericSection1d::recvSelf");
		return res;
	}

	this->setTag(data(0));
	code = data(1);
	int classTag = data(2);

	// Check if the material is null; if so, get a new one
	if (theModel == 0)
		theModel = theBroker.getNewUniaxialMaterial(classTag);

	// Check that the material is of the right type; if not, delete
	// the current one and get a new one of the right type
	else if (theModel->getClassTag() != classTag) {
		delete theModel;
		theModel = theBroker.getNewUniaxialMaterial(classTag);
	}

	// Check if either allocation failed
	if (theModel == 0) {
		g3ErrorHandler->warning("%s -- could not get a UniaxialMaterial",
			"GenericSection1d::recvSelf");
		return -1;
	}

	// Now, receive the material
	theModel->setDbTag(data(3));
	res += theModel->recvSelf(cTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive UniaxialMaterial",
			"GenericSection1d::recvSelf");
		return res;
	}

    return res;
}

void
GenericSection1d::Print (ostream &s, int flag)
{
    s << "GenericSection1d, tag: " << this->getTag() << endl;
    s << "\tResponse code: " << code << endl;
}
