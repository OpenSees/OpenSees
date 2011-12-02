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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-14 23:01:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/GenericSectionNd.cpp,v $
                                                                        
                                                                        
// File: ~/material/GenericSectionNd.C
//
// Written: MHS 
// Created: Apr 2000
// Revision: A
//
// Description: This file contains the class implementation for GenericSectionNd.
//
// What: "@(#) GenericSectionNd.C, revA"

#include <GenericSectionNd.h>
#include <NDMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <G3Globals.h>
#include <classTags.h>

#include <string.h>

GenericSectionNd::GenericSectionNd(int tag, NDMaterial &m, const ID &mCode)
:SectionForceDeformation(tag,SEC_TAG_GenericNd),
otherDbTag(0), theModel(0), code(0)
{
    theModel = m.getCopy();

    if (theModel == 0) {
		g3ErrorHandler->fatal("%s -- failed to get copy of material model",
			"GenericSectionNd::GenericSectionNd");
    }

	order = theModel->getOrder();

	code = new ID(mCode);

	if (code == 0) {
		g3ErrorHandler->fatal("%s -- failed to allocate section ID",
			"GenericSectionNd::GenericSectionNd");
	}

    if (order != code->Size()) {
		g3ErrorHandler->warning("%s -- code size does not match order of material model",
			"GenericSectionNd::GenericSectionNd");
    }
}

GenericSectionNd::GenericSectionNd()
:SectionForceDeformation(0,SEC_TAG_GenericNd), 
otherDbTag(0), theModel(0), code(0), order(0)
{
    
}

GenericSectionNd::~GenericSectionNd()
{
	if (theModel)
		delete theModel;

	if (code)
		delete code;
}

int
GenericSectionNd::setTrialSectionDeformation (const Vector& def)
{
    return theModel->setTrialStrain(def);
}

const Vector&
GenericSectionNd::getSectionDeformation ()
{
    return theModel->getStrain();
}

const Vector&
GenericSectionNd::getStressResultant ()
{
    return theModel->getStress();
}

const Matrix&
GenericSectionNd::getSectionTangent ()
{
    return theModel->getTangent();
}

int
GenericSectionNd::commitState ()
{
    return theModel->commitState();
}

int
GenericSectionNd::revertToLastCommit ()
{
    return theModel->revertToLastCommit();
}

int
GenericSectionNd::revertToStart ()
{
    return theModel->revertToStart();
}

const ID&
GenericSectionNd::getType ()
{
    return *code;
}

int
GenericSectionNd::getOrder () const
{
    return order;
}

SectionForceDeformation*
GenericSectionNd::getCopy ()
{
    GenericSectionNd *theCopy = new GenericSectionNd (this->getTag(),
						      *theModel, *code);

	theCopy->otherDbTag = otherDbTag;

    return theCopy;
}

int
GenericSectionNd::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;

	// Need otherDbTag since code ID and data ID may be the same size
    if (otherDbTag == 0) 
      otherDbTag = theChannel.getDbTag();

	static ID data(5);

	data(0) = this->getTag();
	data(1) = order;
	data(2) = otherDbTag;
	data(3) = theModel->getClassTag();

	int dbTag = theModel->getDbTag();

	if (dbTag == 0) {
		dbTag = theChannel.getDbTag();
		if (dbTag != 0)
			theModel->setDbTag(dbTag);
	}

	data(4) = dbTag;

	// Send the ID vector
	res += theChannel.sendID(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send data ID",
			"GenericSectionNd::sendSelf");
		return res;
	}
    
	// Send the section code
	res += theChannel.sendID(otherDbTag, cTag, *code);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send code ID",
			"GenericSectionNd::sendSelf");
		return res;
	}

	// Ask the NDMaterial to send itself
	res += theModel->sendSelf(cTag, theChannel);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send NDMaterial",
			"GenericSectionNd::sendSelf");
		return res;
	}

    return res;
}

int
GenericSectionNd::recvSelf(int cTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
	int res = 0;

    static ID data(5);

	// Receive the data ID
    res += theChannel.recvID(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive data ID",
			"GenericSectionNd::recvSelf");
		return res;
	}

	this->setTag(data(0));
	order = data(1);
	otherDbTag = data(2);

	// Check if section ID code is null or wrong size, reallocate if so
	if (code == 0)
		code = new ID(order);
	else if (code->Size() != order) {
		delete code;
		code = new ID(order);
	}
	if (code == 0) {
		g3ErrorHandler->warning("%s -- could not allocate new code ID",
			"GenericSectionNd::recvSelf");
		return -1;
	}

	// Receive the code ID
    res += theChannel.recvID(otherDbTag, cTag, *code);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive code ID",
			"GenericSectionNd::recvSelf");
		return res;
	}

	int classTag = data(3);

	// Check if the material is null; if so, get a new one
	if (theModel == 0)
		theModel = theBroker.getNewNDMaterial(classTag);

	// Check that the material is of the right type; if not, delete
	// the current one and get a new one of the right type
	else if (theModel->getClassTag() != classTag) {
		delete theModel;
		theModel = theBroker.getNewNDMaterial(classTag);
	}

	// Check if either allocation failed
	if (theModel == 0) {
		g3ErrorHandler->warning("%s -- could not get an NDMaterial",
			"GenericSectionNd::recvSelf");
		return -1;
	}

	// Now, receive the material
	theModel->setDbTag(data(4));
	res += theModel->recvSelf(cTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive NDMaterial",
			"GenericSectionNd::recvSelf");
		return res;
	}

    return res;
}

void
GenericSectionNd::Print (OPS_Stream &s, int flag)
{
    s << "Generic Section Nd, tag: " << this->getTag() << endln;
    s << "\tsection code: " << code << endln;
}
