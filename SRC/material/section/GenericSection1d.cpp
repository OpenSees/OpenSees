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
                                                                        
// $Revision: 1.6 $
// $Date: 2003-02-14 23:01:34 $
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

#include <classTags.h>

#include <string.h>

Vector GenericSection1d::s(1);
Matrix GenericSection1d::ks(1,1);
ID GenericSection1d::c(1);

GenericSection1d::GenericSection1d(int tag, UniaxialMaterial &m, int type)
:SectionForceDeformation(tag,SEC_TAG_Generic1d), code(type)
{
    theModel = m.getCopy();

    if (!theModel) {
      opserr << "GenericSection1d::GenericSection1d  -- failed to get copy of material model\n";
      exit(-1);
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
  s(0) = theModel->getStress();

  return s;
}

const Matrix&
GenericSection1d::getSectionTangent ()
{
  ks(0,0) = theModel->getTangent();
  
  return ks;
}

const Matrix&
GenericSection1d::getInitialTangent ()
{
  ks(0,0) = theModel->getInitialTangent();
  
  return ks;
}

const Matrix&
GenericSection1d::getSectionFlexibility ()
{
  double tangent = theModel->getTangent();

  if (tangent != 0.0)
    ks(0,0) = 1.0/tangent;
  else
    ks(0,0) = 1.0e12;

  return ks;
}

const Matrix&
GenericSection1d::getInitialFlexibility ()
{
  double tangent = theModel->getInitialTangent();

  ks(0,0) = 1.0/tangent;

  return ks;
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
GenericSection1d::getType ()
{
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
	  opserr << "GenericSection1d::sendSelf -- could not send ID\n";
	  return res;
	}
    
	// Ask the UniaxialMaterial to send itself
	res += theModel->sendSelf(cTag, theChannel);
	if (res < 0) {
	  opserr << "GenericSection1d::sendSelf -- could not send UniaxialMaterial\n";
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
	  opserr << "GenericSection1d::recvSelf -- could not receive ID\n";
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
	  opserr << "GenericSection1d::recvSelf -- could not get a UniaxialMaterial\n";
	  return -1;
	}

	// Now, receive the material
	theModel->setDbTag(data(3));
	res += theModel->recvSelf(cTag, theChannel, theBroker);
	if (res < 0) {
	  opserr << "GenericSection1d::recvSelf -- could not receive UniaxialMaterial\n";
	  return res;
	}

    return res;
}

void
GenericSection1d::Print (OPS_Stream &s, int flag)
{
    s << "GenericSection1d (Uniaxial), tag: " << this->getTag() << endln;
    s << "\tResponse code: " << code << endln;
    s << "\tUniaxialMaterial: " << theModel->getTag() << endln;
}
