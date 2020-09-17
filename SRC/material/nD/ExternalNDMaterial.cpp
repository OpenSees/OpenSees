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

// $Revision: 1.0 $
// $Date: 2019-01-26 00:00:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ExternalNDMaterial.cpp,v $

// Written: M. Salehi
// Created: 26-1 2019
//


#include "ExternalNDMaterial.h"
#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>

ExternalNDMaterial::ExternalNDMaterial (int tag)			
	: NDMaterial(tag, ND_TAG_ExternalNDMaterial)
{
		
}

// Destructor
ExternalNDMaterial::~ExternalNDMaterial()
{
	
}

// get copy
NDMaterial* ExternalNDMaterial::getCopy(void) 
{
	return _NDMGetCopy();
}


// get copy
NDMaterial* ExternalNDMaterial::getCopy(const char *type)
{
	return _NDMGetCopy_Type(type);
}

// Print 
void ExternalNDMaterial::Print(OPS_Stream &s, int flag)
{
	_NDMPrint(flag);
}

int ExternalNDMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;
	return res;
}

int ExternalNDMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	return res;
}

double ExternalNDMaterial::getRho(void)
{
	return _NDMGetRho();
}

int ExternalNDMaterial::setTrialStrain(const Vector &v)
{
	return _NDMSetTrialStrain_V(v);
}

int ExternalNDMaterial::setTrialStrain(const Vector &v, const Vector &r)
{
	return _NDMSetTrialStrain_VR(v, r);
}

int ExternalNDMaterial::setTrialStrainIncr(const Vector &v)
{
	return _NDMSetTrialStrainIncr_V(v);
}

int ExternalNDMaterial::setTrialStrainIncr(const Vector &v, const Vector &r)
{
	return _NDMSetTrialStrainIncr_VR(v,r);
}

const Matrix& ExternalNDMaterial::getTangent (void)
{
	return *_NDMGetTangent();
}

const Matrix& ExternalNDMaterial::getInitialTangent(void)
{
	return *_NDMGetInitialTangent();
}

const Vector& ExternalNDMaterial::getStress(void)
{
	return *_NDMGetStress();
}

const Vector& ExternalNDMaterial::getStrain()
{
	return *_NDMGetStrain();
}

int ExternalNDMaterial::revertToStart(void)
{
	return _NDMRevertToStart();
}


int ExternalNDMaterial::commitState(void)
{
	return _NDMCommitState();
}

int ExternalNDMaterial::revertToLastCommit(void)
{
	return _NDMRevertToLastCommit();
}

const char *ExternalNDMaterial::getType(void) const {
	return _NDMGetType();
}

int ExternalNDMaterial::getOrder(void) const {
	return _NDMGetOrder();
}
