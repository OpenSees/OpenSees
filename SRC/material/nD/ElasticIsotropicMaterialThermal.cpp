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

// $Revision: 1.25 $                                                              
// $Date: 2009-01-29 00:42:03 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicMaterialThermal.cpp,v $                                                                
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for ElasticIsotropicMaterialThermal.
//
// What: "@(#) ElasticIsotropicMaterialThermal.C, revA"
//Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io]

#include <string.h>

#include <ElasticIsotropicMaterialThermal.h>
#include <ElasticIsotropicPlaneStress2D.h>
#include <ElasticIsotropicPlaneStrain2D.h>
#include <ElasticIsotropicAxiSymm.h>
#include <ElasticIsotropicThreeDimensional.h>
#include <ElasticIsotropic3DThermal.h> //added by Liming,UoE
#include <ElasticIsotropicPlateFiber.h>
#include <ElasticIsotropicBeamFiber.h>
#include <ElasticIsotropicBeamFiber2d.h>

#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>


void * OPS_ADD_RUNTIME_VPV(OPS_ElasticIsotropicMaterialThermal)
{
	NDMaterial *theMaterial = 0;
	int softindex=0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	if (numArgs < 3) {
		opserr << "Want: nDMaterial ElasticIsotropic $tag $E $V <$rho> <$alpha> <-cSoft/-sSoft> " << endln;
		return 0;
	}

	int iData[1];
	double dData[4];
	dData[2] = 0.0;
	dData[3] = 0.0;

	int numData = 1;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid integer tag: nDMaterial EasticIsotropic \n";
		return 0;
	}

	if (numArgs > 4)
		numData = 4;
	else
		numData = 2;

	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid data: nDMaterial EasticIsotropic : " << iData[0] << "\n";
		return 0;
	}
	

	//added by Liming, UoE,2014;
	if (numArgs == 6) {
		const char* typeChar = OPS_GetString();
		if ((strcmp(typeChar, "-SteelSoft") == 0) || (strcmp(typeChar, "-SSoft") == 0)||( strcmp(typeChar, "-sSoft") == 0)) {
			softindex = 1;
		}
		else if ((strcmp(typeChar, "-ConcreteSoft") == 0) || (strcmp(typeChar, "-CSoft") == 0) || (strcmp(typeChar, "-cSoft") == 0)) {
			softindex = 2;
		}
	}
	if(numData==4)
		theMaterial = new ElasticIsotropicMaterialThermal(iData[0], dData[0], dData[1], dData[2], dData[3],softindex);		
	else
		theMaterial = new ElasticIsotropicMaterialThermal(iData[0], dData[0], dData[1], dData[2]);

	return theMaterial;
}


//double ElasticIsotropicMaterialThermal::Alpha = 0;

ElasticIsotropicMaterialThermal::ElasticIsotropicMaterialThermal
(int tag, double e, double nu, double r, double alpha, int softindex)
	: NDMaterial(tag, ND_TAG_ElasticIsotropicThermal), E(e), v(nu), rho(r), Alpha(alpha),SoftIndex(softindex)
{

}


ElasticIsotropicMaterialThermal::ElasticIsotropicMaterialThermal
(int tag, int classTag, double e, double nu, double r, double alpha, int softindex)
	:NDMaterial(tag, classTag), E(e), v(nu), rho(r), Alpha(alpha), SoftIndex(softindex)
{

}

ElasticIsotropicMaterialThermal::~ElasticIsotropicMaterialThermal()
{

}

double
ElasticIsotropicMaterialThermal::getRho()
{
	return rho;
}

NDMaterial*
ElasticIsotropicMaterialThermal::getCopy(const char *type)
{
	if (strcmp(type, "ThreeDimensionalThermal") == 0 || strcmp(type, "3DThermal") == 0) {
		ElasticIsotropic3DThermal *theModel;
		theModel = new ElasticIsotropic3DThermal(this->getTag(), E, v, rho, Alpha,SoftIndex);
		return theModel;
	}


	// Handle other cases
	else
		return NDMaterial::getCopy(type);
}

int
ElasticIsotropicMaterialThermal::setTrialStrain(const Vector &v)
{
	opserr << "ElasticIsotropicMaterialThermal::setTrialStrain -- subclass responsibility\n";
	exit(-1);
	return -1;
}

int
ElasticIsotropicMaterialThermal::setTrialStrain(const Vector &v, const Vector &rate)
{
	opserr << "ElasticIsotropicMaterialThermal::setTrialStrain -- subclass responsibility\n";
	exit(-1);
	return -1;
}

int
ElasticIsotropicMaterialThermal::setTrialStrainIncr(const Vector &v)
{
	opserr << "ElasticIsotropicMaterialThermal::setTrialStrainIncr -- subclass responsibility\n";
	exit(-1);
	return -1;
}

int
ElasticIsotropicMaterialThermal::setTrialStrainIncr(const Vector &v, const Vector &rate)
{
	opserr << "ElasticIsotropicMaterialThermal::setTrialStrainIncr -- subclass responsibility\n";
	exit(-1);
	return -1;
}

const Matrix&
ElasticIsotropicMaterialThermal::getTangent(void)
{
	opserr << "ElasticIsotropicMaterialThermal::getTangent -- subclass responsibility\n";
	exit(-1);

	// Just to make it compile
	Matrix *ret = new Matrix();
	return *ret;
}

const Matrix&
ElasticIsotropicMaterialThermal::getInitialTangent(void)
{
	return this->getTangent();
}

const Vector&
ElasticIsotropicMaterialThermal::getStress(void)
{
	opserr << "ElasticIsotropicMaterialThermal::getStress -- subclass responsibility\n";
	exit(-1);

	// Just to make it compile
	Vector *ret = new Vector();
	return *ret;
}

const Vector&
ElasticIsotropicMaterialThermal::getStrain(void)
{
	opserr << "ElasticIsotropicMaterialThermal::getStrain -- subclass responsibility\n";
	exit(-1);

	// Just to make it compile
	Vector *ret = new Vector();
	return *ret;
}

int
ElasticIsotropicMaterialThermal::commitState(void)
{
	opserr << "ElasticIsotropicMaterialThermal::commitState -- subclass responsibility\n";
	exit(-1);
	return -1;
}

int
ElasticIsotropicMaterialThermal::revertToLastCommit(void)
{
	opserr << "ElasticIsotropicMaterialThermal::revertToLastCommit -- subclass responsibility\n";
	exit(-1);

	return -1;
}

int
ElasticIsotropicMaterialThermal::revertToStart(void)
{
	opserr << "ElasticIsotropicMaterialThermal::revertToStart -- subclass responsibility\n";
	exit(-1);
	return -1;
}

NDMaterial*
ElasticIsotropicMaterialThermal::getCopy(void)
{
	opserr << "ElasticIsotropicMaterialThermal::getCopy -- subclass responsibility\n";
	exit(-1);
	return 0;
}

const char*
ElasticIsotropicMaterialThermal::getType(void) const
{
	opserr << "ElasticIsotropicMaterialThermal::getType -- subclass responsibility\n";
	exit(-1);

	return 0;
}

int
ElasticIsotropicMaterialThermal::getOrder(void) const
{
	opserr << "ElasticIsotropicMaterialThermal::getOrder -- subclass responsibility\n";
	exit(-1);
	return -1;
}

int
ElasticIsotropicMaterialThermal::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(4);

	data(0) = this->getTag();
	data(1) = E;
	data(2) = v;
	data(3) = rho;

	res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "ElasticIsotropicMaterialThermal::sendSelf -- could not send Vector\n";
		return res;
	}

	return res;
}

int
ElasticIsotropicMaterialThermal::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;

	static Vector data(4);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "ElasticIsotropicMaterialThermal::recvSelf -- could not recv Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	E = data(1);
	v = data(2);
	rho = data(3);

	return res;
}

void
ElasticIsotropicMaterialThermal::Print(OPS_Stream &s, int flag)
{
	s << "Elastic Isotropic Material Model" << endln;
	s << "\tE:  " << E << endln;
	s << "\tv:  " << v << endln;
	s << "\trho:  " << rho << endln;

	return;
}

int
ElasticIsotropicMaterialThermal::setParameter(const char **argv, int argc,
	Parameter &param)
{
	if (strcmp(argv[0], "E") == 0)
		return param.addObject(1, this);

	else if (strcmp(argv[0], "nu") == 0 || strcmp(argv[0], "v") == 0)
		return param.addObject(2, this);

	else if (strcmp(argv[0], "rho") == 0)
		return param.addObject(3, this);

	return -1;
}

int
ElasticIsotropicMaterialThermal::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		E = info.theDouble;
		return 0;
	case 2:
		v = info.theDouble;
		return 0;
	case 3:
		rho = info.theDouble;
		return 0;
	default:
		return -1;
	}
}
