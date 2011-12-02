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
                                                                        
// $Revision: 1.3 $                                                              
// $Date: 2000-10-10 05:31:36 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicMaterial.cpp,v $                                                                
                                                                        
                                                                        
// File: ~/material/ElasticIsotropicMaterial.C
//
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for ElasticIsotropicMaterial.
//
// What: "@(#) ElasticIsotropicMaterial.C, revA"

#include <string.h>

#include <ElasticIsotropicMaterial.h>
#include <ElasticIsotropicPlaneStress2D.h>
#include <ElasticIsotropicPlaneStrain2D.h>
#include <ElasticIsotropic3D.h>

#include <Tensor.h>

#include <G3Globals.h>

ElasticIsotropicMaterial::ElasticIsotropicMaterial
(int tag, int classTag, double e, double nu)
:NDMaterial(tag, classTag), E(e), v(nu)
{

}

ElasticIsotropicMaterial::ElasticIsotropicMaterial
(int tag, double e, double nu)
:NDMaterial(tag, ND_TAG_ElasticIsotropic), E(e), v(nu)
{

}

ElasticIsotropicMaterial::~ElasticIsotropicMaterial()
{

}

NDMaterial*
ElasticIsotropicMaterial::getCopy (const char *type)
{
    if (strcmp(type,"PlaneStress2D") == 0)
    {
	ElasticIsotropicPlaneStress2D *theModel;
	theModel = new ElasticIsotropicPlaneStress2D (this->getTag(), E, v);
		// DOES NOT COPY sigma, D, and epsilon ...
		// This function should only be called during element instantiation, so
		// no state determination is performed on the material model object
		// prior to copying the material model (calling this function)
	return theModel;
    }

    else if (strcmp(type,"PlaneStrain2D") == 0)
    {
	ElasticIsotropicPlaneStrain2D *theModel;
	theModel = new ElasticIsotropicPlaneStrain2D (this->getTag(), E, v);
		// DOES NOT COPY sigma, D, and epsilon ...
		// This function should only be called during element instantiation, so
		// no state determination is performed on the material model object
		// prior to copying the material model (calling this function)
	return theModel;
    }
    else if (strcmp(type,"ElasticIsotropic3D") == 0)
    {
	ElasticIsotropic3D *theModel;
	theModel = new ElasticIsotropic3D (this->getTag(), E, v);
		// DOES NOT COPY sigma, D, and epsilon ...
		// This function should only be called during element instantiation, so
		// no state determination is performed on the material model object
		// prior to copying the material model (calling this function)
	return theModel;
    }

    // Handle other cases
    else
    {
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::getModel failed to get model %s",
			      type);

	return 0;
    }
}

int
ElasticIsotropicMaterial::setTrialStrain (const Vector &v)
{
    g3ErrorHandler->fatal("ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility");

    return 0;
}

int
ElasticIsotropicMaterial::setTrialStrain (const Vector &v, const Vector &rate)
{
    g3ErrorHandler->fatal("ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility");

    return 0;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Vector &v)
{
    g3ErrorHandler->fatal("ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility");

    return 0;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
    g3ErrorHandler->fatal("ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility");

    return 0;
}

const Matrix&
ElasticIsotropicMaterial::getTangent (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::getTangent -- subclass responsibility");

	// Just to make it compile
	Matrix *ret = new Matrix();
	return *ret;
}

const Vector&
ElasticIsotropicMaterial::getStress (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::getStress -- subclass responsibility");

	// Just to make it compile
	Vector *ret = new Vector();
	return *ret;
}

const Vector&
ElasticIsotropicMaterial::getStrain (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::getStrain -- subclass responsibility");

	// Just to make it compile
	Vector *ret = new Vector();
	return *ret;
}

int
ElasticIsotropicMaterial::setTrialStrain (const Tensor &v)
{
    g3ErrorHandler->fatal("ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility");

    return 0;
}

int
ElasticIsotropicMaterial::setTrialStrain (const Tensor &v, const Tensor &r)
{
    g3ErrorHandler->fatal("ElasticIsotropicMaterial::setTrialStrain -- subclass responsibility");

    return 0;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Tensor &v)
{
    g3ErrorHandler->fatal("ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility");

    return 0;
}

int
ElasticIsotropicMaterial::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    g3ErrorHandler->fatal("ElasticIsotropicMaterial::setTrialStrainIncr -- subclass responsibility");

    return 0;
}

const Tensor&
ElasticIsotropicMaterial::getTangentTensor (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::getTangent -- subclass responsibility");

	// Just to make it compile
	Tensor *t = new Tensor;
	return *t;
}

const Tensor&
ElasticIsotropicMaterial::getStressTensor (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::getStress -- subclass responsibility");

	// Just to make it compile
     Tensor *t = new Tensor;
	 return *t;
}

const Tensor&
ElasticIsotropicMaterial::getStrainTensor (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::getStrain -- subclass responsibility");

	// Just to make it compile
     Tensor *t = new Tensor;
	 return *t;
}

int
ElasticIsotropicMaterial::commitState (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::commitState -- subclass responsibility");

	return 0;
}

int
ElasticIsotropicMaterial::revertToLastCommit (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::revertToLastCommit -- subclass responsibility");

	return 0;
}

int
ElasticIsotropicMaterial::revertToStart (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::revertToStart -- subclass responsibility");

	return 0;
}

NDMaterial*
ElasticIsotropicMaterial::getCopy (void)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::getCopy -- subclass responsibility");

	return 0;
}

const char*
ElasticIsotropicMaterial::getType (void) const
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::getType -- subclass responsibility");

	return 0;
}

int
ElasticIsotropicMaterial::getOrder (void) const
{
	 g3ErrorHandler->fatal("ElasticIsotropicMaterial::getOrder -- subclass responsibility");

	return 0;
}

int
ElasticIsotropicMaterial::sendSelf (int commitTag, Channel &theChannel)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::sendSelf -- subclass responsibility");

	return 0;
}

int
ElasticIsotropicMaterial::recvSelf (int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
	g3ErrorHandler->fatal("ElasticIsotropicMaterial::recvSelf -- subclass responsibility");

	return 0;
}

void
ElasticIsotropicMaterial::Print (ostream &s, int flag)
{
	s << "Elastic Isotropic Material Model" << endl;
	s << "\tE:  " << E << endl;
	s << "\tv:  " << v << endl;

	return;
}

int 
ElasticIsotropicMaterial::setParameter(char **argv, int argc, Information &info)
{
  return -1;
}

int 
ElasticIsotropicMaterial::updateParameter(int parameterID, Information &info)
{ 
  return -1;
}
