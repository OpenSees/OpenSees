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
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/NDMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/NDMaterial.C
//
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for NDMaterial.
//
// What: "@(#) NDMaterial.C, revA"

#include <NDMaterial.h>
#include <Information.h>
#include <G3Globals.h>
#include <Matrix.h>
#include <Vector.h>

NDMaterial::NDMaterial(int tag, int classTag)
:Material(tag,classTag)
{

}


NDMaterial::~NDMaterial()
{

}

int 
NDMaterial::setTrialStrain(const Tensor &v)
{
	g3ErrorHandler->warning("%s -- not yet implemented for Tensors",
		"NDMaterial::setTrialStrain");

	return -1;
}

int 
NDMaterial::setTrialStrain(const Tensor &v, const Tensor &r)
{
	g3ErrorHandler->warning("%s -- not yet implemented for Tensors",
		"NDMaterial::setTrialStrain");

	return -1;
}

int 
NDMaterial::setTrialStrainIncr(const Tensor &v)
{
	g3ErrorHandler->warning("%s -- not yet implemented for Tensors",
		"NDMaterial::setTrialStrainIncr");

	return -1;
}

int 
NDMaterial::setTrialStrainIncr(const Tensor &v, const Tensor &r)
{
	g3ErrorHandler->warning("%s -- not yet implemented for Tensors",
		"NDMaterial::setTrialStrainIncr");

	return -1;
}

const Tensor&
NDMaterial::getTangentTensor(void)
{
	g3ErrorHandler->warning("%s -- not yet implemented",
		"NDMaterial::getTangentTensor");

	Tensor *newTensor = new Tensor;

	return *newTensor;
}

const Tensor&
NDMaterial::getStressTensor(void)
{
	g3ErrorHandler->warning("%s -- not yet implemented",
		"NDMaterial::getStressTensor");

	Tensor *newTensor = new Tensor;

	return *newTensor;

}

const Tensor&
NDMaterial::getStrainTensor(void)
{
	g3ErrorHandler->warning("%s -- not yet implemented",
		"NDMaterial::getStrainTensor");

	Tensor *newTensor = new Tensor;

	return *newTensor;

}

int
NDMaterial::setResponse (char **argv, int argc, Information &matInfo)
{
    if (strcmp(argv[0],"stress") ==0 || strcmp(argv[0],"stresses") == 0) {
		Vector *newVector = new Vector(this->getStress());
		if (newVector == 0) {
			g3ErrorHandler->warning("WARNING NDMaterial::setResponse() - %d out of memory creating vector\n",
				    this->getTag());
			return -1;
		}
		matInfo.theVector = newVector;
		matInfo.theType = VectorType;
		return 1;
    } 

	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
		Vector *newVector = new Vector(this->getStrain());
		if (newVector == 0) {
			g3ErrorHandler->warning("WARNING NDMaterial::setResponse() - %d out of memory creating vector\n",
				    this->getTag());
			return -1;
		}
		matInfo.theVector = newVector;
		matInfo.theType = VectorType;
		return 2;
    } 

    else if (strcmp(argv[0],"tangent") == 0) {
		Matrix *newMatrix = new Matrix(this->getTangent());
		if (newMatrix == 0) {
			g3ErrorHandler->warning("WARNING NDMaterial::setResponse() - %d out of memory creating matrix\n",
				    this->getTag());
			return -1;
		}
		matInfo.theMatrix = newMatrix;
		matInfo.theType = MatrixType;
		return 3;
    } 

    else
		return -1;

}

int 
NDMaterial::getResponse (int responseID, Information &matInfo)
{
	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = this->getStress();
			return 0;
		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = this->getStrain();
			return 0;
		case 3:
			if (matInfo.theMatrix != 0)
				*(matInfo.theMatrix) = this->getTangent();
			return 0;
		default:
			return -1;
	}
}
