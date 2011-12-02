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
// $Date: 2000-12-18 10:48:08 $                                                                  
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
#include <MaterialResponse.h>

Matrix NDMaterial::errMatrix(1,1);
Vector NDMaterial::errVector(1);
Tensor NDMaterial::errTensor(2, def_dim_2, 0.0 );

NDMaterial::NDMaterial(int tag, int classTag)
:Material(tag,classTag)
{

}

NDMaterial::NDMaterial()
:Material(0, 0)
{

}

NDMaterial::~NDMaterial()
{

}

const Vector &
NDMaterial::getCommittedStress(void) 
{
  return this->getStress();
}

const Vector &
NDMaterial::getCommittedStrain(void) 
{
  return this->getStrain();
}

// methods to set and retrieve state.
int 
NDMaterial::setTrialStrain(const Vector &v)
{
   g3ErrorHandler->fatal("NDMaterial::setTrialStrain -- subclass responsibility");
   return 0;    
}

int 
NDMaterial::setTrialStrain(const Vector &v, const Vector &r)
{
   g3ErrorHandler->fatal("NDMaterial::setTrialStrain -- subclass responsibility");
   return 0;    
}

int 
NDMaterial::setTrialStrainIncr(const Vector &v)
{
   g3ErrorHandler->fatal("NDMaterial::setTrialStrainIncr -- subclass responsibility");
   return 0;    
}

int 
NDMaterial::setTrialStrainIncr(const Vector &v, const Vector &r)
{
   g3ErrorHandler->fatal("NDMaterial::setTrialStrainIncr -- subclass responsibility");
   return 0;    
}

const Matrix &
NDMaterial::getTangent(void)
{
   g3ErrorHandler->fatal("NDMaterial::getTangent -- subclass responsibility");
   return errMatrix;    
}

const Vector &
NDMaterial::getStress(void)
{
   g3ErrorHandler->fatal("NDMaterial::getStress -- subclass responsibility");
   return errVector;    
}

const Vector &
NDMaterial::getStrain(void)
{
   g3ErrorHandler->fatal("NDMaterial::getStrain -- subclass responsibility");
   return errVector;    
}

int 
NDMaterial::setTrialStrain(const Tensor &v)
{
   g3ErrorHandler->fatal("NDMaterial::setTrialStrainIncr -- subclass responsibility");
   return 0;    
}

int 
NDMaterial::setTrialStrain(const Tensor &v, const Tensor &r)    
{
   g3ErrorHandler->fatal("NDMaterial::setTrialStrainIncr -- subclass responsibility");
   return 0;    
}

int 
NDMaterial::setTrialStrainIncr(const Tensor &v)
{
   g3ErrorHandler->fatal("NDMaterial::setTrialStrainIncr -- subclass responsibility");
   return 0;    
}

int 
NDMaterial::setTrialStrainIncr(const Tensor &v, const Tensor &r)
{
   g3ErrorHandler->fatal("NDMaterial::setTrialStrainIncr -- subclass responsibility");
   return 0;    
}

const Tensor &
NDMaterial::getTangentTensor(void)
{
   g3ErrorHandler->fatal("NDMaterial::getTangentTensor -- subclass responsibility");
   return errTensor;    
}

const Tensor &
NDMaterial::getStressTensor(void)
{
   g3ErrorHandler->fatal("NDMaterial::getStressTensor -- subclass responsibility");
   return errTensor;    
}

const Tensor &
NDMaterial::getStrainTensor(void)
{
   g3ErrorHandler->fatal("NDMaterial::getStrainTensor -- subclass responsibility");
   return errTensor;    
}

Response*
NDMaterial::setResponse (char **argv, int argc, Information &matInfo)
{
    if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());

    else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
    
	else if (strcmp(argv[0],"tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());
    
	else
		return 0;
}

int 
NDMaterial::getResponse (int responseID, Information &matInfo)
{
	switch (responseID) {
		case 1:
			return matInfo.setVector(this->getStress());

		case 2:
			return matInfo.setVector(this->getStrain());

		case 3:
			return matInfo.setMatrix(this->getTangent());
			
		default:
			return -1;
	}
}
