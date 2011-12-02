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
**                                                                    **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.26 $                                                              
// $Date: 2010-09-13 21:29:28 $                                                                  
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
#include <OPS_Globals.h>
#include <Matrix.h>
#include <Vector.h>
#include <MaterialResponse.h>

#include <PlaneStressMaterial.h>
#include <BeamFiberMaterial.h>
#include <PlateFiberMaterial.h>
#include <string.h>

Matrix NDMaterial::errMatrix(1,1);
Vector NDMaterial::errVector(1);

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

NDMaterial*
NDMaterial::getCopy(const char *type)
{
  if (strcmp(type,"PlaneStress") == 0 ||
      strcmp(type,"PlaneStress2D") == 0) {
    NDMaterial *copy = this->getCopy("ThreeDimensional");
    PlaneStressMaterial *clone = new PlaneStressMaterial(this->getTag(),*copy);
    return clone;
  }
  else if (strcmp(type,"BeamFiber") == 0 ||
	   strcmp(type,"TimoshenkoFiber") == 0) {
    NDMaterial *copy = this->getCopy("ThreeDimensional");
    BeamFiberMaterial *clone = new BeamFiberMaterial(this->getTag(),*copy);
    return clone;
  }
  else if (strcmp(type,"PlateFiber") == 0) {
    NDMaterial *copy = this->getCopy("ThreeDimensional");
    PlateFiberMaterial *clone = new PlateFiberMaterial(this->getTag(),*copy);
    return clone;
  }
  else
    return 0;
}

double
NDMaterial::getRho(void)
{
  return 0.0;
}

// methods to set and retrieve state.
int 
NDMaterial::setTrialStrain(const Vector &v)
{
   opserr << "NDMaterial::setTrialStrain -- subclass responsibility\n";
   return -1;    
}

int 
NDMaterial::setTrialStrain(const Vector &v, const Vector &r)
{
   opserr << "NDMaterial::setTrialStrain -- subclass responsibility\n";
   return -1;    
}

int 
NDMaterial::setTrialStrainIncr(const Vector &v)
{
   opserr << "NDMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;    
}

int 
NDMaterial::setTrialStrainIncr(const Vector &v, const Vector &r)
{
   opserr << "NDMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;    
}

const Matrix &
NDMaterial::getTangent(void)
{
   opserr << "NDMaterial::getTangent -- subclass responsibility\n";
   return errMatrix;    
}

const Vector &
NDMaterial::getStress(void)
{
   opserr << "NDMaterial::getStress -- subclass responsibility\n";
   return errVector;    
}

const Vector &
NDMaterial::getStrain(void)
{
   opserr << "NDMaterial::getStrain -- subclass responsibility\n";
   return errVector;    
}

Response*
NDMaterial::setResponse (const char **argv, int argc, 
			 OPS_Stream &output)
{
  Response *theResponse =0;
  const char *matType = this->getType();

  output.tag("NdMaterialOutput");
  output.attr("matType",this->getClassType());
  output.attr("matTag",this->getTag());

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
    const Vector &res = this->getStress();
    int size = res.Size();
    
    if ( (strcmp(matType,"PlaneStress") == 0 && size == 3) ||
	 (strcmp(matType,"PlaneStrain") == 0 && size == 3)) {
	output.tag("ResponseType","sigma11");
	output.tag("ResponseType","sigma22");
	output.tag("ResponseType","sigma12");
    } else if (strcmp(matType,"ThreeDimensional") == 0 && size == 6) {
	output.tag("ResponseType","sigma11");
	output.tag("ResponseType","sigma22");
	output.tag("ResponseType","sigma33");
	output.tag("ResponseType","sigma12");
	output.tag("ResponseType","sigma13");
	output.tag("ResponseType","sigma23");
    } else {
      for (int i=0; i<size; i++) 
	output.tag("ResponseType","UnknownStress");
    }
    theResponse =  new MaterialResponse(this, 1, this->getStress());

  } else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
    const Vector &res = this->getStrain();
    int size = res.Size();
    if ( (strcmp(matType,"PlaneStress") == 0 && size == 3) ||
	 (strcmp(matType,"PlaneStrain") == 0 && size == 3)) {
	output.tag("ResponseType","eta11");
	output.tag("ResponseType","eta22");
	output.tag("ResponseType","eta12");
    } else if (strcmp(matType,"ThreeDimensional") == 0 && size == 6) {
	output.tag("ResponseType","eps11");
	output.tag("ResponseType","eps22");
	output.tag("ResponseType","eps33");
	output.tag("ResponseType","eps12");
	output.tag("ResponseType","eps13");
	output.tag("ResponseType","eps23");
    } else {
      for (int i=0; i<size; i++) 
	output.tag("ResponseType","UnknownStrain");
    }      
    theResponse =  new MaterialResponse(this, 2, this->getStress());
  }

  output.endTag(); // NdMaterialOutput

  return theResponse;
}

int 
NDMaterial::getResponse (int responseID, Information &matInfo)
{
  switch (responseID) {
  case 1:
    return matInfo.setVector(this->getStress());
    
  case 2:
    return matInfo.setVector(this->getStrain());
    
  default:
    return -1;
  }
}



// AddingSensitivity:BEGIN ////////////////////////////////////////
const Vector &
NDMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
	static Vector dummy(1);
	return dummy;
}

const Vector &
NDMaterial::getStrainSensitivity(int gradIndex)
{
	static Vector dummy(1);
	return dummy;
}

double
NDMaterial::getRhoSensitivity(int gradIndex)
{
	return 0.0;
}

const Matrix &
NDMaterial::getDampTangentSensitivity(int gradIndex)
{
	static Matrix dummy(1,1);
	return dummy;
}
const Matrix &
NDMaterial::getTangentSensitivity(int gradIndex)
{
	static Matrix dummy(1,1);
	return dummy;
}
const Matrix &
NDMaterial::getInitialTangentSensitivity(int gradIndex)
{
	static Matrix dummy(1,1);
	return dummy;
}
int
NDMaterial::commitSensitivity(Vector & strainSensitivity, int gradIndex, int numGrads)
{
	return 0;
}
// AddingSensitivity:END //////////////////////////////////////////


