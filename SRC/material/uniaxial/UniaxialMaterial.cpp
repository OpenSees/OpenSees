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
                                                                        
// $Revision: 1.14 $
// $Date: 2006-08-04 18:17:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/UniaxialMaterial.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 05/98
// Revision: A
//
// Description: This file contains the class implementation for 
// UniaxialMaterial.
//
// What: "@(#) UniaxialMaterial.C, revA"

#include <UniaxialMaterial.h>
#include <string.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <float.h>
#include <Vector.h>
#include <DataOutputHandler.h>


UniaxialMaterial::UniaxialMaterial(int tag, int clasTag)
:Material(tag,clasTag)
{

}

UniaxialMaterial::~UniaxialMaterial()
{
	// does nothing
}

int
UniaxialMaterial::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
  int res = this->setTrialStrain(strain, strainRate);
  if (res == 0) {
    stress = this->getStress();
    tangent = this->getTangent();
  } else {
    opserr << "UniaxialMaterial::setTrial() - material failed in setTrialStrain()\n"; 
  }

  return res;
}

// default operation for strain rate is zero
double
UniaxialMaterial::getStrainRate(void)
{
    return 0.0;
}



// default operation for damping tangent is zero
double
UniaxialMaterial::getDampTangent(void)
{
    return 0.0;
}

// default operation for secant stiffness
double
UniaxialMaterial::getSecant (void)
{
	double strain = this->getStrain();
	double stress = this->getStress();

	if (strain != 0.0)
		return stress/strain;
	else
		return this->getTangent();
}

double 
UniaxialMaterial::getRho(void)
{
	return 0.0;
}

UniaxialMaterial*
UniaxialMaterial::getCopy(SectionForceDeformation *s)
{
	return this->getCopy();
}

Response* 
UniaxialMaterial::setResponse(const char **argv, int argc, Information &matInfo, OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  theOutput.tag("UniaxialMaterialOutput");
  theOutput.attr("matType", this->getClassType());
  theOutput.attr("matTag", this->getTag());

  // stress
  if (strcmp(argv[0],"stress") == 0) {
    theOutput.tag("ResponseType", "sigma11");
    theResponse =  new MaterialResponse(this, 1, this->getStress());
  }  
  // tangent
  else if (strcmp(argv[0],"tangent") == 0) {
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 2, this->getTangent());
  }

  // strain
  else if (strcmp(argv[0],"strain") == 0) {
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 3, this->getStrain());
  }

  // strain
  else if ((strcmp(argv[0],"stressStrain") == 0) || 
	   (strcmp(argv[0],"stressANDstrain") == 0)) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 4, Vector(2));
  }

  theOutput.endTag();
  return theResponse;

}
 
int 
UniaxialMaterial::getResponse(int responseID, Information &matInfo)
{
  static Vector stressStrain(2);
  // each subclass must implement its own stuff    
  switch (responseID) {
    case 1:
      matInfo.setDouble(this->getStress());
      return 0;
      
    case 2:
      matInfo.setDouble(this->getTangent());
      return 0;      

    case 3:
      matInfo.setDouble(this->getStrain());
      return 0;      
    
    case 4:
      stressStrain(0) = this->getStress();
      stressStrain(1) = this->getStrain();
      matInfo.setVector(stressStrain);
      return 0;
      
  default:      
    return -1;
  }
}


// AddingSensitivity:BEGIN ////////////////////////////////////////
int
UniaxialMaterial::setParameter(const char **argv, int argc, Information &info)
{
    return -1;
}

int
UniaxialMaterial::updateParameter(int parameterID, Information &info)
{
    return -1;
}

int
UniaxialMaterial::activateParameter(int parameterID)
{
    return -1;
}

double
UniaxialMaterial::getStressSensitivity(int gradNumber, bool conditional)
{
    return 0.0;
}

double
UniaxialMaterial::getStrainSensitivity(int gradNumber)
{
    return 0.0;
}

double
UniaxialMaterial::getInitialTangentSensitivity(int gradNumber)
{
    return 0.0;
}

double
UniaxialMaterial::getRhoSensitivity(int gradNumber)
{
    return 0.0;
}

double
UniaxialMaterial::getDampTangentSensitivity(int gradNumber)
{
    return 0.0;
}

int
UniaxialMaterial::commitSensitivity(double strainSensitivity, int gradNumber, int numGrads)
{
    return -1;
}

double
UniaxialMaterial::getInitialTangent (void)
{
	opserr << "UniaxialMaterial::getInitialTangent() -- this mehtod " << endln
		<< " is not implemented for the selected material. " << endln;
	return 0.0;
}

// AddingSensitivity:END //////////////////////////////////////////
