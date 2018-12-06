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
                                                                        
// $Revision: 1.22 $
// $Date: 2009-08-25 23:40:17 $
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
#include <stdlib.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>

static MapOfTaggedObjects theUniaxialMaterialObjects;

bool OPS_addUniaxialMaterial(UniaxialMaterial *newComponent) {
  return theUniaxialMaterialObjects.addComponent(newComponent);
}

bool OPS_removeUniaxialMaterial(int tag)
{
    TaggedObject* obj = theUniaxialMaterialObjects.removeComponent(tag);
    if (obj != 0) {
	delete obj;
	return true;
    }
    return false;
}

UniaxialMaterial *OPS_getUniaxialMaterial(int tag) {

  TaggedObject *theResult = theUniaxialMaterialObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "UniaxialMaterial *getUniaxialMaterial(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  UniaxialMaterial *theMat = (UniaxialMaterial *)theResult;

  return theMat;
}

void OPS_clearAllUniaxialMaterial(void) {
  theUniaxialMaterialObjects.clearAll();
}

void OPS_printUniaxialMaterial(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\"uniaxialMaterials\": [\n";        
    MapOfTaggedObjectsIter theObjects = theUniaxialMaterialObjects.getIter();
    theObjects.reset();
    TaggedObject *theObject;
    int count = 0;
    int numComponents = theUniaxialMaterialObjects.getNumComponents();    
    while ((theObject = theObjects()) != 0) {
      UniaxialMaterial *theMaterial = (UniaxialMaterial *)theObject;
      theMaterial->Print(s, flag);
      if (count < numComponents-1)
	s << ",\n";
      count++;      
    }
    s << "\n\t\t]";
  }
}

UniaxialMaterial::UniaxialMaterial(int tag, int clasTag)
:Material(tag,clasTag)
{

}


UniaxialMaterial::UniaxialMaterial()
    :Material(0, 0)
{

}


UniaxialMaterial::~UniaxialMaterial()
{
	// does nothing
}


int
UniaxialMaterial::setTrialStrain(double strain, double temperature, double strainRate)
{
  int res = this->setTrialStrain(strain, strainRate);

  return res;
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


int
UniaxialMaterial::setTrial(double strain, double temperature, double &stress, double &tangent, double &thermalElongation, double strainRate)
{
  int res = this->setTrialStrain(strain, temperature, strainRate);

  if (res == 0) {
    static const char thermal[] = "ThermalElongation";
    const char *thermalPointer = thermal;


    Information info;
    stress = this->getStress();
    tangent = this->getTangent();
    this->getVariable(thermalPointer, info);
    thermalElongation = info.theDouble;

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
/*
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
*/

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
UniaxialMaterial::setResponse(const char **argv, int argc,
			      OPS_Stream &theOutput)
{
  Response *theResponse = 0;

  if ( (strcmp(argv[0],"stress") == 0) ||
       (strcmp(argv[0],"tangent") == 0)||
       (strcmp(argv[0],"strain") == 0) ||
       (strcmp(argv[0],"stressStrain") == 0) ||
       (strcmp(argv[0],"stressANDstrain") == 0) || 
       (strcmp(argv[0],"stressAndStrain") == 0) ||
       (strcmp(argv[0],"stressStrain") == 0) || 
       (strcmp(argv[0],"stressANDstrain") == 0) ||
       (strcmp(argv[0],"stressAndStrain") == 0) ||
       (strcmp(argv[0],"stressStrainTangent") == 0) || 
       (strcmp(argv[0],"stressANDstrainANDtangent") == 0) ||
       (strstr(argv[0],"stressSensitivity") != 0) ||
       (strstr(argv[0],"strainSensitivity") != 0)||
	  (strstr(argv[0], "TempElong") != 0)) {
    
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

    else if (strcmp(argv[0],"plasticStrain") == 0) {
      theOutput.tag("ResponseType", "eps11");
      theResponse =  new MaterialResponse(this, 6, this->getStrain());
    }
    
    // strain
    else if ((strcmp(argv[0],"stressStrain") == 0) || 
	     (strcmp(argv[0],"stressANDstrain") == 0) ||
	     (strcmp(argv[0],"stressAndStrain") == 0)) {
      theOutput.tag("ResponseType", "sig11");
      theOutput.tag("ResponseType", "eps11");
      theResponse =  new MaterialResponse(this, 4, Vector(2));
    }
    
    else if ((strcmp(argv[0],"stressStrainTangent") == 0) || 
	     (strcmp(argv[0],"stressANDstrainANDtangent") == 0)) {
      theOutput.tag("ResponseType", "sig11");
      theOutput.tag("ResponseType", "eps11");
      theOutput.tag("ResponseType", "C11");
      theResponse =  new MaterialResponse(this, 5, Vector(3));
    }

    // stress sensitivity for local sensitivity recorder purpose.  Quan 2009
    // limit:  no more than 10000 random variables/sensitivity parameters
    else if (strstr(argv[0],"stressSensitivity") != 0) {
      char *token = strtok((char *) argv[0], " ");
      if (token != NULL) token = strtok(NULL, " ");
      int gradient = atoi(token);
      theOutput.tag("ResponseType", "sigsens11");
      theResponse =  new MaterialResponse(this, gradient+10000, this->getStress());
    }
    // strain sensivitiy
    else if (strstr(argv[0],"strainSensitivity") != 0) {
      char *token = strtok((char *) argv[0], " ");
      if (token != NULL) token = strtok(NULL, " ");
      int gradient = atoi(token);
      theOutput.tag("ResponseType", "epssens11");
      theResponse =  new MaterialResponse(this, gradient+20000, this->getStrain());
    }
	//Added by Liming, UoE, for temperature and elongation output,[SIF]2017
	else if ((strcmp(argv[0], "TempElong") == 0) ||
		(strcmp(argv[0], "tempANDelong") == 0)) {
		theOutput.tag("ResponseType", "temp11");
		theOutput.tag("ResponseType", "Elong11");
		theResponse = new MaterialResponse(this, 7, Vector(2));
	}
    
    theOutput.endTag();
  }
  
  return theResponse;

}
 
int 
UniaxialMaterial::getResponse(int responseID, Information &matInfo)
{
  static Vector stressStrain(2);
  static Vector stressStrainTangent(3);

  static Vector tempData(2);  //L.jiang [SIF]
  static Information infoData(tempData);  //L.jiang [SIF]

  // each subclass must implement its own stuff   

  // added for sensitivity recorder. Quan 2009
  if ((responseID>10000)&&(responseID<20000)){
      matInfo.setDouble(this->getStressSensitivity(responseID-10000,false));
      return 0;
  }
  else if (responseID>20000){
      matInfo.setDouble(this->getStrainSensitivity(responseID-20000));
      return 0;
  }

  double kInit;
  double stress;
  double strain;

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

  case 6: // an approx to plastic strain
      strain = this->getStrain();
      stress = this->getStress();
      kInit = this->getTangent();
      strain = strain-stress/kInit;
      matInfo.setDouble(strain);
      return 0;      
    
    case 4:
      stressStrain(0) = this->getStress();
      stressStrain(1) = this->getStrain();
      matInfo.setVector(stressStrain);
      return 0;
      
	  case 5:
      stressStrainTangent(0) = this->getStress();
      stressStrainTangent(1) = this->getStrain();
      stressStrainTangent(2) = this->getTangent();
      matInfo.setVector(stressStrainTangent);
      return 0;
	 
	  //Added by Liming, UoE, for temperature and elongation output,[SIF]2017
	  case 7:
		  if ((this->getVariable("TempAndElong", infoData)) != 0) {
			  opserr << "Warning: invalid tag in uniaxialMaterial:getVariable" << endln;
			  return -1;
		  }
		  tempData = infoData.getData();
		  matInfo.setVector(tempData);
		  return 0;
  default:      
    return -1;
  }
}


// AddingSensitivity:BEGIN ////////////////////////////////////////
double
UniaxialMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
    return 0.0;
}

double
UniaxialMaterial::getStrainSensitivity(int gradIndex)
{
    return 0.0;
}

double
UniaxialMaterial::getTangentSensitivity(int gradIndex)
{
    return 0.0;
}

double
UniaxialMaterial::getInitialTangentSensitivity(int gradIndex)
{
    return 0.0;
}

double
UniaxialMaterial::getRhoSensitivity(int gradIndex)
{
    return 0.0;
}

double
UniaxialMaterial::getDampTangentSensitivity(int gradIndex)
{
    return 0.0;
}

int
UniaxialMaterial::commitSensitivity(double strainSensitivity, int gradIndex, int numGrads)
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
