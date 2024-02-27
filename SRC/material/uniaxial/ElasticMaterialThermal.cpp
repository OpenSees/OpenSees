
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
                                                                        
// $Revision: 1.10 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticMaterialThermal.cpp,v $
                                                                        
//
// Description: This file contains the class implementation for 
// ElasticMaterialThermal. 
//
// What: "@(#) ElasticMaterialThermal.C, revA"
// Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 


#include <ElasticMaterialThermal.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>

#include <OPS_Globals.h>

#include <elementAPI.h>

void *
OPS_ElasticMaterialThermal(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  int softindex = 0;
  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "Invalid #args,  want: uniaxialMaterial Elastic tag? E? alpha?<eta?> ... " << endln;
    return 0;
  }
  
  int iData[1];
  double dData1[2];
  double dData2[2];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial Elastic" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData == 1)
	  dData1[1] = 0.0;
  else
	  numData = 2;

  if (OPS_GetDoubleInput(&numData, dData1) != 0) {
	  opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
	  return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();
  if (numData >0) {
	  const char* typeChar = OPS_GetString();
	  if ((strcmp(typeChar, "-SteelSoft") == 0) || (strcmp(typeChar, "-SSoft") == 0) || (strcmp(typeChar, "-sSoft") == 0)) {
		  softindex = 1;
	  }
	  else if ((strcmp(typeChar, "-ConcreteSoft") == 0) || (strcmp(typeChar, "-CSoft") == 0) || (strcmp(typeChar, "-cSoft") == 0)) {
		  softindex = 2;
	  }
  }

  dData2[0] = 0.0; dData2[1] = 0.0;

	numData = numData-1;
	if (numData > 2)
		numData = 2;
	if(numData>0){
	if (OPS_GetDoubleInput(&numData, dData2) != 0) {
      opserr << "Invalid data for uniaxial Elastic " << iData[0] << endln;
      return 0;
	}
	}
 //opserr<< "received alpha"<<dData[1]<<endln;

  // Parsing was successful, allocate the material
  theMaterial = new ElasticMaterialThermal(iData[0], dData1[0], dData1[1], dData2[0], dData2[1], softindex);
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ElasticMaterialThermal\n";
    return 0;
  }

  return theMaterial;
}


double ElasticMaterialThermal::SteelRedFactors[12] = { 1.0, 0.9, 0.8 ,0.7, 0.6 ,0.31, 0.13, 0.09, 0.0675, 0.045, 0.0225 , 0.0 };
double ElasticMaterialThermal::ConcRedFactors[12] = { 0.625, 0.4318 ,0.3036, 0.1875 ,0.1, 0.045, 0.03, 0.015, 0.008 , 0.004,0.001,0.0 };

ElasticMaterialThermal::ElasticMaterialThermal(int tag, double e, double alpha, double et, double eneg, int softindex)
:UniaxialMaterial(tag,MAT_TAG_ElasticMaterialThermal),
 trialStrain(0.0), trialStrainRate(0.0), committedStrain(0.0), committedStrainRate(0.0),
 Epos(e), eta(et),
 Alpha(alpha), E0(e), softIndex(softindex),
 ThermalElongation(0),Temp(0), parameterID(0)
{
	if (eneg <1e-10) {
		Eneg0 = e;
		Eneg = e;
	}
		
	else {
		Eneg0 = eneg;
		Eneg = eneg;
	}

	if (softIndex < 0 || softIndex > 2) {
	  opserr << "ElasticMaterialThermal " << this->getTag() << " receives an invalid softening index, setting softening index to 0" << endln;
	  softIndex = 0;
	}
}

ElasticMaterialThermal::ElasticMaterialThermal()
:UniaxialMaterial(0,MAT_TAG_ElasticMaterialThermal),
 trialStrain(0.0),  trialStrainRate(0.0), committedStrain(0.0), committedStrainRate(0.0),
 Epos(0.0),Eneg(0.0), eta(0.0),
 Alpha(0.0), E0(0.0), Eneg0(0.0), softIndex(0),
 ThermalElongation(0), Temp(0), parameterID(0)
{
	
}

ElasticMaterialThermal::~ElasticMaterialThermal()
{
	// does nothing
}


int
ElasticMaterialThermal::setTrialStrain(double strain, double strainRate)
{
	trialStrain = strain;
	trialStrainRate = strainRate;
	return 0;
}


int
ElasticMaterialThermal::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
	trialStrain = strain;
	trialStrainRate = strainRate;

	if (trialStrain >= 0.0) {
		stress = Epos*trialStrain + eta*trialStrainRate;
		tangent = Epos;
	}
	else {
		stress = Eneg*trialStrain + eta*trialStrainRate;
		tangent = Eneg;
	}

	return 0;
}


double
ElasticMaterialThermal::getStress(void)
{
	if (trialStrain >= 0.0)
		return Epos*trialStrain + eta*trialStrainRate;
	else
		return Eneg*trialStrain + eta*trialStrainRate;
}


double
ElasticMaterialThermal::getTangent(void)
{
	if (trialStrain > 0.0)
		return Epos;
	else if (trialStrain < 0.0)
		return Eneg;
	else
		return (Epos > Eneg) ? Epos : Eneg;
}


double
ElasticMaterialThermal::getInitialTangent(void)
{
	return (Epos > Eneg) ? Epos : Eneg;
}


int
ElasticMaterialThermal::commitState(void)
{
	committedStrain = trialStrain;
	committedStrainRate = trialStrainRate;
	return 0;
}


int
ElasticMaterialThermal::revertToLastCommit(void)
{
	trialStrain = committedStrain;
	trialStrainRate = committedStrainRate;
	return 0;
}


int
ElasticMaterialThermal::revertToStart(void)
{
	trialStrain = 0.0;
	trialStrainRate = 0.0;
	return 0;
}


UniaxialMaterial *
ElasticMaterialThermal::getCopy(void)
{
	ElasticMaterialThermal *theCopy = new ElasticMaterialThermal(this->getTag(), Epos, Alpha, eta, Eneg, softIndex);
	theCopy->trialStrain = trialStrain;
	theCopy->trialStrainRate = trialStrainRate;
	theCopy->committedStrain = committedStrain;
	theCopy->committedStrainRate = committedStrainRate;
	return theCopy;
}


int
ElasticMaterialThermal::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	static Vector data(6);
	data(0) = this->getTag();
	data(1) = Epos;
	data(2) = Eneg;
	data(3) = eta;
	data(4) = committedStrain;
	data(5) = committedStrainRate;
	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0)
		opserr << "ElasticMaterialThermal::sendSelf() - failed to send data\n";

	return res;
}


int
ElasticMaterialThermal::recvSelf(int cTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(6);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);

	if (res < 0) {
		opserr << "ElasticMaterialThermal::recvSelf() - failed to receive data\n";
		Epos = Eneg = 0;
		this->setTag(0);
	}
	else {
		this->setTag(int(data(0)));
		Epos = data(1);
		Eneg = data(2);
		eta = data(3);
		committedStrain = data(4);
		committedStrainRate = data(5);
		this->revertToLastCommit();
	}

	return res;
}


void
ElasticMaterialThermal::Print(OPS_Stream &s, int flag)
{
	s << "Elastic tag: " << this->getTag() << endln;
	s << "  Epos: " << Epos << " Eneg: " << Eneg << " eta: " << eta << endln;
}


int
ElasticMaterialThermal::setParameter(const char **argv, int argc, Parameter &param)
{

	if (strcmp(argv[0], "E") == 0) {
		param.setValue(Epos);
		return param.addObject(1, this);
	}
	if (strcmp(argv[0], "Epos") == 0) {
		param.setValue(Epos);
		return param.addObject(2, this);
	}
	if (strcmp(argv[0], "Eneg") == 0) {
		param.setValue(Eneg);
		return param.addObject(3, this);
	}
	else if (strcmp(argv[0], "eta") == 0) {
		param.setValue(eta);
		return param.addObject(4, this);
	}
	return -1;
}


int
ElasticMaterialThermal::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		Epos = info.theDouble;
		Eneg = info.theDouble;
		return 0;
	case 2:
		Epos = info.theDouble;
		return 0;
	case 3:
		Eneg = info.theDouble;
		return 0;
	case 4:
		eta = info.theDouble;
		return 0;
	default:
		return -1;
	}
}


int
ElasticMaterialThermal::activateParameter(int paramID)
{
	parameterID = paramID;

	return 0;
}


double
ElasticMaterialThermal::getStressSensitivity(int gradIndex, bool conditional)
{
	if (parameterID == 1)
		return trialStrain;
	if (parameterID == 2 && trialStrain > 0.0)
		return trialStrain;
	if (parameterID == 3 && trialStrain < 0.0)
		return trialStrain;
	if (parameterID == 4)
		return trialStrainRate;

	return 0.0;
}


double
ElasticMaterialThermal::getTangentSensitivity(int gradIndex)
{
	if (parameterID == 1)
		return 1.0;
	if (parameterID == 2 && trialStrain >= 0.0)
		return 1.0;
	if (parameterID == 3 && trialStrain <= 0.0)
		return 1.0;

	return 0.0;
}


double
ElasticMaterialThermal::getInitialTangentSensitivity(int gradIndex)
{
	if (parameterID == 1)
		return 1.0;
	if (parameterID == 2)
		return 1.0;
	if (parameterID == 3)
		return 1.0;

	return 0.0;
}


int
ElasticMaterialThermal::commitSensitivity(double strainGradient,
	int gradIndex, int numGrads)
{
	// Nothing to commit ... path independent
	return 0;
}



//Liming for updating the reduction factors////////////start
double
ElasticMaterialThermal::getElongTangent(double TempT, double& ET, double& Elong, double TempTmax)
{
  double ThermalElongation; // Why in the world is this defined as a local variable using
  // the same name as a private class variable??? MHS
  ThermalElongation = 0.0; // Let me go ahead and initialize it too....
  
  Temp = TempT;
  if (softIndex != 0) {

    double *redfactors = 0;
    if (softIndex == 1)
      redfactors = SteelRedFactors;
    if (softIndex == 2)
      redfactors = ConcRedFactors;
    
    for (int i = 0; i < 13; i++) {
      if (Temp <= 80 + 100 * i)
	{
	  if (i == 0) {
	    Epos = E0 * (1.0 - Temp * (1.0 - redfactors[0]) / 80);
	    Eneg = Eneg0 * (1.0 - Temp * (1.0 - redfactors[0]) / 80);
	  }
	  else if (i == 12) {
	    opserr << "Warning:The temperature " << Temp << " for SteelECthermal is out of range\n";
	    return -1;
	  }
	  else {
	    Epos = E0 * (redfactors[i - 1] - (Temp + 20 - 100 * i) * (redfactors[i - 1] - redfactors[i]) / 100);
	    Eneg = Eneg0 * (redfactors[i - 1] - (Temp + 20 - 100 * i) * (redfactors[i - 1] - redfactors[i]) / 100);
	  }
	  break;
	}
      
    }
    
    if (softIndex == 1) {
      if (Temp <= 1) {
	ThermalElongation = Temp * 1.2164e-5;
      }
      else if (Temp <= 730) {
	ThermalElongation = -2.416e-4 + 1.2e-5 * (Temp + 20) + 0.4e-8 * (Temp + 20) * (Temp + 20);
      }
      else if (Temp <= 840) {
	ThermalElongation = 11e-3;
      }
      else if (Temp <= 1180) {
	ThermalElongation = -6.2e-3 + 2e-5 * (Temp + 20);
      }
    }
    else if (softIndex == 2) {
      if (Temp <= 1) {
	ThermalElongation = Temp * 9.213e-6;
      }
      else if (Temp <= 680) {
	ThermalElongation = -1.8e-4 + 9e-6 * (Temp + 20) + 2.3e-11 * (Temp + 20) * (Temp + 20) * (Temp + 20);
      }
      else if (Temp <= 1180) {
	ThermalElongation = 14e-3;
      }
      
    }
    Elong = ThermalElongation;
    
  }
  else {
    ET = E0;
    ThermalElongation = Alpha * TempT;
    Elong = ThermalElongation;
  }
  
  return 0.0; // We have to return something, why is this function not void? ... MHS
}

double
ElasticMaterialThermal::getThermalElongation(void)
{
	return ThermalElongation;
}

int
ElasticMaterialThermal::getVariable(const char *variable, Information &info)
{
	if (strcmp(variable, "ThermalElongation") == 0) {
		info.theDouble = ThermalElongation;
		return 0;
	}
	else if (strcmp(variable, "ElongTangent") == 0) {
		Vector *theVector = info.theVector;
		if (theVector != 0) {
			double tempT, ET, Elong, TempTmax;
			tempT = (*theVector)(0);
			ET = (*theVector)(1);
			Elong = (*theVector)(2);
			TempTmax = (*theVector)(3);
			this->getElongTangent(tempT, ET, Elong, TempTmax);
			(*theVector)(0) = tempT;
			(*theVector)(1) = ET;
			(*theVector)(2) = Elong;
			(*theVector)(3) = TempTmax;
		}
		return 0;
	}
	else if (strcmp(variable, "TempAndElong") == 0) {
		Vector *theVector = info.theVector;
		if (theVector != 0) {
			(*theVector)(0) = Temp;
			(*theVector)(1) = ThermalElongation;
		}
		else {
			opserr << "null Vector in EC" << endln;
		}

		return 0;
	}
	return -1;
}
