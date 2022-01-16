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
                                                                        
// $Revision: 1.11 $
// $Date: 2010-06-01 23:35:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/EPPGapMaterial.cpp,v $

// Written: krm 
// Created: 07/2000
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) EPPGapMaterial.C, revA"

#include <stdlib.h>

#include <EPPGapMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>
#include <OPS_Globals.h>

void * OPS_ADD_RUNTIME_VPV(OPS_EPPGapMaterial)
{
    // Pointer to a uniaxial material that will be returned
    UniaxialMaterial *theMaterial = 0;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 4) {
	opserr << "Invalid #args,  want: uniaxialMaterial ElasticPPGap tag E Fy gap <eta damage>\n";
	return 0;
    }
  
    int tag, damage = 0;
    double dData[4];
    dData[3] = 0.0; // setting default eta to 0.

    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
	opserr << "WARNING invalid tag for uniaxialMaterial EPPGap" << endln;
	return 0;
    }

    numData = OPS_GetNumRemainingInputArgs();
    if(numData > 4) numData = 4;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
	opserr << "Invalid data for uniaxial EPPGap \n";
	return 0;	
    }

    numData = OPS_GetNumRemainingInputArgs();
    if(numData > 0) {
	numData = 1;
	const char* damagestr = OPS_GetString();
	if (strcmp(damagestr, "damage") == 0 || strcmp(damagestr, "Damage") == 0) {
	    damage = 1;
	}
    }

    // Parsing was successful, allocate the material
    theMaterial = new EPPGapMaterial(tag, dData[0], dData[1], dData[2], dData[3], damage);
    if (theMaterial == 0) {
	opserr << "WARNING could not create uniaxialMaterial of type EPPGap\n";
	return 0;
    }

    return theMaterial;
}

EPPGapMaterial::EPPGapMaterial(int tag, double e, double fyl, double gap0, double eta0, int accum)
:UniaxialMaterial(tag,MAT_TAG_EPPGap),
    E(e), fy(fyl), gap(gap0), eta(eta0), minElasticYieldStrain(gap0), damage(accum), 
    parameterID(0), SHVs(0), EnergyP(0.0), trialStrain(0.0), commitStrain(0.0)
{
	if (E == 0.0) {
	  opserr << "EPPGapMaterial::EPPGapMaterial -- E is zero, continuing with E = fy/0.002\n";
	  if (fy != 0.0)
	    E = fabs(fy)/0.002;
	  else {
	    opserr << "EPPGapMaterial::EPPGapMaterial -- E and fy are zero\n";
	    exit(-1);
	  }
	}

	if (fy*gap<0) {
	    opserr << "EPPGapMaterial::EPPGapMaterial -- Alternate signs on fy and gap encountered, continuing anyway\n";
	}
        
    if (eta >= 1.0) {
        opserr << "EPPGapMaterial::EPPGapMaterial -- value of eta must be < 1, setting eta to 0\n";
        eta = 0;
    }

    maxElasticYieldStrain = fy / E + gap;
    //use setTrialStrain method to get trial stress and tangent. Save as committed. Added by ambaker1
    this->setTrialStrain(trialStrain);
    commitStress = trialStress;
    commitTangent = trialTangent;
}

EPPGapMaterial::EPPGapMaterial()
:UniaxialMaterial(0,MAT_TAG_EPPGap),
    E(0.0), fy(0.0), gap(0.0), eta(0.0), minElasticYieldStrain(0.0), 
    maxElasticYieldStrain(0.0), damage(0), parameterID(0), SHVs(0), EnergyP(0.0),
    trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
    commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{

}

EPPGapMaterial::~EPPGapMaterial()
{
  // does nothing
}

int 
EPPGapMaterial::setTrialStrain(double strain, double strainRate)
{
  // set the trial strain
  trialStrain = strain;

  // determine trial stress and tangent
  if (fy >= 0) {
    if (trialStrain > maxElasticYieldStrain) {
      trialStress = fy+(trialStrain-gap-fy/E)*eta*E;
      trialTangent = eta*E;
    } else if (trialStrain < minElasticYieldStrain) {
      trialStress = 0;
      trialTangent = 0;
    } else {
      trialStress =  E*(trialStrain-minElasticYieldStrain);
      trialTangent = E;
    }
  } else {
    if (trialStrain < maxElasticYieldStrain) {
      trialStress =  fy+(trialStrain-gap-fy/E)*eta*E;
      trialTangent = eta*E;
    } else if (trialStrain > minElasticYieldStrain) {
      trialStress =  0;
      trialTangent = 0;
    } else {
      trialStress =  E*(trialStrain-minElasticYieldStrain);
      trialTangent = E;
    }
  }

  return 0;
}

double 
EPPGapMaterial::getStrain(void)
{
    return trialStrain;
}

double 
EPPGapMaterial::getStress(void)
{
  return trialStress;

}

double 
EPPGapMaterial::getTangent(void)
{
  return trialTangent;
}

double 
EPPGapMaterial::getInitialTangent(void)
{
    if ((fy >= 0.0 && gap > 0.0) || (fy < 0.0 && gap < 0.0))
        return 0.0;
    else 
        return E;
}

int 
EPPGapMaterial::commitState(void)
{
    if (fy >= 0) {
       if (trialStrain > maxElasticYieldStrain)  {
           maxElasticYieldStrain = trialStrain;
           minElasticYieldStrain = trialStrain-trialStress/E;
       }
       else if (trialStrain < minElasticYieldStrain && trialStrain > gap
                && damage == 0 )  {
           maxElasticYieldStrain = (trialStrain-eta*gap)/(1-eta)+fy/E;
           minElasticYieldStrain = trialStrain;
       }
    }
    else {
       if (trialStrain < maxElasticYieldStrain)  {
           maxElasticYieldStrain = trialStrain;
           minElasticYieldStrain = trialStrain-trialStress/E;
       }
       else if (trialStrain > minElasticYieldStrain && trialStrain < gap
                && damage == 0 )  {
           maxElasticYieldStrain = (trialStrain-eta*gap)/(1-eta)+fy/E;
           minElasticYieldStrain = trialStrain;
       }
    }

	//added by SAJalali for energy recorder
	EnergyP += 0.5*(commitStress + trialStress)*(trialStrain - commitStrain);

	commitStrain = trialStrain;
	commitStress = trialStress;	//added by SAJalali
    commitTangent = trialTangent;//added by ambaker1

	return 0;
}


int 
EPPGapMaterial::revertToLastCommit(void)
{
    trialStrain = commitStrain;
	trialStress = commitStress;	//added by SAJalali
    trialTangent = commitTangent; //added by ambaker1

    return 0;
}


int 
EPPGapMaterial::revertToStart(void)
{
    trialStrain = 0.0;
    maxElasticYieldStrain = fy/E+gap;
    minElasticYieldStrain = gap;

    //use setTrialStrain method to get trial stress and tangent. Save as committed. Added by ambaker1
    this->setTrialStrain(trialStrain);
    commitStrain = trialStrain;
    commitStress = trialStress;
    commitTangent = trialTangent;

// AddingSensitivity:BEGIN /////////////////////////////////
    if (SHVs != 0) 
      SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

    return 0;
}


UniaxialMaterial *
EPPGapMaterial::getCopy(void)
{
    EPPGapMaterial *theCopy = new EPPGapMaterial(this->getTag(),E,fy,gap,eta,damage);
    theCopy->trialStrain = trialStrain;
    theCopy->trialStress = trialStress;
    theCopy->trialTangent = trialTangent;
    theCopy->commitStrain = commitStrain;
    theCopy->commitStress = commitStress;
    theCopy->commitTangent = commitTangent;
    theCopy->maxElasticYieldStrain = maxElasticYieldStrain;
    theCopy->minElasticYieldStrain = minElasticYieldStrain;
    theCopy->EnergyP = EnergyP;
    theCopy->parameterID = parameterID;

    return theCopy;
}


int 
EPPGapMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(11);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = fy;
  data(3) = gap;
  data(4) = eta;
  data(5) = maxElasticYieldStrain;
  data(6) = minElasticYieldStrain;
  data(7) = damage;
  data(8) = commitStrain;
  data(9) = commitStress;
  data(10) = commitTangent;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "EPPGapMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
EPPGapMaterial::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(11);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "EPPGapMaterial::recvSelf() - failed to recv data\n";
  else {
    this->setTag((int)data(0));
    E = data(1);
    fy = data(2);
    gap = data(3);
    eta = data(4);
    maxElasticYieldStrain = data(5);
    minElasticYieldStrain = data(6);
    damage = (int)data(7);
    commitStrain = data(8);
    commitStress = data(9);
    commitTangent = data(10);

    trialStrain = commitStrain;
    trialStress = commitStress;
    trialTangent = commitTangent;
  }

  return res;
}

void 
EPPGapMaterial::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "EPPGap tag: " << this->getTag() << endln;
		s << "  E: " << E << ", kinematic hardening ratio: " << eta << endln;
		s << "  fy: " << fy << endln;
		s << "  initial gap: " << gap << endln;
		if (damage == 1)
			s << "  damage accumulation specified" << endln;
	}
	
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"EPPGap\", ";
		s << "\"E\": " << E << ", ";
		s << "\"eta\": " << eta << ", ";
		s << "\"fy\": " << fy << ", ";
		s << "\"gap\": " << gap << ", ";
		s << "\"damageFlag\": " << damage << "}";
	}
}

int
EPPGapMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"Fy") == 0 || strcmp(argv[0],"fy") == 0) {
    param.setValue(fy);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"gap") == 0) {
    param.setValue(gap);
    return param.addObject(3, this);
  }

  return 0;
}

int
EPPGapMaterial::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    E = info.theDouble;
    break;
  case 2:
    fy = info.theDouble;
    break;
  case 3:
    gap = info.theDouble;
    break;
  default:
    return -1;
  }

  return 0;
}

int
EPPGapMaterial::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}

double
EPPGapMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  double dsdh = 0.0;

  double dEdh = 0.0;
  double dfydh = 0.0;
  double dgapdh = 0.0;

  if (parameterID == 1)
    dEdh = 1.0;
  if (parameterID == 2)
    dfydh = 1.0;
  if (parameterID == 3)
    dgapdh = 1.0;

  double depsyMindh = 0.0;
  if (SHVs != 0) {
    depsyMindh = (*SHVs)(0,gradIndex);
  }

  if (fy >= 0) {
    if (trialStrain > maxElasticYieldStrain) {
      dsdh = dfydh + (-dgapdh-dfydh/E+fy/(E*E)*dEdh)*eta*E + (trialStrain-gap-fy/E)*eta*dEdh;
    } else if (trialStrain < minElasticYieldStrain) {
      dsdh = 0.0;
    } else {
      dsdh = dEdh*(trialStrain-minElasticYieldStrain) - E*depsyMindh;
    }
  } else {
    if (trialStrain < maxElasticYieldStrain) {
      dsdh = dfydh + (-dgapdh-dfydh/E+fy/(E*E)*dEdh)*eta*E + (trialStrain-gap-fy/E)*eta*dEdh;
    } else if (trialStrain > minElasticYieldStrain) {
      dsdh = 0.0;
    } else {
      dsdh = dEdh*(trialStrain-minElasticYieldStrain) - E*depsyMindh;
    }
  }

  return dsdh;
}

double 
EPPGapMaterial::getTangentSensitivity(int gradIndex)
{
  return 0.0;
}

double
EPPGapMaterial::getInitialTangentSensitivity(int gradIndex)
{
  if (parameterID == 1)
    return 1.0;
  else
    return 0.0;
}

int
EPPGapMaterial::commitSensitivity(double dedh, int gradIndex, int numGrads)
{
  if (SHVs == 0) {
    SHVs = new Matrix(1,numGrads);
  }
  
  //return 0;

  if (gradIndex >= SHVs->noCols()) {
    //opserr << gradIndex << ' ' << SHVs->noCols() << endln;
    return 0;
  }
  
  double dEdh = 0.0;
  double dfydh = 0.0;
  double dgapdh = 0.0;

  if (parameterID == 1)
    dEdh = 1.0;
  if (parameterID == 2)
    dfydh = 1.0;
  if (parameterID == 3)
    dgapdh = 1.0;

  double depsyMindh = (*SHVs)(0,gradIndex);

  if (fy >= 0) {
    if (trialStrain > maxElasticYieldStrain) {
      double dsdh = this->getStressSensitivity(gradIndex, true);
      dsdh += eta*E*dedh; // unconditional
      depsyMindh = dedh - (-trialStress)/(E*E)*dEdh - dsdh/E;
    }
    else if (trialStrain < minElasticYieldStrain && trialStrain > gap && damage == 0)
      depsyMindh = dedh;
  } else {
    if (trialStrain < maxElasticYieldStrain) {
      double dsdh = this->getStressSensitivity(gradIndex, true);
      dsdh += eta*E*dedh; // unconditional
      depsyMindh = dedh - (-trialStress)/(E*E)*dEdh - dsdh/E;
    }
    else if (trialStrain > minElasticYieldStrain && trialStrain < gap && damage == 0)
      depsyMindh = dedh;
  }

  (*SHVs)(0,gradIndex) = depsyMindh;

  return 0;
}
