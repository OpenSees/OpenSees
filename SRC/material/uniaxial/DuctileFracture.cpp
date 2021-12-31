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

// $Revision: 2.0 - adding comments and minor patches in initialization $
// $Date: 2021-12-30 $
// $Revision: 1.0 $
// $Date: 2019-08-22 $

// Written: Kuanshi Zhong 
// Created: 07/2019
//

#include <math.h>
#include <stdlib.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <string.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <elementAPI.h>
#include "DuctileFracture.h"

void* OPS_DuctileFracture()
{
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 5) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: uniaxialMaterial DuctileFracture tag? matTag?";
		opserr << " -c_mono c_mono? -c_cycl c_cycl? -c_symm c_symm?" << endln;
		opserr << " <-E_s E_s> <-esu esu> <-k1 k1> <-k2 k2> " << endln;
		opserr << " <-db db> <-b1 b1> <-b2 b2> <-FImax FImax?> " << endln;
    opserr << " <-c_dete c_dete> <-minStrain minStraing?> <-maxStrain maxStrain?>" << endln;
		return 0;
	}

  // material tags
	int idata[2];
	numdata = 2;
	if (OPS_GetIntInput(&numdata, idata) < 0) {
		opserr << "WARNING invlid int inputs\n";
		return 0;
	}

	// default parameter values
  double FImax = 1.0;
	double E_s = 29000;
  double c_mono;
  double c_cycl;
  double c_symm;
	double epsmin = NEG_INF_STRAIN;
	double epsmax = POS_INF_STRAIN;
	double esu = POS_INF_STRAIN; // necking strain (esu = 0: no necking)
	double k1 = 1; // necking local strain amplification (k1 = 1: no amplification) 
	double k2 = 0; // necking local Triaxiality amplification (k2 = 0: no amplification)
	double db = 0; // bar-diameter for buckle amplification
	double b1 = 0; // buckle initial-imperfection coefficient (b1 = 0: no effect)
	double b2 = 1000.0; // buckle curvature coefficient (b2 = 1000: no local buckle)
	double c_dete = 0.0; // deteriroation coefficient (c_dete = 0: no deterioration)

  // modeling parameter inputs
	numdata = 1;
	while (OPS_GetNumRemainingInputArgs() > 1) {
		const char* type = OPS_GetString();
		if (strcmp(type, "-FImax") == 0) {
			if (OPS_GetDouble(&numdata, &FImax) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-c_mono") == 0) {
			if (OPS_GetDouble(&numdata, &c_mono) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-c_cycl") == 0) {
			if (OPS_GetDouble(&numdata, &c_cycl) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-c_symm") == 0) {
			if (OPS_GetDouble(&numdata, &c_symm) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-E_s") == 0) {
			if (OPS_GetDouble(&numdata, &E_s) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-esu") == 0) {
			if (OPS_GetDouble(&numdata, &esu) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-k1") == 0) {
			if (OPS_GetDouble(&numdata, &k1) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-k2") == 0) {
			if (OPS_GetDouble(&numdata, &k2) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-db") == 0) {
			if (OPS_GetDouble(&numdata, &db) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-b1") == 0) {
			if (OPS_GetDouble(&numdata, &b1) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-b2") == 0) {
			if (OPS_GetDouble(&numdata, &b2) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-c_dete") == 0) {
			if (OPS_GetDouble(&numdata, &c_dete) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
    else if (strcmp(type, "-minStrain") == 0) {
			if (OPS_GetDouble(&numdata, &epsmin) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
		else if (strcmp(type, "-maxStrain") == 0) {
			if (OPS_GetDouble(&numdata, &epsmax) < 0) {
				opserr << "WARNING invalid double inputs\n";
				return 0;
			}
		}
  }
		

		UniaxialMaterial* mat = OPS_getUniaxialMaterial(idata[1]);
		if (mat == 0) {
			opserr << "WARNING component material does not exist\n";
			opserr << "Component material: " << idata[1];
			opserr << "\nuniaxialMaterial DuctileFracture: " << idata[0] << endln;
			return 0;
		}

		UniaxialMaterial* theMat = new DuctileFracture(idata[0], *mat,
			c_mono, c_cycl, c_symm, E_s, esu, k1, k2, 
			db, b1, b2, FImax, c_dete, epsmin, epsmax);
		if (theMat == 0) {
			opserr << "WARNING: failed to create DuctileFracture material\n";
			return 0;
		}

		return theMat;
}

DuctileFracture::DuctileFracture(int tag,UniaxialMaterial &material,
	double c, double lambda, double beta,
	double Es, double e_su, double k_1, double k_2,
	double d_b, double b_1, double b_2, double FI_max, 
  double cdete, double eps_min, double eps_max)
	:UniaxialMaterial(tag,MAT_TAG_DuctileFracture), theMaterial(0), 
	fracFailure(false), trialStrain(0)
{
  // status variable initializations
  FI  = 0; // Fracture index
  FI_VGM	= 0; // Void growth damage component
  FI_MVC	= 0; // Multi-void coalescence damage compoent
  ep_prev	= 0; // Previous plastic strain
  ep_curr	= 0; // Current plastic strain
  dep		= 0; // Incremental plastic strain
  cep_comp	= 0; // Cumulative compressive plastic strain
  es_local = 0; // Local strain
  T = 0; // Triaxiality
  es_max = 0; // The maximum steel strain
  es_min = 0; // The minimum steel strain
  e_memo = 0; // The strain memory factor

  if ( FI_max > 10.0 || FI_max < 0.0 ) {
    opserr << "DuctileFracture::DuctileFracture " <<
      "-FImax must be between 0 and 10, assuming FImax = 1\n" ;
    FImax    = 1.0;
  } else 
    FImax    = FI_max;
  
  c_mono	= c;
  c_cycl	= lambda;
  c_symm	= beta;
  E_s		= Es;
  minStrain = eps_min;
  maxStrain = eps_max;
  esu = e_su;
  k1 = k_1;
  k2 = k_2;
  db = d_b;
  b1 = b_1;
  b2 = b_2;
  c_dete = cdete;
  
  theMaterial = material.getCopy();
  if (theMaterial == 0) {
    opserr <<  "DuctileFracture::DuctileFracture " <<
      " -- failed to get copy of material\n" ;
    exit(-1);
  }
}

DuctileFracture::DuctileFracture()
	:UniaxialMaterial(0, MAT_TAG_DuctileFracture), theMaterial(0),
	fracFailure(false), trialStrain(0)
{
	FI = 0; // Fracture index
	FI_VGM = 0; // Void growth damage component
	FI_MVC = 0; // Multi-void coalescence damage compoent
	ep_prev = 0; // Previous plastic strain
	ep_curr = 0; // Current plastic strain
	dep = 0; // Incremental plastic strain
	cep_comp = 0; // Cumulative compressive plastic strain
	es_local = 0; // Local strain
	T = 0; // Triaxiality
	es_max = 0; // The maximum steel strain
	es_min = 0; // The minimum steel strain
	e_memo = 0; // The strain memory factor

	FImax = 1.0;
	c_mono = 0;
	c_cycl = 0;
	c_symm = 0;
	E_s = 29000;
	minStrain = NEG_INF_STRAIN;
	maxStrain = POS_INF_STRAIN;
	esu = 0;
	k1 = 0;
	k2 = 0;
	db = 0;
	b1 = 0;
	b2 = 1000.0;
	c_dete = 0;
}

DuctileFracture::~DuctileFracture()
{
  if (theMaterial)
    delete theMaterial;
}

int 
DuctileFracture::setTrialStrain(double strain, double strainRate)
{
  if (fracFailure) {
    trialStrain = strain;
    // return 0;
    return theMaterial->setTrialStrain(strain, strainRate);
  } else {
    fracFailure = false;
    trialStrain = strain;
    return theMaterial->setTrialStrain(strain, strainRate);
  }
}

double 
DuctileFracture::getStress(void)
{
  double modifier = 1.0;
  double damageloc = 1.0-FImax+FI;
  if (fracFailure)
	  // Reduce stress to 0.0 
	  return theMaterial->getStress()*1.0e-8;
  else if (FI_MVC > 1.0) {
	  modifier = 1.0 / sqrt(pow(FI_MVC, c_dete));
	  return theMaterial->getStress()*modifier;
  }
  else
    return theMaterial -> getStress();
}

double 
DuctileFracture::getTangent(void)
{
  double modifier = 1.0;
  double damageloc = 1.0-FImax+FI;
  if (fracFailure)
    // Reduce tangent to 0.0 
    return 1.0e-8*theMaterial->getInitialTangent();
  else if (FI_MVC > 1.0) {
	  modifier = 1.0 / sqrt(pow(FI_MVC, c_dete));
	  return theMaterial->getTangent()*modifier;
  }
  else
    return theMaterial->getTangent()*modifier;
}

double 
DuctileFracture::getDampTangent(void)
{
  if (fracFailure)
    return 0.0;
  else
    return theMaterial->getDampTangent();
}

double 
DuctileFracture::getStrain(void)
{
  return theMaterial->getStrain();
}

double 
DuctileFracture::getStrainRate(void)
{
  return theMaterial->getStrainRate();
}

int 
DuctileFracture::commitState(void)
{	
  if (fracFailure) {
    return 0;
  }

  if (trialStrain >= maxStrain || trialStrain <= minStrain) { 
      fracFailure = true;      
      opserr << "DuctileFracture: material tag " << this->getTag() << " failed from excessive strain\n";
      FI = FImax;
      return 0;
  }
  
  // Update es_max
  if (trialStrain > es_max) {
	  es_max = trialStrain;
  }

  // Update es_min
  if (trialStrain < es_min) {
	  es_min = trialStrain;
  }

  // strain memory for c_symm and c_cycl
  e_memo = fabs(es_max - es_min) / 0.05;
  if (e_memo > 1.0) {
	  e_memo = 1.0;
  }

  // Get current stress
  double s_curr = theMaterial->getStress();

  // Necking amplification model
  if (trialStrain > esu) {
	  // Apply necking coefficient
	  es_local = esu + k1 * (trialStrain - esu);
	  T = 0.33 + k2 * (trialStrain - esu);
  }
  else {
	  es_local = trialStrain;
	  T = 0.33;
  }

  // Buckle amplification model
  es_local = es_local - 0.5 * db * (b1 * sinh((es_max - trialStrain) / b2));

  // Compute current plastic strain
  ep_curr = es_local-s_curr/E_s;
  
  // Compute current incremental plastic strain
  dep = ep_curr - ep_prev;  
  
  // Constant triaxiality T = 0.33
 
  // Compute fracture index components
  if (dep > 0) {
	  FI_VGM = FI_VGM+c_mono*(((c_symm-1.0)*e_memo+1.0)*exp(1.3*T)-exp(-1.3*T))*fabs(dep);
  } else if (dep < 0) {
	  FI_VGM = FI_VGM+c_mono*(((c_symm-1.0)*e_memo+1.0)*exp(-1.3*T)-exp(1.3*T))*fabs(dep);
	  if (FI_VGM < 0) {
		  FI_VGM = 0;
	  }
	  cep_comp = cep_comp+fabs(dep);
  }
  FI_MVC = exp(c_cycl*e_memo*cep_comp);
  
  // Compute fracture index
  FI = FI_VGM*FI_MVC;

  // Flag failure if we have reached that point
  if (FI >= FImax ) {
	  fracFailure = true;
	  opserr << "DuctileFracture: material tag " << this->getTag() << " failed\n";
	  } else {
      fracFailure = false;
    }
  
  // Update previous plastic strain
  ep_prev = ep_curr;
  
  // Check if failed at current step
  if (fracFailure) {
    return 0;
  } else {
    return theMaterial->commitState();
  }
    
}

int 
DuctileFracture::revertToLastCommit(void)
{
  // Check if failed at last step
  if (fracFailure) {
    return 0;
  } else {
    return theMaterial->revertToLastCommit();
  }
}

int 
DuctileFracture::revertToStart(void)
{
  fracFailure = false;
  return theMaterial->revertToStart();
}

UniaxialMaterial *
DuctileFracture::getCopy(void)
{
  DuctileFracture *theCopy = 
    new DuctileFracture(this->getTag(), *theMaterial, c_mono, c_cycl, c_symm, E_s, esu, k1, k2, db, b1, b2, FImax, c_dete, minStrain, maxStrain);

  theCopy->fracFailure = fracFailure;
  theCopy->trialStrain = trialStrain;

  return theCopy;
}

int 
DuctileFracture::sendSelf(int cTag, Channel &theChannel)
{
    int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if ( matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "DuctileFracture::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(25);
  dataVec(0)  = FI;
  dataVec(1)  = FI_VGM;
  dataVec(2)  = FI_MVC;
  dataVec(3)  = ep_prev;
  dataVec(4)  = ep_curr;
  dataVec(5)  = dep;
  dataVec(6)  = cep_comp;
  dataVec(7)  = FImax;
  dataVec(8)  = c_mono;
  dataVec(9)  = c_cycl;
  dataVec(10) = c_symm;
  dataVec(11) = E_s;
  dataVec(12) = minStrain;
  dataVec(13) = maxStrain;
  dataVec(14) = esu;
  dataVec(15) = k1;
  dataVec(16) = k2;
  dataVec(17) = db;
  dataVec(18) = b1;
  dataVec(19) = b2;
  dataVec(20) = es_max;
  dataVec(21) = es_min;
  dataVec(22) = e_memo;
  dataVec(23) = c_dete;

  if (fracFailure == true)
    dataVec(24) = 1.0;
  else
    dataVec(24) = 0.0;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "DuctileFracture::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "DuctileFracture::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
DuctileFracture::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "DuctileFracture::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "DuctileFracture::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(25);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "DuctileFracture::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  FI		= dataVec(0);
  FI_VGM	= dataVec(1);
  FI_MVC	= dataVec(2);
  ep_prev	= dataVec(3);
  ep_curr	= dataVec(4);
  dep		= dataVec(5);
  cep_comp	= dataVec(6);
  FImax		= dataVec(7);
  c_mono	= dataVec(8);
  c_cycl	= dataVec(9);
  c_symm	= dataVec(10);
  E_s		= dataVec(11);
  minStrain	= dataVec(12);
  maxStrain	= dataVec(13);
  esu		= dataVec(14);
  k1		= dataVec(15);
  k2		= dataVec(16);
  db		= dataVec(17);
  b1		= dataVec(18);
  b2		= dataVec(19);
  es_max	= dataVec(20);
  es_min	= dataVec(21);
  e_memo	= dataVec(22);
  c_dete	= dataVec(23);

  if (dataVec(24) == 1.0)
    fracFailure = true;
  else
    fracFailure = false;

  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "DuctileFracture::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
DuctileFracture::Print(OPS_Stream &s, int flag)
{
	if (flag == 100) {
		s << FI << endln;
	}
	
	if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
		s << "DuctileFracture tag: " << this->getTag() << endln;
		s << "\tMaterial: " << theMaterial->getTag() << endln;
		s << "\tFI: " << FI << " FImax: " << FImax << endln;
		s << "\tc_mono: " << c_mono << " c_cycl: " << c_cycl << " c_symm: " << c_symm << endln;
	}
		
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": \"" << this->getTag() << "\", ";
		s << "\"type\": \"DuctileFracture\", ";
		s << "\"material\": \"" << theMaterial->getTag() << "\", ";
		s << "\"tFI\": " << FI << ", ";
		s << "\"FImax\": " << FImax << ", ";
		s << "\"tc_mono\": " << c_mono << ", ";
		s << "\"c_cycl\": " << c_cycl << ", ";
		s << "\"c_symm\": " << c_symm << ", ";
	}
}

Response* 
DuctileFracture::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
  if (argc == 0) 
    return 0;

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

 // damage 
  else if (strcmp(argv[0],"damage") == 0) {
    theResponse =  new MaterialResponse(this, 5, FI);
    theOutput.tag("ResponseType", "FI");
  }   

  else if (strcmp(argv[0],"failure") == 0) {
    int res;
    theResponse =  new MaterialResponse(this, 6, res);
    theOutput.tag("ResponseType", "Failure");
  }

  else if (strcmp(argv[0], "vgm") == 0) {
	  theResponse = new MaterialResponse(this, 7, FI_VGM);
	  theOutput.tag("ResponseType", "FI_VGM");
  }

  else if (strcmp(argv[0], "mvc") == 0) {
	  theResponse = new MaterialResponse(this, 8, FI_MVC);
	  theOutput.tag("ResponseType", "FI_MVC");
  }
  // end add


  theOutput.endTag();
  return theResponse;
}

int 
DuctileFracture::getResponse(int responseID, Information &matInfo)
{
  static Vector stressStrain(2);

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
    
  case 5:
    matInfo.setDouble(FI);
    return 0;      

  case 6:
    if (fracFailure == true) 
      matInfo.setInt(1);
    else
      matInfo.setInt(0);
    return 0;   

  case 7:
	  matInfo.setDouble(FI_VGM);
	  return 0;

  case 8:
	  matInfo.setDouble(FI_MVC);
	  return 0;
    
  default:      
    return -1;

  }
}

bool 
DuctileFracture::hasFailed(void) {
  return fracFailure;
}


