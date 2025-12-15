/* ****************************************************************** **
** OpenSees - Open System for Earthquake Engineering Simulation    **
** Pacific Earthquake Engineering Research Center            **
** **
** **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
** **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
** **
** Developed by:                                                      **
** Frank McKenna (fmckenna@ce.berkeley.edu)                         **
** Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
** Filip C. Filippou (filippou@ce.berkeley.edu)                     **
** **
** ****************************************************************** */
                                                                        
 //----------------------------------------------------------------------------------------------------------------------------
 // Developed by:
 // Javad Esmaeelpour (jesmaeel@tennessee.edu)    
 // Mark D. Denavit   (mdenavit@utk.edu)           
 // Michael H. Scott  (michael.scott@oregonstate.edu)
 //
 // Based on TDConcrete implementations by:
 // Adam M. Knaack (adam.knaack@schaefer-inc.com) 
 // Schaefer-Inc, Cincinnati, Ohio, USA
 // Nikola D. Tosic (ntosic@imk.grf.bg.ac.rs)
 // Department for Materials and Structure, Faculty of Civil Engineering, University of Belgrade, Serbia
 // Yahya C. Kurama (ykurama@nd.edu)
 // Department of Civil and Environmental Engineering and Earth Sciences, College of Engineering, University of Notre Dame, Notre Dame, Indiana, USA
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Description: This file contains the source code of CreepShrinkageACI209. 
 // CreepShrinkageACI209 is a wrapper that imposes creep and shrinkage evoluation equations
 // to any uniaxialMaterial.
 //----------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "CreepShrinkageACI209.h"
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <elementAPI.h>
#include <Domain.h>
#include <MaterialResponse.h>
#include <Vector.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>

#include <Concrete02IS.h>
#include <ElasticMaterial.h>

static int numCreepShrinkageACI209 = 0;

void *
OPS_CreepShrinkageACI209() {
  if (numCreepShrinkageACI209 == 0) {
    numCreepShrinkageACI209 = 1;
  }
  
  UniaxialMaterial *theMaterial = 0;
  
  int iData;
  int numData;
  int numArgs;
  
  numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs == 13) {
    double dData[12];
    
    numData = 1;
    if (OPS_GetIntInput(&numData, &iData) != 0) {
      opserr << "WARNING: invalid uniaxialMaterial CreepShrinkageACI209 tag\n";
      return 0;
    }
    
    numData = 12;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "WARNING: invalid material property definition\n";
      return 0;
    }
    
    theMaterial = new CreepShrinkageACI209(iData,dData[0],dData[1],dData[2],dData[3],dData[4],dData[5],dData[6],dData[7],dData[8],dData[9],dData[10],dData[11]);
    if (theMaterial == 0) {
      opserr << "WARNING: could not create uniaxialMaterial of type CreepShrinkageACI209 \n";
      return 0;
    }
    
    return theMaterial;
  }
  if (numArgs == 10) {
    double dData[14];
    
    numData = 1;
    if (OPS_GetIntInput(&numData, &iData) != 0) {
      opserr << "WARNING: invalid uniaxialMaterial CreepShrinkageACI209 tag\n";
      return 0;
    }

    int wrappedMatl;
    if (OPS_GetIntInput(&numData, &wrappedMatl) != 0) {
      opserr << "WARNING: invalid uniaxialMaterial CreepShrinkageACI209 wrapped material tag\n";
      return 0;
    }

    UniaxialMaterial *matl = OPS_getUniaxialMaterial(wrappedMatl);
    if (matl == 0) {
      opserr << "WARNING: CreepShrinkageACI209 - unable to find material with tag " << wrappedMatl << endln;
      return 0;
    }
    
    numData = 8;
    if (OPS_GetDoubleInput(&numData, &dData[6]) != 0) {
      opserr << "WARNING: invalid material property definition\n";
      return 0;
    }
    
    theMaterial = new CreepShrinkageACI209(iData,*matl,dData[6],dData[7],dData[8],dData[9],dData[10],dData[11],dData[12],dData[13]);
    if (theMaterial == 0) {
      opserr << "WARNING: could not create uniaxialMaterial of type CreepShrinkageACI209 \n";
      return 0;
    }
    
    return theMaterial;
  }  
  
  return 0;
}

CreepShrinkageACI209::CreepShrinkageACI209(int tag, double _fc, double _fcu, double _epscu, double _Ec, double _age, double _epsshu, double _epssha, double _tcr, double _epscru, double _epscra, double _epscrd, double _tcast): 
  UniaxialMaterial(tag, MAT_TAG_CreepShrinkageACI209), wrappedMaterial(0),
  Ec(_Ec), age(_age), epsshu(_epsshu), epssha(_epssha), tcr(_tcr), epscru(_epscru), epscra(_epscra), epscrd(_epscrd), tcast(_tcast), maxSize(startSize),
  DSIG_i(0), TIME_i(0), committedCreepStress(0.0)
{
  wrappedMaterial = new Concrete02IS(0,Ec,-fabs(_fc),2*(-fabs(_fc))/Ec,_fcu,_epscu);
  
  epsshu = -fabs(epsshu);
  epscru = fabs(epscru);

  this->revertToStart(); 
  this->expandArrays();
  DSIG_i[0] = 0.0;
  TIME_i[0] = getCurrentTime();
}

void
CreepShrinkageACI209::expandArrays()
{

  if (DSIG_i == 0)
    DSIG_i = new double[maxSize];

  if (TIME_i == 0)
    TIME_i = new double[maxSize];

  if (historyPointCount+1 >= maxSize) {
    maxSize += growSize;

    double *c = new double[maxSize];
    double *e = new double[maxSize];
    for (int i = 0; i <= historyPointCount; i++) {
      c[i] = DSIG_i[i];
      e[i] = TIME_i[i];
    }

    if (DSIG_i != 0)
      delete [] DSIG_i;
    if (TIME_i != 0)
      delete [] TIME_i;
    DSIG_i = c;
    TIME_i = e;
  }
}

CreepShrinkageACI209::CreepShrinkageACI209(int tag, UniaxialMaterial &matl, double _age, double _epsshu, double _epssha, double _tcr, double _epscru, double _epscra, double _epscrd, double _tcast): 
  UniaxialMaterial(tag, MAT_TAG_CreepShrinkageACI209), wrappedMaterial(0),
  age(_age), epsshu(_epsshu), epssha(_epssha), tcr(_tcr), epscru(_epscru), epscra(_epscra), epscrd(_epscrd), tcast(_tcast), maxSize(startSize),
  DSIG_i(0), TIME_i(0), committedCreepStress(0.0)
{
  wrappedMaterial = matl.getCopy();
  if (wrappedMaterial == 0) {
    opserr << "CreepShrinkageACI209::CreepShrinkageACI209 - failed to get copy of material" << endln;
    exit(-1);
  }  

  Ec = wrappedMaterial->getInitialTangent();
    
  epsshu = -fabs(epsshu);
  epscru = fabs(epscru);
  
  this->revertToStart();

  this->expandArrays();
  DSIG_i[0] = 0.0;
  TIME_i[0] = getCurrentTime();
}

CreepShrinkageACI209::CreepShrinkageACI209(void):
  UniaxialMaterial(0, MAT_TAG_CreepShrinkageACI209), wrappedMaterial(0), maxSize(startSize),
  DSIG_i(0), TIME_i(0), committedCreepStress(0.0)
{
}

CreepShrinkageACI209::~CreepShrinkageACI209(void)
{
  if (wrappedMaterial != 0)
    delete wrappedMaterial;
  if (DSIG_i != 0)
    delete [] DSIG_i;
  if (TIME_i != 0)
    delete [] TIME_i;
}

UniaxialMaterial* CreepShrinkageACI209::getCopy(void)
{
  CreepShrinkageACI209 *theCopy = new CreepShrinkageACI209(this->getTag(), *wrappedMaterial, 
                                            age, epsshu, epssha, tcr, epscru, epscra, epscrd, tcast); 

  // Copy state variables
  theCopy->trialStress = trialStress;
  theCopy->trialTangent = trialTangent;
  theCopy->trialTotalStrain = trialTotalStrain;
  theCopy->trialCreepStrain = trialCreepStrain;
  theCopy->trialShrinkageStrain = trialShrinkageStrain;
  theCopy->trialMechanicalStrain = trialMechanicalStrain;
  
  theCopy->committedStress = committedStress;
  theCopy->committedTangent = committedTangent;
  theCopy->committedCreepStrain = committedCreepStrain;
  theCopy->committedCreepStress = committedCreepStress;
  theCopy->committedShrinkageStrain = committedShrinkageStrain;
  theCopy->committedMechanicalStrain = committedMechanicalStrain;
  theCopy->committedTotalStrain = committedTotalStrain;

  theCopy->iterationInStep = iterationInStep;
  theCopy->historyPointCount = historyPointCount;
  
  // Ensure array capacity and copy history
  theCopy->maxSize = maxSize;
  theCopy->DSIG_i = new double[maxSize];
  theCopy->TIME_i = new double[maxSize];
  
  for (int i = 0; i <= historyPointCount; i++) {
    theCopy->DSIG_i[i] = DSIG_i[i];
    theCopy->TIME_i[i] = TIME_i[i];
  }
  
  return theCopy;
}

double
CreepShrinkageACI209::getInitialTangent(void)
{
  return wrappedMaterial->getInitialTangent();
}

double
CreepShrinkageACI209::getCurrentTime(void)
{
  double currentTime = 0.0;
  Domain * theDomain = ops_TheActiveDomain;
  if (theDomain != 0) {
    currentTime = theDomain->getCurrentTime();
  }
  return currentTime;
}	

double
CreepShrinkageACI209::setCreepStrain(double time)
{
  double creep = 0.0, phi = 0.0;
  const double denom_tcr = 1.25 * pow(tcr, -0.118);
  double tmtp = 0.0;
  double a;
  double f2;
  double f3;

  for (int i = 1; i <= historyPointCount; ++i) {
    tmtp = time - TIME_i[i];
    if (tmtp <= 0.0) continue; // Guard against negative time.

    a  = pow(tmtp, epscra);
    f2 = (epscru * a) / (epscrd + a);
    f3 = (1.25 * pow(TIME_i[i] - tcast, -0.118)) / denom_tcr;
    phi = f2 * f3;

    creep += phi * DSIG_i[i] / Ec;
  }
  return creep;
}

double 
CreepShrinkageACI209::setShrink(double time)
{
  double tD = age;
  double shrink = 0.0;
  if (time-(tD) < 0) {
    shrink = 0.0;
  } else {
    shrink = (time-(tD)) / (epssha + (time - (tD))) * epsshu;
  }
  return shrink;
}

int
CreepShrinkageACI209::setTrialStrain(double strain, double strainRate)
{
  double t = getCurrentTime();
  
  
  // Enforce non-negative time only when creep is enabled.
  if (ops_Creep == 1 && t < 0.0) {
    opserr << "CreepShrinkageACI209::setTrialStrain - Negative time (" << t
           << ") not allowed when creep is active.\n";
    return -1;
  }

  trialTotalStrain = strain;
  
  if (ops_Creep == 1) {
    if (fabs(t-TIME_i[historyPointCount]) <= 0.0001) {
      trialCreepStrain = committedCreepStrain;
      trialShrinkageStrain = committedShrinkageStrain;
    } else {
      if (iterationInStep < 1) {
        trialCreepStrain = setCreepStrain(t);
        trialShrinkageStrain = setShrink(t);
      }
    }
  } else {
    trialCreepStrain = committedCreepStrain;
    trialShrinkageStrain = committedShrinkageStrain;
  }
  
  trialMechanicalStrain = trialTotalStrain - trialCreepStrain - trialShrinkageStrain;
  wrappedMaterial->setTrialStrain(trialMechanicalStrain, strainRate);
  trialStress = wrappedMaterial->getStress();
  trialTangent = wrappedMaterial->getTangent();	

  iterationInStep ++;
  return 0;
}

double
CreepShrinkageACI209::getStrain(void)
{
  return trialTotalStrain;
}

double 
CreepShrinkageACI209::getStress(void)
{
  return trialStress;
}

double 
CreepShrinkageACI209::getTangent(void)
{
  return trialTangent;
}


int CreepShrinkageACI209::commitState(void)
{
  iterationInStep = 0;

  if (ops_Creep == 1 && fabs(trialStress - committedCreepStress) > 1e-12) {
    double dSig = trialStress - committedCreepStress;
    this->expandArrays();
    historyPointCount++;
    DSIG_i[historyPointCount] = dSig;
    TIME_i[historyPointCount] = getCurrentTime();

    committedCreepStress = trialStress;
  }

  committedTangent           = trialTangent;
  committedStress            = trialStress;
  committedTotalStrain       = trialTotalStrain;
  committedShrinkageStrain   = trialShrinkageStrain;
  committedCreepStrain       = trialCreepStrain;
  committedMechanicalStrain  = trialMechanicalStrain;

  wrappedMaterial->commitState();
  return 0;
}


int 
  CreepShrinkageACI209::revertToLastCommit(void)
{
  iterationInStep = 0;
  // Restore trial state to last committed state
  trialTotalStrain = committedTotalStrain;
  trialShrinkageStrain = committedShrinkageStrain;
  trialCreepStrain = committedCreepStrain;
  trialMechanicalStrain = committedMechanicalStrain;

  trialTangent = committedTangent;
  trialStress = committedStress;
  
  wrappedMaterial->revertToLastCommit();
  
  return 0;
}

int 
CreepShrinkageACI209::revertToStart(void)
{
  
  committedTangent     = Ec;
  committedStress      = 0.0;
  trialStress          = 0.0;
  committedCreepStress = 0.0;
  trialTangent         = Ec;

 
  historyPointCount = 0;

  
  trialTotalStrain          = 0.0;
  committedTotalStrain      = 0.0;
  trialCreepStrain          = 0.0;
  trialShrinkageStrain      = 0.0;
  trialMechanicalStrain     = 0.0;
  committedCreepStrain      = 0.0;
  committedShrinkageStrain  = 0.0; 
  committedMechanicalStrain = 0.0;

  iterationInStep = 0;

  if (wrappedMaterial) 
    wrappedMaterial->revertToStart();

  return 0;
}

int 
CreepShrinkageACI209::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dbTag = this->getDbTag();
  static ID classTags(4);

  classTags(0) = wrappedMaterial->getClassTag();
  int matDbTag = wrappedMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      wrappedMaterial->setDbTag(matDbTag);
  }
  classTags(1) = matDbTag;
  classTags(2) = this->getTag();
  classTags(3) = maxSize;
  
  res = theChannel.sendID(dbTag, commitTag, classTags);
  if (res < 0) {
    opserr << "CreepShrinkageACI209::sendSelf -- could not send ID" << endln;
    return res;
  }

  Vector data(27 + maxSize*2);
  int i = 0;
  data(i++) = tcr;
  data(i++) = Ec;
  data(i++) = age;
  data(i++) = epsshu;
  data(i++) = epssha;
  data(i++) = epscra;
  data(i++) = epscru;
  data(i++) = epscrd;
  data(i++) = tcast;
  
  data(i++) = committedStress;
  data(i++) = committedTangent;
  data(i++) = committedCreepStress;
  
  data(i++) = historyPointCount;
  data(i++) = trialCreepStrain;
  data(i++) = trialShrinkageStrain;
  data(i++) = trialMechanicalStrain;
  data(i++) = committedMechanicalStrain;
  data(i++) = committedCreepStrain;
  data(i++) = committedShrinkageStrain;
  data(i++) = trialTotalStrain;
  data(i++) = committedTotalStrain;
  data(i++) = iterationInStep;

  for (int j = 0; j < maxSize; j++, i++)
    data(i) = DSIG_i[j];
  for (int j = 0; j < maxSize; j++, i++)
    data(i) = TIME_i[j];
  
  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "CreepShrinkageACI209::sendSelf - failed to send Vector" << endln;
    return res;
  }
  
  res = wrappedMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "CreepShrinkageACI209::sendSelf -- could not send UniaxialMaterial" << endln;
    return res;
  }
	
  return res;
}

int 
CreepShrinkageACI209::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static ID idata(4);
  int dbTag = this->getDbTag();

  res = theChannel.recvID(dbTag, commitTag, idata);
  if (res < 0) {
    opserr << "CreepShrinkageACI209::recvSelf() - failed to receive data\n";
    return res;
  }

  this->setTag(idata(2));  
  maxSize = idata(3);
  
  Vector data(27 + maxSize*2);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "CreepShrinkageACI209::recvSelf() - failed to recvSelf\n";
    return -1;
  }
  
  int i = 0;
  tcr = data(i++);    
  Ec = data(i++);     
  age = data(i++);     
  epsshu = data(i++); 
  epssha = data(i++); 
  epscra = data(i++); 
  epscru = data(i++); 
  epscrd = data(i++); 
  tcast = data(i++); 

  committedStress = data(i++);
  committedTangent = data(i++);
  committedCreepStress = data(i++);
  
  historyPointCount = data(i++);        
  trialCreepStrain = data(i++);         
  trialShrinkageStrain = data(i++);     
  trialMechanicalStrain = data(i++);    
  committedMechanicalStrain = data(i++);
  committedCreepStrain = data(i++);     
  committedShrinkageStrain = data(i++); 
  trialTotalStrain = data(i++);         
  committedTotalStrain = data(i++);     
  iterationInStep = data(i++);          

  if (DSIG_i != 0) delete [] DSIG_i;
  DSIG_i = new double [maxSize];
  if (TIME_i != 0) delete [] TIME_i;
  TIME_i = new double [maxSize];
  
  for (int j = 0; j < maxSize; j++, i++)
    DSIG_i[j]  = data(i);
  for (int j = 0; j < maxSize; j++, i++)
    TIME_i[j]  = data(i);

  trialTangent = committedTangent;
  trialStress = committedStress;
  
  int matClassTag = idata(0);
  
  if (wrappedMaterial == 0) {
    wrappedMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (wrappedMaterial == 0) {
      opserr << "CreepShrinkageACI209::recvSelf -- could not get a UniaxialMaterial" << endln;
      return -1;
    }
  }

  dbTag = idata(1);
  if (wrappedMaterial->getClassTag() != matClassTag) {
    delete wrappedMaterial;
    wrappedMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (wrappedMaterial == 0) {
      opserr << "CreepShrinkageACI209::recvSelf -- could not get a UniaxialMaterial" << endln;
      return -1;
    }
  }
  
  wrappedMaterial->setDbTag(dbTag);
  res = wrappedMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "CreepShrinkageACI209::recvSelf -- count not receive Uniaxialmaterial" << endln;
    return res;
  }
  
  return res;
}

void 
CreepShrinkageACI209::Print(OPS_Stream &s, int flag)
{
  s << "CreepShrinkageACI209:(strain, stress, tangent) " << trialTotalStrain << " " << trialStress << " " << trialTangent << endln;
}


int
CreepShrinkageACI209::getVariable(const char *varName, Information &theInfo)
{
  return -1;
}

Response* CreepShrinkageACI209::setResponse(const char **argv, int argc,
							  OPS_Stream &theOutput)
{	
  Response *theResponse = 0;
  
  theOutput.tag("UniaxialMaterialOutput");
  theOutput.attr("matType", this->getClassType());
  theOutput.attr("matTag", this->getTag());
  
  if (strcmp(argv[0],"stress") == 0) {
    theOutput.tag("ResponseType", "sigma11");
    theResponse =  new MaterialResponse(this, 1, this->getStress());
  } else if (strcmp(argv[0],"tangent") == 0) {
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 2, this->getTangent());
  } else if (strcmp(argv[0],"strain") == 0) {
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 3, this->getStrain());
  } else if ((strcmp(argv[0],"stressStrain") == 0) || (strcmp(argv[0],"stressANDstrain") == 0)) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(this, 4, Vector(2));
  } else if (strcmp(argv[0],"CreepStressStrainTangent")==0) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theOutput.tag("ResponseType", "C11");
    theOutput.tag("ResponseType", "CreepStrain");
    theOutput.tag("ResponseType", "MechStrain");
    theOutput.tag("ResponseType", "ShrinkStrain");
    theResponse = new MaterialResponse(this, 6, Vector(6));
  } else if (strcmp(argv[0],"stressStrainTangent") == 0) {
    theOutput.tag("ResponseType", "sig11");
    theOutput.tag("ResponseType", "eps11");
    theOutput.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(this, 5, Vector(3));
  }
  
  theOutput.endTag();
  return theResponse;
}

int 
CreepShrinkageACI209::getResponse(int responseID, Information &matInfo)
{
  static Vector stressStrain(2);
  static Vector stressStrainTangent(3);
  static Vector CreepStressStrainTangent(6);
  
  switch (responseID) {
  case 1: 
    matInfo.setDouble(trialStress);
    return 0;

  case 2:
    matInfo.setDouble(trialTangent);
    return 0;

  case 3:
    matInfo.setDouble(trialTotalStrain);
    return 0;

  case 4:
    stressStrain(0) = trialStress;
    stressStrain(1) = trialTotalStrain;
    matInfo.setVector(stressStrain);
    return 0;

  case 5:
    stressStrainTangent(0) = trialStress;
    stressStrainTangent(1) = trialTotalStrain;
    stressStrainTangent(2) = trialTangent;
    matInfo.setVector(stressStrainTangent);
    return 0;

  case 6:
    CreepStressStrainTangent(0) = trialStress;
    CreepStressStrainTangent(1) = trialTotalStrain;
    CreepStressStrainTangent(2) = trialTangent;
    CreepStressStrainTangent(3) = trialCreepStrain;
    CreepStressStrainTangent(4) = trialMechanicalStrain;
    CreepStressStrainTangent(5) = trialShrinkageStrain;
    matInfo.setVector(CreepStressStrainTangent);
    return 0;

  default:
    return -1;
  }
}