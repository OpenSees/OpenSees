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

// Documentation: Self-Centering Friction Spring Damper Material
//
// uniaxialMaterial FrictionSpringDamper $tag $Ke $K $Kp $preload <-initialStrain $initStrain>
//
// Required Input Parameters: 
//   $tag         integer tag identifying material
//   $Ke          elastic stiffness
//   $K           loading stiffness
//   $Kp          unloading stiffness
//   $preload     preload
//   $initStrain  initial strain (optional, default = 0.0)
//
// References:
//   1. Wang, W., Fang, C., Zhao, Y., Sause, R., Hu, S., and Ricles, J. (2019). 
//      “Self-centering friction spring dampers for seismic resilience.” 
//      Earthquake Engineering & Structural Dynamics, 48(9), 1045–1065.

#include <elementAPI.h>
#include "FrictionSpringDamper.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

void
localInit()
{
  OPS_Error("FrictionSpringDamper uniaxial material \nWritten by Mark D. Denavit, University of Tennessee, Knoxville\n", 1);
}

void *
OPS_FrictionSpringDamper()
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int numData = OPS_GetNumRemainingInputArgs();
  if (numData < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial FrictionSpringDamper tag? ";
    opserr << "Ke? K? Kp? preload? <-initialStrain initStrain?>" << endln;
    return 0;
  }

  numData = 1;
  int tag;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial FrictionSpringDamper tag \n";
    return 0;
  }
    
  numData = 4;
  double dData[4];
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data for uniaxialMaterial FrictionSpringDamper with tag " << tag << endln;
    return 0;
  }    
  double Ke = dData[0];
  double K = dData[1];
  double Kp = dData[2];
  double preload = dData[3];
  double initStrain = 0.0;

  // Loop through remaining arguments
  while ( OPS_GetNumRemainingInputArgs() > 0 ) {
    const char *sData = OPS_GetString();
    if ( strcmp(sData,"-initialStrain") == 0 ) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &initStrain) != 0) {
        opserr << "WARNING invalid input, want: -initialStrain $initStrain \n";
        return 0;
      }
    } else {
      opserr << "WARNING unknown option " << sData << "\n";
      return 0;
    }
  }

  // Create material
  theMaterial = new FrictionSpringDamper(tag, Ke, K, Kp, preload, initStrain);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type FrictionSpringDamper\n";
    return 0;
  }

  return theMaterial;
}




FrictionSpringDamper::FrictionSpringDamper(int tag, double iKe, double iK, double iKp, double ipreload, double iinitStrain)
:UniaxialMaterial(tag,MAT_TAG_FrictionSpringDamper),
 Ke(iKe), K(iK), Kp(iKp), preload(ipreload), initStrain(iinitStrain)
{
  // Set Remaining Variables
  this->revertToStart();
}

FrictionSpringDamper::FrictionSpringDamper()
:UniaxialMaterial(0,MAT_TAG_FrictionSpringDamper),
 Ke(0.0), K(0.0), Kp(0.0), preload(0.0), initStrain(0.0)
{
  // Set Remaining Variables
  this->revertToStart();
}

FrictionSpringDamper::~FrictionSpringDamper()
{
  // does nothing
}

int FrictionSpringDamper::setTrialStrain(double strain, double strainRate) {

  trialStrain = strain + initStrain;

  double strainIncr = trialStrain - commitStrain;
  double F1 = preload;
  double F2 = preload*Kp/K;
  double Keq1 = 1.0/(1.0/Ke+1.0/K);
  double Keq2 = 1.0/(1.0/Ke+1.0/Kp);
  double backStrain = commitStrain - commitStress/Ke;

  if (strainIncr >= 0.0) {
    if (backStrain < 0.0) {
      if ( trialStrain < (backStrain + (backStrain*Kp-F2)/Ke) ) {
        trialStress = (trialStrain-backStrain)*Ke;
        trialTangent = Ke;
      } else if ( trialStrain < -F2/Ke ) {
        trialStress = -F2 + (trialStrain+F2/Ke)*Keq2;
        trialTangent = Keq2;
      } else if ( trialStrain < F1/Ke ) {
        trialStress = trialStrain*Ke;
        trialTangent = Ke;
      } else {
        trialStress = F1 + (trialStrain-F1/Ke)*Keq1;
        trialTangent = Keq1;
      }
    } else {
      if ( trialStrain < (backStrain + (F1+backStrain*K)/Ke) ) {
        trialStress = (trialStrain-backStrain)*Ke;
        trialTangent = Ke;
      } else {
        trialStress = F1 + (trialStrain-F1/Ke)*Keq1;
        trialTangent = Keq1;
      }
    }

  } else {
    if (backStrain > 0.0) {
      if ( trialStrain > (backStrain + (backStrain*Kp+F2)/Ke) ) {
        trialStress = (trialStrain-backStrain)*Ke;
        trialTangent = Ke;
      } else if ( trialStrain > F2/Ke ) {
        trialStress = F2 + (trialStrain-F2/Ke)*Keq2;
        trialTangent = Keq2;
      } else if ( trialStrain > -F1/Ke ) {
        trialStress = trialStrain*Ke;
        trialTangent = Ke;
      } else {
        trialStress = -F1 + (trialStrain+F1/Ke)*Keq1;
        trialTangent = Keq1;
      }
    } else {
      if ( trialStrain > (backStrain + (-F1+backStrain*K)/Ke) ) {
        trialStress = (trialStrain-backStrain)*Ke;
        trialTangent = Ke;
      } else {
        trialStress = -F1 + (trialStrain+F1/Ke)*Keq1;
        trialTangent = Keq1;
      }
    }

  }

  return 0;
}

double FrictionSpringDamper::getStrain(void) {
  return trialStrain;
}

double FrictionSpringDamper::getStress(void) {
  return trialStress;
}

double FrictionSpringDamper::getTangent(void) {
  return trialTangent;
}

double FrictionSpringDamper::getInitialTangent(void) {
  double Kinit;
  if (initStrain < -preload/Ke) {
    Kinit = K;
  } else if (initStrain <= preload/Ke) {
    Kinit = Ke;
  } else {
    Kinit = K;
  } 
  return Kinit;
}

int FrictionSpringDamper::commitState(void) {
  commitStrain = trialStrain;
  commitStress = trialStress;
  commitTangent = trialTangent;
  return 0;
}

int FrictionSpringDamper::revertToLastCommit(void) {
  trialStrain =  commitStrain;
  trialStress = commitStress;
  trialTangent = commitTangent;
  return 0;
}


int FrictionSpringDamper::revertToStart(void) {
  if (initStrain < -preload/Ke) {
    trialTangent = K;
    trialStress = (initStrain+preload/Ke)*K - preload;
  } else if (initStrain <= preload/Ke) {
    trialTangent = Ke;
    trialStress = initStrain*Ke;
  } else {
    trialTangent = K;
    trialStress = (initStrain-preload/Ke)*K + preload;
  }
  trialStrain = initStrain;
  this->commitState();
  return 0;
}

UniaxialMaterial * FrictionSpringDamper::getCopy(void) {
  FrictionSpringDamper *theCopy = new FrictionSpringDamper(this->getTag(), Ke, K, Kp, preload, initStrain);
  theCopy->trialStrain = this->trialStrain;
  theCopy->trialStress = this->trialStress;
  theCopy->trialTangent = this->trialTangent;
  theCopy->commitStrain = this->commitStrain;
  theCopy->commitStress = this->commitStress;
  theCopy->commitTangent = this->commitTangent;
  return theCopy;
}


int FrictionSpringDamper::sendSelf(int cTag, Channel &theChannel) {
  int res = 0;
  static Vector data(9);
  data(0) = this->getTag();
  data(1) = Ke;
  data(2) = K;
  data(3) = Kp;
  data(4) = preload;
  data(5) = initStrain;
  
  data(6) = commitStrain;
  data(7) = commitStress;
  data(8) = commitTangent;  

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "FrictionSpringDamper::sendSelf() - failed to send data" << endln;

  return res;  
}

int FrictionSpringDamper::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
  int res = 0;
  static Vector data(9);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "FrictionSpringDamper::recvSelf() - failed to recv data" << endln;
  else {
    this->setTag(int(data(0)));
    Ke    = data(1);
    K     = data(2);
    Kp    = data(3);
    preload   = data(4);
    initStrain   = data(5);  
    commitStrain=data(6);
    commitStress=data(7);
    commitTangent=data(8);    
    trialStrain = commitStrain;
    trialStress = commitStress;
    trialTangent = commitTangent;
  }

  return res;
}

void FrictionSpringDamper::Print(OPS_Stream &s, int flag) {
  s << "FrictionSpringDamper tag: " << this->getTag() << endln;
  s << "  Ke: " << Ke << endln;
  s << "  K: " << K << endln;
  s << "  Kp: " << Kp << endln;
  s << "  preload: " << preload << endln;
  s << "  initStrain: " << initStrain << endln;
  s << "  stress: " << trialStress << " tangent: " << trialTangent << endln;
  return;
}

