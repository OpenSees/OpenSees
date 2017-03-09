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
                                                                       
// $Revision: 1.2
// $Date: 03/2017
// $Source:
                                                                     
// Written: R.Rahimi & R.Sepasdar & Dr. M. R. Banan
// Created: 09/2012
// 
//
// Description: This file contains the class implementation of RambergOsgoodSteel.
//-----------------------------------------------------------------------
//              RAMBERG-OSGOOD STEEL MODEL
//      Developed  by REZA RAHIMI, (Reza.Rahimi@Dal.Ca)   (2012)
//
//	Co-Developer: REZA SEPASDAR,
//	Supervisor: Dr. Mo. R. Banan,
//-----------------------------------------------------------------------
 
#include <elementAPI.h>
#include "RambergOsgoodSteel.h"


#include <Vector.h>
#include <float.h>
#include <Channel.h>
#include <math.h>


static int numRambergOsgoodSteel = 0;

void *
OPS_RambergOsgoodSteel(void)
{
  if (numRambergOsgoodSteel == 0) {
    opserr << "RambergOsgoodSteel unaxial material - Written by R.Rahimi & R.Sepasdar & Dr. Mo. R. Banan Shiraz University Copyright 2012; \n";
    numRambergOsgoodSteel++;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[4];
  double mData[2];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial RambergOsgoodSteel tag" << endln;
    return 0;
  }

 
  numData = 4;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid E & ep\n";
    return 0;  
  }
 
  theMaterial = new RambergOsgoodSteel(iData[0], dData[0], dData[1], dData[2], dData[3]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type RambergOsgoodSteel\n";
    return 0;
  }
  return theMaterial;
}


RambergOsgoodSteel::RambergOsgoodSteel(int tag,
                 double _Fy, double _E0, double _rezaA,
                 double _rezaN):
  UniaxialMaterial(tag, MAT_TAG_RambergOsgoodSteel),
  Fy(_Fy), E0(_E0), rezaAA(_rezaA), rezaNN(_rezaN)
{
  konP = 0;
  kon = 0;
  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;

  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;

}

RambergOsgoodSteel::RambergOsgoodSteel(int tag, double _Fy, double _E0):
  UniaxialMaterial(tag, MAT_TAG_RambergOsgoodSteel),
  Fy(_Fy), E0(_E0)
{
  konP = 0;

  rezaAA=0.002;
  rezaNN=3;

  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;

  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;
}



RambergOsgoodSteel::RambergOsgoodSteel(void):
  UniaxialMaterial(0, MAT_TAG_RambergOsgoodSteel)
{
  konP = 0;
}

RambergOsgoodSteel::~RambergOsgoodSteel(void)
{
  // Does nothing
}

UniaxialMaterial*
RambergOsgoodSteel::getCopy(void)
{
  RambergOsgoodSteel *theCopy = new RambergOsgoodSteel(this->getTag(), Fy, E0, rezaAA, rezaNN);
 
  return theCopy;
}

double
RambergOsgoodSteel::getInitialTangent(void)
{
  return E0;
}

int
RambergOsgoodSteel::setTrialStrain(double trialStrain, double strainRate)
{
  double epsy = Fy / E0;
 
 
  eps = trialStrain;

  double deps = eps - epsP;
 
  epsmax = epsmaxP;
  epsmin = epsminP;
  epspl  = epsplP;
  epss0  = epss0P;  
  sigs0  = sigs0P;
  epsr   = epssrP;  
  sigr   = sigsrP;  
  kon = konP;


  if (kon == 0 || kon == 3) {


    if (fabs(deps) < 10.0*DBL_EPSILON) {

      e = E0;
      sig = 0;  
      kon = 3;    
      return 0;

    } else {
      epsmax = epsy;
      epsmin = -epsy;
      if (deps < 0.0) {
        kon = 2;
        epss0 = epsmin;
        sigs0 = Fy;
        epspl = epsmin;
      } else {
        kon = 1;
        epss0 = epsmax;
        sigs0 = Fy;
        epspl = epsmax;
      }
    }
 
  }


  if (kon == 2 && deps > 0.0) {


    kon = 20;
    epsr = epsP;
    sigr = sigP;
    //sigs0 = fabs(sigP) / (pow((1-(fabs(sigP)/(rezaAA * E0))),(1/rezaNN))); // the code will find a new sigma_y based on the yield offset in case of cyclic loading (ADDED to Rev 1.2)
    //epsmin = min(epsP, epsmin);

  } else if (kon == 1 && deps < 0.0) {

    kon = 10;
    epsr = epsP;
    sigr = sigP;
    //sigs0 = fabs(sigP) / (pow((1-(fabs(sigP)/(rezaAA * E0))),(1/rezaNN))); // (ADD to  Rev 1.2)

  } else if (kon == 10 && sigP <= 0 ) {
    kon = 2;
    epsr = epsP;
    sigr = sigP;
  } else if (kon == 20 && sigP >= 0 ) {
    kon = 1; 
    epsr = epsP;
    sigr = sigP;
  }

  /* (REMOVED in Rev 1.2) 
  if (sigr !=0)
  {
          sigs0=2*Fy;
  }
  */ 
  // calculate current stress sig and tangent modulus E
  if (kon == 2 || kon == 1){
  double trialSig[1000], F[1000], dF[1000], o;
  int kk=1;
  trialSig[1]=1;
  double M=10;
  while ( M >= 0.0001)
    {
     
      F[kk]= (trialSig[kk]/E0) + rezaAA*(pow((trialSig[kk]/sigs0),rezaNN)) - fabs(eps-epsr);
      dF[kk]= (1/E0) + rezaAA*(1/sigs0)*rezaNN*(pow((trialSig[kk]/sigs0),(rezaNN-1)));
      trialSig[kk+1]=trialSig[kk]-(F[kk]/dF[kk]);
     
      kk=kk+1;
      sig=trialSig[kk];
      M=fabs(sig-trialSig[kk-1]);
      if (kk == 1000)
        {
          opserr << "NewtonRaphson method does NOT converge at eps=" <<  eps << "\n";
          M=0;
        }
    }
 
  e = 1/((1/E0) + rezaAA*(1/sigs0)*rezaNN*(pow((sig/sigs0),(rezaNN-1))));
  } else if (kon == 20 || kon == 10) {
    e = E0;
    sig = E0 * fabs(eps-epsr);
  }
    if (eps<epsr){
		sig=-1*sig;
	}
  sig=sig+sigr;
  return 0;
}



double
RambergOsgoodSteel::getStrain(void)
{
  return eps;
}

double
RambergOsgoodSteel::getStress(void)
{
  return sig;
}

double
RambergOsgoodSteel::getTangent(void)
{
  return e;
}

int
RambergOsgoodSteel::commitState(void)
{
  epsminP = epsmin;
  epsmaxP = epsmax;
  epsplP = epspl;
  epss0P = epss0;
  sigs0P = sigs0;
  epssrP = epsr;
  sigsrP = sigr;
  konP = kon;
 
  eP = e;
  sigP = sig;
  epsP = eps;

  return 0;
}

int
RambergOsgoodSteel::revertToLastCommit(void)
{
  epsmin = epsminP;
  epsmax = epsmaxP;
  epspl = epsplP;
  epss0 = epss0P;
  sigs0 = sigs0P;
  epsr = epssrP;
  sigr = sigsrP;
  kon = konP;
 
  e = eP;
  sig = sigP;
  eps = epsP;
  return 0;
}

int
RambergOsgoodSteel::revertToStart(void)
{
  eP = E0;
  epsP = 0.0;
  sigP = 0.0;
  sig = 0.0;
  eps = 0.0;
  e = E0;  

  konP = 0;
  epsmaxP = Fy/E0;
  epsminP = -epsmaxP;
  epsplP = 0.0;
  epss0P = 0.0;
  sigs0P = 0.0;
  epssrP = 0.0;
  sigsrP = 0.0;



  return 0;
}

int
RambergOsgoodSteel::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(15);
  data(0) = Fy;
  data(1) = E0;
  data(2) = epsminP;;
  data(3) = epsmaxP;
  data(4) = epsplP;
  data(5) = epss0P;
  data(6) = sigs0P;
  data(7) = epssrP;
  data(8) = sigsrP;
  data(9) = konP;  
  data(10) = epsP;  
  data(11) = sigP;  
  data(12) = eP;    
  data(13) = this->getTag();
  data(14) = sigini;


  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "RambergOsgoodSteel::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int
RambergOsgoodSteel::recvSelf(int commitTag, Channel &theChannel,
             FEM_ObjectBroker &theBroker)
{
  static Vector data(15);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "RambergOsgoodSteel::recvSelf() - failed to recvSelf\n";
    return -1;
  }

  Fy = data(0);
  E0 = data(1);
  epsminP = data(2);
  epsmaxP = data(3);
  epsplP = data(4);
  epss0P = data(5);
  sigs0P = data(6);
  epssrP = data(7);
  sigsrP = data(8);
  konP = data(9);  
  epsP = data(10);  
  sigP = data(11);  
  eP   = data(12);  
  this->setTag(data(13));
  sigini = data(14);

  e = eP;
  sig = sigP;
  eps = epsP;
 
  return 0;
}

void
RambergOsgoodSteel::Print(OPS_Stream &s, int flag)
{
  s << "RambergOsgoodSteel:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}
