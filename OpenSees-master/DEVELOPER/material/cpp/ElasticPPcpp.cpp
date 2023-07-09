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
                                                                        
// $Revision: 1.7 $
// $Date: 2009/03/23 23:17:04 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/ElasticPPcpp.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) ElasticPPcpp.C, revA"

#include <elementAPI.h>
#include "ElasticPPcpp.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numElasticPPcpp = 0;

OPS_Export void *
OPS_ElasticPPcpp()
{
  // print out some KUDO's
  if (numElasticPPcpp == 0) {
    opserr << "ElasticPPcpp unaxial material - Written by fmk UC Berkeley Copyright 2008 - Use at your Own Peril\n";
    numElasticPPcpp =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int    iData[1];
  double dData[2];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ElasticPP tag" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid E & ep\n";
    return 0;	
  }

  // 
  // create a new material
  //

  theMaterial = new ElasticPPcpp(iData[0], dData[0], dData[1]);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ElasticPPCpp\n";
    return 0;
  }

  // return the material
  return theMaterial;
}




ElasticPPcpp::ElasticPPcpp(int tag, double e, double eyp)
:UniaxialMaterial(tag, 0),
 ezero(0.0), E(e), ep(0.0),
 trialStrain(0.0), trialStress(0.0), trialTangent(E),
 commitStrain(0.0), commitStress(0.0), commitTangent(E)
{
  fyp = E*eyp;
  fyn = -fyp;
}

ElasticPPcpp::ElasticPPcpp()
:UniaxialMaterial(0, 0),
 fyp(0.0), fyn(0.0), ezero(0.0), E(0.0), ep(0.0),
 trialStrain(0.0), trialStress(0.0), trialTangent(E),
 commitStrain(0.0), commitStress(0.0), commitTangent(E)
{

}

ElasticPPcpp::~ElasticPPcpp()
{
  // does nothing
}

int 
ElasticPPcpp::setTrialStrain(double strain, double strainRate)
{
    if (fabs(trialStrain - strain) < DBL_EPSILON)
      return 0;

    trialStrain = strain;

    double sigtrial;	// trial stress
    double f;		// yield function

    // compute trial stress
    sigtrial = E * ( trialStrain - ezero - ep );

    //sigtrial  = E * trialStrain;
    //sigtrial -= E * ezero;
    //sigtrial -= E *  ep;

    // evaluate yield function
    if ( sigtrial >= 0.0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    double fYieldSurface = - E * DBL_EPSILON;
    if ( f <= fYieldSurface ) {

      // elastic
      trialStress = sigtrial;
      trialTangent = E;

    } else {

      // plastic
      if ( sigtrial > 0.0 ) {
	trialStress = fyp;
      } else {
	trialStress = fyn;
      }

      trialTangent = 0.0;
    }

    return 0;
}

double 
ElasticPPcpp::getStrain(void)
{
  return trialStrain;
}

double 
ElasticPPcpp::getStress(void)
{
  return trialStress;
}


double 
ElasticPPcpp::getTangent(void)
{
  return trialTangent;
}

int 
ElasticPPcpp::commitState(void)
{
    double sigtrial;	// trial stress
    double f;		// yield function

    // compute trial stress
    sigtrial = E * ( trialStrain - ezero - ep );

    // evaluate yield function
    if ( sigtrial >= 0.0 )
	f =  sigtrial - fyp;
    else
	f = -sigtrial + fyn;

    double fYieldSurface = - E * DBL_EPSILON;
    if ( f > fYieldSurface ) {
      // plastic
      if ( sigtrial > 0.0 ) {
	ep += f / E;
      } else {
	ep -= f / E;
      }
    }

    commitStrain = trialStrain;
    commitTangent=trialTangent;
    commitStress = trialStress;

    return 0;
}	


int 
ElasticPPcpp::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialTangent = commitTangent;
  trialStress = commitStress;

  return 0;
}


int 
ElasticPPcpp::revertToStart(void)
{
  trialStrain = commitStrain = 0.0;
  trialTangent = commitTangent = E;
  trialStress = commitStress = 0.0;

  ep = 0.0;

  return 0;
}


UniaxialMaterial *
ElasticPPcpp::getCopy(void)
{
  ElasticPPcpp *theCopy =
    new ElasticPPcpp(this->getTag(),E,fyp/E);
  theCopy->ep = this->ep;
  
  return theCopy;
}


int 
ElasticPPcpp::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(9);
  data(0) = this->getTag();
  data(1) = ep;
  data(2) = E;
  data(3) = ezero;
  data(4) = fyp;
  data(5) = fyn;
  data(6) = commitStrain;
  data(7) = commitStress;
  data(8) = commitTangent;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticPPcpp::sendSelf() - failed to send data\n";

  return res;
}

int 
ElasticPPcpp::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(9);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticPPcpp::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    ep    = data(1);
    E     = data(2);
    ezero = data(3);
    fyp   = data(4);
    fyn   = data(5);  
    commitStrain=data(6);
    commitStress=data(7);
    commitTangent=data(8);
    trialStrain = commitStrain;
    trialTangent = commitTangent;
    trialStress = commitStress;
  }

  return res;
}

void 
ElasticPPcpp::Print(OPS_Stream &s, int flag)
{
  s << "ElasticPPcpp tag: " << this->getTag() << endln;
  s << "  E: " << E << endln;
  s << "  ep: " << ep << endln;
  s << "  stress: " << trialStress << " tangent: " << trialTangent << endln;
}


