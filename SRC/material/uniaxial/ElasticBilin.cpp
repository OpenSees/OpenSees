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
                                                                        
// $Revision: 1.1 $
// $Date: 2010-08-17 00:23:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticBilin.cpp,v $
                                                                        
// Written: fmk 
// Created: 08/10

// What: "@(#) ElasticBilin.C, revA"


#include <ElasticBilin.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <elementAPI.h>

void *
OPS_NewElasticBilin(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  opserr << argc;

  if (argc != 4 && argc != 6) {
    opserr << "WARNING incoorect num args want: uniaxialMaterial ElasticBilin tag E E2P eps2P <E2N? eps2N?>" << endln;
    return 0;
  }

  int    iData[1];
  double dData[5];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ElasticBilin tag" << endln;
    return 0;
  }
  
  argc--;
  numData = argc;;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: uniaxialMaterial ElasticBilin tag E2P eps2P <E2N? eps2N?>" << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  if (argc == 3) 
    theMaterial = new ElasticBilin(iData[0], dData[0], dData[1], dData[2]);
  else
    theMaterial = new ElasticBilin(iData[0], dData[0], dData[1], dData[3], dData[2], dData[4]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ElasticBilin\n";
    return 0;
  }

  return theMaterial;
}


ElasticBilin::ElasticBilin()
:UniaxialMaterial(0 ,MAT_TAG_ElasticBilin),
 E(0.0), E2P(0.0), E2N(0.0), eps2P(0.0), eps2N(0.0),
 trialStrain(0.0), trialStress(0.0), trialTangent(E)
{
 
}


ElasticBilin::ElasticBilin(int tag, double e, double e2, double eps2)
:UniaxialMaterial(tag, MAT_TAG_ElasticBilin),
 E(e), E2P(e2), E2N(e2), eps2P(eps2), eps2N(-eps2),
 trialStrain(0.0), trialStress(0.0), trialTangent(E)
{
  if (eps2 < 0.0) {
    eps2P = -1*eps2;
    eps2N = eps2;
  }
}

ElasticBilin::ElasticBilin(int tag, double e, double e2p, double e2n, double eps2p, double eps2n)
:UniaxialMaterial(tag, MAT_TAG_ElasticBilin),
 E(e), E2P(e2p), E2N(e2n), eps2P(eps2p), eps2N(eps2n),
 trialStrain(0.0), trialStress(0.0), trialTangent(E)
{
  if (eps2p < 0.0) {
    eps2P = -1*eps2p;
  }
  if (eps2n > 0.0) {
    eps2N = -1*eps2n;
  }
}

ElasticBilin::~ElasticBilin()
{
}

int 
ElasticBilin::setTrialStrain(double strain, double strainRate)
{
    trialStrain = strain;
    
    if (trialStrain >= eps2N && trialStrain <= eps2P) {
      trialTangent = E;
      trialStress = E*trialStrain;
    } else if (trialStrain > eps2P) {
      trialTangent = E2P;
      trialStress = E*eps2P + (trialStrain-eps2P)*E2P;
    } else {
      trialTangent = E2N;
      trialStress = E*eps2N + (trialStrain-eps2N)*E2N;
    }

    this->Print(opserr);
\
    return 0;
}

double 
ElasticBilin::getStrain(void)
{
  return trialStrain;
}

double 
ElasticBilin::getStress(void)
{
  return trialStress;
}


double 
ElasticBilin::getTangent(void)
{
  return trialTangent;
}

int 
ElasticBilin::commitState(void)
{
  return 0;
}	


int 
ElasticBilin::revertToLastCommit(void)
{
  this->setTrialStrain(commitStrain);
  return 0;
}


int 
ElasticBilin::revertToStart(void)
{
  trialStrain = 0;
  trialStress = 0;
  trialTangent = 0;
  commitStrain = 0;

  return 0;
}


UniaxialMaterial *
ElasticBilin::getCopy(void)
{
  ElasticBilin *theCopy =
    new ElasticBilin(this->getTag(), E, E2P, E2N, eps2P, eps2N);
  
  return theCopy;
}


int 
ElasticBilin::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = E2P;
  data(3) = E2N;
  data(4) = eps2P;
  data(5) = eps2N;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticBilin::sendSelf() - failed to send data\n";

  return res;
}

int 
ElasticBilin::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticBilin::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    E       = data(1);
    E2P     = data(2);
    E2N     = data(3);
    eps2P   = data(4);
    eps2N   = data(5);  
  }

  return res;
}

void 
ElasticBilin::Print(OPS_Stream &s, int flag)
{
    s << "ElasticPP tag: " << this->getTag() << endln;
    s << "  E: " << E << endln;
    s << "  E2P: " << E2P << endln;
    s << "  E2N: " << E2N << endln;
    s << "  epsP: " << eps2P << endln;
    s << "  eps2N: " << eps2N << endln;
    s << "strain: "<< trialStrain << " stress: " << trialStress << " tangent: " << trialTangent << endln;
}


