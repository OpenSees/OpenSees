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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticPPMaterial.cpp,v $
                                                                        
                                                                        
// File: ~/material/ElasticPPMaterial.C
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) ElasticPPMaterial.C, revA"


#include <ElasticPPMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


ElasticPPMaterial::ElasticPPMaterial(int tag, double e, double epl)
:UniaxialMaterial(tag,MAT_TAG_ElasticPPMaterial),
 commitStrain(0.0), trialStrain(0.0), E(e), 
 ep(epl), maxElasticYieldStrain(epl), minElasticYieldStrain(-epl)
{

}

ElasticPPMaterial::ElasticPPMaterial()
:UniaxialMaterial(0,MAT_TAG_ElasticPPMaterial),
 commitStrain(0.0), E(0.0), ep(0.0),
 maxElasticYieldStrain(0.0), minElasticYieldStrain(0.0)
{

}

ElasticPPMaterial::~ElasticPPMaterial()
{
  // does nothing
}

int 
ElasticPPMaterial::setTrialStrain(double strain, double strainRate)
{
    trialStrain = strain;
    return 0;
}

double 
ElasticPPMaterial::getStrain(void)
{
    return trialStrain;
}

double 
ElasticPPMaterial::getStress(void)
{
    if (trialStrain >= maxElasticYieldStrain) 
	return E*ep;
    else if (trialStrain <= minElasticYieldStrain)
	return -E*ep;
    else
	return E*(trialStrain-(maxElasticYieldStrain+minElasticYieldStrain)/2);
}



double 
ElasticPPMaterial::getTangent(void)
{
    if (trialStrain >= maxElasticYieldStrain) 
	return 0;
    else if (trialStrain <= minElasticYieldStrain)
	return 0;
    else
	return E;
}

double ElasticPPMaterial::getSecant ()
{
    if (trialStrain != 0.0)
	return this->getStress()/trialStrain;
    else
	return E; 
}

int 
ElasticPPMaterial::commitState(void)
{
    if (trialStrain > maxElasticYieldStrain)  {
	maxElasticYieldStrain = trialStrain;
	minElasticYieldStrain = trialStrain-2*ep;
    }
    else if (trialStrain < minElasticYieldStrain) {
	maxElasticYieldStrain = trialStrain+2*ep;
	minElasticYieldStrain = trialStrain;
    }	

    commitStrain = trialStrain;
    return 0;

}


int 
ElasticPPMaterial::revertToLastCommit(void)
{
    trialStrain = commitStrain;
    return 0;
}


int 
ElasticPPMaterial::revertToStart(void)
{
    commitStrain = 0.0;
    trialStrain = 0.0;
    maxElasticYieldStrain = ep;
    minElasticYieldStrain = -ep;

    return 0;
}


UniaxialMaterial *
ElasticPPMaterial::getCopy(void)
{
    ElasticPPMaterial *theCopy = new ElasticPPMaterial(this->getTag(),E,ep);
    theCopy->trialStrain = trialStrain;
    theCopy->maxElasticYieldStrain = maxElasticYieldStrain;    
    theCopy->minElasticYieldStrain = minElasticYieldStrain;        
    return theCopy;
}


int 
ElasticPPMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = commitStrain;
  data(2) = E;
  data(3) = ep;
  data(4) = maxElasticYieldStrain;
  data(5) = minElasticYieldStrain;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    cerr << "ElasticPPMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
ElasticPPMaterial::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    cerr << "ElasticPPMaterial::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    commitStrain = data(1);
    trialStrain = commitStrain;    
    E = data(2);
    ep = data(3);
    maxElasticYieldStrain = data(4);
    minElasticYieldStrain = data(5);    
  }

  return res;
}

void 
ElasticPPMaterial::Print(ostream &s, int flag)
{
    s << "ElasticPP tag: " << this->getTag() << endl;
    s << "  E: " << E << endl;
    s << "  epsy: " << ep << endl;
}


