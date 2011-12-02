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
// $Date: 2008-12-09 21:23:29 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Neoprene.cpp,v $

// Written: mackie
// Created: 12/2005
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) Neoprene.C, revA"

#include <stdlib.h>

#include "Neoprene.h"
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#include <OPS_Globals.h>

Neoprene::Neoprene(int tag, double e, double gap0)
:UniaxialMaterial(tag,MAT_TAG_EPPGap),
 commitStrain(0.0), trialStrain(0.0), E(e), gap(gap0), minElasticYieldStrain(gap0) 
{
	if (E <= 0.0) {
	  opserr << "Neoprene::Neoprene -- E <= zero\n";
	  exit(-1);
	}
        
        commitStrain = 0.0;

}

Neoprene::Neoprene()
:UniaxialMaterial(0,MAT_TAG_EPPGap),
 E(0.0), gap(0.0), minElasticYieldStrain(0.0)
{

}

Neoprene::~Neoprene()
{
  // does nothing
}

int 
Neoprene::setTrialStrain(double strain, double strainRate)
{
    // set the trial strain
    trialStrain = strain;
    double deps = trialStrain - commitStrain;
    
    // determine trial stress and tangent
    if (gap >= 0) {
        if (trialStrain > gap && deps > 0) {
            // loading with initial stiffness
            trialStress = E*(trialStrain - minElasticYieldStrain);
            trialTangent = E;
        } else if (trialStrain > gap && deps < 0) {
            // unloading along parabolic path
            trialStress = maxElasticYieldStrain*pow(trialStrain - gap,2);
            trialTangent = 2*maxElasticYieldStrain*(trialStrain - gap);
        } else {
            trialStress = 0.0;
            trialTangent = 0.0;
        }
    } else {
        if (trialStrain < gap && deps < 0) {
            // loading with initial stiffness
            trialStress = E*(trialStrain - minElasticYieldStrain);
            trialTangent = E;
        } else if (trialStrain < gap && deps > 0) {
            // unloading along parabolic path
            trialStress = maxElasticYieldStrain*pow(trialStrain - gap,2);
            trialTangent = 2*maxElasticYieldStrain*(trialStrain - gap);
        } else {
            trialStress = 0.0;
            trialTangent = 0.0;
        }
    }
    
    return 0;
}

double 
Neoprene::getStrain(void)
{
    return trialStrain;
}

double 
Neoprene::getStress(void)
{
  return trialStress;

}

double 
Neoprene::getTangent(void)
{
  return trialTangent;
}

double 
Neoprene::getInitialTangent(void)
{
  if (fabs(gap) > 0.0) 
    return 0.0; 
  else 
    return E;
}

int 
Neoprene::commitState(void)
{
    double deps = trialStrain - commitStrain;
    
    if ( gap > 0 && trialStrain > gap ) {
        if ( deps < 0 ) 
            minElasticYieldStrain = trialStrain - trialStress/E;
        else
            maxElasticYieldStrain = trialStress/pow(trialStrain-gap,2);
    } else if ( gap < 0 && trialStrain < gap ) {
        if ( deps > 0 ) 
            minElasticYieldStrain = trialStrain - trialStress/E;
        else
            maxElasticYieldStrain = trialStress/pow(trialStrain-gap,2);
    }

    commitStrain = trialStrain;

    return 0;
}


int 
Neoprene::revertToLastCommit(void)
{
    trialStrain = commitStrain;

    return 0;
}


int 
Neoprene::revertToStart(void)
{
    commitStrain = 0.0;
    trialStrain = 0.0;
    minElasticYieldStrain = gap;

    return 0;
}


UniaxialMaterial *
Neoprene::getCopy(void)
{
    Neoprene *theCopy = new Neoprene(this->getTag(),E,gap);
    theCopy->trialStrain = trialStrain;
    theCopy->maxElasticYieldStrain = maxElasticYieldStrain;
    theCopy->minElasticYieldStrain = minElasticYieldStrain;
    return theCopy;
}


int 
Neoprene::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(6);
  data(0) = this->getTag();
  data(1) = commitStrain;
  data(2) = E;
  data(3) = gap;
  data(4) = maxElasticYieldStrain;
  data(5) = minElasticYieldStrain;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Neoprene::sendSelf() - failed to send data\n";

  return res;
}

int 
Neoprene::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(6);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    opserr << "Neoprene::recvSelf() - failed to recv data\n";
  else {
    this->setTag((int)data(0));
    commitStrain = data(1);
    trialStrain = commitStrain;
    E = data(2);
    gap = data(3);
    maxElasticYieldStrain = data(4);
    minElasticYieldStrain = data(5);
  }

  return res;
}

void 
Neoprene::Print(OPS_Stream &s, int flag)
{
    s << "Neoprene tag: " << this->getTag() << endln;
    s << "  E: " << E << endln;
    s << "  initial gap: " << gap << endln;
}
