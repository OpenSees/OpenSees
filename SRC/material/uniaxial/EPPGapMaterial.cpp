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
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/EPPGapMaterial.cpp,v $

// File: ~/material/EPPGapMaterial.C
//
// Written: krm 
// Created: 07/2000
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) EPPGapMaterial.C, revA"


#include <EPPGapMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


EPPGapMaterial::EPPGapMaterial(int tag, double e, double fyl, double gap0)
:UniaxialMaterial(tag,MAT_TAG_EPPGap),
 commitStrain(0.0), trialStrain(0.0), E(e), fy(fyl), gap(gap0),
 minElasticYieldStrain(gap0)
{
	if (E == 0.0) {
		g3ErrorHandler->warning("%s -- E is zero, continuing with E = fy/0.002",
			"EPPGapMaterial::EPPGapMaterial");
		if (fy != 0.0)
			E = fabs(fy)/0.002;
		else
			g3ErrorHandler->fatal("%s -- E and fy are zero",
				"EPPGapMaterial::EPPGapMaterial");
	}
	else
		maxElasticYieldStrain = fy/E + gap;

	if (fy*gap<0) {
		g3ErrorHandler->warning("%s -- Alternate signs on fy and E encountered, continuing anyway",
			"EPPGapMaterial::EPPGapMaterial");
	}
}

EPPGapMaterial::EPPGapMaterial()
:UniaxialMaterial(0,MAT_TAG_EPPGap),
 E(0.0), fy(0.0), gap(0.0), minElasticYieldStrain(0.0)
{

}

EPPGapMaterial::~EPPGapMaterial()
{
  // does nothing
}

int 
EPPGapMaterial::setTrialStrain(double strain, double strainRate)
{
    trialStrain = strain;
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
    if (fy >= 0) {
       if (trialStrain >= maxElasticYieldStrain)
           return fy;
       else if (trialStrain <= minElasticYieldStrain)
           return 0;
       else
           return E*(trialStrain-minElasticYieldStrain);
    }
    else {
       if (trialStrain <= maxElasticYieldStrain)
           return fy;
       else if (trialStrain >= minElasticYieldStrain)
           return 0;
       else
           return E*(trialStrain-minElasticYieldStrain);
    }
}

double 
EPPGapMaterial::getTangent(void)
{
    if (fy >= 0) {
       if (trialStrain >= maxElasticYieldStrain)
           return 0;
       else if (trialStrain <= minElasticYieldStrain)
           return 0;
       else
           return E;
    }
    else {
       if (trialStrain <= maxElasticYieldStrain)
           return 0;
       else if (trialStrain >= minElasticYieldStrain)
           return 0;
       else
           return E;
    }
}

double 
EPPGapMaterial::getSecant ()
{
    if (fabs(trialStrain-minElasticYieldStrain) > DBL_EPSILON)
        return this->getStress()/(trialStrain-minElasticYieldStrain);
    else
        return 0;
}

int 
EPPGapMaterial::commitState(void)
{
    if (fy >= 0) {
       if (trialStrain > maxElasticYieldStrain)  {
           maxElasticYieldStrain = trialStrain;
           minElasticYieldStrain = trialStrain-fy/E;
       }
       else if (trialStrain < minElasticYieldStrain && trialStrain > gap)  {
           maxElasticYieldStrain = trialStrain+fy/E;
           minElasticYieldStrain = trialStrain;
       }
    }
    else {
       if (trialStrain < maxElasticYieldStrain)  {
           maxElasticYieldStrain = trialStrain;
           minElasticYieldStrain = trialStrain-fy/E;
       }
       else if (trialStrain > minElasticYieldStrain && trialStrain < gap)  {
           maxElasticYieldStrain = trialStrain+fy/E;
           minElasticYieldStrain = trialStrain;
       }
    }

    commitStrain = trialStrain;
    return 0;
}


int 
EPPGapMaterial::revertToLastCommit(void)
{
    trialStrain = commitStrain;
    return 0;
}


int 
EPPGapMaterial::revertToStart(void)
{
    commitStrain = 0.0;
    trialStrain = 0.0;
    maxElasticYieldStrain = fy/E+gap;
    minElasticYieldStrain = gap;

    return 0;
}


UniaxialMaterial *
EPPGapMaterial::getCopy(void)
{
    EPPGapMaterial *theCopy = new EPPGapMaterial(this->getTag(),E,fy,gap);
    theCopy->trialStrain = trialStrain;
    theCopy->maxElasticYieldStrain = maxElasticYieldStrain;
    theCopy->minElasticYieldStrain = minElasticYieldStrain;
    return theCopy;
}


int 
EPPGapMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(7);
  data(0) = this->getTag();
  data(1) = commitStrain;
  data(2) = E;
  data(3) = fy;
  data(4) = gap;
  data(5) = maxElasticYieldStrain;
  data(6) = minElasticYieldStrain;
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    cerr << "EPPGapMaterial::sendSelf() - failed to send data\n";

  return res;
}

int 
EPPGapMaterial::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(7);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0)
    cerr << "EPPGapMaterial::recvSelf() - failed to recv data\n";
  else {
    this->setTag((int)data(0));
    commitStrain = data(1);
    trialStrain = commitStrain;
    E = data(2);
    fy = data(3);
    gap = data(4);
    maxElasticYieldStrain = data(5);
    minElasticYieldStrain = data(6);
  }

  return res;
}

void 
EPPGapMaterial::Print(ostream &s, int flag)
{
    s << "EPPGap tag: " << this->getTag() << endl;
    s << "  E: " << E << endl;
    s << "  fy: " << fy << endl;
    s << "  initial gap: " << gap << endl;
}
