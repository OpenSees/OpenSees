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

// $Revision: 1.6 $
// $Date: 2006-08-15 00:41:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/NewUniaxialMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class implementation for 
// NewUniaxialMaterial.

#include <NewUniaxialMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

NewUniaxialMaterial::NewUniaxialMaterial(int tag)
  :UniaxialMaterial(tag,MAT_TAG_NewUniaxialMaterial),
   trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
   commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{

}

NewUniaxialMaterial::NewUniaxialMaterial()
  :UniaxialMaterial(0,MAT_TAG_NewUniaxialMaterial),
   trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
   commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{

}

NewUniaxialMaterial::~NewUniaxialMaterial()
{

}

int 
NewUniaxialMaterial::setTrialStrain(double strain, double strainRate)
{
  // set the trial strain
  trialStrain = strain;

  // determine trial stress and tangent
  trialStress = 0.0;
  trialTangent = 0.0;

  return 0;
}

double 
NewUniaxialMaterial::getStress(void)
{
  return trialStress;
}

double 
NewUniaxialMaterial::getTangent(void)
{
  return trialTangent;
}

double 
NewUniaxialMaterial::getInitialTangent(void)
{
  // return the initial tangent
  return 0.0;
}

double 
NewUniaxialMaterial::getStrain(void)
{
  return trialStrain;
}

int 
NewUniaxialMaterial::commitState(void)
{
  commitStrain  = trialStrain;
  commitStress  = trialStress;
  commitTangent = trialTangent;

  return 0;
}

int 
NewUniaxialMaterial::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialStress = commitStress;
  trialTangent = commitTangent;

  return 0;
}

int 
NewUniaxialMaterial::revertToStart(void)
{
  trialStrain = 0.;
  trialStress = 0.0;
  trialTangent = 0.0;
  commitStrain = 0.;
  commitStress = 0.0;
  commitTangent = 0.0;

  return 0;
}

UniaxialMaterial *
NewUniaxialMaterial::getCopy(void)
{
  NewUniaxialMaterial *theCopy = new NewUniaxialMaterial(this->getTag());

  return theCopy;
}

int 
NewUniaxialMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
NewUniaxialMaterial::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
NewUniaxialMaterial::Print(OPS_Stream &s, int flag)
{
  s << "NewUniaxialMaterial : " << this->getTag();

  return;
}


