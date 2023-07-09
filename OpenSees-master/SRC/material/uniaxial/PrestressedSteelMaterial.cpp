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
// $Date: 2007-10-26 04:32:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/PrestressedSteelMaterial.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the class implementation for 
// PrestressedSteelMaterial.

#include <PrestressedSteelMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

PrestressedSteelMaterial::PrestressedSteelMaterial(int tag, double e, double f, double epslp, double epsd)
  :UniaxialMaterial(tag,MAT_TAG_PrestressedSteelMaterial),
   trialStrain(0.0), E(e), fy(f), elp(epslp), ed(epsd)
{

}

PrestressedSteelMaterial::PrestressedSteelMaterial()
  :UniaxialMaterial(0,MAT_TAG_PrestressedSteelMaterial),
   trialStrain(0.0), E(0.0), fy(0.0), elp(0.0), ed(0.0)
{

}

PrestressedSteelMaterial::~PrestressedSteelMaterial()
{

}

int 
PrestressedSteelMaterial::setTrialStrain(double strain, double strainRate)
{
  // set the trial strain
  trialStrain = strain;

  return 0;
}

double 
PrestressedSteelMaterial::getStress(void)
{
  if (trialStrain <= elp)
    return E*trialStrain;
  else {
    return fy-0.04/(trialStrain-ed);
  }
}

double 
PrestressedSteelMaterial::getTangent(void)
{
  if (trialStrain <= elp)
    return E;
  else {
    double tmp = (trialStrain-ed);
    return 0.04/(tmp*tmp);
  }
}

double 
PrestressedSteelMaterial::getInitialTangent(void)
{
  return E;
}

double 
PrestressedSteelMaterial::getStrain(void)
{
  return trialStrain;
}

int 
PrestressedSteelMaterial::commitState(void)
{
  return 0;
}

int 
PrestressedSteelMaterial::revertToLastCommit(void)
{
  return 0;
}

int 
PrestressedSteelMaterial::revertToStart(void)
{
  trialStrain = 0.0;

  return 0;
}

UniaxialMaterial *
PrestressedSteelMaterial::getCopy(void)
{
  PrestressedSteelMaterial *theCopy = new PrestressedSteelMaterial(this->getTag(), E, fy, elp, ed);
  
  theCopy->trialStrain = trialStrain;

  return theCopy;
}

int 
PrestressedSteelMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
PrestressedSteelMaterial::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
PrestressedSteelMaterial::Print(OPS_Stream &s, int flag)
{
  s << "PrestressedSteelMaterial : " << this->getTag();

  return;
}


