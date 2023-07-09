/********************************************************************************
(C) Copyright 2001-2022, The Regents of the University of California    
All Rights Reserved.                                               

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************* */

// written: MHS, 2001

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


