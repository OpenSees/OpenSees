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

// $Revision: 1.5 $
// $Date: 2009-03-27 19:19:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/SecantConcrete.cpp,v $

// Written: MHS
// Created: Dec 2001
//
// Description: This file contains the class implementation for 
// SecantConcrete. 

#include <SecantConcrete.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <math.h>
#include <string.h>

SecantConcrete::SecantConcrete(int tag, double f, double ec, double eu)
  :UniaxialMaterial(tag,MAT_TAG_SecantConcrete),
   fc(f), epsc(ec), epsu(eu),
   parameterID(0), SHVs(0)
{
  // Initialize variables
  this->revertToStart();
}

SecantConcrete::SecantConcrete()
  :UniaxialMaterial(0,MAT_TAG_SecantConcrete),
   fc(0.0), epsc(0.0), epsu(0.0),
   parameterID(0), SHVs(0)
{
  // Initialize variables
  this->revertToStart();
}

SecantConcrete::~SecantConcrete()
{
  if (SHVs != 0)
    delete SHVs;
}

void
SecantConcrete::backbone(double strain, double &stress, double &tangent)
{
  if (strain > 0.0 || strain < epsu) {
    stress = 0.0;
    tangent = 0.0;
  }
  else if (strain > epsc) {
    double xi = strain/epsc;
    stress = fc*(2.0*xi-xi*xi);
    tangent = 2.0*fc/epsc*(1.0-xi);
  }
  else {
    tangent = -fc/(epsu-epsc);
    stress = tangent*(strain-epsu);
  }
}

int 
SecantConcrete::setTrialStrain(double strain, double strainRate)
{
  // Set total strain
  Tstrain = strain;
  
  if (Tstrain > 0.0 || Tstrain < epsu) {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
  else if (Tstrain > CminStrain) {
    double sigmaMin, dummy;
    this->backbone(CminStrain, sigmaMin, dummy);
    Ttangent = sigmaMin/CminStrain;
    Tstress = Ttangent*Tstrain;
  }
  else {
    this->backbone(Tstrain, Tstress, Ttangent);
    TminStrain = Tstrain;
  }
  
  return 0;
}

double 
SecantConcrete::getStress(void)
{
  return Tstress;
}

double 
SecantConcrete::getTangent(void)
{
  return Ttangent;
}

double 
SecantConcrete::getStrain(void)
{
  return Tstrain;
}

int 
SecantConcrete::commitState(void)
{
  // Commit trial history variables
  CminStrain = TminStrain;
  
  return 0;
}

int 
SecantConcrete::revertToLastCommit(void)
{
  // Nothing to do here
  return 0;
}

int 
SecantConcrete::revertToStart(void)
{
  // Reset committed history variables
  CminStrain = 0.0;
  
  // Reset trial history variables
  TminStrain = 0.0;
  
  // Initialize state variables
  Tstrain  = 0.0;
  Tstress  = 0.0;
  Ttangent = 0.0;
  
  // Reset sensitivity history variables
  if (SHVs != 0)
    SHVs->Zero();
  
  return 0;
}

UniaxialMaterial *
SecantConcrete::getCopy(void)
{
  SecantConcrete *theCopy = new SecantConcrete(this->getTag(), fc, epsc, epsu);
  
  // Copy committed history variables
  theCopy->CminStrain = CminStrain;
  
  // Copy trial history variables
  theCopy->TminStrain = TminStrain;
  
  // Copy trial state variables
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;
  
  return theCopy;
}

int 
SecantConcrete::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(5);
  
  data(0) = this->getTag();
  data(1) = fc;
  data(2) = epsc;
  data(3) = epsu;
  data(4) = CminStrain;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "SecantConcrete::sendSelf() - failed to send data\n";

  return res;
}

int 
SecantConcrete::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(5);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "SecantConcrete::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    fc = data(1);
    epsc = data(2);
    epsu = data(3);
    CminStrain = data(4);
  }
    
  return res;
}

void 
SecantConcrete::Print(OPS_Stream &s, int flag)
{
  s << "SecantConcrete, tag: " << this->getTag() << endln;
  s << "  fc: " << fc << endln;
  s << "  epsc: " << epsc << endln;
  s << "  epsu: " << epsu << endln;
}


int
SecantConcrete::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"fc") == 0 || strcmp(argv[0],"fpc") == 0) {
    return param.addObject(5, this);
  }
  else if (strcmp(argv[0],"ec") == 0 || strcmp(argv[0],"epsc") == 0) {
    return param.addObject(6, this);
  }
  else if (strcmp(argv[0],"eu") == 0 || strcmp(argv[0],"epsu") == 0) {
    return param.addObject(7, this);
  }

  return -1;
}

int 
SecantConcrete::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case -1:
    return -1;
  case 5:
    fc = info.theDouble;
    return 0;
  case 6:
    epsc = info.theDouble;
    return 0;
  case 7:
    epsu = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
SecantConcrete::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

double
SecantConcrete::getStressSensitivity(int gradIndex, bool conditional)
{
  return this->getStressGradient(gradIndex);
}

double
SecantConcrete::getInitialTangentSensitivity(int gradIndex)
{
  return 0.0;
}

int
SecantConcrete::commitSensitivity(double strainGradient,
				  int gradIndex, int numGrads)
{
  if (SHVs == 0) {
    SHVs = new Matrix(2,numGrads);
    SHVs->Zero();
  }

  return this->setStrainGradient(gradIndex, strainGradient);
}

double
SecantConcrete::backboneCondDeriv(double strain, int gradIndex)
{
  double dfcdh = 0.0;
  double decdh = 0.0;
  double deudh = 0.0;
  
  switch (parameterID) {
  case 5:
    dfcdh = 1.0; break;
  case 6:
    decdh = 1.0; break;
  case 7:
    deudh = 1.0; break;
  default:
    break;
  }
  
  if (strain > 0.0 || strain < epsu)
    return 0.0;
  else if (strain > epsc) {
    double xi = strain/epsc;
    return dfcdh*(2.0*xi-xi*xi) + decdh*2.0*fc/epsc*(xi*xi-xi);
  }
  else {
    double E = -fc/(epsu-epsc);
    double dEdh = -dfcdh/(epsu-epsc) + fc/((epsu-epsc)*(epsu-epsc))*(deudh-decdh);
    return dEdh*(strain-epsu) - E*deudh;
  }
}

double
SecantConcrete::getStressGradient(int gradIndex)
{
  double depsdh = 0.0;
  double dsigdh = 0.0;

  if (SHVs != 0) {
    depsdh = (*SHVs)(0,gradIndex);
    dsigdh = (*SHVs)(1,gradIndex);
  }

  if (Tstrain > 0.0 || Tstrain < epsu)
    return 0.0;
  else if (Tstrain > CminStrain) {
    double sigmaMin, dummy;
    this->backbone(CminStrain, sigmaMin, dummy);
    return Tstrain*(CminStrain*dsigdh-sigmaMin*depsdh)/(CminStrain*CminStrain);
  }
  else
    return this->backboneCondDeriv(Tstrain, gradIndex);
}

double
SecantConcrete::backboneUncondDeriv(double strain, int gradIndex,
				    double depsilondh)
{
  double dfcdh = 0.0;
  double decdh = 0.0;
  double deudh = 0.0;
  
  switch (parameterID) {
  case 5:
    dfcdh = 1.0; break;
  case 6:
    decdh = 1.0; break;
  case 7:
    deudh = 1.0; break;
  default:
    break;
  }
  
  if (strain > 0.0 || strain < epsu)
    return 0.0;
  else if (strain > epsc) {
    double xi = strain/epsc;
    return dfcdh*(2.0*xi-xi*xi) + decdh*2.0*fc/epsc*(xi*xi-xi) + 2.0*fc/epsc*(1.0-xi)*depsilondh;
  }
  else {
    double E = -fc/(epsu-epsc);
    double dEdh = -dfcdh/(epsu-epsc) + fc/((epsu-epsc)*(epsu-epsc))*(deudh-decdh);
    return dEdh*(strain-epsu) + E*(depsilondh-deudh);
  }
  
}

int
SecantConcrete::setStrainGradient(int gradIndex, double depsilondh)
{
  if (SHVs != 0) {

    double &depsdh = (*SHVs)(0,gradIndex);
    double &dsigdh = (*SHVs)(1,gradIndex);

    if (TminStrain < CminStrain) {
      depsdh = depsilondh;
      dsigdh = this->backboneUncondDeriv(TminStrain, gradIndex, depsdh);
    }
  }

  return 0;
}
