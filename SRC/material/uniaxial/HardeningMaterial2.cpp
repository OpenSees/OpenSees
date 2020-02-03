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

// $Revision: 1.8 $
// $Date: 2003/02/14 23:01:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HardeningMaterial2.cpp,v $

// Written: MHS
// Created: May 2000
// Revision: A
//
// Description: This file contains the class implementation for 
// HardeningMaterial2. 

#include <HardeningMaterial2.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>

HardeningMaterial2::HardeningMaterial2(int tag, double e, double s,
				     double hi, double hk, double n)
:UniaxialMaterial(tag,MAT_TAG_Hardening),
 E(e), sigmaY(s), Hiso(hi), Hkin(hk), eta(n),
 parameterID(0), SHVs(0)
{
  // Initialize variables
  this->revertToStart();
}

HardeningMaterial2::HardeningMaterial2()
:UniaxialMaterial(0,MAT_TAG_Hardening),
 E(0.0), sigmaY(0.0), Hiso(0.0), Hkin(0.0), eta(0.0),
 parameterID(0), SHVs(0)
{
  // Initialize variables
  this->revertToStart();
}

HardeningMaterial2::~HardeningMaterial2()
{
  if (SHVs != 0)
    delete SHVs;
}

int 
HardeningMaterial2::setTrialStrain (double strain, double strainRate)
{
  if (fabs(Tstrain-strain) < DBL_EPSILON)
    return 0;
  
  // Set total strain
  Tstrain = strain;
  
  // Elastic trial stress
  Tstress = E * (Tstrain-CplasticStrain);
  
  // Compute trial stress relative to committed back stress
  double xsi = Tstress - CbackStress;
  
  // Compute yield criterion
  double f = fabs(xsi) - (sigmaY + Hiso*Chardening);
  
  // Elastic step ... no updates required
  if (f <= -DBL_EPSILON * E) {
    // Set trial tangent
    Ttangent = E;
  }
  
  // Plastic step ... perform return mapping algorithm
  else {
    double etadt = 0.0;
    
    if (eta != 0.0 || ops_Dt != 0)
      etadt = eta/ops_Dt;
    
    // Compute consistency parameter
    double dGamma = f / (E+Hiso+Hkin+etadt);
    
    // Find sign of xsi
    int sign = (xsi < 0) ? -1 : 1;
    
    // Bring trial stress back to yield surface
    Tstress -= dGamma*E*sign;
    
    // Update plastic strain
    TplasticStrain = CplasticStrain + dGamma*sign;
    
    // Update back stress
    TbackStress = CbackStress + dGamma*Hkin*sign;
    
    // Update internal hardening variable
    Thardening = Chardening + dGamma;
    
    // Set trial tangent
    Ttangent = E*(Hkin+Hiso+etadt) / (E+Hkin+Hiso+etadt);
  }
  
  return 0;
}

double 
HardeningMaterial2::getStress(void)
{
  return Tstress;
}

double 
HardeningMaterial2::getTangent(void)
{
  return Ttangent;
}

double 
HardeningMaterial2::getStrain(void)
{
  return Tstrain;
}

int 
HardeningMaterial2::commitState(void)
{
  // Commit trial history variables
  CplasticStrain = TplasticStrain;
  CbackStress = TbackStress;
  Chardening = Thardening;

  return 0;
}

int 
HardeningMaterial2::revertToLastCommit(void)
{
  return 0;
}

int 
HardeningMaterial2::revertToStart(void)
{
  // Reset committed history variables
  CplasticStrain = 0.0;
  CbackStress = 0.0;
  Chardening = 0.0;
  
  // Reset trial history variables
  TplasticStrain = 0.0;
  TbackStress = 0.0;
  Thardening = 0.0;
  
  // Initialize state variables
  Tstrain = 0.0;
  Tstress = 0.0;
  Ttangent = E;
  
  // Reset sensitivity history variables
  if (SHVs != 0)
    SHVs->Zero();
  
  return 0;
}

UniaxialMaterial *
HardeningMaterial2::getCopy(void)
{
  HardeningMaterial2 *theCopy =
    new HardeningMaterial2(this->getTag(), E, sigmaY, Hiso, Hkin, eta);
  
  // Copy committed history variables
  theCopy->CplasticStrain = CplasticStrain;
  theCopy->CbackStress = CbackStress;
  theCopy->Chardening = Chardening;
  
  // Copy trial history variables
  theCopy->TplasticStrain = TplasticStrain;
  theCopy->TbackStress = TbackStress;
  theCopy->Thardening = Thardening;
  
  // Copy trial state variables
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;
  
  return theCopy;
}

int 
HardeningMaterial2::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(12);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = sigmaY;
  data(3) = Hiso;
  data(4) = Hkin;
  data(5) = eta;
  data(6) = CplasticStrain;
  data(7) = CbackStress;
  data(8) = Chardening;
  data(9) = Tstrain;
  data(10) = Tstress;
  data(11) = Ttangent;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "HardeningMaterial2::sendSelf() - failed to send data\n";

  return res;
}

int 
HardeningMaterial2::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(12);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "HardeningMaterial2::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    E = data(1);
    sigmaY = data(2);
    Hiso = data(3);
    Hkin = data(4);
    eta = data(5);
    CplasticStrain = data(6);
    CbackStress = data(7);
    Chardening = data(8);
    Tstrain = data(9);
    Tstress = data(10);
    Ttangent = data(11);
  }
    
  return res;
}

int
HardeningMaterial2::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return 0;

  if (strcmp(argv[0],"E") == 0 || strcmp(argv[0],"EA") == 0 || strcmp(argv[0],"EI") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0 || strcmp(argv[0],"Fy") == 0) {
    param.setValue(sigmaY);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Hiso") == 0 || strcmp(argv[0],"H_iso") == 0) {
    param.setValue(Hiso);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"Hkin") == 0 || strcmp(argv[0],"H_kin") == 0) {
    param.setValue(Hkin);
    return param.addObject(4, this);
  }

  return -1;
}

int 
HardeningMaterial2::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case -1:
    return -1;
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    sigmaY = info.theDouble;
    return 0;
  case 3:
    Hiso = info.theDouble;
    return 0;
  case 4:
    Hkin = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
HardeningMaterial2::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

double
HardeningMaterial2::getStressGradient(int gradIndex)
{
  // Derivatives of material parameters
  double dEdh = 0.0;
  double dsigmaYdh = 0.0;
  double dHisodh = 0.0;
  double dHkindh = 0.0;

  switch (parameterID) {
  case 1:
    dEdh = 1.0; break;
  case 2:
    dsigmaYdh = 1.0; break;
  case 3:
    dHisodh = 1.0; break;
  case 4:
    dHkindh = 1.0; break;
  default:
    /*return 0.0;*/	break;
  }
  
  double dPlasticStraindh = 0.0;
  double dBackStressdh    = 0.0;
  double dHardeningdh     = 0.0;

  if (SHVs != 0 && gradIndex < SHVs->noCols()) {
    dPlasticStraindh = (*SHVs)(0,gradIndex);
    dBackStressdh    = (*SHVs)(1,gradIndex);
    dHardeningdh     = (*SHVs)(2,gradIndex);
  }

  // Trial stress relative to backstress
  double xi = E*(Tstrain-CplasticStrain) - CbackStress;
  
  // Yield function
  double f = fabs(xi) - sigmaY - Hiso*Chardening;
  
  // Conditional stress derivative
  double dsigmadh = 0.0;
  
  // Elastic step
  if (f <= -DBL_EPSILON * E)
    // CONDITIONAL derivative of stress
    dsigmadh = -E*dPlasticStraindh + dEdh*(Tstrain-CplasticStrain);
  
  // Plastic step
  else {
    // CONDITIONAL derivative of xi
    double dxidh = -E*dPlasticStraindh + dEdh*(Tstrain-CplasticStrain) - dBackStressdh;
    
    // Vector normal to yield surface
    double n = (xi < 0.0) ? -1.0 : 1.0;
    
    // Derivative of yield function
    double dfdh = dxidh*n - dsigmaYdh - Hiso*dHardeningdh - dHisodh*Chardening;
    
    // Change in consistency parameter
    double invEHiHk = 1.0/(E+Hiso+Hkin);
    double dg = f*invEHiHk;
    
    // Derivative of change in consistency parameter
    double ddgdh = -dg*invEHiHk*(dEdh + dHisodh + dHkindh) + dfdh*invEHiHk;
    
    // CONDITIONAL derivative of stress
    dsigmadh = -E*(dPlasticStraindh + ddgdh*n) + dEdh*(Tstrain - CplasticStrain - dg*n);
  }
  
  return dsigmadh;
}

int
HardeningMaterial2::setStrainGradient(int gradIndex, double depsilondh)
{
  if (SHVs == 0)
    return 0;

  double &dPlasticStraindh = (*SHVs)(0,gradIndex);
  double &dBackStressdh    = (*SHVs)(1,gradIndex);
  double &dHardeningdh     = (*SHVs)(2,gradIndex);

  // Derivatives of material parameters
  double dEdh = 0.0;
  double dsigmaYdh = 0.0;
  double dHisodh = 0.0;
  double dHkindh = 0.0;
  
  switch (parameterID) {
  case 1:
    dEdh = 1.0; break;
  case 2:
    dsigmaYdh = 1.0; break;
  case 3:
    dHisodh = 1.0; break;
  case 4:
    dHkindh = 1.0; break;
  default:
    /*return 0;*/ break;
  }
  
  // Trial stress relative to backstress
  double xi = E*(Tstrain-CplasticStrain) - CbackStress;
  
  // Yield function
  double f = fabs(xi) - sigmaY - Hiso*Chardening;
  
  // Elastic step
  if (f <= -DBL_EPSILON * E) {
    // Do not need to update any sensitivity history variables
  }
  
  // Plastic step
  else {
    // UNCONDITIONAL derivative of xi
    double dxidh = E*(depsilondh-dPlasticStraindh) +
      dEdh*(Tstrain-CplasticStrain) - dBackStressdh;
    
    // Vector normal to yield surface
    double n = (xi < 0.0) ? -1.0 : 1.0;
    
    // Derivative of yield function
    double dfdh = dxidh*n - dsigmaYdh - Hiso*dHardeningdh - dHisodh*Chardening;
    
    // Change in consistency parameter
    double invEHiHk = 1.0/(E+Hiso+Hkin);
    double dg = f*invEHiHk;
    
    // Derivative of change in consistency parameter
    double ddgdh = -dg*invEHiHk*(dEdh + dHisodh + dHkindh) + dfdh*invEHiHk;
    
    // Update derivative of plastic strain
    dPlasticStraindh += ddgdh*n;
    
    // Update derivative of internal hardening variable
    dHardeningdh += ddgdh;
    
    // Update derivative of back stress
    dBackStressdh += dg*dHkindh*n + ddgdh*Hkin*n;
  }
  
  return 0;
}

void
HardeningMaterial2::Print(OPS_Stream &s, int flag)
{
  s << "HardeningMaterial2, tag: " << this->getTag() << endln;
  s << "  E: " << E << endln;
  s << "  sigmaY: " << sigmaY << endln;
  s << "  Hiso: " << Hiso << endln;
  s << "  Hkin: " << Hkin << endln;
  s << "  eta: " << eta << endln;
}

double
HardeningMaterial2::getStressSensitivity(int gradIndex, bool conditional)
{
  return this->getStressGradient(gradIndex);
}

double
HardeningMaterial2::getInitialTangentSensitivity(int gradIndex)
{
  if (parameterID == 1)
    return 1.0; 
  else
    return 0.0;
}

int
HardeningMaterial2::commitSensitivity(double strainGradient,
				      int gradIndex, int numGrads)
{
  if (SHVs == 0) {
    SHVs = new Matrix(3,numGrads);
    SHVs->Zero();
  }

  if (gradIndex >= SHVs->noCols())
    return 0;

  return this->setStrainGradient(gradIndex, strainGradient);
}
