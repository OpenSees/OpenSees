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

// $Revision: 1.3 $
// $Date: 2002-12-05 22:49:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2BeamFiber2d.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: Elastic isotropic model where stress 
// components 22, 33, 13, and 23 are condensed out.

#include <J2BeamFiber2d.h>           
#include <Channel.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>

Vector J2BeamFiber2d::sigma(2);
Matrix J2BeamFiber2d::D(2,2);

void *
OPS_J2BeamFiber2dMaterial(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 6) {
    opserr << "Want: nDMaterial J2BeamFiber $tag $E $v $sigmaY $Hiso $Hkin <$rho>" << endln;
    return 0;	
  }
  
  int iData[1];
  double dData[6];
  dData[5] = 0.0;
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial J2BeamFiber \n";
    return 0;
  }
  
  if (numArgs > 6) 
    numData = 6;
  else
    numData = 5;
  
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial J2BeamFiber : " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new J2BeamFiber2d(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]);
  
  return theMaterial;
}

J2BeamFiber2d::J2BeamFiber2d
(int tag, double e, double g, double sy, double hi, double hk):
  NDMaterial(tag, ND_TAG_J2BeamFiber2d),
  E(e), nu(g), sigmaY(sy), Hiso(hi), Hkin(hk),
  parameterID(0), SHVs(0), Tepsilon(2),
  alphan(0.0), alphan1(0.0), dg_n1(0.0)
{
  epsPn[0] = 0.0;
  epsPn[1] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
}

J2BeamFiber2d::J2BeamFiber2d():
  NDMaterial(0, ND_TAG_J2BeamFiber2d),
  E(0.0), nu(0.0), sigmaY(0.0), Hiso(0.0), Hkin(0.0), 
  parameterID(0), SHVs(0), Tepsilon(2), 
  alphan(0.0), alphan1(0.0), dg_n1(0.0)
{
  epsPn[0] = 0.0;
  epsPn[1] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
}

J2BeamFiber2d::~J2BeamFiber2d ()
{
  if (SHVs != 0)
    delete SHVs;
}

int
J2BeamFiber2d::setTrialStrain (const Vector &strain)
{
  Tepsilon = strain;

  return 0;
}

int
J2BeamFiber2d::setTrialStrain (const Vector &strain, const Vector &rate)
{
  Tepsilon = strain;

  return 0;
}

int
J2BeamFiber2d::setTrialStrainIncr (const Vector &strain)
{
  return 0;
}

int
J2BeamFiber2d::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return 0;
}

const Matrix&
J2BeamFiber2d::getTangent (void)
{
  double twoG = E/(1.0+nu);
  double G = 0.5*twoG;

  double sig[2];
  sig[0] = E*(Tepsilon(0)-epsPn[0]);
  sig[1] = G*(Tepsilon(1)-epsPn[1]);

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double two3Hkin = two3*Hkin;

  double xsi[2];
  //xsi[0] = sig[0] - two3*Hkin*1.5*epsPn[0];
  //xsi[1] = sig[1] - two3*Hkin*0.5*epsPn[1];
  xsi[0] = sig[0] -      Hkin*epsPn[0];
  xsi[1] = sig[1] - one3*Hkin*epsPn[1];

  double q = sqrt(two3*xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1]);
  double F = q - root23*(sigmaY + Hiso*alphan);

  if (F < -100*DBL_EPSILON) {
    D(0,0) = E;
    D(1,1) = G;
    D(0,1) = D(1,0) = 0.0;

    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
  }
  else {

    // Solve for dg
    double dg = 0.0;

    static Vector R(3);
    R(0) = 0.0; R(1) = 0.0; R(2) = F;
    static Vector x(3);
    x(0) = xsi[0]; x(1) = xsi[1]; x(2) = dg;

    static Matrix J(3,3);
    static Vector dx(3);

    int iter = 0; int maxIter = 25;
    while (iter < maxIter && R.Norm() > sigmaY*1.0e-14) {
        iter++;

        J(0,0) = 1.0 + dg*two3*(E+Hkin); J(0,1) = 0.0;
        J(1,0) = 0.0; J(1,1) = 1.0 + dg*(twoG+two3Hkin);

        J(0,2) = two3*(E+Hkin)*x(0);
        J(1,2) = (twoG+two3Hkin)*x(1);

        //J(2,0) = x(0)*two3/q; J(2,1) = x(1)*2.0/q;
        J(2,0) = (1.0-two3*Hiso*dg)*x(0)*two3/q;
        J(2,1) = (1.0-two3*Hiso*dg)*x(1)*2.0/q;

        //J(2,2) = -root23*Hiso;
	J(2,2) = -two3*Hiso*q;

        J.Solve(R, dx);
        x.addVector(1.0, dx, -1.0);

        dg = x(2);
        dg_n1 = dg;

        q = sqrt(two3*x(0)*x(0) + 2.0*x(1)*x(1));

        R(0) = x(0) - xsi[0] + dg*two3*(E+Hkin)*x(0);
        R(1) = x(1) - xsi[1] + dg*(twoG+two3Hkin)*x(1);
        R(2) = q - root23*(sigmaY + Hiso*(alphan+dg*root23*q));
    }

    if (iter == maxIter) {
      //opserr << "J2BeamFiber2d::getTangent -- maxIter reached " << R.Norm() << endln;
    }

    alphan1 = alphan + dg*root23*q;

    epsPn1[0] = epsPn[0] + dg*two3*x(0);
    epsPn1[1] = epsPn[1] + dg*2.0*x(1);

    //J(2,0) = (1.0-two3*Hiso*dg)*x(0)*two3/q; J(2,1) = (1.0-two3*Hiso*dg)*x(1)*2.0/q;
    //J(2,2) = -two3*Hiso*q;
    //static Matrix invJ(3,3);
    //J.Invert(invJ);

    J(0,0) = 1.0 + dg*two3*E/(1.0+dg*two3Hkin); J(0,1) = 0.0;
    J(1,0) = 0.0; J(1,1) = 1.0 + dg*twoG/(1.0+dg*two3Hkin);

    J(0,2) = (two3*E-dg*two3*E/(1.0+dg*two3Hkin)*two3Hkin)*x(0);
    J(1,2) = (twoG  -dg*  twoG/(1.0+dg*two3Hkin)*two3Hkin)*x(1);

    //J(2,0) = x(0)/q*two3/(1.0+dg*two3Hkin);
    //J(2,1) = x(1)/q* 2.0/(1.0+dg*two3Hkin);
    J(2,0) = (1.0-two3*Hiso*dg)*x(0)/q*two3/(1.0+dg*two3Hkin);
    J(2,1) = (1.0-two3*Hiso*dg)*x(1)/q* 2.0/(1.0+dg*two3Hkin);

    //J(2,2) = -(x(0)/q*two3/(1.0+dg*two3Hkin)*two3Hkin*x(0))
    //         -(x(1)/q* 2.0/(1.0+dg*two3Hkin)*two3Hkin*x(1));
    //J(2,2) = -q*two3Hkin/(1.0+dg*two3Hkin) - root23*Hiso;
    J(2,2) = -q*two3Hkin/(1.0+dg*two3Hkin) - two3*Hiso*q;

    static Matrix invJ(3,3);
    J.Invert(invJ);

    D(0,0) = invJ(0,0)*E;
    D(1,0) = invJ(1,0)*E;
    D(0,1) = invJ(0,1)*G;
    D(1,1) = invJ(1,1)*G;
  }

  return D;
}

const Matrix&
J2BeamFiber2d::getInitialTangent (void)
{
  double G = 0.5*E/(1.0+nu);

  D(0,0) = E;
  D(1,1) = G;
  D(0,1) = D(1,0) = 0.0;

  return D;
}

const Vector&
J2BeamFiber2d::getStress (void)
{
  double G = 0.5*E/(1.0+nu);

  sigma(0) = E*(Tepsilon(0)-epsPn[0]);
  sigma(1) = G*(Tepsilon(1)-epsPn[1]);

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double xsi[2];
  //xsi[0] = sigma(0) - two3*Hkin*1.5*epsPn[0];
  //xsi[1] = sigma(1) - two3*Hkin*0.5*epsPn[1];
  xsi[0] = sigma(0) -      Hkin*epsPn[0];
  xsi[1] = sigma(1) - one3*Hkin*epsPn[1];

  double q = sqrt(two3*xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1]);
  double F = q - root23*(sigmaY + Hiso*alphan);

  if (F < -100*DBL_EPSILON) {
    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
  }
  else {

    // Solve for dg
    double dg = 0.0;

    static Vector R(3);
    R(0) = 0.0; R(1) = 0.0; R(2) = F;
    static Vector x(3);
    x(0) = xsi[0]; x(1) = xsi[1]; x(2) = dg;

    static Matrix J(3,3);
    static Vector dx(3);

    int iter = 0; int maxIter = 25;
    while (iter < maxIter && R.Norm() > sigmaY*1.0e-14) {
        iter++;

        J(0,0) = 1.0 + dg*two3*(E+Hkin); J(0,1) = 0.0;
        J(1,0) = 0.0; J(1,1) = 1.0 + dg*(2.0*G+two3*Hkin);

        J(0,2) = two3*(E+Hkin)*x(0);
        J(1,2) = (2.0*G+two3*Hkin)*x(1);

        //J(2,0) = x(0)*two3/q; J(2,1) = x(1)*2.0/q;
        J(2,0) = (1.0-two3*Hiso*dg)*x(0)*two3/q;
        J(2,1) = (1.0-two3*Hiso*dg)*x(1)*2.0/q;

        //J(2,2) = -root23*Hiso;
	J(2,2) = -two3*Hiso*q;

        J.Solve(R, dx);
        x = x-dx;

        dg = x(2);
        dg_n1 = dg;

        q = sqrt(two3*x(0)*x(0) + 2.0*x(1)*x(1));

        R(0) = x(0) - xsi[0] + dg*two3*(E+Hkin)*x(0);
        R(1) = x(1) - xsi[1] + dg*(2.0*G+two3*Hkin)*x(1);
        R(2) = q - root23*(sigmaY + Hiso*(alphan+dg*root23*q));
    }

    if (iter == maxIter) {
      //opserr << "J2BeamFiber2d::getStress -- maxIter reached " << R.Norm() << endln;
    }

    alphan1 = alphan + dg*root23*q;

    epsPn1[0] = epsPn[0] + dg*two3*x(0);
    epsPn1[1] = epsPn[1] + dg*2.0*x(1);

    //sigma(0) = x(0) + two3*Hkin*1.5*epsPn1[0];
    //sigma(1) = x(1) + two3*Hkin*0.5*epsPn1[1];
    sigma(0) = x(0) +      Hkin*epsPn1[0];
    sigma(1) = x(1) + one3*Hkin*epsPn1[1];
  }

  return sigma;
}

const Vector&
J2BeamFiber2d::getStrain (void)
{
  return Tepsilon;
}

int
J2BeamFiber2d::commitState (void)
{
  epsPn[0] = epsPn1[0];
  epsPn[1] = epsPn1[1];

  alphan = alphan1;

  return 0;
}

int
J2BeamFiber2d::revertToLastCommit (void)
{
  epsPn1[0] = epsPn[0];
  epsPn1[1] = epsPn[1];

  alphan1 = alphan;

  return 0;
}

int
J2BeamFiber2d::revertToStart (void)
{
  Tepsilon.Zero();

  epsPn[0] = 0.0;
  epsPn[1] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;

  alphan = 0.0;
  alphan1 = 0.0;

  dg_n1 = 0.0;

  if (SHVs != 0)
    SHVs->Zero();

  return 0;
}

NDMaterial*
J2BeamFiber2d::getCopy (void)
{
  J2BeamFiber2d *theCopy =
    new J2BeamFiber2d (this->getTag(), E, nu, sigmaY, Hiso, Hkin);

  return theCopy;
}

NDMaterial*
J2BeamFiber2d::getCopy (const char *type)
{
  if (strcmp(type,this->getType()) == 0)
    return this->getCopy();

  return 0;
}

const char*
J2BeamFiber2d::getType (void) const
{
  return "BeamFiber2d";
}

int
J2BeamFiber2d::getOrder (void) const
{
  return 2;
}

const Vector&
J2BeamFiber2d::getStressSensitivity(int gradIndex, bool conditional)
{
  static Vector sigma(2);

  sigma(0) = 0.0;
  sigma(1) = 0.0;

  double dEdh = 0.0;
  double dsigmaYdh = 0.0;
  double dHkindh = 0.0;
  double dHisodh = 0.0;
  double dGdh = 0.0;

  if (parameterID == 1) { // E
    dEdh = 1.0;
    dGdh = 0.5/(1.0+nu);
  }
  if (parameterID == 2) { // nu
    dGdh = -0.5*E/(1.0 + 2.0*nu + nu*nu);
  }
  if (parameterID == 5) {
    dsigmaYdh = 1.0;
  }
  if (parameterID == 6) {
    dHkindh = 1.0;
  }
  if (parameterID == 7) {
    dHisodh = 1.0;
  }

  double G = 0.5*E/(1.0+nu);

  double depsPdh[2]; depsPdh[0] = 0.0; depsPdh[1] = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    depsPdh[0] = (*SHVs)(0,gradIndex);
    depsPdh[1] = (*SHVs)(1,gradIndex);
    dalphadh = (*SHVs)(2,gradIndex);
  }

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double xsi[2];
  xsi[0] = E*(Tepsilon(0)-epsPn1[0]) -      Hkin*epsPn1[0];
  xsi[1] = G*(Tepsilon(1)-epsPn1[1]) - one3*Hkin*epsPn1[1];

  double q = sqrt(two3*xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1]);
  double F = q - root23*(sigmaY + Hiso*alphan1);

  if (F <= -100*DBL_EPSILON) {
    sigma(0) = dEdh*(Tepsilon(0)-epsPn1[0]) - E*depsPdh[0];
    sigma(1) = dGdh*(Tepsilon(1)-epsPn1[1]) - G*depsPdh[1];
  }
  else {
    static Matrix J(3,3);
    static Vector b(3);
    static Vector dx(3);

    double dg = dg_n1;

    J(0,0) = 1.0 + dg*two3*(E+Hkin); J(0,1) = 0.0;
    J(1,0) = 0.0; J(1,1) = 1.0 + dg*(2.0*G+two3*Hkin);
  
    J(0,2) = two3*(E+Hkin)*xsi[0];
    J(1,2) = (2.0*G+two3*Hkin)*xsi[1];
    
    //J(2,0) = xsi[0]*two3/q; J(2,1) = xsi[1]*2.0/q;
    J(2,0) = (1.0-two3*Hiso*dg)*xsi[0]*two3/q; 
    J(2,1) = (1.0-two3*Hiso*dg)*xsi[1]*2.0/q;

    //J(2,2) = -root23*Hiso;
    J(2,2) = -two3*Hiso*q;

    b(0) = dEdh*Tepsilon(0) - (E+     Hkin)*depsPdh[0] - (dEdh+     dHkindh)*epsPn1[0];
    b(1) = dGdh*Tepsilon(1) - (G+one3*Hkin)*depsPdh[1] - (dGdh+one3*dHkindh)*epsPn1[1];
    b(2) = root23*(dsigmaYdh + dHisodh*alphan1 + Hiso*dalphadh);

    J.Solve(b, dx);

    depsPdh[0] += dx(2)*two3*xsi[0] + dg*two3*dx(0);
    depsPdh[1] += dx(2)* 2.0*xsi[1] + dg* 2.0*dx(1);

    sigma(0) = dx(0) +      Hkin*depsPdh[0] +      dHkindh*epsPn1[0];
    sigma(1) = dx(1) + one3*Hkin*depsPdh[1] + one3*dHkindh*epsPn1[1];
  }

  return sigma;
}

int
J2BeamFiber2d::commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads)
{
  if (SHVs == 0) {
    SHVs = new Matrix(3,numGrads);
  }

  if (gradIndex >= SHVs->noCols()) {
    //opserr << gradIndex << ' ' << SHVs->noCols() << endln;
    return 0;
  }
  
  //return 0;

  double dEdh = 0.0;
  double dsigmaYdh = 0.0;
  double dHkindh = 0.0;
  double dHisodh = 0.0;
  double dGdh = 0.0;

  if (parameterID == 1) { // E
    dEdh = 1.0;
    dGdh = 0.5/(1.0+nu);
  }
  if (parameterID == 2) { // nu
    dGdh = -0.5*E/(1.0 + 2.0*nu + nu*nu);
  }
  if (parameterID == 5) {
    dsigmaYdh = 1.0;
  }
  if (parameterID == 6) {
    dHkindh = 1.0;
  }
  if (parameterID == 7) {
    dHisodh = 1.0;
  }

  double G = 0.5*E/(1.0+nu);

  double depsPdh[2]; depsPdh[0] = 0.0; depsPdh[1] = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    depsPdh[0] = (*SHVs)(0,gradIndex);
    depsPdh[1] = (*SHVs)(1,gradIndex);
    dalphadh = (*SHVs)(2,gradIndex);
  }

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double xsi[2];
  xsi[0] = E*(Tepsilon(0)-epsPn1[0]) -      Hkin*epsPn1[0];
  xsi[1] = G*(Tepsilon(1)-epsPn1[1]) - one3*Hkin*epsPn1[1];

  double q = sqrt(two3*xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1]);
  double F = q - root23*(sigmaY + Hiso*alphan1);

  if (F <= -100*DBL_EPSILON) {
    // Do nothing
  }
  else {
    static Matrix J(3,3);
    static Vector b(3);
    static Vector dx(3);

    double dg = dg_n1;

    J(0,0) = 1.0 + dg*two3*(E+Hkin); J(0,1) = 0.0;
    J(1,0) = 0.0; J(1,1) = 1.0 + dg*(2.0*G+two3*Hkin);
  
    J(0,2) = two3*(E+Hkin)*xsi[0];
    J(1,2) = (2.0*G+two3*Hkin)*xsi[1];
    
    //J(2,0) = xsi[0]*two3/q; J(2,1) = xsi[1]*2.0/q;
    J(2,0) = (1.0-two3*Hiso*dg)*xsi[0]*two3/q; 
    J(2,1) = (1.0-two3*Hiso*dg)*xsi[1]*2.0/q;

    //J(2,2) = -root23*Hiso;
    J(2,2) = -two3*Hiso*q;

    b(0) = E*depsdh(0) + dEdh*Tepsilon(0) - (E+     Hkin)*depsPdh[0] - (dEdh+     dHkindh)*epsPn1[0];
    b(1) = G*depsdh(1) + dGdh*Tepsilon(1) - (G+one3*Hkin)*depsPdh[1] - (dGdh+one3*dHkindh)*epsPn1[1];
    b(2) = root23*(dsigmaYdh + dHisodh*alphan1 + Hiso*dalphadh);

    J.Solve(b, dx);

    dalphadh += dx(2)*root23*q + dg*root23*(xsi[0]*two3*dx(0)+xsi[1]*2.0*dx(1))/q;
    depsPdh[0] += dx(2)*two3*xsi[0] + dg*two3*dx(0);
    depsPdh[1] += dx(2)* 2.0*xsi[1] + dg* 2.0*dx(1);

    (*SHVs)(0,gradIndex) = depsPdh[0];
    (*SHVs)(1,gradIndex) = depsPdh[1];
    (*SHVs)(2,gradIndex) = dalphadh;
  }

  return 0;
}


int
J2BeamFiber2d::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(6+3);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = nu;
  data(3) = sigmaY;
  data(4) = Hiso;
  data(5) = Hkin;
  data(6) = epsPn[0];
  data(7) = epsPn[1];  
  data(8) = alphan;
    
  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "J2BeamFiber2d::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

int
J2BeamFiber2d::recvSelf (int commitTag, Channel &theChannel, 
				    FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(6+3);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
   opserr << "J2BeamFiber2d::recvSelf -- could not recv Vector\n";
   return res;
  }
    
  this->setTag((int)data(0));
  E = data(1);
  nu = data(2);
  sigmaY = data(3);
  Hiso = data(4);
  Hkin = data(5);
  epsPn[0] = data(6);
  epsPn[1] = data(7);
  alphan = data(8);
  
  return res;
}

void
J2BeamFiber2d::Print (OPS_Stream &s, int flag)
{
  s << "J2 Beam Fiber Material Model" << endln;
  s << "\tE:  " << E << endln;
  s << "\tnu:  " << nu << endln;
  s << "\tsigmaY:  " << sigmaY << endln;
  s << "\tHiso:  " << Hiso << endln;
  s << "\tHkin:  " << Hkin << endln;
  
  return;
}

int
J2BeamFiber2d::setParameter(const char **argv, int argc,
			    Parameter &param)
{
  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"nu") == 0) {
    param.setValue(nu);
    return param.addObject(2, this);  
  }
  else if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0  || strcmp(argv[0],"Fy") == 0) {
    param.setValue(sigmaY);
    return param.addObject(5, this);
  }
  else if (strcmp(argv[0],"Hkin") == 0) {
    param.setValue(Hkin);
    return param.addObject(6, this);
  }
  else if (strcmp(argv[0],"Hiso") == 0) {
    param.setValue(Hiso);
    return param.addObject(7, this);
  }

  return -1;
}

int 
J2BeamFiber2d::updateParameter(int parameterID, Information &info)
{ 
  switch(parameterID) {
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    nu = info.theDouble;
    return 0;
  case 5:
    sigmaY = info.theDouble;
    return 0;
  case 6:
    Hkin = info.theDouble;
    return 0;
  case 7:
    Hiso = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
J2BeamFiber2d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}
