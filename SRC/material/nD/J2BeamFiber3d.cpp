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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2BeamFiber3d.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: Elastic isotropic model where stress 
// components 22, 33, 13, and 23 are condensed out.

#include <J2BeamFiber3d.h>           
#include <Channel.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>

Vector J2BeamFiber3d::sigma(3);
Matrix J2BeamFiber3d::D(3,3);

void * OPS_ADD_RUNTIME_VPV(OPS_J2BeamFiber3dMaterial)
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
  
  theMaterial = new J2BeamFiber3d(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]);
  
  return theMaterial;
}

J2BeamFiber3d::J2BeamFiber3d
(int tag, double e, double g, double sy, double hi, double hk):
  NDMaterial(tag, ND_TAG_J2BeamFiber3d),
  E(e), nu(g), sigmaY(sy), Hiso(hi), Hkin(hk),
  parameterID(0), SHVs(0), Tepsilon(3),
  alphan(0.0), alphan1(0.0), dg_n1(0.0)
{
  epsPn[0] = 0.0;
  epsPn[1] = 0.0;
  epsPn[2] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
  epsPn1[2] = 0.0;
}

J2BeamFiber3d::J2BeamFiber3d():
  NDMaterial(0, ND_TAG_J2BeamFiber3d),
  E(0.0), nu(0.0), sigmaY(0.0), Hkin(0.0), 
  parameterID(0), SHVs(0), Tepsilon(3), 
  alphan(0.0), alphan1(0.0), dg_n1(0.0)
{
  epsPn[0] = 0.0;
  epsPn[1] = 0.0;
  epsPn[2] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
  epsPn1[2] = 0.0;
}

J2BeamFiber3d::~J2BeamFiber3d ()
{
  if (SHVs != 0)
    delete SHVs;
}

int
J2BeamFiber3d::setTrialStrain (const Vector &strain)
{
  Tepsilon = strain;

  return 0;
}

int
J2BeamFiber3d::setTrialStrain (const Vector &strain, const Vector &rate)
{
  Tepsilon = strain;

  return 0;
}

int
J2BeamFiber3d::setTrialStrainIncr (const Vector &strain)
{
  return 0;
}

int
J2BeamFiber3d::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return 0;
}

const Matrix&
J2BeamFiber3d::getTangent (void)
{
  double twoG = E/(1.0+nu);
  double G = 0.5*twoG;

  double sig[3];
  sig[0] = E*(Tepsilon(0)-epsPn[0]);
  sig[1] = G*(Tepsilon(1)-epsPn[1]);
  sig[2] = G*(Tepsilon(2)-epsPn[2]);

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double two3Hkin = two3*Hkin;

  double xsi[3];
  //xsi[0] = sig[0] - two3*Hkin*1.5*epsPn[0];
  //xsi[1] = sig[1] - two3*Hkin*0.5*epsPn[1];
  xsi[0] = sig[0] -      Hkin*epsPn[0];
  xsi[1] = sig[1] - one3*Hkin*epsPn[1];
  xsi[2] = sig[2] - one3*Hkin*epsPn[2];

  double q = sqrt(two3*xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1] + 2.0*xsi[2]*xsi[2]);
  double F = q - root23*(sigmaY + Hiso*alphan);

  if (F < -100*DBL_EPSILON) {
    D(0,0) = E;
    D(1,1) = G;
    D(2,2) = G;
    D(0,1) = D(1,0) = 0.0;
    D(0,2) = D(2,0) = 0.0;
    D(1,2) = D(2,1) = 0.0;

    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
    epsPn1[2] = epsPn[2];
    alphan1 = alphan;
  }
  else {

    // Solve for dg
    double dg = 0.0;

    static Vector R(4);
    R(0) = 0.0; R(1) = 0.0; R(2) = 0.0; R(3) = F;
    static Vector x(4);
    x(0) = xsi[0]; x(1) = xsi[1]; x(2) = xsi[2]; x(3) = dg;

    static Matrix J(4,4);
    static Vector dx(4);

    int iter = 0; int maxIter = 25;
    while (iter < maxIter && R.Norm() > sigmaY*1.0e-14) {
        iter++;

        J(0,0) = 1.0 + dg*two3*(E+Hkin); J(0,1) = 0.0; J(0,2) = 0.0;
        J(1,0) = 0.0; J(1,1) = 1.0 + dg*(twoG+two3Hkin); J(1,2) = 0.0;
        J(2,0) = 0.0; J(2,1) = 0.0; J(2,2) = 1.0 + dg*(twoG+two3Hkin);

        J(0,3) = two3*(E+Hkin)*x(0);
        J(1,3) = (twoG+two3Hkin)*x(1);
        J(2,3) = (twoG+two3Hkin)*x(2);

        //J(2,0) = x(0)*two3/q; J(2,1) = x(1)*2.0/q;
        J(3,0) = (1.0-two3*Hiso*dg)*x(0)*two3/q;
        J(3,1) = (1.0-two3*Hiso*dg)*x(1)*2.0/q;
        J(3,2) = (1.0-two3*Hiso*dg)*x(2)*2.0/q;

        //J(2,2) = -root23*Hiso;
	J(3,3) = -two3*Hiso*q;

        J.Solve(R, dx);
        x.addVector(1.0, dx, -1.0);

        dg = x(3);
        dg_n1 = dg;

        q = sqrt(two3*x(0)*x(0) + 2.0*x(1)*x(1) + 2.0*x(2)*x(2));

        R(0) = x(0) - xsi[0] + dg*two3*(E+Hkin)*x(0);
        R(1) = x(1) - xsi[1] + dg*(twoG+two3Hkin)*x(1);
        R(2) = x(2) - xsi[2] + dg*(twoG+two3Hkin)*x(2);
        R(3) = q - root23*(sigmaY + Hiso*(alphan+dg*root23*q));
    }

    if (iter == maxIter) {
      //opserr << "J2BeamFiber3d::getTangent -- maxIter reached " << R.Norm() << endln;
    }

    alphan1 = alphan + dg*root23*q;

    epsPn1[0] = epsPn[0] + dg*two3*x(0);
    epsPn1[1] = epsPn[1] + dg*2.0*x(1);
    epsPn1[2] = epsPn[2] + dg*2.0*x(2);

    //J(2,0) = (1.0-two3*Hiso*dg)*x(0)*two3/q; J(2,1) = (1.0-two3*Hiso*dg)*x(1)*2.0/q;
    //J(2,2) = -two3*Hiso*q;
    //static Matrix invJ(3,3);
    //J.Invert(invJ);

    J(0,0) = 1.0 + dg*two3*E/(1.0+dg*two3Hkin); J(0,1) = 0.0; J(0,2) = 0.0;
    J(1,0) = 0.0; J(1,1) = 1.0 + dg*twoG/(1.0+dg*two3Hkin); J(1,2) = 0.0;
    J(2,0) = 0.0; J(2,1) = 0.0; J(2,2) = 1.0 + dg*twoG/(1.0+dg*two3Hkin);

    J(0,3) = (two3*E-dg*two3*E/(1.0+dg*two3Hkin)*two3Hkin)*x(0);
    J(1,3) = (twoG  -dg*  twoG/(1.0+dg*two3Hkin)*two3Hkin)*x(1);
    J(2,3) = (twoG  -dg*  twoG/(1.0+dg*two3Hkin)*two3Hkin)*x(2);

    //J(2,0) = x(0)/q*two3/(1.0+dg*two3Hkin);
    //J(2,1) = x(1)/q* 2.0/(1.0+dg*two3Hkin);
    J(3,0) = (1.0-two3*Hiso*dg)*x(0)/q*two3/(1.0+dg*two3Hkin);
    J(3,1) = (1.0-two3*Hiso*dg)*x(1)/q* 2.0/(1.0+dg*two3Hkin);
    J(3,2) = (1.0-two3*Hiso*dg)*x(2)/q* 2.0/(1.0+dg*two3Hkin);

    //J(2,2) = -(x(0)/q*two3/(1.0+dg*two3Hkin)*two3Hkin*x(0))
    //         -(x(1)/q* 2.0/(1.0+dg*two3Hkin)*two3Hkin*x(1));
    //J(2,2) = -q*two3Hkin/(1.0+dg*two3Hkin) - root23*Hiso;
    J(3,3) = -q*two3Hkin/(1.0+dg*two3Hkin) - two3*Hiso*q;

    static Matrix invJ(4,4);
    J.Invert(invJ);

    D(0,0) = invJ(0,0)*E;
    D(1,0) = invJ(1,0)*E;
    D(2,0) = invJ(2,0)*E;
    D(0,1) = invJ(0,1)*G;
    D(1,1) = invJ(1,1)*G;
    D(2,1) = invJ(2,1)*G;
    D(0,2) = invJ(0,2)*G;
    D(1,2) = invJ(1,2)*G;
    D(2,2) = invJ(2,2)*G;
  }

  return D;
}

const Matrix&
J2BeamFiber3d::getInitialTangent (void)
{
  double G = 0.5*E/(1.0+nu);

  D(0,0) = E;
  D(1,1) = G;
  D(2,2) = G;
  D(0,1) = D(1,0) = 0.0;
  D(0,2) = D(2,0) = 0.0;
  D(1,2) = D(2,1) = 0.0;

  return D;
}

const Vector&
J2BeamFiber3d::getStress (void)
{
  double G = 0.5*E/(1.0+nu);

  sigma(0) = E*(Tepsilon(0)-epsPn[0]);
  sigma(1) = G*(Tepsilon(1)-epsPn[1]);
  sigma(2) = G*(Tepsilon(2)-epsPn[2]);

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double xsi[3];
  //xsi[0] = sigma(0) - two3*Hkin*1.5*epsPn[0];
  //xsi[1] = sigma(1) - two3*Hkin*0.5*epsPn[1];
  xsi[0] = sigma(0) -      Hkin*epsPn[0];
  xsi[1] = sigma(1) - one3*Hkin*epsPn[1];
  xsi[2] = sigma(2) - one3*Hkin*epsPn[2];

  double q = sqrt(two3*xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1] + 2.0*xsi[2]*xsi[2]);
  double F = q - root23*(sigmaY + Hiso*alphan);

  if (F < -100*DBL_EPSILON) {
    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
    epsPn1[2] = epsPn[2];
    alphan1 = alphan;
  }
  else {

    // Solve for dg
    double dg = 0.0;

    static Vector R(4);
    R(0) = 0.0; R(1) = 0.0; R(2) = 0.0; R(3) = F;
    static Vector x(4);
    x(0) = xsi[0]; x(1) = xsi[1]; x(2) = xsi[2]; x(3) = dg;

    static Matrix J(4,4);
    static Vector dx(4);

    int iter = 0; int maxIter = 25;
    while (iter < maxIter && R.Norm() > sigmaY*1.0e-14) {
        iter++;

        J(0,0) = 1.0 + dg*two3*(E+Hkin); J(0,1) = 0.0; J(0,2) = 0.0;
        J(1,0) = 0.0; J(1,1) = 1.0 + dg*(2.0*G+two3*Hkin); J(1,2) = 0.0;
        J(2,0) = 0.0; J(2,1) = 0.0; J(2,2) = 1.0 + dg*(2.0*G+two3*Hkin);

        J(0,3) = two3*(E+Hkin)*x(0);
        J(1,3) = (2.0*G+two3*Hkin)*x(1);
        J(2,3) = (2.0*G+two3*Hkin)*x(2);

        //J(2,0) = x(0)*two3/q; J(2,1) = x(1)*2.0/q;
        J(3,0) = (1.0-two3*Hiso*dg)*x(0)*two3/q;
        J(3,1) = (1.0-two3*Hiso*dg)*x(1)*2.0/q;
        J(3,2) = (1.0-two3*Hiso*dg)*x(2)*2.0/q;

        //J(2,2) = -root23*Hiso;
	J(3,3) = -two3*Hiso*q;

        J.Solve(R, dx);
        x = x-dx;

        dg = x(3);
        dg_n1 = dg;

        q = sqrt(two3*x(0)*x(0) + 2.0*x(1)*x(1) + 2.0*x(2)*x(2));

        R(0) = x(0) - xsi[0] + dg*two3*(E+Hkin)*x(0);
        R(1) = x(1) - xsi[1] + dg*(2.0*G+two3*Hkin)*x(1);
        R(2) = x(2) - xsi[2] + dg*(2.0*G+two3*Hkin)*x(2);
        R(3) = q - root23*(sigmaY + Hiso*(alphan+dg*root23*q));
    }

    if (iter == maxIter) {
      //opserr << "J2BeamFiber3d::getStress -- maxIter reached " << R.Norm() << endln;
    }

    alphan1 = alphan + dg*root23*q;

    epsPn1[0] = epsPn[0] + dg*two3*x(0);
    epsPn1[1] = epsPn[1] + dg*2.0*x(1);
    epsPn1[2] = epsPn[2] + dg*2.0*x(2);

    //sigma(0) = x(0) + two3*Hkin*1.5*epsPn1[0];
    //sigma(1) = x(1) + two3*Hkin*0.5*epsPn1[1];
    sigma(0) = x(0) +      Hkin*epsPn1[0];
    sigma(1) = x(1) + one3*Hkin*epsPn1[1];
    sigma(2) = x(2) + one3*Hkin*epsPn1[2];
  }

  return sigma;
}

const Vector&
J2BeamFiber3d::getStrain (void)
{
  return Tepsilon;
}

int
J2BeamFiber3d::commitState (void)
{
  epsPn[0] = epsPn1[0];
  epsPn[1] = epsPn1[1];
  epsPn[2] = epsPn1[2];

  alphan = alphan1;

  return 0;
}

int
J2BeamFiber3d::revertToLastCommit (void)
{
  epsPn1[0] = epsPn[0];
  epsPn1[1] = epsPn[1];
  epsPn1[2] = epsPn[2];

  alphan1 = alphan;

  return 0;
}

int
J2BeamFiber3d::revertToStart (void)
{
  Tepsilon.Zero();

  epsPn[0] = 0.0;
  epsPn[1] = 0.0;
  epsPn[2] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
  epsPn1[2] = 0.0;

  alphan = 0.0;
  alphan1 = 0.0;

  dg_n1 = 0.0;

  if (SHVs != 0)
    SHVs->Zero();

  return 0;
}

NDMaterial*
J2BeamFiber3d::getCopy (void)
{
  J2BeamFiber3d *theCopy =
    new J2BeamFiber3d (this->getTag(), E, nu, sigmaY, Hiso, Hkin);

  return theCopy;
}

NDMaterial*
J2BeamFiber3d::getCopy (const char *type)
{
  if (strcmp(type,this->getType()) == 0)
    return this->getCopy();

  return 0;
}

const char*
J2BeamFiber3d::getType (void) const
{
  return "BeamFiber";
}

int
J2BeamFiber3d::getOrder (void) const
{
  return 3;
}

const Vector&
J2BeamFiber3d::getStressSensitivity(int gradIndex, bool conditional)
{
  static Vector sigma(3);

  sigma(0) = 0.0;
  sigma(1) = 0.0;
  sigma(2) = 0.0;

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

  double depsPdh[3]; depsPdh[0] = 0.0; depsPdh[1] = 0.0; depsPdh[2] = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    depsPdh[0] = (*SHVs)(0,gradIndex);
    depsPdh[1] = (*SHVs)(1,gradIndex);
    depsPdh[2] = (*SHVs)(2,gradIndex);
    dalphadh = (*SHVs)(3,gradIndex);
  }

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double xsi[3];
  xsi[0] = E*(Tepsilon(0)-epsPn1[0]) -      Hkin*epsPn1[0];
  xsi[1] = G*(Tepsilon(1)-epsPn1[1]) - one3*Hkin*epsPn1[1];
  xsi[2] = G*(Tepsilon(2)-epsPn1[2]) - one3*Hkin*epsPn1[2];

  double q = sqrt(two3*xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1] + 2.0*xsi[2]*xsi[2]);
  double F = q - root23*(sigmaY + Hiso*alphan1);

  if (F <= -100*DBL_EPSILON) {
    sigma(0) = dEdh*(Tepsilon(0)-epsPn1[0]) - E*depsPdh[0];
    sigma(1) = dGdh*(Tepsilon(1)-epsPn1[1]) - G*depsPdh[1];
    sigma(2) = dGdh*(Tepsilon(2)-epsPn1[2]) - G*depsPdh[2];
  }
  else {
    static Matrix J(4,4);
    static Vector b(4);
    static Vector dx(4);

    double dg = dg_n1;

    J(0,0) = 1.0 + dg*two3*(E+Hkin); J(0,1) = 0.0; J(0,2) = 0.0;
    J(1,0) = 0.0; J(1,1) = 1.0 + dg*(2.0*G+two3*Hkin); J(1,2) = 0.0;
    J(2,0) = 0.0; J(2,1) = 0.0; J(2,2) = 1.0 + dg*(2.0*G+two3*Hkin);
  
    J(0,3) = two3*(E+Hkin)*xsi[0];
    J(1,3) = (2.0*G+two3*Hkin)*xsi[1];
    J(2,3) = (2.0*G+two3*Hkin)*xsi[2];
    
    //J(2,0) = xsi[0]*two3/q; J(2,1) = xsi[1]*2.0/q;
    J(3,0) = (1.0-two3*Hiso*dg)*xsi[0]*two3/q; 
    J(3,1) = (1.0-two3*Hiso*dg)*xsi[1]*2.0/q;
    J(3,2) = (1.0-two3*Hiso*dg)*xsi[2]*2.0/q;

    //J(2,2) = -root23*Hiso;
    J(3,3) = -two3*Hiso*q;

    b(0) = dEdh*Tepsilon(0) - (E+     Hkin)*depsPdh[0] - (dEdh+     dHkindh)*epsPn1[0];
    b(1) = dGdh*Tepsilon(1) - (G+one3*Hkin)*depsPdh[1] - (dGdh+one3*dHkindh)*epsPn1[1];
    b(2) = dGdh*Tepsilon(2) - (G+one3*Hkin)*depsPdh[2] - (dGdh+one3*dHkindh)*epsPn1[2];
    b(3) = root23*(dsigmaYdh + dHisodh*alphan1 + Hiso*dalphadh);

    J.Solve(b, dx);

    depsPdh[0] += dx(3)*two3*xsi[0] + dg*two3*dx(0);
    depsPdh[1] += dx(3)* 2.0*xsi[1] + dg* 2.0*dx(1);
    depsPdh[2] += dx(3)* 2.0*xsi[2] + dg* 2.0*dx(2);

    sigma(0) = dx(0) +      Hkin*depsPdh[0] +      dHkindh*epsPn1[0];
    sigma(1) = dx(1) + one3*Hkin*depsPdh[1] + one3*dHkindh*epsPn1[1];
    sigma(2) = dx(2) + one3*Hkin*depsPdh[2] + one3*dHkindh*epsPn1[2];
  }

  return sigma;
}

int
J2BeamFiber3d::commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads)
{
  if (SHVs == 0) {
    SHVs = new Matrix(4,numGrads);
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

  double depsPdh[3]; depsPdh[0] = 0.0; depsPdh[1] = 0.0; depsPdh[2] = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    depsPdh[0] = (*SHVs)(0,gradIndex);
    depsPdh[1] = (*SHVs)(1,gradIndex);
    depsPdh[2] = (*SHVs)(2,gradIndex);
    dalphadh = (*SHVs)(3,gradIndex);
  }

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double xsi[3];
  xsi[0] = E*(Tepsilon(0)-epsPn1[0]) -      Hkin*epsPn1[0];
  xsi[1] = G*(Tepsilon(1)-epsPn1[1]) - one3*Hkin*epsPn1[1];
  xsi[2] = G*(Tepsilon(2)-epsPn1[2]) - one3*Hkin*epsPn1[2];

  double q = sqrt(two3*xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1] + 2.0*xsi[2]*xsi[2]);
  double F = q - root23*(sigmaY + Hiso*alphan1);

  if (F <= -100*DBL_EPSILON) {
    // Do nothing
  }
  else {
    static Matrix J(4,4);
    static Vector b(4);
    static Vector dx(4);

    double dg = dg_n1;

    J(0,0) = 1.0 + dg*two3*(E+Hkin); J(0,1) = 0.0; J(0,2) = 0.0;
    J(1,0) = 0.0; J(1,1) = 1.0 + dg*(2.0*G+two3*Hkin); J(1,2) = 0.0;
    J(2,0) = 0.0; J(2,1) = 0.0; J(2,2) = 1.0 + dg*(2.0*G+two3*Hkin);
  
    J(0,3) = two3*(E+Hkin)*xsi[0];
    J(1,3) = (2.0*G+two3*Hkin)*xsi[1];
    J(2,3) = (2.0*G+two3*Hkin)*xsi[2];
    
    //J(2,0) = xsi[0]*two3/q; J(2,1) = xsi[1]*2.0/q;
    J(3,0) = (1.0-two3*Hiso*dg)*xsi[0]*two3/q; 
    J(3,1) = (1.0-two3*Hiso*dg)*xsi[1]*2.0/q;
    J(3,2) = (1.0-two3*Hiso*dg)*xsi[2]*2.0/q;

    //J(2,2) = -root23*Hiso;
    J(3,3) = -two3*Hiso*q;

    b(0) = E*depsdh(0) + dEdh*Tepsilon(0) - (E+     Hkin)*depsPdh[0] - (dEdh+     dHkindh)*epsPn1[0];
    b(1) = G*depsdh(1) + dGdh*Tepsilon(1) - (G+one3*Hkin)*depsPdh[1] - (dGdh+one3*dHkindh)*epsPn1[1];
    b(2) = G*depsdh(2) + dGdh*Tepsilon(2) - (G+one3*Hkin)*depsPdh[2] - (dGdh+one3*dHkindh)*epsPn1[2];
    b(3) = root23*(dsigmaYdh + dHisodh*alphan1 + Hiso*dalphadh);

    J.Solve(b, dx);

    dalphadh += dx(3)*root23*q + dg*root23*(xsi[0]*two3*dx(0) + xsi[1]*2.0*dx(1) + xsi[2]*2.0*dx(2))/q;
    depsPdh[0] += dx(3)*two3*xsi[0] + dg*two3*dx(0);
    depsPdh[1] += dx(3)* 2.0*xsi[1] + dg* 2.0*dx(1);
    depsPdh[2] += dx(3)* 2.0*xsi[2] + dg* 2.0*dx(2);

    (*SHVs)(0,gradIndex) = depsPdh[0];
    (*SHVs)(1,gradIndex) = depsPdh[1];
    (*SHVs)(2,gradIndex) = depsPdh[2];
    (*SHVs)(3,gradIndex) = dalphadh;
  }

  return 0;
}


int
J2BeamFiber3d::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(6);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = nu;
  data(3) = sigmaY;
  data(4) = Hiso;
  data(5) = Hkin;
  
  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "J2BeamFiber3d::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

int
J2BeamFiber3d::recvSelf (int commitTag, Channel &theChannel, 
				    FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(6);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
   opserr << "J2BeamFiber3d::recvSelf -- could not recv Vector\n";
   return res;
  }
    
  this->setTag((int)data(0));
  E = data(1);
  nu = data(2);
  sigmaY = data(3);
  Hiso = data(4);
  Hkin = data(5);
  
  return res;
}

void
J2BeamFiber3d::Print (OPS_Stream &s, int flag)
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
J2BeamFiber3d::setParameter(const char **argv, int argc,
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
J2BeamFiber3d::updateParameter(int parameterID, Information &info)
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
J2BeamFiber3d::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}
