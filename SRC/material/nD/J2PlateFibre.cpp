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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2PlateFibre.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: Elastic isotropic model where stress 
// components 22, 33, 13, and 23 are condensed out.

#include <J2PlateFibre.h>           
#include <Channel.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <string.h>
#include <stdlib.h>

Vector J2PlateFibre::sigma(5);
Matrix J2PlateFibre::D(5,5);

void *
OPS_J2PlateFibreMaterial(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 6) {
    opserr << "Want: nDMaterial J2PlateFibre $tag $E $v $sigmaY $Hiso $Hkin <$rho>" << endln;
    return 0;	
  }
  
  int iData[1];
  double dData[6];
  dData[5] = 0.0;
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial J2PlateFibre \n";
    return 0;
  }
  
  if (numArgs > 6) 
    numData = 6;
  else
    numData = 5;
  
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial J2PlateFibre : " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new J2PlateFibre(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]);
  
  return theMaterial;
}

J2PlateFibre::J2PlateFibre
(int tag, double e, double g, double sy, double hi, double hk):
  NDMaterial(tag, ND_TAG_J2PlateFibre),
  E(e), nu(g), sigmaY(sy), Hiso(hi), Hkin(hk),
  parameterID(0), SHVs(0), Tepsilon(5), dg_n1(0.0)
{
  epsPn[0] = 0.0;
  epsPn[1] = 0.0;
  epsPn[2] = 0.0;
  epsPn[3] = 0.0;
  epsPn[4] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
  epsPn1[2] = 0.0;
  epsPn1[3] = 0.0;
  epsPn1[4] = 0.0;

  alphan = 0.0;
  alphan1 = 0.0;
}

J2PlateFibre::J2PlateFibre():
  NDMaterial(0, ND_TAG_J2PlateFibre),
  E(0.0), nu(0.0), sigmaY(0.0), Hiso(0.0), Hkin(0.0), 
  parameterID(0), SHVs(0), Tepsilon(5), dg_n1(0.0)
{
  epsPn[0] = 0.0;
  epsPn[1] = 0.0;
  epsPn[2] = 0.0;
  epsPn[3] = 0.0;
  epsPn[4] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
  epsPn1[2] = 0.0;
  epsPn1[3] = 0.0;
  epsPn1[4] = 0.0;

  alphan = 0.0;
  alphan1 = 0.0;
}

J2PlateFibre::~J2PlateFibre ()
{
  if (SHVs != 0)
    delete SHVs;
}

int
J2PlateFibre::setTrialStrain (const Vector &strain)
{
  Tepsilon = strain;

  return 0;
}

int
J2PlateFibre::setTrialStrain (const Vector &strain, const Vector &rate)
{
  Tepsilon = strain;

  return 0;
}

int
J2PlateFibre::setTrialStrainIncr (const Vector &strain)
{
  return 0;
}

int
J2PlateFibre::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return 0;
}

const Matrix&
J2PlateFibre::getTangent (void)
{
  double twoG = E/(1+nu);
  double G = 0.5*twoG;
  double C00 = E/(1-nu*nu); double C11 = C00;
  double C01 = nu*C00; double C10 = C01;

  double sig[5];
  sig[0] = C00*(Tepsilon(0)-epsPn[0]) + C01*(Tepsilon(1)-epsPn[1]);
  sig[1] = C01*(Tepsilon(0)-epsPn[0]) + C00*(Tepsilon(1)-epsPn[1]);
  sig[2] = G*(Tepsilon(2)-epsPn[2]);
  sig[3] = G*(Tepsilon(3)-epsPn[3]);
  sig[4] = G*(Tepsilon(4)-epsPn[4]);

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double two3Hkin = two3*Hkin;

  double xsi[5];
  xsi[0] = sig[0] - two3Hkin*(2*epsPn[0]+epsPn[1]);
  xsi[1] = sig[1] - two3Hkin*(2*epsPn[1]+epsPn[0]);
  xsi[2] = sig[2] - one3*Hkin*epsPn[2];
  xsi[3] = sig[3] - one3*Hkin*epsPn[3];
  xsi[4] = sig[4] - one3*Hkin*epsPn[4];

  double q = sqrt(two3*(xsi[0]*xsi[0] + xsi[1]*xsi[1] - xsi[0]*xsi[1]) +
		  2.0*(xsi[2]*xsi[2] + xsi[3]*xsi[3] + xsi[4]*xsi[4]));
  double F = q - root23*(sigmaY + Hiso*alphan);

  if (F < -100*DBL_EPSILON) {
    D.Zero();
    D(0,0) = C00; D(0,1) = C01;
    D(1,0) = C01; D(1,1) = C00;
    D(2,2) = G;
    D(3,3) = G;
    D(4,4) = G;

    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
    epsPn1[2] = epsPn[2];
    epsPn1[3] = epsPn[3];
    epsPn1[4] = epsPn[4];
  }
  else {

    // Solve for dg
    double dg = 0.0;

    static Vector R(6);
    static Vector x(6);
    x(0) = xsi[0]; R(0) = 0.0;
    x(1) = xsi[1]; R(1) = 0.0;
    x(2) = xsi[2]; R(2) = 0.0;
    x(3) = xsi[3]; R(3) = 0.0;
    x(4) = xsi[4]; R(4) = 0.0;
    x(5) = dg;     R(5) = F;

    static Matrix J(6,6);
    static Vector dx(6);

    int iter = 0; int maxIter = 25;
    while (iter < maxIter && R.Norm() > 1.0e-14) {
      iter++;

      J(0,0) = 1.0 + dg*(two3*C00-one3*C01+two3Hkin); J(0,1) = dg*(two3*C01-one3*C00);
      J(1,0) = dg*(two3*C10-one3*C11); J(1,1) = 1.0 + dg*(two3*C11-one3*C10+two3Hkin);
      J(2,2) = 1.0 + dg*(twoG+two3Hkin);
      J(3,3) = 1.0 + dg*(twoG+two3Hkin);
      J(4,4) = 1.0 + dg*(twoG+two3Hkin);

      J(0,5) = (two3*C00-one3*C01+two3Hkin)*x(0) + (two3*C01-one3*C00)*x(1);
      J(1,5) = (two3*C10-one3*C11)*x(0) + (two3*C11-one3*C10+two3Hkin)*x(1);
      J(2,5) = (twoG+two3Hkin)*x(2);
      J(3,5) = (twoG+two3Hkin)*x(3);
      J(4,5) = (twoG+two3Hkin)*x(4);

      J(5,0) = (1.0-two3*Hiso*dg)*(two3*x(0)-one3*x(1))/q;
      J(5,1) = (1.0-two3*Hiso*dg)*(two3*x(1)-one3*x(0))/q;
      J(5,2) = (1.0-two3*Hiso*dg)*2.0*x(2)/q;
      J(5,3) = (1.0-two3*Hiso*dg)*2.0*x(3)/q;
      J(5,4) = (1.0-two3*Hiso*dg)*2.0*x(4)/q;

      J(5,5) = -two3*Hiso*q;

      J.Solve(R, dx);
      x.addVector(1.0, dx, -1.0);

      dg = x(5);
      dg_n1 = dg;

      q = sqrt(two3*(x(0)*x(0) + x(1)*x(1) - x(0)*x(1)) + 2.0*(x(2)*x(2) + x(3)*x(3) + x(4)*x(4)));

      R(0) = x(0) - xsi[0] + dg*((two3*C00-one3*C01+two3Hkin)*x(0) + (two3*C01-one3*C00)*x(1));
      R(1) = x(1) - xsi[1] + dg*((two3*C10-one3*C11)*x(0) + (two3*C11-one3*C10+two3Hkin)*x(1));
      R(2) = x(2) - xsi[2] + dg*(twoG+two3Hkin)*x(2);
      R(3) = x(3) - xsi[3] + dg*(twoG+two3Hkin)*x(3);
      R(4) = x(4) - xsi[4] + dg*(twoG+two3Hkin)*x(4);
      R(5) = q - root23*(sigmaY + Hiso*(alphan+dg*root23*q));
    }

    if (iter == maxIter) {
      //opserr << "J2PlateFibre::getTangent -- maxIter reached " << R.Norm() << endln;
    }

    alphan1 = alphan + dg*root23*q;

    epsPn1[0] = epsPn[0] + dg*(two3*x(0)-one3*x(1));
    epsPn1[1] = epsPn[1] + dg*(two3*x(1)-one3*x(0));
    epsPn1[2] = epsPn[2] + dg*2.0*x(2);
    epsPn1[3] = epsPn[3] + dg*2.0*x(3);
    epsPn1[4] = epsPn[4] + dg*2.0*x(4);

    double beta = 1.0+dg*two3Hkin;

    //J.Zero();
    J(0,0) = 1.0 + dg*(two3*C00-one3*C01)/beta; J(0,1) = dg*(two3*C01-one3*C00)/beta;
    J(1,0) = dg*(two3*C10-one3*C11)/beta; J(1,1) = 1.0 + dg*(two3*C11-one3*C10)/beta;
    J(2,2) = 1.0 + dg*twoG/beta;
    J(3,3) = 1.0 + dg*twoG/beta;
    J(4,4) = 1.0 + dg*twoG/beta;

    //J(0,5) = ((two3*C00-one3*C01) - dg*(two3*C00-one3*C01)/beta*two3Hkin)*x(0) + 
    //  ((two3*C01-one3*C00) - dg*(two3*C01-one3*C00)/beta*two3Hkin)*x(1);
    //J(1,5) = ((two3*C10-one3*C11) - dg*(two3*C10-one3*C11)/beta*two3Hkin)*x(0) + 
    //  ((two3*C11-one3*C10) - dg*(two3*C11-one3*C10)/beta*two3Hkin)*x(1);
    J(0,5) = (two3*C00-one3*C01)*(1.0-dg/beta*two3Hkin)*x(0) + (two3*C01-one3*C00)*(1.0-dg/beta*two3Hkin)*x(1);
    J(1,5) = (two3*C10-one3*C11)*(1.0-dg/beta*two3Hkin)*x(0) + (two3*C11-one3*C10)*(1.0-dg/beta*two3Hkin)*x(1);
    J(2,5) = (twoG - dg*twoG/beta*two3Hkin)*x(2);
    J(3,5) = (twoG - dg*twoG/beta*two3Hkin)*x(3);
    J(4,5) = (twoG - dg*twoG/beta*two3Hkin)*x(4);

    J(5,0) = (1.0-two3*Hiso*dg)*(two3*x(0)-one3*x(1))/q/beta;
    J(5,1) = (1.0-two3*Hiso*dg)*(two3*x(1)-one3*x(0))/q/beta;
    J(5,2) = (1.0-two3*Hiso*dg)*x(2)/q*2.0/beta;
    J(5,3) = (1.0-two3*Hiso*dg)*x(3)/q*2.0/beta;
    J(5,4) = (1.0-two3*Hiso*dg)*x(4)/q*2.0/beta;

    J(5,5) = -q*two3Hkin/beta - two3*Hiso*q;

    static Matrix invJ(6,6);
    J.Invert(invJ);

    D(0,0) = invJ(0,0)*C00 + invJ(0,1)*C10;
    D(1,0) = invJ(1,0)*C00 + invJ(1,1)*C10;
    D(2,0) = invJ(2,0)*C00 + invJ(2,1)*C10;
    D(3,0) = invJ(3,0)*C00 + invJ(3,1)*C10;
    D(4,0) = invJ(4,0)*C00 + invJ(4,1)*C10;

    D(0,1) = invJ(0,0)*C01 + invJ(0,1)*C11;
    D(1,1) = invJ(1,0)*C01 + invJ(1,1)*C11;
    D(2,1) = invJ(2,0)*C01 + invJ(2,1)*C11;
    D(3,1) = invJ(3,0)*C01 + invJ(3,1)*C11;
    D(4,1) = invJ(4,0)*C01 + invJ(4,1)*C11;

    D(0,2) = invJ(0,2)*G;
    D(1,2) = invJ(1,2)*G;
    D(2,2) = invJ(2,2)*G;
    D(3,2) = invJ(3,2)*G;
    D(4,2) = invJ(4,2)*G;

    D(0,3) = invJ(0,3)*G;
    D(1,3) = invJ(1,3)*G;
    D(2,3) = invJ(2,3)*G;
    D(3,3) = invJ(3,3)*G;
    D(4,3) = invJ(4,3)*G;

    D(0,4) = invJ(0,4)*G;
    D(1,4) = invJ(1,4)*G;
    D(2,4) = invJ(2,4)*G;
    D(3,4) = invJ(3,4)*G;
    D(4,4) = invJ(4,4)*G;
  }

  return D;
}

const Matrix&
J2PlateFibre::getInitialTangent (void)
{
  double G = 0.5*E/(1+nu);
  double C00 = E/(1-nu*nu);
  double C01 = nu*C00;

  D.Zero();
  D(0,0) = C00; D(0,1) = C01;
  D(1,0) = C01; D(1,1) = C00;
  D(2,2) = G;
  D(3,3) = G;
  D(4,4) = G;

  return D;
}

const Vector&
J2PlateFibre::getStress (void)
{
  double twoG = E/(1+nu);
  double G = 0.5*twoG;
  double C00 = E/(1-nu*nu); double C11 = C00;
  double C01 = nu*C00; double C10 = C01;

  sigma(0) = C00*(Tepsilon(0)-epsPn[0]) + C01*(Tepsilon(1)-epsPn[1]);
  sigma(1) = C01*(Tepsilon(0)-epsPn[0]) + C00*(Tepsilon(1)-epsPn[1]);
  sigma(2) = G*(Tepsilon(2)-epsPn[2]);
  sigma(3) = G*(Tepsilon(3)-epsPn[3]);
  sigma(4) = G*(Tepsilon(4)-epsPn[4]);

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double two3Hkin = two3*Hkin;

  double xsi[5];
  xsi[0] = sigma(0) - two3Hkin*(2*epsPn[0]+epsPn[1]);
  xsi[1] = sigma(1) - two3Hkin*(2*epsPn[1]+epsPn[0]);
  xsi[2] = sigma(2) - one3*Hkin*epsPn[2];
  xsi[3] = sigma(3) - one3*Hkin*epsPn[3];
  xsi[4] = sigma(4) - one3*Hkin*epsPn[4];

  double q = sqrt(two3*(xsi[0]*xsi[0] + xsi[1]*xsi[1] - xsi[0]*xsi[1]) +
		  2.0*(xsi[2]*xsi[2] + xsi[3]*xsi[3] + xsi[4]*xsi[4]));
  double F = q - root23*(sigmaY + Hiso*alphan);

  if (F < -100*DBL_EPSILON) {
    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
    epsPn1[2] = epsPn[2];
    epsPn1[3] = epsPn[3];
    epsPn1[4] = epsPn[4];
  }
  else {
    // Solve for dg
    double dg = 0.0;

    static Vector R(6);
    static Vector x(6);
    x(0) = xsi[0]; R(0) = 0.0;
    x(1) = xsi[1]; R(1) = 0.0;
    x(2) = xsi[2]; R(2) = 0.0;
    x(3) = xsi[3]; R(3) = 0.0;
    x(4) = xsi[4]; R(4) = 0.0;
    x(5) = dg;     R(5) = F;

    static Matrix J(6,6);
    static Vector dx(6);

    int iter = 0; int maxIter = 25;
    while (iter < maxIter && R.Norm() > 1.0e-14) {
      iter++;

      J(0,0) = 1.0 + dg*(two3*C00-one3*C01+two3Hkin); J(0,1) = dg*(two3*C01-one3*C00);
      J(1,0) = dg*(two3*C10-one3*C11); J(1,1) = 1.0 + dg*(two3*C11-one3*C10+two3Hkin);
      J(2,2) = 1.0 + dg*(twoG+two3Hkin);
      J(3,3) = 1.0 + dg*(twoG+two3Hkin);
      J(4,4) = 1.0 + dg*(twoG+two3Hkin);

      J(0,5) = (two3*C00-one3*C01+two3Hkin)*x(0) + (two3*C01-one3*C00)*x(1);
      J(1,5) = (two3*C10-one3*C11)*x(0) + (two3*C11-one3*C10+two3Hkin)*x(1);
      J(2,5) = (twoG+two3Hkin)*x(2);
      J(3,5) = (twoG+two3Hkin)*x(3);
      J(4,5) = (twoG+two3Hkin)*x(4);

      J(5,0) = (1.0-two3*Hiso*dg)*(two3*x(0)-one3*x(1))/q;
      J(5,1) = (1.0-two3*Hiso*dg)*(two3*x(1)-one3*x(0))/q;
      J(5,2) = (1.0-two3*Hiso*dg)*2.0*x(2)/q;
      J(5,3) = (1.0-two3*Hiso*dg)*2.0*x(3)/q;
      J(5,4) = (1.0-two3*Hiso*dg)*2.0*x(4)/q;

      J(5,5) = -two3*Hiso*q;

      J.Solve(R, dx);
      x.addVector(1.0, dx, -1.0);

      dg = x(5);
      dg_n1 = dg;

      q = sqrt(two3*(x(0)*x(0) + x(1)*x(1) - x(0)*x(1)) + 2.0*(x(2)*x(2) + x(3)*x(3) + x(4)*x(4)));

      R(0) = x(0) - xsi[0] + dg*((two3*C00-one3*C01+two3Hkin)*x(0) + (two3*C01-one3*C00)*x(1));
      R(1) = x(1) - xsi[1] + dg*((two3*C10-one3*C11)*x(0) + (two3*C11-one3*C10+two3Hkin)*x(1));
      R(2) = x(2) - xsi[2] + dg*(twoG+two3Hkin)*x(2);
      R(3) = x(3) - xsi[3] + dg*(twoG+two3Hkin)*x(3);
      R(4) = x(4) - xsi[4] + dg*(twoG+two3Hkin)*x(4);
      R(5) = q - root23*(sigmaY + Hiso*(alphan+dg*root23*q));
    }
    if (iter == maxIter) {
      //opserr << "J2PlateFibre::getStress -- maxIter reached " << R.Norm() << endln;
    }

    alphan1 = alphan + dg*root23*q;

    epsPn1[0] = epsPn[0] + dg*(two3*x(0)-one3*x(1));
    epsPn1[1] = epsPn[1] + dg*(two3*x(1)-one3*x(0));
    epsPn1[2] = epsPn[2] + dg*2.0*x(2);
    epsPn1[3] = epsPn[3] + dg*2.0*x(3);
    epsPn1[4] = epsPn[4] + dg*2.0*x(4);

    sigma(0) = x(0) + two3Hkin*(2.0*epsPn1[0]+epsPn1[1]);
    sigma(1) = x(1) + two3Hkin*(2.0*epsPn1[1]+epsPn1[0]);
    sigma(2) = x(2) + one3*Hkin*epsPn1[2];
    sigma(3) = x(3) + one3*Hkin*epsPn1[3];
    sigma(4) = x(4) + one3*Hkin*epsPn1[4];
  }

  return sigma;
}

const Vector&
J2PlateFibre::getStrain (void)
{
  return Tepsilon;
}

int
J2PlateFibre::commitState (void)
{
  epsPn[0] = epsPn1[0];
  epsPn[1] = epsPn1[1];
  epsPn[2] = epsPn1[2];
  epsPn[3] = epsPn1[3];
  epsPn[4] = epsPn1[4];

  alphan = alphan1;

  return 0;
}

int
J2PlateFibre::revertToLastCommit (void)
{
  epsPn1[0] = epsPn[0];
  epsPn1[1] = epsPn[1];
  epsPn1[2] = epsPn[2];
  epsPn1[3] = epsPn[3];
  epsPn1[4] = epsPn[4];

  alphan1 = alphan;

  return 0;
}

int
J2PlateFibre::revertToStart (void)
{
  Tepsilon.Zero();

  epsPn[0] = 0.0;
  epsPn[1] = 0.0;
  epsPn[2] = 0.0;
  epsPn[3] = 0.0;
  epsPn[4] = 0.0;

  epsPn1[0] = 0.0;
  epsPn1[1] = 0.0;
  epsPn1[2] = 0.0;
  epsPn1[3] = 0.0;
  epsPn1[4] = 0.0;

  alphan = 0.0;
  alphan1 = 0.0;

  dg_n1 = 0.0;

  if (SHVs != 0)
    SHVs->Zero();

  return 0;
}

NDMaterial*
J2PlateFibre::getCopy (void)
{
  J2PlateFibre *theCopy =
    new J2PlateFibre (this->getTag(), E, nu, sigmaY, Hiso, Hkin);

  return theCopy;
}

NDMaterial*
J2PlateFibre::getCopy (const char *type)
{
  if (strcmp(type,this->getType()) == 0)
    return this->getCopy();

  return 0;
}

const char*
J2PlateFibre::getType (void) const
{
  return "PlateFiber";
}

int
J2PlateFibre::getOrder (void) const
{
  return 5;
}

const Vector&
J2PlateFibre::getStressSensitivity(int gradIndex, bool conditional)
{
  sigma(0) = 0.0;
  sigma(1) = 0.0;
  sigma(2) = 0.0;
  sigma(3) = 0.0;
  sigma(4) = 0.0;

  double twoG = E/(1+nu);
  double G = 0.5*twoG;
  double C00 = E/(1-nu*nu); double C11 = C00;
  double C01 = nu*C00; double C10 = C01;

  double dEdh = 0.0;
  double dsigmaYdh = 0.0;
  double dHkindh = 0.0;
  double dHisodh = 0.0;
  double dGdh = 0.0;
  double dC00dh = 0.0;
  double dC01dh = 0.0;

  if (parameterID == 1) { // E
    dEdh = 1.0;
    dGdh = 0.5/(1+nu);
    //double C00 = E/(1-nu*nu);
    //double C01 = nu*C00;
    dC00dh = 1.0/(1-nu*nu);
    dC01dh = nu*dC00dh;
  }
  if (parameterID == 2) { // G (nu)
    dGdh = 1.0;
    dGdh = -0.5*E/(1.0 + 2*nu + nu*nu);
    //double C00 = E/(1-nu*nu);
    //double C01 = nu*C00;
    dC00dh = -E/((1-nu*nu)*(1-nu*nu))*(-2*nu);
    dC01dh = nu*dC00dh + C00;
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

  double depsPdh[5]; 
  depsPdh[0] = 0.0; 
  depsPdh[1] = 0.0;
  depsPdh[2] = 0.0;
  depsPdh[3] = 0.0;
  depsPdh[4] = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    depsPdh[0] = (*SHVs)(0,gradIndex);
    depsPdh[1] = (*SHVs)(1,gradIndex);
    depsPdh[2] = (*SHVs)(2,gradIndex);
    depsPdh[3] = (*SHVs)(3,gradIndex);
    depsPdh[4] = (*SHVs)(4,gradIndex);
    dalphadh = (*SHVs)(5,gradIndex);
  }

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double two3Hkin = two3*Hkin;

  double xsi[5];
  xsi[0] = C00*(Tepsilon(0)-epsPn[0]) + C01*(Tepsilon(1)-epsPn[1]);
  xsi[1] = C01*(Tepsilon(0)-epsPn[0]) + C00*(Tepsilon(1)-epsPn[1]);
  xsi[2] = G*(Tepsilon(2)-epsPn[2]);
  xsi[3] = G*(Tepsilon(3)-epsPn[3]);
  xsi[4] = G*(Tepsilon(4)-epsPn[4]);

  xsi[0] -= two3Hkin*(2*epsPn[0]+epsPn[1]); 
  xsi[1] -= two3Hkin*(2*epsPn[1]+epsPn[0]); 
  xsi[2] -= one3*Hkin*epsPn[2]; 
  xsi[3] -= one3*Hkin*epsPn[3]; 
  xsi[4] -= one3*Hkin*epsPn[4];

  double dbetadh[5]; 
  dbetadh[0] = two3*dHkindh*(2*epsPn[0]+epsPn[1]) + two3Hkin*(2*depsPdh[0]+depsPdh[1]);
  dbetadh[1] = two3*dHkindh*(2*epsPn[1]+epsPn[0]) + two3Hkin*(2*depsPdh[1]+depsPdh[0]);
  dbetadh[2] = one3*(dHkindh*epsPn[2] + Hkin*depsPdh[2]);
  dbetadh[3] = one3*(dHkindh*epsPn[3] + Hkin*depsPdh[3]);
  dbetadh[4] = one3*(dHkindh*epsPn[4] + Hkin*depsPdh[4]);

  double dxsidh[5];
  dxsidh[0] = -C00*depsPdh[0] - C01*depsPdh[1] + dC00dh*(Tepsilon(0)-epsPn[0]) + dC01dh*(Tepsilon(1)-epsPn[1]) - dbetadh[0];
  dxsidh[1] = -C01*depsPdh[0] - C00*depsPdh[1] + dC01dh*(Tepsilon(0)-epsPn[0]) + dC00dh*(Tepsilon(1)-epsPn[1]) - dbetadh[1];
  dxsidh[2] = -G*depsPdh[2] + dGdh*(Tepsilon(2)-epsPn[2]) - dbetadh[2];
  dxsidh[3] = -G*depsPdh[3] + dGdh*(Tepsilon(3)-epsPn[3]) - dbetadh[3];
  dxsidh[4] = -G*depsPdh[4] + dGdh*(Tepsilon(4)-epsPn[4]) - dbetadh[4];

  double q = two3*(xsi[0]*xsi[0] + xsi[1]*xsi[1] - xsi[0]*xsi[1]) + 
    2.0*(xsi[2]*xsi[2] + xsi[3]*xsi[3] + xsi[4]*xsi[4]);
  double F = q - root23*(sigmaY + Hiso*alphan1);

  if (F <= -100*DBL_EPSILON) {
    sigma(0) = dC00dh*(Tepsilon(0)-epsPn[0]) + dC01dh*(Tepsilon(1)-epsPn[1]) - C00*depsPdh[0] - C01*depsPdh[1];
    sigma(1) = dC01dh*(Tepsilon(0)-epsPn[0]) + dC00dh*(Tepsilon(1)-epsPn[1]) - C01*depsPdh[0] - C00*depsPdh[1];
    sigma(2) = dGdh*(Tepsilon(2)-epsPn1[2]) - G*depsPdh[2];
    sigma(3) = dGdh*(Tepsilon(3)-epsPn1[3]) - G*depsPdh[3];
    sigma(4) = dGdh*(Tepsilon(4)-epsPn1[4]) - G*depsPdh[4];
  }
  else {
    static Matrix J(6,6);
    static Vector b(6);
    static Vector dx(6);

    double dg = dg_n1;

    J(0,0) = 1.0 + dg*(two3*C00-one3*C01+two3Hkin); J(0,1) = dg*(two3*C01-one3*C00);
    J(1,0) = dg*(two3*C10-one3*C11); J(1,1) = 1.0 + dg*(two3*C11-one3*C10+two3Hkin);
    J(2,2) = 1.0 + dg*(twoG+two3Hkin);
    J(3,3) = 1.0 + dg*(twoG+two3Hkin);
    J(4,4) = 1.0 + dg*(twoG+two3Hkin);

    J(0,5) = (two3*C00-one3*C01+two3Hkin)*xsi[0] + (two3*C01-one3*C00)*xsi[1];
    J(1,5) = (two3*C10-one3*C11)*xsi[0] + (two3*C11-one3*C10+two3Hkin)*xsi[1];
    J(2,5) = (twoG+two3Hkin)*xsi[2];
    J(3,5) = (twoG+two3Hkin)*xsi[3];
    J(4,5) = (twoG+two3Hkin)*xsi[4];

    J(5,0) = (1.0-two3*Hiso*dg)*(two3*xsi[0]-one3*xsi[1])/q;
    J(5,1) = (1.0-two3*Hiso*dg)*(two3*xsi[1]-one3*xsi[0])/q;
    J(5,2) = (1.0-two3*Hiso*dg)*2.0*xsi[2]/q;
    J(5,3) = (1.0-two3*Hiso*dg)*2.0*xsi[3]/q;
    J(5,4) = (1.0-two3*Hiso*dg)*2.0*xsi[4]/q;

    J(5,5) = -two3*Hiso*q;

    J.Solve(b, dx);
  }

  return sigma;
}

int
J2PlateFibre::commitSensitivity(const Vector &depsdh, int gradIndex, int numGrads)
{
  if (SHVs == 0) {
    SHVs = new Matrix(6,numGrads);
  }

  if (gradIndex >= SHVs->noCols()) {
    //opserr << gradIndex << ' ' << SHVs->noCols() << endln;
    return 0;
  }

  double twoG = E/(1+nu);
  double G = 0.5*twoG;
  double C00 = E/(1-nu*nu); double C11 = C00;
  double C01 = nu*C00; double C10 = C01;

  double dEdh = 0.0;
  double dsigmaYdh = 0.0;
  double dHkindh = 0.0;
  double dHisodh = 0.0;
  double dGdh = 0.0;
  double dC00dh = 0.0;
  double dC01dh = 0.0;

  if (parameterID == 1) { // E
    dEdh = 1.0;
    dGdh = 0.5/(1+nu);
    //double C00 = E/(1-nu*nu);
    //double C01 = nu*C00;
    dC00dh = 1.0/(1-nu*nu);
    dC01dh = nu*dC00dh;
  }
  if (parameterID == 2) { // G (nu)
    dGdh = 1.0;
    dGdh = -0.5*E/(1.0 + 2*nu + nu*nu);
    //double C00 = E/(1-nu*nu);
    //double C01 = nu*C00;
    dC00dh = -E/((1-nu*nu)*(1-nu*nu))*(-2*nu);
    dC01dh = nu*dC00dh + C00;
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

  double depsPdh[5];
  depsPdh[0] = 0.0;
  depsPdh[1] = 0.0;
  depsPdh[2] = 0.0;
  depsPdh[3] = 0.0;
  depsPdh[4] = 0.0;
  double dalphadh = 0.0;
  if (SHVs != 0) {
    depsPdh[0] = (*SHVs)(0,gradIndex);
    depsPdh[1] = (*SHVs)(1,gradIndex);
    depsPdh[2] = (*SHVs)(2,gradIndex);
    depsPdh[3] = (*SHVs)(3,gradIndex);
    depsPdh[4] = (*SHVs)(4,gradIndex);
    dalphadh = (*SHVs)(5,gradIndex);
  }

  static const double one3 = 1.0/3;
  static const double two3 = 2.0*one3;
  static const double root23 = sqrt(two3);

  double two3Hkin = two3*Hkin;

  double xsi[5];
  xsi[0] = C00*(Tepsilon(0)-epsPn[0]) + C01*(Tepsilon(1)-epsPn[1]);
  xsi[1] = C01*(Tepsilon(0)-epsPn[0]) + C00*(Tepsilon(1)-epsPn[1]);
  xsi[2] = G*(Tepsilon(2)-epsPn[2]);
  xsi[3] = G*(Tepsilon(3)-epsPn[3]);
  xsi[4] = G*(Tepsilon(4)-epsPn[4]);

  xsi[0] -= two3Hkin*(2*epsPn[0]+epsPn[1]);
  xsi[1] -= two3Hkin*(2*epsPn[1]+epsPn[0]);
  xsi[2] -= one3*Hkin*epsPn[2];
  xsi[3] -= one3*Hkin*epsPn[3];
  xsi[4] -= one3*Hkin*epsPn[4];

  double q = two3*(xsi[0]*xsi[0] + xsi[1]*xsi[1] - xsi[0]*xsi[1]) +
    2.0*(xsi[2]*xsi[2] + xsi[3]*xsi[3] + xsi[4]*xsi[4]); 
  double F = q - root23*(sigmaY + Hiso*alphan1);

  if (F <= -100*DBL_EPSILON) {
    // Do nothing
  }
  else {
    static Matrix J(6,6);
    static Vector b(6);
    static Vector dx(6);
    
    double dg = dg_n1;

    J(0,0) = 1.0 + dg*(two3*C00-one3*C01+two3Hkin); J(0,1) = dg*(two3*C01-one3*C00);
    J(1,0) = dg*(two3*C10-one3*C11); J(1,1) = 1.0 + dg*(two3*C11-one3*C10+two3Hkin);
    J(2,2) = 1.0 + dg*(twoG+two3Hkin);
    J(3,3) = 1.0 + dg*(twoG+two3Hkin);
    J(4,4) = 1.0 + dg*(twoG+two3Hkin);

    J(0,5) = (two3*C00-one3*C01+two3Hkin)*xsi[0] + (two3*C01-one3*C00)*xsi[1];
    J(1,5) = (two3*C10-one3*C11)*xsi[0] + (two3*C11-one3*C10+two3Hkin)*xsi[1];
    J(2,5) = (twoG+two3Hkin)*xsi[2];
    J(3,5) = (twoG+two3Hkin)*xsi[3];
    J(4,5) = (twoG+two3Hkin)*xsi[4];

    J(5,0) = (1.0-two3*Hiso*dg)*(two3*xsi[0]-one3*xsi[1])/q;
    J(5,1) = (1.0-two3*Hiso*dg)*(two3*xsi[1]-one3*xsi[0])/q;
    J(5,2) = (1.0-two3*Hiso*dg)*2.0*xsi[2]/q;
    J(5,3) = (1.0-two3*Hiso*dg)*2.0*xsi[3]/q;
    J(5,4) = (1.0-two3*Hiso*dg)*2.0*xsi[4]/q;

    J(5,5) = -two3*Hiso*q;
    
    J.Solve(b, dx);
  }

  return 0;
}


int
J2PlateFibre::sendSelf (int commitTag, Channel &theChannel)
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
    opserr << "J2PlateFibre::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

int
J2PlateFibre::recvSelf (int commitTag, Channel &theChannel, 
				    FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(6);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
   opserr << "J2PlateFibre::recvSelf -- could not recv Vector\n";
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
J2PlateFibre::Print (OPS_Stream &s, int flag)
{
  s << "J2 Plate Fibre Material Model" << endln;
  s << "\tE:  " << E << endln;
  s << "\tnu:  " << nu << endln;
  s << "\tsigmaY:  " << sigmaY << endln;
  s << "\tHiso:  " << Hiso << endln;
  s << "\tHkin:  " << Hkin << endln;
  
  return;
}

int
J2PlateFibre::setParameter(const char **argv, int argc,
			    Parameter &param)
{
  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);
  
  else if (strcmp(argv[0],"nu") == 0)
    return param.addObject(2, this);  
  
  else if (strcmp(argv[0],"sigmaY") == 0 || strcmp(argv[0],"fy") == 0)
    return param.addObject(5, this);

  else if (strcmp(argv[0],"Hkin") == 0)
    return param.addObject(6, this);

  else if (strcmp(argv[0],"Hiso") == 0)
    return param.addObject(7, this);

  return -1;
}

int 
J2PlateFibre::updateParameter(int parameterID, Information &info)
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
J2PlateFibre::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}
