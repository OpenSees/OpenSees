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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicBeamFiber2d.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: Elastic isotropic model where stress 
// components 22, 33, 13, and 23 are condensed out.

#include <ElasticIsotropicBeamFiber2d.h>           
#include <Channel.h>
#include <string.h>

Vector ElasticIsotropicBeamFiber2d::sigma(2);
Matrix ElasticIsotropicBeamFiber2d::D(2,2);

ElasticIsotropicBeamFiber2d::ElasticIsotropicBeamFiber2d
(int tag, double E, double nu, double rho):
  ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicBeamFiber2d, E, nu, rho),
  Tepsilon(2)
{

}

ElasticIsotropicBeamFiber2d::ElasticIsotropicBeamFiber2d():
  ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicBeamFiber2d, 0.0, 0.0),
  Tepsilon(2)
{

}

ElasticIsotropicBeamFiber2d::~ElasticIsotropicBeamFiber2d ()
{

}

int
ElasticIsotropicBeamFiber2d::setTrialStrain (const Vector &strain)
{
  Tepsilon = strain;

  return 0;
}

int
ElasticIsotropicBeamFiber2d::setTrialStrain (const Vector &strain, const Vector &rate)
{
  Tepsilon = strain;

  return 0;
}

int
ElasticIsotropicBeamFiber2d::setTrialStrainIncr (const Vector &strain)
{
  return 0;
}

int
ElasticIsotropicBeamFiber2d::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return 0;
}

const Matrix&
ElasticIsotropicBeamFiber2d::getTangent (void)
{
  double mu = 0.50*E/(1.0+v);

  D(0,0) = E;
  D(1,1) = mu;

  return D;
}

const Matrix&
ElasticIsotropicBeamFiber2d::getInitialTangent (void)
{
  double mu = 0.50*E/(1.0+v);

  D(0,0) = E;
  D(1,1) = mu;

  return D;
}

const Vector&
ElasticIsotropicBeamFiber2d::getStress (void)
{
  double mu = 0.50*E/(1.0+v);

  sigma(0) =  E*Tepsilon(0);
  sigma(1) = mu*Tepsilon(1);

  return sigma;
}

const Vector&
ElasticIsotropicBeamFiber2d::getStrain (void)
{
  return Tepsilon;
}

int
ElasticIsotropicBeamFiber2d::commitState (void)
{
  return 0;
}

int
ElasticIsotropicBeamFiber2d::revertToLastCommit (void)
{
  return 0;
}

int
ElasticIsotropicBeamFiber2d::revertToStart (void)
{
  return 0;
}

NDMaterial*
ElasticIsotropicBeamFiber2d::getCopy (void)
{
  ElasticIsotropicBeamFiber2d *theCopy =
    new ElasticIsotropicBeamFiber2d (this->getTag(), E, v, rho);

  return theCopy;
}

const char*
ElasticIsotropicBeamFiber2d::getType (void) const
{
  return "BeamFiber2d";
}

int
ElasticIsotropicBeamFiber2d::getOrder (void) const
{
  return 2;
}

const Vector&
ElasticIsotropicBeamFiber2d::getStressSensitivity(int gradIndex,
						  bool conditional)
{
  sigma(0) = 0.0;
  sigma(1) = 0.0;

  if (parameterID == 1) { // E
    //double mu = 0.5*E/(1.0+v);
    double dmudE = 0.5/(1.0+v);

    sigma(0) = Tepsilon(0);
    sigma(1) = dmudE*Tepsilon(1);
  }

  if (parameterID == 2) { // nu
    //double mu = 0.5*E/(1.0+v);
    double dmudnu = -0.5*E/(1.0 + 2*v + v*v);

    sigma(0) = 0.0;
    sigma(1) = dmudnu*Tepsilon(1);
  }

  return sigma;
}
