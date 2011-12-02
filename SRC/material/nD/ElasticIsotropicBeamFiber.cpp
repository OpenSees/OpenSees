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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicBeamFiber.cpp,v $

// Written: MHS
// Created: Aug 2001
//
// Description: Elastic isotropic model where stress components 22, 33, and 23
// are condensed out.

#include <ElasticIsotropicBeamFiber.h>           
#include <Channel.h>
#include <Tensor.h>

Vector ElasticIsotropicBeamFiber::sigma(3);
Matrix ElasticIsotropicBeamFiber::D(3,3);

ElasticIsotropicBeamFiber::ElasticIsotropicBeamFiber
(int tag, double E, double nu, double rho):
  ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicBeamFiber, E, nu, rho),
  Tepsilon(3), parameterID(0)
{

}

ElasticIsotropicBeamFiber::ElasticIsotropicBeamFiber():
  ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicBeamFiber, 0.0, 0.0),
  Tepsilon(3), parameterID(0)
{

}

ElasticIsotropicBeamFiber::~ElasticIsotropicBeamFiber ()
{

}

int
ElasticIsotropicBeamFiber::setTrialStrain (const Vector &strain)
{
  Tepsilon = strain;

  return 0;
}

int
ElasticIsotropicBeamFiber::setTrialStrain (const Vector &strain, const Vector &rate)
{
  Tepsilon = strain;

  return 0;
}

int
ElasticIsotropicBeamFiber::setTrialStrainIncr (const Vector &strain)
{
  return 0;
}

int
ElasticIsotropicBeamFiber::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return 0;
}

const Matrix&
ElasticIsotropicBeamFiber::getTangent (void)
{
  double mu = 0.5*E/(1.0+v);

  D(0,0) = E;
  D(1,1) = mu;
  D(2,2) = mu;
  
  return D;
}

const Matrix&
ElasticIsotropicBeamFiber::getInitialTangent (void)
{
  double mu = 0.5*E/(1.0+v);

  D(0,0) = E;
  D(1,1) = mu;
  D(2,2) = mu;
  
  return D;
}

const Vector&
ElasticIsotropicBeamFiber::getStress (void)
{
  double mu = 0.5*E/(1.0+v);

  sigma(0) =  E*Tepsilon(0);
  sigma(1) = mu*Tepsilon(1);
  sigma(2) = mu*Tepsilon(2);
  
  return sigma;
}

const Vector&
ElasticIsotropicBeamFiber::getStrain (void)
{
  return Tepsilon;
}

int
ElasticIsotropicBeamFiber::commitState (void)
{
  return 0;
}

int
ElasticIsotropicBeamFiber::revertToLastCommit (void)
{
  return 0;
}

int
ElasticIsotropicBeamFiber::revertToStart (void)
{
  return 0;
}

NDMaterial*
ElasticIsotropicBeamFiber::getCopy (void)
{
  ElasticIsotropicBeamFiber *theCopy =
    new ElasticIsotropicBeamFiber (this->getTag(), E, v, rho);

  return theCopy;
}

const char*
ElasticIsotropicBeamFiber::getType (void) const
{
  return "BeamFiber";
}

int
ElasticIsotropicBeamFiber::getOrder (void) const
{
  return 3;
}

int
ElasticIsotropicBeamFiber::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);

  return -1;
}

int
ElasticIsotropicBeamFiber::updateParameter(int paramID, Information &info)
{
  if (paramID == 1)
    E = info.theDouble;

  return 0;
}

int
ElasticIsotropicBeamFiber::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
ElasticIsotropicBeamFiber::getStressSensitivity(int gradIndex,
						bool conditional)
{
  sigma.Zero();

  if (parameterID == 1) { // E
    //double mu = 0.5*E/(1.0+v);
    double dmudE = 0.5/(1.0+v);

    sigma(0) = Tepsilon(0);
    sigma(1) = dmudE*Tepsilon(1);
    sigma(2) = dmudE*Tepsilon(2);
  }

  return sigma;
}
