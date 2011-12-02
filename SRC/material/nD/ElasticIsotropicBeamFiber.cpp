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
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicBeamFiber, E, nu, rho),
 Tepsilon(3)
{

}

ElasticIsotropicBeamFiber::ElasticIsotropicBeamFiber():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicBeamFiber, 0.0, 0.0),
 Tepsilon(3)
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
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

	D(0,0) = mu*(mu2+3.0*lam)/(mu+lam);
	D(1,1) = mu;
	D(2,2) = mu;

	return D;
}

const Matrix&
ElasticIsotropicBeamFiber::getInitialTangent (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

	D(0,0) = mu*(mu2+3.0*lam)/(mu+lam);
	D(1,1) = mu;
	D(2,2) = mu;

	return D;
}

const Vector&
ElasticIsotropicBeamFiber::getStress (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

	sigma(0) = mu*(mu2+3.0*lam)/(mu+lam)*Tepsilon(0);
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

