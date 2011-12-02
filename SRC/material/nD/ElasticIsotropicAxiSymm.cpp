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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicAxiSymm.cpp,v $
                                                                 
#include <ElasticIsotropicAxiSymm.h>                                                                        
#include <Channel.h>

Vector ElasticIsotropicAxiSymm::sigma(4);
Matrix ElasticIsotropicAxiSymm::D(4,4);

ElasticIsotropicAxiSymm::ElasticIsotropicAxiSymm
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicAxiSymm, E, nu, rho),
 epsilon(4)
{

}

ElasticIsotropicAxiSymm::ElasticIsotropicAxiSymm():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicAxiSymm, 0.0, 0.0),
 epsilon(4)
{

}

ElasticIsotropicAxiSymm::~ElasticIsotropicAxiSymm ()
{

}

int
ElasticIsotropicAxiSymm::setTrialStrain (const Vector &strain)
{
	epsilon = strain;

	return 0;
}

int
ElasticIsotropicAxiSymm::setTrialStrain (const Vector &strain, const Vector &rate)
{
	epsilon = strain;

	return 0;
}

int
ElasticIsotropicAxiSymm::setTrialStrainIncr (const Vector &strain)
{
	epsilon+=strain;

	return 0;
}

int
ElasticIsotropicAxiSymm::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
	epsilon+=strain;

	return 0;
}

const Matrix&
ElasticIsotropicAxiSymm::getTangent (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

	D(0,0) = D(1,1) = D(2,2) = mu2+lam;
	D(0,1) = D(1,0) = lam;
	D(0,2) = D(2,0) = lam;
	D(1,2) = D(2,1) = lam;
	D(3,3) = mu;

	return D;
}

const Matrix&
ElasticIsotropicAxiSymm::getInitialTangent (void)
{
	double mu2 = E/(1.0+v);
	double lam = v*mu2/(1.0-2.0*v);
	double mu = 0.50*mu2;

	D(0,0) = D(1,1) = D(2,2) = mu2+lam;
	D(0,1) = D(1,0) = lam;
	D(0,2) = D(2,0) = lam;
	D(1,2) = D(2,1) = lam;
	D(3,3) = mu;

	return D;
}

const Vector&
ElasticIsotropicAxiSymm::getStress (void)
{
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;

  double eps0 = epsilon(0);
  double eps1 = epsilon(1);
  double eps2 = epsilon(2);

  mu2 += lam;

  //sigma = D*epsilon;
  sigma(0) = mu2*eps0 + lam*(eps1+eps2);
  sigma(1) = mu2*eps1 + lam*(eps0+eps2);
  sigma(2) = mu2*eps2 + lam*(eps0+eps1);
  sigma(3) = mu*epsilon(3);
	
  return sigma;
}

const Vector&
ElasticIsotropicAxiSymm::getStrain (void)
{
	return epsilon;
}

int
ElasticIsotropicAxiSymm::commitState (void)
{
  return 0;
}

int
ElasticIsotropicAxiSymm::revertToLastCommit (void)
{
  return 0;
}

int
ElasticIsotropicAxiSymm::revertToStart (void)
{
  epsilon.Zero();
  return 0;
}

NDMaterial*
ElasticIsotropicAxiSymm::getCopy (void)
{
	ElasticIsotropicAxiSymm *theCopy =
		new ElasticIsotropicAxiSymm (this->getTag(), E, v, rho);

	theCopy->epsilon = epsilon;

	return theCopy;
}

const char*
ElasticIsotropicAxiSymm::getType (void) const
{
	return "AxiSymmetric";
}

int
ElasticIsotropicAxiSymm::getOrder (void) const
{
	return 4;
}

