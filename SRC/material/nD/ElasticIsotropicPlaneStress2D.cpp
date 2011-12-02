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
// $Date: 2002-12-05 22:49:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStress2D.cpp,v $
                                                                        
                                                                        

#include <ElasticIsotropicPlaneStress2D.h>           
#include <Channel.h>
#include <Tensor.h>

Vector ElasticIsotropicPlaneStress2D::sigma(3);
Matrix ElasticIsotropicPlaneStress2D::D(3,3);

ElasticIsotropicPlaneStress2D::ElasticIsotropicPlaneStress2D
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicPlaneStress2d, E, nu, rho),
 epsilon(3)
{

}

ElasticIsotropicPlaneStress2D::ElasticIsotropicPlaneStress2D():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicPlaneStress2d, 0.0, 0.0),
 epsilon(3)
{

}

ElasticIsotropicPlaneStress2D::~ElasticIsotropicPlaneStress2D ()
{

}

int
ElasticIsotropicPlaneStress2D::setTrialStrain (const Vector &strain)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropicPlaneStress2D::setTrialStrain (const Vector &strain, const Vector &rate)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropicPlaneStress2D::setTrialStrainIncr (const Vector &strain)
{
  epsilon += strain;
  return 0;
}

int
ElasticIsotropicPlaneStress2D::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsilon += strain;
  return 0;
}

const Matrix&
ElasticIsotropicPlaneStress2D::getTangent (void)
{
    double d00 = E/(1.0-v*v);
    double d01 = v*d00;
    double d22 = 0.5*(d00-d01);

    D(0,0) = D(1,1) = d00;
    D(1,0) = D(0,1) = d01;
    D(2,2) = d22;

    return D;
}

const Matrix&
ElasticIsotropicPlaneStress2D::getInitialTangent (void)
{
    double d00 = E/(1.0-v*v);
    double d01 = v*d00;
    double d22 = 0.5*(d00-d01);

    D(0,0) = D(1,1) = d00;
    D(1,0) = D(0,1) = d01;
    D(2,2) = d22;

    return D;
}

const Vector&
ElasticIsotropicPlaneStress2D::getStress (void)
{
    double d00 = E/(1.0-v*v);
    double d01 = v*d00;
    double d22 = 0.5*(d00-d01);

    double eps0 = epsilon(0);
    double eps1 = epsilon(1);

    //sigma = D*epsilon;
    sigma(0) = d00*eps0 + d01*eps1;
    sigma(1) = d01*eps0 + d00*eps1;
    sigma(2) = d22*epsilon(2);
	
    return sigma;
}

const Vector&
ElasticIsotropicPlaneStress2D::getStrain (void)
{
  return epsilon;
}

int
ElasticIsotropicPlaneStress2D::commitState (void)
{
  return 0;
}

int
ElasticIsotropicPlaneStress2D::revertToLastCommit (void)
{
  return 0;
}

int
ElasticIsotropicPlaneStress2D::revertToStart (void)
{
  epsilon.Zero();
  return 0;
}

NDMaterial*
ElasticIsotropicPlaneStress2D::getCopy (void)
{
	ElasticIsotropicPlaneStress2D *theCopy =
		new ElasticIsotropicPlaneStress2D (this->getTag(), E, v, rho);

	theCopy->epsilon = epsilon;

	return theCopy;
}

const char*
ElasticIsotropicPlaneStress2D::getType (void) const
{
	return "PlaneStress";
}

int
ElasticIsotropicPlaneStress2D::getOrder (void) const
{
	return 3;
}

