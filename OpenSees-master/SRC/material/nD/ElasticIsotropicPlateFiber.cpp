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
                                                                        
// $Revision: 1.6 $
// $Date: 2002-12-05 22:49:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlateFiber.cpp,v $
                                                                        
                                                                        

#include <ElasticIsotropicPlateFiber.h>           
#include <Channel.h>

Vector ElasticIsotropicPlateFiber::sigma(5);
Matrix ElasticIsotropicPlateFiber::D(5,5);

ElasticIsotropicPlateFiber::ElasticIsotropicPlateFiber
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicPlateFiber, E, nu, rho),
 epsilon(5)
{

}

ElasticIsotropicPlateFiber::ElasticIsotropicPlateFiber():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicPlateFiber, 0.0, 0.0),
 epsilon(5)
{

}

ElasticIsotropicPlateFiber::~ElasticIsotropicPlateFiber ()
{

}

int
ElasticIsotropicPlateFiber::setTrialStrain (const Vector &strain)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropicPlateFiber::setTrialStrain (const Vector &strain, const Vector &rate)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropicPlateFiber::setTrialStrainIncr (const Vector &strain)
{
  epsilon += strain;
  return 0;
}

int
ElasticIsotropicPlateFiber::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsilon += strain;
  return 0;
}

const Matrix&
ElasticIsotropicPlateFiber::getTangent (void)
{
    double d00 = E/(1.0-v*v);
    double d01 = v*d00;
    double d22 = 0.5*(d00-d01);
    
    D(0,0) = D(1,1) = d00;
    D(0,1) = D(1,0) = d01;
    D(2,2) = d22;
    D(3,3) = d22;
    D(4,4) = d22;

    return D;
}

const Matrix&
ElasticIsotropicPlateFiber::getInitialTangent (void)
{
    double d00 = E/(1.0-v*v);
    double d01 = v*d00;
    double d22 = 0.5*(d00-d01);
    
    D(0,0) = D(1,1) = d00;
    D(0,1) = D(1,0) = d01;
    D(2,2) = d22;
    D(3,3) = d22;
    D(4,4) = d22;

    return D;
}

const Vector&
ElasticIsotropicPlateFiber::getStress (void)
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
    sigma(3) = d22*epsilon(3);
    sigma(4) = d22*epsilon(4);
	
    return sigma;
}

const Vector&
ElasticIsotropicPlateFiber::getStrain (void)
{
  return epsilon;
}

int
ElasticIsotropicPlateFiber::commitState (void)
{
  return 0;
}

int
ElasticIsotropicPlateFiber::revertToLastCommit (void)
{
  return 0;
}

int
ElasticIsotropicPlateFiber::revertToStart (void)
{
  epsilon.Zero();
  return 0;
}

NDMaterial*
ElasticIsotropicPlateFiber::getCopy (void)
{
	ElasticIsotropicPlateFiber *theCopy =
		new ElasticIsotropicPlateFiber (this->getTag(), E, v, rho);

	theCopy->epsilon = epsilon;

	return theCopy;
}

const char*
ElasticIsotropicPlateFiber::getType (void) const
{
	return "PlateFiber";
}

int
ElasticIsotropicPlateFiber::getOrder (void) const
{
	return 5;
}
