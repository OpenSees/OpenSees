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

Vector ElasticIsotropicPlaneStress2D::sigma(3);
Matrix ElasticIsotropicPlaneStress2D::D(3,3);

ElasticIsotropicPlaneStress2D::ElasticIsotropicPlaneStress2D
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicPlaneStress2d, E, nu, rho),
 epsilon(3), Cepsilon(3)
{
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticIsotropicPlaneStress2D::ElasticIsotropicPlaneStress2D():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicPlaneStress2d, 0.0, 0.0),
 epsilon(3), Cepsilon(3)
{
  epsilon.Zero();
  Cepsilon.Zero();
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
  Cepsilon=epsilon;
  return 0;
}

int
ElasticIsotropicPlaneStress2D::revertToLastCommit (void)
{
  epsilon=Cepsilon;
  return 0;
}

int
ElasticIsotropicPlaneStress2D::revertToStart (void)
{
  epsilon.Zero();
  Cepsilon.Zero();
  return 0;
}

NDMaterial*
ElasticIsotropicPlaneStress2D::getCopy (void)
{
  ElasticIsotropicPlaneStress2D *theCopy =
    new ElasticIsotropicPlaneStress2D (this->getTag(), E, v, rho);
  
  theCopy->epsilon = epsilon;
  theCopy->Cepsilon = Cepsilon;
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

int 
ElasticIsotropicPlaneStress2D::sendSelf(int commitTag, Channel &theChannel)
{
  
  static Vector data(7);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = v;
  data(3) = rho;
  data(4) = Cepsilon(0);
  data(5) = Cepsilon(1);
  data(6) = Cepsilon(2);
  
 int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
 if (res < 0) {
   opserr << "ElasticIsotropicPlaneStress2D::sendSelf -- could not send Vector\n";
   return res;
 }

 return res;
}

int 
ElasticIsotropicPlaneStress2D::recvSelf(int commitTag, Channel &theChannel, 
				      FEM_ObjectBroker &theBroker)
{
  static Vector data(7);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticIsotropicPlaneStress2D::sendSelf -- could not send Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  E = data(1);
  v = data(2);
  rho = data(3);
  epsilon(0)=data(4);
  epsilon(1)=data(5);
  epsilon(2)=data(6);
  Cepsilon = epsilon;

  return res;
}
