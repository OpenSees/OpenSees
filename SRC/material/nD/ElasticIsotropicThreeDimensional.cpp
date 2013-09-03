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
                                                                        
#include <ElasticIsotropicThreeDimensional.h>           
#include <Channel.h>

Vector ElasticIsotropicThreeDimensional::sigma(6);
Matrix ElasticIsotropicThreeDimensional::D(6,6);

ElasticIsotropicThreeDimensional::ElasticIsotropicThreeDimensional
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicThreeDimensional, E, nu, rho),
 epsilon(6), Cepsilon(6)
{
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticIsotropicThreeDimensional::ElasticIsotropicThreeDimensional():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicThreeDimensional, 0.0, 0.0),
 epsilon(6), Cepsilon(6)
{
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticIsotropicThreeDimensional::~ElasticIsotropicThreeDimensional ()
{

}

int
ElasticIsotropicThreeDimensional::setTrialStrain (const Vector &strain)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropicThreeDimensional::setTrialStrain (const Vector &strain, const Vector &rate)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropicThreeDimensional::setTrialStrainIncr (const Vector &strain)
{
  epsilon += strain;
  return 0;
}

int
ElasticIsotropicThreeDimensional::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsilon += strain;
  return 0;
}

const Matrix&
ElasticIsotropicThreeDimensional::getTangent (void)
{
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  return D;
}

const Matrix&
ElasticIsotropicThreeDimensional::getInitialTangent (void)
{
  //  return this->getTangent();
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  return D;
}

const Vector&
ElasticIsotropicThreeDimensional::getStress (void)
{
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  double eps0 = epsilon(0);
  double eps1 = epsilon(1);
  double eps2 = epsilon(2);

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  //sigma = D*epsilon;
  sigma(0) = mu2*eps0 + lam*(eps1+eps2);
  sigma(1) = mu2*eps1 + lam*(eps0+eps2);
  sigma(2) = mu2*eps2 + lam*(eps0+eps1);
  sigma(3) = mu*epsilon(3);
  sigma(4) = mu*epsilon(4);
  sigma(5) = mu*epsilon(5);
	
  return sigma;
}

const Vector&
ElasticIsotropicThreeDimensional::getStrain (void)
{
  return epsilon;
}

int
ElasticIsotropicThreeDimensional::commitState (void)
{
  Cepsilon=epsilon;
  return 0;
}

int
ElasticIsotropicThreeDimensional::revertToLastCommit (void)
{
  epsilon=Cepsilon;
  return 0;
}

int
ElasticIsotropicThreeDimensional::revertToStart (void)
{
  epsilon.Zero();
  Cepsilon.Zero();
  return 0;
}

NDMaterial*
ElasticIsotropicThreeDimensional::getCopy (void)
{
  ElasticIsotropicThreeDimensional *theCopy =
    new ElasticIsotropicThreeDimensional (this->getTag(), E, v, rho);
  
  theCopy->epsilon = epsilon;
  theCopy->Cepsilon = Cepsilon;
  
  return theCopy;
}

const char*
ElasticIsotropicThreeDimensional::getType (void) const
{
  return "ThreeDimensional";
}

int
ElasticIsotropicThreeDimensional::getOrder (void) const
{
  return 6;
}

int 
ElasticIsotropicThreeDimensional::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(10);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = v;
  data(3) = rho;
  data(4) = Cepsilon(0);
  data(5) = Cepsilon(1);
  data(6) = Cepsilon(2);
  data(7) = Cepsilon(3);
  data(8) = Cepsilon(4);
  data(9) = Cepsilon(5);
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticIsotropicThreeDimensional::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
ElasticIsotropicThreeDimensional::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(10);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticIsotropicThreeDimensional::sendSelf -- could not send Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  E = data(1);
  v = data(2);
  rho = data(3);
  epsilon(0)=data(4);
  epsilon(1)=data(5);
  epsilon(2)=data(6);
  epsilon(3)=data(7);
  epsilon(4)=data(8);
  epsilon(5)=data(9);

  Cepsilon = epsilon;

  return res;
}


const Vector&
ElasticIsotropicThreeDimensional::getStressSensitivity(int gradIndex,
						       bool conditional)
{
  if (parameterID < 1 || parameterID > 2) {
    sigma.Zero();
    return sigma;
  }

  double dmu2dh = 0.0;
  double dlamdh = 0.0;

  //double lam = v*mu2/(1.0-2.0*v);

  if (parameterID == 1) { // E
    //double mu2 = E/(1.0+v);
    dmu2dh = 1.0/(1.0+v);
    dlamdh = dmu2dh*v/(1-2*v);
  }

  if (parameterID == 2) { // nu
    double mu2 = E/(1.0+v);
    dmu2dh = -E/(1.0 + 2*v + v*v);
    dlamdh = mu2/(1.0 - 4*v + 4*v*v) + dmu2dh*v/(1-2*v);
  }

  //double mu = 0.50*mu2;
  double dmudh  = 0.5*dmu2dh;

  //mu2 += lam;
  dmu2dh += dlamdh;

  double eps0 = epsilon(0);
  double eps1 = epsilon(1);
  double eps2 = epsilon(2);

  //sigma = D*epsilon;
  sigma(0) = dmu2dh*eps0 + dlamdh*(eps1+eps2);
  sigma(1) = dmu2dh*eps1 + dlamdh*(eps0+eps2);
  sigma(2) = dmu2dh*eps2 + dlamdh*(eps0+eps1);
  sigma(3) = dmudh*epsilon(3);
  sigma(4) = dmudh*epsilon(4);
  sigma(5) = dmudh*epsilon(5);

  return sigma;
}
