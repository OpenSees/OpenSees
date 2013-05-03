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
                                                                        
#include <ElasticOrthotropicThreeDimensional.h>           
#include <Channel.h>

Vector ElasticOrthotropicThreeDimensional::sigma(6);
Matrix ElasticOrthotropicThreeDimensional::D(6,6);

ElasticOrthotropicThreeDimensional::ElasticOrthotropicThreeDimensional
(int tag, double Ex, double Ey, double Ez,
double vxy, double vyz, double vzx,
double Gxy, double Gyz, double Gzx, double rho) :
 ElasticOrthotropicMaterial (tag, 
ND_TAG_ElasticOrthotropicThreeDimensional, Ex, Ey, Ez, vxy, vyz, vzx, Gxy, Gyz, Gzx, rho),
 epsilon(6), Cepsilon(6)
{
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticOrthotropicThreeDimensional::ElasticOrthotropicThreeDimensional():
 ElasticOrthotropicMaterial (0, ND_TAG_ElasticOrthotropicThreeDimensional, 
0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0),
 epsilon(6), Cepsilon(6)
{
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticOrthotropicThreeDimensional::~ElasticOrthotropicThreeDimensional ()
{

}

int
ElasticOrthotropicThreeDimensional::setTrialStrain (const Vector &strain)
{
  epsilon = strain;
  return 0;
}

int
ElasticOrthotropicThreeDimensional::setTrialStrain (const Vector &strain, const Vector &rate)
{
  epsilon = strain;
  return 0;
}

int
ElasticOrthotropicThreeDimensional::setTrialStrainIncr (const Vector &strain)
{
  epsilon += strain;
  return 0;
}

int
ElasticOrthotropicThreeDimensional::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsilon += strain;
  return 0;
}

const Matrix&
ElasticOrthotropicThreeDimensional::getTangent (void)
{
  double vyx = vxy*Ey/Ex;
  double vzy = vyz*Ez/Ey;
  double vxz = vzx*Ex/Ez;

  double d = (1.0 - vxy*vyx-vyz*vzy-vzx*vxz - 2.0*vxy*vyz*vzx)/(Ex*Ey*Ez);

  D(0,0) = (1.0-vyz*vzy)/(Ey*Ez*d);
  D(1,1) = (1.0-vzx*vxz)/(Ez*Ex*d);
  D(2,2) = (1.0-vxy*vyx)/(Ex*Ey*d);

  D(1,0) = (vxy + vxz*vzy)/(Ez*Ex*d);
  D(0,1) = D(1,0);

  D(2,0) = (vxz + vxy*vyz)/(Ex*Ey*d);
  D(0,2) = D(2,0);

  D(2,1) = (vyz + vxz*vyx)/(Ex*Ey*d);
  D(1,2) = D(2,1);

  D(3,3) = Gxy;
  D(4,4) = Gyz;
  D(5,5) = Gzx;

  return D;
}

const Matrix&
ElasticOrthotropicThreeDimensional::getInitialTangent (void)
{
  return this->getTangent();
}

const Vector&
ElasticOrthotropicThreeDimensional::getStress (void)
{
  double eps0 = epsilon(0);
  double eps1 = epsilon(1);
  double eps2 = epsilon(2);

  double vyx = vxy*Ey/Ex;
  double vzy = vyz*Ez/Ey;
  double vxz = vzx*Ex/Ez;

  double d = (1.0 - vxy*vyx-vyz*vzy-vzx*vxz - 2.0*vxy*vyz*vzx)/(Ex*Ey*Ez);

  D(0,0) = (1.0-vyz*vzy)/(Ey*Ez*d);
  D(1,1) = (1.0-vzx*vxz)/(Ez*Ex*d);
  D(2,2) = (1.0-vxy*vyx)/(Ex*Ey*d);

  D(1,0) = (vxy + vxz*vzy)/(Ez*Ex*d);
  D(0,1) = D(1,0);

  D(2,0) = (vxz + vxy*vyz)/(Ex*Ey*d);
  D(0,2) = D(2,0);

  D(2,1) = (vyz + vxz*vyx)/(Ex*Ey*d);
  D(1,2) = D(2,1);

  sigma(0) = D(0,0)*eps0 + D(0,1)*eps1 + D(0,2)*eps2;
  sigma(1) = D(1,0)*eps0 + D(1,1)*eps1 + D(1,2)*eps2;
  sigma(2) = D(2,0)*eps0 + D(2,1)*eps1 + D(2,2)*eps2;

  sigma(3) = Gxy*epsilon(3);
  sigma(4) = Gyz*epsilon(4);
  sigma(5) = Gzx*epsilon(5);
	
  return sigma;
}

const Vector&
ElasticOrthotropicThreeDimensional::getStrain (void)
{
  return epsilon;
}

int
ElasticOrthotropicThreeDimensional::commitState (void)
{
  Cepsilon=epsilon;
  return 0;
}

int
ElasticOrthotropicThreeDimensional::revertToLastCommit (void)
{
  epsilon=Cepsilon;
  return 0;
}

int
ElasticOrthotropicThreeDimensional::revertToStart (void)
{
  epsilon.Zero();
  Cepsilon.Zero();
  return 0;
}

NDMaterial*
ElasticOrthotropicThreeDimensional::getCopy (void)
{
  ElasticOrthotropicThreeDimensional *theCopy =
    new ElasticOrthotropicThreeDimensional (this->getTag(), Ex, Ey, Ez, 
vxy, vyz, vzx, Gxy, Gyz, Gzx, rho);
  
  theCopy->epsilon = epsilon;
  theCopy->Cepsilon = Cepsilon;
  
  return theCopy;
}

const char*
ElasticOrthotropicThreeDimensional::getType (void) const
{
  return "ThreeDimensional";
}

int
ElasticOrthotropicThreeDimensional::getOrder (void) const
{
  return 6;
}

int 
ElasticOrthotropicThreeDimensional::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(17);
  
  data(0) = this->getTag();
  data(1) = Ex;
  data(2) = Ey;
  data(3) = Ez;
  data(4) = vxy;
  data(5) = vyz;
  data(6) = vzx;
  data(7) = Gxy;
  data(8) = Gyz;
  data(9) = Gzx;
  data(10) = rho;
  data(11) = Cepsilon(0);
  data(12) = Cepsilon(1);
  data(13) = Cepsilon(2);
  data(14) = Cepsilon(3);
  data(15) = Cepsilon(4);
  data(16) = Cepsilon(5);
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticOrthotropicThreeDimensional::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
ElasticOrthotropicThreeDimensional::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(17);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticOrthotropicThreeDimensional::sendSelf -- could not send Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  Ex = data(1);
  Ey = data(2);
  Ez = data(3);
  vxy = data(4);
  vyz = data(5);
  vzx = data(6);
  Gxy = data(7);
  Gyz = data(8);
  Gzx = data(9);
  rho = data(10);
  epsilon(0)=data(11);
  epsilon(1)=data(12);
  epsilon(2)=data(13);
  epsilon(3)=data(14);
  epsilon(4)=data(15);
  epsilon(5)=data(16);

  Cepsilon = epsilon;

  return res;
}


const Vector&
ElasticOrthotropicThreeDimensional::getStressSensitivity(int gradIndex,
						       bool conditional)
{
  if (parameterID < 1 || parameterID > 2) {
    sigma.Zero();
    return sigma;
  }

  sigma.Zero();
  return sigma;
}
