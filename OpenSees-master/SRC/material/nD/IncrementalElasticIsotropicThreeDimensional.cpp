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
                                                                        
#include <Channel.h>
#include <IncrementalElasticIsotropicThreeDimensional.h>           
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <MaterialResponse.h>

bool IncrementalElasticIsotropicThreeDimensional::printnow = true;
// Vector IncrementalElasticIsotropicThreeDimensional::sigma(6);
Matrix IncrementalElasticIsotropicThreeDimensional::D(6,6);

void *OPS_IncrementalElasticIsotropicThreeDimensional(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  

  if (numArgs < 3) {
    opserr << "Want: nDMaterial IncrementalElasticIsotropic3D $tag $E $V <$rho>" << endln;
    return 0; 
  }
  
  int iData[1];
  double dData[3];
  dData[2] = 0.0;
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: nDMaterial IncrementalElasticIsotropic3D \n";
    return 0;
  }
  
  if (numArgs > 3) 
    numData = 3;
  else
    numData = 2;
  
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid data: nDMaterial IncrementalElasticIsotropic3D : " << iData[0] <<"\n";
    return 0;
  }  
  
  theMaterial = new IncrementalElasticIsotropicThreeDimensional(iData[0], dData[0], dData[1], dData[2]);
  
  return theMaterial;
}



IncrementalElasticIsotropicThreeDimensional::IncrementalElasticIsotropicThreeDimensional(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_IncrementalElasticIsotropicThreeDimensional, E, nu, rho),
 epsilon(6), epsilon_n(6), sigma(6), sigma_n(6)
{
  epsilon.Zero();
  sigma.Zero();
  sigma_n.Zero();
  epsilon_n.Zero();
}

IncrementalElasticIsotropicThreeDimensional::IncrementalElasticIsotropicThreeDimensional():
 ElasticIsotropicMaterial (0, ND_TAG_IncrementalElasticIsotropicThreeDimensional, 0.0, 0.0),
 epsilon(6), epsilon_n(6), sigma(6), sigma_n(6)
{
  epsilon.Zero();
  sigma.Zero();
  sigma_n.Zero();
  epsilon_n.Zero();
}

IncrementalElasticIsotropicThreeDimensional::~IncrementalElasticIsotropicThreeDimensional ()
{

}

int
IncrementalElasticIsotropicThreeDimensional::setTrialStrain (const Vector &strain)
{
  epsilon = strain;

  // static Vector depsilon(6);
  // depsilon.Zero();
  
  // sigma = sigma_n;


  // depsilon = epsilon - epsilon_n;

  // double mu2 = E/(1.0+v);
  // double lam = v*mu2/(1.0-2.0*v);
  // double mu = 0.50*mu2;
  // mu2 += lam;

  // double deps0 = depsilon(0);
  // double deps1 = depsilon(1);
  // double deps2 = depsilon(2);

  // D(0,0) = D(1,1) = D(2,2) = mu2;
  // D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  // D(3,3) = mu;
  // D(4,4) = mu;
  // D(5,5) = mu;

  // //dsigma = D*depsilon;
  // sigma(0) = sigma(0) + mu2*deps0 + lam*(deps1 + deps2);
  // sigma(1) = sigma(1) + mu2*deps1 + lam*(deps0 + deps2);
  // sigma(2) = sigma(2) + mu2*deps2 + lam*(deps0 + deps1);
  // sigma(3) = sigma(3) + mu*depsilon(3);
  // sigma(4) = sigma(4) + mu*depsilon(4);
  // sigma(5) = sigma(5) + mu*depsilon(5);

  return 0;
}

int
IncrementalElasticIsotropicThreeDimensional::setTrialStrain (const Vector &strain, const Vector &rate)
{
  epsilon = strain;

  // static Vector depsilon(6);
  // depsilon.Zero();
  
  // sigma = sigma_n;


  // depsilon = epsilon - epsilon_n;

  // double mu2 = E/(1.0+v);
  // double lam = v*mu2/(1.0-2.0*v);
  // double mu = 0.50*mu2;
  // mu2 += lam;

  // double deps0 = depsilon(0);
  // double deps1 = depsilon(1);
  // double deps2 = depsilon(2);

  // D(0,0) = D(1,1) = D(2,2) = mu2;
  // D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  // D(3,3) = mu;
  // D(4,4) = mu;
  // D(5,5) = mu;

  // //dsigma = D*depsilon;
  // sigma(0) = sigma(0) + mu2*deps0 + lam*(deps1 + deps2);
  // sigma(1) = sigma(1) + mu2*deps1 + lam*(deps0 + deps2);
  // sigma(2) = sigma(2) + mu2*deps2 + lam*(deps0 + deps1);
  // sigma(3) = sigma(3) + mu*depsilon(3);
  // sigma(4) = sigma(4) + mu*depsilon(4);
  // sigma(5) = sigma(5) + mu*depsilon(5);

  return 0;
}

int
IncrementalElasticIsotropicThreeDimensional::setTrialStrainIncr (const Vector &strain)
{
  epsilon += strain;
  return 0;
}

int
IncrementalElasticIsotropicThreeDimensional::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsilon += strain;
  return 0;
}

const Matrix&
IncrementalElasticIsotropicThreeDimensional::getTangent (void)
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
IncrementalElasticIsotropicThreeDimensional::getInitialTangent (void)
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
IncrementalElasticIsotropicThreeDimensional::getStress (void)
{	
  static Vector depsilon(6);
  depsilon.Zero();
  
  sigma = sigma_n;

  depsilon = epsilon - epsilon_n;

  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  mu2 += lam;

  double deps0 = depsilon(0);
  double deps1 = depsilon(1);
  double deps2 = depsilon(2);

  D(0,0) = D(1,1) = D(2,2) = mu2;
  D(0,1) = D(1,0) = D(0,2) = D(2,0) = D(1,2) = D(2,1) = lam;
  D(3,3) = mu;
  D(4,4) = mu;
  D(5,5) = mu;

  //dsigma = D*depsilon;
  sigma(0) = sigma(0) + mu2*deps0 + lam*(deps1 + deps2);
  sigma(1) = sigma(1) + mu2*deps1 + lam*(deps0 + deps2);
  sigma(2) = sigma(2) + mu2*deps2 + lam*(deps0 + deps1);
  sigma(3) = sigma(3) + mu*depsilon(3);
  sigma(4) = sigma(4) + mu*depsilon(4);
  sigma(5) = sigma(5) + mu*depsilon(5);

  // if(printnow)
  // {
  //   // opserr << "   [[[[[[[[[  sigma = " << sigma << endln;
  //   printnow = false;
  // }

  return sigma;
}

const Vector&
IncrementalElasticIsotropicThreeDimensional::getStrain (void)
{
  return epsilon_n;
}

int
IncrementalElasticIsotropicThreeDimensional::commitState (void)
{
  epsilon_n=epsilon;
  sigma_n=sigma;
  // printnow = true;
  return 0;
}

int
IncrementalElasticIsotropicThreeDimensional::revertToLastCommit (void)
{
  epsilon=epsilon_n;
  sigma=sigma_n;
  return 0;
}

int
IncrementalElasticIsotropicThreeDimensional::revertToStart (void)
{
  epsilon.Zero();
  epsilon_n.Zero();
  sigma.Zero();
  sigma_n.Zero();
  return 0;
}

NDMaterial*
IncrementalElasticIsotropicThreeDimensional::getCopy (void)
{
  IncrementalElasticIsotropicThreeDimensional *theCopy =
    new IncrementalElasticIsotropicThreeDimensional (this->getTag(), E, v, rho);
  
  theCopy->epsilon = epsilon;
  theCopy->epsilon_n = epsilon_n;
  theCopy->sigma = sigma;
  theCopy->sigma_n = sigma_n;

  return theCopy;
}


NDMaterial*
IncrementalElasticIsotropicThreeDimensional::getCopy (const char *type)
{
  return this->getCopy();
}

const char*
IncrementalElasticIsotropicThreeDimensional::getType (void) const
{
  return "ThreeDimensional";
}

int
IncrementalElasticIsotropicThreeDimensional::getOrder (void) const
{
  return 6;
}

int 
IncrementalElasticIsotropicThreeDimensional::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(28);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = v;
  data(3) = rho;
  data(4) = epsilon(0);
  data(5) = epsilon(1);
  data(6) = epsilon(2);
  data(7) = epsilon(3);
  data(8) = epsilon(4);
  data(9) = epsilon(5);
  data(10) = epsilon_n(0);
  data(11) = epsilon_n(1);
  data(12) = epsilon_n(2);
  data(13) = epsilon_n(3);
  data(14) = epsilon_n(4);
  data(15) = epsilon_n(5);
  data(16) = sigma(0);
  data(17) = sigma(1);
  data(18) = sigma(2);
  data(19) = sigma(3);
  data(20) = sigma(4);
  data(21) = sigma(5);
  data(22) = sigma_n(0);
  data(23) = sigma_n(1);
  data(24) = sigma_n(2);
  data(25) = sigma_n(3);
  data(26) = sigma_n(4);
  data(27) = sigma_n(5);
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "IncrementalElasticIsotropicThreeDimensional::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
IncrementalElasticIsotropicThreeDimensional::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(28);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "IncrementalElasticIsotropicThreeDimensional::sendSelf -- could not send Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  E = data(1);
  v = data(2);
  rho = data(3);

  epsilon(0) = data(4) ;
  epsilon(1) = data(5) ;
  epsilon(2) = data(6) ;
  epsilon(3) = data(7) ;
  epsilon(4) = data(8) ;
  epsilon(5) = data(9) ;
  epsilon_n(0) = data(10) ;
  epsilon_n(1) = data(11) ;
  epsilon_n(2) = data(12) ;
  epsilon_n(3) = data(13) ;
  epsilon_n(4) = data(14) ;
  epsilon_n(5) = data(15) ;
  sigma(0) = data(16) ;
  sigma(1) = data(17) ;
  sigma(2) = data(18) ;
  sigma(3) = data(19) ;
  sigma(4) = data(20) ;
  sigma(5) = data(21) ;
  sigma_n(0) = data(22) ;
  sigma_n(1) = data(23) ;
  sigma_n(2) = data(24) ;
  sigma_n(3) = data(25) ;
  sigma_n(4) = data(26) ;
  sigma_n(5) = data(27) ;

  return res;
}


const Vector&
IncrementalElasticIsotropicThreeDimensional::getStressSensitivity(int gradIndex,
						       bool conditional)
{
  sigma.Zero();

  return sigma;
}


Response*
IncrementalElasticIsotropicThreeDimensional::setResponse (const char **argv, int argc, OPS_Stream &output)
{
    if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
    {
        return new MaterialResponse(this, 1, this->getStress());
    }
    else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
        return new MaterialResponse(this, 2, this->getStrain());
    else
        return 0;
}

int
IncrementalElasticIsotropicThreeDimensional::getResponse(int responseID, Information &matInfo)
{
    switch (responseID) {
        case -1:
            return -1;
        case 1:
            if (matInfo.theVector != 0)
                *(matInfo.theVector) = getStress();
            return 0;
        case 2:
            if (matInfo.theVector != 0)
                *(matInfo.theVector) = getStrain();
            return 0;
        default:
            return -1;
    }
}