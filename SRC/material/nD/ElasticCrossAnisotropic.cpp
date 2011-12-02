//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Elastic Cross Anisotropic Material implementation:
//# CLASS:             ElasticIsotropic3D
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhaohui Yang, Yi Bian, Boris Jeremic
//#
//#
//# DATE:              10Oct2002 Yi Bian
//# UPDATE HISTORY:    March 20, 2003 Revised by Joey Yang & Boris Jeremic, UC Davis
//#
//#
//===============================================================================
//


#include <ElasticCrossAnisotropic.h>

//Tensor ElasticCrossAnisotropic :: rank2 (2, def_dim_2, 0.0 ) ;
//Tensor ElasticCrossAnisotropic :: rank4 (2, def_dim_2, 0.0 ) ;

Vector ElasticCrossAnisotropic::sigma(6);
Matrix ElasticCrossAnisotropic::D(6, 6);

///////////////////////////////////////////////////////////////////////////////
ElasticCrossAnisotropic::ElasticCrossAnisotropic(int tag, 
                                                 double Ehp, 
                                                 double Evp, 
                                                 double nuhvp, 
                                                 double nuhhp, 
                                                 double Ghvp, 
                                                 double rhop):
NDMaterial(tag, ND_TAG_ElasticCrossAnisotropic3D),
Tepsilon(6), 
Cepsilon(6),
Eh(Ehp), 
Ev(Evp), 
nuhv(nuhvp),
nuhh(nuhhp),
Ghv(Ghvp),
rho(rhop)
{
   D.Zero();
   Dt = tensor( 4, def_dim_4, 0.0 );
   this->convertD2TensorEijkl();
}

///////////////////////////////////////////////////////////////////////////////
ElasticCrossAnisotropic::ElasticCrossAnisotropic():
Tepsilon(6), Cepsilon(6)
{
  D.Zero();
}

///////////////////////////////////////////////////////////////////////////////
ElasticCrossAnisotropic::~ElasticCrossAnisotropic ()
{

}

///////////////////////////////////////////////////////////////////////////////
double ElasticCrossAnisotropic::getrho()
{
  return rho;
}

///////////////////////////////////////////////////////////////////////////////
NDMaterial* ElasticCrossAnisotropic::getCopy (const char *type)
{
    if (strcmp(type,"ThreeDimensional") == 0)
    {
  ElasticCrossAnisotropic *theModel;
  theModel = new ElasticCrossAnisotropic (this->getTag(), Eh, Ev, nuhv, nuhh, Ghv, rho);
          // This function should only be called during element instantiation, so
    // no state determination is performed on the material model object
    // prior to copying the material model (calling this function)
  return theModel;
    }

    else
    {
  opserr <<"ElasticCrossAnisotropic::getModel failed to get model " << type << "\n";
  return 0;
    }
}

///////////////////////////////////////////////////////////////////////////////
int ElasticCrossAnisotropic::setTrialStrain (const Vector &strain)
{
  Tepsilon = strain;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
int ElasticCrossAnisotropic::setTrialStrain (const Vector &strain, const Vector &rate)
{
  Tepsilon = strain;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
int ElasticCrossAnisotropic::setTrialStrainIncr (const Vector &strain)
{
  //epsilon += strain;
  Tepsilon = Cepsilon;
  Tepsilon.addVector (1.0, strain, 1.0);

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
int ElasticCrossAnisotropic::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  //epsilon += strain;
  Tepsilon = Cepsilon;
  Tepsilon.addVector (1.0, strain, 1.0);

  return 0;
}

const Matrix&
ElasticCrossAnisotropic::getTangent (void)
{
    this->setInitElasticStiffness();
    return D;
}

const Vector&
ElasticCrossAnisotropic::getStress (void)
{
    double eps0 = Tepsilon(0);
    double eps1 = Tepsilon(1);
    double eps2 = Tepsilon(2);
    double eps3 = Tepsilon(3);
    double eps4 = Tepsilon(4);
    double eps5 = Tepsilon(5);

    //sigma = D*epsilon;
    sigma(0) = D(0,0) * eps0 + D(0,1) * eps1 + D(0,2) * eps2;
    sigma(1) = D(1,0) * eps0 + D(1,1) * eps1 + D(1,2) * eps2;
    sigma(2) = D(2,0) * eps0 + D(2,1) * eps1 + D(2,2) * eps2;
    sigma(3) = D(3,3) * eps3;
    sigma(4) = D(4,4) * eps4;
    sigma(5) = D(5,5) * eps5;

    return sigma;
}

const Vector&
ElasticCrossAnisotropic::getStrain (void)
{
    return Tepsilon;
}

int
ElasticCrossAnisotropic::setTrialStrain (const Tensor &v)
{
    Strain = v;
    return 0;
}

int
ElasticCrossAnisotropic::setTrialStrain (const Tensor &v, const Tensor &r)
{
    Strain = v;
    return 0;
}

int
ElasticCrossAnisotropic::setTrialStrainIncr (const Tensor &v)
{
    Strain = Strain + v;
    return 0;
}

int
ElasticCrossAnisotropic::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    Strain = Strain + v;
    return 0;
}

const Tensor&
ElasticCrossAnisotropic::getTangentTensor (void)
{
   //Dt.print();
   return Dt;
}

const stresstensor ElasticCrossAnisotropic::getStressTensor (void)
{
        Stress = Dt("ijkl") * Strain("kl");
        return Stress;
}

const straintensor ElasticCrossAnisotropic::getStrainTensor (void)
{
  return Strain;
}

const straintensor
ElasticCrossAnisotropic::getPlasticStrainTensor (void)
{
     //Return zero straintensor
     straintensor t;
     return t;
}

int
ElasticCrossAnisotropic::commitState (void)
{
  //sigma = D * epsilon;
  Cepsilon = Tepsilon;

        return 0;
}

int
ElasticCrossAnisotropic::revertToLastCommit (void)
{
  Tepsilon = Cepsilon;
  return 0;
}

int
ElasticCrossAnisotropic::revertToStart (void)
{
  Cepsilon.Zero();
  return 0;
}

NDMaterial*
ElasticCrossAnisotropic::getCopy (void)
{
  ElasticCrossAnisotropic *theCopy =
  new ElasticCrossAnisotropic (this->getTag(), Eh, Ev, nuhv, nuhh, Ghv, rho);

  theCopy->Cepsilon = this -> Cepsilon;
  theCopy->sigma = this->sigma;
  //theCopy->Strain = this->Strain;
  //theCopy->Stress = this->Stress;
  return theCopy;
}

const char*
ElasticCrossAnisotropic::getType (void) const
{
  return "ThreeDimensional";
}

int
ElasticCrossAnisotropic::getOrder (void) const
{
  return 6;
}

int
ElasticCrossAnisotropic::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(13);

  data(0) = this->getTag();
  data(1) = Eh;
  data(2) = Ev;
  data(3) = nuhv;
  data(4) = nuhh;
  data(5) = Ghv;
        data(6) = rho;
  data(7) = Cepsilon(0);
  data(8) = Cepsilon(1);
  data(9) = Cepsilon(2);
        data(10) = Cepsilon(3);
  data(11) = Cepsilon(4);
  data(12) = Cepsilon(5);

        res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << " ElasticCrossAnisotropic::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}

int
ElasticCrossAnisotropic::recvSelf(int commitTag, Channel &theChannel,
     FEM_ObjectBroker &theBroker)
{
  int res = 0;

        static Vector data(13);

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "ElasticCrossAnisotropic::recvSelf -- could not receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));
        Eh = data(1);
  Ev = data(2);
  nuhv = data(3);
  nuhh = data(4);
  Ghv = data(5);
  rho = data(6);

  Cepsilon(0) = data(7);
  Cepsilon(1) = data(8);
  Cepsilon(2) = data(9);
        Cepsilon(3) = data(10);
  Cepsilon(4) = data(11);
  Cepsilon(5) = data(12);

  D.Zero();
  return res;
}

void
ElasticCrossAnisotropic::Print (OPS_Stream &s, int flag)
{
  s << "Elastic Cross-Anisotropic Material Model\n";
  s << "\tEh:  " << Eh << "\tEv:  " << Ev << "\n";
  s << "\tnuhv:  " << nuhv << "\tnuhh:  " << nuhh << "\n";
  s << "\tGhv:  " << Ghv << "\trho:  " << rho << "\n";
  return;
}

//int
//ElasticCrossAnisotropic::setParameter(char **argv, int argc, Information &info)
//{
//  return -1;
//}

/*int
ElasticCrossAnisotropic::updateParameter(int parameterID, Information &info)
{
  return -1;
}
*/

void ElasticCrossAnisotropic::setInitElasticStiffness(void)
{
   //Old codes from Yi Bian
   //double A = 1/((1 + nuhv) * (1 - 2 * nuhv));
   //D(0, 0) = A * (1 - nuhv) * Ev;
   //D(1, 1) = D(2, 2) = A * (1 - nuhv) * Eh;
   //D(0, 1) = D(0, 2) = D(1, 0) = D(2, 0) = A * Eh * nuhh;
   //D(3, 3) = Eh/(1 + nuhv);
   //D(4, 4) = D(5, 5) = 2 * Ghv;

//  Compliance matrix C Refer to Thor C. Heimdahl and Andrew Drescher
//  "Elastic Anisotropy of Tire Shreds", {ASCE} Journal of Geotechnical
//  and Geoenvironmental Engineering, 1999, pp. 383-389
//  |Sxx|    | 1/Eh     -nuhh/Eh -nuhv/Ev     0       0     0    |
//  |Syy|    |-nuhh/Eh    1/Eh   -nuhv/Ev     0       0     0    |
//  |Szz|    |-nuhv/Ev  -nuhv/Ev   1/Ev       0       0     0    |
//  |Sxy| C= |   0         0       0  2(1+nuhh)/Eh     0     0    |
//  |Sxz|    |   0         0       0         0     1/(Ghv)   0    |
//  |Syz|    |   0         0       0         0         0   1/(Ghv)|
//
   // Form compliance matrix D
   double A = 1.0/Eh;
   double B = 1.0/Ev;
   //opserr << "A " << A << " B " << B;
   D(0,0) = D(1,1) = A;
   D(2,2) = B;
   D(0,1) = D(1,0) = -nuhh*A;
   D(0,2) = D(2,0) = D(1,2) = D(2,1) = -nuhv*B;
   D(3,3) = 2*(1.0+nuhh)*A;
   //D(4,4) = D(5,5) = 0.5/Ghv;
   D(4,4) = D(5,5) = 1.0/Ghv;
   //opserr << " C " << D;

   // Do invertion once to get Elastic matrix D
   D.Invert( D );
   //opserr << " D " << D;

   return;

}

void ElasticCrossAnisotropic::convertD2TensorEijkl(void)
{
   this->setInitElasticStiffness();

   //Convert Matric D to 4th order Elastic constants tensor Dt;
   Dt.val(1,1,1,1) = D(0,0); 
   Dt.val(1,1,2,2) = D(0,1); 
   Dt.val(1,1,3,3) = D(0,2); // --> Sigma_xx

   Dt.val(1,2,1,2) = D(3,3); 
   Dt.val(1,2,2,1) = D(3,3); // --> Sigma_xy

   Dt.val(1,3,1,3) = D(4,4); 
   Dt.val(1,3,3,1) = D(4,4); // --> Sigma_xz

   Dt.val(2,1,1,2) = D(3,3); 
   Dt.val(2,1,2,1) = D(3,3); // --> Sigma_yx

   Dt.val(2,2,1,1) = D(1,0); 
   Dt.val(2,2,2,2) = D(1,1); 
   Dt.val(2,2,3,3) = D(1,2); // --> Sigma_yy

   Dt.val(2,3,2,3) = D(5,5); 
   Dt.val(2,3,3,2) = D(5,5); // --> Sigma_yz

   Dt.val(3,1,1,3) = D(4,4); 
   Dt.val(3,1,3,1) = D(4,4); // --> Sigma_zx

   Dt.val(3,2,2,3) = D(5,5); 
   Dt.val(3,2,3,2) = D(5,5); // --> Sigma_zy

   Dt.val(3,3,1,1) = D(2,0); 
   Dt.val(3,3,2,2) = D(2,1); 
   Dt.val(3,3,3,3) = D(2,2); // --> Sigma_zz

   return;
}
