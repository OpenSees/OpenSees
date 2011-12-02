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
//# PURPOSE:           Pressure dependent elastic isotropic material implementation:
//# CLASS:             PressureDependentElastic3D
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhaohui Yang, Boris Jeremic
//#
//#
//# DATE:              07July2001
//# UPDATE HISTORY:    22Nov2002 small fixes, formating...
//#
//#
//===============================================================================



#include <PressureDependentElastic3D.h>



Matrix PressureDependentElastic3D::D(6,6);   // global for ElasticIsotropic3D only
Vector PressureDependentElastic3D::sigma(6); // global for ElasticIsotropic3D only


PressureDependentElastic3D::PressureDependentElastic3D
(int tag, double E, double nu, double rhop, double expp, double pr, double pop):
 ElasticIsotropicMaterial (tag, ND_TAG_PressureDependentElastic3D, E, nu, rhop),
 epsilon(6), exp(expp), p_ref(pr), p_cutoff(pop)
{
  // Set up the elastic constant matrix for 3D elastic isotropic
  Dt = tensor( 4, def_dim_4, 0.0 );
  ComputeElasticStiffness();

}

PressureDependentElastic3D::PressureDependentElastic3D():
 ElasticIsotropicMaterial (0, ND_TAG_PressureDependentElastic3D, 0.0, 0.0, 0.0),
 epsilon(6)
{
  Dt = tensor( 4, def_dim_4, 0.0 );
}

PressureDependentElastic3D::~PressureDependentElastic3D ()
{

}

int PressureDependentElastic3D::setTrialStrain (const Vector &v)
  {
    epsilon = v;
    return 0;
  }

int PressureDependentElastic3D::setTrialStrain (const Vector &v, const Vector &r)
  {
    epsilon = v;
    return 0;
  }

int PressureDependentElastic3D::setTrialStrainIncr (const Vector &v)
  {
    epsilon += v;
    return 0;
  }

int PressureDependentElastic3D::setTrialStrainIncr (const Vector &v, const Vector &r)
  {
    epsilon += v;
    return 0;
  }

const Matrix& PressureDependentElastic3D::getTangent (void)
  {
    //Update E
    Stress = getStressTensor();
    double p = Stress.p_hydrostatic();

    if (p <= p_cutoff)
      p = p_cutoff;
    double Ec = E * pow(p/p_ref, exp);

    double mu2 = Ec/(1.0+v);
    double lam = v*mu2/(1.0-2.0*v);
    double mu  = 0.50*mu2;

    mu2 += lam;

    D(0,0) = D(1,1) = D(2,2) = mu2;
    D(0,1) = D(1,0) = lam;
    D(0,2) = D(2,0) = lam;
    D(1,2) = D(2,1) = lam;
    D(3,3) = mu;
    D(4,4) = mu;
    D(5,5) = mu;

    return D;
  }

const Vector& PressureDependentElastic3D::getStress (void)
  {
    //Update E
    Stress = getStressTensor();
    double p = Stress.p_hydrostatic();

    if (p <= p_cutoff)
    p = p_cutoff;
    double Ec = E * pow(p/p_ref, exp);

    double mu2 = Ec/(1.0+v);
    double lam = v*mu2/(1.0-2.0*v);
    double mu  = 0.50*mu2;

    mu2 += lam;

    double eps0 = epsilon(0);
    double eps1 = epsilon(1);
    double eps2 = epsilon(2);

    sigma(0) = mu2*eps0 + lam*(eps1+eps2);
    sigma(1) = mu2*eps1 + lam*(eps2+eps0);
    sigma(2) = mu2*eps2 + lam*(eps0+eps1);

    sigma(3) = mu*epsilon(3);
    sigma(4) = mu*epsilon(4);
    sigma(5) = mu*epsilon(5);

    return sigma;
    //sigma = D*epsilon;
    //return sigma;
  }

const Vector& PressureDependentElastic3D::getStrain (void)
  {
    return epsilon;
  }

int PressureDependentElastic3D::setTrialStrain (const Tensor &v)
  {
    Strain = v;
    return 0;
  }

int PressureDependentElastic3D::setTrialStrain (const Tensor &v, const Tensor &r)
  {
    Strain = v;
    return 0;
  }

int PressureDependentElastic3D::setTrialStrainIncr (const Tensor &v)
  {
    //opserr << " before set Tri St Incr " << Strain;
    //opserr << " Strain Incr " << v << endlnn;
    Strain = Strain + v;
    //opserr << " after setTrialStrainIncr  " << Strain << endlnn;
    return 0;
  }

int PressureDependentElastic3D::setTrialStrainIncr (const Tensor &v, const Tensor &r)
  {
    Strain = Strain + v;
    return 0;
  }

const Tensor& PressureDependentElastic3D::getTangentTensor (void)
  {
    ComputeElasticStiffness();
    return Dt;
  }

const stresstensor PressureDependentElastic3D::getStressTensor (void)
  {
    //First time? then return zero stress
    if ( Dt.rank() == 0 )
      return Stress;
    else
      {
        Stress = Dt("ijkl") * Strain("kl");
        return Stress;
      }
  }

const straintensor PressureDependentElastic3D::getStrainTensor (void)
  {
    return Strain;
  }

const straintensor PressureDependentElastic3D::getPlasticStrainTensor (void)
  {
    //Return zero straintensor
    straintensor t;
    return t;
  }

int PressureDependentElastic3D::commitState (void)
  {
    //Set the new Elastic constants
    tensor ret( 4, def_dim_4, 0.0 );

    // Kronecker delta tensor
    tensor I2("I", 2, def_dim_2);

    tensor I_ijkl = I2("ij")*I2("kl");


    //I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;

    //Update E according to hydrostatic stress
    Stress = this->getStressTensor();
    //Dt("ijkl") * Strain("kl");
    double p = Stress.p_hydrostatic();
    //opserr << " p = " <<  p;

    //Cut-off pressure
    if (p <= p_cutoff)
      p = p_cutoff;

    double Ec = E * pow(p/p_ref, exp);
    //opserr << " Eo = " << E << " Ec = " << Ec << " Exp:" << exp<< " p_ref:" << p_ref << " po: " << po<< endln;
    //opserr << " coef = " << Ec/E << endln;

    // Building elasticity tensor
    ret = I_ijkl*( Ec*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( Ec / (1.0 + v) );

    //ret.print();
    Dt = ret;

    return 0;
  }

int PressureDependentElastic3D::revertToLastCommit (void)
  {
    // Nothing was previously committed
    return 0;
  }

int PressureDependentElastic3D::revertToStart (void)
  {
    // Nothing was previously committed
    return 0;
  }

NDMaterial* PressureDependentElastic3D::getCopy (void)
  {
    PressureDependentElastic3D *theCopy =
    new PressureDependentElastic3D (this->getTag(), E, v, rho, exp, p_ref, p_cutoff);
    //opserr << "In Get copy" <<  *theCopy << endln;
    theCopy->epsilon = this->epsilon;
    theCopy->sigma = this->sigma;
    theCopy->Strain = this->Strain;
    theCopy->Stress = this->Stress;
    // D and Dt are created in the constructor call
    return theCopy;
  }

const char* PressureDependentElastic3D::getType (void) const
  {
    return "ThreeDimensional";
  }

int PressureDependentElastic3D::getOrder (void) const
  {
    return 6;
  }

int PressureDependentElastic3D::sendSelf(int commitTag, Channel &theChannel)
  {
    int res = 0;

    static Vector data(6);

    data(0) = this->getTag();
    data(1) = E;
    data(2) = v;
    data(3) = exp;
    data(4) = p_ref;
    data(5) = p_cutoff;

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
    if (res < 0)
      {
        opserr << "PressureDependentElastic3D::sendSelf -- could not send Vector\n";
        return res;
      }

    return res;
  }

int PressureDependentElastic3D::recvSelf(int commitTag,
                                         Channel &theChannel,
                                         FEM_ObjectBroker &theBroker)
  {
    int res = 0;

    static Vector data(6);

    res += theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0)
      {
        opserr << "PressureDependentElastic3D::recvSelf -- could not recv Vector\n";
        return res;
      }

    this->setTag((int)data(0));
    E = data(1);
    v = data(2);
    exp = data(3);
    p_ref = data(4);
    p_cutoff = data(5);

    // Set up the elastic constant matrix for 3D elastic isotropic
    this->ComputeElasticStiffness();

    return res;
  }

void PressureDependentElastic3D::Print(OPS_Stream &s, int flag)
  {
    s << "PressureDependentElastic3D" << endln;
    s << "\ttag: " << this->getTag() << endln;
    s << "\tE: " << E << endln;
    s << "\tv: " << v << endln;
    s << "\texp: " << exp << endln;
    s << "\tp_ref: " << p_ref << endln;
    s << "\tp_cutoff: " << p_cutoff << endln;
    //s << "\tD: " << D << endln;
  }


//================================================================================
void PressureDependentElastic3D::ComputeElasticStiffness(void)
  {
    tensor ret( 4, def_dim_4, 0.0 );
    // Kronecker delta tensor
    tensor I2("I", 2, def_dim_2);
    tensor I_ijkl = I2("ij")*I2("kl");
    //I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
    //Initialize E according to initial pressure in the gauss point
    Stress = getStressTensor();
    //Dt("ijkl") * Strain("kl");
    double p = Stress.p_hydrostatic();
    if (p <= p_cutoff)
      p = p_cutoff;
    double Eo = E * pow(p/p_ref, exp);
    //opserr << " p_ref = " <<  p_ref << " p = " << p << endln;
    //opserr << " coef = " << Eo/E << endln;
    //opserr << " E@ref = " << E << " Eo = " << Eo << endln;

    // Building elasticity tensor
    ret = I_ijkl*( Eo*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( Eo / (1.0 + v) );
    //ret.print();
    Dt = ret;
    return;
  }
