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

Matrix PressureDependentElastic3D::D(6,6);
Vector PressureDependentElastic3D::sigma(6);

Tensor PressureDependentElastic3D::Dt(4, def_dim_4, 0.0 );
stresstensor PressureDependentElastic3D::Stress;

PressureDependentElastic3D::PressureDependentElastic3D
(int tag, double E, double nu, double rhop, double expp, double pr, double pop):
 ElasticIsotropicMaterial(tag, ND_TAG_PressureDependentElastic3D, E, nu, rhop),
 exp0(expp), p_ref(pr), p_cutoff(pop), epsilon(6), p_n(0.0), p_n1(0.0)
{

}

PressureDependentElastic3D::PressureDependentElastic3D():
 ElasticIsotropicMaterial(0, ND_TAG_PressureDependentElastic3D, 0.0, 0.0, 0.0),
 exp0(0.0), p_ref(1.0), p_cutoff(1.0), epsilon(6), p_n(0.0), p_n1(0.0)
{

}

PressureDependentElastic3D::~PressureDependentElastic3D ()
{

}

int
PressureDependentElastic3D::setTrialStrain(const Vector &v)
{
  epsilon = v;
  
  return 0;
}

int
PressureDependentElastic3D::setTrialStrain(const Vector &v, const Vector &r)
{
  epsilon = v;
  
  return 0;
}

int
PressureDependentElastic3D::setTrialStrainIncr(const Vector &v)
{
  epsilon += v;
  
  return 0;
}

int
PressureDependentElastic3D::setTrialStrainIncr(const Vector &v, const Vector &r)
{
  epsilon += v;
  
  return 0;
}

const Matrix&
PressureDependentElastic3D::getTangent (void)
{
  double p = p_n;
  if (p <= p_cutoff)
    p = p_cutoff;
  double Eo = E * pow(p/p_ref, exp0);

  double mu2 = Eo/(1.0+v);
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

const Matrix&
PressureDependentElastic3D::getInitialTangent (void)
{
  double Eo = E;

  double mu2 = Eo/(1.0+v);
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

const Vector&
PressureDependentElastic3D::getStress (void)
{
  double p = p_n;
  if (p <= p_cutoff)
    p = p_cutoff;
  double Eo = E * pow(p/p_ref, exp0);

  double mu2 = Eo/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu  = 0.50*mu2;
  
  mu2 += lam;
  
  double eps0 = epsilon(0);
  double eps1 = epsilon(1);
  double eps2 = epsilon(2);
  
  sigma(0) = mu2*eps0 + lam*(eps1+eps2);
  sigma(1) = mu2*eps1 + lam*(eps2+eps0);
  sigma(2) = mu2*eps2 + lam*(eps0+eps1);
  
  p_n1 = (sigma(0)+sigma(1)+sigma(2))/3.0;

  sigma(3) = mu*epsilon(3);
  sigma(4) = mu*epsilon(4);
  sigma(5) = mu*epsilon(5);
  
  return sigma;
}

const Vector&
PressureDependentElastic3D::getStrain (void)
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
    Strain += v;
    return 0;
}

int PressureDependentElastic3D::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    Strain += v;
    return 0;
}

const Tensor& PressureDependentElastic3D::getTangentTensor (void)
{
    return  ComputeElasticStiffness();
}

const stresstensor& PressureDependentElastic3D::getStressTensor (void)
{
    tensor Dt0 = this->getTangentTensor();
	straintensor Da = Strain - CStrain;
	stresstensor De = Dt0("ijkl") * Da("kl");
		De.null_indices();
	Stress = CStress + De;
	
    return Stress;
}

const straintensor& PressureDependentElastic3D::getStrainTensor (void)
{
    return Strain;
}

int PressureDependentElastic3D::commitState (void)
{
  CStress = Stress;
  CStrain = Strain;

  p_n = p_n1;

  return 0;
}

int PressureDependentElastic3D::revertToLastCommit (void)
{
  Stress = CStress;
  Strain = CStrain;

  return 0;
}

int PressureDependentElastic3D::revertToStart (void)
{
  // added: C.McGann, U.Washington for InitialStateAnalysis
  if (ops_InitialStateAnalysis) {
	// do nothing, keep state variables from last step
  } else {
	// normal call for revertToStart (not initialStateAnalysis)
	Stress = Stress*0.0;
    Strain = Strain*0.0;

    p_n = 0.0;
  }

  return 0;
}

NDMaterial*
PressureDependentElastic3D::getCopy (void)
{
    PressureDependentElastic3D *theCopy =
		new PressureDependentElastic3D (this->getTag(), E, v, rho, exp0, p_ref, p_cutoff);

    theCopy->p_n  = p_n;
    theCopy->p_n1 = p_n1;

    return theCopy;
}

NDMaterial*
PressureDependentElastic3D::getCopy(const char *type)
{
  if (strcmp(type,"ThreeDimensional") == 0)
    return this->getCopy();
  else {
    opserr << "PressureDependentElastic3D::getCopy " << type << " not supported" << endln;
    return 0;
  }
}

const char*
PressureDependentElastic3D::getType(void) const
{
    return "ThreeDimensional";
}

int
PressureDependentElastic3D::getOrder(void) const
{
    return 6;
}


int PressureDependentElastic3D::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(7);

    data(0) = this->getTag();
    data(1) = E;
    data(2) = v;
    data(3) = exp0;
    data(4) = p_ref;
    data(5) = p_cutoff;
    data(6) = p_n;

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

    static Vector data(7);

    res += theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0)
      {
        opserr << "PressureDependentElastic3D::recvSelf -- could not recv Vector\n";
        return res;
      }

    this->setTag((int)data(0));
    E = data(1);
    v = data(2);
    exp0 = data(3);
    p_ref = data(4);
    p_cutoff = data(5);
    p_n = data(6);

    return res;
  }

void PressureDependentElastic3D::Print(OPS_Stream &s, int flag)
{
    s << "PressureDependentElastic3D" << "\n";
    s << "tag: " << this->getTag() << "\n";
    s << "E: " << E << "\n";
    s << "v: " << v << "\n";
    s << "exp: " << exp0 << "\n";
    s << "p_ref: " << p_ref << "\n";
    s << "p_cutoff: " << p_cutoff << "\n";
}


//================================================================================
const Tensor& PressureDependentElastic3D::ComputeElasticStiffness(void)
{
    tensor I21("I", 2, def_dim_2);
    tensor I22("I", 2, def_dim_2);
    tensor I_ijkl = I21("ij")*I22("kl");
      I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
   
	double p = CStress.p_hydrostatic();
    if (p <= p_cutoff)
		p = p_cutoff;
    double Eo = E * pow(p/p_ref, exp0);

    Dt = I_ijkl*( Eo*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( Eo / (1.0 + v) );
    
	return Dt;
}

//================================================================================
int PressureDependentElastic3D::setStressTensor(const Tensor& stressIn)
{
	CStress = stressIn;

	return 0;
}

