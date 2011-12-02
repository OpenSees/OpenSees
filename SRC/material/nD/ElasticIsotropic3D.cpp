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
//# PURPOSE:           Elastic Isotropic Material implementation:
//# CLASS:             ElasticIsotropic3D
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhaohui Yang, Boris Jeremic
//#
//#
//# DATE:              10Oct2000
//# UPDATE HISTORY:    22Nov2002 small fixes 
//#                    Aug2006   Z.Cheng
//#
//===============================================================================
                                                                        
#include <ElasticIsotropic3D.h>


Matrix ElasticIsotropic3D::D(6,6);	  // global for ElasticIsotropic3D only
Vector ElasticIsotropic3D::sigma(6);	 // global for ElasticIsotropic3D onyl

Tensor ElasticIsotropic3D::Dt(4, def_dim_4, 0.0);
stresstensor ElasticIsotropic3D::Stress;

ElasticIsotropic3D::ElasticIsotropic3D
(int tag, double E, double nu, double rho):
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropic3D, E, nu, rho),
 epsilon(6) 
{
	// Set up the elastic constant matrix for 3D elastic isotropic 
	D.Zero();
}

ElasticIsotropic3D::ElasticIsotropic3D():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropic3D, 0.0, 0.0, 0.0),
 epsilon(6)
{

}

ElasticIsotropic3D::~ElasticIsotropic3D ()
{

}

int
ElasticIsotropic3D::setTrialStrain (const Vector &v)
{
	epsilon = v;

	return 0;
}

int
ElasticIsotropic3D::setTrialStrain (const Vector &v, const Vector &r)
{
	epsilon = v;

	return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Vector &v)
{
	epsilon += v;

	return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Vector &v, const Vector &r)
{
	epsilon += v;

	return 0;
}

const Matrix&
ElasticIsotropic3D::getTangent (void)
{
   double mu2 = E/(1.0+v);
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
ElasticIsotropic3D::getStress (void)
{
       double mu2 = E/(1.0+v);
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
}

const Vector&
ElasticIsotropic3D::getStrain (void)
{
  return epsilon;
}

int
ElasticIsotropic3D::setTrialStrain (const Tensor &v)
{
    Strain = v;
    return 0;
}

int
ElasticIsotropic3D::setTrialStrain (const Tensor &v, const Tensor &r)
{
    Strain = v;
    return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Tensor &v)
{
    //opserr << " before set Tri St Incr " << Strain;
    //opserr << " Strain Incr " << v << endlnn;
    Strain = Strain + v;
    //opserr << " after setTrialStrainIncr  " << Strain << endlnn;
    return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    Strain = Strain + v;
    return 0;
}

const Tensor&
ElasticIsotropic3D::getTangentTensor (void)
{
  setInitElasticStiffness();  
  return Dt;
}

const stresstensor&
ElasticIsotropic3D::getStressTensor (void)
{
  setInitElasticStiffness();
  Stress = Dt("ijkl") * Strain("kl");
  return Stress;
}

const straintensor&
ElasticIsotropic3D::getStrainTensor (void)
{
    return Strain;
}

int
ElasticIsotropic3D::commitState (void)
{
    return 0;
}

int
ElasticIsotropic3D::revertToLastCommit (void)
{
	return 0;
}

int
ElasticIsotropic3D::revertToStart (void)
{
	return 0;
}

NDMaterial*
ElasticIsotropic3D::getCopy (void)
{
	ElasticIsotropic3D *theCopy =
		new ElasticIsotropic3D (this->getTag(), E, v, rho);
	theCopy->epsilon = this->epsilon;
	theCopy->Strain = this->Strain;

	return theCopy;
}

const char*
ElasticIsotropic3D::getType (void) const
{
	return "ThreeDimensional";
}

int
ElasticIsotropic3D::getOrder (void) const
{
	return 6;
}

void
ElasticIsotropic3D::Print(OPS_Stream &s, int flag)
{
	s << "ElasticIsotropic3D" << endln;
	s << "\ttag: " << this->getTag() << endln;
	s << "\tE: " << E << endln;
	s << "\tv: " << v << endln;
	s << "\trho: " << rho << endln;
}


//================================================================================
void ElasticIsotropic3D::setInitElasticStiffness(void)
{        				       
    // Kronecker delta tensor
    tensor I2("I", 2, def_dim_2);

    tensor I_ijkl = I2("ij")*I2("kl");

    //I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;
    
    // Building elasticity tensor
    Dt = I_ijkl*( E*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( E / (1.0 + v) );

    return;

}


