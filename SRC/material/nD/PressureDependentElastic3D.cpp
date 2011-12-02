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


Tensor PressureDependentElastic3D::Dt(4, def_dim_4, 0.0 );
stresstensor PressureDependentElastic3D::Stress;

PressureDependentElastic3D::PressureDependentElastic3D
(int tag, double E, double nu, double rhop, double expp, double pr, double pop):
 ElasticIsotropicMaterial (tag, ND_TAG_PressureDependentElastic3D, E, nu, rhop),
 exp0(expp), p_ref(pr), p_cutoff(pop)
{

}

PressureDependentElastic3D::PressureDependentElastic3D():
 ElasticIsotropicMaterial (0, ND_TAG_PressureDependentElastic3D, 0.0, 0.0, 0.0)
{

}

PressureDependentElastic3D::~PressureDependentElastic3D ()
{

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
	Stress = Stress*0.0;
	Strain = Strain*0.0;

    return 0;
}

NDMaterial* PressureDependentElastic3D::getCopy (void)
{
    PressureDependentElastic3D *theCopy =
		new PressureDependentElastic3D (this->getTag(), E, v, rho, exp0, p_ref, p_cutoff);

    return theCopy;
}

const char* PressureDependentElastic3D::getType (void) const
{
    return "ThreeDimensional";
}


int PressureDependentElastic3D::sendSelf(int commitTag, Channel &theChannel)
{
    int res = 0;

    static Vector data(6);

    data(0) = this->getTag();
    data(1) = E;
    data(2) = v;
    data(3) = exp0;
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
    exp0 = data(3);
    p_ref = data(4);
    p_cutoff = data(5);

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

