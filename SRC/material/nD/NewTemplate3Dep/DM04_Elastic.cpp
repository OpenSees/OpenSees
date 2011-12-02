///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              
// CLASS:             
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, 
// DATE:              Fall 2005
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//

// This is based on G = G0*Pat*[(2.97-e)^2/(1+e)]*(p/Pat)^0.5
// Richart et al 1970, Li & Dafalias 2000, Dafalias & Mazari 2004
// Parameters:
// 1: G0:    reference shear mudulus factor (no unit);
// 2: v:     Poissoin's ratio;
// 3: Pat:   Atmospheric pressure;
// 4: k_c;   cut-off factor, for p < k_c*Pat, let p = k_c*Pat to calculate G;
// 5: e0;    initial void ratio 

#ifndef DM04_Elastic_CPP
#define DM04_Elastic_CPP

#include "DM04_Elastic.h"

//BJtensor DM04_Elastic::StiffnessH(4, def_dim_4, 0.0);

DM04_Elastic::DM04_Elastic(int G0_in, 
                           int v_in,
                           int Pat_in,
                           int k_c_in,
                           int e0_in, 
                           const stresstensor& initialStress, 
                           const straintensor& initialStrain)
: ElasticState(initialStress, initialStrain),
  G0_index(G0_in),
  v_index(v_in),
  Pat_index(Pat_in),
  k_c_index(k_c_in),
  e0_index(e0_in)
   
{

}

// Create a new 
ElasticState* DM04_Elastic::newObj() 
{
	ElasticState * Els = new  DM04_Elastic(this->G0_index, 
	                                       this->v_index,
	                                       this->Pat_index,
	                                       this->k_c_index,
	                                       this->e0_index,
	                                       this->Stress,
	                                       this->Strain);
     return Els;
}

// Get Stiffness Tensor
const BJtensor& DM04_Elastic::getElasticStiffness(const MaterialParameter &MaterialParameter_in) const
{
    // Kronecker delta tensor
    BJtensor I2("I", 2, def_dim_2);

    BJtensor I_ijkl = I2("ij")*I2("kl");
    I_ijkl.null_indices();
    BJtensor I_ikjl = I_ijkl.transpose0110();
    BJtensor I_iljk = I_ijkl.transpose0111();
    BJtensor I4s = (I_ikjl+I_iljk)*0.5;
    
    double G0 = getG0(MaterialParameter_in);
    double v = getv(MaterialParameter_in);
    double Pat = getPat(MaterialParameter_in); 
    double k_c = getk_c(MaterialParameter_in);
    double e0 = gete0(MaterialParameter_in);
    
    double epsilon_v = this->getStrain().Iinvariant1();
    double e = e0 + (1 + e0) *epsilon_v;
    double ef = (2.97-e)*(2.97-e)/(1.0+e);
    double p_cal = this->getStress().p_hydrostatic();
    double p_cut = Pat *k_c;

    double p = (p_cal > p_cut) ? p_cal : p_cut;
            
    double G = G0 *Pat *ef *sqrt(p/Pat);   
    double K = G * (2.0*(1.0+v)/(3.0*(1.0-2.0*v)));
       
    // Building elasticity tensor
    ElasticState::ElasticStiffness = I_ijkl *(K - 2.0*G/3.0) + I4s *(2.0*G);

    return ElasticState::ElasticStiffness;
}

// Get G0
double DM04_Elastic::getG0(const MaterialParameter &MaterialParameter_in) const
{
	if ( G0_index > MaterialParameter_in.getNum_Material_Parameter() || G0_index < 2) { 
		cout << "DM04_Elastic: Invalid Input. " << endl;
		exit (1);
	}
	else
		return MaterialParameter_in.getMaterial_Parameter(G0_index - 1); 
}

// Get v
double DM04_Elastic::getv(const MaterialParameter &MaterialParameter_in) const
{
	double v = 0.0;
	if ( v_index > MaterialParameter_in.getNum_Material_Parameter() || v_index < 2) { 
		cout << "DM04_Elastic: Invalid Input. " << endl;
		exit (1);
	}
	else {
		v = MaterialParameter_in.getMaterial_Parameter(v_index - 1);
		if (v >= 0.5 || v <= -1.0) {
			cout << "Warning!! DM04_Elastic: Invalid possoin's ratio. " << endl;
			exit (1);
		}
		return v;
	}
}

// Get Pat
double DM04_Elastic::getPat(const MaterialParameter &MaterialParameter_in) const
{
	if ( Pat_index > MaterialParameter_in.getNum_Material_Parameter() || Pat_index < 1) { 
		cout << "DM04_Elastic: Invalid Input. " << endl;
		exit (1);
	}
	else
		return MaterialParameter_in.getMaterial_Parameter(Pat_index - 1); 
}

// Get k_cut
double DM04_Elastic::getk_c(const MaterialParameter &MaterialParameter_in) const
{
	if ( k_c_index > MaterialParameter_in.getNum_Material_Parameter() || k_c_index < 1) { 
		cout << "DM04_Elastic: Invalid Input. " << endl;
		exit (1);
	}
	else
		return MaterialParameter_in.getMaterial_Parameter(k_c_index - 1); 
}

// Get e0
double DM04_Elastic::gete0(const MaterialParameter &MaterialParameter_in) const
{
	if ( e0_index > MaterialParameter_in.getNum_Material_Parameter() || e0_index < 2) { 
		cout << "DM04_Elastic: Invalid Input. " << endl;
		exit (1);
	}
	else
		return MaterialParameter_in.getMaterial_Parameter(e0_index - 1); 
}

#endif

