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

// This is based on "e = e_ref - log(p/p_ref)"
// So K = (1+e)*p/kappa (Wood, 1990, Soil Behavior and Critical State Soil Mechanics)
// Parameters:
// 1: kappa: unliading slope in e-log(p) relation curve
// 2: v:     poissoin's ratio
// 3: Kc:    cut-off bulk modulus, when K<Kc, let K = Kc
// 4: e0:    initial void ratio

#ifndef elnp_Elastic_CPP
#define elnp_Elastic_CPP

#include "elnp_Elastic.h"

//BJtensor elnp_Elastic::StiffnessH(4, def_dim_4, 0.0);

elnp_Elastic::elnp_Elastic(int kappa_in, 
                           int v_in,
                           int K_c_in,
                           int e0_in,
                           const stresstensor& initialStress,
                           const straintensor& initialStrain)
: ElasticState(initialStress, initialStrain),
  kappa_index(kappa_in),
  v_index(v_in),
  K_c_index(K_c_in),
  e0_index(e0_in)
{

}

// Create a new 
ElasticState* elnp_Elastic::newObj() 
{
    ElasticState * Els = new  elnp_Elastic(this->kappa_index, 
                                           this->v_index, 
                                           this->K_c_index,
                                           this->e0_index,
                                           this->Stress,
                                           this->Strain);
     return Els;
}

// Get Stiffness Tensor
const BJtensor& elnp_Elastic::getElasticStiffness(const MaterialParameter &MaterialParameter_in) const
{
    // Kronecker delta tensor
    BJtensor I2("I", 2, def_dim_2);

    BJtensor I_ijkl = I2("ij")*I2("kl");
    I_ijkl.null_indices();
    BJtensor I_ikjl = I_ijkl.transpose0110();
    BJtensor I_iljk = I_ijkl.transpose0111();
    BJtensor I4s = (I_ikjl+I_iljk)*0.5;
    
    double kappa = getkappa(MaterialParameter_in);
    double v = getv(MaterialParameter_in);
    if (v >= 0.5 || v < -1.0) {
      cout << "Warning!! elnp_Elastic: Invalid Possoin's ratio. " << endl;
      exit (1);
    }
    double Kc = getK_c(MaterialParameter_in);
    double e0 = gete0(MaterialParameter_in);
    
    double p = this->getStress().p_hydrostatic();
    double epsilon_v = this->getStrain().Iinvariant1();
    double e = e0 + (1 + e0) *epsilon_v;
    double Kk = (1.0 + e) *p /kappa;
    
    double K = (Kk > Kc) ? Kk : Kc ;
    double G = K *1.5*(1.0-2.0*v)/(1.0+v);
    
    // To avoid numerical problems
    //if (G < (1.0e-3) *K)  G = 1.0e-3 *K;
       
    // Building elasticity tensor
    ElasticState::ElasticStiffness = I_ijkl *(K - 2.0*G/3.0) + I4s *(2.0*G);

    return ElasticState::ElasticStiffness;
}

// Get kappa
double elnp_Elastic::getkappa(const MaterialParameter &MaterialParameter_in) const
{
    if ( kappa_index > MaterialParameter_in.getNum_Material_Parameter() || kappa_index < 2) { 
        cout << "elnp_Elastic: Invalid Input. " << endl;
        exit (1);
    }
    else
        return MaterialParameter_in.getMaterial_Parameter(kappa_index - 1); 
}

// Get v
double elnp_Elastic::getv(const MaterialParameter &MaterialParameter_in) const
{
    if ( v_index > MaterialParameter_in.getNum_Material_Parameter() || v_index < 2) { 
        cout << "elnp_Elastic: Invalid Input. " << endl;
        exit (1);
    }
    else
      return MaterialParameter_in.getMaterial_Parameter(v_index - 1); 
}

// Get K_c
double elnp_Elastic::getK_c(const MaterialParameter &MaterialParameter_in) const
{
    if ( K_c_index > MaterialParameter_in.getNum_Material_Parameter() || K_c_index < 2) { 
        cout << "elnp_Elastic: Invalid Input. " << endl;
        exit (1);
    }
    else
        return MaterialParameter_in.getMaterial_Parameter(K_c_index - 1); 
}

// Get e0
double elnp_Elastic::gete0(const MaterialParameter &MaterialParameter_in) const
{
    if ( e0_index > MaterialParameter_in.getNum_Material_Parameter() || e0_index < 2) { 
        cout << "elnp_Elastic: Invalid Input. " << endl;
        exit (1);
    }
    else
        return MaterialParameter_in.getMaterial_Parameter(e0_index - 1); 
}

#endif
