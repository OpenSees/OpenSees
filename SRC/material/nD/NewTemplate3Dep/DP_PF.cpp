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

#ifndef DP_PF_CPP
#define DP_PF_CPP

#include "DP_PF.h"

straintensor DP_PF::DPm;

//================================================================================
DP_PF::DP_PF(int dilatant_which_in, int index_dilatant_in, int alpha_which_in, int index_alpha_in)
: dilatant_which(dilatant_which_in), index_dilatant(index_dilatant_in), 
  alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
DP_PF::~DP_PF() 
{  

}

//================================================================================
PlasticFlow* DP_PF::newObj() 
{  
     PlasticFlow  *new_PF = new DP_PF(dilatant_which, index_dilatant, 
                                      alpha_which, index_alpha);
     
     return new_PF;
}

//================================================================================
const straintensor& DP_PF::PlasticFlowTensor(const stresstensor &Stre, 
                                                    const straintensor &Stra, 
                                                    const MaterialParameter &MaterialParameter_in) const
{
	BJtensor KroneckerI("I", 2, def_dim_2);
	if (alpha_which == -1) {
		double temp0 = Stre.Jinvariant2();
        temp0 = sqrt(0.5*temp0);
		if (fabs(temp0) < 1.0e-7) {
			DPm = KroneckerI *getdilatant(MaterialParameter_in);
			return DPm;
		}
		DPm = KroneckerI *getdilatant(MaterialParameter_in) + ( Stre.deviator() *(0.5/temp0) );
		return DPm;
	}
	else {
		double p = Stre.p_hydrostatic();
		stresstensor s_back = getalpha(MaterialParameter_in);
		stresstensor s_bar = Stre.deviator() - (s_back * p);
		double temp1 = ( s_bar("ij") * s_bar("ij") ).trace();
        temp1 = sqrt(0.5*temp1);
		if (fabs(temp1) < 1.0e-7) {
			DPm = KroneckerI *getdilatant(MaterialParameter_in);
			return DPm;
		}
		stresstensor I_back = KroneckerI - s_back;
		double temp2 = ( s_bar("ij") * I_back("ij") ).trace();
		DPm = s_bar + ( KroneckerI * (temp2/3.0) );
		DPm = DPm*(0.5/temp1);
        DPm += ( KroneckerI *getdilatant(MaterialParameter_in) );
		return DPm;
	}
}

//================================================================================
double DP_PF::getdilatant(const MaterialParameter &MaterialParameter_in) const
{
	// to dilatant
	if ( dilatant_which == 0) {
		if ( index_dilatant <= MaterialParameter_in.getNum_Material_Parameter() && index_dilatant > 0)
			return MaterialParameter_in.getMaterial_Parameter(index_dilatant-1); 
		else {
			cout << "DP_PF: Invalid Input. " << endl;
			exit (1);
		}
	}
	else if ( dilatant_which == 1) {
		if ( index_dilatant <= MaterialParameter_in.getNum_Internal_Scalar() && index_dilatant > 0)
			return MaterialParameter_in.getInternal_Scalar(index_dilatant-1); 
		else {
			cout << "DP_PF: Invalid Input. " << endl;
			exit (1);
		}
    }
	else {
		cout << "DP_PF: Invalid Input. " << endl;
		exit(1);
	}
}

//================================================================================ 
const stresstensor& DP_PF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	//to get alpha(backstress)
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
		DPm = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return DPm;
	}
	else {
		cout << "DP_PF: Invalid Input. " << endl;
		exit (1);
	}
}

#endif

