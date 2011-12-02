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

#ifndef DP_YF_CPP
#define DP_YF_CPP

#include "DP_YF.h"

stresstensor DP_YF::DPst;

//================================================================================
DP_YF::DP_YF(int a_which_in, int index_a_in, 
             int k_which_in, int index_k_in, 
             int alpha_which_in, int index_alpha_in)
: a_which(a_which_in), index_a(index_a_in), 
  k_which(k_which_in), index_k(index_k_in),
  alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
DP_YF::~DP_YF() 
{

}

//================================================================================
YieldFunction* DP_YF::newObj() 
{
	YieldFunction  *new_YF = new DP_YF(a_which, index_a, 
	                                   k_which, index_k, 
	                                   alpha_which, index_alpha);

	return new_YF;
}

//================================================================================
double DP_YF::YieldFunctionValue( const stresstensor& Stre, 
                                 const MaterialParameter &MaterialParameter_in ) const
{
	// f = a*I1 + [0.5(sij-p*aij)(sij-p*aij)]^0.5 - k = 0
    if (a_which == -1) {
		cout << "DP_YF: Invalid Input Parameter. " << endl;
		exit(1);
	}
	if (alpha_which == -1)
		return Stre.Iinvariant1() * geta(MaterialParameter_in) 
		     + sqrt(Stre.Jinvariant2()) 
		     - getk(MaterialParameter_in);
	else {
		double temp1 = Stre.Iinvariant1() * geta(MaterialParameter_in) - getk(MaterialParameter_in);
		double p = Stre.p_hydrostatic();
		stresstensor s_back = getalpha(MaterialParameter_in);
		stresstensor s_bar = Stre.deviator() - (s_back * p);
		double temp2 = ( s_bar("ij") * s_bar("ij") ).trace();
		return temp1 + sqrt(0.5*temp2);
	}
}

//================================================================================
const stresstensor& DP_YF::StressDerivative(const stresstensor& Stre, 
                                            const MaterialParameter &MaterialParameter_in) const
{
	BJtensor KroneckerI("I", 2, def_dim_2);
	if (alpha_which == -1) {
		double temp0 = Stre.Jinvariant2();
        temp0 = sqrt(0.5*temp0);
		if (fabs(temp0) < 1.0e-7) {
			DPst = KroneckerI *geta(MaterialParameter_in);
			return DPst;
		}
		DPst = KroneckerI *geta(MaterialParameter_in) + ( Stre.deviator() *(0.5/temp0) );
		return DPst;
	}
	else {
		double p = Stre.p_hydrostatic();
		stresstensor s_back = getalpha(MaterialParameter_in);
		stresstensor s_bar = Stre.deviator() - (s_back * p);
		double temp1 = ( s_bar("ij") * s_bar("ij") ).trace();
        temp1 = sqrt(0.5*temp1);
		if (fabs(temp1) < 1.0e-7) {
			DPst = KroneckerI *geta(MaterialParameter_in);
			return DPst;
		}
		stresstensor I_back = KroneckerI - s_back;
		double temp2 = ( s_bar("ij") * I_back("ij") ).trace();
		DPst = s_bar + ( KroneckerI * (temp2/3.0) );
		DPst = DPst*(0.5/temp1);
        DPst += ( KroneckerI *geta(MaterialParameter_in) );
		return DPst;
	}
}

//================================================================================
double DP_YF::InScalarDerivative(const stresstensor& Stre, 
                                 const MaterialParameter &MaterialParameter_in, 
                                 int which) const
{
	if (a_which == 1 && which == index_a)
		return Stre.Iinvariant1();
	
	if (k_which == 1 && which == index_k)
		return -1.0;
		
	return 0.0;
}

//================================================================================
const stresstensor& DP_YF::InTensorDerivative(const stresstensor& Stre, 
                                              const MaterialParameter &MaterialParameter_in, 
                                              int which) const
{
	if (alpha_which != 2 || which != 1){
		cout << "DP_YF: Invalid Input Parameter. " << endl;
		exit (1);
	}
	
	double p = Stre.p_hydrostatic();
	stresstensor s_back = getalpha(MaterialParameter_in);
	stresstensor s_bar = Stre.deviator() - (s_back * p);
	double temp1 = ( s_bar("ij") * s_bar("ij") ).trace();
	temp1 = sqrt(0.5*temp1);
	if (temp1 < 1.0e-7) {
		DPst = DPst*0.0;
		return DPst;
	}
	DPst = s_bar *(-0.5*p/temp1);
	return DPst;
}

//================================================================================
int DP_YF::getNumInternalScalar() const
{
	int Numyf = 0;
	
	if ( a_which == 1)
		Numyf++;
	
	if ( k_which == 1)
		Numyf++;
	
	return Numyf;
}

//================================================================================
int DP_YF::getNumInternalTensor() const
{
	if (alpha_which == 2)
		return 1;
	else
		return 0;
}

//================================================================================   
int DP_YF::getYieldFunctionRank() const
{
	return 1;
}

//================================================================================   
double DP_YF::geta(const MaterialParameter &MaterialParameter_in) const
{
	// to get a
	if ( a_which == 0) {
		if ( index_a <= MaterialParameter_in.getNum_Material_Parameter() && index_a > 0)
			return MaterialParameter_in.getMaterial_Parameter(index_a-1); 
		else {
			cout << "DP_YF: Invalid Input. " << endl;
			exit (1);
		}
	}
	else if ( a_which == 1) {
		if ( index_a <= MaterialParameter_in.getNum_Internal_Scalar() && index_a > 0)
			return MaterialParameter_in.getInternal_Scalar(index_a-1); 
		else {
			cout << "DP_YF: Invalid Input. " << endl;
			exit (1);
		}
    }
	else {
		cout << "DP_YF: Invalid Input. " << endl;
		exit(1);
	}
}

//================================================================================   
double DP_YF::getk(const MaterialParameter &MaterialParameter_in) const
{
	// to get k
	if ( k_which == 0) {
		if ( index_k <= MaterialParameter_in.getNum_Material_Parameter() && index_k > 0)
			return MaterialParameter_in.getMaterial_Parameter(index_k-1); 
		else {
			cout << "DP_YF: Invalid Input. " << endl;
			exit (1);
		}
	}
	else if ( k_which == 1) {
		if ( index_k <= MaterialParameter_in.getNum_Internal_Scalar() && index_k > 0)
			return MaterialParameter_in.getInternal_Scalar(index_k-1); 
		else {
			cout << "DP_YF: Invalid Input. " << endl;
			exit (1);
		}
    }
	else {
		cout << "DP_YF: Invalid Input. " << endl;
		exit(1);
	}
}

//================================================================================ 
const stresstensor& DP_YF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	//to get alpha
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
		DPst = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return DPst;
	}
	else {
		cout << "DP_YF: Invalid Input. " << endl;
		exit (1);
	}
}


#endif

