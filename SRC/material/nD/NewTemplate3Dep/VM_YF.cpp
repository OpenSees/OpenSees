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

#ifndef VM_YF_CPP
#define VM_YF_CPP

#include "VM_YF.h"

stresstensor VM_YF::VMst;

//================================================================================
VM_YF::VM_YF(int k_which_in, int index_k_in, int alpha_which_in, int index_alpha_in)
: k_which(k_which_in), index_k(index_k_in), 
  alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
VM_YF::~VM_YF() 
{

}

//================================================================================
YieldFunction* VM_YF::newObj() 
{
	YieldFunction  *new_YF = new VM_YF(k_which, index_k, alpha_which, index_alpha);

	return new_YF;
}

//================================================================================
double VM_YF::YieldFunctionValue( const stresstensor &Stre, 
                                  const MaterialParameter &MaterialParameter_in ) const
{
	// f = 1.5*(sij-aij)*(sij-aij) - k*k = 0
	// or f = 3.0*J2^2 - k*k = 0
	
	if (alpha_which == -1)
		return Stre.Jinvariant2()*3.0 - pow(getk(MaterialParameter_in), 2);
	if (alpha_which == 2) {
		stresstensor s_back = getbackstress(MaterialParameter_in);
		stresstensor s_bar = Stre.deviator() - s_back;
		double temp2 = ( s_bar("ij") * s_bar("ij") ).trace();
		return temp2*1.5 - pow(getk(MaterialParameter_in), 2);
	}
	else {
		cout << "Warning!! VM_YF: Invalid Input Parameter. " << endl;
		exit (1);
	}
}

//================================================================================
const stresstensor& VM_YF::StressDerivative(const stresstensor &Stre, 
                                            const MaterialParameter &MaterialParameter_in) const
{
	if (alpha_which == -1) {
		VMst = Stre.deviator() *3.0;
		return VMst;
	}
	if (alpha_which == 2) {
		stresstensor s_back = getbackstress(MaterialParameter_in);
		stresstensor s_bar = Stre.deviator() - s_back;
		VMst = s_bar *3.0;
		return VMst;
	}
	else {
		cout << "Warning!! VM_YF: Invalid Input Parameter. " << endl;
		exit (1);
	}
}

//================================================================================
double VM_YF::InScalarDerivative(const stresstensor &Stre, 
                                 const MaterialParameter &MaterialParameter_in, 
                                 int index_) const
{
	if (k_which == 1 && index_ == index_k)
		return -2.0*getk(MaterialParameter_in);
	else {
		cout << "Warning!! VM_YF: Invalid Input Parameter. " << endl;
		exit (1);
	}
}

//================================================================================
const stresstensor& VM_YF::InTensorDerivative(const stresstensor &Stre, 
                                              const MaterialParameter &MaterialParameter_in, 
                                              int index_) const
{
	if (alpha_which == 2 || index_ == index_alpha) {
		stresstensor s_back = getbackstress(MaterialParameter_in);
		stresstensor s_bar = Stre.deviator() - s_back;
		VMst = s_bar *(-3.0);
		return VMst;
	}
	else {
		cout << "Warning!! VM_YF: Invalid Input Parameter. " << endl;
		exit (1);
	}
}

//================================================================================
int VM_YF::getNumInternalScalar() const
{
	if ( k_which == 1)
		return 1;
	else
		return 0;
}

//================================================================================
int VM_YF::getNumInternalTensor() const
{
	if (alpha_which == 2)
		return 1;
	else
		return 0;
}

//================================================================================   
int VM_YF::getYieldFunctionRank() const
{
	return 2;
}

//================================================================================   
double VM_YF::getk(const MaterialParameter &MaterialParameter_in) const
{
	// to get k
	if ( k_which == 0 && index_k <= MaterialParameter_in.getNum_Material_Parameter() && index_k > 2)
		return MaterialParameter_in.getMaterial_Parameter(index_k-1);
	else if( k_which == 1 && index_k <= MaterialParameter_in.getNum_Internal_Scalar() && index_k > 0)
		return MaterialParameter_in.getInternal_Scalar(index_k-1);
	else {
		cout << "Warning!! VM_YF: Invalid Input. " << endl;
		exit (1);
	}
}

//================================================================================ 
const stresstensor& VM_YF::getbackstress(const MaterialParameter &MaterialParameter_in) const
{
	//to get backstress
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
		VMst = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return VMst;
	}
	else {
		cout << "Warning!! VM_YF: Invalid Input. " << endl;
		exit (1);
	}
}


#endif

