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

#ifndef VM_PF_CPP
#define VM_PF_CPP

#include "VM_PF.h"

straintensor VM_PF::VMm;
stresstensor VM_PF::VMb;

//================================================================================
VM_PF::VM_PF(int alpha_which_in, int index_alpha_in)
: alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
VM_PF::~VM_PF() 
{  

}

//================================================================================
PlasticFlow* VM_PF::newObj() 
{  
     PlasticFlow  *new_PF = new VM_PF(alpha_which, index_alpha);
     
     return new_PF;
}

//================================================================================
const straintensor& VM_PF::PlasticFlowTensor(const stresstensor &Stre, 
                                             const straintensor &Stra, 
                                             const MaterialParameter &MaterialParameter_in) const
{
	if (alpha_which == -1) {
		VMm = Stre.deviator() *3.0;
		return VMm;
	}
	if (alpha_which == 2) {
		stresstensor s_back = getalpha(MaterialParameter_in);
		stresstensor s_bar = Stre.deviator() - s_back;
		VMm = s_bar *3.0;
		return VMm;
	}
	else {
		cout << "Warning!! VM_PF: Invalid Input Parameter. " << endl;
		exit (1);
	}
}

//================================================================================ 
const stresstensor& VM_PF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	//to get alpha (backstress)
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
		VMb = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return VMb;
	}
	
	cout << "Warning!! VM_PF: Invalid Input. " << endl;
	exit (1);
}

#endif

