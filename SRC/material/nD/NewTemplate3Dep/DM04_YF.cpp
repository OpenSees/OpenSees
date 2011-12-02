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
// PROGRAMMER:        Zhao Cheng
// DATE:              Fall 2005
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//

// Ref: Dafalias and Manzari 2004: J. Eng. Mech. 130(6), pp 622-634
// Parameters:
// 1- m:      f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*p*m = 0;
// 2- alpha:  f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*p*m = 0;

#ifndef DM04_YF_CPP
#define DM04_YF_CPP

#include "DM04_YF.h"

stresstensor DM04_YF::DM04st;

//================================================================================
DM04_YF::DM04_YF(int m_which_in, int index_m_in, 
                 int alpha_which_in, int index_alpha_in)
: m_which(m_which_in), index_m(index_m_in), 
  alpha_which(alpha_which_in), index_alpha(index_alpha_in)
{

}

//================================================================================
DM04_YF::~DM04_YF() 
{

}

//================================================================================
YieldFunction* DM04_YF::newObj() 
{
	YieldFunction  *new_YF = new DM04_YF(m_which, index_m, 
                                         alpha_which, index_alpha);

	return new_YF;
}

//================================================================================
double DM04_YF::YieldFunctionValue( const stresstensor& Stre, 
                                    const MaterialParameter &MaterialParameter_in ) const
{
	//f = [(sij-p*aij)(sij-p*aij)] - (2/3)*(p*m)^2 = 0
	
	double p = Stre.p_hydrostatic();

	double m = getm(MaterialParameter_in);
	stresstensor alpha = getalpha(MaterialParameter_in);
	
    stresstensor s_bar = Stre.deviator() - (alpha * p);
	double temp1 = ( s_bar("ij") * s_bar("ij") ).trace();
	
    return temp1 - (2.0/3.0)*m*m*p*p - 1.0e-6;
}

//================================================================================
const stresstensor& DM04_YF::StressDerivative(const stresstensor& Stre, 
                                              const MaterialParameter &MaterialParameter_in) const
{
     BJtensor KroneckerI("I", 2, def_dim_2);

     double p = Stre.p_hydrostatic();
     
     double m = getm(MaterialParameter_in);
     stresstensor alpha = getalpha(MaterialParameter_in);
     
	 stresstensor s_bar = Stre.deviator() - (alpha * p);

     double s_bar_alpha = (s_bar("ij")*alpha("ij")).trace();
     
     DM04st = s_bar + ( KroneckerI * ( s_bar_alpha/3.0 + (4.0/9.0)*m*m*p ));
     
     return DM04st;
}

//================================================================================
const stresstensor& DM04_YF::InTensorDerivative(const stresstensor& Stre, 
                                                const MaterialParameter &MaterialParameter_in, 
                                                int which) const
{
    if (which == index_alpha) {
    
        stresstensor alpha = getalpha(MaterialParameter_in);

	    double p = Stre.p_hydrostatic();	        
        stresstensor s_bar = Stre.deviator() - (alpha *p);	    
		
		DM04_YF::DM04st = s_bar *(-p);	    
    }
	else {
		opserr << "DM04_YF: Invalid Input. " << endln;
		exit (1);
	}
    
	return DM04_YF::DM04st;
}

//================================================================================
int DM04_YF::getNumInternalScalar() const
{
	return 0;
}

//================================================================================
int DM04_YF::getNumInternalTensor() const
{
	return 1;
}

//================================================================================   
int DM04_YF::getYieldFunctionRank() const
{
	return 1;
}

//================================================================================   
double DM04_YF::getm(const MaterialParameter &MaterialParameter_in) const
{
	// to get m
	double m = 0.0;
    if ( m_which == 0 && index_m <= MaterialParameter_in.getNum_Material_Parameter() && index_m > 0 ) {
        m = MaterialParameter_in.getMaterial_Parameter(index_m-1);
	    if (m <= 0.0) {
		    opserr << "DM04_YF: Invalid Input, m <= 0.0. " << endln;
		    exit (1);
	    }
		return m;
    }
	else {
		opserr << "Warning!! DM04_YF: Invalid Input. " << endln;
		exit (1);
	}
}

//================================================================================ 
const stresstensor& DM04_YF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	//to get alpha
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0 ) {
		DM04_YF::DM04st = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return DM04_YF::DM04st;
	}
	else {
		opserr << "Warning!! DM04_YF: Invalid Input. " << endln;
		exit (1);
	}
}

#endif

