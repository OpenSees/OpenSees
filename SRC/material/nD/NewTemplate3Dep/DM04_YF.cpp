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
	// f = [(sij-p*aij)(sij-p*aij)]^0.5 - sqrt(2/3)*p*m = 0

	double sqrt23 = sqrt(2.0/3.0);
	double p = Stre.p_hydrostatic();
	double m = getm(MaterialParameter_in);
	stresstensor alpha = getalpha(MaterialParameter_in);
	stresstensor s_bar = Stre.deviator() - (alpha *p);
	double temp1 = ( s_bar("ij") * s_bar("ij") ).trace();
	
	return sqrt(temp1) - sqrt23*m*p; // a little moving of cone point
}

//================================================================================
const stresstensor& DM04_YF::StressDerivative(const stresstensor& Stre, 
                                              const MaterialParameter &MaterialParameter_in) const
{
//    //using f = 0
//    BJtensor KroneckerI("I", 2, def_dim_2);
//    stresstensor r;
//    stresstensor n;
//    double nr = 0.0;
//    double sqrt23 = sqrt(2.0/3.0);
//    double p = Stre.p_hydrostatic();
//    double m = getm(MaterialParameter_in);
//    stresstensor alpha = getalpha(MaterialParameter_in);
//    if (p != 0.0) {
//	    r = Stre.deviator() *(1.0/p);
//	    n = (r - alpha) *(1.0/(sqrt23*m));
//	    nr = ( n("ij") * r("ij") ).trace();
//	}
//	
//    DM04st = n + ( KroneckerI *(nr/3.0) );
//    return DM04st;

 	// no using f = 0
     BJtensor KroneckerI("I", 2, def_dim_2);
     stresstensor n;
     double n_alpha = 0.0;
     double sqrt23 = sqrt(2.0/3.0);
     double p = Stre.p_hydrostatic();
     double m = getm(MaterialParameter_in);
     stresstensor alpha = getalpha(MaterialParameter_in);
     stresstensor s_bar = Stre.deviator() - (alpha *p);
     double _s_bar_ = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
     if (p > 0.0 && _s_bar_ > 0.0) {
 	    n = s_bar * (1.0/_s_bar_);
 	    n_alpha = ( n("ij") * alpha("ij") ).trace();
     }
     
     DM04st = n + ( KroneckerI * ( (n_alpha + sqrt23*m)/3.0 ) );
     
     return DM04st;

}

//================================================================================
const stresstensor& DM04_YF::InTensorDerivative(const stresstensor& Stre, 
                                                const MaterialParameter &MaterialParameter_in, 
                                                int which) const
{
	// using f =0
    //if (which == index_alpha) {
	//	double sqrt23 = sqrt(2.0/3.0);
	//	double m = getm(MaterialParameter_in);
	//	double p = Stre.p_hydrostatic();
	//	stresstensor alpha = getalpha(MaterialParameter_in);
	//	stresstensor s_bar = Stre.deviator() - (alpha *p);
	//	DM04st = s_bar *(-1.0/(sqrt23*m));
	//}
	//else {
	//	cout << "DM04_YF: Invalid Input. " << endl;
	//	exit (1);
	//}
    //
	//return DM04st;

	// no using f =0
    if (which == index_alpha) {
	    stresstensor n;
	    double p = Stre.p_hydrostatic();
	    stresstensor alpha = getalpha(MaterialParameter_in);
        stresstensor s_bar = Stre.deviator() - (alpha *p);
        double _s_bar_ = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
	    if (p > 0.0 && _s_bar_ > 0.0)
		    n = s_bar * (1.0/_s_bar_);
		DM04_YF::DM04st = n *(-p);
	}
	else {
		cout << "DM04_YF: Invalid Input. " << endl;
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
		    cout << "DM04_YF: Invalid Input, m <= 0.0. " << endl;
		    exit (1);
	    }
		return m;
    }
	else {
		cout << "Warning!! DM04_YF: Invalid Input. " << endl;
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
		cout << "Warning!! DM04_YF: Invalid Input. " << endl;
		exit (1);
	}
}

#endif

