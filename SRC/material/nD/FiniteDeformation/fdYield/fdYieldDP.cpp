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
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              July 2004
//# UPDATE HISTORY:    change the input parameter from friction angle and cohesion
//#                    to alpha and k, and it is now no cone case index.
//#
//===============================================================================

#ifndef fdYieldDP_CPP
#define fdYieldDP_CPP

#include "fdYieldDP.h"

//--------------------------------------------------------------------
fdYieldDP::fdYieldDP(double alpha_in, double k_in)
 : alpha(alpha_in), k(k_in)
{   

}

//--------------------------------------------------------------------
fdYield * fdYieldDP::newObj()
{
    fdYield *newfdyd = new fdYieldDP(alpha, k);
    
    return newfdyd;
}

int fdYieldDP::getNumRank()
{
    return 1;
}

double fdYieldDP::getTolerance()
{
    double tol = 1.0e-6 * k;
    
    return (tol*tol) > 1.0e-7?  (tol*tol) : 1.0e-7;
}

//--------------------------------------------------------------------------------------------
// Yd =  alpha*I1 + sqrt(0.5*Sij*Sij) - (k+q) = 0, 
// Note here NumRank = 1: No Kinematic hardening
//--------------------------------------------------------------------------------------------

double fdYieldDP::Yd(const stresstensor &sts, const FDEPState &fdepstate ) const
{    
    // NumRank=1, No Ki Hardeing

    stresstensor dev = sts.deviator();
    double I1 = sts.Iinvariant1();
    tensor st = dev("ij")*dev("ij");
      st.null_indices();
    double x = st.trace();
    
    return alpha*I1 + sqrt(0.5*x) - k;
}

//--------------------------------------------------------------------
stresstensor fdYieldDP::dYods(const stresstensor &sts, const FDEPState &fdepstate ) const
{   
    // NumRank=1, No Ki Hardening
    
    double y = 0.0;

    tensor tI2("I", 2, def_dim_2);
    stresstensor dev = sts.deviator();
    tensor st = dev("ij")*dev("ij");  
      st.null_indices();
    double x = st.trace();

    if (fabs(x) > 1.0e-6)
      y = 1.0/sqrt(2.0*x);
    
    return tI2 * alpha + dev * y;
}					      

//--------------------------------------------------------------------
double fdYieldDP::dYodq(const stresstensor &sts, const FDEPState &fdepstate ) const
{      
    // NumRank=1, No Ki Hardening
    
    return -1.0;
}

////--------------------------------------------------------------------
//stresstensor fdYieldDP::dYoda(const stresstensor &sts, const FDEPState &fdepstate ) const
//{
//}
	
//--------------------------------------------------------------------
OPS_Stream& operator<<(OPS_Stream& os, const fdYieldDP &fdydDP)
{
    os << "fdYieldDP Parameters: " << "\n";
    os << "alpha: " << fdydDP.alpha << "\n";
    os << "k: " << fdydDP.k << "\n";
    
    return os;
}


#endif

