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
//# UPDATE HISTORY:
//#
//===============================================================================

#ifndef fdFlowDP_CPP
#define fdFlowDP_CPP

#include "fdFlowDP.h"

//--------------------------------------------------------------------
fdFlowDP::fdFlowDP(double dilation_in, double k_in)
:dilation(dilation_in), k(k_in)
{

}

//--------------------------------------------------------------------
fdFlow * fdFlowDP::newObj()
{
     fdFlow *newfdyd = new fdFlowDP(dilation, k);
     
     return newfdyd;
}

//-------------------------------------------------------------------
// Q = dilation*I1 + sqrt(0.5*Sij*Sji) = 0, 
// Note here NumRank = 1, No Kinematic Hardening
//-------------------------------------------------------------------

//--------------------------------------------------------------------
stresstensor fdFlowDP::dFods(const stresstensor &sts, const FDEPState &fdepstate) const
{    
    // NumRank=1, No Ki Hardening
    
    double y = 0.0;

    tensor tI2("I", 2 , def_dim_2);
    stresstensor dev = sts.deviator();
    tensor st = dev("ij")*dev("ij");  
      st.null_indices();
    double x = st.trace();

    if (fabs(x) > 1.0e-6)
      y = 1.0/sqrt(2.0*x);

    return tI2 *dilation + dev *y;
}

//--------------------------------------------------------------------
double fdFlowDP::dFodq(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    //NumRank=1, No Ki Hardening
    
    return -1.0;
}

////--------------------------------------------------------------------
//stresstensor fdFlowDP::dFoda(const stresstensor &sts, const FDEPState &fdepstate) const
//{
//}

//--------------------------------------------------------------------
tensor fdFlowDP::d2Fodsds(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    // NumRank=1, No Ki Hardeing
    
    tensor z4(4, def_dim_4, 0.0);
    
    tensor tI2("I", 2, def_dim_2);
    stresstensor dev = sts.deviator();
    tensor st = dev("ij")*dev("ij")*2.0;  
      st.null_indices();
    double x = st.trace();
    tensor I4 = tI2("ij")*tI2("kl"); 
      I4.null_indices();
    I4 = I4 - I4*(1.0/3.0);  
    tensor st4 = dev("ij")*dev("kl");
      st4.null_indices();

    if (fabs(x) > 1.0e-6)
      z4 = I4 *(1.0/sqrt(x)) - st4 *(2.0/x/sqrt(x));
    
    return z4;   
}

////--------------------------------------------------------------------
//stresstensor fdFlowDP::d2Fodsdq(const stresstensor &sts, const FDEPState &fdepstate ) const
//{
//}

////--------------------------------------------------------------------
//tensor fdFlowDP::d2Fodsda(const stresstensor &sts, const FDEPState &fdepstate ) const
//{
//}
    
////--------------------------------------------------------------------
//double fdFlowDP::d2Fodqdq(const stresstensor &sts, const FDEPState &fdepstate ) const
//{
//}

////--------------------------------------------------------------------
//stresstensor fdFlowDP::d2Fodqda(const stresstensor &sts, const FDEPState &fdepstate ) const
//{
//}

//tensor fdFlowDP::d2Fodada(const stresstensor &sts, const FDEPState &fdepstate ) const
//{
//}

//--------------------------------------------------------------------
OPS_Stream& operator<<(OPS_Stream& os, const fdFlowDP &fdfdDP)
{
    os << "fdFlowDP Parameters: " << "\n";
    os << "DilatedAngl: " << fdfdDP.dilation << "\n";
    os << "k: " << fdfdDP.k<< "\n";
    
    return os;
}


#endif

