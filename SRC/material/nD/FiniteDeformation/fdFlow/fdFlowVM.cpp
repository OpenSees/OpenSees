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

#ifndef fdFlowVM_CPP
#define fdFlowVM_CPP

#include "fdFlowVM.h"

//--------------------------------------------------------------------
fdFlowVM::fdFlowVM(double Y0_in) :Y0(Y0_in)
{

}

//--------------------------------------------------------------------
fdFlow * fdFlowVM::newObj()
{
     fdFlow *newfdyd = new fdFlowVM(Y0);
     return newfdyd;
}

//-------------------------------------------------------------------
// Q = 1.5*(S_ij-a_ij)*(S_ij-a_ij)  - (Y0+q)*(Y0+q) = 0, Note here NumRank = 2
//-------------------------------------------------------------------

//--------------------------------------------------------------------
stresstensor fdFlowVM::dFods(const stresstensor &sts, const FDEPState &fdepstate) const
{    
    return sts.deviator() * 3.0;
}

//--------------------------------------------------------------------
double fdFlowVM::dFodq(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    double q = fdepstate.getStressLikeInVar();
    return -2.0 * ( Y0+q );
}

//--------------------------------------------------------------------
stresstensor fdFlowVM::dFoda(const stresstensor &sts, const FDEPState &fdepstate) const
{    
    return sts.deviator() * (-3.0);
}

//--------------------------------------------------------------------
tensor fdFlowVM::d2Fodsds(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    tensor I2("I", 2 , def_dim_2);
    tensor I4 = I2("ij")*I2("kl"); I4.null_indices();
    //I4 = (I4.transpose0110()+I4.transpose0111())*1.5 - I4; //For symmetric tensor
    I4 = I4*3.0 - I4;	//For general tensor
    return I4;
}

//--------------------------------------------------------------------
tensor fdFlowVM::d2Fodsda(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    tensor I2("I", 2 , def_dim_2);
    tensor I4 = I2("ij")*I2("kl"); I4.null_indices();
    //I4 = (I4.transpose0110()+I4.transpose0111())*1.5 - I4; //For symmetric tensor
    I4 = I4*3.0 - I4;	//For general tensor
    return I4 *(-1.0);
}

//--------------------------------------------------------------------
double fdFlowVM::d2Fodqdq(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    return -2.0;
}

//--------------------------------------------------------------------
tensor fdFlowVM::d2Fodada(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    tensor I2("I", 2 , def_dim_2);
    tensor I4 = I2("ij")*I2("kl"); I4.null_indices();
    //I4 = (I4.transpose0110()+I4.transpose0111())*1.5; //For symmetric tensor
    I4 = I4 *3.0;       //For general tensor
    return I4;
}

//--------------------------------------------------------------------
OPS_Stream& operator<<(OPS_Stream& os, const fdFlowVM &fdflVM)
{
    os << "fdFlowVM Parameters: " << "\n";
    os << "Y0: " << fdflVM.Y0 << "\n";
    return os;
}


#endif

