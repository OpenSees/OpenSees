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

#ifndef fdYieldVM_CPP
#define fdYieldVM_CPP

#include "fdYieldVM.h"

//--------------------------------------------------------------------
fdYieldVM::fdYieldVM(double Y0_in) :Y0(Y0_in)
{

}

//--------------------------------------------------------------------
fdYield * fdYieldVM::newObj()
{
    fdYield *newfdyd = new fdYieldVM(Y0);
    return newfdyd;
}

int fdYieldVM::getNumRank()
{
    //return 2;
    return 1;
}

double fdYieldVM::getTolerance()
{
    double tol = (Y0*(1.0e-8)) * (Y0*(1.0e-8));
    return tol > 1.0e-8?  tol : 1.0e-8;
    //return Y0*1.0e-8 > 1.0e-8? Y0*1.0e-8 : 1.0e-8;
}

//--------------------------------------------------------------------------------------
// Yd =  3.0*(J2) - (Y0+q)*(Y0+q) = 0, Note here NumRank = 2: No Kinematic hardening
// Yd =  |S_ij| - sqrt(2/3)*(Y0+q) = 0, Note here NumRank = 1: No Kinematic hardening
// Yd =  1.5 * (S_ij - a_ij)*(S_ij-a_ij) - (Y0+q)*(Y0+q) = 0, Note here NumRank = 2
//--------------------------------------------------------------------------------------

double fdYieldVM::Yd(const stresstensor &sts, const FDEPState &fdepstate ) const
{
    //// NumRank=2, No Ki Hardeing
    //double J2 = sts.Jinvariant2();
    //double q = fdepstate.getStressLikeInVar();
    //return 3.0*J2 - (Y0+q)*(Y0+q);
    
    //// NumRank=1, No Ki Hardeing
    //return sqrt(2.0*J2) - sqrt(2.0/3.0);
    
    // NumRank=2, With Ki Hardeing
    stresstensor a = fdepstate.getStressLikeKiVar();
    double q = fdepstate.getStressLikeInVar();
    stresstensor dev = sts.deviator() - a;
    tensor st = dev("ij")*dev("ij");   
      st.null_indices();
    double x = st.trace();
    return 1.5*x - (Y0+q)*(Y0+q);
}

//--------------------------------------------------------------------
stresstensor fdYieldVM::dYods(const stresstensor &sts, const FDEPState &fdepstate ) const
{   
    //// NumRank=2, No Ki Hardeing
    //return sts.deviator() * 3.0;

    //// NumRank=1, No Ki Hardeing
    //double J2 = sts.Jinvariant2();
    //return sts.deviator()/(sqrt(8.0*J2);

    // NumRank=2, With Ki Hardeing
    stresstensor a = fdepstate.getStressLikeKiVar();
    stresstensor dev = sts.deviator() - a;
    return dev *3.0;
}					      

//--------------------------------------------------------------------
double fdYieldVM::dYodq(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    //// NumRank=2, No Ki Hardeing
    //double q = fdepstate.getStressLikeInVar();
    //return -2.0 * (Y0+q);
    
    //// NumRank=1, No Ki Hardeing
    //return sqrt(2.0/3.0);

    // NumRank=2, With Ki Hardeing
    double q = fdepstate.getStressLikeInVar();
    return -2.0 * (Y0+q);
}

//--------------------------------------------------------------------
stresstensor fdYieldVM::dYoda(const stresstensor &sts, const FDEPState &fdepstate ) const
{   
    //// NumRank=2, No Ki Hardeing
    //return sts.deviator() * 3.0;

    //// NumRank=1, No Ki Hardeing
    //double J2 = sts.Jinvariant2();
    //return sts.deviator()/(sqrt(8.0*J2);

    // NumRank=2, With Ki Hardeing
    stresstensor a = fdepstate.getStressLikeKiVar();
    stresstensor dev = sts.deviator() - a;
    return dev *(-3.0);
}	

//--------------------------------------------------------------------
OPS_Stream& operator<<(OPS_Stream& os, const fdYieldVM &fdydVM)
{
    os << "fdYieldVM Parameters: " << "\n";
    os << "Y0: " << fdydVM.Y0 << "\n";
    return os;
}


#endif

