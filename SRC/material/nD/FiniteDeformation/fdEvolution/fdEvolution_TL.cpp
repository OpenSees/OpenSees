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

#ifndef fdEvolution_TL_CPP
#define fdEvolution_TL_CPP

#include "fdEvolution_TL.h"

//------------------------------------------------------------------------
fdEvolution_TL::fdEvolution_TL(double H_linear_in)
:H_linear(H_linear_in)
{
    
}

//------------------------------------------------------------------------
fdEvolution_T * fdEvolution_TL::newObj() 
{   
    fdEvolution_T *newEL = new fdEvolution_TL( *this );    
    return newEL;
}

//------------------------------------------------------------------------
tensor fdEvolution_TL::HModulus(const stresstensor &sts, const FDEPState &fdepstate) const 
{
    tensor eta = fdepstate.getStrainLikeKiVar();
    tensor I2("I", 2 , def_dim_2);
    tensor I4 = I2("ij")*I2("kl"); I4.null_indices();
    //I4 = (I4.transpose0110()+I4.transpose0111())*0.5; //For symmetric tensor
    I4 = I4.transpose0110();  //For general tensor
    return I4*H_linear;
}

//------------------------------------------------------------------------
void fdEvolution_TL::print()
{
    opserr << (*this);
}

//------------------------------------------------------------------------
OPS_Stream& operator<< (OPS_Stream& os, const fdEvolution_TL & fdetl)
{
   os.precision(5);
   os.width(10);
   os << "Tensor Linear Evolution Law's Modulus: " << "\n";
   
   return os;
}

#endif
