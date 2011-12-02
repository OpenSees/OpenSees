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

#ifndef fdEvolution_T_CPP
#define fdEvolution_T_CPP

#include "fdEvolution_T.h"

fdEvolution_T * fdEvolution_T::newObj() 
{   
    fdEvolution_T *newEL = new fdEvolution_T( *this );    
    return newEL;
}

tensor fdEvolution_T::HModulus(const stresstensor &sts, const FDEPState &fdepstate) const
{
    tensor Z400(4, def_dim_4, 0.0);
    return Z400;
}

void fdEvolution_T::print()
{
    opserr << (*this);
}

OPS_Stream& operator<< (OPS_Stream& os, const fdEvolution_T & ev)
{
   os << "Tensor Evolution Law's Parameters: " << "\n";
   return os;
}

#endif

