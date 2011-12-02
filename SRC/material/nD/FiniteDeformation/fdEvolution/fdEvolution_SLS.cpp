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

#ifndef fdEvolution_SLS_CPP
#define fdEvolution_SLS_CPP

#include "fdEvolution_SLS.h"

// Linear and / or Saturated isotropic hardening

//------------------------------------------------------------------------
fdEvolution_SLS::fdEvolution_SLS(double H_linear_in,
		                 double q_saturated_in,
		                 double delta_in)
:H_linear(H_linear_in), q_saturated(q_saturated_in), delta(delta_in)
{
    
}

//------------------------------------------------------------------------
fdEvolution_S * fdEvolution_SLS::newObj() 
{   
    fdEvolution_S *newEL = new fdEvolution_SLS( *this );    
    return newEL;
}

//------------------------------------------------------------------------
double fdEvolution_SLS::HModulus(const stresstensor &sts, const FDEPState &fdepstate) const 
{
    double xi = fdepstate.getStrainLikeInVar();
    return H_linear + delta*q_saturated*exp(-delta*xi);
}

//------------------------------------------------------------------------
void fdEvolution_SLS::print()
{
    opserr << (*this);
}

//------------------------------------------------------------------------
OPS_Stream& operator<< (OPS_Stream& os, const fdEvolution_SLS & fdesl)
{
   os.precision(5);
   os.width(10);
   os << "Linear & Saturate Scalar Linear Evolution Law's Modulus: " << "\n";
   os << "H_linear = " << fdesl.H_linear << "; " << "\n";
   os << "H_initial = " << fdesl.q_saturated << "; " << "\n";
   os << "Delta= " << fdesl.delta << "; " << "\n";
   
   return os;
}

#endif
