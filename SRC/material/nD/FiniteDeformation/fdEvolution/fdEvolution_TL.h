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

#ifndef fdEvolution_TL_H
#define fdEvolution_TL_H

#include "fdEvolution_T.h"

class fdEvolution_TL : public fdEvolution_T
{
 private:
    double H_linear;

 public:
    
    fdEvolution_TL(double H_linear_in = 0.0);
    fdEvolution_TL(const stresstensor &sts, const FDEPState &fdepstate);
    //virtual ~fdEvolution_TL() {};

    fdEvolution_T *newObj();
    
    // Derivative of stress like hardening variable to strain like variable
    tensor HModulus(const stresstensor &sts, const FDEPState &fdepstate) const;  


    void print();

    friend OPS_Stream& operator<< (OPS_Stream& os, const fdEvolution_TL & fdetl);

};


#endif
