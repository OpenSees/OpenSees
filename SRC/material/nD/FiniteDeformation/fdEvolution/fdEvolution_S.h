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

#ifndef fdEvolution_S_H
#define fdEvolution_S_H 

#include <stresst.h>
#include <straint.h>
#include <math.h>

#include <FDEPState.h>

class fdEvolution_S
{
 public:
    
    fdEvolution_S() {};
    virtual ~fdEvolution_S() {};

    virtual fdEvolution_S *newObj();
    
    // Derivative of stress like hardening variable to strain like variable
    virtual double HModulus(const stresstensor &sts, const FDEPState &fdepstate) const;  

    void print();

    friend OPS_Stream& operator<< (OPS_Stream& os, const fdEvolution_S & ev);

};


#endif
