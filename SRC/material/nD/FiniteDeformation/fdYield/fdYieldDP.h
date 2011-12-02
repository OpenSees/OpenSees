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

#ifndef fdYieldDP_H
#define fdYieldDP_H

#include "fdYield.h"

class fdYieldDP : public fdYield
{
  private:
    double alpha;
    double k;
  
  public:
    fdYieldDP(double alpha_in, double k_in);
    // virtual ~fdYieldDP() {}; 
    
    fdYield *newObj();   

    int getNumRank();
    double getTolerance();
    
    double Yd(const stresstensor &sts, const FDEPState &fdepstate ) const;	

    stresstensor dYods(const stresstensor &sts, const FDEPState &fdepstate ) const; 
    double dYodq(const stresstensor &sts, const FDEPState &fdepstate ) const;	    
    //stresstensor dYoda(const stresstensor &sts, const FDEPState &fdepstate ) const;

    void print() { opserr << *this; };   

    friend OPS_Stream& operator<< (OPS_Stream& os, const fdYieldDP & fdydDP);
};


#endif
