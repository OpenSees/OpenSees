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

#ifndef fdFlowDP_H
#define fdFlowDP_H

#include "fdFlow.h"

class fdFlowDP : public fdFlow
{
  private:
    double DilatedAngle;
    double Cohesion;
    int ConeIndex;
    double k1;
    double k2;
  public:
    fdFlowDP(double DilatedAngle_in, double Cohesion_in, int ConeIndex_in = 0);
    // virtual ~fdFlowDP() {}; 
    
    fdFlow *newObj();   

    stresstensor dFods(const stresstensor &sts, const FDEPState &fdepstate ) const ; 
    double dFodq(const stresstensor &sts, const FDEPState &fdepstate ) const;	    
    stresstensor dFoda(const stresstensor &sts, const FDEPState &fdepstate ) const ; 
    
    tensor d2Fodsds(const stresstensor &sts, const FDEPState &fdepstate ) const ;     	 
    stresstensor d2Fodsdq(const stresstensor &sts, const FDEPState &fdepstate ) const ;
    tensor d2Fodsda(const stresstensor &sts, const FDEPState &fdepstate ) const ;
     
    double d2Fodqdq(const stresstensor &sts, const FDEPState &fdepstate ) const; 
    stresstensor d2Fodqda(const stresstensor &sts, const FDEPState &fdepstate ) const ;

    tensor d2Fodada(const stresstensor &sts, const FDEPState &fdepstate ) const ;

    void print() { opserr << *this; };   

    friend OPS_Stream& operator<< (OPS_Stream& os, const fdFlowDP &fdflDP);
};


#endif
