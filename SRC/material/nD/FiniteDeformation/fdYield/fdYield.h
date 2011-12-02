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

#ifndef fdYield_H
#define fdYield_H

#include <stresst.h>
#include <straint.h>

#include <FDEPState.h>

class fdYield
{
  public:
    fdYield();
    virtual ~fdYield() {}; 
    
    virtual fdYield *newObj() = 0;   

    virtual int getNumRank();
    virtual double getTolerance();
    virtual double Yd(const stresstensor &sts, const FDEPState &fdepstate ) const = 0;	
    
    virtual stresstensor dYods(const stresstensor &sts, const FDEPState &fdepstate ) const = 0; 
    
    virtual double dYodq(const stresstensor &sts, const FDEPState &fdepstate ) const;	 
    virtual stresstensor dYoda(const stresstensor &sts, const FDEPState &fdepstate ) const; 
    
    virtual void print() = 0; 

    friend OPS_Stream& operator<< (OPS_Stream& os, const fdYield & fdyd);
};


#endif
