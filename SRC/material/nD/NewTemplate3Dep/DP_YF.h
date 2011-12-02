///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              
// CLASS:             
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, 
// DATE:              Fall 2005
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//
#ifndef DP_YF_H
#define DP_YF_H

#include "YieldFunction.h"
#include <math.h>
#include <iostream.h>

class DP_YF : public YieldFunction
{
  public:
    DP_YF(int a_which_in = -1, int index_a_in = 0, 
          int k_which_in = -1, int index_k_in = 0,
          int alpha_which_in = -1, int index_alpha_in = 0);
    ~DP_YF();
      
    YieldFunction *newObj();
  
    double YieldFunctionValue(const stresstensor &Stre, 
                              const MaterialParameter &MaterialParameter_in) const;
           
    const stresstensor& StressDerivative(const stresstensor &Stre, 
                                         const MaterialParameter &MaterialParameter_in) const;
    
    double InScalarDerivative(const stresstensor &Stre, 
                              const MaterialParameter &MaterialParameter_in, 
                              int which) const; 
                         
    const stresstensor& InTensorDerivative(const stresstensor &Stre, 
                                           const MaterialParameter &MaterialParameter_in, 
                                           int which) const;
    
    int getNumInternalScalar() const;
    int getNumInternalTensor() const;
    int getYieldFunctionRank() const;

  private:
    double geta(const MaterialParameter &MaterialParameter_in) const;
    double getk(const MaterialParameter &MaterialParameter_in) const;
    const stresstensor& getalpha(const MaterialParameter &MaterialParameter_in) const;
    
  private:    
    int a_which;
    int index_a;
    int k_which;
    int index_k;
    int alpha_which;
    int index_alpha;
        
    static stresstensor DPst;
};


#endif

