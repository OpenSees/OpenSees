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
#ifndef CC_YF_H
#define CC_YF_H

#include "YieldFunction.h"
#include <math.h>
#include <iostream.h>

class CC_YF : public YieldFunction
{
  public:
    CC_YF(int M_which_in = -1, int index_M_in = 0, 
          int p0_which_in = -1, int index_p0_in = 0);
    ~CC_YF();
      
    YieldFunction *newObj();
  
    double YieldFunctionValue(const stresstensor& Stre, 
                              const MaterialParameter &MaterialParameter_in) const;
                              
    const stresstensor& StressDerivative(const stresstensor& Stre, 
                                         const MaterialParameter &MaterialParameter_in) const;
    
    double InScalarDerivative(const stresstensor& Stre, 
                              const MaterialParameter &MaterialParameter_in, 
                              int which) const;
    
    int getNumInternalScalar() const;
    int getNumInternalTensor() const;
    int getYieldFunctionRank() const;

  private:
    double getM(const MaterialParameter &MaterialParameter_in) const;
    double getP0(const MaterialParameter &MaterialParameter_in) const;
        
  private:    
    int M_which;
    int index_M;
    int p0_which;
    int index_p0;
    
  private:
    static stresstensor CCst;
};

#endif

