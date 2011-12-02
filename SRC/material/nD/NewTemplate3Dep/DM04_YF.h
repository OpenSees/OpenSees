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
#ifndef DM04_YF_H
#define DM04_YF_H

#include "YieldFunction.h"

#define YIELDFUNCTION_TAGS_DM04_YF 5

#include <math.h>

class DM04_YF : public YieldFunction
{
  public:
    DM04_YF(int m_which_in =0, int index_m_in =0, 
            int alpha_which_in =0, int index_alpha_in =0);
    ~DM04_YF();
      
    YieldFunction *newObj();
  
    double YieldFunctionValue(const stresstensor& Stre, 
                              const MaterialParameter &MaterialParameter_in) const;
    
    const stresstensor& StressDerivative(const stresstensor& Stre, 
                                         const MaterialParameter &MaterialParameter_in) const;
       
    const stresstensor& InTensorDerivative(const stresstensor& Stre, 
                                           const MaterialParameter &MaterialParameter_in, 
                                           int which) const;
    
    int getNumInternalScalar() const;
    int getNumInternalTensor() const;
    int getYieldFunctionRank() const;

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  private:
    double getm(const MaterialParameter &MaterialParameter_in) const;
    const stresstensor& getalpha(const MaterialParameter &MaterialParameter_in) const;
    
  private:    
    int m_which; 
    int index_m;
    int alpha_which; 
    int index_alpha;
        
    static stresstensor DM04st;
};


#endif

