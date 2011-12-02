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
#ifndef DP_PF_H
#define DP_PF_H

#include "PlasticFlow.h"
#include <math.h>
#include <iostream.h>

class DP_PF : public PlasticFlow
{
  public:   
    DP_PF(int dilatant_which_in = -1, int index_dilatant_in = 0, 
          int alpha_which_in = -1, int index_alpha_in = 0);
    
    ~DP_PF();     
    
    PlasticFlow* newObj();
   
    const straintensor& PlasticFlowTensor(const stresstensor &Stre, 
                                          const straintensor &Stra, 
                                          const MaterialParameter &MaterialParameter_in) const;

  private:
     
    double getdilatant(const MaterialParameter &MaterialParameter_in) const;
    const stresstensor& getalpha(const MaterialParameter &MaterialParameter_in) const; 

  private:
    
    int dilatant_which;
    int index_dilatant;
    int alpha_which;
    int index_alpha;    
    
    static straintensor DPm;
};


#endif

