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

#ifndef DM04_Elastic_H
#define DM04_Elastic_H

#include "ElasticState.h"
#include <math.h>
#include <iostream.h>

//stresstensor zerostress;
//straintensor zerostrain;

class DM04_Elastic : public ElasticState
{  
  public:
  
    DM04_Elastic(int G0_in, 
                 int v_in,
                 int Pat_in,
                 int k_c_in,
                 int e0_in,
                 const stresstensor& initialStress = zerostress, 
                 const straintensor& initialStrain = zerostrain);                    
    
    ElasticState* newObj();
    
    const BJtensor& getElasticStiffness (const MaterialParameter &MaterialParameter_in) const;
    
  private:
  
    double getG0(const MaterialParameter &MaterialParameter_in) const;
    double getv(const MaterialParameter &MaterialParameter_in) const;
    double getPat(const MaterialParameter &MaterialParameter_in) const;
    double getk_c(const MaterialParameter &MaterialParameter_in) const;
    double gete0(const MaterialParameter &MaterialParameter_in) const;    
     
  private:
  
    int G0_index;   
    int v_index;
    int Pat_index;
    int k_c_index;
    int e0_index;    
};

#endif

