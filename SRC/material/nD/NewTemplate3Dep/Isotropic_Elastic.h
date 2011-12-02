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

#ifndef Isotropic_Elastic_H
#define Isotropic_Elastic_H

#include "ElasticState.h"
#include <iostream.h>

//stresstensor zerostress;
//straintensor zerostrain;

class Isotropic_Elastic : public ElasticState
{  
  public:
                    
    Isotropic_Elastic(int E_in, 
                   int v_in,
                   const stresstensor& initialStress = zerostress, 
                   const straintensor& initialStrain = zerostrain);
    
    ElasticState *newObj();
    
    const BJtensor& getElasticStiffness(const MaterialParameter &MaterialParameter_in) const;
    
  private:
    
    double getE(const MaterialParameter &MaterialParameter_in) const;
    double getv(const MaterialParameter &MaterialParameter_in) const;
  
  private:
    
    int E_index;   
    int v_index;
};


#endif

