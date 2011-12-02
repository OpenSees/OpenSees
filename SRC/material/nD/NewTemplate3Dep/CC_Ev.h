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

#ifndef CC_Ev_H
#define CC_Ev_H 

#include "ScalarEvolution.h"
#include "ElasticState.h"
#include <iostream.h>

class CC_Ev : public ScalarEvolution
{
  public:
  
   CC_Ev(int lambda_index_in,
         int kappa_index_in,
         int e0_index_in,
         int p0_index_in);

    ScalarEvolution* newObj();

    double H(const straintensor& plastic_flow, const stresstensor& Stre, 
             const straintensor& Stra, const MaterialParameter& material_parameter);
       
  private:
  
    double getlambda(const MaterialParameter& material_parameter) const;
    double getkappa(const MaterialParameter& material_parameter) const;
    double gete0(const MaterialParameter& material_parameter) const;
    double getp0(const MaterialParameter& material_parameter) const; 
    
  private:
  
    int lambda_index;
    int kappa_index;
    int e0_index;
    int p0_index;
};


#endif

