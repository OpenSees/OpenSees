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
#ifndef RMC_PF_H
#define RMC_PF_H

#define PLASTICFLOW_TAGS_RMC_PF 2

#include "PlasticFlow.h"

#include <math.h>

class RMC_PF : public PlasticFlow
{
  public:   
    RMC_PF(int dilatant_which_in = -1, int index_dilatant_in = 0, 
          int r_which_in = -1, int index_r_in = 0);
    
    ~RMC_PF();     
    
    PlasticFlow* newObj();
   
    const straintensor& PlasticFlowTensor(const stresstensor &Stre, 
                                          const straintensor &Stra, 
                                          const MaterialParameter &MaterialParameter_in) const;


    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

  private:
     
    double getdilatant(const MaterialParameter &MaterialParameter_in) const;
    double getr(const MaterialParameter &MaterialParameter_in) const; 
    
    double RoundedFunctionf1(double s, double r) const;
    double RoundedFunctionf2(double s, double r) const;
    double RoundedFunctiondf1(double s, double r) const;
    double RoundedFunctiondf2(double s, double r) const;    
    
 private:
    
    int dilatant_which;
    int index_dilatant;
    int r_which;
    int index_r;    
    
    static straintensor RMCm;
};


#endif

