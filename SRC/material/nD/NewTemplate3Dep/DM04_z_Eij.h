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

#ifndef DM04_z_Eij_H
#define DM04_z_Eij_H 

#include "TensorEvolution.h"
#include <iostream.h>

class DM04_z_Eij : public TensorEvolution
{
  public:
  
    DM04_z_Eij(int m_index_in,
               int c_z_index_in,
               int z_max_index_in,
               int alpha_index_in,
               int z_index_in);

    TensorEvolution* newObj();

    const straintensor& Hij(const straintensor& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter);
       
  private:
  
    double getm(const MaterialParameter& material_parameter) const;
    double getc_z(const MaterialParameter& material_parameter) const;
    double getz_max(const MaterialParameter& material_parameter) const;
        
    const stresstensor& getalpha(const MaterialParameter& material_parameter) const;
    const stresstensor& getz(const MaterialParameter& material_parameter) const;   

  private:
  
    int m_index;
    int c_z_index;
    int z_max_index;
    int alpha_index;
    int z_index;

    static stresstensor DM04_z_t;
};


#endif

