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

#ifndef DM04_alpha_Eij_H
#define DM04_alpha_Eij_H 

#define TENSOR_EVOLUTION_TAGS_DM04_alpha_Eij 3

#include "TensorEvolution.h"

class DM04_alpha_Eij : public TensorEvolution
{
  public:
  
    DM04_alpha_Eij(int e0_index_in =0,
                   int e_r_index_in =0,
                   int lambda_c_index_in =0,
                   int xi_index_in =0,
                   int Pat_index_in =0,
                   int m_index_in =0,
                   int M_cal_index_in =0,
                   int cc_index_in =0,
                   int nb_index_in =0,
                   int h0_index_in =0,
                   int ch_index_in =0,
                   int G0_index_in =0,
                   int alpha_index_in =0,
                   int z_index_in =0);

    TensorEvolution* newObj();

    const straintensor& Hij(const straintensor& plastic_flow, const stresstensor& Stre, 
                            const straintensor& Stra, const MaterialParameter& material_parameter);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
              
  private:
    
    double gete0(const MaterialParameter& material_parameter) const; 
    double gete_r(const MaterialParameter& material_parameter) const;      
    double getlambda_c(const MaterialParameter& material_parameter) const;
    double getxi(const MaterialParameter& material_parameter) const;
    double getPat(const MaterialParameter& material_parameter) const;
    double getm(const MaterialParameter& material_parameter) const;
    double getM_cal(const MaterialParameter& material_parameter) const;
    double getcc(const MaterialParameter& material_parameter) const;
    double getnb(const MaterialParameter& material_parameter) const;
    double geth0(const MaterialParameter& material_parameter) const;
    double getch(const MaterialParameter& material_parameter) const;
    double getG0(const MaterialParameter& material_parameter) const;
        
    const stresstensor& getalpha(const MaterialParameter& material_parameter) const;
    const stresstensor& getz(const MaterialParameter& material_parameter) const;
    
    inline double getParameters(const MaterialParameter& material_parameter, int which_in) const;    
    inline double getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const;
    inline double getg(double c, double cos3theta) const;
    
  private:
  
    // things need to be memorized
    int a_index;
    stresstensor alpha_in;

  private:
  
    int e0_index;
    int e_r_index;
    int lambda_c_index;
    int xi_index;
    int Pat_index;
    int m_index;
    int M_cal_index;
    int cc_index;
    int nb_index;
    int h0_index;    
    int ch_index;
    int G0_index;
    int alpha_index;
    int z_index;   
    
    static stresstensor DM04_alpha_t;

};


#endif

