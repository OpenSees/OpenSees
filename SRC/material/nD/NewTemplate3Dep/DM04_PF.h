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
#ifndef DM04_PF_H
#define DM04_PF_H

#include "PlasticFlow.h"
#include <math.h>
#include <iostream.h>

class DM04_PF : public PlasticFlow
{
  public:   
    DM04_PF(int e0_which, int index_e0,
            int e_r_which_in, int index_e_r_in,
            int lambda_c_which_in, int index_lambda_c_in,
            int xi_which_in, int index_xi_in,
            int Pat_which_in, int index_Pat_in,
            int m_which_in, int index_m_in,
            int M_cal_which_in, int index_M_cal_in,
            int cc_which_in, int index_cc_in,        
            int A0_which_in, int index_A0_in,
            int nd_which_in, int index_nd_in,
            int alpha_which_in, int index_alpha_in,
            int z_which_in, int index_z_in);
    ~DM04_PF();     
    PlasticFlow* newObj();
   
    const straintensor& PlasticFlowTensor(const stresstensor &Stre, 
                                          const straintensor &Stra, 
                                          const MaterialParameter &MaterialParameter_in) const;

  private: 

    double gete0(const MaterialParameter &MaterialParameter_in) const;
    double gete_r(const MaterialParameter &MaterialParameter_in) const;      
    double getlambda_c(const MaterialParameter &MaterialParameter_in) const;
    double getxi(const MaterialParameter &MaterialParameter_in) const;
    double getPat(const MaterialParameter &MaterialParameter_in) const;
    double getm(const MaterialParameter &MaterialParameter_in) const;
    double getM_cal(const MaterialParameter &MaterialParameter_in) const;
    double getcc(const MaterialParameter &MaterialParameter_in) const;
    double getA0(const MaterialParameter &MaterialParameter_in) const;
    double getnd(const MaterialParameter &MaterialParameter_in) const;
    
    const stresstensor& getalpha(const MaterialParameter &MaterialParameter_in) const;
    const stresstensor& getz(const MaterialParameter &MaterialParameter_in) const;

    inline double getParameters(const MaterialParameter &MaterialParameter_in, int parIndex_in, int which_in) const;    
    inline double getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const;
    inline double getg(double c, double cos3theta) const;

  private:
    
    int e0_which;         int index_e0;
    int e_r_which;        int index_e_r;
    int lambda_c_which;   int index_lambda_c;    
    int xi_which;         int index_xi;
    int Pat_which;        int index_Pat; 
    int m_which;          int index_m;
    int M_cal_which;      int index_M_cal; 
    int cc_which;         int index_cc;
    int A0_which;         int index_A0;   
    int nd_which;         int index_nd;

    int alpha_which;      int index_alpha;
    int z_which;          int index_z;
            
    static straintensor DM04m;
    static stresstensor DM04temp;
};


#endif

