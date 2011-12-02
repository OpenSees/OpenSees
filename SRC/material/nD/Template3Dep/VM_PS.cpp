///*
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Von Mises  yield criterion                                #
//# CLASS:             VMPotentialSurface                                        #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 31 '00                                             #
//# UPDATE HISTORY:                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//================================================================================
//*/

#ifndef VM_PS_CPP
#define VM_PS_CPP

#include "VM_PS.h"


//================================================================================
// Copy constrstructor
//================================================================================
//
//VMPotentialSurface::VMPotentialSurface(const VMPotentialSurface &VMYS ) {
//
//}

//================================================================================
//create a colne of itself
//================================================================================

PotentialSurface * VMPotentialSurface::newObj() {  

     PotentialSurface  *new_YS = new VMPotentialSurface();
     return new_YS;

}


//================================================================================
// tensor dQ/dsigma_ij = 3*( S_ij )        c.f. pp.274 W.F.Chen
//================================================================================

tensor VMPotentialSurface::dQods(const EPState *EPS) const {
    
    // Deviatoric stress tensor of sigma - alpha
    stresstensor sigma = EPS->getStress();
    int nod = EPS->getNTensorVar();
    stresstensor alpha;

    if ( nod >=1 ) //May not have kinematic hardening
      alpha = EPS->getTensorVar(1);

    stresstensor sigma_bar = sigma - alpha;   
    stresstensor s_bar = sigma_bar.deviator();
    //s_bar.null_indices();
    
    // Tensor dQ/dsigma_ij 
    tensor dQods = s_bar * 3.0;
    dQods.null_indices();
    
    return dQods;
}

//================================================================================
// tensor d2Qods2 = d[ 3*(S_ij - alpha_ij) ] /dsigma_ij 
//================================================================================

tensor VMPotentialSurface::d2Qods2(const EPState *EPS) const {

    tensor I("I", 2, def_dim_2);
    tensor temp1 = I("im") * I("jn");
    temp1.null_indices();

    tensor I2("I", 2, def_dim_2);
    tensor temp2 = I2("mn") * I2("ij") * (1.0/3.0);
    temp2.null_indices();
    
    //// Deviatoric stress tensor
    //stresstensor s = EPS->CurrentStress.deviator();
    //s.null_indices();
    
    tensor d2Qods2 = (temp1 - temp2) * 3.0;
    d2Qods2.null_indices();
    
    return d2Qods2;
}

// For Consistent Algorithm, Z Cheng, Jan 2004
tensor VMPotentialSurface::d2Qodsdt1(const EPState *EPS) const 
{  
  tensor d2Qoverdsdt1(4, def_dim_4, 0.0);
  tensor I2("I", 2, def_dim_2);
  tensor I4 = I2("ij") * I2("mn");
  I2.null_indices();
  d2Qoverdsdt1 = ( I4.transpose0110() + I4.transpose0111() ) * (-1.5);
  return d2Qoverdsdt1;
}

//================================================================================
// friend OPS_Stream functions for output
//================================================================================

OPS_Stream& operator<< (OPS_Stream& os, const VMPotentialSurface & YS)
{
    os << "Von Mises Potential Surface Parameters: " << endln;
    return os;
}

#endif

