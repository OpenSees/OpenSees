///*
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Von Mises  yield criterion                                #
//# CLASS:             VMYieldSurface                                            #
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

#ifndef VM_YS_CPP
#define VM_YS_CPP

#include "VM_YS.h"


//================================================================================
// Copy constrstructor
//================================================================================
//
//VMYieldSurface::VMYieldSurface(const VMYieldSurface &VMYS ) {
//
//}

//================================================================================
//create a colne of itself
//================================================================================

YieldSurface * VMYieldSurface::newObj() {  

     VMYieldSurface  *new_YS = new VMYieldSurface();
     return new_YS;

}

//================================================================================
//  Yield criterion evaluation function f = 3/2 *Sij_bar Sij_bar- k^2
//  c.f. pp.274 W.F.Chen Plasticity for Structure Engineers
//================================================================================

double VMYieldSurface::f(const EPState *EPS) const {
    //deviatoric stress tensor
    stresstensor sigma = EPS->getStress();
    stresstensor alpha = EPS->getTensorVar(1);

    stresstensor sigma_bar = sigma - alpha;   
    stresstensor s_bar = sigma_bar.deviator();


    double k = EPS->getScalarVar(1);
    double k2 = k * k;
    
    stresstensor temp1 = s_bar("ij") * s_bar("ij");
    double temp = temp1.trace();
    temp = temp * 3.0 / 2.0;

    double f   = temp - k2;

    return f;
}


//================================================================================
// tensor dF/dsigma_ij = 3*( S_ij - alpha_ij )
//================================================================================

tensor VMYieldSurface::dFods(const EPState *EPS) const {

    stresstensor sigma = EPS->getStress();
    stresstensor alpha = EPS->getTensorVar(1);

    stresstensor sigma_bar = sigma - alpha;   
    stresstensor s_bar = sigma_bar.deviator();
    tensor dFods = 3.0 * s_bar;
    
    return dFods;
}

//================================================================================
// double xi1 = dF/dk = 2.0*k  Derivative in terms of first scalar var
//================================================================================

double VMYieldSurface::xi_s1(const EPState *EPS) const {

    double k = EPS->getScalarVar( 1 );
    
    return 2.0 * k;
}


//================================================================================
// tensor xi_k1=dF/d_alpha_ij=-3.0*S_bar_pq  Derivative in terms of 1st tensor var
//================================================================================

tensor VMYieldSurface::xi_t1(const EPState *EPS) const {

    stresstensor sigma = EPS->getStress();
    stresstensor alpha = EPS->getTensorVar(1);

    stresstensor sigma_bar = sigma - alpha;   
    stresstensor s_bar = sigma_bar.deviator();
    tensor xi = -3.0 * s_bar;
    
    return xi;
}



#endif

