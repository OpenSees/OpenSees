///*
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Drucker - Prager  yield criterion                         #
//# CLASS:             DPYieldSurface                                            #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 08 '00                                             #
//# UPDATE HISTORY:                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//================================================================================
//*/

#ifndef DP_YS_CPP
#define DP_YS_CPP

#include "DP_YS.h"


//================================================================================
// Normal constructor
//================================================================================

//DPYieldSurface::DPYieldSurface( double a1d = 0.5) : alfa1(a1d) { } 

//================================================================================
// Normal constructor
//================================================================================

DPYieldSurface::~DPYieldSurface( )
{
} 

//================================================================================
// Copy constrstructor
//================================================================================
//
//DPYieldSurface::DPYieldSurface(const DPYieldSurface &DPYS ) {
//
//}

//================================================================================
//create a colne of itself
//================================================================================

//alpha machines complains on this
//DPYieldSurface * DPYieldSurface::newObj() {  

YieldSurface * DPYieldSurface::newObj() {  

     YieldSurface  *new_YS = new DPYieldSurface();
     return new_YS;

}

//================================================================================
//  Yield criterion evaluation function F(EPState)
//================================================================================

double DPYieldSurface::f(const EPState *EPS) const {
    double temp1 = EPS->getStress().Iinvariant1();
    double temp2 = temp1 * EPS->getScalarVar(1);
    double temp3 = EPS->CurrentStress.Jinvariant2();
    double eps = pow( d_macheps(), (1./2.) );
    if ( temp3 < 0.0 && fabs(temp3) < eps )
      {
        temp3 = 0.0;
      }
    double temp4 = sqrt(temp3);
    double k = EPS->getScalarVar(2);
    double f   = temp2 + temp4 - k;

    return f;
}

//================================================================================
// tensor dF/dsigma_ij  ( eq. 5.206 in Chen )
//================================================================================

tensor DPYieldSurface::dFods(const EPState *EPS) const {
    tensor KroneckerI("I", 2, def_dim_2);
    
    // Second invariant of deviatoric stress tensor
    double stressJ2 = EPS->CurrentStress.Jinvariant2();
    
    // Deviatoric stress tensor
    stresstensor s = EPS->CurrentStress.deviator();
    s.null_indices();
    
    // Tensor dF/dsigma_ij  ( eq. 5.206 in Chen )
    tensor dFods = KroneckerI * EPS->getScalarVar(1);
    dFods.null_indices();
    
    if ( stressJ2 > 0.0 )
      {
        dFods = dFods + s * ( 1.0 / (2.0 * sqrt(stressJ2) ) );
      }
    return dFods;
}

//================================================================================
// double xi_s1 = dF/dS1 = dF/dalfa1 = I1  Derivative in terms of first scalar var 
//================================================================================

double DPYieldSurface::xi_s1( const EPState *EPS ) const {

    double temp = EPS->getStress().Iinvariant1();
    return temp;

}

//================================================================================
// double xi_s2 = dF/dS2 = dF/k = -1.0  Derivative in terms of second scalar var 
//================================================================================

double DPYieldSurface::xi_s2( const EPState *EPS ) const {

    return -1.0;

}


//================================================================================
// double xi_t1 = dF/dt1 = dF/dalpha  Derivative in terms of first tensorial var 
// Need more work
//================================================================================

tensor DPYieldSurface::xi_t1( const EPState *EPS) const {
    stresstensor xi;
    return xi;
}

//================================================================================
OPS_Stream& operator<< (OPS_Stream& os, const DPYieldSurface & YS)
{
   os << "Drucker-Prager Yield Surface Parameters: " << endln;
   os.precision(4);
   //os << "alfa1 = " << YS.getalfa1() << endlnn;
   return os;
}

#endif

