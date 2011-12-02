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
//# UPDATE HISTORY:    20Aug2004 ZC added kinematic hardening part               #
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
//  f = a*I1 + 0.5*sqrt(2.0)*||s_ij - p*alpha_ij|| - k =0
//================================================================================

double DPYieldSurface::f(const EPState *EPS) const {
    double temp1 = EPS->getStress().Iinvariant1();
    double temp2 = temp1 * EPS->getScalarVar(1);

    stresstensor alpha;
    stresstensor s_bar;
    stresstensor sigma = EPS->getStress();
    double p = sigma.p_hydrostatic();
    stresstensor sdev = sigma.deviator();
    double halfRt2 = 0.5 * sqrt(2.0);
    int nod = EPS->getNTensorVar();
    if ( nod >=1 )  { //May not have kinematic hardening
      alpha = EPS->getTensorVar(1);   
      s_bar = sdev - (alpha*p);
    }
    else {
	  s_bar = sdev;
    }  
    stresstensor temp3 = s_bar("ij") * s_bar("ij");
    temp3.null_indices();
    double temp4 = temp3.trace();
    temp4 = sqrt(temp4);
    double temp5 = halfRt2 * temp4;  
    
    double k = EPS->getScalarVar(2);
    
    double f   = temp2 + temp5 - k;

    return f;
}

//================================================================================
// tensor dF/dsigma_ij
//================================================================================

tensor DPYieldSurface::dFods(const EPState *EPS) const {
    stresstensor dFods;
    
	tensor KroneckerI("I", 2, def_dim_2);
    double temp1 =  EPS->getScalarVar(1);    
    
    stresstensor alpha;
    stresstensor s_bar;
    stresstensor Tnsr0;
    double temp6 = 0.0;
    stresstensor temp5;
    stresstensor sigma = EPS->getStress();   
    double p = sigma.p_hydrostatic();
    stresstensor sdev = sigma.deviator();
    double halfRt2 = 0.5 * sqrt(2.0);
    int nod = EPS->getNTensorVar();
    if ( nod >=1 )  { //May not have kinematic hardening
      alpha = EPS->getTensorVar(1);  
      s_bar = sdev - (alpha*p);
      temp5 = alpha("ij") * s_bar("ij");
      temp5.null_indices();
      temp6 = temp5.trace();
      Tnsr0 = KroneckerI*(temp6/3.0);
    }
    else {
	  s_bar = sdev;
    }
    stresstensor temp3 = s_bar("ij") * s_bar("ij");
    temp3.null_indices();
    double temp4 = temp3.trace();
    temp4 = sqrt(temp4);
    Tnsr0 += s_bar;
    double eps = pow( d_macheps(), 0.5 );
    if ( fabs(temp4) > eps )  {
      Tnsr0 = Tnsr0 * (halfRt2/temp4);
    }
        	    
    dFods = (KroneckerI*temp1) + Tnsr0; 
        
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
// double xi_t1 = dF/dt1 = dF/dalpha 
//================================================================================

tensor DPYieldSurface::xi_t1( const EPState *EPS) const {
    stresstensor xi;
    
    stresstensor alpha;
    stresstensor sigma = EPS->getStress();
    double p = sigma.p_hydrostatic();
    double halfRt2 = 0.5 * sqrt(2.0);
    int nod = EPS->getNTensorVar();
    if ( nod >=1 )  { //May not have kinematic hardening
      alpha = EPS->getTensorVar(1); 
    }
    stresstensor sigma_bar = sigma - alpha*p;   
    stresstensor s_bar = sigma_bar.deviator();
    stresstensor temp3 = s_bar("ij") * s_bar("ij");
    temp3.null_indices();
    double temp4 = temp3.trace();
    temp4 = sqrt(temp4);
    double eps = pow( d_macheps(), 0.5 );
    if ( fabs(temp4) > eps )  {
      xi = s_bar * (-halfRt2*p/temp4); //Note the Negative 1 here
    }
        
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
