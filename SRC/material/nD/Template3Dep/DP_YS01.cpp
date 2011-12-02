
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Drucker-Prager yield criterion 01 (with Pc)	         #
//#                      (Ref. Geotechnique                    			 #
//#                      V.47 No.2 255-272, 1997)                                #
//# CLASS:             DPYieldSurface01                                          #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 03 '00                                             #
//# UPDATE HISTORY:    December 13, '00                                             #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//================================================================================


#ifndef DP_YS01_CPP
#define DP_YS01_CPP

#include "DP_YS01.h"
#include <basics.h>

//================================================================================
// Normal constructor
//================================================================================

DPYieldSurface01::DPYieldSurface01(double pc)  
{ 
     Pc = pc;
} 


//================================================================================
//create a colne of itself
//================================================================================

YieldSurface * DPYieldSurface01::newObj() {  

     YieldSurface  *new_YS = new DPYieldSurface01( Pc );
     return new_YS;

}

//================================================================================
// Copy constrstructor
//================================================================================
//
//DPYieldSurface01::DPYieldSurface01(DPYieldSurface01 &MDYS ) { }
//


//================================================================================
//  Yield criterion evaluation function F(EPState)
//================================================================================

double DPYieldSurface01::f(const EPState *EPS) const
{
  stresstensor S = EPS->getStress().deviator();
  //cout << "s " << S;
  //S.print("sigma", "S");

  double p =EPS->getStress().p_hydrostatic();
  p = p - Pc;
  opserr << "p " << p;

  stresstensor alpha = EPS->getTensorVar( 1 );
  //alpha.print("alpha", " alpha");
  //opserr << "alpha " << alpha;
  
  double m = EPS->getScalarVar(1);

  //stresstensor temp1 = S - p*alpha;
  stresstensor temp1 = S - p*alpha; // S +p * alpha
  temp1.null_indices();
  //temp1.print("temp1 ", " temp1");

  stresstensor temp2 = temp1("ij") * temp1("ij");  

  double temp3 = sqrt( temp2.trace() );

  temp3 = temp3 - sqrt(2.0/3.0) * m * p;        
  //printf("\n========Inside f  temp3 = %.4f x = %.4f\n ", temp3, x);

  return temp3;
}


//================================================================================
// tensor dF/dsigma_ij
//================================================================================

tensor DPYieldSurface01::dFods(const EPState *EPS) const {
  
  tensor dFoverds( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  stresstensor S = EPS->getStress().deviator();
  //S.reportshort("S");

  double p = EPS->getStress().p_hydrostatic();
  p = p - Pc;
  //printf("Here we go!  p %f\n", p);
  	    
  stresstensor alpha = EPS->getTensorVar( 1 ); // getting alpha_ij from EPState
  //alpha.reportshort("alpha");
  
  //stresstensor n = EPS->getTensorVar( 3 );     // getting n_ij from EPState
  //n.reportshort("n");
  
  //-------------------------------------------------
  // might be moved to Evolution Law
    stresstensor r = S * (1.0 / p);
    //r.reportshort("r");
    stresstensor r_bar = r - alpha;
    //r_bar.reportshort("r_bar"); 
    stresstensor norm2 = r_bar("ij") * r_bar("ij");
    double norm = sqrt( norm2.trace() );
    
    //opserr << "d_macheps " << d_macheps() << endlnn;

    stresstensor n;
    if ( norm >= d_macheps() ){ 
      n =  r_bar*(1.0 / norm );
    }
    else {
      opserr << "DPYieldSurface01::dFods  |n_ij| = 0, divide by zero! Program exits.\n";
      exit(-1);
    }
    //EPS->setTensorVar( 3, n); //update n_ij//
  //-------------------------------------------------


  double m = EPS->getScalarVar( 1 );

  
  //tensorial multiplication produces 1st-order tensor
  //tensor temp = n("ij") * n("ij");
  //double temp1 = temp.trace();
  //printf("==== n_ij*n_ij %e\n", temp1);

  //!!Very important:  N = n_pq * alpha_pq +sqrt(2/3)*m (always) = n_pq*r_pq(not hold when not on the yield surface)
  tensor temp = n("ij") * alpha("ij");
  double N = temp.trace() + sqrt(2.0/3.0)*m; 
  //printf("    N =  %e\n", N);

  dFoverds =  n - N *I2 *(1.0/3.0);

  return dFoverds;

}


//================================================================================
// double xi_s1 = dF/dm = -(2/3)^0.5 p  Derivative in terms of first scalar var 
//================================================================================

double DPYieldSurface01::xi_s1( const EPState *EPS ) const {

    double p = EPS->getStress().p_hydrostatic();
    //double p = EPS->getStress().Iinvariant1()/3.0;
    p = p - Pc;
    return -1.0*sqrt(2.0/3.0) * p;

}

//================================================================================
// double xi_t1 = dF/dt1 = dF/dalpha= -p*n  Derivative in terms of first tensorial var 
//================================================================================

tensor DPYieldSurface01::xi_t1( const EPState *EPS) const {
    tensor dFoverds( 2, def_dim_2, 0.0);
    tensor I2("I", 2, def_dim_2);

    stresstensor S = EPS->getStress().deviator();
    
    double p = EPS->getStress().p_hydrostatic();
    p = p - Pc;
    	    
    stresstensor alpha = EPS->getTensorVar( 1 ); // getting alpha_ij from EPState
    
    stresstensor r = S * (1.0 / p); //for p = sig_kk/3
    stresstensor r_bar = r - alpha;
    stresstensor norm2 = r_bar("ij") * r_bar("ij");
    double norm = sqrt( norm2.trace() );
    
    stresstensor n;
    if ( norm >= d_macheps() ){ 
      n = r_bar *(1.0 / norm );
    }
    else {
      opserr << "DPYieldSurface01::dFods  |n_ij| = 0, divide by zero! Program exits.\n";
      exit(-1);
    }
      
    return (-1.0) * n * p;
}


OPS_Stream& operator<< (OPS_Stream& os, const DPYieldSurface01 & YS)
{
   os << "Drucker-Prager Yield Surface 01 Parameters: " << endln;
   os << "Pc = " << YS.Pc << endln;
   return os;
}

#endif

