
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Manzari-Dafalias  potential surface                       #
//#                    (Ref. Geotechnique v.47 No.2 255-272, 1997)               #
//# CLASS:             MDPotentialSurface                                        #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 03 '93                                             #
//# UPDATE HISTORY:    August 08 '00                                             #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//================================================================================


#ifndef MD_PS_CPP
#define MD_PS_CPP

#include "MD_PS.h"
#include <basics.h>


//================================================================================
// Normal constructor
//================================================================================

MDPotentialSurface::MDPotentialSurface( ) { } 


//================================================================================
//create a colne of itself
//================================================================================

PotentialSurface * MDPotentialSurface::newObj() {  

     MDPotentialSurface  *new_PS = new MDPotentialSurface();
     return new_PS;

}

//================================================================================
//  tensor dQ/dsigma_ij: the normal to the potential surface
//================================================================================

tensor MDPotentialSurface::dQods(const EPState *EPS) const { 

  tensor dQoverds( 2, def_dim_2, 0.0);
  tensor I2("I", 2, def_dim_2);

  tensor S = EPS->getStress().deviator();
  double p = EPS->getStress().p_hydrostatic();
  
  tensor alpha = EPS->getTensorVar( 1 ); //the first tensor hardening var is alpha_ij

  //double m = EPS->getScalarVar( 1 );	 //the first scalar hardening var is m
  double D = EPS->getScalarVar( 2 );	 //the second scalar hardening var is D
  //printf("Here we go!  D= %e\n", D);

  tensor r = S *(1.0/p);
  tensor temp_f1 = r - alpha;
  tensor temp_f2 = temp_f1("ij") * temp_f1("ij");
  double temp_f3 = sqrt( temp_f2.trace() );

  stresstensor n;
  if ( temp_f3 >= d_macheps() ){ 
    n = ( r  - alpha ) * (1.0 / temp_f3 );
  }
  else {
    ::printf(" \n\n n_ij not defined!!!! Program exits\n");
    exit(1);
  }
  
  dQoverds =  n + D * I2 *(1.0/3.0);

  return dQoverds;
}


//================================================================================
// tensor d2Q/dsigma_ij2: the second derivatives to the potential surface
//================================================================================

tensor MDPotentialSurface::d2Qods2(const EPState *EPS) const {

  tensor d2Qoverds2( 2, def_dim_2, 0.0); // dummy second derivatives. To be redefined.

  return d2Qoverds2;

}



#endif

