///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Drucker - Prager  yield criterion                         #
//# CLASS:             DPPotentialSurface                                        #
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
//################################################################################
//*/

#ifndef DP_PS_CPP
#define DP_PS_CPP

#include "DP_PS.h"



//================================================================================
// Normal constructor
//================================================================================

//DPPotentialSurface::DPPotentialSurface( double a2d = 0.3) : alfa2(a2d) 
//{
//
//} 


//================================================================================
// Copy constrstructor  === not necessary for no pointer member
//================================================================================

DPPotentialSurface::DPPotentialSurface(const DPPotentialSurface &DPPS ) {

    alfa2 =  DPPS.getalfa2();

}


//================================================================================
//create a colne of itself
//================================================================================

PotentialSurface * DPPotentialSurface::newObj() {  

     DPPotentialSurface  *new_PS = new DPPotentialSurface( this->getalfa2() );
     return new_PS;

}



//================================================================================
// tensor dQ/dsigma_ij  ( eq. 5.206 in Chen )
//================================================================================

tensor DPPotentialSurface::dQods(const EPState *EPS) const {

    tensor KroneckerI("I", 2, def_dim_2);
    
    // Second invariant of deviatoric stress tensor
    double stressJ2 = EPS->getStress().Jinvariant2();
    // Deviatoric stress tensor
    stresstensor s = EPS->getStress().deviator();
    s.null_indices();
    // Tensor dQ/dsigma_ij  ( eq. 5.207 in Chen )
    tensor dQods = KroneckerI * getalfa2();
    //tensor dQods = KroneckerI * EPS->getScalarVar(1);
    dQods.null_indices();
    if ( stressJ2 > 0.0 )
      {
        dQods = dQods + s*(1.0/(2.0*sqrt(stressJ2)));
      }
    return dQods;
}


//================================================================================
// tensor d2Q/dsigma_ij_2  ( eq. 5.206 in Chen )
//================================================================================

tensor DPPotentialSurface::d2Qods2(const EPState *EPS) const {

    double stressJ2 = EPS->getStress().Jinvariant2();
    double stressJ2sqrt =pow( stressJ2, 0.5);

    tensor I("I", 2, def_dim_2);
    tensor temp1 = I("im") * I("jn");
    temp1.null_indices();

    tensor I2("I", 2, def_dim_2);
    tensor temp2 = I2("mn") * I2("ij") * (1.0/3.0);
    temp2.null_indices();

    tensor temp3 = temp1 - temp2;
    temp3.null_indices();

    stresstensor s = EPS->getStress().deviator();
    s.null_indices();

    tensor temp4 = s("ij") * s("mn");
    temp4.null_indices();

    tensor d2Qods2 = temp3* (1.0/2.0/stressJ2sqrt) - temp4 * (1.0/4.0/stressJ2/stressJ2sqrt);

    return d2Qods2;
}

//================================================================================

double DPPotentialSurface::getalfa2() const {

    return alfa2; 

}


#endif

