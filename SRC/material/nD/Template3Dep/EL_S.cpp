/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
#                                                                                #
# CLASS:             EvolutionLaw_S(base Evolution Law class for scalar var      #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              09-02-2000                                                  #
# UPDATE HISTORY:    09-12-2000                                                  #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: Here are some initial definitions of some virtual funs,     #
#                    some of which will be redefined in derived classes.         #
#                                                                                #
//================================================================================
*/

#ifndef EL_CPP
#define EL_CPP

#include "EL_S.h"
#include <basics.h>
    

//================================================================================
//  Create a clone of itself 
//================================================================================
EvolutionLaw_S * EvolutionLaw_S::newObj() {
    
    EvolutionLaw_S *newEL = new EvolutionLaw_S( *this );
    
    return newEL;

}


//================================================================================
// Evaluating h_s ( for the evaluation of Kp )
//================================================================================

double EvolutionLaw_S::h_s( EPState *EPS, PotentialSurface *PS){

    return 0.0;

}


//================================================================================
//  Print content of the base class, will be overwritten! 
//================================================================================
void EvolutionLaw_S::print()
{
    cout << (*this);
}




#endif

