/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
#                                                                                #
# CLASS:             EvolutionLaw_T                                                #
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

#ifndef EL_T_CPP
#define EL_T_CPP

#include "EL_T.h"
#include <basics.h>
    

//================================================================================
//  Create a clone of itself 
//================================================================================
EvolutionLaw_T * EvolutionLaw_T::newObj() {
    
    EvolutionLaw_T *newEL = new EvolutionLaw_T( *this );
    
    return newEL;

}

//================================================================================
// Evaluating h_ ( for the evaluation of Kp )
//================================================================================

tensor EvolutionLaw_T::h_t( EPState *EPS, PotentialSurface *PS)
{
    // Return zero valued tensor
    stresstensor temp;
    return temp;
}

//================================================================================
// updating E, e D and m for Manzari-Dafalias model
//================================================================================

int EvolutionLaw_T::updateEeDm( EPState *EPS, double st_vol, double dLamda)
{
   // do nothing
   return 0;
}  



//================================================================================
//  Print vars defined in Linear Evolution Law_T
//================================================================================
void EvolutionLaw_T::print()
{
    opserr << (*this);
}

//================================================================================
// Overloaded Insertion Operator
// prints base Evolution Law_T's contents 
//================================================================================
OPS_Stream& operator<< (OPS_Stream& os, const EvolutionLaw_T & EL)
{
   os << "Base of Tensorial Evolution Law's Parameters: Nothing" << endln;
   return os;
}




#endif

