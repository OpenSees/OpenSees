//
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Mazari - Dafalias  potential criterion                    #
//# CLASS:             MDPotentialSurface                                        #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 08 '00                                             #
//# UPDATE HISTORY:                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# SHORT EXPLANATION:                                                           #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//================================================================================
//

#ifndef MD_PS_H    
#define MD_PS_H

#include "EPState.h"
#include "PS.h"
#include <BJtensor.h>


class MDPotentialSurface : public PotentialSurface
{

  public:
    PotentialSurface *newObj();  //create a colne of itself
    MDPotentialSurface();          // Default constructor

    tensor dQods(const EPState *EPS) const; 
    tensor d2Qods2(const EPState *EPS) const ;   
    void print() { cout << *this; };

    //================================================================================
    // Overloaded Insertion Operator
    // prints an PotentialSurface's contents 
    //================================================================================
    friend ostream& operator<< (ostream& os, const MDPotentialSurface &PS)
    {
       os << "Manzari-Dafalias Potential Surface Parameters: " << endln;
       return os;
    }
};

#endif

