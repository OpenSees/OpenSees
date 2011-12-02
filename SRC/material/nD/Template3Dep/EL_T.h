
///================================================================================
// COPYRIGHT (C):     :-))                                                        |
// PROJECT:           Object Oriented Finite Element Program                      |
// PURPOSE:           General platform for elaso-plastic efficient and easy       |
//                    constitutive model implementation                           |
//                                                                                |
// CLASS:             EvolutionLaw_T (the base class for tensorial evolution law) |
//                                                                                |
// VERSION:                                                                       |
// LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    |
// TARGET OS:         DOS || UNIX || . . .                                        |
// DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 |
// PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 |
//                                                                                |
//                                                                                |
// DATE:              09-02-2000                                                  |
// UPDATE HISTORY:                                                                |
//                                                                                |
//                                                                                |
//                                                                                |
//                                                                                |
// SHORT EXPLANATION: The base class is necessary for we need to use to runtime   |
//                    polymorphism of C++ to achieve the fexibility of the 	  |
//                    platform.                                                   |
//                                                                                |
//=================================================================================
//

#ifndef EL_T_H
#define EL_T_H 

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>

#include "EPState.h"
#include "PS.h"

class EvolutionLaw_T
{
 public:
    
    EvolutionLaw_T() {};			 //Normal Constructor

    virtual EvolutionLaw_T *newObj();  //create a colne of itself

    // Not necessary since the increment of internal var can be evalulated in constitutive driver!
    //virtual void UpdateVar( EPState *EPS, double dlamda ) = 0; // Evolve only one internal variable
    
    virtual void print(); 	//Print the contents of the evolution law
    
    //virtual void InitVars( EPState *EPS) = 0;  // Initializing eo and E for Manzari-Dafalias model, 
    //                                  // other model might not need it!

    //virtual void setInitD( EPState *EPS) = 0;  // Initializing D once current st hits the y.s. for Manzari-Dafalias model , 
    //                                  // other model might not need it
    
    //Why can't i put const before EPState???????????????????????? Aug. 16, 2000
    //virtual double h( EPState *EPS, double d ) = 0;    // Evaluating hardening function h
    virtual tensor h_t( EPState *EPS, PotentialSurface *PS);    // Evaluating hardening function h

    //================================================================================
    // Overloaded Insertion Operator
    // prints an Evolution Law_T's contents 
    //================================================================================
    friend ostream& operator<< (ostream& os, const EvolutionLaw_T & EL)
    {
       os << "Base of Tensorial Evolution Law's Parameters: Nothing" << endln;
       return os;
    }
};


#endif

