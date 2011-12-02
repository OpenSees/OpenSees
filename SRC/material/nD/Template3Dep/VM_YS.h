///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Von Mises         yield criterion                         #
//# CLASS:             VMYieldSurface                                            #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 03 '93                                             #
//# UPDATE HISTORY:    August 08 '00                                             #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# SHORT EXPLANATION: Von Mises yield surface                                   #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/

#ifndef VM_YS_H
#define VM_YS_H

#include <stresst.h>
#include <BJtensor.h>

#include "EPState.h"
#include "YS.h"


class VMYieldSurface : public YieldSurface
{
  //Private vars to define the Von Mises Yield Surface
  private:		  
  
  public:
    // Create a colne of itself
    YieldSurface *newObj();  

    // Default constructor
    VMYieldSurface ( ) {}     

    // Copy constructor
    //VMYieldSurface (const VMYieldSurface & );   

    // Evaluation of f
    double f(const EPState *EPS) const;

    //First derivative of F over sigma
    tensor dFods(const EPState *EPS) const;

    // Redefine 1st derivative of F over first scalar internal variable
    double xi_s1(const EPState *EPS) const;   

    // Redefine 1st derivative of F over first tensorial internal variable
    tensor xi_t1(const EPState *EPS) const;

    void print() { cout << *this; }; 
  
    //================================================================================
    // Overloaded Insertion Operator
    // prints an VM YieldSurface's contents 
    //================================================================================
    friend ostream& operator<< (ostream& os, const VMYieldSurface & YS)
    {
       os << "Von Mises Yield Surface Parameters: " << endln;
       return os;
    }

};

#endif

