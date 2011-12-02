///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Von Mises         potential criterion                     #
//# CLASS:             VMpotentialSurface                                        #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 03 '93                                             #
//# UPDATE HISTORY:    September 01 '00                                          #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# SHORT EXPLANATION: Von Mises potential surface                               #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/

#ifndef VM_PS_H
#define VM_PS_H

#include <stresst.h>
#include <BJtensor.h>

#include "EPState.h"
#include "PS.h"


class VMPotentialSurface : public PotentialSurface
{
  // Private vars to define the Von Mises Potential Surface
  private:		  
  
  public:
    PotentialSurface *newObj();  //create a colne of itself
    VMPotentialSurface ( ) {}      // Default constructor
    //VMPotentialSurface (const VMPotentialSurface & );   // copy constructor

    tensor dQods(const EPState *EPS) const;
    tensor d2Qods2(const EPState *EPS) const;

    void print() { opserr << *this; }; 
  
    //================================================================================
    // Overloaded Insertion Operator
    // prints an VM PotentialSurface's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const VMPotentialSurface & YS);

};

#endif

