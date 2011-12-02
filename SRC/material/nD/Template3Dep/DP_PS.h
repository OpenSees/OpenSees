///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Drucker - Prager  potential criterion                     #
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
//# SHORT EXPLANATION:                                                           #
//#                                                                              #
//# if alfa1#=0 && alfa2#=0 && alfa1#=alfa2 =>                                   #
//#              Drucker-Prager with non-associated flow rule                    #
//#                                                                              #
//# if alfa1#=0 && alfa2#=0 && alfa1==alfa2 =>                                   #
//#              Drucker-Prager with associated flow rule                        #
//#                                                                              #
//# if alfa1==0 && alfa2#=0 =>                                                   #
//#              Von Mises with non-associated Drucker-Prager flow rule          #
//#                                                                              #
//# if alfa1#=0 && alfa2==0 =>                                                   #
//#              Drucker-Prager with non-associated Von Mises flow rule          #
//#                                                                              #
//# if alfa1==0 && alfa2==0 =>                                                   #
//#              Von Mises with associated flow rule                             #
//#                                                                              #
//################################################################################
//*/

#ifndef DP_PS_H
#define DP_PS_H

#include "PS.h"
#include "EPState.h"
#include <BJtensor.h>


class DPPotentialSurface : public PotentialSurface
{
  private:		  // Private vars to define the Drucker-Prager Potential Surface
    double alfa2;	  // potential surface orientation angle  
  
  public:
    DPPotentialSurface( double a2d = 0.0 ) : alfa2(a2d) {}   // Default constructor
    virtual ~DPPotentialSurface() { }; //Virtual functions not all redefined
    DPPotentialSurface( const DPPotentialSurface &DPPS ); //Copy constructor
    PotentialSurface *newObj();              //create a colne of itself
    double getalfa2() const;

    tensor dQods(const EPState *EPS) const;
    tensor d2Qods2(const EPState *EPS) const;
    
    tensor d2Qodsds1(const EPState *EPS) const; // For Consistent Algorithm, Z Cheng, Jan 2004        
    
    void print() { opserr << *this; };

    //================================================================================
    // Overloaded Insertion Operator
    // prints an DP-PotentialSurface's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const DPPotentialSurface &PS);

};

#endif

