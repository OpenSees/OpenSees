///*
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Drucker-Prager yield criterion 01(with Pc)                #
//# CLASS:             DPYieldSurface01                                          #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 03 '00                                             #
//# UPDATE HISTORY:    December 13, 00                                           #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# SHORT EXPLANATION:                                                           #
//#                                                                              #
//#                                                                              #
//================================================================================
//*/

#ifndef DP_YS01_H
#define DP_YS01_H

#include <stresst.h>
#include "EPState.h"
#include "YS.h"
#include <BJtensor.h>


class DPYieldSurface01 : public YieldSurface
{
  // Private vars to define the Mazari-Dafalias Yield Surface
  private:
    double Pc;

  public:
    YieldSurface *newObj();                  //create a colne of itself
    DPYieldSurface01(double pc);             // Default constructor
    //DPYieldSurface01(const DPYieldSurface01 &);  // Default constructor

    double f(const EPState *EPS) const;
    tensor dFods(const EPState *EPS) const;

    // Redefine 1st derivative of F over scalar internal variables
    double xi_s1( const EPState *EPS ) const;  // df/dm
    //double xi_s2( const EPState *EPS ) const;

    // Redefine 1st derivative of F over tensorial internal variables
    tensor xi_t1(const EPState *EPS) const; // dF / d alpha_ij

    void print() {opserr << *this; };
  
    //================================================================================
    // Overloaded Insertion Operator
    friend OPS_Stream& operator<< (OPS_Stream& os, const DPYieldSurface01 & YS);

};

#endif

