///*
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Mazari - Dafalias  yield criterion                        #
//# CLASS:             MDYieldSurface                                            #
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
//#                                                                              #
//================================================================================
//*/

#ifndef MD_YS_H
#define MD_YS_H

#include <stresst.h>
#include "EPState.h"
#include "YS.h"
#include <BJtensor.h>


class MDYieldSurface : public YieldSurface
{
  private:		  // Private vars to define the Mazari-Dafalias Yield Surface

  public:
    YieldSurface *newObj();                  //create a colne of itself
    MDYieldSurface();                          // Default constructor
    //MDYieldSurface(const MDYieldSurface &);  // Default constructor

    double f(const EPState *EPS) const;
    tensor dFods(const EPState *EPS) const;

    // Redefine 1st derivative of F over scalar internal variables
    double xi_s1( const EPState *EPS ) const;  // df/dm
    //double xi_s2( const EPState *EPS ) const;

    // Redefine 1st derivative of F over tensorial internal variables
    tensor xi_t1(const EPState *EPS) const; // dF / d alpha_ij

    void print() { opserr << *this; };
  
    //================================================================================
    // Overloaded Insertion Operator
    friend OPS_Stream& operator<< (OPS_Stream& os, const MDYieldSurface & YS);

};

#endif

