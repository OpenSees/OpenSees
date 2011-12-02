///*
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           CAM CLAY yield criterion                                  #
//# CLASS:             CAMYieldSurface01(for monotonic loading)                  #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, ZHaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              Mar. 28, 2001                                             #
//# UPDATE HISTORY:                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# SHORT EXPLANATION:                                                           #
//#                                                                              #
//#                                                                              #
//================================================================================
//*/

#ifndef CAM_YS_H
#define CAM_YS_H

#include <stresst.h>
#include "EPState.h"
#include "YS.h"
#include <BJtensor.h>


class CAMYieldSurface : public YieldSurface
{
  // Private vars to define the Mazari-Dafalias Yield Surface
  private:
    double M;	// the slope of critical state line

  public:
    YieldSurface *newObj();                  //create a colne of itself
    CAMYieldSurface(double Mp = 1.2);                          // Default constructor

    double f(const EPState *EPS) const;
    tensor dFods(const EPState *EPS) const;

    // Redefine 1st derivative of F over scalar internal variables
    double xi_s1( const EPState *EPS ) const;  // df/dm

    // Redefine 1st derivative of F over tensorial internal variables

    double getM() const;
    void print() { opserr << *this; };
    
// moved to stresstensor
//   private:
//     tensor dpoverds( ) const;
//     tensor dqoverds(const EPState *EPS) const;
//     tensor dthetaoverds(const EPState *EPS) const;
		         
    //================================================================================
    // Overloaded Insertion Operator
    friend OPS_Stream& operator<< (OPS_Stream& os, const CAMYieldSurface & YS);

};

#endif

