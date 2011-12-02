///*
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           CAM CLAY potential criterion                              #
//# CLASS:             CAMPotentialSurface                                           #
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

#ifndef CAM_PS_H
#define CAM_PS_H

#include <stresst.h>
#include "EPState.h"
#include "PS.h"
#include <BJtensor.h>


class CAMPotentialSurface : public PotentialSurface
{
  // Private vars to define the Mazari-Dafalias Potential Surface
  private:
    double M;	// the slope of critical state line

  public:
    PotentialSurface *newObj();                  //create a colne of itself
    CAMPotentialSurface(double Mp = 1.2);        // Default constructor

    tensor dQods(const EPState *EPS) const;
    tensor d2Qods2(const EPState *EPS) const;

    void print() { opserr << *this; };
    double getM() const;
    

// moved to stresstensor Boris Jeremic@ucdavis.edu 21Aug2001
//  private:
//    tensor dpoverds( ) const;
//    tensor dqoverds(const EPState *EPS) const;
//    tensor dthetaoverds(const EPState *EPS) const;
//    tensor d2poverds2( void ) const;
//    tensor d2qoverds2(const EPState *EPS) const;
//    tensor d2thetaoverds2(const EPState *EPS) const;
//	     	          
    //================================================================================
    // Overloaded Insertion Operator
    friend OPS_Stream& operator<< (OPS_Stream& os, const CAMPotentialSurface & YS);

};

#endif

