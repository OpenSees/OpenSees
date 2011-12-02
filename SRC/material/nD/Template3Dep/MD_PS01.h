//
//================================================================================
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Manzari-Dafalias potential criterion 01 (with Pc)         #
//# CLASS:             MDPotentialSurface01                                      #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
//#                                                                              #
//#                                                                              #
//# DATE:              August 08 '00                                             #
//# UPDATE HISTORY:    December 13 '00                                           #
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

#ifndef MD_PS01_H    
#define MD_PS01_H

#include "EPState.h"
#include "PS.h"
#include <BJtensor.h>


class MDPotentialSurface01 : public PotentialSurface
{
  private:
    double Pc;

  public:
    PotentialSurface *newObj();  //create a colne of itself
    MDPotentialSurface01(double pc);        // Default constructor

    tensor dQods(const EPState *EPS) const; 
    tensor d2Qods2(const EPState *EPS) const ;   

    //aux. functions for d2Qods2
    tensor dnods(const EPState *EPS) const;
    tensor dthetaoverds(const EPState *EPS) const;
    double dgoverdt(double theta, double c) const;
    tensor apqdnods(const EPState *EPS) const;
    
    void print() { opserr << *this; };

    //================================================================================
    // Overloaded Insertion Operator
    // prints an PotentialSurface's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const MDPotentialSurface01 &PS);

};

#endif

