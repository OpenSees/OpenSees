
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           OpenSees                                                  #
//# PURPOSE:           Triaxial Failure Criterion for Concrete - yield criterion #
//# CLASS:             TriFCYieldSurface                                         #
//#                                                                              #
//# VERSION:           1.0                                                       #
//# LANGUAGE:          C++ (ili tako nesto)                                      #
//# TARGET OS:                                                                   #
// DESIGNER(S):       Boris Jeremic and Zhaohui Yang [jeremic,zhyang]@ucdavis.edu|
// PROGRAMMER(S):     Vlado Vukadin                                              |
//#                                                                              #
//#                                                                              #
//# DATE:             June 01, 2002                                              #
//# UPDATE HISTORY:    bice tako dobr da nece biti potreban update :)            #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# SHORT EXPLANATION:                                                           #
//#                                                                              #
//# Yield surface is based on article by Menetrey, P. and William, K.J.          #
//# published in 1995 in  ACI Structural Journal pp 311-318. Purpose of the     #
//# Yield surface is to model triaxial strenght of concrete                      #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/

#ifndef Tri_a_fail_crit_YS_H
#define Tri_a_fail_crit_YS_H

#include <stresst.h>
#include <BJtensor.h>
#include "EPState.h"
#include "YS.h"

class TriFCYieldSurface : public YieldSurface
{
  // Private vars to define the TriFCYieldSurface Yield Surface
  private:
     double fcomp;
     double ftens;
     double el; 
     double c;
  public:
    YieldSurface *newObj();  //create a clone of itself
    
    TriFCYieldSurface (double fc, double ft, double e, double coh );   // Default constructor

     

    virtual ~TriFCYieldSurface  ( );     // Destructor

    double f(const EPState *EPS) const;
    tensor dFods(const EPState *EPS) const;

    // Redefine 1st derivative of F over scalar internal variables
    //double xi_s1( const EPState *EPS ) const;
    //double xi_s2( const EPState *EPS ) const;

    // Redefine 1st derivative of F over tensorial internal variables
    //tensor xi_t1(const EPState *EPS) const
    //{
    
    //}
    
    double getfcomp() const;
    double getftens() const;
    double getel() const;
    double get_c() const;

    void print() {opserr << *this; }; 
  
    //================================================================================
    // Overloaded Insertion Operator
    // prints an XX YieldSurface's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const TriFCYieldSurface & YS);

};

#endif

