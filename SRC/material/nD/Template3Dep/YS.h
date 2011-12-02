/*
################################################################################
# COPYRIGHT (C):     :-))                                                      #
# PROJECT:           Object Oriented Finite Element Program                    #
# PURPOSE:           General platform for elaso-plastic constitutive model     #
#                    implementation                                            #
# CLASS:             YieldSurface (the base class for all yield surfaces)      #
#                                                                              #
# VERSION:                                                                     #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
# TARGET OS:         DOS || UNIX || . . .                                      #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                               #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
#                                                                              #
#                                                                              #
# DATE:              08-03-2000                                                #
# UPDATE HISTORY:                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
# SHORT EXPLANATION: The goal is to create a platform for efficient and easy   #
#                    implemetation of any elasto-plastic constitutive model!   #
#                                                                              #
################################################################################*/

#ifndef YS_H
#define YS_H

#include <stresst.h>
#include <straint.h>
#include <BJtensor.h>

#include "EPState.h"

class YieldSurface
{
  public:
    YieldSurface() {} ;			 //Normal Constructor
    virtual YieldSurface *newObj() = 0;  //create a colne of itself
  
    virtual double f( const EPState *EPS ) const = 0;	 //pure virtual func
    virtual tensor dFods( const EPState *EPS ) const = 0;  //pure virtual func
    virtual void print() = 0; //pure virtual func

    // 1st derivative of F over scalar internal variables (at most 4 scalar internal vars allowed currently)
    virtual double xi_s1( const EPState *EPS ) const;	 
    virtual double xi_s2( const EPState *EPS ) const;
    virtual double xi_s3( const EPState *EPS ) const;
    virtual double xi_s4( const EPState *EPS ) const;

    // 1st derivative of F over scalar internal variables (at most 4 tensor internal vars allowed currently)
    virtual tensor xi_t1( const EPState *EPS ) const;	 
    virtual tensor xi_t2( const EPState *EPS ) const;
    virtual tensor xi_t3( const EPState *EPS ) const;
    virtual tensor xi_t4( const EPState *EPS ) const;


    //================================================================================
    // Overloaded Insertion Operator
    // prints an YieldSurface's contents 
    //================================================================================
    friend ostream& operator<< (ostream& os, const YieldSurface & YS);
};


#endif

