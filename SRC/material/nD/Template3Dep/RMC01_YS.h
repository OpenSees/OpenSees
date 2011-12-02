//================================================================================
//# COPY LEFT and RIGHT:                                                         #
//# Commercial    use    of    this  program without express permission of the   #
//# University  of  California, is strictly encouraged. Copyright and Copyleft   #
//# are covered by the following clause:                                         #
//#                                                                              #
//# Woody's license:                                                             #
//# ``This    source    code is Copyrighted in U.S., by the The Regents of the   #
//# University  of  California,  for  an indefinite period, and anybody caught   #
//# using  it  without  our  permission,  will be mighty good friends of ourn,   #
//# cause  we  don't give a darn. Hack it. Compile it. Debug it. Run it. Yodel   #
//# it. Enjoy it. We wrote it, that's all we wanted to do.'' bj                  #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           Rounded Mohr Coulomb Potential Surface                    #
//# CLASS:             RMC01YieldSurface                                         #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++                                                       #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic jeremic@ucdavis.edu                         #
//#                    Zhao Cheng,                                               #
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic                                 #
//#                                                                              #
//#                                                                              #
//# DATE:              12 Feb. 2003                                              #
//# UPDATE HISTORY:                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//# SHORT EXPLANATION: Functions for rounded Mohr-Coulomb yield function         #
//#                                                                              #
//================================================================================

#ifndef RMC01_YS_H
#define RMC01_YS_H

#include "RMC01.h"
#include <stresst.h>
#include <BJtensor.h>
#include "EPState.h"
#include "YS.h"


class RMC01YieldSurface : public YieldSurface
{
  // Private vars to define the RMC01 Yield Surface
  private:
    
  
  public:
    YieldSurface *newObj();  //create a clone of itself
    
    RMC01YieldSurface ( ) {}    // Default constructor
    virtual ~RMC01YieldSurface() {}     // Destructor

    double f(const EPState *EPS) const;
    tensor dFods(const EPState *EPS) const;

    // Redefine 1st derivative of F over scalar internal variables
    double xi_s1( const EPState *EPS ) const;
    double xi_s2( const EPState *EPS ) const;

    // Redefine 1st derivative of F over tensorial internal variables
//    tensor xi_t1(const EPState *EPS) const;

    void print() { opserr << *this; }; 
  
    //================================================================================
    // Overloaded Insertion Operator
    // prints an RMC01 YieldSurface's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const RMC01YieldSurface & YS);

};

#endif

