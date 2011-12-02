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
//# CLASS:                                                                       #
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
//# SHORT EXPLANATION: Functions for rounded Mohr-Coulomb potential function     #
//#                                                                              #
//================================================================================

#ifndef RMC01_PS_H
#define RMC01_PS_H

#include "RMC01.h"
#include <stresst.h>
#include <BJtensor.h>
#include "EPState.h"
#include "PS.h"


class RMC01PotentialSurface : public PotentialSurface
{
  // Private vars to define the RMC01 Potential Surface
  private:
  		  
  
  public:
    RMC01PotentialSurface( ){ };   // Default constructor
    ~RMC01PotentialSurface() { }; //Not all virtual functions  redefined
    PotentialSurface *newObj(); //create a colne of itself

    tensor dQods(const EPState *EPS  ) const;
    tensor d2Qods2(const EPState *EPS) const;
    void print() { opserr << *this; };

    //================================================================================
    // Overloaded Insertion Operator
    // prints an RMC01-PotentialSurface's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const RMC01PotentialSurface &PS);

};

#endif

