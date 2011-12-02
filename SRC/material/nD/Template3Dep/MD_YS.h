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

#include "EPState.h"
#include "YS.h"


class MDYieldSurface : public YieldSurface
{
  private:		  // Private vars to define the Mazari-Dafalias Yield Surface

  public:
    YieldSurface *newObj();                  //create a colne of itself
    MDYieldSurface();                          // Default constructor
    //MDYieldSurface(const MDYieldSurface &);  // Default constructor

    double f(const EPState *EPS) const;
    tensor dFods(const EPState *EPS) const;
    void print() { cout << *this; }; //pure virtual func
  
    //================================================================================
    // Overloaded Insertion Operator
    // prints an YieldSurface's contents 
    //================================================================================
    friend ostream& operator<< (ostream& os, const MDYieldSurface & YS)
    {
       os << "Manzari-Dafalias Yield Surface Parameters: " << endln;
       return os;
    }

};

#endif

