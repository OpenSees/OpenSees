
//================================================================================
// COPYRIGHT (C):     :-))                                                       |
// PROJECT:           Object Oriented Finite Element Program                     |
// PURPOSE:           General platform for elaso-plastic constitutive model      |
//                    implementation                                             |
//                                                                               |
// CLASS:             EvolutionLaw_L_Eeq (linear Evolution law)                  |
//                                                                               |
//                                                                               |
// VERSION:                                                                      |
// LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )   |
// TARGET OS:         DOS || UNIX || . . .                                       |
// DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                |
// PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                |
//                                                                               |
//                                                                               |
// DATE:              09-02-2000                                                 |
// UPDATE HISTORY:                                                               |
//                                                                               |
//                                                                               |
//                                                                               |
//                                                                               |
// SHORT EXPLANATION: This is a linear evolution law for the evoltion of a       |
//                    scalar variable k which depends on plastic equi. strain    |
//                    i.e. dk = a*de_eq_p                                        |
//                                                                               |
//================================================================================

#ifndef EL_LEeq_H
#define EL_LEeq_H

#include <math.h>

#include <iostream.h>
#include <iomanip.h>

#include "EL_S.h"

class EvolutionLaw_L_Eeq : public EvolutionLaw_S
{
  // Private vars to define the evolution law

    double  a;  //coefficient to define the linear hardening rule of a scalar hardening var

  public:
    //EvolutionLaw_L_Eeq( );    // default constructor---no parameters
    
    EvolutionLaw_L_Eeq( double ad = 10.0) : a(ad) {}
                         
    EvolutionLaw_L_Eeq(const EvolutionLaw_L_Eeq &LEL );   // Copy constructor
    
    EvolutionLaw_S *newObj();                     //create a colne of itself
    
    //void InitVars(EPState *EPS);    // Initialize all hardening vars called only once 
    //                                // after material point is formed if necessary.
    
    //void setInitD(EPState  *EPS);   // set initial D once current stress hits the y.s.
    //                                // was primarily for Manzari-Dafalias model

    //double h( EPState *EPS,  double norm_dQods);     // Evaluating hardening function h
    double h_s( EPState *EPS, PotentialSurface *PS);    // Evaluating hardening function h
    
    //void UpdateVar( EPState *EPS, double dlamda );  // Evolve corresponding var linearly using de_eq_p
    //Moved to CDriver.cpp

    void print();

    // some accessor functions
    double geta() const;      // Linear coefficient used to evolve internal var
    void   seta( double ad);

    //================================================================================
    // Overloaded Insertion Operator	  Zhaohui Added Aug. 13, 2000
    // prints Linear EvolutionLaw's contents 
    //================================================================================
    friend ostream& operator<< (ostream& os, const EvolutionLaw_L_Eeq & LEL)
    {
        os.unsetf( ios::scientific );
        os.precision(5);

        os.width(10);       
        os << endln << "Linear Scalar Evolution Law's parameters:" << endln;
	os << "a = " << LEL.geta() << "; " << endln;
               
        return os;
    }  

    
};


#endif




