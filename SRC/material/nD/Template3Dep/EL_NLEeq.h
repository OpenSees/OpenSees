
//================================================================================
// COPYRIGHT (C):     :-))                                                       |
// PROJECT:           Object Oriented Finite Element Program                     |
// PURPOSE:           General platform for elaso-plastic constitutive model      |
//                    implementation                                             |
//                                                                               |
// CLASS:             NonLinear EvolutionLaw (on plastic equivalent strain)      |
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
// SHORT EXPLANATION: This is a nonlinear evolution law for the evoltion of a    |
//                    scalar variable k which depends on plastic equi. strain    |
//                    i.e. dk = f( de_eq_p )                                     |
//                                                                               |
//================================================================================

#ifndef EL_NLEeq_H
#define EL_NLEeq_H

#include <math.h>

#include <iostream.h>
#include <iomanip.h>

#include "EL_S.h"

class EvolutionLaw_NL_Eeq : public EvolutionLaw_S
{
  // Private vars to define the evolution law

    //Coefficients to define the nonlinear hardening rule of a scalar var eta
    double eeqEtaPeak, etaResidual, etaStart, etaPeak, e, d;

  public:
    // default constructor
    EvolutionLaw_NL_Eeq( double eeqEtaPeakd = 0.003,
                         double etaResiduald = .2,
                         double etaStartd = 0.3,
                         double etaPeakd = 0.5,
                         double ed = 0.5,
                         double dd = 500.0) :
                       eeqEtaPeak(eeqEtaPeakd), etaResidual(etaResiduald),
                       etaStart(etaStartd), etaPeak(etaPeakd), e(ed), d(dd) {}
                         
    EvolutionLaw_NL_Eeq(const EvolutionLaw_NL_Eeq &NLEL );   // Copy constructor
    
    EvolutionLaw_S *newObj();                     //create a colne of itself
    
    //void InitVars(EPState *EPS);    // Initialize all hardening vars called only once 
    //                                // after material point is formed if necessary.
    
    //void setInitD(EPState  *EPS);   // set initial D once current stress hits the y.s.
    //                                // was primarily for Manzari-Dafalias model

    //double h( EPState *EPS, double d );     // Evaluating hardening function h
    double h_s( EPState *EPS, PotentialSurface *PS);    // Evaluating hardening function h
    
    //void UpdateVar( EPState *EPS, double dlamda );  // Evolve corresponding var linearly using de_eq_p

    void print();

    // some accessor functions
    // some accessor functions
    double geteeqEtaPeak() const;
    double getetaResidual() const;
    double getetaStart() const;
    double getetaPeak() const;
    double gete() const;
    double getd() const;

    //================================================================================
    // Overloaded Insertion Operator	  Zhaohui Added Aug. 13, 2000
    // prints nonlinear EvolutionLaw's contents 
    //================================================================================
    friend ostream& operator<< (ostream& os, const EvolutionLaw_NL_Eeq & NLEL)
    {
        os.unsetf( ios::scientific );
        os.precision(5);
    
        //os.width(10);       
        os << endln << "Nonlinear Evolution(plastic equivalent strain ) Law's parameters:" << endln;
        os << "eeqEtaPeak = " << NLEL.geteeqEtaPeak() << "; ";
        os << "etaResidual = " << NLEL.getetaResidual() << "; " << endln;
        os << "etaStart = " << NLEL.getetaStart() << "; ";
        os << "etaPeak = " << NLEL.getetaPeak() << "; ";
        os << "e = " << NLEL.gete() << "; ";
        os << "d = " << NLEL.getd() << "; " << endln;
        //os.width(10);       
               
        return os;
    }  

    
};


#endif




