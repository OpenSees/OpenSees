/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.1 $
// $Date: 2009-03-20 18:36:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/Houbolt.h,v $

#ifndef Houbolt_h
#define Houbolt_h

// Written : krm
// Created : 11/2012
//
// Description: This file contains the class definition for Houbolt.
// Houbolt is an algorithmic class for performing a transient analysis
// using Houbolt linear multistep method.
// ref: J.C.Houbolt "A recurrence matrix solution for the dynamic response of 
//      elastic aircraft". J. Aeronaut. Sci. 17 (1950): 540–550.
// It is 2nd order accurate and unconditionally stable, but may dissipate
// too much energy in low frequency range for earthquake problems.  Useful 
// for other dynamic problems, however.
//
// Since it is not necessarily selfstarting, we bootstrap using steps of 
// the trapezoid rule. Note this should not be used with variable step size 
// otherwise frequent step size changes will render this a trapezoidal integrator.
//
// What: "@(#) Houbolt.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class Houbolt : public TransientIntegrator
{
public:
    Houbolt();
    ~Houbolt();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);        
    
    int domainChanged(void);    
    int newStep(double deltaT);    
    int revertToLastStep(void);        
    int update(const Vector &deltaU);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);        
    
 protected:

 private:
    int step;       // keep track of previous points performed
    double dt;      // last dt, if not same as previous we do trapezoidal step
    
    double c1, c2, c3;                  // some constants we need to keep
    Vector *Utm2;                       // response quantities at time t-2
    Vector *Utm1;                       // response quantities at time t-1
    Vector *Ut, *Utdot, *Utdotdot;      // response quantities at time t
    Vector *U, *Udot, *Udotdot;         // response quantities at time t+deltaT
};

#endif
