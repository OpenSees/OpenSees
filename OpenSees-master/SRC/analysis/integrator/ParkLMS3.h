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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/ParkLMS3.h,v $

#ifndef ParkLMS3_h
#define ParkLMS3_h

// Written : krm
// Created : 11/2012
//
// Description: This file contains the class definition for ParkLMS3.
// ParkLMS3 is an algorithmic class for performing a transient analysis
// using Park 1975 linear multistep (3 step) method.
// ref: K.C.Park "Evaluating time integration methods for nonlinear dynamic 
//      analysis" in: T. Belytschko, T.L. Geers (Eds.), Finite Element Analysis 
//      of Transient Nonlinear Behavior, AMD 14, ASME, New York (1975): 35–58.
// It is 2nd order accurate and unconditionally stable.
//
// Since it is not necessarily selfstarting, we bootstrap using steps of 
// the trapezoid rule. Note this should not be used with variable step size 
// otherwise frequent step size changes will render this a trapezoidal integrator.
//
// What: "@(#) ParkLMS3.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class ParkLMS3 : public TransientIntegrator
{
public:
    ParkLMS3();
    ~ParkLMS3();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);        
    
    int domainChanged(void);    
    int newStep(double deltaT);    
    int revertToLastStep(void);        
    int update(const Vector &deltaU);

    const Vector &getVel(void);
    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);        
    
 protected:

 private:
    int step;       // keep track of previous points performed
    double dt;      // last dt, if not same as previous we do trapezoidal step
    
    double c1, c2, c3;                  // some constants we need to keep
    Vector *Utm2, *Utm2dot;             // response quantities at time t-2
    Vector *Utm1, *Utm1dot;             // response quantities at time t-1
    Vector *Ut, *Utdot, *Utdotdot;      // response quantities at time t
    Vector *U, *Udot, *Udotdot;         // response quantities at time t+deltaT
};

#endif
