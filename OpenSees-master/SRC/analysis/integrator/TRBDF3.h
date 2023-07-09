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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/TRBDF3.h,v $

#ifndef TRBDF3_h
#define TRBDF3_h

// Written : krm
// Created : 11/2012
//
// Description: This file contains the class definition for TRBDF3.
// TRBDF3 is an algorithmic class for performing a transient analysis
// using the so-called TRBDF3 integration scheme. It is based in part on the 
// two part algorithm of Bathe that is now coded as TRBDF2. Note that TRBDF2 
// is actually described in an earlier paper by Bathe, as well as the basis 
// for TRBDF3 used here:
// ref: K.J.Bathe and M.M.I. Baige, "On a composite implicit time integration procedure
//      for nonlinear dynamics", Computers and Structures 83(2005),2513-2524
//
// But here we further subdivide the step into a total of 3 sub-steps with the first being
// trapezoid rule, second three point backward Euler, and third is Houbolt (although 
// technically he notes in the paper that he uses a 3 sub-step procedure with 2 trapezoid
// steps followed by Houbolt).
//
// note: the implementation does not do sub-step, it just cycles between trapezoidal, 
// euler methods, and Houbolt. So to compare to paper or to Newmark or some other analysis, 
// user should specify dt/3 step size.
//
// What: "@(#) TRBDF3.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class TRBDF3 : public TransientIntegrator
{
public:
    TRBDF3();
    ~TRBDF3();
    
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
    int step;      // a flag indicating whether trap or euler step
    double dt;     // last dt, if not same as previous we do trapezoidal step
    
    double c1, c2, c3;                  // some constants we need to keep
    Vector *Utm2, *Utm2dot;             // response quantities at time t-2
    Vector *Utm1, *Utm1dot;             // response quantities at time t-1
    Vector *Ut, *Utdot, *Utdotdot;      // response quantities at time t
    Vector *U, *Udot, *Udotdot;         // response quantities at time t+deltaT
};

#endif
