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
// $Date: 2009/03/20 22:37:09 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewIntegrator/Trapezoidal.h,v $

#ifndef Trapezoidal_h
#define Trapezoidal_h

// Written : fmk 
// Created : 02/09
//
// Description: This file contains the class definition for Trapezoidal.
// Trapezoidal is a transient integrator which uses trapezoidal rule for
// Ut+deltaT and Vt+deltaT (U,V,A the displ, vel and accel)
//   Ut+deltaT = Ut + [deltaT/2.0] * (Vt + Vt+deltaT) 
//   Vt+deltaT = Vt + [deltaT/2.0] * (At + At+deltaT) 
// and solves transient equation at time t+deltaT.
//
// It is the numerically the same as Newmark with gamma=0.5, beta=0.25 (average acel method)

// What: "@(#) Trapezoidal.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class Trapezoidal : public TransientIntegrator
{
public:
    Trapezoidal();
    ~Trapezoidal();
    
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
    int revertToStart();
    
 protected:

 private:
    double c1, c2, c3;              // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot;         // response quantities at time t
    Vector *U, *Udot, *Udotdot;            // response quantities at time t+deltaT
};

#endif
