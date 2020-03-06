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

// $Revision$
// $Date$
// $URL$

#ifndef CentralDifference_h
#define CentralDifference_h

// Written: fmk
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for CentralDifference.
// CentralDifference is an algorithmic class for performing a transient analysis
// using the central difference integration scheme.

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class CentralDifference : public TransientIntegrator
{
public:
    // constructors
    CentralDifference();
    CentralDifference(double alphaM, double betaK, double betaKi, double betaKc);
    
    // destructor
    ~CentralDifference();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    
    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    const Vector &getVel(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:
    
private:
    double deltaT;
    
    // rayleigh damping factors
    double alphaM;
    double betaK;
    double betaKi;
    double betaKc;
    
    int updateCount;                // method should only have one update per step
    double c2, c3;                  // some constants we need to keep
    Vector *Utm1;                   // disp response quantity at time t-deltaT
    Vector *Ut, *Utdot, *Utdotdot;  // response quantities at time t
    Vector *Udot, *Udotdot;         // response quantities at time t+deltaT
};

#endif
