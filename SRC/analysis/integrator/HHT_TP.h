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

#ifndef HHT_TP_h
#define HHT_TP_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 08/15
// Revision: A
//
// Description: This file contains the class definition for HHT_TP.
// HHT_TP is an algorithmic class for performing a transient analysis
// using the HHT integration scheme based on the trapezoidal rule.

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class HHT_TP : public TransientIntegrator
{
public:
    // constructors
    HHT_TP();
    HHT_TP(double alpha);
    HHT_TP(double alpha, double beta, double gamma);
    
    // destructor
    ~HHT_TP();
    
    // method to set up the system of equations
    int formUnbalance(void);
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    int formEleResidual(FE_Element *theEle);
    int formNodUnbalance(DOF_Group *theDof);
    
    // methods to update the domain
    int domainChanged(void);
    int newStep(double deltaT);
    int revertToLastStep(void);
    int update(const Vector &deltaU);
    int commit(void);

    const Vector &getVel(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:
    
private:
    double alpha;
    double beta;
    double gamma;
    double deltaT;
    
    double c1, c2, c3;                      // some constants we need to keep
    double alphaM, alphaD, alphaR, alphaP;  // weighting factors we need to keep
    Vector *Ut, *Utdot, *Utdotdot;          // response quantities at time t
    Vector *U, *Udot, *Udotdot;             // response quantities at time t + deltaT
    Vector *Put;                            // unbalance at time t
};

#endif
