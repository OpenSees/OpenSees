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

// $Revision: 1.2 $
// $Date: 2005-12-21 00:31:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/AlphaOS.h,v $


#ifndef AlphaOS_h
#define AlphaOS_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/05
// Revision: A
//
// Description: This file contains the class definition for AlphaOS.
// AlphaOS is an algorithmic class for performing a transient analysis
// using the Alpha-Operator-Splitting integration scheme.
// The parameter alpha corresponds to 1+alpha_{HHT}.
//
// What: "@(#) AlphaOS.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class AlphaOS : public TransientIntegrator
{
public:
    // constructors
    AlphaOS();
    AlphaOS(double alpha);
    AlphaOS(double alpha, 
        double alphaM, double betaK, double betaKi, double betaKc);
    AlphaOS(double alpha, double beta, double gamma);
    AlphaOS(double alpha, double beta, double gamma,
        double alphaM, double betaK, double betaKi, double betaKc);

    // destructor
    ~AlphaOS();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);       
    
    int domainChanged(void);
    int newStep(double deltaT);    
    int revertToLastStep(void);
    int update(const Vector &deltaU);
    int commit(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);        
    
protected:
    virtual int formElementResidual(void);
    
private:
    double alpha;
    double beta;
    double gamma;
    double deltaT;
    
    // rayleigh damping factors
    double alphaM;
    double betaK;
    double betaKi;
    double betaKc;
    
    int updateCount;                // method should only have one update per step
    double c1, c2, c3;              // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot;  // response quantities at time t
    Vector *U, *Udot, *Udotdot;     // response quantities at time t+deltaT
    Vector *Ualpha, *Ualphadot;     // response quantities at time t+alpha*deltaT
    Vector *Upt, *Uptdot;           // predictor quantities at time t
};

#endif
