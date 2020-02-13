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

#ifndef KRAlphaExplicit_h
#define KRAlphaExplicit_h

// Developed: Chinmoy Kolay (chk311@lehigh.edu)
// Implemented: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 08/14
// Revision: A
//
// Description: This file contains the class definition for KRAlphaExplicit.
// KRAlphaExplicit is an algorithmic class for performing a transient analysis
// using the explicit Kolay-Ricles integration scheme based on the midpoint rule.
//
// Reference: Kolay, C. and J. Ricles (2014). "Development of a family of
// unconditionally stable explicit direct integration algorithms with
// controllable numerical energy dissipation." Earthquake Engineering and
// Structural Dynamics, 43(9):1361–1380.

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;
class Matrix;

class KRAlphaExplicit : public TransientIntegrator
{
public:
    // constructors
    KRAlphaExplicit();
    KRAlphaExplicit(double rhoInf,
        bool updElemDisp = false);
    
    // destructor
    ~KRAlphaExplicit();
    
    // method to set up the system of equations
    int formTangent(int statFlag);
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    
    // methods to update the domain
    int domainChanged(void);
    int newStep(double deltaT);
    int revertToLastStep(void);
    int update(const Vector &aiPlusOne);
    int commit(void);

    const Vector &getVel(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:
    
private:
    double alphaM;
    double alphaF;
    double beta;
    double gamma;
    bool updElemDisp;  // a flag indicating if element displacements are updated during commit
    double deltaT;
    
    Matrix *alpha1, *alpha3;  // integration parameter matrices, alpha2 = (0.5 + gamma)*alpha1
    Matrix *Mhat;             // effective mass matrix for linear SOE
    
    int updateCount;                            // method should only have one update per step
    int initAlphaMatrices;                      // a flag to initialize the alpha matrices
    double c1, c2, c3;                          // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot;              // response quantities at time t
    Vector *U, *Udot, *Udotdot;                 // response quantities at time t + deltaT
    Vector *Ualpha, *Ualphadot, *Ualphadotdot;  // response quantities at time t + alpha*deltaT
    Vector *Utdothat;                           // velocity-like vector at time t
};

#endif
