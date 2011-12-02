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
// $Date: 2000-12-13 08:27:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/HHT.h,v $
                                                                        
                                                                        
#ifndef HHT_h
#define HHT_h

// File: ~/analysis/integrator/HHT.h
// 
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for HHT.
// HHT is an algorithmic class for performing a transient analysis
// using the HHT integration scheme.
//
// What: "@(#) HHT.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class HHT : public TransientIntegrator
{
  public:
    HHT();
    HHT(double alpha);
    HHT(double alpha, double alphaM, double betaKcurrent, 
	double betaKinit, double betaKlastCommit);        

    ~HHT();

    // Override residual calculation needed for Rayleigh damping
    virtual int formEleResidual(FE_Element *theEle);
    virtual int formNodUnbalance(DOF_Group *theDof);        
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);        

    int domainChanged(void);    
    int initialize(void);    
    int newStep(double deltaT);    
    int revertToLastStep(void);        
    int update(const Vector &deltaU);

    int commit(void);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(ostream &s, int flag =0);        
    
  protected:
    
  private:
    double alpha;
    double gamma;
    double beta;
    bool rayleighDamping; // flag indicating wheter rayleigh damping being used

    // D = alphaM * M + betaK * K(current) + betaKi * Kinit + betaKc * Kcommit
    double alphaM, betaK, betaKi, betaKc;      // raleigh factors 

    double c1, c2, c3;  // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot; // response quantities at time t
    Vector *U, *Udot, *Udotdot; // response quantities at time t + delta t
    Vector *Ualpha, *Udotalpha; // response quantities at time t+alpha delta t
};

#endif

