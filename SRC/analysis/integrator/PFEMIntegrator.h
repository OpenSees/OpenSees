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

// $Revision: 1.0 $
// $Date: 2013-01-16 12:51:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/PFEMIntegrator.h,v $


#ifndef PFEMIntegrator_h
#define PFEMIntegrator_h

// Written : Minjie Zhu
// Created : Jan 2013
//
// Description: This file contains the class definition for PFEMIntegrator.
// PFEMIntegrator is an algorithmic class for performing a PFEM analysis
// using the Newmark or BackEuler integration scheme.
//
// What: "@(#) PFEMIntegrator.h, revA"

#include <TransientIntegrator.h>
#include <Vector.h>

class DOF_Group;
class FE_Element;

class PFEMIntegrator : public TransientIntegrator
{
public:
    // constructors
    PFEMIntegrator();
    PFEMIntegrator(double gamma, double beta, int disp, int init);

    // destructor
    ~PFEMIntegrator();

    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    int formEleResidual(FE_Element* theEle);
    int formNodUnbalance(DOF_Group* theDof);
    int formTangent(int statFlag);

    int domainChanged(void);
    int newStep(double deltaT);
    int revertToLastStep(void);
    int update(const Vector &deltaU);
    int commit();

    const Vector &getVel(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag = 0);

    // AddingSensitivity:BEGIN //////////////////////////////////
    int revertToStart();
    int formSensitivityRHS(int gradNum);
    int formIndependentSensitivityRHS();
    int formIndependentSensitivityLHS(int statusFlag = CURRENT_TANGENT);
    int saveSensitivity   (const Vector &v, int gradNum, int numGrads);
    int commitSensitivity (int gradNum, int numGrads);
    // AddingSensitivity:END ////////////////////////////////////

protected:
    int displ;      // a flag indicating whether displ(1), vel(2) or accel(3) increments
    int init;  // 1-disp, 2-vel, 3-accel
    double gamma;
    double beta;

    double c1, c2, c3;              // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot;  // response quantities at time t
    Vector *U, *Udot, *Udotdot;     // response quantities at time t+deltaT
    bool determiningMass;           // flag to check if just want the mass contribution

    // Adding Sensitivity
    int sensitivityFlag;
    int gradNumber;
    Vector *massMatrixMultiplicator;
    Vector *dampingMatrixMultiplicator;
    int assemblyFlag;
    Vector independentRHS;
    Vector dUn, dVn, dAn;
    //////////////////////

private:
    int populateU();
    int populateUn();
};

#endif
