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

// $Revision: 1.9 $
// $Date: 2010-02-04 01:06:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/Newmark.h,v $


#ifndef Newmark_h
#define Newmark_h

// Written : fmk 
// Created : 11/98
// Modified: 02/05 ahs
// Revision: A
//
// Description: This file contains the class definition for Newmark.
// Newmark is an algorithmic class for performing a transient analysis
// using the Newmark integration scheme.
//
// What: "@(#) Newmark.h, revA"

#include <TransientIntegrator.h>
#include <Vector.h>

class DOF_Group;
class FE_Element;

class Newmark : public TransientIntegrator
{
public:
    // constructors
    Newmark();
    Newmark(double gamma, double beta, bool disp = true, bool aflag=false);

    // destructor
    ~Newmark();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    int formEleResidual(FE_Element* theEle);
    int formNodUnbalance(DOF_Group* theDof);
    
    int domainChanged(void);    
    int newStep(double deltaT);    
    int revertToLastStep(void);        
    int update(const Vector &deltaU);

    double getCFactor(void);

    const Vector &getVel(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);        
    
    // AddingSensitivity:BEGIN //////////////////////////////////
    int revertToStart();
    int formSensitivityRHS(int gradNum);
    int formIndependentSensitivityRHS();
    int saveSensitivity   (const Vector &v, int gradNum, int numGrads);
    int commitSensitivity (int gradNum, int numGrads);  
   int computeSensitivities( );
    // AddingSensitivity:END ////////////////////////////////////
    
protected:
    bool displ;      // a flag indicating whether displ or accel increments
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
};

#endif
