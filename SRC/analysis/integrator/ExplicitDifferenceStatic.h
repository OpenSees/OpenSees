/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** ****************************************************************** */

// Written: Jose A. Abell (UANDES)
// Created: 2024
//
// Description: This file contains the class definition for ExplicitDifferenceStatic.
// ExplicitDifferenceStatic is an algorithmic class for performing a transient analysis
// using an explicit central difference scheme with local non-viscous damping.
//
// The scheme uses a leap-frog integration approach with velocity defined at half time steps.
// It incorporates FLAC-style local non-viscous damping that is proportional to the 
// unbalanced force magnitude and opposes velocity. This damping is adaptive and uses
// velocity sign memory to provide stable energy dissipation.
//
// Stability: For an undamped system, dt <= 2/omega_max
//            For a damped system, dt <= 2/(omega_max * sqrt(1 + xi^2) - xi)
//
// References:
// - FLAC3D Manual, Section on Dynamic Analysis
// - Cundall, P. A. (1987). "Distinct element models of rock and soil structure"

#ifndef ExplicitDifferenceStatic_h
#define ExplicitDifferenceStatic_h

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

// Default threshold for velocity sign memory deadband
#define ED_VSIGN_EPS 1e-4

class ExplicitDifferenceStatic : public TransientIntegrator
{
public:
    // Constructors
    ExplicitDifferenceStatic();
    ExplicitDifferenceStatic(double alphaM, double betaK, double betaKi, double betaKc);
    
    // Destructor
    ~ExplicitDifferenceStatic();

    // Methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    
    // Method to obtain current velocity for modal damping
    const Vector &getVel(void);
    
    // Methods to update mesh and define state
    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);
    
    // Method to form nodal unbalance with local non-viscous damping
    int formNodalUnbalance(void);
    
    // Methods for parallel processing
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    // Method to print information
    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    // Time step parameters
    double deltaT;                  // Current time step
    static double deltaT1;          // Previous time step (for variable dt schemes)
    
    // Rayleigh damping coefficients
    double alphaM;                  // Mass proportional damping
    double betaK;                   // Current stiffness proportional damping
    double betaKi;                  // Initial stiffness proportional damping
    double betaKc;                  // Committed stiffness proportional damping
    
    // Integration parameters
    int updateCount;                // Counter for update calls per step
    double c2, c3;                  // Integration constants
    
    // State vectors at different time levels
    Vector *U;                      // Displacement at t+dt (working)
    Vector *Ut;                     // Displacement at t
    Vector *Utdotdot;               // Acceleration at t
    Vector *Utdotdot1;              // Acceleration at t+dt (for output)
    Vector *Udot;                   // Velocity (working)
    Vector *Utdot;                  // Velocity at t+0.5*dt (leap-frog)
    Vector *Utdot1;                 // Velocity at t+dt (for output)
    
    // Local non-viscous damping state variables
    Vector *velSignMem;             // Velocity sign memory per DOF: {-1, 0, +1}
    Vector *prevUnbal;              // Previous step unbalanced force per DOF
    double vSignEps;                // Deadband threshold for velocity sign change
};

#endif
