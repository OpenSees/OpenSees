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

// Written: Jose A. Abell (UANDES) & Massimo Petracca (ASDEA)
// Created: 3 December 2024
// Revision: A
//
// Description: This file contains the class definition for ExplicitBathe.
// ExplicitBathe is an algorithmic class for performing a transient analysis
// using the explicit Bathe time integration scheme. 
//
// This algorithm is similar to a central difference scheme but perturbed to 
// introduce numerical damping. It is a second-order accurate explicit scheme. 
// Unlike the CentralDifference class, this one only assembles the mass matrix 
// on the right hand side, making it much easier to use with diagonal mass matrices. 
//
// The time-step required for stability is approximately twice that of 
// standard explicit central difference schemes, making it more efficient for 
// many applications.
//
// Key features:
// - Second-order accurate
// - Conditionally stable (dt <= 2/omega_max for undamped systems)
// - Built-in numerical damping (controlled by parameter p, p=0.54 is a good choice)
// - Only mass matrix on RHS (no tangent matrix assembly required)
// - Optional automatic critical time step calculation
//
// Reference:
// Gunwoo Noh, Klaus-Jürgen Bathe,
// "An explicit time integration scheme for the analysis of wave propagations",
// Computers & Structures, Volume 129, 2013, Pages 178-193,
// ISSN 0045-7949,
// https://doi.org/10.1016/j.compstruc.2013.06.007.

#ifndef ExplicitBathe_h
#define ExplicitBathe_h

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class ExplicitBathe : public TransientIntegrator
{
public:
    // Constructors
    ExplicitBathe();
    
    ExplicitBathe(double p, int compute_critical_timestep_ = 0);
    
    // Destructor
    ~ExplicitBathe();
    
    // Methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    
    // Methods to update mesh and define state
    int domainChanged(void);
    int newStep(double deltaT);
    int update(const Vector &U);
    int commit(void);

    // Method to obtain current velocity
    const Vector &getVel(void);
    
    // Methods for parallel processing
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    // Method to print information
    void Print(OPS_Stream &s, int flag = 0);

protected:

private:
    // Time step
    double deltaT;

    // Explicit Bathe parameters
    double p;           // Sub-step parameter (0 < p < 1), controls damping
    double q0;          // Integration coefficient q0
    double q1;          // Integration coefficient q1
    double q2;          // Integration coefficient q2

    // State vectors at different time levels
    // The scheme has two sub-steps per time step:
    // Step 1: t -> t + p*dt (using acceleration A_tpdt)
    // Step 2: t + p*dt -> t + dt (using accelerations A_tpdt and A_tdt)
    
    Vector *U_t;        // Displacement at time t
    Vector *V_t;        // Velocity at time t 
    Vector *A_t;        // Acceleration at time t 
    
    Vector *U_tpdt;     // Displacement at time t + p*dt
    Vector *V_tpdt;     // Velocity at time t + p*dt
    Vector *V_fake;     // Temporary velocity for setting response
    Vector *A_tpdt;     // Acceleration at time t + p*dt
    
    Vector *U_tdt;      // Displacement at time t + dt
    Vector *V_tdt;      // Velocity at time t + dt
    Vector *A_tdt;      // Acceleration at time t + dt
    Vector *R_tdt;      // Forces at time t + dt (unused, kept for compatibility)

    // Update counter
    int updateCount;    // Counts updates per step (should be exactly 2)

    // Integration coefficients (computed from p)
    double a0;          // = p * dt
    double a1;          // = (p * dt)^2 / 2
    double a2;          // = a0 / 2
    double a3;          // = (1 - p) * dt
    double a4;          // = ((1 - p) * dt)^2 / 2
    double a5;          // = q0 * a3
    double a6;          // = (0.5 + q1) * a3
    double a7;          // = q2 * a3

    // Critical time step computation
    int compute_critical_timestep;          // Control flag for critical dt computation
    double damped_minimum_critical_timestep;    // Minimum critical dt (damped)
    double undamped_minimum_critical_timestep;  // Minimum critical dt (undamped)
    int damped_critical_element_tag;            // Element with minimum damped dt
    int undamped_critical_element_tag;          // Element with minimum undamped dt
};

#endif
