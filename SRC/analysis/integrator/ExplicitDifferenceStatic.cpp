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

#include <ExplicitDifferenceStatic.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <cmath>
#include <cstring>

#define OPS_Export 

// OPS interface function for creating ExplicitDifferenceStatic integrator
void* OPS_ExplicitDifferenceStatic(void)
{
    ExplicitDifferenceStatic *theIntegrator = new ExplicitDifferenceStatic();

    // integrator ExplicitDifferenceStatic <-alpha $alpha> <-simple>
    double alpha = 0.59;
    bool combined = true;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *opt = OPS_GetString();
        if (strcmp(opt, "-alpha") == 0) {
            int numData = 1;
            if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &alpha) < 0 ||
                alpha < 0.0 || alpha >= 1.0) {
                opserr << "WARNING ExplicitDifferenceStatic -alpha needs a value in [0, 1)\n";
                delete theIntegrator;
                return 0;
            }
        } else if (strcmp(opt, "-simple") == 0) {
            combined = false;
        } else {
            opserr << "WARNING ExplicitDifferenceStatic - unknown option " << opt << "\n";
            delete theIntegrator;
            return 0;
        }
    }
    theIntegrator->setLocalDamping(alpha, combined);

    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitDifferenceStatic integrator\n";

    return theIntegrator;
}

// Default constructor
ExplicitDifferenceStatic::ExplicitDifferenceStatic()
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitDifferenceStatic),
    deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    updateCount(0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0),
    Udot(0), Utdotdot1(0), U(0), Utdot1(0), 
    velSignMem(0), prevUnbal(0), vSignEps(ED_VSIGN_EPS), alphaLNVD(0.59), useCombined(true)
{
}

// Constructor with Rayleigh damping parameters
ExplicitDifferenceStatic::ExplicitDifferenceStatic(
    double _alphaM, double _betaK, double _betaKi, double _betaKc)
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitDifferenceStatic),
    deltaT(0.0),
    alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    updateCount(0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0),
    Udot(0), Utdotdot1(0), U(0), Utdot1(0), 
    velSignMem(0), prevUnbal(0), vSignEps(ED_VSIGN_EPS), alphaLNVD(0.59), useCombined(true)
{
}

// Destructor - clean up allocated memory
ExplicitDifferenceStatic::~ExplicitDifferenceStatic()
{
    if (Ut != 0)
        delete Ut;
    if (Utdot != 0)
        delete Utdot;
    if (Utdotdot != 0)
        delete Utdotdot;
    if (Udot != 0)
        delete Udot;
    if (Utdotdot1 != 0)
        delete Utdotdot1;
    if (U != 0)
        delete U;
    if (Utdot1 != 0)
        delete Utdot1;
    if (velSignMem != 0)
        delete velSignMem;
    if (prevUnbal != 0)
        delete prevUnbal;
}

// Advance to a new time step
// This method implements the leap-frog scheme:
// - Velocity is stored at half time steps: v_{t+0.5*dt}
// - Displacement is at full time steps: u_t, u_{t+dt}
// - Acceleration is at full time steps: a_t, a_{t+dt}
int ExplicitDifferenceStatic::newStep(double _deltaT)
{
    updateCount = 0;
    deltaT = _deltaT;

    if (deltaT <= 0.0)  {
        opserr << "ExplicitDifferenceStatic::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();

    // Leap-frog update: advance velocity to t+0.5*dt and displacement to t+dt
    Utdot->addVector(1.0, *Utdotdot, deltaT);           // v_{t+0.5*dt} += a_t * dt
    Ut->addVector(1.0, *Utdot, deltaT);                 // u_{t+dt} = u_t + v_{t+0.5*dt} * dt

    if (Ut == 0)  {
        opserr << "ExplicitDifferenceStatic::newStep() - domainChange() failed or hasn't been called\n";
        return -2;
    }

    // Zero acceleration for explicit method (Ma = RHS, no Ma on LHS)
    (*Utdotdot) *= 0.0;

    // Set trial response quantities
    theModel->setVel(*Utdot);
    theModel->setAccel(*Utdotdot);
    theModel->setDisp(*Ut);

    // Increment time and apply loads
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "ExplicitDifferenceStatic::newStep() - failed to update the domain\n";
        return -3;
    }

    // Set response at t to be that at t+deltaT of previous step
    (*Utdotdot) = (*Utdotdot1);
    
    return 0;
}

// Form element tangent matrix (mass matrix for explicit scheme)
int ExplicitDifferenceStatic::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    theEle->addMtoTang();
    return 0;
}

// Form nodal tangent matrix (mass matrix for explicit scheme)
int ExplicitDifferenceStatic::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    theDof->addMtoTang();
    return 0;
}

// Handle domain changes (mesh modifications, initial conditions, etc.)
int ExplicitDifferenceStatic::domainChanged()
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();

    opserr << "ExplicitDifferenceStatic::domainChanged()" << endln;

    // Set Rayleigh damping factors if specified
    if (alphaM != 0.0 || betaK != 0.0 || betaKi != 0.0 || betaKc != 0.0)
        theModel->setRayleighDampingFactors(alphaM, betaK, betaKi, betaKc);

    // Allocate or reallocate vectors if size has changed
    if (Ut == 0 || Ut->Size() != size)  {
        // Delete old vectors
        if (Ut != 0)
            delete Ut;
        if (Utdot != 0)
            delete Utdot;
        if (Utdotdot != 0)
            delete Utdotdot;
        if (Udot != 0)
            delete Udot;
        if (Utdotdot1 != 0)
            delete Utdotdot1;
        if (U != 0)
            delete U;
        if (Utdot1 != 0)
            delete Utdot1;
        if (velSignMem != 0)
            delete velSignMem;
        if (prevUnbal != 0)
            delete prevUnbal;

        // Create new vectors
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        Udot = new Vector(size);
        U = new Vector(size);
        Utdotdot1 = new Vector(size);
        Utdot1 = new Vector(size);
        velSignMem = new Vector(size);
        prevUnbal = new Vector(size);

        // Verify allocation was successful
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            U == 0 || U->Size() != size ||
            Utdotdot1 == 0 || Utdotdot1->Size() != size ||
            Utdot1 == 0 || Utdot1->Size() != size  || 
            velSignMem->Size() != size || 
            prevUnbal->Size() != size)  {

            opserr << "ExplicitDifferenceStatic::domainChanged - ran out of memory\n";

            // Clean up on failure
            if (Ut != 0)
                delete Ut;
            if (Utdot != 0)
                delete Utdot;
            if (Utdotdot != 0)
                delete Utdotdot;
            if (Udot != 0)
                delete Udot;
            if (U != 0)
                delete U;
            if (Utdotdot1 != 0)
                delete Utdotdot1;
            if (Utdot1 != 0)
                delete Utdot1;
            if (velSignMem != 0)
                delete velSignMem;
            if (prevUnbal != 0)
                delete prevUnbal;

            Ut = 0; Utdot = 0; Utdotdot = 0;
            Udot = 0; U = 0, Utdotdot1 = 0;
            Utdot1 = 0; velSignMem = 0; prevUnbal = 0;
        
            return -1;
        }
    }

    // Initialize state vectors from committed DOF values
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0)  {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        // Initialize displacements
        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {            
                (*Ut)(loc) = disp(i);
            }
        }

        // Initialize velocities
        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Utdot)(loc) = vel(i);
                (*Utdot1)(loc) = vel(i);
            }
        }

        // Initialize accelerations
        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Utdotdot)(loc) = accel(i);
                (*Utdotdot1)(loc) = accel(i);
            }
        }
    }

    return 0;
}

// Form nodal unbalance with local non-viscous damping
// 
// This method implements FLAC-style adaptive local non-viscous damping:
// F_damped = F_unbalanced + F_d
// 
// where F_d can be either:
// 1. Simple form: F_d = -alpha * |F_unbalanced| * sign(v)
// 2. Combined form: F_d = 0.5 * alpha * |F_unbalanced| * (sign(dF/dt) - sign(v))
// 
// The combined form provides better energy dissipation and stability.
// Velocity sign memory with deadband prevents sign chatter near zero velocity.
void ExplicitDifferenceStatic::setLocalDamping(double alpha, bool combined)
{
    alphaLNVD = alpha;
    useCombined = combined;
}

// Local non-viscous damping (Cundall 1987; FLAC), applied equation by equation to
// the fully assembled unbalance B = P - R(u). (DOF_Group::getUnbalance only holds
// applied nodal loads, so damping must act on the assembled vector.)
//   simple:   F_d = -alpha |F| sign(v)
//   combined: F_d = 0.5 alpha |F| (sign(dF) - sign(v))
// v is the leap-frog velocity at t + dt/2; sign(v) is held inside a small deadband.
int ExplicitDifferenceStatic::formUnbalance(void)
{
    int res = this->TransientIntegrator::formUnbalance();
    if (res < 0 || alphaLNVD <= 0.0 || Utdot == 0 || velSignMem == 0 || prevUnbal == 0)
        return res;

    LinearSOE *theLinSOE = this->getLinearSOE();
    static Vector B;
    B = theLinSOE->getB();
    const int n = B.Size();
    if (Utdot->Size() != n || velSignMem->Size() != n || prevUnbal->Size() != n)
        return res;

    for (int eq = 0; eq < n; ++eq) {
        const double v = (*Utdot)(eq);
        const double F = B(eq);

        double s_v = (*velSignMem)(eq);
        if (std::abs(v) > vSignEps) {
            s_v = (v > 0.0) ? 1.0 : -1.0;
            (*velSignMem)(eq) = s_v;
        }

        double Fd;
        if (useCombined) {
            const double dF = F - (*prevUnbal)(eq);
            const double s_fdot = (dF > 0.0) ? 1.0 : ((dF < 0.0) ? -1.0 : 0.0);
            Fd = 0.5 * alphaLNVD * std::abs(F) * (s_fdot - s_v);
        } else {
            Fd = -alphaLNVD * std::abs(F) * s_v;
        }
        (*prevUnbal)(eq) = F;
        B(eq) = F + Fd;
    }
    return theLinSOE->setB(B);
}

// Update the response quantities
int ExplicitDifferenceStatic::update(const Vector &Udotdot)
{
    updateCount++;
    if (updateCount > 2)  {
        opserr << "WARNING ExplicitDifferenceStatic::update() - called more than once -";
        opserr << " ExplicitDifferenceStatic integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING ExplicitDifferenceStatic::update() - no AnalysisModel set\n";
        return -2;
    }

    if (Ut == 0)  {
        opserr << "WARNING ExplicitDifferenceStatic::update() - domainChange() failed or not called\n";
        return -3;
    }

    if (Udotdot.Size() != Utdotdot->Size()) {
        opserr << "WARNING ExplicitDifferenceStatic::update() - Vectors of incompatible size ";
        opserr << " expecting " << Utdotdot->Size() << " obtained " << Udotdot.Size() << endln;
        return -4;
    }

    // Update acceleration: weighted average for stability
    // a_{t+dt} = (3*a_new + a_old) / 4
    double halfT = deltaT * 0.125;

    Utdotdot1->addVector(0.0, Udotdot, 3.0);
    Utdotdot1->addVector(1.0, *Utdotdot, 1.0);

    // Update velocity for output (v at t+dt from leap-frog v at t+0.5*dt)
    Utdot1->addVector(0.0, *Utdot, 1.0);
    Utdot1->addVector(1.0, *Utdotdot1, halfT);

    // Set response in model
    theModel->setResponse(*Ut, *Utdot1, Udotdot);

    if (theModel->updateDomain() < 0)  {
        opserr << "ExplicitDifferenceStatic::update() - failed to update the domain\n";
        return -5;
    }

    // Store acceleration for next step
    (*Utdotdot) = Udotdot;
    (*Utdotdot1) = Udotdot;

    return 0;
}

// Commit the state
int ExplicitDifferenceStatic::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING ExplicitDifferenceStatic::commit() - no AnalysisModel set\n";
        return -1;
    }

    return theModel->commitDomain();
}

// Send object state for parallel processing
int ExplicitDifferenceStatic::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(6);
    data(0) = alphaM;
    data(1) = betaK;
    data(2) = betaKi;
    data(3) = betaKc;
    data(4) = alphaLNVD;
    data(5) = useCombined ? 1.0 : 0.0;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING ExplicitDifferenceStatic::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}

// Receive object state for parallel processing
int ExplicitDifferenceStatic::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(6);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING ExplicitDifferenceStatic::recvSelf() - could not receive data\n";
        return -1;
    }

    alphaM = data(0);
    betaK = data(1);
    betaKi = data(2);
    betaKc = data(3);
    alphaLNVD = data(4);
    useCombined = (data(5) != 0.0);

    return 0;
}

// Print integrator information
void ExplicitDifferenceStatic::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "ExplicitDifferenceStatic - currentTime: " << currentTime << endln;
        s << "  Rayleigh Damping - alphaM: " << alphaM << "  betaK: " << betaK;
        s << "  betaKi: " << betaKi << "  betaKc: " << betaKc << endln;
        s << "  Local non-viscous damping coefficient: 0.59" << endln;
    }
    else
        s << "ExplicitDifferenceStatic - no associated AnalysisModel\n";
}

// Get current velocity (for modal damping interface)
const Vector &ExplicitDifferenceStatic::getVel()
{
    return *Utdot;
}
