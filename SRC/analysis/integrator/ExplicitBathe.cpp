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

// Written: Jose A. Abell (UANDES) & Massimo Petracca (ASDEA)
// Created: 3 December 2024

#include <ExplicitBathe.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <cmath>
#include <limits>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>

#define OPS_Export

// LAPACK eigenvalue solver declarations
#ifdef _WIN32
extern "C" int DGGEV(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                     double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                     double *BETA, double *VL, int *LDVL, double *VR,
                     int *LDVR, double *WORK, int *LWORK, int *INFO);
#else
extern "C" int dggev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA,
                      double *B, int *LDB, double *ALPHAR, double *ALPHAI,
                      double *BETA, double *VL, int *LDVL, double *VR,
                      int *LDVR, double *WORK, int *LWORK, int *INFO);
#endif

// OPS interface function for creating ExplicitBathe integrator
void *OPS_ExplicitBathe(void) {
    TransientIntegrator *theIntegrator = nullptr;

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "WARNING: Insufficient arguments for ExplicitBathe integrator.\n";
        opserr << "Usage: integrator ExplicitBathe p <compute_critical_timestep>\n";
        opserr << "  p = damping parameter (0 < p < 1, typically 0.5-0.95)\n";
        opserr << "  compute_critical_timestep = optional flag (0 or 1)\n";
        return nullptr;
    }

    // Read damping parameter p
    double p;
    numArgs = 1;
    if (OPS_GetDoubleInput(&numArgs, &p) < 0) {
        opserr << "WARNING: Invalid input for ExplicitBathe integrator. p parameter\n";
        return nullptr;
    }

    // Validate p parameter
    if (p <= 0.0 || p >= 1.0) {
        opserr << "WARNING: ExplicitBathe parameter p must be in range (0, 1).\n";
        opserr << "  Typical values: 0.5-0.95 (higher = more damping)\n";
        return nullptr;
    }

    // Read optional critical timestep computation flag
    int compute_critical_timestep = 0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
        if (OPS_GetIntInput(&numArgs, &compute_critical_timestep) < 0) {
            opserr << "WARNING: Invalid compute_critical_timestep parameter\n";
            return nullptr;
        }
    }

    // Create integrator
    theIntegrator = new ExplicitBathe(p, compute_critical_timestep);

    if (theIntegrator == nullptr) {
        opserr << "WARNING - out of memory creating ExplicitBathe integrator\n";
    }
    
    return theIntegrator;
}

// Default constructor
ExplicitBathe::ExplicitBathe()
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(0.0), q0(0.0), q1(0.0), q2(0.0), 
      U_t(0), V_t(0), A_t(0),
      U_tpdt(0), V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0), V_tdt(0), A_tdt(0), R_tdt(0),
      updateCount(0),
      a0(0.), a1(0.), a2(0.), a3(0.), a4(0.), a5(0.), a6(0.), a7(0.),
      compute_critical_timestep(0),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0)
{}

// Main constructor with parameters
//
// Integration coefficients q0, q1, q2 are computed from p:
// q1 = (1 - 2p) / (2p(1-p))
// q2 = 0.5 - p*q1
// q0 = -q1 - q2 + 0.5
ExplicitBathe::ExplicitBathe(double _p, int compute_critical_timestep_)
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(_p), q0(0.0), q1(0.0), q2(0.0), 
      U_t(0), V_t(0), A_t(0),
      U_tpdt(0), V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0), V_tdt(0), A_tdt(0), R_tdt(0),
      updateCount(0),
      a0(0.), a1(0.), a2(0.), a3(0.), a4(0.), a5(0.), a6(0.), a7(0.),
      compute_critical_timestep(compute_critical_timestep_),
      damped_minimum_critical_timestep(0.0),
      undamped_minimum_critical_timestep(0.0),
      damped_critical_element_tag(0),
      undamped_critical_element_tag(0)
{
    // Calculate integration coefficients from p parameter
    q1 = (1.0 - 2.0*p) / (2.0*p*(1.0 - p));
    q2 = 0.5 - p * q1;
    q0 = -q1 - q2 + 0.5;

    opserr << "ExplicitBathe: p = " << p 
           << ", compute_critical_timestep = " << compute_critical_timestep << endln;
}

// Destructor - clean up allocated memory
ExplicitBathe::~ExplicitBathe() {
    if (U_t) delete U_t;
    if (V_t) delete V_t;
    if (A_t) delete A_t;
    if (U_tpdt) delete U_tpdt;
    if (V_tpdt) delete V_tpdt;
    if (V_fake) delete V_fake;
    if (A_tpdt) delete A_tpdt;
    if (U_tdt) delete U_tdt;
    if (V_tdt) delete V_tdt;
    if (A_tdt) delete A_tdt;
    if (R_tdt) delete R_tdt;
}

// Handle domain changes (mesh modifications, initial conditions, etc.)
int ExplicitBathe::domainChanged() {
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    if (!theModel || !theLinSOE) {
        opserr << "ExplicitBathe::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    int size = x.Size();

    if (size == 0) {
        opserr << "ExplicitBathe::domainChanged - invalid size\n";
        return -1;
    }

    // Allocate or reallocate state vectors if size has changed
    if (!U_t || U_t->Size() != size) {
        // Delete old vectors
        if (U_t) delete U_t;
        if (V_t) delete V_t;
        if (A_t) delete A_t;
        if (U_tpdt) delete U_tpdt;
        if (V_tpdt) delete V_tpdt;
        if (V_fake) delete V_fake;
        if (A_tpdt) delete A_tpdt;
        if (U_tdt) delete U_tdt;
        if (V_tdt) delete V_tdt;
        if (A_tdt) delete A_tdt;

        // Create new vectors
        U_t = new Vector(size);
        V_t = new Vector(size);
        A_t = new Vector(size);
        U_tpdt = new Vector(size);
        V_tpdt = new Vector(size);
        V_fake = new Vector(size);
        A_tpdt = new Vector(size);
        U_tdt = new Vector(size);
        V_tdt = new Vector(size);
        A_tdt = new Vector(size);

        // Verify allocation was successful
        if (!U_t || !V_t || !A_t || !U_tpdt || !V_tpdt || !V_fake ||
            !A_tpdt || !U_tdt || !V_tdt || !A_tdt) {
            opserr << "ExplicitBathe::domainChanged - out of memory\n";
            return -1;
        }
    }

    // Initialize state vectors from committed DOF values
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    
    while ((dofPtr = theDOFs()) != nullptr) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        // Initialize displacements
        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*U_t)(loc) = disp(i);
                (*U_tpdt)(loc) = disp(i);
                (*U_tdt)(loc) = disp(i);
            }
        }

        // Initialize velocities
        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*V_t)(loc) = vel(i);
            }
        }

        // Initialize accelerations
        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*A_t)(loc) = accel(i);
            }
        }
    }

    // Initialize critical timestep values
    damped_minimum_critical_timestep = std::numeric_limits<double>::infinity();
    undamped_minimum_critical_timestep = std::numeric_limits<double>::infinity();

    // Reset computation flag if it was already computed
    if (compute_critical_timestep == 2) {
        compute_critical_timestep = 1;
    }

    return 0;
}

// Compute critical time step for all elements
// 
// This method computes the critical time step for each element by solving
// the generalized eigenvalue problem: K*v = lambda*M*v
// 
// The critical time step is then: dt_crit = 2/omega_max
// For damped systems: dt_crit = 2/(omega_max * (sqrt(1 + xi^2) - xi))
void computeCriticalTimestep(AnalysisModel *theModel, 
                             double &damped_min_dt, 
                             double &undamped_min_dt,
                             int &damped_elem_tag,
                             int &undamped_elem_tag)
{
    Domain* theDomain = theModel->getDomainPtr();
    Element *ele;
    ElementIter &elements = theDomain->getElements();
    
    while ((ele = elements()) != 0) {
        const Matrix &M = ele->getMass();
        const Matrix &K = ele->getInitialStiff();

        int n = M.noRows();
        if (n == 0 || K.noRows() != n) {
            continue;  // Skip elements without proper matrices
        }

        // Create lumped mass matrix (row-sum lumping)
        Matrix Mlumped(n, n);
        Mlumped.Zero();
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                sum += M(i, j);
            }
            Mlumped(i, i) = sum;
        }

        // Prepare matrices for LAPACK (column-major format)
        double *K_data = new double[n * n];
        double *M_data = new double[n * n];
        
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                K_data[j * n + i] = K(i, j);
                M_data[j * n + i] = Mlumped(i, i == j ? i : 0);
            }
        }

        // Eigenvalue problem arrays
        double *alphar = new double[n];
        double *alphai = new double[n];
        double *beta = new double[n];
        double *vl = nullptr;
        double *vr = nullptr;
        
        char jobvl = 'N';  // Don't compute left eigenvectors
        char jobvr = 'N';  // Don't compute right eigenvectors
        int lda = n, ldb = n;
        int info;
        
        // Workspace query
        int lwork = -1;
        double wkopt;
        
#ifdef _WIN32
        DGGEV(&jobvl, &jobvr, &n, K_data, &lda, M_data, &ldb,
              alphar, alphai, beta, vl, &lda, vr, &ldb, &wkopt, &lwork, &info);
#else
        dggev_(&jobvl, &jobvr, &n, K_data, &lda, M_data, &ldb,
               alphar, alphai, beta, vl, &lda, vr, &ldb, &wkopt, &lwork, &info);
#endif
        
        lwork = (int)wkopt;
        double *work = new double[lwork];
        
        // Actual computation
#ifdef _WIN32
        DGGEV(&jobvl, &jobvr, &n, K_data, &lda, M_data, &ldb,
              alphar, alphai, beta, vl, &lda, vr, &ldb, work, &lwork, &info);
#else
        dggev_(&jobvl, &jobvr, &n, K_data, &lda, M_data, &ldb,
               alphar, alphai, beta, vl, &lda, vr, &ldb, work, &lwork, &info);
#endif

        if (info > 0) {
            opserr << "WARNING: Eigenvalue computation failed for element " 
                   << ele->getTag() << "\n";
        }

        // Find maximum eigenvalue
        double maxEigenvalue = 0.0;
        for (int i = 0; i < n; ++i) {
            if (beta[i] != 0.0) {
                double lambda = alphar[i] / beta[i];
                if (lambda > maxEigenvalue) {
                    maxEigenvalue = lambda;
                }
            }
        }

        // Compute critical timesteps
        if (maxEigenvalue > 0.0) {
            double w_max = std::sqrt(maxEigenvalue);
            
            // Get Rayleigh damping coefficients
            Vector coefs = ele->getRayleighDampingFactors();
            double alphaM = coefs(0);
            double betaK = coefs(1);
            
            // Compute damping ratio
            double xi = 0.5 * (alphaM / w_max + betaK * w_max);

            // Critical timesteps
            double undamped_dt = 2.0 / w_max;
            double damped_dt = 2.0 / w_max * (std::sqrt(1.0 + xi*xi) - xi);

            // Update minimums
            if (damped_dt < damped_min_dt) {
                damped_min_dt = damped_dt;
                damped_elem_tag = ele->getTag();
            }
            if (undamped_dt < undamped_min_dt) {
                undamped_min_dt = undamped_dt;
                undamped_elem_tag = ele->getTag();
            }
        }

        // Clean up
        delete[] K_data;
        delete[] M_data;
        delete[] alphar;
        delete[] alphai;
        delete[] beta;
        delete[] work;
    }
}

// Advance to a new time step
//
// This method implements the first sub-step of the Bathe scheme:
// From time t to time t + p*dt
// 
// u_{t+p*dt} = u_t + p*dt*v_t + (p*dt)^2/2 * a_t
// v_{t+p*dt} = v_t + p*dt*a_t
// 
// Then solve: M*a_{t+p*dt} = R_{t+p*dt} - C*v_{t+p*dt} - K*u_{t+p*dt}
int ExplicitBathe::newStep(double _deltaT) {
    deltaT = _deltaT;

    if (!U_t || !V_t || !A_t) {
        opserr << "ExplicitBathe::newStep() - state variables not initialized\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();

    // Compute critical timestep if requested (only once)
    if (compute_critical_timestep == 1) {
        computeCriticalTimestep(theModel, 
                                damped_minimum_critical_timestep,
                                undamped_minimum_critical_timestep,
                                damped_critical_element_tag,
                                undamped_critical_element_tag);
        
        compute_critical_timestep = 2;  // Mark as computed
        
        opserr << "ExplicitBathe: Critical timestep analysis\n";
        opserr << "  UNDAMPED: dt_crit = " << undamped_minimum_critical_timestep 
               << " @ element #" << undamped_critical_element_tag << "\n";
        opserr << "  DAMPED:   dt_crit = " << damped_minimum_critical_timestep 
               << " @ element #" << damped_critical_element_tag << "\n";
    }

    // Report timestep stability factor
    if (compute_critical_timestep > 0) {
        double dT_factor = deltaT / damped_minimum_critical_timestep;
        opserr << "ExplicitBathe::newStep() - dt = " << deltaT 
               << ", dt_crit = " << damped_minimum_critical_timestep 
               << ", factor = " << dT_factor;
        if (dT_factor < 1.0) {
            opserr << " [OK]";
        } else {
            opserr << " [WARNING: dt > dt_crit!]";
        }
        opserr << endln;
    }

    // Compute integration coefficients for this time step
    a0 = p * deltaT;
    a1 = std::pow(p * deltaT, 2) / 2.0;
    a2 = a0 / 2.0;
    a3 = (1.0 - p) * deltaT;
    a4 = std::pow((1.0 - p) * deltaT, 2) / 2.0;
    a5 = q0 * a3;
    a6 = (0.5 + q1) * a3;
    a7 = q2 * a3;

    // Predict displacement and velocity at t + p*dt
    *U_tpdt = *U_t;
    U_tpdt->addVector(1.0, *V_t, a0);      // += v_t * p*dt
    U_tpdt->addVector(1.0, *A_t, a1);      // += a_t * (p*dt)^2/2
    
    *V_fake = *V_t;
    V_fake->addVector(1.0, *A_t, a0);      // += a_t * p*dt
    
    A_tpdt->Zero();  // Will be computed by solving system

    // Set response in model
    theModel->setResponse(*U_tpdt, *V_fake, *A_tpdt);

    // Update domain time
    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime + p * deltaT;

    if (theModel->updateDomain(newtime, p * deltaT) < 0) {
        opserr << "ExplicitBathe::newStep() - failed to update the domain\n";
        return -3;
    }

    return 0;
}

// Update the response quantities
//
// This method is called twice per time step:
// 1st call: After solving at t + p*dt, updates to prepare for t + dt
// 2nd call: After solving at t + dt, finalizes the time step
int ExplicitBathe::update(const Vector &U) {
    updateCount++;
    if (updateCount > 2) {
        opserr << "WARNING ExplicitBathe::update() - called more than twice\n";
        opserr << "  ExplicitBathe requires exactly 2 updates per time step\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING ExplicitBathe::update() - no AnalysisModel set\n";
        return -2;
    }

    if (U_t == 0) {
        opserr << "WARNING ExplicitBathe::update() - domainChange() failed or not called\n";
        return -3;
    }

    if (U.Size() != A_t->Size()) {
        opserr << "WARNING ExplicitBathe::update() - Vector size mismatch\n";
        return -4;
    }

    LinearSOE *theLinSOE = this->getLinearSOE();
    
    // Store acceleration at t + p*dt
    *A_tpdt = U;

    // Update velocity at t + p*dt (corrected)
    // v_{t+p*dt} = v_t + (a_t + a_{t+p*dt}) * p*dt/2
    *V_tpdt = *V_t;
    V_tpdt->addVector(1.0, *A_t, a2);      // += a_t * p*dt/2
    V_tpdt->addVector(1.0, *A_tpdt, a2);   // += a_{t+p*dt} * p*dt/2

    // Prepare fake velocity for response setting
    *V_fake = *V_tpdt;
    V_fake->addVector(1.0, *A_tpdt, a3);   // += a_{t+p*dt} * (1-p)*dt

    // Predict displacement at t + dt
    // u_{t+dt} = u_{t+p*dt} + v_{t+p*dt} * (1-p)*dt + a_{t+p*dt} * ((1-p)*dt)^2/2
    *U_tdt = *U_tpdt;
    U_tdt->addVector(1.0, *V_tpdt, a3);    // += v_{t+p*dt} * (1-p)*dt
    U_tdt->addVector(1.0, *A_tpdt, a4);    // += a_{t+p*dt} * ((1-p)*dt)^2/2

    A_tdt->Zero();  // Will be computed by solving system

    // Set response in model
    theModel->setResponse(*U_tdt, *V_fake, *A_tdt);

    // Update domain time
    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime + (1.0 - p) * deltaT;

    if (theModel->updateDomain(newtime, (1.0 - p) * deltaT) < 0) {
        opserr << "ExplicitBathe::update() - failed to update the domain\n";
        return -3;
    }

    // Solve for acceleration at t + dt
    this->formUnbalance();
    theLinSOE->solve();
    *A_tdt = theLinSOE->getX();

    // Report maximum acceleration for monitoring
    double A_max = A_tdt->pNorm(0);
    opserr << "ExplicitBathe::update() - max acceleration = " << A_max << endln;

    // Final velocity update at t + dt using all three accelerations
    // v_{t+dt} = v_{t+p*dt} + q0*a_t*(1-p)*dt + (0.5+q1)*a_{t+p*dt}*(1-p)*dt + q2*a_{t+dt}*(1-p)*dt
    *V_tdt = *V_tpdt;
    V_tdt->addVector(1.0, *A_t, a5);       // += q0 * a_t * (1-p)*dt
    V_tdt->addVector(1.0, *A_tpdt, a6);    // += (0.5+q1) * a_{t+p*dt} * (1-p)*dt
    V_tdt->addVector(1.0, *A_tdt, a7);     // += q2 * a_{t+dt} * (1-p)*dt

    // Set final response
    theModel->setResponse(*U_tdt, *V_tdt, *A_tdt);
    
    if (theModel->updateDomain() < 0) {
        opserr << "ExplicitBathe::update() - failed to update the domain\n";
        return -4;
    }

    return 0;
}

// Form element tangent matrix (mass matrix for explicit scheme)
int ExplicitBathe::formEleTangent(FE_Element *theEle) {
    theEle->zeroTangent();
    theEle->addMtoTang();
    return 0;
}

// Form nodal tangent matrix (mass matrix for explicit scheme)
int ExplicitBathe::formNodTangent(DOF_Group *theDof) {
    theDof->zeroTangent();
    theDof->addMtoTang();
    return 0;
}

// Commit the state for this time step
//
// This method is called after both sub-steps are complete.
// It updates the state vectors for the next time step.
int ExplicitBathe::commit() {
    updateCount = 0;  // Reset update counter for next step

    // Update state vectors: t+dt becomes new t
    *U_t = *U_tdt;
    *V_t = *V_tdt;
    *A_t = *A_tdt;

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == nullptr) {
        opserr << "ExplicitBathe::commit() - no AnalysisModel set\n";
        return -1;
    }

    return theModel->commitDomain();
}

// Get current velocity (for modal damping interface)
const Vector &ExplicitBathe::getVel() {
    return *V_t;
}

// Send object state for parallel processing
int ExplicitBathe::sendSelf(int cTag, Channel &theChannel) {
    Vector data(1);
    data(0) = p;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}

// Receive object state for parallel processing
int ExplicitBathe::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    Vector data(1);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::recvSelf() - could not receive data\n";
        return -1;
    }

    p = data(0);

    // Recalculate integration coefficients from received p
    q1 = (1.0 - 2.0*p) / (2.0*p*(1.0 - p));
    q2 = 0.5 - p * q1;
    q0 = -q1 - q2 + 0.5;

    return 0;
}

// Print integrator information
void ExplicitBathe::Print(OPS_Stream &stream, int flag) {
    stream << "Explicit Bathe Method\n";
    stream << "  Time Step: " << deltaT << "\n";
    stream << "  Damping parameter p: " << p << "\n";
    stream << "  Integration coefficients: q0 = " << q0 
           << ", q1 = " << q1 << ", q2 = " << q2 << "\n";
    if (compute_critical_timestep > 0) {
        stream << "  Critical timestep (damped): " << damped_minimum_critical_timestep << "\n";
        stream << "  Critical timestep (undamped): " << undamped_minimum_critical_timestep << "\n";
    }
}
