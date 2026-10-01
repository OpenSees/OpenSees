#include <ExplicitBatheLNVD.h>
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
#define OPS_Export

#include <NodeIter.h>
#include <LoadPatternIter.h>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>

#include <limits>

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

void *OPS_ExplicitBatheLNVD(void) {
    // Pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = nullptr;

    // Check if enough arguments are provided
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 2) {
        opserr << "WARNING: Insufficient arguments for ExplicitBatheLNVD integrator. Expected at least 2 arguments (p, alpha_flac).\n";
        return nullptr;
    }

    // Read input parameters
    double p, alpha_flac;
    numArgs = 2;
    double input[2];
    if (OPS_GetDoubleInput(&numArgs, input) < 0) {
        opserr << "WARNING: Invalid input for ExplicitBatheLNVD integrator. p and alpha_flac parameters\n";
        return nullptr;
    }
    p = input[0];
    alpha_flac = input[1];

    int compute_critical_timestep = 0;
    if (OPS_GetNumRemainingInputArgs() > 0) {
        numArgs = 1;
        if (OPS_GetIntInput(&numArgs, &compute_critical_timestep) < 0) {
            opserr << "WARNING: Invalid input for ExplicitBatheLNVD integrator. compute_critical_timestep parameter\n";
            return nullptr;
        }
    }

    // Create the ExplicitBatheLNVD integrator with the provided parameters
    theIntegrator = new ExplicitBatheLNVD(p, alpha_flac, compute_critical_timestep);

    if (theIntegrator == nullptr) {
        opserr << "WARNING - out of memory creating ExplicitBatheLNVD integrator\n";
    }
    return theIntegrator;
}

ExplicitBatheLNVD::ExplicitBatheLNVD()
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBatheLNVD),
      deltaT(0.0), p(0.0), q0(0.0), q1(0.0), q2(0.0), alpha_flac(0.0),
      U_t(0), V_t(0), A_t(0),
      U_tpdt(0), V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0), V_tdt(0), A_tdt(0), F_damp(0), F_unbal(0), updateCount(0),
      a0(0.), a1(0.), a2(0.), a3(0.), a4(0.), a5(0.), a6(0.), a7(0.), compute_critical_timestep(0)
{}

ExplicitBatheLNVD::ExplicitBatheLNVD(double _p, double _alpha_flac, int compute_critical_timestep_)
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBatheLNVD),
      deltaT(0.0), p(_p), q0(0.0), q1(0.0), q2(0.0), alpha_flac(_alpha_flac),
      U_t(0), V_t(0), A_t(0),
      U_tpdt(0), V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0), V_tdt(0), A_tdt(0), F_damp(0), F_unbal(0), updateCount(0),
      a0(0.), a1(0.), a2(0.), a3(0.), a4(0.), a5(0.), a6(0.), a7(0.), compute_critical_timestep(compute_critical_timestep_)
{
    // Calculate the integration constants based on p
    q1 = (1 - 2*p)/(2*p*(1-p));
    q2 = 0.5 - p * q1;
    q0 = -q1 -q2 + 0.5;

    opserr << "ExplicitBatheLNVD - p = " << p << " alpha_flac = " << alpha_flac 
           << " compute_critical_timestep = " << compute_critical_timestep << endln;
}

ExplicitBatheLNVD::~ExplicitBatheLNVD() {
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
    if (F_damp) delete F_damp;
    if (F_unbal) delete F_unbal;
}

int ExplicitBatheLNVD::domainChanged() {
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();

    if (!theModel || !theLinSOE) {
        opserr << "ExplicitBatheLNVD::domainChanged - missing model or linear system\n";
        return -1;
    }

    const Vector &x = theLinSOE->getX();
    int size = x.Size();

    if (size == 0) {
        opserr << "ExplicitBatheLNVD::domainChanged - invalid size\n";
        return -1;
    }

    // Allocate memory for state variables
    if (!U_t || U_t->Size() != size) {
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
        if (F_damp) delete F_damp;
        if (F_unbal) delete F_unbal;

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
        F_damp = new Vector(size);
        F_unbal = new Vector(size);

        if (!U_t || !V_t || !A_t || !U_tpdt || !V_tpdt || !A_tpdt || 
            !U_tdt || !V_tdt || !A_tdt || !F_damp || !F_unbal) {
            opserr << "ExplicitBatheLNVD::domainChanged - out of memory\n";
            return -1;
        }
    }

    // Initialize state variables
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != nullptr) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        const Vector &disp = dofPtr->getCommittedDisp();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*U_t)(loc) = disp(i);
                (*U_tpdt)(loc) = disp(i);
                (*U_tdt)(loc) = disp(i);
            }
        }

        const Vector &vel = dofPtr->getCommittedVel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*V_t)(loc) = vel(i);
            }
        }

        const Vector &accel = dofPtr->getCommittedAccel();
        for (int i = 0; i < idSize; ++i) {
            int loc = id(i);
            if (loc >= 0) {
                (*A_t)(loc) = accel(i);
            }
        }
    }

    // Initialize damping and unbalanced force vectors
    F_damp->Zero();
    F_unbal->Zero();

    damped_minimum_critical_timestep = std::numeric_limits<double>::infinity();
    undamped_minimum_critical_timestep = std::numeric_limits<double>::infinity();

    if (compute_critical_timestep == 2) {
        compute_critical_timestep = 1;
    }

    return 0;
}

int ExplicitBatheLNVD::newStep(double _deltaT) {
    deltaT = _deltaT;

    // Ensure state variables are initialized
    if (!U_t || !V_t || !A_t) {
        opserr << "ExplicitBatheLNVD::newStep() - state variables not initialized\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();

    // Critical timestep computation (same as ExplicitBathe)
    if (compute_critical_timestep == 1) {
        Domain* theDomain = theModel->getDomainPtr();
        Element * ele;
        ElementIter &elements = theDomain->getElements();
        while ((ele = elements()) != 0) {
            const Matrix &M = ele->getMass();
            const Matrix &K = ele->getInitialStiff();

            int n = M.noRows();
            if (n == 0 || K.noRows() != n) {
                continue;
            }

            Matrix Mlumped(n, n);
            Mlumped.Zero();
            for (int i = 0; i < n; ++i) {
                double sum = 0.0;
                for (int j = 0; j < n; ++j) {
                    sum += M(i, j);
                }
                Mlumped(i, i) = sum;
            }

            // Convert matrices to column-major format for LAPACK
            double *M_data = new double[n * n];
            double *K_data = new double[n * n];
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    M_data[j * n + i] = Mlumped(i, j);
                    K_data[j * n + i] = K(i, j);
                }
            }

            // Perform generalized eigenvalue analysis
            char jobvl = 'N';
            char jobvr = 'N';
            double *alphar = new double[n];
            double *alphai = new double[n];
            double *beta = new double[n];
            double *vl = nullptr;
            double *vr = nullptr;
            int lda = n;
            int ldb = n;
            int info;
            int lwork = 8 * n;
            double *work = new double[lwork];

#ifdef _WIN32
            DGGEV(&jobvl, &jobvr, &n, K_data, &lda, M_data, &ldb, alphar, alphai, beta, vl, &lda, vr, &ldb, work, &lwork, &info);
#else
            dggev_(&jobvl, &jobvr, &n, K_data, &lda, M_data, &ldb, alphar, alphai, beta, vl, &lda, vr, &ldb, work, &lwork, &info);
#endif

            if (info > 0) {
                opserr << "WARNING: Eigenvalue computation did not converge for element " << ele->getTag() << "\n";
            }

            // Extract the largest eigenvalue
            double maxEigenvalue = 0.0;
            for (int i = 0; i < n; ++i) {
                if (beta[i] != 0) {
                    double lambda = alphar[i] / beta[i];
                    if (lambda > maxEigenvalue) {
                        maxEigenvalue = lambda;
                    }
                }
            }

            double w_max = std::sqrt(maxEigenvalue);
            double alphaM = 0., betaK = 0., betaK0 = 0., betaKc = 0.;
            Vector coefs = ele->getRayleighDampingFactors();
            alphaM = coefs(0);
            betaK = coefs(1);
            betaK0 = coefs(2);
            betaKc = coefs(3);

            double xi = 0.5*(alphaM / w_max + betaK*w_max);

            // Compute critical timestep
            double undamped_critical_timestep = 2.0 / w_max; 
            double damped_critical_timestep = 2.0 / w_max * (std::sqrt(1 + xi*xi) - xi);

            if (damped_minimum_critical_timestep > damped_critical_timestep) {
                damped_minimum_critical_timestep = damped_critical_timestep;
                damped_critical_element_tag = ele->getTag();
            }
            if (undamped_minimum_critical_timestep > undamped_critical_timestep) {
                undamped_minimum_critical_timestep = undamped_critical_timestep;
                undamped_critical_element_tag = ele->getTag();
            }

            // Clean up allocated memory
            delete[] M_data;
            delete[] K_data;
            delete[] alphar;
            delete[] alphai;
            delete[] beta;
            delete[] work;
        }
        compute_critical_timestep = 2;
        opserr << " Overall UNDAMPED Critical timestep = " << undamped_minimum_critical_timestep << " @ element # " << undamped_critical_element_tag << "\n";
        opserr << " Overall  DAMPED  Critical timestep = " << damped_minimum_critical_timestep << " @ element # " << damped_critical_element_tag << "\n";
    }

    if (compute_critical_timestep > 0) {
        double dT_factor = deltaT / damped_minimum_critical_timestep;
        opserr << "    ExplicitBatheLNVD::newStep()  dt =  " << deltaT << "   dt_crit = " <<  damped_minimum_critical_timestep 
               << " dT_factor = " << dT_factor 
               << (dT_factor < 1. ? " Ok!" : " WARNING! dt > dt_crit") << endln;
    }

    // A. Initial Calculations
    a0 = p * deltaT;
    a1 = std::pow(p * deltaT, 2) / 2;
    a2 = a0 / 2;
    a3 = (1 - p) * deltaT;
    a4 = std::pow((1 - p) * deltaT, 2) / 2;
    a5 = q0 * a3;
    a6 = (0.5 + q1) * a3;
    a7 = q2 * a3;

    // Prepare for first matrix inversion
    *U_tpdt = *U_t;
    U_tpdt->addVector(1.0, *V_t, a0);
    U_tpdt->addVector(1.0, *A_t, a1);
    *V_fake = *V_t;
    V_fake->addVector(1.0, *A_t, a0);
    A_tpdt->Zero();

    theModel->setResponse(*U_tpdt, *V_fake, *A_tpdt);



    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime + p*deltaT;

    if (theModel->updateDomain(newtime, p*deltaT) < 0) {
        opserr << "ExplicitBatheLNVD - failed to update the domain\n";
        return -3;
    }

    // // Form unbalanced forces and solve for acceleration
    // this->formUnbalance();
    
    // // Store unbalanced forces before solving
    // LinearSOE *theLinSOE = this->getLinearSOE();
    // *F_unbal = theLinSOE->getB();
    
    // // Compute FLAC damping forces
    // this->computeFLACDamping(V_fake);
    // // this->computeFLACDamping(V_t);
    
    // // Add damping forces to the system
    // F_unbal->addVector(1.0, *F_damp, 1.0);
    // theLinSOE->setB(*F_unbal);

    return 0;
}

int ExplicitBatheLNVD::update(const Vector &U) {
    updateCount++;
    if (updateCount > 2) {
        opserr << "WARNING ExplicitBatheLNVD::update() - called more than once -";
        opserr << " ExplicitBatheLNVD integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING ExplicitBatheLNVD::update() - no AnalysisModel set\n";
        return -2;
    }

    // check domainChanged() has been called
    if (U_t == 0) {
        opserr << "WARNING ExplicitBatheLNVD::update() - domainChange() failed or not called\n";
        return -3;
    }

    // check vectors are of correct size
    if (A_t->Size() != A_tdt->Size()) {
        opserr << "WARNING ExplicitBatheLNVD::update() - Vectors of incompatible size ";
        opserr << " expecting " << A_t->Size() << " obtained " << A_tdt->Size() << endln;
        return -4;
    }

    
    *A_tpdt = U;

    // determine the response at t + p*deltaT
    *V_tpdt = *V_t;
    V_tpdt->addVector(1.0, *A_t, a2);
    V_tpdt->addVector(1.0, *A_tpdt, a2);

    *V_fake = *V_tpdt;
    V_fake->addVector(1.0, *A_tpdt, a3);

    // Get displacement at t + deltaT
    *U_tdt = *U_tpdt;
    U_tdt->addVector(1.0, *V_tpdt, a3);
    U_tdt->addVector(1.0, *A_tpdt, a4);

    A_tdt->Zero();

    theModel->setResponse(*U_tdt, *V_fake, *A_tdt);

    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime + (1-p)*deltaT;

    if (theModel->updateDomain(newtime, (1-p)*deltaT) < 0) {
        opserr << "ExplicitBatheLNVD::update() - failed to update the domain\n";
        return -3;
    }

    // Form unbalanced forces and solve for acceleration
    
    // // Store unbalanced forces before solving
    LinearSOE *theLinSOE = this->getLinearSOE();
    // theLinSOE->zeroB();
    // *F_unbal = theLinSOE->getB();
    this->formUnbalance();
    
    // // Compute FLAC damping forces
    // this->computeFLACDamping(V_fake);
    
    // // Add damping forces to the system
    // F_unbal->addVector(1.0, *F_damp, 1.0);
    // theLinSOE->setB(*F_unbal);
    
    theLinSOE->solve();
    *A_tdt = theLinSOE->getX();

    double A_max = A_tdt->pNorm(0);
    opserr << "    ExplicitBatheLNVD::update()  A_max =  " << A_max << endln;

    *V_tdt = *V_tpdt;
    V_tdt->addVector(1.0, *A_t, a5);
    V_tdt->addVector(1.0, *A_tpdt, a6);
    V_tdt->addVector(1.0, *A_tdt, a7);

    // set response at t to be that at t+deltaT of previous step
    theModel->setResponse(*U_tdt, *V_tdt, *A_tdt);
    if (theModel->updateDomain() < 0) {
        opserr << "ExplicitBatheLNVD::update() - failed to update the domain\n";
        return -4;
    }

    return 0;
}

int 
ExplicitBatheLNVD::formNodalUnbalance(void)
{
    // opserr << "my formNodalUnbalance" << endln;

    // loop through the DOF_Groups and add the unbalance
    DOF_GrpIter &theDOFs = (this->getAnalysisModel())->getDOFs();
    DOF_Group *dofPtr;
    int res = 0;

    static Vector Fdamping(10);

    while ((dofPtr = theDOFs()) != 0) { 

        const Vector &F_unbalanced = dofPtr->getUnbalance(this);
        const Vector &Vtrial = dofPtr->getTrialVel();
        // const Vector &Vtrial = dofPtr->getCommittedVel();

        Fdamping.resize(F_unbalanced.Size());
        Fdamping.Zero();

        for (int i = 0; i < F_unbalanced.Size(); ++i)
        {
            double f_unbal_i = F_unbalanced(i);
            double v_i = Vtrial(i);

        
            // Compute damping force: alpha_flac * |f_unbal| * sign(v)
            double sign_v = (v_i > 0.0) ? 1.0 : ((v_i < 0.0) ? -1.0 : 0.0);
            // double sign_v = (v_i > 0.0) ? -1.0 : ((v_i < 0.0) ? 1.0 : 0.0);
            Fdamping(i) = -alpha_flac * std::abs(f_unbal_i) * sign_v;
            // F_unbalanced += Fdamping(i);
            Fdamping(i) += F_unbalanced(i);
        }

        LinearSOE *theLinSOE = this->getLinearSOE();
        if (theLinSOE->addB(Fdamping, dofPtr->getID()) <0) {
            opserr << "WARNING IncrementalIntegrator::formNodalUnbalance -";
            opserr << " failed in addB for ID " << dofPtr->getID();
            res = -2;
        }
    }
    
    return res;
}

void ExplicitBatheLNVD::computeFLACDamping(Vector * vel) {
    // FLAC local non-viscous damping: f_damp = alpha_flac * |f_unbal| * sign(v)
    F_damp->Zero();
    
    int size = F_unbal->Size();
    for (int i = 0; i < size; i++) {
        double f_unbal_i = (*F_unbal)(i);
        double v_i = (*vel)(i);

        
        // Compute damping force: alpha_flac * |f_unbal| * sign(v)
        double sign_v = (v_i > 0.0) ? 1.0 : ((v_i < 0.0) ? -1.0 : 0.0);
        // double sign_v = (v_i > 0.0) ? -1.0 : ((v_i < 0.0) ? 1.0 : 0.0);
        (*F_damp)(i) = -alpha_flac * std::abs(f_unbal_i) * sign_v;
    }
}

int ExplicitBatheLNVD::formEleTangent(FE_Element *theEle) {
    theEle->zeroTangent();
    theEle->addMtoTang();
    return 0;
}

int ExplicitBatheLNVD::formNodTangent(DOF_Group *theDof) {
    theDof->zeroTangent();
    theDof->addMtoTang();
    return 0;
}

int ExplicitBatheLNVD::commit() {
    updateCount = 0;

    *U_t = *U_tdt;
    *V_t = *V_tdt;
    *A_t = *A_tdt;

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == nullptr) {
        opserr << "ExplicitBatheLNVD::commit() - no AnalysisModel set\n";
        return -1;
    }

    return theModel->commitDomain();
}

const Vector &ExplicitBatheLNVD::getVel() {
    return *V_t;
}

int ExplicitBatheLNVD::sendSelf(int cTag, Channel &theChannel) {
    Vector data(2);
    data(0) = p;
    data(1) = alpha_flac;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheLNVD::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}

int ExplicitBatheLNVD::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    Vector data(2);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBatheLNVD::recvSelf() - could not receive data\n";
        return -1;
    }

    p = data(0);
    alpha_flac = data(1);

    // Recalculate integration constants based on received p
    q1 = (1 - 2*p)/(2*p*(1-p));
    q2 = 0.5 - p * q1;
    q0 = -q1 -q2 + 0.5;

    return 0;
}

void ExplicitBatheLNVD::Print(OPS_Stream &stream, int flag) {
    stream << "Explicit Bathe FLAC Method:\n";
    stream << "  Time Step: " << deltaT << "\n";
    stream << "  p: " << p << ", q0: " << q0 << ", q1: " << q1 << ", q2: " << q2 << "\n";
    stream << "  FLAC damping coefficient (alpha_flac): " << alpha_flac << "\n";
}