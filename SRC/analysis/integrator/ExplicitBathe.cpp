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

void *OPS_ExplicitBathe(void) {
    // Pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = nullptr;

    // Check if enough arguments are provided
    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 1) {
        opserr << "WARNING: Insufficient arguments for ExplicitBathe integrator. Expected 1 argument (p).\n";
        return nullptr;
    }

    // Read input parameters
    double p;
    numArgs = 1;
    if (OPS_GetDoubleInput(&numArgs, &p) < 0) {
        opserr << "WARNING: Invalid input for ExplicitBathe integrator. p parameter\n";
        return nullptr;
    }

    int compute_critical_timestep = 0;
    if (OPS_GetNumRemainingInputArgs() > 0)
    {
    	if (OPS_GetIntInput(&numArgs, &compute_critical_timestep) < 0) {
        	opserr << "WARNING: Invalid input for ExplicitBathe integrator. compute_critical_timestep parameter\n";
        	return nullptr;
    	}
    }


    // Create the ExplicitBathe integrator with the provided parameters
    theIntegrator = new ExplicitBathe(p, compute_critical_timestep);

    if (theIntegrator == nullptr) {
        opserr << "WARNING - out of memory creating ExplicitBathe integrator\n";
    }
    return theIntegrator;
}

ExplicitBathe::ExplicitBathe()
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(0.0), q0(0.0), q1(0.0), q2(0.0), 
      U_t(0),V_t(0),A_t(0),
      U_tpdt(0),V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0),V_tdt(0),A_tdt(0), updateCount(0),
      a0(0.),a1(0.),a2(0.),a3(0.),a4(0.),a5(0.),a6(0.),a7(0.),compute_critical_timestep(0)
{}

ExplicitBathe::ExplicitBathe(double _p, int compute_critical_timestep_)//);, double _q0, double _q1, double _q2, double _s)
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitBathe),
      deltaT(0.0), p(_p), q0(0.0), q1(0.0), q2(0.0), 
      U_t(0),V_t(0),A_t(0),
      U_tpdt(0),V_tpdt(0), V_fake(0), A_tpdt(0),
      U_tdt(0),V_tdt(0),A_tdt(0), updateCount(0),
      a0(0.),a1(0.),a2(0.),a3(0.),a4(0.),a5(0.),a6(0.),a7(0.),compute_critical_timestep(compute_critical_timestep_)
{
    // Calculate the integration constants based on p
    q1 = (1 - 2*p)/(2*p*(1-p));
    q2 = 0.5 - p * q1;
    q0 = -q1 -q2 + 0.5;
    // s = -1;

    opserr << "ExplicitBathe - @jaabell p =" << p << " compute_critical_timestep = " << compute_critical_timestep << endln;

}

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
}



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

        if (!U_t || 
            !V_t || 
            !A_t || 
            !U_tpdt || 
            !V_tpdt || 
            !A_tpdt || 
            !U_tdt || 
            !V_tdt || 
            !A_tdt 
            ) {
            opserr << "ExplicitBathe::domainChanged - out of memory\n";
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
                (*U_tpdt)(loc) = disp(i);  // Assume Ut + pdt = Ut initially
                (*U_tdt)(loc) = disp(i);  // Assume  Ut +  dt = Ut initially
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

    // *U_t = *U_tdt;
    // *V_t = *V_tdt;
    // *A_t = *A_tdt;
    damped_minimum_critical_timestep = std::numeric_limits<double>::infinity();
    undamped_minimum_critical_timestep = std::numeric_limits<double>::infinity();

    if (compute_critical_timestep ==2)
    {
    	compute_critical_timestep = 1;
    }

    return 0;
}




int ExplicitBathe::newStep(double _deltaT) {
    deltaT = _deltaT;


    // Ensure state variables are initialized
    if (!U_t || !V_t || !A_t) {
        opserr << "ExplicitBathe::newStep() - state variables not initialized\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();


    if(compute_critical_timestep == 1)
    {
    	Domain* theDomain = theModel->getDomainPtr();
    	Element * ele;
        ElementIter &elements = theDomain->getElements();
        while ((ele = elements()) != 0)
        {
            const Matrix &M = ele->getMass();
            const Matrix &K = ele->getInitialStiff();

            int n = M.noRows();
            if (n == 0 || K.noRows() != n) {
                continue;
            }

            Matrix Mlumped(n, n);
            Mlumped.Zero();
            for (int i = 0; i < n; ++i)
            {
                double sum = 0.0;
                for (int j = 0; j < n; ++j)
                {
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

            // Perform generalized eigenvalue analysis using DGGEV function from LAPACK
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

            // Extract the largest eigenvalue (critical timestep calculation)
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
            double alphaM=0.,  betaK=0.,  betaK0=0.,  betaKc=0.;
            Vector coefs = ele->getRayleighDampingFactors();
            alphaM = coefs(0);
            betaK = coefs(1);
            betaK0 = coefs(2);
            betaKc = coefs(3);

            double xi = 0.5*(alphaM / w_max + betaK*w_max);


            // Compute critical timestep using the largest eigenvalue
            double undamped_critical_timestep = 2.0 / w_max; 
            double damped_critical_timestep = 2.0 / w_max * (std::sqrt(1 + xi*xi) - xi); // We can use up to 1.9 times this value!
            // opserr << "Element " << ele->getTag() << ": Critical timestep = " << critical_timestep << "\n";

            if (damped_minimum_critical_timestep > damped_critical_timestep)
            {
            	damped_minimum_critical_timestep = damped_critical_timestep;
            	damped_critical_element_tag = ele->getTag();
            }
            if (undamped_minimum_critical_timestep > undamped_critical_timestep)
            {
            	undamped_minimum_critical_timestep = undamped_critical_timestep;
            	undamped_critical_element_tag = ele->getTag();
            }
            // minimum_critical_timestep = std::min(minimum_critical_timestep, critical_timestep);

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

    if (compute_critical_timestep > 0)
    {
        double dT_factor =deltaT / damped_minimum_critical_timestep;
        opserr << "    ExplicitBathe::newStep()  dt =  " << deltaT << "   dt_crit = " <<  damped_minimum_critical_timestep 
        << " dT_factor = " << dT_factor 
        << (dT_factor < 1. ? " Ok!" : " WARNING! dt > dt_crit") << endln;
    }




    // A. Initial Calculations
    a0 = p * deltaT;
    a1 = std::pow(p * deltaT,2)/2;
    a2 = a0/2;
    a3 = (1 - p) * deltaT;
    a4 = std::pow((1 - p) * deltaT,2)/2;
    a5 = q0 * a3;
    a6 = (0.5 + q1) * a3;
    a7 = q2 * a3;


    // Prepare for first matrix inversion, which occurs after newstep. 
    *U_tpdt = *U_t;// + a0 * (*V_t) + a1 * (*A_t);
    U_tpdt->addVector(1.0, *V_t, a0);
    U_tpdt->addVector(1.0, *A_t, a1);
    *V_fake = *V_t; // + a0 * (*A_t);
    V_fake->addVector(1.0, *A_t, a0);
    A_tpdt->Zero(); //*A_tpdt = *A_t;

    theModel->setResponse(*U_tpdt, *V_fake, *A_tpdt);

    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime + p*deltaT;
    // theModel->setCurrentDomainTime(newtime);

    if (theModel->updateDomain(newtime, p*deltaT) < 0)  {
        opserr << "ExplicitDifference - failed to update the domain\n";
        return -3;
    }

    // Now the SOE gets solved and we get update called


    return 0;
}




int ExplicitBathe::update(const Vector &U) {
  
   updateCount++;
    if (updateCount > 2)  {
        opserr << "WARNING ExplicitBathe::update() - called more than once -";
        opserr << " ExplicitBathe integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING ExplicitBathe::update() - no souAnalysisModel set\n";
        return -2;
    }

    // check domainChanged() has been called, i.e. Ut will not be zero
    if (U_t == 0)  {
        opserr << "WARNING ExplicitBathe::update() - domainChange() failed or not called\n";
        return -3;
    }

    // check Udotdot is of correct size
    if (A_t->Size() != A_tdt->Size()) {
        opserr << "WARNING ExplicitBathe::update() - Vectors of incompatible size ";
        opserr << " expecting " << A_t->Size() << " obtained " << A_tdt->Size() << endln;
        return -4;
    }

    int size = A_t->Size();

    LinearSOE *theLinSOE = this->getLinearSOE();
    *A_tpdt = U;//theLinSOE->getX();

    // determine the response at t  + p *deltaT
    *V_tpdt = *V_t;
    V_tpdt -> addVector(1.0, *A_t, a2);
    V_tpdt -> addVector(1.0, *A_tpdt, a2);

    *V_fake = *V_tpdt;
    V_fake -> addVector(1.0, *A_tpdt, a3);

    // Get displacement at t  + deltaT
    *U_tdt = *U_tpdt;
    U_tdt -> addVector(1.0, *V_tpdt, a3);
    U_tdt -> addVector(1.0, *A_tpdt, a4);

    A_tdt->Zero();

    theModel->setResponse(*U_tdt, *V_fake, *A_tdt);

    double oldtime = theModel->getCurrentDomainTime();
    double newtime = oldtime+(1-p)*deltaT;

    if (theModel->updateDomain(newtime, (1-p)*deltaT) < 0)  {
        opserr << "ExplicitDifference::newStep() - failed to update the domain\n";
        return -3;
    }

    this->formUnbalance();
    theLinSOE->solve();
    *A_tdt = theLinSOE->getX();

    double A_max = A_tdt->pNorm(0);

    opserr << "    ExplicitBathe::update()  A_max =  " << A_max << endln;


    *V_tdt = *V_tpdt;
    V_tdt->addVector(1.0, *A_t, a5);
    V_tdt->addVector(1.0, *A_tpdt, a6);
    V_tdt->addVector(1.0, *A_tdt, a7);

    // set response at t to be that at t+deltaT of previous step
    theModel->setResponse(*U_tdt,*V_tdt,*A_tdt);
    if (theModel->updateDomain() < 0)  {
        opserr << "Newmark::update() - failed to update the domain\n";
        return -4;
    }

    return 0;
}

int ExplicitBathe::formEleTangent(FE_Element *theEle) {

    theEle->zeroTangent();
    theEle->addMtoTang();

    return 0;
}

int ExplicitBathe::formNodTangent(DOF_Group *theDof) {
    theDof->zeroTangent();
    theDof->addMtoTang();

    return 0;
}

int ExplicitBathe::commit() {

    updateCount = 0;

    *U_t = *U_tdt;
    *V_t = *V_tdt;
    *A_t = *A_tdt;

    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == nullptr) {
        opserr << "ExplicitBathe::commit() - no AnalysisModel set\n";
        return -1;
    }

    // double time = theModel->getCurrentDomainTime();
    // time += deltaT;
    // theModel->setCurrentDomainTime(time);

    return theModel->commitDomain();
}

const Vector &ExplicitBathe::getVel() {
    return *V_t;
}

int ExplicitBathe::sendSelf(int cTag, Channel &theChannel) {
    Vector data(1);
    data(0) = p;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}

int ExplicitBathe::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
    Vector data(1);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
        opserr << "ExplicitBathe::recvSelf() - could not receive data\n";
        return -1;
    }

    p = data(0);

    // Recalculate integration constants based on received p
    q1 = (1 - 2*p)/(2*p*(1-p));
    q2 = 0.5 - p * q1;
    q0 = -q1 -q2 + 0.5;
    // s = -1;

    return 0;
}

void ExplicitBathe::Print(OPS_Stream &stream, int flag) {
    stream << "Explicit Bathe Method:\n";
    stream << "  Time Step: " << deltaT << "\n";
    stream << "  p: " << p << ", q0: " << q0 << ", q1: " << q1 << ", q2: " << q2 << "\n";
}
