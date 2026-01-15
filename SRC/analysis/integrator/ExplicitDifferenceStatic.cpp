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
#define OPS_Export 


void* OPS_ExplicitDifferenceStatic(void)
{
    TransientIntegrator *theIntegrator = 0;
    theIntegrator = new ExplicitDifferenceStatic();

    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating ExplicitDifferenceStatic integrator\n";

    return theIntegrator;
}


ExplicitDifferenceStatic::ExplicitDifferenceStatic()
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitDifferenceStatic),
    deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    updateCount(0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0),
    Udot(0), Utdotdot1(0), U(0), Utdot1(0), velSignMem(0), prevUnbal(0), vSignEps(ED_VSIGN_EPS)
{

}


ExplicitDifferenceStatic::ExplicitDifferenceStatic(
    double _alphaM, double _betaK, double _betaKi, double _betaKc)
    : TransientIntegrator(INTEGRATOR_TAGS_ExplicitDifferenceStatic),
    deltaT(0.0),
    alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    updateCount(0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0),
    Udot(0), Utdotdot1(0), U(0), Utdot1(0), velSignMem(0), prevUnbal(0), vSignEps(ED_VSIGN_EPS)
{
}


ExplicitDifferenceStatic::~ExplicitDifferenceStatic()
{
    // clean up the memory created

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


int ExplicitDifferenceStatic::newStep(double _deltaT)
{
    updateCount = 0;

    deltaT = _deltaT;

    if (deltaT <= 0.0)  {
        opserr << "ExplicitDifferenceStatic::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -1;
    }

    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();

    //calculate vel at t+0.5deltaT and U at t+delatT
    Utdot->addVector(1.0, *Utdotdot, deltaT);
    Ut->addVector(1.0, *Utdot, deltaT);

    int size = Utdotdot->Size();

    if (Ut == 0)  {
        opserr << "ExplicitDifferenceStatic::newStep() - domainChange() failed or hasn't been called\n";
        return -2;
    }

    // for leap-frog method Ma=f-ku-cv, on the right side there is no Ma
    (*Utdotdot) *= 0;

    // set the garbage response quantities for the nodes
    theModel->setVel(*Utdot);
    theModel->setAccel(*Utdotdot);
    theModel->setDisp(*Ut);

    // increment the time to t and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "ExplicitDifferenceStatic::newStep() - failed to update the domain\n";
        return -3;
    }

    // set response at t to be that at t+deltaT of previous step
    (*Utdotdot) = (*Utdotdot1);
    
    return 0;
}


int ExplicitDifferenceStatic::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();

    theEle->addMtoTang();

    return 0;
}


int ExplicitDifferenceStatic::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();

    theDof->addMtoTang();

    return(0);
}


int ExplicitDifferenceStatic::domainChanged()
{

    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();


    opserr << "ExplicitDifferenceStatic::domainChanged()" << endln;

    // if damping factors exist set them in the element & node of the domain
    if (alphaM != 0.0 || betaK != 0.0 || betaKi != 0.0 || betaKc != 0.0)
        theModel->setRayleighDampingFactors(alphaM, betaK, betaKi, betaKc);


    // create the new Vector objects
    if (Ut == 0 || Ut->Size() != size)  {

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

        // create the new

        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        Udot = new Vector(size);
        U = new Vector(size);
        Utdotdot1 = new Vector(size);
        Utdot1 = new Vector(size);
        velSignMem = new Vector(size);
        prevUnbal = new Vector(size);

        // check we obtained the new
        if ( Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            U == 0 || U->Size() != size ||
            Utdotdot1 == 0 || Utdotdot1->Size() != size ||
            Utdot1 == 0 || Utdot1->Size() != size  || velSignMem->Size() != size || prevUnbal->Size() != size
        )  {

            opserr << "ExplicitDifferenceStatic::domainChanged - ran out of memory\n";

            // delete the old
    
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

    // now go through and populate U, Udot and Udotdot by iterating through
    // the DOF_Groups and getting the last committed velocity and accel
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0)  {

        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        int i;
        const Vector &disp = dofPtr->getCommittedDisp();
        for (i = 0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {            
                (*Ut)(loc) = disp(i);
            }
        }

        const Vector &vel = dofPtr->getCommittedVel();
        for (i = 0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Utdot)(loc) = vel(i);
                (*Utdot1)(loc) = vel(i);
            }
        }

        const Vector &accel = dofPtr->getCommittedAccel();
        for (i = 0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Utdotdot)(loc) = accel(i);
                (*Utdotdot1)(loc) = accel(i);
            }
        }
    }

    velSignMem->Zero();   // no previous sign yet

    opserr << "WARNING: ExplicitDifferenceStatic::domainChanged() - assuming Ut-1 = Ut\n";

    return 0;
}


// int 
// ExplicitDifferenceStatic::formNodalUnbalance(void)
// {
//     // opserr << "my formNodalUnbalance" << endln;

//     const double alpha_flac = 0.8;

//     // loop through the DOF_Groups and add the unbalance
//     DOF_GrpIter &theDOFs = (this->getAnalysisModel())->getDOFs();
//     DOF_Group *dofPtr;
//     int res = 0;

//     static Vector Fdamping(10);

//     while ((dofPtr = theDOFs()) != 0) { 

//         const Vector &F_unbalanced = dofPtr->getUnbalance(this);
//         const Vector &Vtrial = dofPtr->getTrialVel();
//         // const Vector &Vtrial = dofPtr->getCommittedVel();

//         Fdamping.resize(F_unbalanced.Size());
//         Fdamping.Zero();

//         for (int i = 0; i < F_unbalanced.Size(); ++i)
//         {
//             double f_unbal_i = F_unbalanced(i);
//             double v_i = Vtrial(i);

        
//             // Compute damping force: alpha_flac * |f_unbal| * sign(v)
//             double sign_v = (v_i > 0.0) ? 1.0 : ((v_i < 0.0) ? -1.0 : 0.0);
//             // double sign_v = (v_i > 0.0) ? -1.0 : ((v_i < 0.0) ? 1.0 : 0.0);
//             Fdamping(i) = -alpha_flac * std::abs(f_unbal_i) * sign_v;
//             // F_unbalanced += Fdamping(i);
//             Fdamping(i) += F_unbalanced(i);
//         }

//         LinearSOE *theLinSOE = this->getLinearSOE();
//         if (theLinSOE->addB(Fdamping, dofPtr->getID()) <0) {
//             opserr << "WARNING IncrementalIntegrator::formNodalUnbalance -";
//             opserr << " failed in addB for ID " << dofPtr->getID();
//             res = -2;
//         }
//     }
    
//     return res;
// }

int
ExplicitDifferenceStatic::formNodalUnbalance(void)
{
    const double alpha_flac = 0.59;
    const bool useCombined = true; 
    
    DOF_GrpIter &theDOFs = (this->getAnalysisModel())->getDOFs();
    DOF_Group *dofPtr;
    int res = 0;

    static Vector Fdamping(10);

    while ((dofPtr = theDOFs()) != 0) {

        const Vector &F_unbalanced = dofPtr->getUnbalance(this);
        const Vector &Vtrial       = dofPtr->getTrialVel();
        const ID     &id           = dofPtr->getID();     // global equation numbers  // NEW

        Fdamping.resize(F_unbalanced.Size());
        Fdamping.Zero();

        for (int i = 0; i < F_unbalanced.Size(); ++i) {

            // const double f_unbal_i = F_unbalanced(i);
            // const double v_i       = Vtrial(i);

            // // --- velocity sign memory (per global equation) ------------------ // NEW
            // double s = 0.0;
            // const int eq = (i < id.Size()) ? id(i) : -1;  // global equation index
            // if (eq >= 0 && eq < velSignMem->Size()) {
            //     s = (*velSignMem)(eq);                       // last stored sign (-1,0,+1)
            //     if (std::abs(v_i) > vSignEps) {
            //         s = (v_i > 0.0) ? 1.0 : -1.0;        // update sign if outside deadband
            //         (*velSignMem)(eq) = s;                   // persist it
            //     }
            // } else {
            //     // constrained or out-of-range DOF → no damping (s = 0)
            // }

            // // Local non-viscous damping:  Fd = -alpha * |Funbal| * sign(v_mem)
            // const double Fd_i = -alpha_flac * std::abs(f_unbal_i) * s;

            // // Accumulate RHS contribution for this DOF_Group entry:
            // Fdamping(i) = f_unbal_i + Fd_i;              // i.e., Funbal - alpha*|Funbal|*sign(v)

            const int eq = id(i);
            double v = Vtrial(i);
            double F = F_unbalanced(i);

            // --- velocity sign with memory + deadband
            double s_v = (*velSignMem)(eq);
            if (std::abs(v) > vSignEps) {
                s_v = (v > 0.0) ? 1.0 : -1.0;
                (*velSignMem)(eq) = s_v;
            }

            // --- "combined" part: sign of dF/dt ~ sign(ΔF)
            double dF = F - (*prevUnbal)(eq);
            double s_fdot = (dF > 0.0) ? 1.0 : ((dF < 0.0) ? -1.0 : 0.0);

            // --- choose damping
            double Fd;
            if (useCombined) {
                Fd = 0.5 * alpha_flac * std::abs(F) * (s_fdot - s_v);
            } else {
                Fd = -alpha_flac * std::abs(F) * s_v;
            }
            (*prevUnbal)(eq) = F;

            // add to RHS: Funbal + Fd
            Fdamping(i) = F + Fd;
        }

        LinearSOE *theLinSOE = this->getLinearSOE();
        if (theLinSOE->addB(Fdamping, id) < 0) {
            opserr << "WARNING ExplicitDifferenceStatic::formNodalUnbalance -"
                   << " failed in addB for ID " << id;
            res = -2;
        }
    }

    return res;
}

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
        opserr << "WARNING ExplicitDifferenceStatic::update() - no souAnalysisModel set\n";
        return -2;
    }

    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING ExplicitDifferenceStatic::update() - domainChange() failed or not called\n";
        return -3;
    }

    // check Udotdot is of correct size
    if (Udotdot.Size() != Utdotdot->Size()) {
        opserr << "WARNING ExplicitDifferenceStatic::update() - Vectors of incompatible size ";
        opserr << " expecting " << Utdotdot->Size() << " obtained " << Udotdot.Size() << endln;
        return -4;
    }

    int size = Udotdot.Size();


    // determine the response at t+deltaT
    double halfT = deltaT *0.125;

    Utdotdot1->addVector(0.0, Udotdot, 3.0);
    Utdotdot1->addVector(1.0, *Utdotdot, 1.0);

    //Velosity to output, because Utdot is velosity is defined at t+0.5deltaT
    Utdot1->addVector(0.0, *Utdot, 1.0);
    Utdot1->addVector(1.0, *Utdotdot1, halfT);


    theModel->setResponse(*Ut, *Utdot1, Udotdot);

    if (theModel->updateDomain() < 0)  {
        opserr << "ExplicitDifferenceStatic::update() - failed to update the domain\n";
        return -5;
    }



    // set response at t to be that at t+deltaT of previous step

    (*Utdotdot) = Udotdot;
    (*Utdotdot1) = Udotdot;




    return 0;
}


int ExplicitDifferenceStatic::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0) {
        opserr << "WARNING ExplicitDifferenceStatic::commit() - no AnalysisModel set\n";
        return -1;
    }

    return theModel->commitDomain();
}


int ExplicitDifferenceStatic::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = alphaM;
    data(1) = betaK;
    data(2) = betaKi;
    data(3) = betaKc;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING ExplicitDifferenceStatic::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int ExplicitDifferenceStatic::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING ExplicitDifferenceStatic::recvSelf() - could not receive data\n";
        return -1;
    }

    alphaM = data(0);
    betaK = data(1);
    betaKi = data(2);
    betaKc = data(3);

    return 0;
}


void ExplicitDifferenceStatic::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "ExplicitDifferenceStatic - currentTime: " << currentTime << endln;
        s << "  Rayleigh Damping - alphaM: " << alphaM << "  betaK: " << betaK;
        s << "  betaKi: " << betaKi << "  betaKc: " << betaKc << endln;
    }
    else
        s << "ExplicitDifferenceStatic - no associated AnalysisModel\n";
}



//a interface to get velosity for modal damping
const Vector &
ExplicitDifferenceStatic::getVel()
{
    return *Utdot;
}

