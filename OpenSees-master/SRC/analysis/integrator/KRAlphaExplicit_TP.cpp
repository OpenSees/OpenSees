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

#include <KRAlphaExplicit_TP.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include "fullGEN/FullGenLinSOE.h"
#include "fullGEN/FullGenLinLapackSolver.h"
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#define OPS_Export


void *    OPS_KRAlphaExplicit_TP(void)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 1) {
        opserr << "WARNING - incorrect number of args want KRAlphaExplicit_TP $rhoInf\n";
        return 0;
    }
    
    double rhoInf;
    if (OPS_GetDouble(&argc, &rhoInf) != 0) {
        opserr << "WARNING - invalid args want KRAlphaExplicit_TP $rhoInf\n";
        return 0;
    }
    
    theIntegrator = new KRAlphaExplicit_TP(rhoInf);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating KRAlphaExplicit_TP integrator\n";
    
    return theIntegrator;
}


KRAlphaExplicit_TP::KRAlphaExplicit_TP()
    : TransientIntegrator(INTEGRATOR_TAGS_KRAlphaExplicit_TP),
    alphaI(0.5), alphaF(0.5), beta(0.25), gamma(0.5),
    deltaT(0.0), alpha1(0), alpha3(0), Mhat(0),
    updateCount(0), initAlphaMatrices(1),
    c1(0.0), c2(0.0), c3(0.0),
    alphaM(0.0), alphaD(0.5), alphaR(0.5), alphaP(0.5),
    Ut(0), Utdot(0), Utdotdot(0),
    U(0), Udot(0), Udotdot(0),
    Utdothat(0), Put(0)
{
    
}


KRAlphaExplicit_TP::KRAlphaExplicit_TP(double _rhoInf)
    : TransientIntegrator(INTEGRATOR_TAGS_KRAlphaExplicit_TP),
    alphaI((2.0-_rhoInf)/(1.0+_rhoInf)), alphaF(1.0/(1.0+_rhoInf)),
    beta(1.0/(1.0+_rhoInf)/(1.0+_rhoInf)), gamma(0.5*(3.0-_rhoInf)/(1.0+_rhoInf)),
    deltaT(0.0), alpha1(0), alpha3(0), Mhat(0),
    updateCount(0),  initAlphaMatrices(1),
    c1(0.0), c2(0.0), c3(0.0),
    alphaM(0.0), alphaD(alphaF), alphaR(alphaF), alphaP(alphaF),
    Ut(0), Utdot(0), Utdotdot(0),
    U(0), Udot(0), Udotdot(0),
    Utdothat(0), Put(0)
{
    
}


KRAlphaExplicit_TP::~KRAlphaExplicit_TP()
{
    // clean up the memory created
    if (alpha1 != 0)
        delete alpha1;
    if (alpha3 != 0)
        delete alpha3;
    if (Mhat != 0)
        delete Mhat;
    if (Ut != 0)
        delete Ut;
    if (Utdot != 0)
        delete Utdot;
    if (Utdotdot != 0)
        delete Utdotdot;
    if (U != 0)
        delete U;
    if (Udot != 0)
        delete Udot;
    if (Udotdot != 0)
        delete Udotdot;
    if (Utdothat != 0)
        delete Utdothat;
    if (Put != 0)
        delete Put;
}


int KRAlphaExplicit_TP::newStep(double _deltaT)
{
    updateCount = 0;
    
    if (beta == 0 || gamma == 0 )  {
        opserr << "WARNING KRAlphaExplicit_TP::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::newStep() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -2;
    }
    
    // initialize the integration matrices
    if (initAlphaMatrices || _deltaT != deltaT)  {
        
        // update time step increment
        deltaT = _deltaT;
        if (deltaT <= 0.0)  {
            opserr << "WARNING KRAlphaExplicit_TP::newStep() - error in variable\n";
            opserr << "dT = " << deltaT << endln;
            return -3;
        }
        
        // get the ConvergenceTest so we can switch back later
        ConvergenceTest *theTest = this->getConvergenceTest();
        
        // set up the FullLinearSOE (needed to compute the alpha matrices)
        int size = theLinSOE->getNumEqn();
        FullGenLinSolver *theFullLinSolver = new FullGenLinLapackSolver();
        LinearSOE *theFullLinSOE = new FullGenLinSOE(size, *theFullLinSolver);
        if (theFullLinSOE == 0)  {
            opserr << "WARNING KRAlphaExplicit_TP::newStep() - failed to create FullLinearSOE\n";
            return -4;
        }
        theFullLinSOE->setLinks(*theModel);
        
        // now switch the SOE to the FullLinearSOE
        this->IncrementalIntegrator::setLinks(*theModel, *theFullLinSOE, theTest);
        
        // get a pointer to the A matrix of the FullLinearSOE
        const Matrix *tmp = theFullLinSOE->getA();
        if (tmp == 0)  {
            opserr << "WARNING KRAlphaExplicit_TP::newStep() - ";
            opserr << "failed to get A matrix of FullGeneral LinearSOE\n";
            return -5;
        }
        
        // calculate the integration parameter matrices
        c1 = beta*deltaT*deltaT;
        c2 = gamma*deltaT;
        c3 = 1.0;
        this->TransientIntegrator::formTangent(INITIAL_TANGENT);
        Matrix A(*tmp);
        
        c1 *= (1.0 - alphaF);
        c2 *= (1.0 - alphaF);
        c3 = (1.0 - alphaI);
        this->TransientIntegrator::formTangent(INITIAL_TANGENT);
        Matrix B3(*tmp);
        
        // solve [M + gamma*deltaT*C + beta*deltaT^2*K]*[alpha3] = 
        // [alphaI*M + alphaF*gamma*deltaT*C + alphaF*beta*deltaT^2*K] for alpha3
        A.Solve(B3, *alpha3);
        
        c1 = 0.0;
        c2 = 0.0;
        c3 = 1.0;
        this->TransientIntegrator::formTangent(INITIAL_TANGENT);
        Matrix B1(*tmp);
        
        // solve [M + gamma*deltaT*C + beta*deltaT^2*K]*[alpha1] = [M] for alpha1
        A.Solve(B1, *alpha1);
        
        // calculate the effective mass matrix Mhat
        Mhat->addMatrix(0.0, B1, 1.0);
        Mhat->addMatrixProduct(1.0, B1, *alpha3, -1.0);
        
        // switch the SOE back to the user specified one
        this->IncrementalIntegrator::setLinks(*theModel, *theLinSOE, theTest);
        
        // get initial unbalance at time t (in case domain changed)
        // warning: this will use committed stiffness prop. damping
        // from current step instead of previous step
        (*Utdotdot) = *Udotdot;
        alphaM = 1.0;
        alphaD = alphaR = alphaP = (1.0 - alphaF);
        Udotdot->addMatrixVector(0.0, *alpha3, *Utdotdot, 1.0);
        theModel->setAccel(*Udotdot);
        this->TransientIntegrator::formUnbalance();
        (*Put) = theLinSOE->getB();
        // reset accelerations at t+deltaT
        (*Udotdot) = *Utdotdot;
        theModel->setAccel(*Udotdot);
        
        // set flag that initializations are done
        initAlphaMatrices = 0;
    }
    
    if (U == 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::newStep() - domainChange() failed or hasn't been called\n";
        return -6;
    }
    
    // set weighting factors for subsequent iterations
    alphaM = 0.0;
    alphaD = alphaR = alphaP = alphaF;
    
    // determine new response at time t+deltaT
    Utdothat->addMatrixVector(0.0, *alpha1, *Utdotdot, deltaT);
    
    U->addVector(1.0, *Utdot, deltaT);
    double a1 = (0.5 + gamma)*deltaT;
    U->addVector(1.0, *Utdothat, a1);
    
    Udot->addVector(1.0, *Utdothat, 1.0);
    
    // no need to zero accelerations because alphaM = 0
    //Udotdot->Zero();
    //theModel->setResponse(*U, *Udot, *Udotdot);
    
    // set the trial response quantities
    theModel->setDisp(*U);
    theModel->setVel(*Udot);
    
    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::newStep() - failed to update the domain\n";
        return -7;
    }
    
    return 0;
}


int KRAlphaExplicit_TP::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int KRAlphaExplicit_TP::formTangent(int statFlag)
{
    statusFlag = statFlag;
    
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::formTangent() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -1;
    }
    
    theLinSOE->zeroA();
    
    int size = theLinSOE->getNumEqn();
    ID id(size);
    for (int i=1; i<size; i++)  {
        id(i) = id(i-1) + 1;
    }
    if (theLinSOE->addA(*Mhat, id) < 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::formTangent() - ";
        opserr << "failed to add Mhat to A\n";
        return -2;
    }
    
    return 0;
}


int KRAlphaExplicit_TP::formUnbalance()
{
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::formUnbalance() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -1;
    }
    
    theLinSOE->setB(*Put);
    
    // do modal damping
    const Vector *modalValues = theModel->getModalDampingFactors();
    if (modalValues != 0)  {
        this->addModalDampingForce(modalValues);
    }
    
    if (this->formElementResidual() < 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::formUnbalance() ";
        opserr << " - this->formElementResidual failed\n";
        return -2;
    }
    
    if (this->formNodalUnbalance() < 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::formUnbalance() ";
        opserr << " - this->formNodalUnbalance failed\n";
        return -3;
    }
    
    return 0;
}


int KRAlphaExplicit_TP::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    if (statusFlag == CURRENT_TANGENT)
        theEle->addKtToTang(c1);
    else if (statusFlag == INITIAL_TANGENT)
        theEle->addKiToTang(c1);
    
    theEle->addCtoTang(c2);
    theEle->addMtoTang(c3);
    
    return 0;
}


int KRAlphaExplicit_TP::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(c2);
    theDof->addMtoTang(c3);
    
    return 0;
}


int KRAlphaExplicit_TP::formEleResidual(FE_Element *theEle)
{
    theEle->zeroResidual();
    
    // this does not work because for some elements damping is returned
    // with the residual as well as the damping tangent 
    //theEle->addRtoResidual(alphaR);
    //theEle->addD_Force(*Udot, -alphaD);
    //theEle->addM_Force(*Udotdot, -alphaM);
    
    // instead use residual including the inertia terms and then correct
    // the mass contribution (only works because alphaR = alphaD) 
    theEle->addRIncInertiaToResidual(alphaR);
    theEle->addM_Force(*Udotdot, alphaR-alphaM);
    
    return 0;
}


int KRAlphaExplicit_TP::formNodUnbalance(DOF_Group *theDof)
{
    theDof->zeroUnbalance();
    
    theDof->addPtoUnbalance(alphaP);
    theDof->addD_Force(*Udot, -alphaD);
    theDof->addM_Force(*Udotdot, -alphaM);
    
    return 0;
}


int KRAlphaExplicit_TP::domainChanged()
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    
    // create the new Matrix and Vector objects
    if (Ut == 0 || Ut->Size() != size)  {
        
        // delete the old
        if (alpha1 != 0)
            delete alpha1;
        if (alpha3 != 0)
            delete alpha3;
        if (Mhat != 0)
            delete Mhat;
        if (Ut != 0)
            delete Ut;
        if (Utdot != 0)
            delete Utdot;
        if (Utdotdot != 0)
            delete Utdotdot;
        if (U != 0)
            delete U;
        if (Udot != 0)
            delete Udot;
        if (Udotdot != 0)
            delete Udotdot;
        if (Utdothat != 0)
            delete Utdothat;
        if (Put != 0)
            delete Put;
        
        // create the new
        alpha1 = new Matrix(size,size);
        alpha3 = new Matrix(size,size);
        Mhat = new Matrix(size,size);
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        Utdothat = new Vector(size);
        Put = new Vector(size);
        
        // check we obtained the new
        if (alpha1 == 0 || alpha1->noRows() != size || alpha1->noCols() != size ||
            alpha3 == 0 || alpha3->noRows() != size || alpha3->noCols() != size ||
            Mhat == 0 || Mhat->noRows() != size || Mhat->noCols() != size ||
            Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Utdothat == 0 || Utdothat->Size() != size ||
            Put == 0 || Put->Size() != size)  {
            
            opserr << "WARNING KRAlphaExplicit_TP::domainChanged() - ";
            opserr << "ran out of memory\n";
            
            // delete the old
            if (alpha1 != 0)
                delete alpha1;
            if (alpha3 != 0)
                delete alpha3;
            if (Mhat != 0)
                delete Mhat;
            if (Ut != 0)
                delete Ut;
            if (Utdot != 0)
                delete Utdot;
            if (Utdotdot != 0)
                delete Utdotdot;
            if (U != 0)
                delete U;
            if (Udot != 0)
                delete Udot;
            if (Udotdot != 0)
                delete Udotdot;
            if (Utdothat != 0)
                delete Utdothat;
            if (Put != 0)
                delete Put;
            
            alpha1 = 0; alpha3 = 0; Mhat = 0;
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Utdothat = 0; Put = 0;
            
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
        for (i=0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*U)(loc) = disp(i);
            }
        }
        
        const Vector &vel = dofPtr->getCommittedVel();
        for (i=0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Udot)(loc) = vel(i);
            }
        }
        
        const Vector &accel = dofPtr->getCommittedAccel();
        for (i=0; i < idSize; i++)  {
            int loc = id(i);
            if (loc >= 0)  {
                (*Udotdot)(loc) = accel(i);
            }
        }
    }
    
    // recalculate integration parameter matrices b/c domain changed
    initAlphaMatrices = 1;
    
    return 0;
}


int KRAlphaExplicit_TP::update(const Vector &aiPlusOne)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING KRAlphaExplicit_TP::update() - called more than once -";
        opserr << " KRAlphaExplicit_TP integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::update() - no AnalysisModel set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check aiPlusOne is of correct size
    if (aiPlusOne.Size() != U->Size())  {
        opserr << "WARNING KRAlphaExplicit_TP::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << aiPlusOne.Size() << endln;
        return -4;
    }
    
    //  determine the response at t+deltaT
    //U->addVector(1.0, aiPlusOne, c1);  // c1 = 0.0
    
    //Udot->addVector(1.0, aiPlusOne, c2);  // c2 = 0.0
    
    Udotdot->addVector(0.0, aiPlusOne, c3);
    
    // update the response at the DOFs
    theModel->setAccel(*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::update() - failed to update the domain\n";
        return -5;
    }
    
    return 0;
}


int KRAlphaExplicit_TP::commit(void)
{
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::commit() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -1;
    }
    
    // set response at t of next step to be that at t+deltaT
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // get unbalance Put and store it for next step
    alphaM = 1.0;
    alphaD = alphaR = alphaP = (1.0 - alphaF);
    Udotdot->addMatrixVector(0.0, *alpha3, *Utdotdot, 1.0);
    theModel->setAccel(*Udotdot);
    this->TransientIntegrator::formUnbalance();
    (*Put) = theLinSOE->getB();
    // reset accelerations at t+deltaT
    (*Udotdot) = *Utdotdot;
    theModel->setAccel(*Udotdot);
    
    return theModel->commitDomain();
}


const Vector &
KRAlphaExplicit_TP::getVel()
{
  return *Udot;
}

int KRAlphaExplicit_TP::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = alphaI;
    data(1) = alphaF;
    data(2) = beta;
    data(3) = gamma;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int KRAlphaExplicit_TP::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING KRAlphaExplicit_TP::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alphaI = data(0);
    alphaF = data(1);
    beta   = data(2);
    gamma  = data(3);
    
    alphaM = 0.0;
    alphaD = alphaF;
    alphaR = alphaF;
    alphaP = alphaF;
    
    return 0;
}


void KRAlphaExplicit_TP::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "KRAlphaExplicit_TP - currentTime: " << currentTime << endln ;
        s << "  alphaI: " << alphaI << "  alphaF: " << alphaF  << "  beta: " << beta  << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else
        s << "KRAlphaExplicit_TP - no associated AnalysisModel\n";
}
