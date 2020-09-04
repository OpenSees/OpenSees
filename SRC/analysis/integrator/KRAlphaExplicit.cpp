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

#include <KRAlphaExplicit.h>
#include <FE_Element.h>
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


void *    OPS_KRAlphaExplicit(void)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 1 && argc != 2) {
        opserr << "WARNING - incorrect number of args want KRAlphaExplicit $rhoInf <-updateElemDisp>\n";
        return 0;
    }
    
    bool updElemDisp = false;
    double rhoInf;
    int numData = 1;
    if (OPS_GetDouble(&numData, &rhoInf) != 0) {
        opserr << "WARNING - invalid args want KRAlphaExplicit $rhoInf <-updateElemDisp>\n";
        return 0;
    }
    
    if (argc == 2) {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-updateElemDisp") == 0)
            updElemDisp = true;
    }
    
    theIntegrator = new KRAlphaExplicit(rhoInf, updElemDisp);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating KRAlphaExplicit integrator\n";
    
    return theIntegrator;
}


KRAlphaExplicit::KRAlphaExplicit()
    : TransientIntegrator(INTEGRATOR_TAGS_KRAlphaExplicit),
    alphaM(0.5), alphaF(0.5), beta(0.25), gamma(0.5),
    updElemDisp(0), deltaT(0.0),
    alpha1(0), alpha3(0), Mhat(0),
    updateCount(0), initAlphaMatrices(1),
    c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0),
    U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0),
    Utdothat(0)
{
    
}


KRAlphaExplicit::KRAlphaExplicit(double _rhoInf,
    bool updelemdisp)
    : TransientIntegrator(INTEGRATOR_TAGS_KRAlphaExplicit),
    alphaM((2.0-_rhoInf)/(1.0+_rhoInf)), alphaF(1.0/(1.0+_rhoInf)),
    beta(1.0/(1.0+_rhoInf)/(1.0+_rhoInf)), gamma(0.5*(3.0-_rhoInf)/(1.0+_rhoInf)),
    updElemDisp(updelemdisp), deltaT(0.0),
    alpha1(0), alpha3(0), Mhat(0),
    updateCount(0),  initAlphaMatrices(1),
    c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0),
    U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0),
    Utdothat(0)
{
    
}


KRAlphaExplicit::~KRAlphaExplicit()
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
    if (Ualpha != 0)
        delete Ualpha;
    if (Ualphadot != 0)
        delete Ualphadot;
    if (Ualphadotdot != 0)
        delete Ualphadotdot;
    if (Utdothat != 0)
        delete Utdothat;
}


int KRAlphaExplicit::newStep(double _deltaT)
{
    updateCount = 0;
    
    if (beta == 0 || gamma == 0 )  {
        opserr << "WARNING KRAlphaExplicit::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING KRAlphaExplicit::newStep() - no AnalysisModel set\n";
        return -2;
    }
    
    if (initAlphaMatrices || _deltaT != deltaT)  {
        
        // update time step increment
        deltaT = _deltaT;
        if (deltaT <= 0.0)  {
            opserr << "WARNING KRAlphaExplicit::newStep() - error in variable\n";
            opserr << "dT = " << deltaT << endln;
            return -3;
        }
        
        // get the LinearSOE and the ConvergenceTest so we can switch back later
        LinearSOE *theLinSOE = this->getLinearSOE();
        ConvergenceTest *theTest = this->getConvergenceTest();
        
        // set up the FullLinearSOE (needed to compute the alpha matrices)
        int size = theLinSOE->getNumEqn();
        FullGenLinSolver *theFullLinSolver = new FullGenLinLapackSolver();
        LinearSOE *theFullLinSOE = new FullGenLinSOE(size, *theFullLinSolver);
        if (theFullLinSOE == 0)  {
            opserr << "WARNING KRAlphaExplicit::newStep() - failed to create FullLinearSOE\n";
            return -4;
        }
        theFullLinSOE->setLinks(*theModel);
        
        // now switch the SOE to the FullLinearSOE
        this->IncrementalIntegrator::setLinks(*theModel, *theFullLinSOE, theTest);
        
        // get a pointer to the A matrix of the FullLinearSOE
        const Matrix *tmp = theFullLinSOE->getA();
        if (tmp == 0)  {
            opserr << "WARNING KRAlphaExplicit::newStep() - ";
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
        c3 = (1.0 - alphaM);
        this->TransientIntegrator::formTangent(INITIAL_TANGENT);
        Matrix B3(*tmp);
        
        // solve [M + gamma*deltaT*C + beta*deltaT^2*K]*[alpha3] = 
        // [alphaM*M + alphaF*gamma*deltaT*C + alphaF*beta*deltaT^2*K] for alpha3
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
        
        initAlphaMatrices = 0;
    }
    
    if (U == 0)  {
        opserr << "WARNING KRAlphaExplicit::newStep() - domainChange() failed or hasn't been called\n";
        return -6;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // determine new response at time t+deltaT
    Utdothat->addMatrixVector(0.0, *alpha1, *Utdotdot, deltaT);
    
    U->addVector(1.0, *Utdot, deltaT);
    double a1 = (0.5 + gamma)*deltaT;
    U->addVector(1.0, *Utdothat, a1);
    
    Udot->addVector(1.0, *Utdothat, 1.0);
    
    // determine the response at t+alpha*deltaT
    Ualpha->addVector(0.0, *Ut, (1.0-alphaF));
    Ualpha->addVector(1.0, *U, alphaF);
    
    Ualphadot->addVector(0.0, *Utdot, (1.0-alphaF));
    Ualphadot->addVector(1.0, *Udot, alphaF);
    
    Ualphadotdot->addMatrixVector(0.0, *alpha3, *Utdotdot, 1.0);
    
    // set the trial response quantities
    theModel->setResponse(*Ualpha, *Ualphadot, *Ualphadotdot);
    
    // increment the time to t+alpha*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += alphaF*deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "WARNING KRAlphaExplicit::newStep() - failed to update the domain\n";
        return -7;
    }
    
    return 0;
}


int KRAlphaExplicit::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int KRAlphaExplicit::formTangent(int statFlag)
{
    statusFlag = statFlag;
    
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING KRAlphaExplicit::formTangent() - ";
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
        opserr << "WARNING KRAlphaExplicit::formTangent() - ";
        opserr << "failed to add Mhat to A\n";
        return -2;
    }
    
    return 0;
}


int KRAlphaExplicit::formEleTangent(FE_Element *theEle)
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


int KRAlphaExplicit::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(c2);
    theDof->addMtoTang(c3);
    
    return 0;
}


int KRAlphaExplicit::domainChanged()
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
        if (Ualpha != 0)
            delete Ualpha;
        if (Ualphadot != 0)
            delete Ualphadot;
        if (Ualphadotdot != 0)
            delete Ualphadotdot;
        if (Utdothat != 0)
            delete Utdothat;
        
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
        Ualpha = new Vector(size);
        Ualphadot = new Vector(size);
        Ualphadotdot = new Vector(size);
        Utdothat = new Vector(size);
        
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
            Ualpha == 0 || Ualpha->Size() != size ||
            Ualphadot == 0 || Ualphadot->Size() != size ||
            Ualphadotdot == 0 || Ualphadotdot->Size() != size ||
            Utdothat == 0 || Utdothat->Size() != size)  {
            
            opserr << "WARNING KRAlphaExplicit::domainChanged() - ";
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
            if (Ualpha != 0)
                delete Ualpha;
            if (Ualphadot != 0)
                delete Ualphadot;
            if (Ualphadotdot != 0)
                delete Ualphadotdot;
            if (Utdothat != 0)
                delete Utdothat;
            
            alpha1 = 0; alpha3 = 0; Mhat = 0;
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Ualpha = 0; Ualphadot = 0; Ualphadotdot = 0;
            Utdothat = 0;
            
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


int KRAlphaExplicit::update(const Vector &aiPlusOne)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING KRAlphaExplicit::update() - called more than once -";
        opserr << " KRAlphaExplicit integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING KRAlphaExplicit::update() - no AnalysisModel set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING KRAlphaExplicit::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check aiPlusOne is of correct size
    if (aiPlusOne.Size() != U->Size())  {
        opserr << "WARNING KRAlphaExplicit::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << aiPlusOne.Size() << endln;
        return -4;
    }
    
    //  determine the response at t+deltaT
    //U->addVector(1.0, aiPlusOne, c1);  // c1 = 0.0
    
    //Udot->addVector(1.0, aiPlusOne, c2);  // c2 = 0.0
    
    Udotdot->addVector(0.0, aiPlusOne, c3);
    
    // update the response at the DOFs
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "KRAlphaExplicit::update() - failed to update the domain\n";
        return -5;
    }
    // do not update displacements in elements only at nodes
    theModel->setDisp(*U);
    
    return 0;
}


int KRAlphaExplicit::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING KRAlphaExplicit::commit() - no AnalysisModel set\n";
        return -1;
    }
    
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-alphaF)*deltaT;
    theModel->setCurrentDomainTime(time);
    
    // update the displacements in the elements
    if (updElemDisp == true)
        theModel->updateDomain();
    
    return theModel->commitDomain();
}

const Vector &
KRAlphaExplicit::getVel()
{
  return *Udot;
}

int KRAlphaExplicit::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(5);
    data(0) = alphaM;
    data(1) = alphaF;
    data(2) = beta;
    data(3) = gamma;
    if (updElemDisp == false) 
        data(4) = 0.0;
    else
        data(4) = 1.0;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING KRAlphaExplicit::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int KRAlphaExplicit::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(5);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING KRAlphaExplicit::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alphaM = data(0);
    alphaF = data(1);
    beta   = data(2);
    gamma  = data(3);
    if (data(4) == 0.0)
        updElemDisp = false;
    else
        updElemDisp = true;
    
    return 0;
}


void KRAlphaExplicit::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "KRAlphaExplicit - currentTime: " << currentTime << endln ;
        s << "  alphaM: " << alphaM << "  alphaF: " << alphaF  << "  beta: " << beta  << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
        if (updElemDisp)
            s << "  updateElemDisp: yes\n";
        else
            s << "  updateElemDisp: no\n";
    } else
        s << "KRAlphaExplicit - no associated AnalysisModel\n";
}
