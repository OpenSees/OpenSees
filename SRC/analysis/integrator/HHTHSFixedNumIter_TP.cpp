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

// $Revision$
// $Date$
// $URL$

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 10/05
// Revision: A
//
// Description: This file contains the implementation of the HHTHSFixedNumIter_TP class.

#include <HHTHSFixedNumIter_TP.h>
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
#include <ConvergenceTest.h>
#include <elementAPI.h>
#define OPS_Export 


void *
OPS_ADD_RUNTIME_VPV(OPS_HHTHSFixedNumIter_TP)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 1 && argc != 3 && argc != 4 && argc != 6) {
        opserr << "WARNING - incorrect number of args want HHTHSFixedNumIter_TP $rhoInf <-polyOrder $O>\n";
        opserr << "          or HHTHSFixedNumIter_TP $alphaI $alphaF $beta $gamma <-polyOrder $O>\n";
        return 0;
    }
    
    double dData[4];
    int polyOrder = 2;
    bool updDomFlag = true;
    int numData;
    if (argc < 4)
        numData = 1;
    else
        numData = 4;
    
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING - invalid args want HHTHSFixedNumIter_TP $rhoInf <-polyOrder $O>\n";
        opserr << "          or HHTHSFixedNumIter_TP $alphaI $alphaF $beta $gamma <-polyOrder $O>\n";
        return 0;
    }
    
    if (argc == 3 || argc == 6) {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-polyOrder") == 0) {
            numData = 1;
            if (OPS_GetInt(&numData, &polyOrder) != 0) {
                opserr << "WARNING - invalid polyOrder want HHTHSFixedNumIter_TP $rhoInf <-polyOrder $O>\n";
                opserr << "          or HHTHSFixedNumIter_TP $alphaI $alphaF $beta $gamma <-polyOrder $O>\n";
            }
        }
    }
    
    if (argc < 4)
        theIntegrator = new HHTHSFixedNumIter_TP(dData[0], polyOrder, updDomFlag);
    else
        theIntegrator = new HHTHSFixedNumIter_TP(dData[0], dData[1], dData[2], dData[3], polyOrder, updDomFlag);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating HHTHSFixedNumIter_TP integrator\n";
    
    return theIntegrator;
}


HHTHSFixedNumIter_TP::HHTHSFixedNumIter_TP()
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHSFixedNumIter_TP),
    alphaI(0.5), alphaF(0.5), beta(0.25), gamma(0.5),
    polyOrder(2), updDomFlag(true),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0), x(1.0),
    alphaM(0.5), alphaD(0.5), alphaR(0.5), alphaP(0.5),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0), Put(0)
{
    
}


HHTHSFixedNumIter_TP::HHTHSFixedNumIter_TP(double _rhoInf,
    int polyorder, bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHSFixedNumIter_TP),
    alphaI((2.0-_rhoInf)/(1.0+_rhoInf)), alphaF(1.0/(1.0+_rhoInf)),
    beta(1.0/(1.0+_rhoInf)/(1.0+_rhoInf)), gamma(0.5*(3.0-_rhoInf)/(1.0+_rhoInf)),
    polyOrder(polyorder), updDomFlag(upddomflag),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0), x(1.0),
    alphaM(alphaI), alphaD(alphaF), alphaR(alphaF), alphaP(alphaF),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0), Put(0)
{
    
}


HHTHSFixedNumIter_TP::HHTHSFixedNumIter_TP(double _alphaI, double _alphaF,
    double _beta, double _gamma, int polyorder, bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHSFixedNumIter_TP),
    alphaI(_alphaI), alphaF(_alphaF), beta(_beta), gamma(_gamma),
    polyOrder(polyorder), updDomFlag(upddomflag),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0), x(1.0),
    alphaM(alphaI), alphaD(alphaF), alphaR(alphaF), alphaP(alphaF),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0), Put(0)
{
    
}


HHTHSFixedNumIter_TP::~HHTHSFixedNumIter_TP()
{
    // clean up the memory created
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
    if (Utm1 != 0)
        delete Utm1;
    if (Utm2 != 0)
        delete Utm2;
    if (scaledDeltaU != 0)
        delete scaledDeltaU;
    if (Put != 0)
        delete Put;
}


int HHTHSFixedNumIter_TP::newStep(double _deltaT)
{
    if (beta == 0 || gamma == 0 )  {
        opserr << "HHTHSFixedNumIter_TP::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    deltaT = _deltaT;
    if (deltaT <= 0.0)  {
        opserr << "HHTHSFixedNumIter_TP::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;
    }
    
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING HHTHSFixedNumIter_TP::newStep() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -3;
    }
    
    // set the constants
    c1 = 1.0;
    c2 = gamma/(beta*deltaT);
    c3 = 1.0/(beta*deltaT*deltaT);
    
    if (U == 0)  {
        opserr << "HHTHSFixedNumIter_TP::newStep() - domainChange() failed or hasn't been called\n";
        return -4;
    }
    
    // set weighting factors for subsequent iterations
    alphaM = alphaI;
    alphaD = alphaR = alphaP = alphaF;
    
    // determine new velocities and accelerations at t+deltaT
    double a1 = (1.0 - gamma/beta);
    double a2 = deltaT*(1.0 - 0.5*gamma/beta);
    Udot->addVector(a1, *Utdotdot, a2);
    
    double a3 = -1.0/(beta*deltaT);
    double a4 = 1.0 - 0.5/beta;
    Udotdot->addVector(a4, *Utdot, a3);
    
    // set the trial response quantities
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    
    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    theModel->applyLoadDomain(time);
    
    return 0;
}


int HHTHSFixedNumIter_TP::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
        (*Ut) = *Utm1;
        (*Utm1) = *Utm2;
    }
    
    return 0;
}


int HHTHSFixedNumIter_TP::formUnbalance()
{
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING HHTHSFixedNumIter_TP::formUnbalance() - ";
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
        opserr << "WARNING HHTHSFixedNumIter_TP::formUnbalance() ";
        opserr << " - this->formElementResidual failed\n";
        return -2;
    }
    
    if (this->formNodalUnbalance() < 0)  {
        opserr << "WARNING HHTHSFixedNumIter_TP::formUnbalance() ";
        opserr << " - this->formNodalUnbalance failed\n";
        return -3;
    }
    
    return 0;
}


int HHTHSFixedNumIter_TP::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    if (statusFlag == CURRENT_TANGENT)
        theEle->addKtToTang(alphaF*c1);
    else if (statusFlag == INITIAL_TANGENT)
        theEle->addKiToTang(alphaF*c1);
    
    theEle->addCtoTang(alphaF*c2);
    theEle->addMtoTang(alphaI*c3);
    
    return 0;
}


int HHTHSFixedNumIter_TP::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(alphaI*c3);
    
    return 0;
}


int HHTHSFixedNumIter_TP::formEleResidual(FE_Element *theEle)
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


int HHTHSFixedNumIter_TP::formNodUnbalance(DOF_Group *theDof)
{
    theDof->zeroUnbalance();
    
    theDof->addPtoUnbalance(alphaP);
    theDof->addD_Force(*Udot, -alphaD);
    theDof->addM_Force(*Udotdot, -alphaM);
    
    return 0;
}


int HHTHSFixedNumIter_TP::domainChanged()
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    
    // create the new Vector objects
    if (Ut == 0 || Ut->Size() != size)  {
        
        // delete the old
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
        if (Utm1 != 0)
            delete Utm1;
        if (Utm2 != 0)
            delete Utm2;
        if (scaledDeltaU != 0)
            delete scaledDeltaU;
        if (Put != 0)
            delete Put;
        
        // create the new
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        Utm1 = new Vector(size);
        Utm2 = new Vector(size);
        scaledDeltaU = new Vector(size);
        Put = new Vector(size);
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Utm1 == 0 || Utm1->Size() != size ||
            Utm2 == 0 || Utm2->Size() != size ||
            scaledDeltaU == 0 || scaledDeltaU->Size() != size ||
            Put == 0 || Put->Size() != size)  {
            
            opserr << "HHTHSFixedNumIter_TP::domainChanged() - ran out of memory\n";
            
            // delete the old
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
            if (Utm1 != 0)
                delete Utm1;
            if (Utm2 != 0)
                delete Utm2;
            if (scaledDeltaU != 0)
                delete scaledDeltaU;
            if (Put != 0)
                delete Put;
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Utm1 = 0; Utm2 = 0; scaledDeltaU = 0; Put = 0;
            
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
                (*Utm1)(loc) = disp(i);
                (*Ut)(loc) = disp(i);
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
    
    // now get unbalance at last commit and store it
    // warning: this will use committed stiffness prop. damping
    // from current step instead of previous step
    alphaM = (1.0 - alphaI);
    alphaD = alphaR = alphaP = (1.0 - alphaF);
    this->TransientIntegrator::formUnbalance();
    (*Put) = theLinSOE->getB();
    
    if (polyOrder == 2)
        opserr << "\nWARNING: HHTHSFixedNumIter_TP::domainChanged() - assuming Ut-1 = Ut\n";
    else if (polyOrder == 3)
        opserr << "\nWARNING: HHTHSFixedNumIter_TP::domainChanged() - assuming Ut-2 = Ut-1 = Ut\n";
    
    return 0;
}


int HHTHSFixedNumIter_TP::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING HHTHSFixedNumIter_TP::update() - no AnalysisModel set\n";
        return -1;
    }
    ConvergenceTest *theTest = this->getConvergenceTest();
    if (theTest == 0)  {
        opserr << "WARNING HHTHSFixedNumIter_TP::update() - no ConvergenceTest set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING HHTHSFixedNumIter_TP::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING HHTHSFixedNumIter_TP::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
        return -4;
    }
    
    // get interpolation location and scale displacement increment 
    x = (double) theTest->getNumTests()/theTest->getMaxNumTests();
    if (polyOrder == 1)  {
        (*scaledDeltaU) = x*((*U)+deltaU) - (x-1.0)*(*Ut)  - (*U);
    }
    else if (polyOrder == 2)  {
        (*scaledDeltaU) = x*(x+1.0)/2.0*((*U)+deltaU) - (x-1.0)*(x+1.0)*(*Ut) 
                        + (x-1.0)*x/2.0*(*Utm1) - (*U);
    }
    else if (polyOrder == 3)  {
        (*scaledDeltaU) = x*(x+1.0)*(x+2.0)/6.0*((*U)+deltaU) - (x-1.0)*(x+1.0)*(x+2.0)/2.0*(*Ut)
                        + (x-1.0)*x*(x+2.0)/2.0*(*Utm1) - (x-1.0)*x*(x+1.0)/6.0*(*Utm2) - (*U);
    }
    else  {
        opserr << "WARNING HHTHSFixedNumIter_TP::update() - polyOrder > 3 not supported\n";
        return -5;
    }
    
    // determine the response at t+deltaT
    U->addVector(1.0, *scaledDeltaU, c1);
    
    Udot->addVector(1.0, *scaledDeltaU, c2);
    
    Udotdot->addVector(1.0, *scaledDeltaU, c3);
    
    // update the response at the DOFs
    theModel->setResponse(*U, *Udot, *Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "HHTHSFixedNumIter_TP::update() - failed to update the domain\n";
        return -6;
    }
    
    return 0;
}


int HHTHSFixedNumIter_TP::commit(void)
{
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING HHTHSFixedNumIter_TP::commit() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -1;
    }
    
    if (updDomFlag == true)  {
        if (this->formTangent(statusFlag) < 0)  {
            opserr << "WARNING HHTHSFixedNumIter_TP::commit() - "
                << "the Integrator failed in formTangent()\n";
            return -2;
        }
        
        if (theLinSOE->solve() < 0)  {
            opserr << "WARNING HHTHSFixedNumIter_TP::commit() - "
                << "the LinearSysOfEqn failed in solve()\n";
            return -3;
        }
        const Vector &deltaU = theLinSOE->getX();
        
        //  determine the response at t+deltaT
        U->addVector(1.0, deltaU, c1);
        
        Udot->addVector(1.0, deltaU, c2);
        
        Udotdot->addVector(1.0, deltaU, c3);
        
        // update the response at the DOFs
        theModel->setResponse(*U, *Udot, *Udotdot);
    }
    
    // set response at t of next step to be that at t+deltaT
    (*Utm2) = *Utm1;
    (*Utm1) = *Ut;
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // get unbalance Put and store it for next step
    alphaM = (1.0 - alphaI);
    alphaD = alphaR = alphaP = (1.0 - alphaF);
    this->TransientIntegrator::formUnbalance();
    (*Put) = theLinSOE->getB();
    
    return theModel->commitDomain();
}

const Vector &
HHTHSFixedNumIter_TP::getVel()
{
  return *Udot;
}

int HHTHSFixedNumIter_TP::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(6);
    data(0) = alphaI;
    data(1) = alphaF;
    data(2) = beta;
    data(3) = gamma;
    data(4) = polyOrder;
    if (updDomFlag == true) 
        data(5) = 1.0;
    else
        data(5) = 0.0;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTHSFixedNumIter_TP::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int HHTHSFixedNumIter_TP::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(6);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTHSFixedNumIter_TP::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alphaI    = data(0);
    alphaF    = data(1);
    beta      = data(2);
    gamma     = data(3);
    polyOrder = int(data(4));
    if (data(5) == 1.0)
        updDomFlag = true;
    else
        updDomFlag = false;
    
    alphaM = alphaI;
    alphaD = alphaF;
    alphaR = alphaF;
    alphaP = alphaF;
    
    return 0;
}


void HHTHSFixedNumIter_TP::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "HHTHSFixedNumIter_TP - currentTime: " << currentTime << endln;
        s << "  alphaI: " << alphaI << "  alphaF: " << alphaF;
        s << "  beta: " << beta  << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
        s << "  polyOrder: " << polyOrder << endln;
        if (updDomFlag)
            s << "  update Domain: yes\n";
        else
            s << "  update Domain: no\n";
    } else
        s << "HHTHSFixedNumIter_TP - no associated AnalysisModel\n";
}
