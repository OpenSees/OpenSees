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
// Description: This file contains the implementation of the HHTGeneralizedExplicit_TP class.

#include <HHTGeneralizedExplicit_TP.h>
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
#include <math.h>
#include <elementAPI.h>
#define OPS_Export


void *
OPS_ADD_RUNTIME_VPV(OPS_HHTGeneralizedExplicit_TP)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 2 && argc != 4) {
        opserr << "WARNING - incorrect number of args want HHTGeneralizedExplicit_TP $rhoB $alphaF\n";
        opserr << "          or HHTGeneralizedExplicit_TP $alphaI $alphaF $beta $gamma\n";
        return 0;
    }
    
    double dData[4];
    if (OPS_GetDouble(&argc, dData) != 0) {
        opserr << "WARNING - invalid args want HHTGeneralizedExplicit_TP $rhoB $alphaF\n";
        opserr << "          or HHTGeneralizedExplicit_TP $alphaI $alphaF $beta $gamma\n";
        return 0;
    }
    
    if (argc == 2)
        theIntegrator = new HHTGeneralizedExplicit_TP(dData[0], dData[1]);
    else if (argc == 4)
        theIntegrator = new HHTGeneralizedExplicit_TP(dData[0], dData[1], dData[2], dData[3]);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating HHTGeneralizedExplicit_TP integrator\n";
    
    return theIntegrator;
}


HHTGeneralizedExplicit_TP::HHTGeneralizedExplicit_TP()
    : TransientIntegrator(INTEGRATOR_TAGS_HHTGeneralizedExplicit_TP),
    alphaI(0.5), alphaF(0.5), beta(0.25), gamma(0.5),
    deltaT(0.0), updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    alphaM(0.0), alphaD(0.5), alphaR(0.5), alphaP(0.5),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Put(0)
{
    
}


HHTGeneralizedExplicit_TP::HHTGeneralizedExplicit_TP(double _rhoB, double _alphaF)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTGeneralizedExplicit_TP),
    alphaI((2.0-_rhoB)/(1.0+_rhoB)), alphaF(_alphaF),
    beta((5.0-3*_rhoB+3*_alphaF*(-2.0-_rhoB+pow(_rhoB,2))
    +pow(_alphaF,2)*(2.0+3*_rhoB-pow(_rhoB,3)))
    /((-1.0+_alphaF)*(-2.0+_rhoB)*pow(1.0+_rhoB,2))),
    gamma(0.5+alphaI-_alphaF),
    deltaT(0.0), updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    alphaM(0.0), alphaD(alphaF), alphaR(alphaF), alphaP(alphaF),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Put(0)
{
    
}


HHTGeneralizedExplicit_TP::HHTGeneralizedExplicit_TP(double _alphaI, double _alphaF,
    double _beta, double _gamma)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTGeneralizedExplicit_TP),
    alphaI(_alphaI), alphaF(_alphaF), beta(_beta), gamma(_gamma),
    deltaT(0.0), updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    alphaM(0.0), alphaD(alphaF), alphaR(alphaF), alphaP(alphaF),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Put(0)
{
    
}


HHTGeneralizedExplicit_TP::~HHTGeneralizedExplicit_TP()
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
    if (Put != 0)
        delete Put;
}


int HHTGeneralizedExplicit_TP::newStep(double _deltaT)
{
    updateCount = 0;
    
    if (gamma == 0)  {
        opserr << "HHTExplicit::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << endln;
        return -1;
    }
    
    deltaT = _deltaT;
    if (deltaT <= 0.0)  {
        opserr << "HHTGeneralizedExplicit_TP::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;
    }
    
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::newStep() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -3;
    }
    
    // set the constants
    c1 = beta*deltaT*deltaT;
    c2 = gamma*deltaT;
    c3 = 1.0;
    
    if (U == 0)  {
        opserr << "HHTGeneralizedExplicit_TP::newStep() - domainChange() failed or hasn't been called\n";
        return -4;
    }
    
    // set weighting factors for subsequent iterations
    alphaM = 0.0;
    alphaD = alphaR = alphaP = alphaF;
    
    // determine new displacements and velocities at time t+deltaT
    U->addVector(1.0, *Utdot, deltaT);
    double a1 = (0.5 - beta)*deltaT*deltaT;
    U->addVector(1.0, *Utdotdot, a1);
    
    double a2 = deltaT*(1.0 - gamma);
    Udot->addVector(1.0, *Utdotdot, a2);
    
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
        opserr << "HHTGeneralizedExplicit_TP::newStep() - failed to update the domain\n";
        return -5;
    }
    
    return 0;
}


int HHTGeneralizedExplicit_TP::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int HHTGeneralizedExplicit_TP::formUnbalance()
{
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::formUnbalance() - ";
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
        opserr << "WARNING HHTGeneralizedExplicit_TP::formUnbalance() ";
        opserr << " - this->formElementResidual failed\n";
        return -2;
    }
    
    if (this->formNodalUnbalance() < 0)  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::formUnbalance() ";
        opserr << " - this->formNodalUnbalance failed\n";
        return -3;
    }
    
    return 0;
}


int HHTGeneralizedExplicit_TP::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    theEle->addCtoTang(alphaF*c2);
    theEle->addMtoTang(alphaI*c3);
    
    return 0;
}


int HHTGeneralizedExplicit_TP::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(alphaI*c3);
    
    return 0;
}


int HHTGeneralizedExplicit_TP::formEleResidual(FE_Element *theEle)
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


int HHTGeneralizedExplicit_TP::formNodUnbalance(DOF_Group *theDof)
{
    theDof->zeroUnbalance();
    
    theDof->addPtoUnbalance(alphaP);
    theDof->addD_Force(*Udot, -alphaD);
    theDof->addM_Force(*Udotdot, -alphaM);
    
    return 0;
}


int HHTGeneralizedExplicit_TP::domainChanged()
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
        if (Put != 0)
            delete Put;
        
        // create the new
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        Put = new Vector(size);
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Put == 0 || Put->Size() != size)  {
            
            opserr << "HHTGeneralizedExplicit_TP::domainChanged() - ran out of memory\n";
            
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
            if (Put != 0)
                delete Put;
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Put = 0;
            
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
        for (i=0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                (*U)(loc) = disp(i);
            }
        }
        
        const Vector &vel = dofPtr->getCommittedVel();
        for (i=0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                (*Udot)(loc) = vel(i);
            }
        }
        
        const Vector &accel = dofPtr->getCommittedAccel();
        for (i=0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
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
    
    return 0;
}


int HHTGeneralizedExplicit_TP::update(const Vector &aiPlusOne)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::update() - called more than once -";
        opserr << " HHTGeneralizedExplicit_TP integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::update() - no AnalysisModel set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check aiPlusOne is of correct size
    if (aiPlusOne.Size() != U->Size())  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << aiPlusOne.Size() << endln;
        return -4;
    }
    
    //  determine the response at t+deltaT
    U->addVector(1.0, aiPlusOne, c1);
    
    Udot->addVector(1.0, aiPlusOne, c2);
    
    Udotdot->addVector(0.0, aiPlusOne, c3);
    
    // update the response at the DOFs
    // unfortunately we need to update the displacements in the elements
    // here, so we can not use this method for hybrid simulations
    theModel->setResponse(*U, *Udot, *Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "HHTGeneralizedExplicit_TP::update() - failed to update the domain\n";
        return -5;
    }
    
    return 0;
}


int HHTGeneralizedExplicit_TP::commit(void)
{
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::commit() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -1;
    }
    
    // set response at t of next step to be that at t+deltaT
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
HHTGeneralizedExplicit_TP::getVel()
{
  return *Udot;
}

int HHTGeneralizedExplicit_TP::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = alphaI;
    data(1) = alphaF;
    data(2) = beta;
    data(3) = gamma;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int HHTGeneralizedExplicit_TP::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTGeneralizedExplicit_TP::recvSelf() - could not receive data\n";
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


void HHTGeneralizedExplicit_TP::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "HHTGeneralizedExplicit_TP - currentTime: " << currentTime << endln;
        s << "  alphaI: " << alphaI << "  alphaF: " << alphaF;
        s << "  beta: " << beta  << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else
        s << "HHTGeneralizedExplicit_TP - no associated AnalysisModel\n";
}
