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
// Description: This file contains the implementation of the HHTHSIncrReduct class.

#include <HHTHSIncrReduct.h>
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


void *
OPS_ADD_RUNTIME_VPV(OPS_HHTHSIncrReduct)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 2 && argc != 5) {
        opserr << "WARNING - incorrect number of args want HHTHSIncrReduct $rhoInf $reduct\n";
        opserr << "          or HHTHSIncrReduct $alphaI $alphaF $beta $gamma $reduct\n";
        return 0;
    }
    
    double dData[5];
    if (OPS_GetDouble(&argc, dData) != 0) {
        opserr << "WARNING - invalid args want HHTHSIncrReduct $rhoInf $reduct\n";
        opserr << "          or HHTHSIncrReduct $alphaI $alphaF $beta $gamma $reduct\n";
        return 0;
    }
    
    if (argc == 2)
        theIntegrator = new HHTHSIncrReduct(dData[0], dData[1]);
    else
        theIntegrator = new HHTHSIncrReduct(dData[0], dData[1], dData[2], dData[3], dData[4]);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating HHTHSIncrReduct integrator\n";
    
    return theIntegrator;
}


HHTHSIncrReduct::HHTHSIncrReduct()
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHSIncrReduct),
    alphaI(0.5), alphaF(0.5), beta(0.25), gamma(0.5), reduct(1.0),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0),
    scaledDeltaU(0)
{
    
}


HHTHSIncrReduct::HHTHSIncrReduct(double _rhoInf, double _reduct)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHSIncrReduct),
    alphaI((2.0-_rhoInf)/(1.0+_rhoInf)), alphaF(1.0/(1.0+_rhoInf)),
    beta(1.0/(1.0+_rhoInf)/(1.0+_rhoInf)), gamma(0.5*(3.0-_rhoInf)/(1.0+_rhoInf)),
    reduct(_reduct), deltaT(0.0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0),
    scaledDeltaU(0)
{
    
}


HHTHSIncrReduct::HHTHSIncrReduct(double _alphaI, double _alphaF,
    double _beta, double _gamma, double _reduct)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHSIncrReduct),
    alphaI(_alphaI), alphaF(_alphaF), beta(_beta), gamma(_gamma),
    reduct(_reduct), deltaT(0.0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0),
    scaledDeltaU(0)
{
    
}


HHTHSIncrReduct::~HHTHSIncrReduct()
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
    if (Ualpha != 0)
        delete Ualpha;
    if (Ualphadot != 0)
        delete Ualphadot;
    if (Ualphadotdot != 0)
        delete Ualphadotdot;
    if (scaledDeltaU != 0)
        delete scaledDeltaU;
}


int HHTHSIncrReduct::newStep(double _deltaT)
{
    if (beta == 0 || gamma == 0 )  {
        opserr << "HHTHSIncrReduct::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    deltaT = _deltaT;
    if (deltaT <= 0.0)  {
        opserr << "HHTHSIncrReduct::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;
    }
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    
    // set the constants
    c1 = 1.0;
    c2 = gamma/(beta*deltaT);
    c3 = 1.0/(beta*deltaT*deltaT);
    
    if (U == 0)  {
        opserr << "HHTHSIncrReduct::newStep() - domainChange() failed or hasn't been called\n";
        return -3;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // determine new velocities and accelerations at t+deltaT
    double a1 = (1.0 - gamma/beta);
    double a2 = deltaT*(1.0 - 0.5*gamma/beta);
    Udot->addVector(a1, *Utdotdot, a2);
    
    double a3 = -1.0/(beta*deltaT);
    double a4 = 1.0 - 0.5/beta;
    Udotdot->addVector(a4, *Utdot, a3);
    
    // determine the velocities and accelerations at t+alpha*deltaT
    (*Ualphadot) = *Utdot;
    Ualphadot->addVector((1.0-alphaF), *Udot, alphaF);
    
    (*Ualphadotdot) = *Utdotdot;
    Ualphadotdot->addVector((1.0-alphaI), *Udotdot, alphaI);
    
    // set the trial response quantities
    theModel->setVel(*Ualphadot);
    theModel->setAccel(*Ualphadotdot);
    
    // increment the time to t+alpha*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += alphaF*deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "HHTHSIncrReduct::newStep() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int HHTHSIncrReduct::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int HHTHSIncrReduct::formEleTangent(FE_Element *theEle)
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


int HHTHSIncrReduct::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(alphaI*c3);
    
    return 0;
}


int HHTHSIncrReduct::domainChanged()
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
        if (Ualpha != 0)
            delete Ualpha;
        if (Ualphadot != 0)
            delete Ualphadot;
        if (Ualphadotdot != 0)
            delete Ualphadotdot;
        if (scaledDeltaU != 0)
            delete scaledDeltaU;
        
        // create the new
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        Ualpha = new Vector(size);
        Ualphadot = new Vector(size);
        Ualphadotdot = new Vector(size);
        scaledDeltaU = new Vector(size);
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Ualpha == 0 || Ualpha->Size() != size ||
            Ualphadot == 0 || Ualphadot->Size() != size ||
            Ualphadotdot == 0 || Ualphadotdot->Size() != size ||
            scaledDeltaU == 0 || scaledDeltaU->Size() != size)  {
            
            opserr << "HHTHSIncrReduct::domainChanged() - ran out of memory\n";
            
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
            if (Ualpha != 0)
                delete Ualpha;
            if (Ualphadot != 0)
                delete Ualphadot;
            if (Ualphadotdot != 0)
                delete Ualphadotdot;
            if (scaledDeltaU != 0)
                delete scaledDeltaU;
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Ualpha = 0; Ualphadot = 0; Ualphadotdot = 0;
            scaledDeltaU = 0;
            
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
    
    return 0;
}


int HHTHSIncrReduct::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING HHTHSIncrReduct::update() - no AnalysisModel set\n";
        return -1;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING HHTHSIncrReduct::update() - domainChange() failed or not called\n";
        return -2;
    }
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING HHTHSIncrReduct::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
        return -3;
    }
    
    // get scaled increment
    (*scaledDeltaU) = reduct*deltaU;
    
    // determine the response at t+deltaT
    U->addVector(1.0, *scaledDeltaU, c1);
    
    Udot->addVector(1.0, *scaledDeltaU, c2);
    
    Udotdot->addVector(1.0, *scaledDeltaU, c3);
    
    // determine response at t+alpha*deltaT
    (*Ualpha) = *Ut;
    Ualpha->addVector((1.0-alphaF), *U, alphaF);
    
    (*Ualphadot) = *Utdot;
    Ualphadot->addVector((1.0-alphaF), *Udot, alphaF);
    
    (*Ualphadotdot) = *Utdotdot;
    Ualphadotdot->addVector((1.0-alphaI), *Udotdot, alphaI);
    
    // update the response at the DOFs
    theModel->setResponse(*Ualpha, *Ualphadot, *Ualphadotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "HHTHSIncrReduct::update() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int HHTHSIncrReduct::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING HHTHSIncrReduct::commit() - no AnalysisModel set\n";
        return -1;
    }
    
    // update the response at the DOFs
    theModel->setResponse(*U, *Udot, *Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "HHTHSIncrReduct::commit() - failed to update the domain\n";
        return -2;
    }
    
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-alphaF)*deltaT;
    theModel->setCurrentDomainTime(time);
    
    return theModel->commitDomain();
}


const Vector &
HHTHSIncrReduct::getVel()
{
  return *Udot;
}

int HHTHSIncrReduct::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(5);
    data(0) = alphaI;
    data(1) = alphaF;
    data(2) = beta;
    data(3) = gamma;
    data(4) = reduct;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTHSIncrReduct::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int HHTHSIncrReduct::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(5);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTHSIncrReduct::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alphaI = data(0);
    alphaF = data(1);
    beta   = data(2);
    gamma  = data(3);
    reduct = data(4);
    
    return 0;
}


void HHTHSIncrReduct::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "HHTHSIncrReduct - currentTime: " << currentTime << endln;
        s << "  alphaI: " << alphaI << "  alphaF: " << alphaF;
        s << "  beta: " << beta  << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
        s << "  reductionFactor: " << reduct << endln;
    } else
        s << "HHTHSIncrReduct - no associated AnalysisModel\n";
}
