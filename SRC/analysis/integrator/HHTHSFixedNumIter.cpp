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
// Description: This file contains the implementation of the HHTHSFixedNumIter class.

#include <HHTHSFixedNumIter.h>
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
OPS_ADD_RUNTIME_VPV(OPS_HHTHSFixedNumIter)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 1 && argc != 3 && argc != 4 && argc != 6) {
        opserr << "WARNING - incorrect number of args want HHTHSFixedNumIter $rhoInf <-polyOrder $O>\n";
        opserr << "          or HHTHSFixedNumIter $alphaI $alphaF $beta $gamma <-polyOrder $O>\n";
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
        opserr << "WARNING - invalid args want HHTHSFixedNumIter $rhoInf <-polyOrder $O>\n";
        opserr << "          or HHTHSFixedNumIter $alphaI $alphaF $beta $gamma <-polyOrder $O>\n";
        return 0;
    }
    
    if (argc == 3 || argc == 6) {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-polyOrder") == 0) {
            numData = 1;
            if (OPS_GetInt(&numData, &polyOrder) != 0) {
                opserr << "WARNING - invalid polyOrder want HHTHSFixedNumIter $rhoInf <-polyOrder $O>\n";
                opserr << "          or HHTHSFixedNumIter $alphaI $alphaF $beta $gamma <-polyOrder $O>\n";
            }
        }
    }
    
    if (argc < 4)
        theIntegrator = new HHTHSFixedNumIter(dData[0], polyOrder, updDomFlag);
    else
        theIntegrator = new HHTHSFixedNumIter(dData[0], dData[1], dData[2], dData[3], polyOrder, updDomFlag);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating HHTHSFixedNumIter integrator\n";
    
    return theIntegrator;
}


HHTHSFixedNumIter::HHTHSFixedNumIter()
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHSFixedNumIter),
    alphaI(0.5), alphaF(0.5), beta(0.25), gamma(0.5),
    polyOrder(2), updDomFlag(true),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0), x(1.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0)
{
    
}


HHTHSFixedNumIter::HHTHSFixedNumIter(double _rhoInf,
    int polyorder, bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHSFixedNumIter),
    alphaI((2.0-_rhoInf)/(1.0+_rhoInf)), alphaF(1.0/(1.0+_rhoInf)),
    beta(1.0/(1.0+_rhoInf)/(1.0+_rhoInf)), gamma(0.5*(3.0-_rhoInf)/(1.0+_rhoInf)),
    polyOrder(polyorder), updDomFlag(upddomflag),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0), x(1.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0)
{
    
}


HHTHSFixedNumIter::HHTHSFixedNumIter(double _alphaI, double _alphaF,
    double _beta, double _gamma, int polyorder, bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHSFixedNumIter),
    alphaI(_alphaI), alphaF(_alphaF), beta(_beta), gamma(_gamma),
    polyOrder(polyorder), updDomFlag(upddomflag),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0), x(1.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0)
{
    
}


HHTHSFixedNumIter::~HHTHSFixedNumIter()
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
    if (Utm1 != 0)
        delete Utm1;
    if (Utm2 != 0)
        delete Utm2;
    if (scaledDeltaU != 0)
        delete scaledDeltaU;
}


int HHTHSFixedNumIter::newStep(double _deltaT)
{
    deltaT = _deltaT;
    if (beta == 0 || gamma == 0 )  {
        opserr << "HHTHSFixedNumIter::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "HHTHSFixedNumIter::newStep() - error in variable\n";
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
        opserr << "HHTHSFixedNumIter::newStep() - domainChange() failed or hasn't been called\n";
        return -3;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Utm2) = *Utm1;
    (*Utm1) = *Ut;
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
    theModel->applyLoadDomain(time);
    
    return 0;
}


int HHTHSFixedNumIter::revertToLastStep()
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


int HHTHSFixedNumIter::formEleTangent(FE_Element *theEle)
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


int HHTHSFixedNumIter::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(alphaI*c3);
    
    return 0;
}


int HHTHSFixedNumIter::domainChanged()
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
        if (Utm1 != 0)
            delete Utm1;
        if (Utm2 != 0)
            delete Utm2;
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
        Utm1 = new Vector(size);
        Utm2 = new Vector(size);
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
            Utm1 == 0 || Utm1->Size() != size ||
            Utm2 == 0 || Utm2->Size() != size ||
            scaledDeltaU == 0 || scaledDeltaU->Size() != size)  {
            
            opserr << "HHTHSFixedNumIter::domainChanged() - ran out of memory\n";
            
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
            if (Utm1 != 0)
                delete Utm1;
            if (Utm2 != 0)
                delete Utm2;
            if (scaledDeltaU != 0)
                delete scaledDeltaU;
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Ualpha = 0; Ualphadot = 0; Ualphadotdot = 0;
            Utm1 = 0; Utm2 = 0; scaledDeltaU = 0;
            
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
    
    if (polyOrder == 2)
        opserr << "\nWARNING: HHTHSFixedNumIter::domainChanged() - assuming Ut-1 = Ut\n";
    else if (polyOrder == 3)
        opserr << "\nWARNING: HHTHSFixedNumIter::domainChanged() - assuming Ut-2 = Ut-1 = Ut\n";
    
    return 0;
}


int HHTHSFixedNumIter::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING HHTHSFixedNumIter::update() - no AnalysisModel set\n";
        return -1;
    }
    ConvergenceTest *theTest = this->getConvergenceTest();
    if (theTest == 0)  {
        opserr << "WARNING HHTHSFixedNumIter::update() - no ConvergenceTest set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING HHTHSFixedNumIter::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING HHTHSFixedNumIter::update() - Vectors of incompatible size ";
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
        opserr << "WARNING HHTHSFixedNumIter::update() - polyOrder > 3 not supported\n";
        return -5;
    }
    
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
        opserr << "HHTHSFixedNumIter::update() - failed to update the domain\n";
        return -6;
    }
    
    return 0;
}


int HHTHSFixedNumIter::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING HHTHSFixedNumIter::commit() - no AnalysisModel set\n";
        return -1;
    }
    
    if (updDomFlag == true)  {
        LinearSOE *theSOE = this->getLinearSOE();
        if (theSOE == 0)  {
            opserr << "WARNING HHTHSFixedNumIter::commit() - no LinearSOE set\n";
            return -2;
        }
        
        if (this->formTangent(statusFlag) < 0)  {
            opserr << "WARNING HHTHSFixedNumIter::commit() - "
                << "the Integrator failed in formTangent()\n";
            return -3;
        }
        
        if (theSOE->solve() < 0)  {
            opserr << "WARNING HHTHSFixedNumIter::commit() - "
                << "the LinearSysOfEqn failed in solve()\n";
            return -4;
        }
        const Vector &deltaU = theSOE->getX();
        
        //  determine the response at t+deltaT
        U->addVector(1.0, deltaU, c1);
        
        Udot->addVector(1.0, deltaU, c2);
        
        Udotdot->addVector(1.0, deltaU, c3);
    }
    
    // update the response at the DOFs
    theModel->setResponse(*U, *Udot, *Udotdot);
    
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-alphaF)*deltaT;
    theModel->setCurrentDomainTime(time);
    
    return theModel->commitDomain();
}

const Vector &
HHTHSFixedNumIter::getVel()
{
  return *Udot;
}

int HHTHSFixedNumIter::sendSelf(int cTag, Channel &theChannel)
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
        opserr << "WARNING HHTHSFixedNumIter::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int HHTHSFixedNumIter::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(6);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTHSFixedNumIter::recvSelf() - could not receive data\n";
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
    
    return 0;
}


void HHTHSFixedNumIter::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "HHTHSFixedNumIter - currentTime: " << currentTime << endln;
        s << "  alphaI: " << alphaI << "  alphaF: " << alphaF;
        s << "  beta: " << beta  << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
        s << "  polyOrder: " << polyOrder << endln;
        if (updDomFlag)
            s << "  update Domain: yes\n";
        else
            s << "  update Domain: no\n";
    } else
        s << "HHTHSFixedNumIter - no associated AnalysisModel\n";
}
