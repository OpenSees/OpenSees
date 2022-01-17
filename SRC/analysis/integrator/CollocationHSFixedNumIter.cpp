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
// Description: This file contains the implementation of CollocationHSFixedNumIter.

#include <CollocationHSFixedNumIter.h>
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
#include <math.h>
#include <elementAPI.h>
#define OPS_Export 


void *
OPS_ADD_RUNTIME_VPV(OPS_CollocationHSFixedNumIter)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 1 && argc != 3 && argc != 5) {
        opserr << "WARNING - incorrect number of args want CollocationHSFixedNumIter $theta <-polyOrder $O>\n";
        opserr << "          or CollocationHSFixedNumIter $theta $beta $gamma <-polyOrder $O>\n";
        return 0;
    }
    
    double dData[3];
    int polyOrder = 2;
    int numData = 0;
    
    // count number of numeric parameters
    while (OPS_GetNumRemainingInputArgs() > 0)  {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-polyOrder") == 0) {
            break;
        }
        numData++;
    }
    // reset to read from beginning
    OPS_ResetCurrentInputArg(2);
    
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING - invalid args want CollocationHSFixedNumIter $theta <-polyOrder $O>\n";
        opserr << "          or CollocationHSFixedNumIter $theta $beta $gamma <-polyOrder $O>\n";
        return 0;
    }
    
    if (numData+2 == argc) {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-polyOrder") == 0) {
            int numData2 = 1;
            if (OPS_GetInt(&numData2, &polyOrder) != 0) {
                opserr << "WARNING - invalid polyOrder want CollocationHSFixedNumIter $rhoInf <-polyOrder $O>\n";
                opserr << "          or CollocationHSFixedNumIter $alphaI $alphaF $beta $gamma <-polyOrder $O>\n";
            }
        }
    }
    
    if (numData == 1)
        theIntegrator = new CollocationHSFixedNumIter(dData[0], polyOrder);
    else if (numData == 3)
        theIntegrator = new CollocationHSFixedNumIter(dData[0], dData[1], dData[2], polyOrder);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating CollocationHSFixedNumIter integrator\n";
    
    return theIntegrator;
}


CollocationHSFixedNumIter::CollocationHSFixedNumIter()
    : TransientIntegrator(INTEGRATOR_TAGS_CollocationHSFixedNumIter),
    theta(1.0), beta(0.25), gamma(0.5), polyOrder(2),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0), x(1.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0)
{
    
}


CollocationHSFixedNumIter::CollocationHSFixedNumIter(
    double _theta, int polyorder)
    : TransientIntegrator(INTEGRATOR_TAGS_CollocationHSFixedNumIter),
    theta(_theta), beta(0.0), gamma(0.5), polyOrder(polyorder),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0), x(1.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0)
{
    beta = -6.018722044382699e+002 * pow(theta,9) +
            6.618777151634235e+003 * pow(theta,8) +
           -3.231561059595987e+004 * pow(theta,7) +
            9.195359004558867e+004 * pow(theta,6) +
           -1.680788908312227e+005 * pow(theta,5) +
            2.047005794710718e+005 * pow(theta,4) +
           -1.661421563528177e+005 * pow(theta,3) +
            8.667950092619179e+004 * pow(theta,2) +
           -2.638652989051994e+004 * theta +
            3.572862280471971e+003;
}


CollocationHSFixedNumIter::CollocationHSFixedNumIter(
    double _theta, double _beta, double _gamma, int polyorder)
    : TransientIntegrator(INTEGRATOR_TAGS_CollocationHSFixedNumIter),
    theta(_theta), beta(_beta), gamma(_gamma), polyOrder(polyorder),
    deltaT(0.0), c1(0.0), c2(0.0), c3(0.0), x(1.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0)
{
    
}


CollocationHSFixedNumIter::~CollocationHSFixedNumIter()
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
}


int CollocationHSFixedNumIter::newStep(double _deltaT)
{
    if (theta <= 0.0 )  {
        opserr << "CollocationHSFixedNumIter::newStep() - error in variable\n";
        opserr << "theta: " << theta << " <= 0.0\n";
        return -1;
    }
    
    deltaT = _deltaT;
    if (deltaT <= 0.0)  {
        opserr << "CollocationHSFixedNumIter::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;
    }
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    
    // set the constants
    c1 = 1.0;
    c2 = gamma/(beta*theta*deltaT);
    c3 = 1.0/(beta*theta*theta*deltaT*deltaT);
    
    if (U == 0)  {
        opserr << "CollocationHSFixedNumIter::newStep() - domainChange() failed or hasn't been called\n"; 
        return -3;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Utm2) = *Utm1;
    (*Utm1) = *Ut;
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // determine new velocities and accelerations at t+theta*deltaT
    double a1 = (1.0 - gamma/beta);
    double a2 = theta*deltaT*(1.0 - 0.5*gamma/beta);
    Udot->addVector(a1, *Utdotdot, a2);
    
    double a3 = -1.0/(beta*theta*deltaT);
    double a4 = 1.0 - 0.5/beta;
    Udotdot->addVector(a4, *Utdot, a3);
    
    // set the trial response quantities
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    
    // increment the time to t+theta*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += theta*deltaT;
    theModel->applyLoadDomain(time);
    
    return 0;
}


int CollocationHSFixedNumIter::revertToLastStep()
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


int CollocationHSFixedNumIter::formEleTangent(FE_Element *theEle)
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


int CollocationHSFixedNumIter::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(c2);
    theDof->addMtoTang(c3);
    
    return 0;
}


int CollocationHSFixedNumIter::domainChanged()
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
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Utm1 == 0 || Utm1->Size() != size ||
            Utm2 == 0 || Utm2->Size() != size ||
            scaledDeltaU == 0 || scaledDeltaU->Size() != size)  {
            
            opserr << "CollocationHSFixedNumIter::domainChanged() - ran out of memory\n";
            
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
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
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
        opserr << "\nWARNING: CollocationHSFixedNumIter::domainChanged() - assuming Ut-1 = Ut\n";
    else if (polyOrder == 3)
        opserr << "\nWARNING: CollocationHSFixedNumIter::domainChanged() - assuming Ut-2 = Ut-1 = Ut\n";
    
    return 0;
}


int CollocationHSFixedNumIter::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING CollocationHSFixedNumIter::update() - no AnalysisModel set\n";
        return -1;
    }
    ConvergenceTest *theTest = this->getConvergenceTest();
    if (theTest == 0)  {
        opserr << "WARNING CollocationHSFixedNumIter::update() - no ConvergenceTest set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING CollocationHSFixedNumIter::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING CollocationHSFixedNumIter::update() - Vectors of incompatible size ";
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
        opserr << "WARNING CollocationHSFixedNumIter::update() - polyOrder > 3 not supported\n";
        return -5;
    }
    
    // determine the response at t+theta*deltaT
    U->addVector(1.0, *scaledDeltaU, c1);
    
    Udot->addVector(1.0, *scaledDeltaU, c2);
    
    Udotdot->addVector(1.0, *scaledDeltaU, c3);
    
    // update the response at the DOFs
    theModel->setResponse(*U, *Udot, *Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "CollocationHSFixedNumIter::update() - failed to update the domain\n";
        return -5;
    }
    
    return 0;
}


int CollocationHSFixedNumIter::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING CollocationHSFixedNumIter::commit() - no AnalysisModel set\n";
        return -1;
    }
    
    LinearSOE *theSOE = this->getLinearSOE();
    if (theSOE == 0)  {
        opserr << "WARNING CollocationHSFixedNumIter::commit() - no LinearSOE set\n";
        return -2;
    }
    
    if (theSOE->solve() < 0)  {
        opserr << "WARNING CollocationHSFixedNumIter::commit() - "
            << "the LinearSysOfEqn failed in solve()\n";
        return -3;
    }
    const Vector &deltaU = theSOE->getX();
    
    // determine the response at t+theta*deltaT
    U->addVector(1.0, deltaU, c1);
    
    Udot->addVector(1.0, deltaU, c2);
    
    Udotdot->addVector(1.0, deltaU, c3);
    
    // determine response quantities at t+deltaT
    Udotdot->addVector(1.0/theta, *Utdotdot, (theta-1.0)/theta);
    
    (*Udot) = *Utdot;
    double a1 = deltaT*(1.0 - gamma);
    double a2 = deltaT*gamma;
    Udot->addVector(1.0, *Utdotdot, a1);
    Udot->addVector(1.0, *Udotdot, a2);
    
    (*U) = *Ut;
    U->addVector(1.0, *Utdot, deltaT);
    double a3 = deltaT*deltaT*(0.5 - beta);
    double a4 = deltaT*deltaT*beta;
    U->addVector(1.0, *Utdotdot, a3);
    U->addVector(1.0, *Udotdot, a4);
    
    // update the response at the DOFs
    theModel->setResponse(*U, *Udot, *Udotdot);
    
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-theta)*deltaT;
    theModel->setCurrentDomainTime(time);
    
    return theModel->commitDomain();
}


const Vector &
CollocationHSFixedNumIter::getVel()
{
  return *Udot;
}

int CollocationHSFixedNumIter::sendSelf(int cTag, Channel &theChannel)
{
    static Vector data(4);
    data(0) = theta;
    data(1) = beta;
    data(2) = gamma;
    data(3) = polyOrder;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING CollocationHSFixedNumIter::sendSelf() - failed to send the data\n";
        return -1;
    }
    
    return 0;
}


int CollocationHSFixedNumIter::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING CollocationHSFixedNumIter::recvSelf() - could not receive data\n";
        return -1;
    }
    
    theta     = data(0);
    beta      = data(1);
    gamma     = data(2);
    polyOrder = int(data(3));
    
    return 0;
}


void CollocationHSFixedNumIter::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "CollocationHSFixedNumIter - currentTime: " << currentTime << endln;
        s << "  theta: " << theta << endln;
        s << "  polyOrder: " << polyOrder << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else 
        s << "CollocationHSFixedNumIter - no associated AnalysisModel\n";
}
