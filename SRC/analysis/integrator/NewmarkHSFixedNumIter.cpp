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
// Created: 09/05
// Revision: A
//
// Description: This file contains the implementation of the NewmarkHSFixedNumIter class.

#include <NewmarkHSFixedNumIter.h>
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
OPS_ADD_RUNTIME_VPV(OPS_NewmarkHSFixedNumIter)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 2 && argc != 4) {
        opserr << "WARNING - incorrect number of args want NewmarkHSFixedNumIter $gamma $beta <-polyOrder $O>\n";
        return 0;
    }
    
    double dData[2];
    int polyOrder = 2;
    bool updDomFlag = true;
    int numData = 2;
    
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING - invalid args want NewmarkHSFixedNumIter $gamma $beta <-polyOrder $O>\n";
        return 0;
    }
    
    if (argc == 4) {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-polyOrder") == 0) {
            numData = 1;
            if (OPS_GetInt(&numData, &polyOrder) != 0) {
                opserr << "WARNING - invalid polyOrder want NewmarkHSFixedNumIter $gamma $beta <-polyOrder $O>\n";
            }
        }
    }
    
    theIntegrator = new NewmarkHSFixedNumIter(dData[0], dData[1], polyOrder, updDomFlag);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating NewmarkHSFixedNumIter integrator\n";
    
    return theIntegrator;
}


NewmarkHSFixedNumIter::NewmarkHSFixedNumIter()
    : TransientIntegrator(INTEGRATOR_TAGS_NewmarkHSFixedNumIter),
    gamma(0.5), beta(0.25), polyOrder(2), updDomFlag(true),
    c1(0.0), c2(0.0), c3(0.0), x(1.0),
    Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0)
{
    
}


NewmarkHSFixedNumIter::NewmarkHSFixedNumIter(double _gamma,
    double _beta, int polyorder, bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_NewmarkHSFixedNumIter),
    gamma(_gamma), beta(_beta), polyOrder(polyorder), updDomFlag(upddomflag),
    c1(0.0), c2(0.0), c3(0.0), x(1.0),
    Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
    Utm1(0), Utm2(0), scaledDeltaU(0)
{
    
}


NewmarkHSFixedNumIter::~NewmarkHSFixedNumIter()
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


int NewmarkHSFixedNumIter::newStep(double deltaT)
{
    if (beta == 0 || gamma == 0)  {
        opserr << "NewmarkHSFixedNumIter::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "NewmarkHSFixedNumIter::newStep() - error in variable\n";
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
        opserr << "NewmarkHSFixedNumIter::newStep() - domainChange() failed or hasn't been called\n";
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
    
    // set the trial response quantities
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    
    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    theModel->applyLoadDomain(time);
    
    //correctForce = true;
    
    return 0;
}


int NewmarkHSFixedNumIter::revertToLastStep()
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


int NewmarkHSFixedNumIter::formEleTangent(FE_Element *theEle)
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


int NewmarkHSFixedNumIter::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(c2);
    theDof->addMtoTang(c3);
    
    return 0;
}


int NewmarkHSFixedNumIter::domainChanged()
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
            
            opserr << "NewmarkHSFixedNumIter::domainChanged() - ran out of memory\n";
            
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
        opserr << "\nWARNING: NewmarkHSFixedNumIter::domainChanged() - assuming Ut-1 = Ut\n";
    else if (polyOrder == 3)
        opserr << "\nWARNING: NewmarkHSFixedNumIter::domainChanged() - assuming Ut-2 = Ut-1 = Ut\n";
    
    return 0;
}


int NewmarkHSFixedNumIter::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING NewmarkHSFixedNumIter::update() - no AnalysisModel set\n";
        return -1;
    }
    ConvergenceTest *theTest = this->getConvergenceTest();
    if (theTest == 0)  {
        opserr << "WARNING NewmarkHSFixedNumIter::update() - no ConvergenceTest set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING NewmarkHSFixedNumIter::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING NewmarkHSFixedNumIter::update() - Vectors of incompatible size";
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
        opserr << "WARNING NewmarkHSFixedNumIter::update() - polyOrder > 3 not supported\n";
        return -5;
    }
    
    // determine the response at t+deltaT
    U->addVector(1.0, *scaledDeltaU, c1);
    
    Udot->addVector(1.0, *scaledDeltaU, c2);
    
    Udotdot->addVector(1.0, *scaledDeltaU, c3);
    
    // update the response at the DOFs
    theModel->setResponse(*U, *Udot, *Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "NewmarkHSFixedNumIter::update() - failed to update the domain\n";
        return -6;
    }
    
    return 0;
}


int NewmarkHSFixedNumIter::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING NewmarkHSFixedNumIter::commit() - no AnalysisModel set\n";
        return -1;
    }
    
    if (updDomFlag == true)  {
        LinearSOE *theSOE = this->getLinearSOE();
        if (theSOE == 0)  {
            opserr << "WARNING NewmarkHSFixedNumIter::commit() - no LinearSOE set\n";
            return -2;
        }
        
        if (this->formTangent(statusFlag) < 0)  {
            opserr << "WARNING NewmarkHSFixedNumIter::commit() - "
                << "the Integrator failed in formTangent()\n";
            return -3;
        }
        
        if (theSOE->solve() < 0)  {
            opserr << "WARNING NewmarkHSFixedNumIter::commit() - "
                << "the LinearSysOfEqn failed in solve()\n";
            return -4;
        }
        const Vector &deltaU = theSOE->getX();
        
        //  determine the response at t+deltaT
        U->addVector(1.0, deltaU, c1);
        
        Udot->addVector(1.0, deltaU, c2);
        
        Udotdot->addVector(1.0, deltaU, c3);
        
        // update the response at the DOFs
        theModel->setResponse(*U, *Udot, *Udotdot);
    }
    
    return theModel->commitDomain();
}


const Vector &
NewmarkHSFixedNumIter::getVel()
{
  return *Udot;
}

int NewmarkHSFixedNumIter::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = gamma;
    data(1) = beta;
    data(2) = polyOrder;
    if (updDomFlag == true) 
        data(3) = 1.0;
    else
        data(3) = 0.0;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING NewmarkHSFixedNumIter::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int NewmarkHSFixedNumIter::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING NewmarkHSFixedNumIter::recvSelf() - could not receive data\n";
        return -1;
    }
    
    gamma     = data(0);
    beta      = data(1);
    polyOrder = int(data(2));
    if (data(3) == 1.0)
        updDomFlag = true;
    else
        updDomFlag = false;
    
    return 0;
}


void NewmarkHSFixedNumIter::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "NewmarkHSFixedNumIter - currentTime: " << currentTime << endln;
        s << "  gamma: " << gamma << "  beta: " << beta << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
        s << "  polyOrder: " << polyOrder << endln;
        if (updDomFlag)
            s << "  update Domain: yes\n";
        else
            s << "  update Domain: no\n";
    } else 
        s << "NewmarkHSFixedNumIter - no associated AnalysisModel\n";
}


/* force correction does not work if numIter = 1, yields slightly
// improved results if numIter = 2 and has no effect if numIter >= 3
int NewmarkHSFixedNumIter::formElementResidual(void)
{
    // calculate Residual Force
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theSOE = this->getLinearSOE();
    
    // loop through the FE_Elements and add the residual
    FE_Element *elePtr;
    FE_EleIter &theEles = theModel->getFEs();
    while ((elePtr = theEles()) != 0)  {
        if (theSOE->addB(elePtr->getResidual(this),elePtr->getID()) < 0)  {
            opserr << "WARNING NewmarkHSFixedNumIter::formElementResidual() -";
            opserr << " failed in addB for ID " << elePtr->getID();
            return -1;
        }
        if (correctForce)  {
            const Vector &deltaU = theSOE->getX();
            if (statusFlag == CURRENT_TANGENT)  {
                if (theSOE->addB(elePtr->getK_Force(deltaU), elePtr->getID()) < 0)  {
                    opserr << "WARNING NewmarkHSFixedNumIter::formElementResidual() -";
                    opserr << " failed in addB for ID " << elePtr->getID();
                    return -2;
                }
            } else if (statusFlag == INITIAL_TANGENT)  {
                if (theSOE->addB(elePtr->getKi_Force(deltaU), elePtr->getID()) < 0)  {
                    opserr << "WARNING NewmarkHSFixedNumIter::formElementResidual() -";
                    opserr << " failed in addB for ID " << elePtr->getID();
                    return -2;
                }
            }
        }
    }
    
    correctForce = false;
    
    return 0;
}*/
