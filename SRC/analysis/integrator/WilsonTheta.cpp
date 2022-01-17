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

// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the implementation of WilsonTheta.

#include <WilsonTheta.h>
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
OPS_ADD_RUNTIME_VPV(OPS_WilsonTheta)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc != 1) {
        opserr << "WARNING - incorrect number of args want WilsonTheta $theta\n";
        return 0;
    }
    
    double theta;    
    if (OPS_GetDouble(&argc, &theta) != 0) {
        opserr << "WARNING - invalid args want WilsonTheta $theta\n";
        return 0;
    }
    
    theIntegrator = new WilsonTheta(theta);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating WilsonTheta integrator\n";
    
    return theIntegrator;
}


WilsonTheta::WilsonTheta()
    : TransientIntegrator(INTEGRATOR_TAGS_WilsonTheta),
    theta(0), deltaT(0),
    c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0)
{
    
}


WilsonTheta::WilsonTheta(double _theta)
    : TransientIntegrator(INTEGRATOR_TAGS_WilsonTheta),
    theta(_theta), deltaT(0.0),
    c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0)
{
    
}


WilsonTheta::~WilsonTheta()
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
}


int WilsonTheta::newStep(double _deltaT)
{
    deltaT = _deltaT;
    if (theta <= 0.0 )  {
        opserr << "WilsonTheta::newStep() - error in variable\n";
        opserr << "theta: " << theta << " <= 0.0\n";
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "WilsonTheta::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;
    }
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    
    // set the constants
    c1 = 1.0;
    c2 = 3.0/(theta*deltaT);
    c3 = 2.0*c2/(theta*deltaT);
    
    if (U == 0)  {
        opserr << "WilsonTheta::newStep() - domainChange() failed or hasn't been called\n"; 
        return -3;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // determine new velocities and accelerations at t+theta*deltaT
    double a1 = -0.5*theta*deltaT;
    Udot->addVector(-2.0, *Utdotdot, a1);
    
    double a2 = -6.0/theta/deltaT; 
    Udotdot->addVector(-2.0, *Utdot, a2);
    
    // set the trial response quantities
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    
    // increment the time to t+theta*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += theta*deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "WilsonTheta::newStep() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int WilsonTheta::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int WilsonTheta::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    if (statusFlag == CURRENT_TANGENT)
        theEle->addKtToTang(c1);
    else if (statusFlag == INITIAL_TANGENT)
        theEle->addKiToTang(c1);
    else if (statusFlag == HALL_TANGENT)  {
        theEle->addKtToTang(c1*cFactor);
        theEle->addKiToTang(c1*iFactor);   
    }
 
    theEle->addCtoTang(c2);
    theEle->addMtoTang(c3);
    
    return 0;
}


int WilsonTheta::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(c2);
    theDof->addMtoTang(c3);
    
    return 0;
}


int WilsonTheta::domainChanged()
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
        
        // create the new
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size) {
            
            opserr << "WilsonTheta::domainChanged() - ran out of memory\n";
            
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
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            
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


int WilsonTheta::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING WilsonTheta::update() - no AnalysisModel set\n";
        return -1;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING WilsonTheta::update() - domainChange() failed or not called\n";
        return -2;
    }
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING WilsonTheta::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
        return -3;
    }
    
    // determine the response at t+theta*deltaT
    U->addVector(1.0, deltaU, c1);
    
    Udot->addVector(1.0, deltaU, c2);
    
    Udotdot->addVector(1.0, deltaU, c3);
    
    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "WilsonTheta::update() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int WilsonTheta::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING WilsonTheta::commit() - no AnalysisModel set\n";
        return -1;
    }
        
    // determine response quantities at t+deltaT
    Udotdot->addVector(1.0/theta, *Utdotdot, (theta-1.0)/theta);
    
    (*Udot) = *Utdot;
    double a1 = 0.5*deltaT;
    Udot->addVector(1.0, *Udotdot, a1);
    Udot->addVector(1.0, *Utdotdot, a1);
    
    (*U) = *Ut;
    U->addVector(1.0, *Utdot, deltaT);
    double a2 = deltaT*deltaT/6.0;
    U->addVector(1.0, *Udotdot, a2);
    U->addVector(1.0, *Utdotdot, 2.0*a2);
    
    // update the response at the DOFs 
    theModel->setResponse(*U,*Udot,*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "WilsonTheta::commit() - failed to update the domain\n";
        return -2;
    }
    
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-theta)*deltaT;
    theModel->setCurrentDomainTime(time);
    
    return theModel->commitDomain();
}


const Vector &
WilsonTheta::getVel()
{
  return *Udot;
}

int WilsonTheta::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(1);
    data(0) = theta;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WilsonTheta::sendSelf() - failed to send the data\n";
        return -1;
    }
    
    return 0;
}


int WilsonTheta::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(1);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING WilsonTheta::recvSelf() - could not receive data\n";
        return -1;
    }
    
    theta  = data(0);
    
    return 0;
}


void WilsonTheta::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t WilsonTheta - currentTime: " << currentTime << endln;
        s << "  theta: " << theta << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else
        s << "\t WilsonTheta - no associated AnalysisModel\n";
}
