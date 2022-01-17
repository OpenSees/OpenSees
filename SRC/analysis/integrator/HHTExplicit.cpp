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
// Description: This file contains the implementation of the HHTExplicit class.

#include <HHTExplicit.h>
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
OPS_ADD_RUNTIME_VPV(OPS_HHTExplicit)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc < 1 || argc > 3) {
        opserr << "WARNING - incorrect number of args want HHTExplicit $alpha <-updateElemDisp>\n";
        opserr << "          or HHTExplicit $alpha $gamma <-updateElemDisp>\n";
        return 0;
    }
    
    bool updElemDisp = false;
    double dData[2];
    int numData = 0;
    
    // count number of numeric parameters
    while (OPS_GetNumRemainingInputArgs() > 0)  {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-updateElemDisp") == 0) {
            break;
        }
        numData++;
    }
    // reset to read from beginning
    OPS_ResetCurrentInputArg(2);
    
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING - invalid args want HHTExplicit $alpha <-updateElemDisp>\n";
        opserr << "          or HHTExplicit $alpha $gamma <-updateElemDisp>\n";
        return 0;
    }
    
    if (numData+1 == argc) {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-updateElemDisp") == 0)
            updElemDisp = true;
    }
    
    if (numData == 1)
        theIntegrator = new HHTExplicit(dData[0], updElemDisp);
    else if (numData == 2)
        theIntegrator = new HHTExplicit(dData[0], dData[1], updElemDisp);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating HHTExplicit integrator\n";
    
    return theIntegrator;
}


HHTExplicit::HHTExplicit()
    : TransientIntegrator(INTEGRATOR_TAGS_HHTExplicit),
    alpha(1.0), gamma(0.5), updElemDisp(false),
    deltaT(0.0), updateCount(0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0)
{
    
}


HHTExplicit::HHTExplicit(double _alpha,
    bool updelemdisp)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTExplicit),
    alpha(_alpha), gamma(0.5), updElemDisp(updelemdisp),
    deltaT(0.0), updateCount(0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0)
{
    
}


HHTExplicit::HHTExplicit(double _alpha, double _gamma,
    bool updelemdisp)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTExplicit),
    alpha(_alpha), gamma(_gamma), updElemDisp(updelemdisp),
    deltaT(0.0), updateCount(0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0)
{
    
}


HHTExplicit::~HHTExplicit()
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
}


int HHTExplicit::newStep(double _deltaT)
{
    updateCount = 0;
    
    if (gamma == 0)  {
        opserr << "HHTExplicit::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << endln;
        return -1;
    }
    
    deltaT = _deltaT;
    if (deltaT <= 0.0)  {
        opserr << "HHTExplicit::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;
    }
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    
    // set the constants
    c2 = gamma*deltaT;
    c3 = 1.0;
    
    if (U == 0)  {
        opserr << "HHTExplicit::newStep() - domainChange() failed or hasn't been called\n";
        return -3;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // determine new displacements and velocities at time t+deltaT
    U->addVector(1.0, *Utdot, deltaT);
    double a1 = 0.5*deltaT*deltaT;
    U->addVector(1.0, *Utdotdot, a1);
    
    double a2 = deltaT*(1.0 - gamma);
    Udot->addVector(1.0, *Utdotdot, a2);
    
    // determine the response at t+alpha*deltaT
    (*Ualpha) = *Ut;
    Ualpha->addVector((1.0-alpha), *U, alpha);
    
    (*Ualphadot) = *Utdot;
    Ualphadot->addVector((1.0-alpha), *Udot, alpha);
    
    Udotdot->Zero();
    
    // set the trial response quantities for the elements
    theModel->setResponse(*Ualpha, *Ualphadot, *Udotdot);
    
    // increment the time to t+alpha*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += alpha*deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "HHTExplicit::newStep() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int HHTExplicit::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int HHTExplicit::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    theEle->addCtoTang(alpha*c2);
    theEle->addMtoTang(c3);
    
    return 0;
}


int HHTExplicit::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(alpha*c2);
    theDof->addMtoTang(c3);
    
    return 0;
}


int HHTExplicit::domainChanged()
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
        
        // create the new
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        Ualpha = new Vector(size);
        Ualphadot = new Vector(size);
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Ualpha == 0 || Ualpha->Size() != size ||
            Ualphadot == 0 || Ualphadot->Size() != size)  {
            
            opserr << "HHTExplicit::domainChanged() - ran out of memory\n";
            
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
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Ualpha = 0; Ualphadot = 0;
            
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
    
    return 0;
}


int HHTExplicit::update(const Vector &aiPlusOne)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING HHTExplicit::update() - called more than once -";
        opserr << " HHTExplicit integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING HHTExplicit::update() - no AnalysisModel set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING HHTExplicit::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check aiPlusOne is of correct size
    if (aiPlusOne.Size() != U->Size())  {
        opserr << "WARNING HHTExplicit::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << aiPlusOne.Size() << endln;
        return -4;
    }
    
    //  determine the response at t+deltaT
    Udot->addVector(1.0, aiPlusOne, c2);
    
    Udotdot->addVector(0.0, aiPlusOne, c3);
    
    // update the response at the DOFs
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "HHTExplicit::update() - failed to update the domain\n";
        return -5;
    }
    // do not update displacements in elements only at nodes
    theModel->setDisp(*U);
    
    return 0;
}


int HHTExplicit::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING HHTExplicit::commit() - no AnalysisModel set\n";
        return -1;
    }
    
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-alpha)*deltaT;
    theModel->setCurrentDomainTime(time);
    
    // update the displacements in the elements
    if (updElemDisp == true)
        theModel->updateDomain();
    
    return theModel->commitDomain();
}

const Vector &
HHTExplicit::getVel()
{
  return *Udot;
}

int HHTExplicit::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(3);
    data(0) = alpha;
    data(1) = gamma;
    if (updElemDisp == false)
        data(2) = 0.0;
    else
        data(2) = 1.0;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTExplicit::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int HHTExplicit::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(3);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTExplicit::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alpha = data(0);
    gamma = data(1);
    if (data(2) == 0.0)
        updElemDisp = false;
    else
        updElemDisp = true;
    
    return 0;
}


void HHTExplicit::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "HHTExplicit - currentTime: " << currentTime << endln;
        s << "  alpha: " << alpha << " gamma: " << gamma << endln;
        s << "  c2: " << c2 << " c3: " << c3 << endln;
        if (updElemDisp)
            s << "  updateElemDisp: yes\n";
        else
            s << "  updateElemDisp: no\n";
    } else
        s << "HHTExplicit - no associated AnalysisModel\n";
}
