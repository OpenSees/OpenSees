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

// $Revision: 1.2 $
// $Date: 2005-12-21 00:32:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/HHTExplicit.cpp,v $


// File: ~/analysis/integrator/HHTExplicit.cpp
// 
// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 10/05
// Revision: A
//
// Description: This file contains the implementation of the HHTExplicit class.
//
// What: "@(#) HHTExplicit.cpp, revA"

#include <HHTExplicit.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


HHTExplicit::HHTExplicit()
    : TransientIntegrator(INTEGRATOR_TAGS_HHTExplicit),
    alpha(1.0), gamma(0.0), deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0)
{
    
}


HHTExplicit::HHTExplicit(double _alpha)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTExplicit),
    alpha(_alpha), gamma(0.5), deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0)
{
    
}


HHTExplicit::HHTExplicit(double _alpha,
    double _alphaM, double _betaK, double _betaKi, double _betaKc)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTExplicit),
    alpha(_alpha), gamma(0.5), deltaT(0.0),
    alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0)
{
    
}


HHTExplicit::HHTExplicit(double _alpha, double _gamma)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTExplicit),
    alpha(_alpha), gamma(_gamma), deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0)
{
    
}


HHTExplicit::HHTExplicit(double _alpha, double _gamma,
    double _alphaM, double _betaK, double _betaKi, double _betaKc)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTExplicit),
    alpha(_alpha), gamma(_gamma), deltaT(0.0),
    alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    c2(0.0), c3(0.0), 
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

    deltaT = _deltaT;
    if (gamma == 0 )  {
        opserr << "HHTExplicit::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << endln;
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "HHTExplicit::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;	
    }
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModelPtr();

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

    // determine the displacements and velocities at t+alpha*deltaT
    (*Ualpha) = *Ut;
    Ualpha->addVector((1.0-alpha), *U, alpha);
    
    (*Ualphadot) = *Utdot;
    Ualphadot->addVector((1.0-alpha), *Udot, alpha);
        
    // set the trial response quantities for the elements
    theModel->setDisp(*Ualpha);
    theModel->setVel(*Ualphadot);

    // increment the time to t+alpha*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += alpha*deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "HHTExplicit::newStep() - failed to update the domain\n";
        return -4;
    }
    
    // determine the accelerations at t+alpha*deltaT
    Udotdot->Zero();
    
    // set the new trial response quantities for the nodes
    theModel->setAccel(*Udotdot);
    
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
    AnalysisModel *myModel = this->getAnalysisModelPtr();
    LinearSOE *theLinSOE = this->getLinearSOEPtr();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    
    // if damping factors exist set them in the ele & node of the domain
    if (alphaM != 0.0 || betaK != 0.0 || betaKi != 0.0 || betaKc != 0.0)
        myModel->setRayleighDampingFactors(alphaM, betaK, betaKi, betaKc);
    
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
            
            opserr << "HHTExplicit::domainChanged - ran out of memory\n";
            
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
    DOF_GrpIter &theDOFs = myModel->getDOFs();
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

    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel == 0)  {
        opserr << "WARNING HHTExplicit::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING HHTExplicit::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check aiPlusOne is of correct size
    if (aiPlusOne.Size() != U->Size())  {
        opserr << "WARNING HHTExplicit::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << aiPlusOne.Size() << endln;
        return -3;
    }
    
    //  determine the response at t+deltaT
    Udot->addVector(1.0, aiPlusOne, c2);

    (*Udotdot) = aiPlusOne;
    
    // update the response at the DOFs
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    
    return 0;
}


int HHTExplicit::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel == 0)  {
        opserr << "WARNING HHTExplicit::commit() - no AnalysisModel set\n";
        return -1;
    }	  
        
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-alpha)*deltaT;
    theModel->setCurrentDomainTime(time);

    return theModel->commitDomain();
}


int HHTExplicit::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(6);
    data(0) = alpha;
    data(1) = gamma;
    data(2) = alphaM;
    data(3) = betaK;
    data(4) = betaKi;
    data(5) = betaKc;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTExplicit::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int HHTExplicit::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(6);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTExplicit::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alpha  = data(0);
    gamma  = data(1);
    alphaM = data(2);
    betaK  = data(3);
    betaKi = data(4);
    betaKc = data(5);
    
    return 0;
}


void HHTExplicit::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t HHTExplicit - currentTime: " << currentTime << endln ;
        s << "  alpha: " << alpha << " gamma: " << gamma << endln;
        s << "  c2: " << c2 << " c3: " << c3 << endln;
        s << "  Rayleigh Damping - alphaM: " << alphaM;
        s << "  betaK: " << betaK << "   betaKi: " << betaKi << endln;	    
    } else 
        s << "\t HHTExplicit - no associated AnalysisModel\n";
}



