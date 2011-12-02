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

// $Revision: 1.1 $
// $Date: 2005-12-19 22:39:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/NewmarkHybridSimulation.cpp,v $


// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 09/05
// Revision: A
//
// Description: This file contains the implementation of the NewmarkHybridSimulation class.
//
// What: "@(#) NewmarkHybridSimulation.cpp, revA"

#include <NewmarkHybridSimulation.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ConvergenceTest.h>


NewmarkHybridSimulation::NewmarkHybridSimulation()
    : TransientIntegrator(INTEGRATOR_TAGS_NewmarkHybridSimulation),
    gamma(0), beta(0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    theTest(0), rFact(1.0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0)
{
    
}

NewmarkHybridSimulation::NewmarkHybridSimulation(double _gamma, double _beta,
    ConvergenceTest &theT)
    : TransientIntegrator(INTEGRATOR_TAGS_NewmarkHybridSimulation),
    gamma(_gamma), beta(_beta),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    theTest(&theT), rFact(1.0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0)
{

}

NewmarkHybridSimulation::NewmarkHybridSimulation(double _gamma, double _beta,
    ConvergenceTest &theT,
    double _alphaM, double _betaK, double _betaKi, double _betaKc)
    : TransientIntegrator(INTEGRATOR_TAGS_NewmarkHybridSimulation),
    gamma(_gamma), beta(_beta), 
    alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    theTest(&theT), rFact(1.0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0)
{

}

NewmarkHybridSimulation::~NewmarkHybridSimulation()
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

int NewmarkHybridSimulation::newStep(double deltaT)
{
    if (beta == 0 || gamma == 0)  {
        opserr << "NewmarkHybridSimulation::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "NewmarkHybridSimulation::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;	
    }

    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    
    // set the constants
    c1 = 1.0;
    c2 = gamma/(beta*deltaT);
    c3 = 1.0/(beta*deltaT*deltaT);
    
    if (U == 0)  {
        opserr << "NewmarkHybridSimulation::newStep() - domainChange() failed or hasn't been called\n";
        return -3;	
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;        
    (*Utdot) = *Udot;  
    (*Utdotdot) = *Udotdot;
    
    // increment the time and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    theModel->applyLoadDomain(time);
    
    // determine new velocities and accelerations at t+deltaT
    double a1 = (1.0 - gamma/beta); 
    double a2 = (deltaT)*(1.0 - 0.5*gamma/beta);
    Udot->addVector(a1, *Utdotdot, a2);
    
    double a3 = -1.0/(beta*deltaT);
    double a4 = 1.0 - 0.5/beta;
    Udotdot->addVector(a4, *Utdot, a3);
    
    // set the trial response quantities for the nodes
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    
    return 0;
}

int NewmarkHybridSimulation::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next newStep
    if (U != 0)  {
        (*U) = *Ut;        
        (*Udot) = *Utdot;  
        (*Udotdot) = *Utdotdot;  
    }

    return 0;
}

int NewmarkHybridSimulation::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    if (statusFlag == CURRENT_TANGENT)  {
        theEle->addKtToTang(c1);
        theEle->addCtoTang(c2);
        theEle->addMtoTang(c3);
    } else if (statusFlag == INITIAL_TANGENT)  {
        theEle->addKiToTang(c1);
        theEle->addCtoTang(c2);
        theEle->addMtoTang(c3);
    }
    
    return 0;
}    

int NewmarkHybridSimulation::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();

    theDof->addCtoTang(c2);
    theDof->addMtoTang(c3);  
    
    return 0;
}    

int NewmarkHybridSimulation::domainChanged()
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
            Udotdot == 0 || Udotdot->Size() != size)  {
            
            opserr << "NewmarkHybridSimulation::domainChanged - ran out of memory\n";
            
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
    DOF_GrpIter &theDOFs = myModel->getDOFs();
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

int NewmarkHybridSimulation::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel == 0)  {
        opserr << "WARNING NewmarkHybridSimulation::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING NewmarkHybridSimulation::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING NewmarkHybridSimulation::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
        return -3;
    }
    
    // determine the displacement increment reduction factor
    rFact = 1.0/(theTest->getMaxNumTests() - theTest->getNumTests() + 1.0);

    //  determine the response at t+deltaT
    (*U) += rFact*deltaU;
    
    Udot->addVector(1.0, deltaU, rFact*c2);
    
    Udotdot->addVector(1.0, deltaU, rFact*c3);
    
    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);        
    if (theModel->updateDomain() < 0)  {
        opserr << "NewmarkHybridSimulation::update() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}    

int NewmarkHybridSimulation::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(8);
    data(0) = gamma;
    data(1) = beta;
    data(2) = theTest->getClassTag();
    data(3) = theTest->getDbTag();
    data(4) = alphaM;
    data(5) = betaK;
    data(6) = betaKi;
    data(7) = betaKc;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING NewmarkHybridSimulation::sendSelf() - could not send data\n";
        return -1;
    }

    if (theTest->sendSelf(cTag, theChannel) < 0)  {
        opserr << "WARNING NewmarkHybridSimulation::sendSelf() - failed to send CTest object\n";
        return -1;
    }

    return 0;
}

int NewmarkHybridSimulation::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(8);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING NewmarkHybridSimulation::recvSelf() - could not receive data\n";
        gamma = 0.5; beta = 0.25; 
        return -1;
    }
    
    gamma  = data(0);
    beta   = data(1);
    int ctType = int(data(2));
    int ctDb   = int(data(3));
    alphaM = data(4);
    betaK  = data(5);
    betaKi = data(6);
    betaKc = data(7);
    
    theTest = theBroker.getNewConvergenceTest(ctType);
    theTest->setDbTag(ctDb);
    if (theTest->recvSelf(cTag, theChannel, theBroker) < 0) {
        opserr << "WARNING NewmarkHybridSimulation::recvSelf() - failed to recv CTest object\n";
        return -1;
    }
    
    return 0;
}

void NewmarkHybridSimulation::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t NewmarkHybridSimulation - currentTime: " << currentTime;
        s << "  gamma: " << gamma << "  beta: " << beta << endln;
        s << " c1: " << c1 << " c2: " << c2 << " c3: " << c3 << endln;
        s << "  Rayleigh Damping - alphaM: " << alphaM;
        s << "  betaK: " << betaK << "   betaKi: " << betaKi << endln;	    
    } else 
        s << "\t NewmarkHybridSimulation - no associated AnalysisModel\n";
}
