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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/HHTHybridSimulation.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 10/05
// Revision: A
//
// Description: This file contains the implementation of the HHTHybridSimulation class.
//
// What: "@(#) HHTHybridSimulation.cpp, revA"

#include <HHTHybridSimulation.h>
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


HHTHybridSimulation::HHTHybridSimulation()
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHybridSimulation),
    alphaI(1.0), alphaF(1.0),
    beta(0.0), gamma(0.0), deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    theTest(0), rFact(1.0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0)
{

}


HHTHybridSimulation::HHTHybridSimulation(double _rhoInf, ConvergenceTest &theT)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHybridSimulation),
    alphaI((2.0-_rhoInf)/(1.0+_rhoInf)), alphaF(1.0/(1.0+_rhoInf)),
    beta(1.0/(1.0+_rhoInf)/(1.0+_rhoInf)), gamma(0.5*(3.0-_rhoInf)/(1.0+_rhoInf)),
    deltaT(0.0), alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    theTest(&theT), rFact(1.0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0)
{

}


HHTHybridSimulation::HHTHybridSimulation(double _rhoInf, ConvergenceTest &theT,
    double _alphaM, double _betaK, double _betaKi, double _betaKc)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHybridSimulation),
    alphaI((2.0-_rhoInf)/(1.0+_rhoInf)), alphaF(1.0/(1.0+_rhoInf)),
    beta(1.0/(1.0+_rhoInf)/(1.0+_rhoInf)), gamma(0.5*(3.0-_rhoInf)/(1.0+_rhoInf)),
    deltaT(0.0), alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    theTest(&theT), rFact(1.0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0)
{

}


HHTHybridSimulation::HHTHybridSimulation(double _alphaI, double _alphaF,
    double _beta, double _gamma, ConvergenceTest &theT)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHybridSimulation),
    alphaI(_alphaI), alphaF(_alphaF),
    beta(_beta), gamma(_gamma), deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    theTest(&theT), rFact(1.0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0)
{

}


HHTHybridSimulation::HHTHybridSimulation(double _alphaI, double _alphaF,
    double _beta, double _gamma, ConvergenceTest &theT,
    double _alphaM, double _betaK, double _betaKi, double _betaKc)
    : TransientIntegrator(INTEGRATOR_TAGS_HHTHybridSimulation),
    alphaI(_alphaI), alphaF(_alphaF),
    beta(_beta), gamma(_gamma), deltaT(0.0),
    alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    theTest(&theT), rFact(1.0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0)
{

}


HHTHybridSimulation::~HHTHybridSimulation()
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
}


int HHTHybridSimulation::newStep(double _deltaT)
{
    deltaT = _deltaT;
    if (beta == 0 || gamma == 0 )  {
        opserr << "HHTHybridSimulation::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "HHTHybridSimulation::newStep() - error in variable\n";
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
        opserr << "HHTHybridSimulation::newStep() - domainChange() failed or hasn't been called\n";
        return -3;	
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;        
    (*Utdot) = *Udot;  
    (*Utdotdot) = *Udotdot;

    // increment the time to t+alpha*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += alphaF*deltaT;
//    theModel->applyLoadDomain(time);
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "HHTHybridSimulation::newStep() - failed to update the domain\n";
        return -4;
    }

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

    // set the trial response quantities for the nodes
    theModel->setVel(*Ualphadot);
    theModel->setAccel(*Ualphadotdot);
        
    return 0;
}


int HHTHybridSimulation::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }

    return 0;
}


int HHTHybridSimulation::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    if (statusFlag == CURRENT_TANGENT)  {
        theEle->addKtToTang(alphaF*c1);
        theEle->addCtoTang(alphaF*c2);
        theEle->addMtoTang(alphaI*c3);
    } else if (statusFlag == INITIAL_TANGENT)  {
        theEle->addKiToTang(alphaF*c1);
        theEle->addCtoTang(alphaF*c2);
        theEle->addMtoTang(alphaI*c3);
    }
    
    return 0;
}    


int HHTHybridSimulation::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();

    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(alphaI*c3);
    
    return 0;
}


int HHTHybridSimulation::domainChanged()
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
        if (Ualphadotdot != 0)
            delete Ualphadotdot;
        
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
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Ualpha == 0 || Ualpha->Size() != size ||
            Ualphadot == 0 || Ualphadot->Size() != size ||
            Ualphadotdot == 0 || Ualphadotdot->Size() != size)  {
            
            opserr << "HHTHybridSimulation::domainChanged - ran out of memory\n";
            
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
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Ualpha = 0; Ualphadot = 0; Ualphadotdot = 0;

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


int HHTHybridSimulation::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel == 0)  {
        opserr << "WARNING HHTHybridSimulation::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING HHTHybridSimulation::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING HHTHybridSimulation::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
        return -3;
    }
    
    // determine the displacement increment reduction factor
    rFact = 1.0/(theTest->getMaxNumTests() - theTest->getNumTests() + 1.0);

    //  determine the response at t+deltaT
    (*U) += rFact*deltaU;

    Udot->addVector(1.0, deltaU, rFact*c2);
    
    Udotdot->addVector(1.0, deltaU, rFact*c3);

    // determine displacement and velocity at t+alpha*deltaT
    (*Ualpha) = *Ut;
    Ualpha->addVector((1.0-alphaF), *U, alphaF);

    (*Ualphadot) = *Utdot;
    Ualphadot->addVector((1.0-alphaF), *Udot, alphaF);
    
    (*Ualphadotdot) = *Utdotdot;
    Ualphadotdot->addVector((1.0-alphaI), *Udotdot, alphaI);

    // update the response at the DOFs
    theModel->setResponse(*Ualpha,*Ualphadot,*Ualphadotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "HHTHybridSimulation::update() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int HHTHybridSimulation::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel == 0)  {
        opserr << "WARNING HHTHybridSimulation::commit() - no AnalysisModel set\n";
        return -1;
    }	  
    
    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);
//    if (theModel->updateDomain() < 0)  {
//        opserr << "HHTHybridSimulation::commit() - failed to update the domain\n";
//        return -4;
//    }
    
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-alphaF)*deltaT;
    theModel->setCurrentDomainTime(time);

    return theModel->commitDomain();
}


int HHTHybridSimulation::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(10);
    data(0) = alphaI;
    data(1) = alphaF;
    data(2) = beta;
    data(3) = gamma;
    data(4) = theTest->getClassTag();
    data(5) = theTest->getDbTag();
    data(6) = alphaM;
    data(7) = betaK;
    data(8) = betaKi;
    data(9) = betaKc;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTHybridSimulation::sendSelf() - could not send data\n";
        return -1;
    }

    if (theTest->sendSelf(cTag, theChannel) < 0)  {
        opserr << "WARNING HHTHybridSimulation::sendSelf() - failed to send CTest object\n";
        return -1;
    }

    return 0;
}


int HHTHybridSimulation::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(10);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING HHTHybridSimulation::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alphaI = data(0);
    alphaF = data(1);
    beta   = data(2);
    gamma  = data(3);
    int ctType = int(data(4));
    int ctDb   = int(data(5));
    alphaM = data(6);
    betaK  = data(7);
    betaKi = data(8);
    betaKc = data(9);
    
    theTest = theBroker.getNewConvergenceTest(ctType);
    theTest->setDbTag(ctDb);
    if (theTest->recvSelf(cTag, theChannel, theBroker) < 0) {
        opserr << "WARNING HHTHybridSimulation::recvSelf() - failed to recv CTest object\n";
        return -1;
    }

    return 0;
}


void HHTHybridSimulation::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t HHTHybridSimulation - currentTime: " << currentTime << endln;
        s << "  alphaI: " << alphaI << " alphaF: " << alphaF  << " beta: " << beta  << " gamma: " << gamma << endln;
        s << "  c1: " << c1 << " c2: " << c2 << " c3: " << c3 << endln;
        s << "  Rayleigh Damping - alphaM: " << alphaM;
        s << "  betaK: " << betaK << "   betaKi: " << betaKi << endln;	    
    } else 
        s << "\t HHTHybridSimulation - no associated AnalysisModel\n";
}

