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

// $Revision: 1.5 $
// $Date: 2005-12-19 22:43:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/CentralDifference.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/05
// Revision: A
//
// Description: This file contains the implementation of the CentralDifference class.
//
// What: "@(#) CentralDifference.cpp, revA"

#include <CentralDifference.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


CentralDifference::CentralDifference()
    : TransientIntegrator(INTEGRATOR_TAGS_CentralDifference),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    c2(0.0), c3(0.0),  deltaT(0.0),
    Utm1(0), Ut(0), Utdot(0), Utdotdot(0),
    Udot(0), Udotdot(0)
{
    
}


CentralDifference::CentralDifference(
    double _alphaM, double _betaK, double _betaKi , double _betaKc)
    : TransientIntegrator(INTEGRATOR_TAGS_CentralDifference),
    alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    c2(0.0), c3(0.0),  deltaT(0.0),
    Utm1(0), Ut(0), Utdot(0), Utdotdot(0),
    Udot(0), Udotdot(0)
{
    
}


CentralDifference::~CentralDifference()
{
    // clean up the memory created
    if (Utm1 != 0)
        delete Utm1;
    if (Ut != 0)
        delete Ut;
    if (Utdot != 0)
        delete Utdot;
    if (Utdotdot != 0)
        delete Utdotdot;
    if (Udot != 0)
        delete Udot;
    if (Udotdot != 0)
        delete Udotdot;
}


int CentralDifference::newStep(double _deltaT)
{
    updateCount = 0;
    
    deltaT = _deltaT;
    if (deltaT <= 0.0)  {
        opserr << "CentralDifference::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;	
    }
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    
    // set the constants
    c2 = 0.5/deltaT;
    c3 = 1.0/(deltaT*deltaT);
    
    if (Ut == 0)  {
        opserr << "CentralDifference::newStep() - domainChange() failed or hasn't been called\n";
        return -3;	
    }
    
    // increment the time and apply the load
    double time = theModel->getCurrentDomainTime();
    theModel->applyLoadDomain(time);
    
    // determine the garbage velocities and accelerations at t
    Utdot->addVector(0.0, *Utm1, -c2);
    
    Utdotdot->addVector(0.0, *Ut, -2.0*c3);
    Utdotdot->addVector(1.0, *Utm1, c3);
    
    // set the garbage response quantities for the nodes
    theModel->setVel(*Utdot);
    theModel->setAccel(*Utdotdot);
    
    // set response at t to be that at t+deltaT of previous step
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    return 0;
}


int CentralDifference::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    theEle->addCtoTang(c2);
    theEle->addMtoTang(c3);
    
    return 0;
}    


int CentralDifference::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(c2);
    theDof->addMtoTang(c3);

    return(0);
}    


int CentralDifference::domainChanged()
{
    AnalysisModel *myModel = this->getAnalysisModelPtr();
    LinearSOE *theLinSOE = this->getLinearSOEPtr();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    
    // if damping factors exist set them in the element & node of the domain
    if (alphaM != 0.0 || betaK != 0.0 || betaKi != 0.0 || betaKc != 0.0)
        myModel->setRayleighDampingFactors(alphaM, betaK, betaKi, betaKc);

    // create the new Vector objects
    if (Ut == 0 || Ut->Size() != size)  {
        
        if (Utm1 != 0)
            delete Utm1;
        if (Ut != 0)
            delete Ut;
        if (Utdot != 0)
            delete Utdot;
        if (Utdotdot != 0)
            delete Utdotdot;
        if (Udot != 0)
            delete Udot;
        if (Udotdot != 0)
            delete Udotdot;
        
        // create the new
        Utm1 = new Vector(size);
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        
        // check we obtained the new
        if (Utm1 == 0 || Utm1->Size() != size ||
            Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size)  {
            
            opserr << "CentralDifference::domainChanged - ran out of memory\n";
            
            // delete the old
            if (Utm1 != 0)
                delete Utm1;
            if (Ut != 0)
                delete Ut;
            if (Utdot != 0)
                delete Utdot;
            if (Utdotdot != 0)
                delete Utdotdot;
            if (Udot != 0)
                delete Udot;
            if (Udotdot != 0)
                delete Udotdot;
            
            Utm1 = 0;
            Ut = 0; Utdot = 0; Utdotdot = 0;
            Udot = 0; Udotdot = 0;

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
                (*Utm1)(loc) = disp(i);
                (*Ut)(loc) = disp(i);		
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
    
    opserr << "WARNING: CentralDifference::domainChanged() - assuming Ut-1 = Ut\n";
    
    return 0;
}


int CentralDifference::update(const Vector &U)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING CentralDifference::update() - called more than once -";
        opserr << " CentralDifference integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel == 0)  {
        opserr << "WARNING CentralDifference::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING CentralDifference::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check U is of correct size
    if (U.Size() != Ut->Size()) {
        opserr << "WARNING CentralDifference::update() - Vectors of incompatible size ";
        opserr << " expecting " << Ut->Size() << " obtained " << U.Size() << endln;
        return -3;
    }
    
    //  determine the response at t+deltaT
    Udot->addVector(0.0, U, 3.0);
    Udot->addVector(1.0, *Ut, -4.0);
    Udot->addVector(1.0, *Utm1, 1.0);
    (*Udot) *= c2;
    
    Udotdot->addVector(0.0, *Udot, 1.0);
    Udotdot->addVector(1.0, *Utdot, -1.0);
    (*Udotdot) /= deltaT;
    
    // set response at t to be that at t+deltaT of previous step
    (*Utm1) = *Ut;
    (*Ut) = U;
   
    // update the response at the DOFs
    theModel->setResponse(U, *Udot, *Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "CentralDifference::update() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}    


int CentralDifference::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel == 0) {
        opserr << "WARNING CentralDifference::commit() - no AnalysisModel set\n";
        return -1;
    }	  
    
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    theModel->setCurrentDomainTime(time);
    
    return theModel->commitDomain();
}


int CentralDifference::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = alphaM;
    data(1) = betaK;
    data(2) = betaKi;
    data(3) = betaKc;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING CentralDifference::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int CentralDifference::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING CentralDifference::recvSelf() - could not receive data\n"; 
        return -1;
    }
    
    alphaM = data(0);
    betaK  = data(1);
    betaKi = data(2);
    betaKc = data(3);
    
    return 0;
}


void CentralDifference::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t CentralDifference - currentTime: " << currentTime << endln;
        s << "  Rayleigh Damping - alphaM: " << alphaM;
        s << "  betaK: " << betaK << "  betaKi: " << betaKi << endln;	    
    } else 
        s << "\t CentralDifference - no associated AnalysisModel\n";
}
