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
// $Date: 2009-03-20 18:36:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/TRBDF3.cpp,v $

// Written : krm
// Created : 11/2012
//
// Description: This file contains the implementation of the TRBDF3 class.
// See comments in the header file for what it's doing and relation to TRBDF2.
//
// What: "@(#) TRBDF3.C, revA"

#include <TRBDF3.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <api/runtimeAPI.h>

void *
OPS_ADD_RUNTIME_VPV(OPS_TRBDF3)
{
    return new TRBDF3();
}

TRBDF3::TRBDF3()
  : TransientIntegrator(INTEGRATOR_TAGS_TRBDF3),
    step(0), dt(0.0),
    c1(0.0), c2(0.0), c3(0.0), 
    Utm2(0), Utm2dot(0),
    Utm1(0), Utm1dot(0), 
    Ut(0), Utdot(0), Utdotdot(0), 
    U(0), Udot(0), Udotdot(0)
{
    
}
TRBDF3::~TRBDF3()
{
  // clean up the memory created
    if (Utm2 != 0)
        delete Utm2;
    if (Utm2dot != 0)
        delete Utm2dot;
    
  if (Utm1 != 0)
    delete Utm1;
  if (Utm1dot != 0)
    delete Utm1dot;

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


int TRBDF3::newStep(double deltaT)
{
    // check the vectors have been created
    if (U == 0)  {
        opserr << "TRBDF3::newStep() - domainChange() failed or hasn't been called\n";
        return -3;	
    }

    // mark step as Trapezoidal (=0), Backward Euler (=1), or Houbolt (=2)
    if (deltaT != dt || step == 2) {
        step = 0;
    } else if (step == 0) {
        step = 1;
    } else {
        step = 2;
    }

    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();

    // set response at t to be that at t+deltaT of previous step
    dt = deltaT;

    (*Utm2) = *Utm1;
    (*Utm2dot) = *Utm1dot;

    (*Utm1) = *Ut;
    (*Utm1dot) = *Utdot;

    (*Ut) = *U;        
    (*Utdot) = *Udot;  
    (*Utdotdot) = *Udotdot;

    // set the constants
    if (step == 0)  { // trapezoidal
        c1 = 1.0;
        c2 = 2.0/deltaT;
        c3 = 4.0/(deltaT*deltaT);

        (*Udot) *= -1.0;

        double a3 = -4.0/deltaT;
        double a4 = -1;
        Udotdot->addVector(a4, *Utdot, a3);

    } else if (step == 1)  {  // backward euler
        c1 = 1.0;
        c2 = 1.5/deltaT;
        c3 = 2.25/(deltaT*deltaT);

        (*Udot) = *Utm1;
        Udot->addVector(0.5/deltaT, *Ut, -1/(2.0*deltaT));

        (*Udotdot) = *Utm1dot;
        Udotdot->addVector(0.5/deltaT, *Utdot, -4.0/(2.0*deltaT));    
        Udotdot->addVector(1.0, *Udot, 3.0/(2.0*deltaT));    

    } else {   // Houbolt
        c1 = 1.0;
        c2 = 11/(6.0*deltaT);
        c3 = 2/(deltaT*deltaT);

        (*Udot) = *Utm2;
        Udot->addVector(-1/(3.0*deltaT), *Utm1, 3/(2.0*deltaT));
        Udot->addVector(1.0, *Ut, -7/(6.0*deltaT));

        (*Udotdot) = *Utm2;
        Udotdot->addVector(-1/(deltaT*deltaT), *Utm1, 4/(deltaT*deltaT));
        Udotdot->addVector(1.0, *Ut, -3/(deltaT*deltaT));

    }

    // set the trial response quantities
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);

    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "TRBDF3::newStep() - failed to update the domain\n";
        return -4;
    }

    return 0;
}


int TRBDF3::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next newStep
    if (U != 0)  {
        (*U) = *Ut;        
        (*Udot) = *Utdot;  
        (*Udotdot) = *Utdotdot;  
        // this is crude but restarts the cycle of steps
        step = 2;
    }

    return 0;
}


int TRBDF3::formEleTangent(FE_Element *theEle)
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
    } else if (statusFlag == HALL_TANGENT)  {
        theEle->addKtToTang(c1*cFactor);
        theEle->addKiToTang(c1*iFactor);
        theEle->addCtoTang(c2);
        theEle->addMtoTang(c3);
    } else {
      opserr << "TRBDF3::formEleTangent - unknown FLAG\n";
    }    
    
    return 0;
}    


int TRBDF3::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    theDof->addCtoTang(c2);      
    theDof->addMtoTang(c3);
    
    return 0;
}    


int TRBDF3::domainChanged()
{
    AnalysisModel *myModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    
    // create the new Vector objects
    if (Ut == 0 || Ut->Size() != size)  {
        
        // delete the old
        if (Utm2 != 0)
            delete Utm2;
        if (Utm2dot != 0)
            delete Utm2dot;
        if (Utm1 != 0)
            delete Utm1;
        if (Utm1dot != 0)
            delete Utm1dot;

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
        Utm2 = new Vector(size);
        Utm2dot = new Vector(size);
        Utm1 = new Vector(size);
        Utm1dot = new Vector(size);
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        
        // check we obtained the new
        if (Utm2 == 0 || Utm2->Size() != size ||
            Utm2dot == 0 || Utm2dot->Size() != size ||
            Utm1 == 0 || Utm1->Size() != size ||
            Utm1dot == 0 || Utm1dot->Size() != size ||
            Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size)  {
            
            // delete the old
            if (Utm2 != 0)
                delete Utm2;
            if (Utm2dot != 0)
                delete Utm2dot;
            if (Utm1 != 0)
                delete Utm1;
            if (Utm1dot != 0)
                delete Utm1dot;
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

            Utm2 = 0; Utm2dot = 0; 
            Utm1 = 0; Utm1dot = 0;  
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


int TRBDF3::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING TRBDF3::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING TRBDF3::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING TRBDF3::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
        return -3;
    }
    
    //  determine the response at t+deltaT
    (*U) += deltaU;
    Udot->addVector(1.0, deltaU, c2);
    Udotdot->addVector(1.0, deltaU, c3);

    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);
    if (theModel->updateDomain() < 0)  {
      opserr << "TRBDF3::update() - failed to update the domain\n";
      return -4;
    }
    
    return 0;
}    


const Vector &
TRBDF3::getVel()
{
  return *Udot;
}

int TRBDF3::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to send
    return 0;
}


int TRBDF3::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // nothing to receive
    return 0;
}


void TRBDF3::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t TRBDF3 - currentTime: " << currentTime;
    } else 
        s << "\t TRBDF3 - no associated AnalysisModel\n";
}
