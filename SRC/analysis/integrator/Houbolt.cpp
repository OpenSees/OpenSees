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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/Houbolt.cpp,v $

// Written : krm
// Created : 11/2012
//
// Description: This file contains the implementation of the Houbolt class.
// See comments in the header file for description of algorithm.
//
// What: "@(#) Houbolt.C, revA"

#include <Houbolt.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h> // OPS_ADD_RUNTIME_VPV

void *
OPS_ADD_RUNTIME_VPV(OPS_Houbolt)
{
    return new Houbolt();
}

Houbolt::Houbolt()
  : TransientIntegrator(INTEGRATOR_TAGS_Houbolt),
    step(0), dt(0.0), 
    c1(0.0), c2(0.0), c3(0.0), 
    Utm2(0), Utm1(0),  
    Ut(0), Utdot(0), Utdotdot(0), 
    U(0), Udot(0), Udotdot(0)
{
    
    
}
Houbolt::~Houbolt()
{
  // clean up the memory created
    if (Utm2 != 0)
        delete Utm2;
    if (Utm1 != 0)
        delete Utm1;

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


int Houbolt::newStep(double deltaT)
{
  // check the vectors have been created
  if (U == 0)  {
    opserr << "Houbolt::newStep() - domainChange() failed or hasn't been called\n";
    return -3;	
  }

  // mark step as bootstrap or not
  if ( deltaT != dt )
      step = 0;
  else
      step++;

  // get a pointer to the AnalysisModel
  AnalysisModel *theModel = this->getAnalysisModel();

  // set response at t to be that at t+deltaT of previous step
  dt = deltaT;

  // save storage quantites before updating
  (*Utm2) = *Utm1;
  (*Utm1) = *Ut;

  (*Ut) = *U;        
  (*Utdot) = *Udot;  
  (*Utdotdot) = *Udotdot;

  // set the constants
  if (step < 2)  { // trapezoidal
    c1 = 1.0;
    c2 = 2.0/deltaT;
    c3 = 4.0/(deltaT*deltaT);

    (*Udot) *= -1.0;

    double a3 = -4.0/deltaT;
    double a4 = -1;
    Udotdot->addVector(a4, *Utdot, a3);
    
  } else {   // Houbolt
      c1 = 1.0;
      c2 = 11.0/(6.0*deltaT);
      c3 = 2.0/(deltaT*deltaT);
      
      (*Udot) = *Utm2;
      Udot->addVector(-1.0/(3.0*deltaT), *Utm1, 3.0/(2.0*deltaT));
      Udot->addVector(1.0, *Ut, -7.0/(6.0*deltaT));
      
      (*Udotdot) = *Utm2;
      Udotdot->addVector(-1.0/(deltaT*deltaT), *Utm1, 4.0/(deltaT*deltaT));
      Udotdot->addVector(1.0, *Ut, -3.0/(deltaT*deltaT));
      
  }
    
    // set the trial response quantities
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);

    // increment the time to t+deltaT and apply the load
  double time = theModel->getCurrentDomainTime();
  time += deltaT;
  if (theModel->updateDomain(time, deltaT) < 0)  {
    opserr << "Houbolt::newStep() - failed to update the domain\n";
    return -4;
  }
  
  return 0;
}


int Houbolt::revertToLastStep()
{
  // set response at t+deltaT to be that at t .. for next newStep
  if (U != 0)  {
    (*U) = *Ut;        
    (*Udot) = *Utdot;  
    (*Udotdot) = *Utdotdot;
      
      // this is crude, should be based on previous steps taken
    step = 0;
  }

  return 0;
}


int Houbolt::formEleTangent(FE_Element *theEle)
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
      opserr << "Houbold::formEleTangent - unknown FLAG\n";
    }   
    
    return 0;
}    


int Houbolt::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    theDof->addCtoTang(c2);      
    theDof->addMtoTang(c3);
    
    return 0;
}    


int Houbolt::domainChanged()
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
        if (Utm1 != 0)
            delete Utm1;

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
        Utm1 = new Vector(size);
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        
        // check we obtained the new
        if (Utm2 == 0 || Utm2->Size() != size ||
            Utm1 == 0 || Utm1->Size() != size ||
            Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size)  {
            
            // delete the old
            if (Utm2 != 0)
                delete Utm2;
            if (Utm1 != 0)
                delete Utm1;
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

            Utm2 = 0; Utm1 = 0;
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


int Houbolt::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING Houbolt::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING Houbolt::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING Houbolt::update() - Vectors of incompatible size ";
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
      opserr << "Houbolt::update() - failed to update the domain\n";
      return -4;
    }
    
    return 0;
}    

const Vector &
Houbolt::getVel()
{
  return *Udot;
}

int Houbolt::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to send
    return 0;
}


int Houbolt::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // nothing to receive
    return 0;
}


void Houbolt::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t Houbolt - currentTime: " << currentTime;
    } else 
        s << "\t Houbolt - no associated AnalysisModel\n";
}

