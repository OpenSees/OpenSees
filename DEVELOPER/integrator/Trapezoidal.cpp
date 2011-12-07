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
// $Date: 2009/03/20 22:37:09 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewIntegrator/Trapezoidal.cpp,v $

// Written : fmk
// Created : 02/09
//
// Description: This file contains the implementation of the Trapezoidal class.
// ref: K.J.Bathe, "Conserving Energy and Momentum in Nonlinear Dynamics: A Simple
//      Implicit Time Integration Scheme", Computers ans Structures 85(2007),437-445
//
// NOTE: In the implementation we use deltaT as opposed to deltaT/2 in the paper.
//
// What: "@(#) Trapezoidal.C, revA"

#include <elementAPI.h>
#include "Trapezoidal.h"
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

// init function called once on first loading of package

OPS_Export void
localInit() 
{
  OPS_Error("Trapezoidal - Written by fmk UC Berkeley Copyright 2009\n", 1);
}

OPS_Export void *
OPS_Trapezoidal()
{
  //
  // parse args
  //

  // none to parse for this soe

  //
  // create the SOE
  //

  TransientIntegrator *theIntegrator = new Trapezoidal();

  // return it

  return theIntegrator;

}




Trapezoidal::Trapezoidal()
  : TransientIntegrator(0),
    c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), 
    U(0), Udot(0), Udotdot(0)
{
    
}
Trapezoidal::~Trapezoidal()
{
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


int Trapezoidal::newStep(double deltaT)
{
  // check the vectors have been created
  if (U == 0)  {
    opserr << "Trapezoidal::newStep() - domainChange() failed or hasn't been called\n";
    return -3;	
  }

  // get a pointer to the AnalysisModel
  AnalysisModel *theModel = this->getAnalysisModel();

  // set response at t to be that at t+deltaT of previous step
  (*Ut) = *U;        
  (*Utdot) = *Udot;  
  (*Utdotdot) = *Udotdot;

  c1 = 1.0;
  c2 = 2.0/deltaT;
  c3 = 4.0/(deltaT*deltaT);

  (*Udot) *= -1.0;
  
  double a3 = -4.0/deltaT;
  double a4 = -1;
  Udotdot->addVector(a4, *Utdot, a3);
    
  // set the trial response quantities
  theModel->setVel(*Udot);
  theModel->setAccel(*Udotdot);    

  // increment the time to t+deltaT and apply the load
  double time = theModel->getCurrentDomainTime();
  time += deltaT;
  if (theModel->updateDomain(time, deltaT) < 0)  {
    opserr << "Trapezoidal::newStep() - failed to update the domain\n";
    return -4;
  }
  
  return 0;
}


int Trapezoidal::revertToLastStep()
{
  // set response at t+deltaT to be that at t .. for next newStep
  if (U != 0)  {
    (*U) = *Ut;        
    (*Udot) = *Utdot;  
    (*Udotdot) = *Utdotdot;  
  }

  return 0;
}


int Trapezoidal::formEleTangent(FE_Element *theEle)
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


int Trapezoidal::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();

    theDof->addCtoTang(c2);      
    theDof->addMtoTang(c3);
    
    return 0;
}    


int Trapezoidal::domainChanged()
{
    AnalysisModel *myModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    
    // create the new Vector objects
    if (Ut == 0 || Ut->Size() != size)  {
        
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
            
	  // delete any we obtained

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


int Trapezoidal::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING Trapezoidal::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING Trapezoidal::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING Trapezoidal::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
        return -3;
    }
    
    (*U) += deltaU;
    Udot->addVector(1.0, deltaU, c2);
    Udotdot->addVector(1.0, deltaU, c3);

    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);
    if (theModel->updateDomain() < 0)  {
      opserr << "Trapezoidal::update() - failed to update the domain\n";
      return -4;
    }
    
    return 0;
}    


int Trapezoidal::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}


int Trapezoidal::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}


void Trapezoidal::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t Trapezoidal - currentTime: " << currentTime;
    } else 
        s << "\t Trapezoidal - no associated AnalysisModel\n";
}


// AddingSensitivity:BEGIN //////////////////////////////
int Trapezoidal::revertToStart()
{
    if (Ut != 0) 
        Ut->Zero();
    if (Utdot != 0) 
        Utdot->Zero();
    if (Utdotdot != 0) 
        Utdotdot->Zero();

    if (U != 0) 
        U->Zero();
    if (Udot != 0) 
        Udot->Zero();
    if (Udotdot != 0) 
        Udotdot->Zero();
    
    return 0;
}
// AddingSensitivity:END ////////////////////////////////

