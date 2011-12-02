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

// $Revision: 1.19 $
// $Date: 2010-02-04 01:06:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/Newmark.cpp,v $

// Written : fmk
// Created : 11/98
// Modified: 02/05 ahs
// Revision: A
//
// Description: This file contains the implementation of the Newmark class.
//
// What: "@(#) Newmark.C, revA"

#include <Newmark.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <elementAPI.h>

TransientIntegrator *
OPS_NewNewmark(void)
{
  // Pointer to a uniaxial material that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want Newmark $gamma $beta <-form $typeUnknown>\n";
    return 0;
  }

  bool dispFlag = true;
  double dData[2];
  int numData = 2;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want Newmark $gamma $beta <-form $typeUnknown>\n";
    return 0;
  }
  
  if (argc == 2)
    theIntegrator = new Newmark(dData[0], dData[1]);
  else {
    char nextString[10];
    OPS_GetString(nextString, 10);
    if (strcmp(nextString,"-form") == 0) {
      OPS_GetString(nextString, 10);
      if ((nextString[0] == 'D') || (nextString[0] == 'd')) 
	dispFlag = true;
      else if ((nextString[0] == 'A') || (nextString[0] == 'a')) 
	dispFlag = false;      
    }    
    theIntegrator = new Newmark(dData[0], dData[1], dispFlag);
  }

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating Newmark integrator\n";

  return theIntegrator;
}


Newmark::Newmark()
    : TransientIntegrator(INTEGRATOR_TAGS_Newmark),
    displ(true), gamma(0), beta(0), 
    c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    determiningMass(false)
{
    
}


Newmark::Newmark(double _gamma, double _beta, bool dispFlag)
    : TransientIntegrator(INTEGRATOR_TAGS_Newmark),
    displ(dispFlag), gamma(_gamma), beta(_beta), 
    c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    determiningMass(false)
{
    
}

Newmark::~Newmark()
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


int Newmark::newStep(double deltaT)
{
    if (beta == 0 || gamma == 0)  {
        opserr << "Newmark::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "Newmark::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;	
    }

    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    
    // set the constants
    if (displ == true)  {
        c1 = 1.0;
        c2 = gamma/(beta*deltaT);
        c3 = 1.0/(beta*deltaT*deltaT);
    } else  {
        c1 = beta*deltaT*deltaT;
        c2 = gamma*deltaT;
        c3 = 1.0;
    }
    
    if (U == 0)  {
        opserr << "Newmark::newStep() - domainChange() failed or hasn't been called\n";
        return -3;	
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;        
    (*Utdot) = *Udot;  
    (*Utdotdot) = *Udotdot;
    
    if (displ == true)  {    
        // determine new velocities and accelerations at t+deltaT
        double a1 = (1.0 - gamma/beta); 
        double a2 = (deltaT)*(1.0 - 0.5*gamma/beta);
        Udot->addVector(a1, *Utdotdot, a2);
        
        double a3 = -1.0/(beta*deltaT);
        double a4 = 1.0 - 0.5/beta;
        Udotdot->addVector(a4, *Utdot, a3);

        // set the trial response quantities
        theModel->setVel(*Udot);
        theModel->setAccel(*Udotdot);
    } else  {
        // determine new displacements and velocities at t+deltaT      
        double a1 = (deltaT*deltaT/2.0);
        U->addVector(1.0, *Utdot, deltaT);
        U->addVector(1.0, *Utdotdot, a1);
        
        Udot->addVector(1.0, *Utdotdot, deltaT);

        // set the trial response quantities
        theModel->setDisp(*U);
        theModel->setVel(*Udot);
    }

    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "Newmark::newStep() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int Newmark::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next newStep
    if (U != 0)  {
        (*U) = *Ut;        
        (*Udot) = *Utdot;  
        (*Udotdot) = *Utdotdot;  
    }

    return 0;
}


int Newmark::formEleTangent(FE_Element *theEle)
{
    if (determiningMass == true)
        return 0;
    
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


int Newmark::formNodTangent(DOF_Group *theDof)
{
    if (determiningMass == true)
        return 0;
    
    theDof->zeroTangent();

    theDof->addCtoTang(c2);      
    theDof->addMtoTang(c3);
    
    return 0;
}    


int Newmark::domainChanged()
{
    AnalysisModel *myModel = this->getAnalysisModel();
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
            Udotdot == 0 || Udotdot->Size() != size)  {
            
            opserr << "Newmark::domainChanged - ran out of memory\n";
            
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


int Newmark::update(const Vector &deltaU)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING Newmark::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING Newmark::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING Newmark::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
        return -3;
    }
    
    //  determine the response at t+deltaT
    if (displ == true)  {
        (*U) += deltaU;

        Udot->addVector(1.0, deltaU, c2);

        Udotdot->addVector(1.0, deltaU, c3);
    } else  {
        U->addVector(1.0, deltaU, c1);
        
        Udot->addVector(1.0, deltaU, c2);
        
        (*Udotdot) += deltaU;
    }
    
    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "Newmark::update() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}    


int Newmark::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(3);
    data(0) = gamma;
    data(1) = beta;
    if (displ == true) 
        data(2) = 1.0;
    else
        data(2) = 0.0;

    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING Newmark::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int Newmark::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(3);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING Newmark::recvSelf() - could not receive data\n";
        gamma = 0.5; beta = 0.25; 
        return -1;
    }
    
    gamma  = data(0);
    beta   = data(1);
    if (data(2) == 1.0)
        displ = true;
    else
        displ = false;

    return 0;
}


void Newmark::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t Newmark - currentTime: " << currentTime;
        s << "  gamma: " << gamma << "  beta: " << beta << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else 
        s << "\t Newmark - no associated AnalysisModel\n";
}


// AddingSensitivity:BEGIN //////////////////////////////
int Newmark::revertToStart()
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

