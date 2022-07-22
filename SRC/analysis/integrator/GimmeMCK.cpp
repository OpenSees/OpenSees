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

// Written: MHS
// Created: May 2020
// Revision: A
//
// Description: This file contains the implementation of the GimmeMCK class.

#include <GimmeMCK.h>
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


void*  OPS_GimmeMCK(void)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc < 3) {
        opserr << "WARNING - incorrect number of args want GimmeMCK $m $c $k <$ki>\n";
        return 0;
    }

    int numdata = 3;
    double ddata[3];    
    if (OPS_GetDouble(&numdata, ddata) != 0) {
        opserr << "WARNING - invalid args want GimmeMCK $m $c $k <$ki>\n";
        return 0;
    }
    numdata = 1;
    double ki = 0.0;
    if (argc > 3) {
      if (OPS_GetDouble(&numdata, &ki) != 0) {
        opserr << "WARNING - invalid args want GimmeMCK $m $c $k <$ki>\n";
        return 0;
      }
    }
    
    theIntegrator = new GimmeMCK(ddata[0], ddata[1], ddata[2], ki);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating GimmeMCK integrator\n";
    
    return theIntegrator;
}


GimmeMCK::GimmeMCK()
    : TransientIntegrator(INTEGRATOR_TAGS_GimmeMCK),
      m(0.0), c(0.0), k(0.0), ki(0.0), updateCount(0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0)
{
    
}


GimmeMCK::GimmeMCK(double mm, double cc, double kk, double kki)
    : TransientIntegrator(INTEGRATOR_TAGS_GimmeMCK),
      m(mm), c(cc), k(kk), ki(kki), updateCount(0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0)
{
    
}


GimmeMCK::~GimmeMCK()
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


int GimmeMCK::newStep(double deltaT)
{
    updateCount = 0;

    // Allow deltaT == 0.0
    if (deltaT > 0.0)  {
      opserr << "GimmeMCK::newStep() - dT will be ignored\n";
      opserr << "  will use dT=0 and not update the domain" << endln;
      //return -2;
    }
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    
    // set the constants
    //c1 = m;
    //c2 = c;
    //c3 = k;
    
    if (U == 0)  {
        opserr << "GimmeMCK::newStep() - domainChange() failed or hasn't been called\n";
        return -3;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;

    // Do not update response
    // determine new response at time t+deltaT
    //U->addVector(1.0, *Utdot, deltaT);
    //double a1 = 0.5*deltaT*deltaT;
    //U->addVector(1.0, *Utdotdot, a1);
    
    //double a2 = deltaT*(1.0 - gamma);
    //Udot->addVector(1.0, *Utdotdot, a2);
    
    //Udotdot->Zero();
    
    // set the trial response quantities
    theModel->setResponse(*U, *Udot, *Udotdot);

    // Do not update time
    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    //time += deltaT;
    if (theModel->updateDomain(time, 0.0*deltaT) < 0)  {
        opserr << "GimmeMCK::newStep() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int GimmeMCK::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int GimmeMCK::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();

    if (k != 0.0)
      theEle->addKtToTang(k);
    if (ki != 0.0)
      theEle->addKiToTang(ki);
    if (c != 0.0)
      theEle->addCtoTang(c);
    if (m != 0.0)
      theEle->addMtoTang(m);
    
    return 0;
}


int GimmeMCK::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();

    if (c != 0.0)
      theDof->addCtoTang(c);
    if (m != 0.0)
      theDof->addMtoTang(m);
    
    return 0;
}


int GimmeMCK::domainChanged()
{
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theLinSOE = this->getLinearSOE();
    const Vector &x = theLinSOE->getX();
    int size = x.Size();
    
    // create the new Vector objects
    if (U == 0 || U->Size() != size)  {
        
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
            
            opserr << "GimmeMCK::domainChanged() - ran out of memory\n";
            
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


int GimmeMCK::update(const Vector &aiPlusOne)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING GimmeMCK::update() - called more than once -";
        opserr << " GimmeMCK integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING GimmeMCK::update() - no AnalysisModel set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING GimmeMCK::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check aiPlusOne is of correct size
    if (aiPlusOne.Size() != U->Size())  {
        opserr << "WARNING GimmeMCK::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << aiPlusOne.Size() << endln;
        return -4;
    }

    // Don't update the response
    //  determine the response at t+deltaT
    //Udot->addVector(1.0, aiPlusOne, c2);
    
    //Udotdot->addVector(0.0, aiPlusOne, c3);
    
    // update the response at the DOFs
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "GimmeMCK::update() - failed to update the domain\n";
        return -5;
    }
    
    return 0;
}


const Vector &
GimmeMCK::getVel()
{
  return *Udot;
}

int GimmeMCK::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = m;
    data(1) = c;
    data(2) = k;
    data(3) = ki;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING GimmeMCK::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int GimmeMCK::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING GimmeMCK::recvSelf() - could not receive data\n";
        return -1;
    }
    
    m = data(0);
    c = data(1);
    k = data(2);
    ki = data(3);
    
    return 0;
}


void GimmeMCK::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "GimmeMCK - currentTime: " << currentTime << endln;
        s << "  m: " << m << endln;
        s << "  c: " << c << endln;
	s << "  k: " << k << endln;
	s << "  ki: " << ki << endln;
    } else
        s << "GimmeMCK - no associated AnalysisModel\n";
}
