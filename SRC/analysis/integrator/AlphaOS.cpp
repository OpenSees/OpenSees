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
// Created: 02/05
// Revision: A
//
// Description: This file contains the implementation of the AlphaOS class.

#include <AlphaOS.h>
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
OPS_ADD_RUNTIME_VPV(OPS_AlphaOS)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc < 1 || argc > 4) {
        opserr << "WARNING - incorrect number of args want AlphaOS $alpha <-updateElemDisp>\n";
        opserr << "          or AlphaOS $alpha $beta $gamma <-updateElemDisp>\n";
        return 0;
    }
    
    bool updElemDisp = false;
    double dData[3];
    int numData;
    if (argc < 3)
        numData = 1;
    else
        numData = 3;
    
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING - invalid args want AlphaOS $alpha <-updateElemDisp>\n";
        opserr << "          or AlphaOS $alpha $beta $gamma <-updateElemDisp>\n";
        return 0;
    }
    
    if (argc == 2 || argc == 4) {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-updateElemDisp") == 0)
            updElemDisp = true;
    }
    
    if (argc < 3)
        theIntegrator = new AlphaOS(dData[0], updElemDisp);
    else
        theIntegrator = new AlphaOS(dData[0], dData[1], dData[2], updElemDisp);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating AlphaOS integrator\n";
    
    return theIntegrator;
}


AlphaOS::AlphaOS()
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS),
    alpha(1.0), beta(0.0), gamma(0.0),
    updElemDisp(0), deltaT(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Upt(0)
{
    
}


AlphaOS::AlphaOS(double _alpha,
    bool updelemdisp)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS),
    alpha(_alpha), beta((2-_alpha)*(2-_alpha)*0.25), gamma(1.5-_alpha),
    updElemDisp(updelemdisp), deltaT(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Upt(0)
{
    
}


AlphaOS::AlphaOS(double _alpha, double _beta, double _gamma,
    bool updelemdisp)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS),
    alpha(_alpha), beta(_beta), gamma(_gamma),
    updElemDisp(updelemdisp), deltaT(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Upt(0)
{
    
}


AlphaOS::~AlphaOS()
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
    if (Upt != 0)
        delete Upt;
}


int AlphaOS::newStep(double _deltaT)
{
    updateCount = 0;
    
    deltaT = _deltaT;
    if (beta == 0 || gamma == 0 )  {
        opserr << "AlphaOS::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "AlphaOS::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << "\n";
        return -2;
    }
    
    // get a pointer to the AnalysisModel
    AnalysisModel *theModel = this->getAnalysisModel();
    
    // set the constants
    c1 = 1.0;
    c2 = gamma/(beta*deltaT);
    c3 = 1.0/(beta*deltaT*deltaT);
    
    if (U == 0)  {
        opserr << "AlphaOS::newStep() - domainChange() failed or hasn't been called\n";
        return -3;
    }
    
    // set response at t to be that at t+deltaT of previous step
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // determine new displacements and velocities at time t+deltaT
    U->addVector(1.0, *Utdot, deltaT);
    double a1 = (0.5 - beta)*deltaT*deltaT;
    U->addVector(1.0, *Utdotdot, a1);
    
    double a2 = deltaT*(1.0 - gamma);
    Udot->addVector(1.0, *Utdotdot, a2);
    
    // determine the response at t+alpha*deltaT
    (*Ualpha) = *Upt;
    Ualpha->addVector((1.0-alpha), *U, alpha);
    
    (*Ualphadot) = *Utdot;
    Ualphadot->addVector((1.0-alpha), *Udot, alpha);
    
    Udotdot->Zero();
    
    // set the trial response quantities
    theModel->setResponse(*Ualpha, *Ualphadot, *Udotdot);
    
    // increment the time to t+alpha*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += alpha*deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "AlphaOS::newStep() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}

const Vector &
AlphaOS::getVel()
{
  return *Udot;
}

int AlphaOS::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int AlphaOS::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    if (statusFlag == CURRENT_TANGENT)
        theEle->addKtToTang(alpha*c1);
    else if (statusFlag == INITIAL_TANGENT)
        theEle->addKiToTang(alpha*c1);
    else if (statusFlag == HALL_TANGENT)  {
        theEle->addKtToTang(alpha*c1*cFactor);
        theEle->addKiToTang(alpha*c1*iFactor);
    }
    
    theEle->addCtoTang(alpha*c2);
    theEle->addMtoTang(c3);
    
    return 0;
}


int AlphaOS::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(alpha*c2);
    theDof->addMtoTang(c3);
    
    return 0;
}


int AlphaOS::domainChanged()
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
        if (Ualpha != 0)
            delete Ualpha;
        if (Ualphadot != 0)
            delete Ualphadot;
        if (Upt != 0)
            delete Upt;
        
        // create the new
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        Ualpha = new Vector(size);
        Ualphadot = new Vector(size);
        Upt = new Vector(size);
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Ualpha == 0 || Ualpha->Size() != size ||
            Ualphadot == 0 || Ualphadot->Size() != size ||
            Upt == 0 || Upt->Size() != size)  {
            
            opserr << "AlphaOS::domainChanged() - ran out of memory\n";
            
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
            if (Upt != 0)
                delete Upt;
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Ualpha = 0; Ualphadot = 0;
            Upt = 0;
            
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
                (*Upt)(loc) = disp(i);
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


int AlphaOS::update(const Vector &deltaU)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING AlphaOS::update() - called more than once -";
        opserr << " AlphaOS integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING AlphaOS::update() - no AnalysisModel set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING AlphaOS::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING AlphaOS::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << "\n";
        return -4;
    }
    
    // save the predictor displacements
    (*Upt) = *U;
    
    //  determine the response at t+deltaT
    U->addVector(1.0, deltaU, c1);
    
    Udot->addVector(1.0, deltaU, c2);
    
    Udotdot->addVector(0.0, deltaU, c3);
    
    // update the response at the DOFs
    theModel->setVel(*Udot);
    theModel->setAccel(*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "AlphaOS::update() - failed to update the domain\n";
        return -5;
    }
    // do not update displacements in elements only at nodes
    theModel->setDisp(*U);
    
    return 0;
}


int AlphaOS::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING AlphaOS::commit() - no AnalysisModel set\n";
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


int AlphaOS::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = alpha;
    data(1) = beta;
    data(2) = gamma;
    if (updElemDisp == false) 
        data(3) = 0.0;
    else
        data(3) = 1.0;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING AlphaOS::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int AlphaOS::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING AlphaOS::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alpha  = data(0);
    beta   = data(1);
    gamma  = data(2);
    if (data(3) == 0.0)
        updElemDisp = false;
    else
        updElemDisp = true;
    
    return 0;
}


void AlphaOS::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "AlphaOS - currentTime: " << currentTime << endln;
        s << "  alpha: " << alpha << "  beta: " << beta  << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
        if (updElemDisp)
            s << "  updateElemDisp: yes\n";
        else
            s << "  updateElemDisp: no\n";
    } else
        s << "AlphaOS - no associated AnalysisModel\n";
}


int AlphaOS::formElementResidual(void)
{
    // calculate Residual Force
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theSOE = this->getLinearSOE();
    
    // loop through the FE_Elements and add the residual
    FE_Element *elePtr;
    FE_EleIter &theEles = theModel->getFEs();
    while((elePtr = theEles()) != 0)  {
        if (theSOE->addB(elePtr->getResidual(this), elePtr->getID()) < 0)  {
            opserr << "WARNING AlphaOS::formElementResidual() -";
            opserr << " failed in addB for ID " << elePtr->getID();
            return -1;
        }
        if (alpha < 1.0)  {
            if (statusFlag == CURRENT_TANGENT)  {
                if (theSOE->addB(elePtr->getK_Force(*Ut-*Upt), elePtr->getID(), alpha-1.0) < 0)  {
                    opserr << "WARNING AlphaOS::formElementResidual() -";
                    opserr << " failed in addB for ID " << elePtr->getID();
                    return -2;
                }
            } else if (statusFlag == INITIAL_TANGENT)  {
                if (theSOE->addB(elePtr->getKi_Force(*Ut-*Upt), elePtr->getID(), alpha-1.0) < 0)  {
                    opserr << "WARNING AlphaOS::formElementResidual() -";
                    opserr << " failed in addB for ID " << elePtr->getID();
                    return -2;
                }
            }
        }
    }
    
    return 0;
}
