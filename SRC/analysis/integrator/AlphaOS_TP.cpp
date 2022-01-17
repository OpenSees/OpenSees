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
// Description: This file contains the implementation of the AlphaOS_TP class.

#include <AlphaOS_TP.h>
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
OPS_ADD_RUNTIME_VPV(OPS_AlphaOS_TP)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    int argc = OPS_GetNumRemainingInputArgs();
    if (argc < 1 || argc > 4) {
        opserr << "WARNING - incorrect number of args want AlphaOS_TP $alpha <-updateElemDisp>\n";
        opserr << "          or AlphaOS_TP $alpha $beta $gamma <-updateElemDisp>\n";
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
        opserr << "WARNING - invalid args want AlphaOS_TP $alpha <-updateElemDisp>\n";
        opserr << "          or AlphaOS_TP $alpha $beta $gamma <-updateElemDisp>\n";
        return 0;
    }
    
    if (argc == 2 || argc == 4) {
        const char *argvLoc = OPS_GetString();
        if (strcmp(argvLoc, "-updateElemDisp") == 0)
            updElemDisp = true;
    }
    
    if (argc < 3)
        theIntegrator = new AlphaOS_TP(dData[0], updElemDisp);
    else
        theIntegrator = new AlphaOS_TP(dData[0], dData[1], dData[2], updElemDisp);
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating AlphaOS_TP integrator\n";
    
    return theIntegrator;
}


AlphaOS_TP::AlphaOS_TP()
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS_TP),
    alpha(1.0), beta(0.0), gamma(0.0),
    updElemDisp(0), deltaT(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    alphaD(1.0), alphaR(1.0), alphaKU(0.0), alphaP(1.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Upt(0), Put(0)
{
    
}


AlphaOS_TP::AlphaOS_TP(double _alpha,
    bool updelemdisp)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS_TP),
    alpha(_alpha), beta((2-_alpha)*(2-_alpha)*0.25), gamma(1.5-_alpha),
    updElemDisp(updelemdisp), deltaT(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    alphaD(alpha), alphaR(alpha), alphaKU(0.0), alphaP(alpha),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Upt(0), Put(0)
{
    
}


AlphaOS_TP::AlphaOS_TP(double _alpha, double _beta, double _gamma,
    bool updelemdisp)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS_TP),
    alpha(_alpha), beta(_beta), gamma(_gamma),
    updElemDisp(updelemdisp), deltaT(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    alphaD(alpha), alphaR(alpha), alphaKU(0.0), alphaP(alpha),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Upt(0), Put(0)
{
    
}


AlphaOS_TP::~AlphaOS_TP()
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
    if (Upt != 0)
        delete Upt;
    if (Put != 0)
        delete Put;
}


int AlphaOS_TP::newStep(double _deltaT)
{
    updateCount = 0;
    
    if (beta == 0 || gamma == 0 )  {
        opserr << "AlphaOS_TP::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    deltaT = _deltaT;
    if (deltaT <= 0.0)  {
        opserr << "AlphaOS_TP::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << "\n";
        return -2;
    }
    
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING AlphaOS_TP::newStep() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -3;
    }
    
    // set the constants
    c1 = 1.0;
    c2 = gamma/(beta*deltaT);
    c3 = 1.0/(beta*deltaT*deltaT);
    
    if (U == 0)  {
        opserr << "AlphaOS_TP::newStep() - domainChange() failed or hasn't been called\n";
        return -4;
    }
    
    // set weighting factors for subsequent iterations
    alphaD = alphaR = alphaP = alpha;
    alphaKU = 0.0;
    
    // determine new response at time t+deltaT
    U->addVector(1.0, *Utdot, deltaT);
    double a1 = (0.5 - beta)*deltaT*deltaT;
    U->addVector(1.0, *Utdotdot, a1);
    
    double a2 = deltaT*(1.0 - gamma);
    Udot->addVector(1.0, *Utdotdot, a2);
    
    // set the trial response quantities
    theModel->setDisp(*U);
    theModel->setVel(*Udot);
    
    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "AlphaOS_TP::newStep() - failed to update the domain\n";
        return -5;
    }
    
    return 0;
}

const Vector &
AlphaOS_TP::getVel()
{
  return *Udot;
}

int AlphaOS_TP::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }
    
    return 0;
}


int AlphaOS_TP::formUnbalance()
{
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING AlphaOS_TP::formUnbalance() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -1;
    }
    
    theLinSOE->setB(*Put);
    
    // do modal damping
    const Vector *modalValues = theModel->getModalDampingFactors();
    if (modalValues != 0)  {
        this->addModalDampingForce(modalValues);
    }
    
    if (this->formElementResidual() < 0)  {
        opserr << "WARNING AlphaOS_TP::formUnbalance() ";
        opserr << " - this->formElementResidual failed\n";
        return -2;
    }
    
    if (this->formNodalUnbalance() < 0)  {
        opserr << "WARNING AlphaOS_TP::formUnbalance() ";
        opserr << " - this->formNodalUnbalance failed\n";
        return -3;
    }
    
    return 0;
}


int AlphaOS_TP::formEleTangent(FE_Element *theEle)
{
    theEle->zeroTangent();
    
    if (statusFlag == CURRENT_TANGENT)
        theEle->addKtToTang(alpha*c1);
    else if (statusFlag == INITIAL_TANGENT)
        theEle->addKiToTang(alpha*c1);
    
    theEle->addCtoTang(alpha*c2);
    theEle->addMtoTang(c3);
    
    return 0;
}


int AlphaOS_TP::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(alpha*c2);
    theDof->addMtoTang(c3);
    
    return 0;
}


int AlphaOS_TP::formEleResidual(FE_Element *theEle)
{
    theEle->zeroResidual();
    
    // this does not work because for some elements damping is returned
    // with the residual as well as the damping tangent 
    //theEle->addRtoResidual(alphaR);
    //theEle->addD_Force(*Udot, -alphaD);
    //theEle->addM_Force(*Udotdot, -alphaM);
    
    // instead use residual including the inertia terms and then correct
    // the mass contribution (only works because alphaR = alphaD) 
    theEle->addRIncInertiaToResidual(alphaR);
    theEle->addM_Force(*Udotdot, alphaR);
    
    return 0;
}


int AlphaOS_TP::formNodUnbalance(DOF_Group *theDof)
{
    theDof->zeroUnbalance();
    
    theDof->addPtoUnbalance(alphaP);
    theDof->addD_Force(*Udot, -alphaD);
    
    return 0;
}


int AlphaOS_TP::domainChanged()
{
    AnalysisModel *theModel = this->getAnalysisModel();
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
        if (Upt != 0)
            delete Upt;
        if (Put != 0)
            delete Put;
        
        // create the new
        Ut = new Vector(size);
        Utdot = new Vector(size);
        Utdotdot = new Vector(size);
        U = new Vector(size);
        Udot = new Vector(size);
        Udotdot = new Vector(size);
        Upt = new Vector(size);
        Put = new Vector(size);
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Upt == 0 || Upt->Size() != size ||
            Put == 0 || Put->Size() != size)  {
            
            opserr << "AlphaOS_TP::domainChanged() - ran out of memory\n";
            
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
            if (Upt != 0)
                delete Upt;
            if (Put != 0)
                delete Put;
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Upt = 0; Put = 0;
            
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
    
    // now get unbalance at last commit and store it
    // warning: this will use committed stiffness prop. damping
    // from current step instead of previous step
    alphaD = alphaR = alphaKU = alphaP = (1.0 - alpha);
    if (alpha < 1.0)  {
        this->TransientIntegrator::formUnbalance();
        (*Put) = theLinSOE->getB();
    } else {
        Put->Zero();
    }
    
    return 0;
}


int AlphaOS_TP::update(const Vector &deltaU)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING AlphaOS_TP::update() - called more than once -";
        opserr << " AlphaOS_TP integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING AlphaOS_TP::update() - no AnalysisModel set\n";
        return -2;
    }
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING AlphaOS_TP::update() - domainChange() failed or not called\n";
        return -3;
    }
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING AlphaOS_TP::update() - Vectors of incompatible size ";
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
        opserr << "AlphaOS_TP::update() - failed to update the domain\n";
        return -5;
    }
    // do not update displacements in elements only at nodes
    theModel->setDisp(*U);
    
    return 0;
}


int AlphaOS_TP::commit(void)
{
    // get a pointer to the LinearSOE and the AnalysisModel
    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0)  {
        opserr << "WARNING AlphaOS_TP::commit() - ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -1;
    }
    
    // set response at t of next step to be that at t+deltaT
    (*Ut) = *U;
    (*Utdot) = *Udot;
    (*Utdotdot) = *Udotdot;
    
    // get unbalance Put and store it for next step
    alphaD = alphaR = alphaKU = alphaP = (1.0 - alpha);
    if (alpha < 1.0)  {
        this->TransientIntegrator::formUnbalance();
        (*Put) = theLinSOE->getB();
    } else {
        Put->Zero();
    }
    
    // update the displacements in the elements
    if (updElemDisp == true)
        theModel->updateDomain();
    
    return theModel->commitDomain();
}


int AlphaOS_TP::sendSelf(int cTag, Channel &theChannel)
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
        opserr << "WARNING AlphaOS_TP::sendSelf() - could not send data\n";
        return -1;
    }
    
    return 0;
}


int AlphaOS_TP::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING AlphaOS_TP::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alpha  = data(0);
    beta   = data(1);
    gamma  = data(2);
    if (data(3) == 0.0)
        updElemDisp = false;
    else
        updElemDisp = true;
    
    alphaD  = alpha;
    alphaR  = alpha;
    alphaKU = 0.0;
    alphaP  = alpha;
    
    return 0;
}


void AlphaOS_TP::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "AlphaOS_TP - currentTime: " << currentTime << endln;
        s << "  alpha: " << alpha << "  beta: " << beta  << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
        if (updElemDisp)
            s << "  updateElemDisp: yes\n";
        else
            s << "  updateElemDisp: no\n";
    } else 
        s << "AlphaOS_TP - no associated AnalysisModel\n";
}


int AlphaOS_TP::formElementResidual(void)
{
    // calculate Residual Force
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theSOE = this->getLinearSOE();
    
    // loop through the FE_Elements and add the residual
    FE_Element *elePtr;
    FE_EleIter &theEles = theModel->getFEs();
    while((elePtr = theEles()) != 0)  {
        if (theSOE->addB(elePtr->getResidual(this), elePtr->getID()) < 0)  {
            opserr << "WARNING AlphaOS_TP::formElementResidual() -";
            opserr << " failed in addB for ID " << elePtr->getID();
            return -1;
        }
        if (alphaKU > 0.0)  {
            if (statusFlag == CURRENT_TANGENT)  {
                if (theSOE->addB(elePtr->getK_Force(*Ut-*Upt), elePtr->getID(), -alphaKU) < 0)  {
                    opserr << "WARNING AlphaOS_TP::formElementResidual() -";
                    opserr << " failed in addB for ID " << elePtr->getID();
                    return -2;
                }
            } else if (statusFlag == INITIAL_TANGENT)  {
                if (theSOE->addB(elePtr->getKi_Force(*Ut-*Upt), elePtr->getID(), -alphaKU) < 0)  {
                    opserr << "WARNING AlphaOS_TP::formElementResidual() -";
                    opserr << " failed in addB for ID " << elePtr->getID();
                    return -2;
                }
            }
        }
    }
    
    return 0;
}
