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
// $Date: 2009-05-19 22:10:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/AlphaOSGeneralized.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 10/05
// Revision: A
//
// Description: This file contains the implementation of the AlphaOSGeneralized class.
//
// What: "@(#) AlphaOSGeneralized.cpp, revA"

#include <AlphaOSGeneralized.h>
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


AlphaOSGeneralized::AlphaOSGeneralized()
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOSGeneralized),
    alphaI(1.0), alphaF(1.0), beta(0.0), gamma(0.0),
    updDomFlag(0), deltaT(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0), Upt(0), Uptdot(0)
{
    
}


AlphaOSGeneralized::AlphaOSGeneralized(double _rhoInf,
    bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOSGeneralized),
    alphaI((2.0-_rhoInf)/(1.0+_rhoInf)), alphaF(1.0/(1.0+_rhoInf)),
    beta(1.0/(1.0+_rhoInf)/(1.0+_rhoInf)), gamma(0.5*(3.0-_rhoInf)/(1.0+_rhoInf)),
    updDomFlag(upddomflag), deltaT(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0), Upt(0), Uptdot(0)
{
    
}


AlphaOSGeneralized::AlphaOSGeneralized(double _alphaI, double _alphaF,
    double _beta, double _gamma,
    bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOSGeneralized),
    alphaI(_alphaI), alphaF(_alphaF), beta(_beta), gamma(_gamma),
    updDomFlag(upddomflag), deltaT(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Ualphadotdot(0), Upt(0), Uptdot(0)
{
    
}


AlphaOSGeneralized::~AlphaOSGeneralized()
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
    if (Upt != 0)
        delete Upt;
    if (Uptdot != 0)
        delete Uptdot;
}


int AlphaOSGeneralized::newStep(double _deltaT)
{
    updateCount = 0;

    deltaT = _deltaT;
    if (beta == 0 || gamma == 0 )  {
        opserr << "AlphaOSGeneralized::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }
    
    if (deltaT <= 0.0)  {
        opserr << "AlphaOSGeneralized::newStep() - error in variable\n";
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
        opserr << "AlphaOSGeneralized::newStep() - domainChange() failed or hasn't been called\n";
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
            
    // determine the response at t+alphaF*deltaT
    (*Ualpha) = *Upt;
    Ualpha->addVector((1.0-alphaF), *U, alphaF);
    
    (*Ualphadot) = *Utdot;
    Ualphadot->addVector((1.0-alphaF), *Udot, alphaF);

    (*Ualphadotdot) = (1.0-alphaI)*(*Utdotdot);
        
    // set the trial response quantities
    theModel->setResponse(*Ualpha,*Ualphadot,*Ualphadotdot);
    
    // increment the time to t+alphaF*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += alphaF*deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "AlphaOSGeneralized::newStep() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}


int AlphaOSGeneralized::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next step
    if (U != 0)  {
        (*U) = *Ut;
        (*Udot) = *Utdot;
        (*Udotdot) = *Utdotdot;
    }

    return 0;
}


int AlphaOSGeneralized::formEleTangent(FE_Element *theEle)
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


int AlphaOSGeneralized::formNodTangent(DOF_Group *theDof)
{
    theDof->zeroTangent();
    
    theDof->addCtoTang(alphaF*c2);
    theDof->addMtoTang(alphaI*c3);
    
    return 0;
}    


int AlphaOSGeneralized::domainChanged()
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
        if (Ualphadotdot != 0)
            delete Ualphadotdot;
        if (Upt != 0)
            delete Upt;
        if (Uptdot != 0)
            delete Uptdot;
        
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
        Upt = new Vector(size);
        Uptdot = new Vector(size);
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size ||
            Ualpha == 0 || Ualpha->Size() != size ||
            Ualphadot == 0 || Ualphadot->Size() != size ||
            Ualphadotdot == 0 || Ualphadotdot->Size() != size ||
            Upt == 0 || Upt->Size() != size ||
            Uptdot == 0 || Uptdot->Size() != size)  {
            
            opserr << "AlphaOSGeneralized::domainChanged - ran out of memory\n";
            
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
            if (Upt != 0)
                delete Upt;
            if (Uptdot != 0)
                delete Uptdot;
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Ualpha = 0; Ualphadot = 0; Ualphadotdot = 0;
            Upt = 0; Uptdot = 0;

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


int AlphaOSGeneralized::update(const Vector &deltaU)
{
    updateCount++;
    if (updateCount > 1)  {
        opserr << "WARNING AlphaOSGeneralized::update() - called more than once -";
        opserr << " AlphaOSGeneralized integration scheme requires a LINEAR solution algorithm\n";
        return -1;
    }
    
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING AlphaOSGeneralized::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING AlphaOSGeneralized::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING AlphaOSGeneralized::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << "\n";
        return -3;
    }
    
    // save the predictor displacements and velocities
    (*Upt) = *U;
    (*Uptdot) = *Udot;

    //  determine the response at t+deltaT
    (*U) += deltaU;

    Udot->addVector(1.0, deltaU, c2);

    (*Udotdot) = c3*deltaU;
    
    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);
    if (updDomFlag == true)  {
        if (theModel->updateDomain() < 0)  {
            opserr << "AlphaOSGeneralized::update() - failed to update the domain\n";
            return -4;
        }
    }
    
    return 0;
}    


int AlphaOSGeneralized::commit(void)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING AlphaOSGeneralized::commit() - no AnalysisModel set\n";
        return -1;
    }
        
    // set the time to be t+deltaT
    double time = theModel->getCurrentDomainTime();
    time += (1.0-alphaF)*deltaT;
    theModel->setCurrentDomainTime(time);
    
    return theModel->commitDomain();
}    


int AlphaOSGeneralized::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(5);
    data(0) = alphaI;
    data(1) = alphaF;
    data(2) = beta;
    data(3) = gamma;
    if (updDomFlag == false) 
        data(4) = 0.0;
    else
        data(4) = 1.0;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING AlphaOSGeneralized::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int AlphaOSGeneralized::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(5);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING AlphaOSGeneralized::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alphaI = data(0);
    alphaF = data(1);
    beta   = data(2);
    gamma  = data(3);
    if (data(4) == 0.0)
        updDomFlag = false;
    else
        updDomFlag = true;
    
    return 0;
}


void AlphaOSGeneralized::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0)  {
        double currentTime = theModel->getCurrentDomainTime();
        s << "AlphaOSGeneralized - currentTime: " << currentTime << endln;
        s << "  alphaI: " << alphaI << "  alphaF: " << alphaF  << "  beta: " << beta  << "  gamma: " << gamma << endln;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else 
        s << "AlphaOSGeneralized - no associated AnalysisModel\n";
}


int AlphaOSGeneralized::formElementResidual(void)
{
    // calculate Residual Force     
    AnalysisModel *theModel = this->getAnalysisModel();
    LinearSOE *theSOE = this->getLinearSOE();
    
    // loop through the FE_Elements and add the residual
    FE_Element *elePtr;
    int res = 0;
    FE_EleIter &theEles = theModel->getFEs();
    while((elePtr = theEles()) != 0)  {
        if (theSOE->addB(elePtr->getResidual(this),elePtr->getID()) < 0)  {
            opserr << "WARNING AlphaOSGeneralized::formElementResidual -";
            opserr << " failed in addB for ID " << elePtr->getID();
            res = -2;
        }        
        if (statusFlag == CURRENT_TANGENT)  {
            if (theSOE->addB(elePtr->getK_Force(*Ut-*Upt), elePtr->getID(), alphaF-1.0) < 0)  {
                opserr << "WARNING AlphaOSGeneralized::formElementResidual -";
                opserr << " failed in addB for ID " << elePtr->getID();
                res = -2;
            }
        } else if (statusFlag == INITIAL_TANGENT)  {
            if (theSOE->addB(elePtr->getKi_Force(*Ut-*Upt), elePtr->getID(), alphaF-1.0) < 0)  {
                opserr << "WARNING AlphaOSGeneralized::formElementResidual -";
                opserr << " failed in addB for ID " << elePtr->getID();
                res = -2;
            }
        }
    }

    return res;
}
