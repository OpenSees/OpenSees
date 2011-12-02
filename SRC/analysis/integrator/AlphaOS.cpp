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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/AlphaOS.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/05
// Revision: A
//
// Description: This file contains the implementation of the AlphaOS class.
//
// What: "@(#) AlphaOS.cpp, revA"

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


AlphaOS::AlphaOS()
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS),
    alpha(1.0), beta(0.0), gamma(0.0),
    updDomFlag(0), deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0), 
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Upt(0), Uptdot(0)
{
    
}


AlphaOS::AlphaOS(double _alpha,
    bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS),
    alpha(_alpha), beta((2-_alpha)*(2-_alpha)*0.25), gamma(1.5-_alpha),
    updDomFlag(upddomflag), deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Upt(0), Uptdot(0)
{
    
}


AlphaOS::AlphaOS(double _alpha, 
    double _alphaM, double _betaK, double _betaKi, double _betaKc,
    bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS),
    alpha(_alpha), beta((2-_alpha)*(2-_alpha)*0.25), gamma(1.5-_alpha),  
    updDomFlag(upddomflag), deltaT(0.0),
    alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Upt(0), Uptdot(0)
{
    
}


AlphaOS::AlphaOS(double _alpha, double _beta, double _gamma,
    bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS),
    alpha(_alpha), beta(_beta), gamma(_gamma),
    updDomFlag(upddomflag), deltaT(0.0),
    alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Upt(0), Uptdot(0)
{
    
}


AlphaOS::AlphaOS(double _alpha, double _beta, double _gamma,
    double _alphaM, double _betaK, double _betaKi, double _betaKc,
    bool upddomflag)
    : TransientIntegrator(INTEGRATOR_TAGS_AlphaOS),
    alpha(_alpha), beta(_beta), gamma(_gamma),
    updDomFlag(upddomflag), deltaT(0.0),
    alphaM(_alphaM), betaK(_betaK), betaKi(_betaKi), betaKc(_betaKc),
    updateCount(0), c1(0.0), c2(0.0), c3(0.0),
    Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
    Ualpha(0), Ualphadot(0), Upt(0), Uptdot(0)
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
    if (Uptdot != 0)
        delete Uptdot;
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
    theModel->setResponse(*Ualpha,*Ualphadot,*Udotdot);
    
    // increment the time to t+alpha*deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += alpha*deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "AlphaOS::newStep() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
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
    
    if (statusFlag == CURRENT_TANGENT)  {
        theEle->addKtToTang(alpha*c1);
        theEle->addCtoTang(alpha*c2);
        theEle->addMtoTang(c3);
    } else if (statusFlag == INITIAL_TANGENT)  {
        theEle->addKiToTang(alpha*c1);
        theEle->addCtoTang(alpha*c2);
        theEle->addMtoTang(c3);
    }
    
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
    
    // if damping factors exist set them in the element & node of the domain
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
            Upt == 0 || Upt->Size() != size ||
            Uptdot == 0 || Uptdot->Size() != size)  {
            
            opserr << "AlphaOS::domainChanged - ran out of memory\n";
            
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
            if (Uptdot != 0)
                delete Uptdot;
            
            Ut = 0; Utdot = 0; Utdotdot = 0;
            U = 0; Udot = 0; Udotdot = 0;
            Ualpha = 0; Ualphadot = 0;
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
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING AlphaOS::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check deltaU is of correct size
    if (deltaU.Size() != U->Size())  {
        opserr << "WARNING AlphaOS::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << "\n";
        return -3;
    }
    
    // save the predictor displacements and velocities
    (*Upt) = *U;
    (*Uptdot) = *Udot;

    //  determine the response at t+deltaT
    (*U) += deltaU;

    Udot->addVector(1.0, deltaU, c2);

    Udotdot->addVector(1.0, deltaU, c3);
    
    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);
    if (updDomFlag == true)  {
        if (theModel->updateDomain() < 0)  {
            opserr << "AlphaOS::update() - failed to update the domain\n";
            return -4;
        }
    }
    
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
    
    return theModel->commitDomain();
}    


int AlphaOS::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(8);
    data(0) = alpha;
    data(1) = beta;
    data(2) = gamma;
    data(3) = alphaM;
    data(4) = betaK;
    data(5) = betaKi;
    data(6) = betaKc;
    if (updDomFlag == false) 
        data(7) = 0.0;
    else
        data(7) = 1.0;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING AlphaOS::sendSelf() - could not send data\n";
        return -1;
    }

    return 0;
}


int AlphaOS::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(8);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING AlphaOS::recvSelf() - could not receive data\n";
        return -1;
    }
    
    alpha  = data(0);
    beta   = data(1);
    gamma  = data(2);
    alphaM = data(3);
    betaK  = data(4);
    betaKi = data(5);
    betaKc = data(6);
    if (data(7) == 0.0)
        updDomFlag = false;
    else
        updDomFlag = true;
    
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
        s << "  Rayleigh Damping - alphaM: " << alphaM << "  betaK: " << betaK;
        s << "  betaKi: " << betaKi << "  betaKc: " << betaKc << endln;	    
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
    int res = 0;    
    FE_EleIter &theEles = theModel->getFEs();
    while((elePtr = theEles()) != 0)  {
        if (theSOE->addB(elePtr->getResidual(this),elePtr->getID()) < 0)  {
            opserr << "WARNING AlphaOS::formElementResidual -";
            opserr << " failed in addB for ID " << elePtr->getID();
            res = -2;
        }        
        if (statusFlag == CURRENT_TANGENT)  {
            if (theSOE->addB(elePtr->getK_Force(*Ut-*Upt), elePtr->getID(), alpha-1.0) < 0)  {
                opserr << "WARNING AlphaOS::formElementResidual -";
                opserr << " failed in addB for ID " << elePtr->getID();
                res = -2;
            }
        } else if (statusFlag == INITIAL_TANGENT)  {
            if (theSOE->addB(elePtr->getKi_Force(*Ut-*Upt), elePtr->getID(), alpha-1.0) < 0)  {
                opserr << "WARNING AlphaOS::formElementResidual -";
                opserr << " failed in addB for ID " << elePtr->getID();
                res = -2;
            }
        }
    }

    return res;
}
