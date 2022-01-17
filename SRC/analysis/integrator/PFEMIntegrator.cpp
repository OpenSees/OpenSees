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

// $Revision: 1.0 $
// $Date: 2013-01-16 12:51:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/PFEMIntegrator.h,v $


// Written : Minjie Zhu
// Created : Jan 2013
//
// Description: This file contains the class definition for PFEMIntegrator.
// PFEMIntegrator is an algorithmic class for performing a PFEM analysis
// using the BackEuler integration scheme.
//
// What: "@(#) PFEMIntegrator.h, revA"

#include <PFEMIntegrator.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <Domain.h>
#include <Pressure_Constraint.h>
#include <Pressure_ConstraintIter.h>
#include <NodeIter.h>
#include <Node.h>
#include <LoadPatternIter.h>
#include <LoadPattern.h>
#include <FE_EleIter.h>
#include <elementAPI.h>
#include "sparseGEN/PFEMLinSOE.h"

void *
OPS_ADD_RUNTIME_VPV(OPS_PFEMIntegrator)
{
    TransientIntegrator *theIntegrator = 0;

    int dispFlag = 2;
    int init = 2;
    double dData[2] = {-1.0, -1.0};
    int numData = 2;
    if (OPS_GetNumRemainingInputArgs() > 1) {
        if (OPS_GetDouble(&numData, dData) < 0) {
            OPS_ResetCurrentInputArg(-2);
        }
    }

    if (OPS_GetNumRemainingInputArgs() > 1) {
        const char *nextString = OPS_GetString();
        if (strcmp(nextString,"-form") == 0) {
            nextString = OPS_GetString();
            if ((nextString[0] == 'D') || (nextString[0] == 'd')) {
                dispFlag = 1;
                init = 1;
            } else if ((nextString[0] == 'A') || (nextString[0] == 'a')) {
                dispFlag = 3;
                init = 3;
            } else if ((nextString[0] == 'V') || (nextString[0] == 'v')) {
                dispFlag = 2;
                init = 2;
            }
        } else {
            opserr << "WARNING: first option must be -form\n";
            return 0;
        }
        if (OPS_GetNumRemainingInputArgs() > 1) {
            nextString = OPS_GetString();
            if (strcmp(nextString, "-init") == 0) {
                nextString = OPS_GetString();
                if ((nextString[0] == 'D') || (nextString[0] == 'd')) {
                    init = 1;
                } else if ((nextString[0] == 'A') || (nextString[0] == 'a')) {
                    init = 3;
                } else if ((nextString[0] == 'V') || (nextString[0] == 'v')) {
                    init = 2;
                }
            }
        } else {
            opserr << "WARNING: second option must be -init\n";
            return 0;
        }
    }

    theIntegrator = new PFEMIntegrator(dData[0], dData[1], dispFlag, init);

    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating Newmark integrator\n";

    return theIntegrator;
}

PFEMIntegrator::PFEMIntegrator()
        : TransientIntegrator(INTEGRATOR_TAGS_PFEMIntegrator),
          displ(2), init(1), gamma(-1), beta(-1),
          c1(0.0), c2(0.0), c3(0.0),
          Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
          determiningMass(false),
          sensitivityFlag(0), gradNumber(0), massMatrixMultiplicator(0),
          dampingMatrixMultiplicator(0), assemblyFlag(false), independentRHS(),
          dUn(), dVn(), dAn()
{

}

PFEMIntegrator::PFEMIntegrator(double _gamma, double _beta, int _disp, int _init)
    : TransientIntegrator(INTEGRATOR_TAGS_PFEMIntegrator),
      displ(_disp), init(_init), gamma(_gamma), beta(_beta),
      c1(0.0), c2(0.0), c3(0.0),
      Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
      determiningMass(false),
      sensitivityFlag(0), gradNumber(0), massMatrixMultiplicator(0),
      dampingMatrixMultiplicator(0), assemblyFlag(false), independentRHS(),
      dUn(), dVn(), dAn()
{
    
}

PFEMIntegrator::~PFEMIntegrator()
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

    // clean up sensitivity
    if (massMatrixMultiplicator!=0)
        delete massMatrixMultiplicator;

    if (dampingMatrixMultiplicator!=0)
        delete dampingMatrixMultiplicator;
}


int PFEMIntegrator::newStep(double deltaT)
{
    if (beta == 0 || gamma == 0)  {
        opserr << "Newmark::newStep() - error in variable\n";
        opserr << "gamma = " << gamma << " beta = " << beta << endln;
        return -1;
    }

    if (deltaT <= 0.0)  {
        opserr << "PFEMIntegrator::newStep() - error in variable\n";
        opserr << "dT = " << deltaT << endln;
        return -2;	
    }

    // get a pointer to the AnalysisModel and Domain
    AnalysisModel *theModel = this->getAnalysisModel();
    if(theModel == 0) {
        opserr << "Analysis model has not been linked - PFEMIntegrator::newStep()\n";
        return -1;
    }
    Domain* theDomain = theModel->getDomainPtr();
    if(theDomain == 0) {
        opserr<<"WARNING: no domain is set for the model";
        opserr<<" -- PFEMIntegrator::newStep()\n";
        return -1;
    }

    // set the constants
    if (displ == 1)  {
        if (gamma>0 && beta>0) {
            c1 = 1.0;
            c2 = gamma / (beta * deltaT);
            c3 = 1.0 / (beta * deltaT * deltaT);
        } else {
            c1 = 1.0;
            c2 = 1.0 / deltaT;
            c3 = 1.0 / (deltaT * deltaT);
        }
    } else if (displ == 2) {
        if (gamma>0 && beta>0) {
            c1 = deltaT * beta / gamma;
            c2 = 1.0;
            c3 = 1.0 / (gamma * deltaT);
        } else {
            c1 = deltaT;
            c2 = 1.0;
            c3 = 1.0/deltaT;
        }
    } else if (displ == 3) {
        if (gamma>0 && beta>0) {
            c1 = beta * deltaT * deltaT;
            c2 = gamma * deltaT;
            c3 = 1.0;
        } else {
            c1 = deltaT * deltaT;
            c2 = deltaT;
            c3 = 1.0;
        }
    }

    // check if domainchange() is called
    if (U == 0)  {
        opserr << "PFEMIntegrator::newStep() - domainChange() failed or hasn't been called\n";
        return -3;	
    }

    populateUn();
    populateU();

    if (init == 1) {

        // determine new velocities and accelerations at t+deltaT
        *Udot = *Utdot;
        *Udotdot = *Utdotdot;
        if (gamma>0 && beta>0) {
            Udot->addVector(1.0-gamma/beta, *Utdotdot, deltaT*(1.0-0.5*gamma/beta));
            Udot->addVector(1.0, *U, gamma / (deltaT * beta));
            Udot->addVector(1.0, *Ut, -gamma / (deltaT * beta));
            Udotdot->addVector(1.0-0.5/beta, *Utdot, -1.0 / (beta * deltaT));
            Udotdot->addVector(1.0, *U, 1.0/(deltaT * deltaT * beta));
            Udotdot->addVector(1.0, *Ut, -1.0/(deltaT * deltaT * beta));
        } else {
            Udot->addVector(0.0, *U, 1.0/deltaT);
            Udot->addVector(1.0, *Ut, -1.0/deltaT);
            Udotdot->addVector(0.0, *U, 1.0/(deltaT*deltaT));
            Udotdot->addVector(1.0, *Ut, -1.0/(deltaT*deltaT));
            Udotdot->addVector(1.0, *Utdot, -1.0/deltaT);
        }

        // set the trial response quantities
        theModel->setVel(*Udot);
        theModel->setAccel(*Udotdot);

    } else if (init == 2) {
        // determine new displacements and accelerations at t+deltaT
        *U = *Ut;
        *Udotdot = *Utdotdot;
        if (gamma>0 && beta>0) {
            U->addVector(1.0, *Utdot, deltaT*(1.0-beta/gamma));
            U->addVector(1.0, *Udot, deltaT*beta/gamma);
            U->addVector(1.0, *Utdotdot, deltaT*deltaT*(0.5-beta/gamma));
            Udotdot->addVector((1.0-1.0/gamma), *Udot, 1.0/(gamma*deltaT));
            Udotdot->addVector(1.0, *Utdot, -1.0/(gamma*deltaT));
        } else {
            U->addVector(1.0, *Udot, deltaT);
            Udotdot->addVector(0.0, *Udot, 1.0/deltaT);
            Udotdot->addVector(1.0, *Utdot, -1.0/deltaT);
        }

        // set the trial response quantities
        theModel->setDisp(*U);
        theModel->setAccel(*Udotdot);

    } else  {
        // determine new displacements and velocities at t+deltaT
        *U = *Ut;
        *Udot = *Utdot;

        if (gamma>0 && beta>0) {
            U->addVector(1.0, *Utdot, deltaT);
            U->addVector(1.0, *Utdotdot, deltaT * deltaT * (0.5 - beta));
            U->addVector(1.0, *Udotdot, deltaT * deltaT * beta);
            Udot->addVector(1.0, *Utdotdot, deltaT * (1 - gamma));
            Udot->addVector(1.0, *Udotdot, deltaT * gamma);
        } else {
            U->addVector(1.0, *Utdot, deltaT);
            U->addVector(1.0, *Udotdot, deltaT*deltaT);
            Udot->addVector(1.0, *Udotdot, deltaT);
        }

        // set the trial response quantities
        theModel->setDisp(*U);
        theModel->setVel(*Udot);
    }
    
    // increment the time to t+deltaT and apply the load
    double time = theModel->getCurrentDomainTime();
    time += deltaT;
    if (theModel->updateDomain(time, deltaT) < 0)  {
        opserr << "PFEMIntegrator::newStep() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}

int PFEMIntegrator::commit()
{
    // get a pointer to the AnalysisModel and Domain
    AnalysisModel *theModel = this->getAnalysisModel();
    if(theModel == 0) {
        opserr << "Analysis model has not been linked - PFEMIntegrator::newStep()\n";
        return -1;
    }

    return theModel->commitDomain();
}


int PFEMIntegrator::revertToLastStep()
{
    // set response at t+deltaT to be that at t .. for next newStep
    if (U != 0)  {
        (*U) = *Ut;        
        (*Udot) = *Utdot;  
        (*Udotdot) = *Utdotdot;  
    }

    return 0;
}

int
PFEMIntegrator::formTangent(int statFlag)
{
    int result = 0;
    statusFlag = statFlag;

    LinearSOE *theLinSOE = this->getLinearSOE();
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theLinSOE == 0 || theModel == 0) {
        opserr << "WARNING TransientIntegrator::formTangent() ";
        opserr << "no LinearSOE or AnalysisModel has been set\n";
        return -1;
    }

    // the loops to form and add the tangents are broken into two for
    // efficiency when performing parallel computations

    theLinSOE->zeroA();

    // do modal damping
    bool inclModalMatrix=theModel->inclModalDampingMatrix();
    if (inclModalMatrix == true) {
        const Vector *modalValues = theModel->getModalDampingFactors();
        if (modalValues != 0) {
            this->addModalDampingMatrix(modalValues);
        }
    }


    // loop through the DOF_Groups and add the unbalance
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    DOF_Group *dofPtr;

    while ((dofPtr = theDOFs()) != 0) {
        if (theLinSOE->addA(dofPtr->getTangent(this),dofPtr->getID()) <0) {
            opserr << "TransientIntegrator::formTangent() - failed to addA:dof\n";
            result = -1;
        }
    }

    // loop through the FE_Elements getting them to add the tangent
    FE_EleIter &theEles2 = theModel->getFEs();
    FE_Element *elePtr;
    while((elePtr = theEles2()) != 0)     {
        if (theLinSOE->addA(elePtr->getTangent(this),elePtr->getID()) < 0) {
            opserr << "TransientIntegrator::formTangent() - failed to addA:ele\n";
            result = -2;
        }
    }
    return result;
}


int PFEMIntegrator::formEleTangent(FE_Element *theEle)
{
    if (determiningMass == true)
        return 0;
    
    theEle->zeroTangent();
    
    if (statusFlag == CURRENT_TANGENT)  {
        theEle->addKtToTang(c1);
        theEle->addCtoTang(c2);
        theEle->addMtoTang(c3);
        if (sensitivityFlag == 1) {
            theEle->addKgToTang(c1);
        }
    } else if (statusFlag == INITIAL_TANGENT)  {
        theEle->addKiToTang(c1);
        theEle->addCtoTang(c2);
        theEle->addMtoTang(c3);
    } else {
        opserr << "Newmark::formEleTangent - unknown FLAG\n";
    }
    
    return 0;
}    


int PFEMIntegrator::formNodTangent(DOF_Group *theDof)
{
    if (determiningMass == true)
        return 0;
    
    theDof->zeroTangent();

    theDof->addCtoTang(c2);      
    theDof->addMtoTang(c3);
    
    return 0;
}    

int
PFEMIntegrator::formEleResidual(FE_Element *theEle)
{
    if (sensitivityFlag == 0) {  // NO SENSITIVITY ANALYSIS

        this->TransientIntegrator::formEleResidual(theEle);

    }
    else {  // (ASSEMBLE ALL TERMS)

        theEle->zeroResidual();

        // Compute the time-stepping parameters on the form
        // udotdot = 1/dt*vn+1 - 1/dt*vn
        // u       = un + dt*vn+1


        // Obtain sensitivity vectors from previous step
        // dVn.resize(U->Size()); dVn.Zero();
        // Vector dUn(U->Size());

        // AnalysisModel *myModel = this->getAnalysisModel();
        // DOF_GrpIter &theDOFs = myModel->getDOFs();
        // DOF_Group *dofPtr = 0;
        // while ((dofPtr = theDOFs()) != 0) {

        //     const ID &id = dofPtr->getID();
        //     int idSize = id.Size();

        //     //const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
        //     for (int i=0; i < idSize; i++) {
        //         int loc = id(i);
        //         if (loc >= 0) {
        //             //dUn(loc) = dispSens(i);
	// 	    dUn(loc) = 0;
        //         }
        //     }

        //     //const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
        //     for (int i=0; i < idSize; i++) {
        //         int loc = id(i);
        //         if (loc >= 0) {
        //             //dVn(loc) = velSens(i);
	// 	    dVn(loc) = 0;
        //         }
        //     }
        // }

        // Now we're ready to make calls to the FE Element:

        // The term -dPint/dh|u fixed
        theEle->addResistingForceSensitivity(gradNumber); 

        // The term -dM/dh*acc
        theEle->addM_ForceSensitivity(gradNumber, *Udotdot, -1.0);

        // The term -M*(-1/dt*dvn)
        theEle->addM_Force(dVn, c3);

        // The term -K*(dun)
        theEle->addK_Force(dUn, -1.0);

	// The term -Kg*(dun)
        theEle->addKg_Force(dUn, -1.0);

        // The term -dC/dh*vel
        theEle->addD_ForceSensitivity(gradNumber, *Udot,-1.0);
		
    }

    return 0;
}    

int
PFEMIntegrator::formNodUnbalance(DOF_Group *theDof)
{
    if (sensitivityFlag == 0) {  // NO SENSITIVITY ANALYSIS

        this->TransientIntegrator::formNodUnbalance(theDof);

    } else {  // ASSEMBLE ALL TERMS

        theDof->zeroUnbalance();

        // The term -M*(-1/dt*dvn)
        theDof->addM_Force(dVn, c3);

        // The term -dM/dh*acc
        theDof->addM_ForceSensitivity(*Udotdot, -1.0);

        // The term -C*(dvn)
        //theDof->addD_Force(dVn ,-1.0);

	
        // The term -dC/dh*vel
        theDof->addD_ForceSensitivity(*Udot,-1.0);


        // In case of random loads (have already been formed by 'applyLoadSensitivity')
        theDof->addPtoUnbalance();

    }

    return 0;
}    


int PFEMIntegrator::domainChanged()
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
	dVn.resize(size); dVn.Zero();
	dUn.resize(size); dUn.Zero();
	dAn.resize(size); dAn.Zero();
        
        // check we obtained the new
        if (Ut == 0 || Ut->Size() != size ||
            Utdot == 0 || Utdot->Size() != size ||
            Utdotdot == 0 || Utdotdot->Size() != size ||
            U == 0 || U->Size() != size ||
            Udot == 0 || Udot->Size() != size ||
            Udotdot == 0 || Udotdot->Size() != size)  {
            
            opserr << "PFEMIntegrator::domainChanged - ran out of memory\n";
            
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

        const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
        for (i=0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                dUn(loc) = dispSens(i);
            }
        }

        const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
        for (i=0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                dVn(loc) = velSens(i);
            }
        }

        const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
        for (i=0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                dAn(loc) = accelSens(i);
            }
        }
    }

    return 0;
}


int PFEMIntegrator::update(const Vector &delta)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel == 0)  {
        opserr << "WARNING PFEMIntegrator::update() - no AnalysisModel set\n";
        return -1;
    }	
    
    // check domainChanged() has been called, i.e. Ut will not be zero
    if (Ut == 0)  {
        opserr << "WARNING PFEMIntegrator::update() - domainChange() failed or not called\n";
        return -2;
    }	
    
    // check deltaU is of correct size
    if (delta.Size() != U->Size())  {
        opserr << "WARNING PFEMIntegrator::update() - Vectors of incompatible size ";
        opserr << " expecting " << U->Size() << " obtained " << delta.Size() << endln;
        return -3;
    }
    
    //  determine the response at t+deltaT
    if (displ == 1)  {
        (*U) += delta;
        Udot->addVector(1.0, delta, c2);
        Udotdot->addVector(1.0, delta, c3);
    } else if (displ == 2) {
        U->addVector(1.0, delta, c1);
        (*Udot) += delta;
        Udotdot->addVector(1.0, delta, c3);
    } else  {
        U->addVector(1.0, delta, c1);
        Udot->addVector(1.0, delta, c2);
        (*Udotdot) += delta;
    }
    
    // update the response at the DOFs
    theModel->setResponse(*U,*Udot,*Udotdot);
    if (theModel->updateDomain() < 0)  {
        opserr << "PFEMIntegrator::update() - failed to update the domain\n";
        return -4;
    }
    
    return 0;
}    


const Vector &
PFEMIntegrator::getVel()
{
  return *Udot;
}

int PFEMIntegrator::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = gamma;
    data(1) = beta;
    data(2) = displ;
    data(3) = init;


    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING Newmark::sendSelf() - could not send data\n";
        return -1;
    }
    return 0;
}


int PFEMIntegrator::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0)  {
        opserr << "WARNING Newmark::recvSelf() - could not receive data\n";
        gamma = 0.5; beta = 0.25; displ = 2; init = 1;
        return -1;
    }

    gamma  = data(0);
    beta   = data(1);
    displ  = data(2);
    init   = data(3);
    return 0;
}


void PFEMIntegrator::Print(OPS_Stream &s, int flag)
{
    if (flag != 0) return;
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
        double currentTime = theModel->getCurrentDomainTime();
        s << "\t PFEMIntegrator - currentTime: " << currentTime;
        s << "  c1: " << c1 << "  c2: " << c2 << "  c3: " << c3 << endln;
    } else 
        s << "\t PFEMIntegrator - no associated AnalysisModel\n";
}

int
PFEMIntegrator::populateUn()
{
    AnalysisModel *myModel = this->getAnalysisModel();

    // now go through and populate U, Udot and Udotdot by iterating through
    // the DOF_Groups and getting the last committed velocity and accel
    DOF_GrpIter &theDOFs = myModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();

        int i;
        const Vector &disp = dofPtr->getCommittedDisp();
        for (i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                (*Ut)(loc) = disp(i);
            }
        }

        const Vector &vel = dofPtr->getCommittedVel();
        for (i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                (*Utdot)(loc) = vel(i);
            }
        }

        const Vector &accel = dofPtr->getCommittedAccel();
        for (i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                (*Utdotdot)(loc) = accel(i);
            }
        }
    }

    return 0;
}

int
PFEMIntegrator::populateU()
{
    AnalysisModel *myModel = this->getAnalysisModel();
    Domain* domain = myModel->getDomainPtr();
    if (domain == 0) return -1;

    // now go through and populate U, Udot and Udotdot by iterating through
    // the DOF_Groups and getting the last committed velocity and accel
    DOF_GrpIter &theDOFs = myModel->getDOFs();
    DOF_Group *dofPtr;
    while ((dofPtr = theDOFs()) != 0) {
        const ID &id = dofPtr->getID();
        int idSize = id.Size();
        int nodetag = dofPtr->getNodeTag();

        int i;
        const Vector &disp = dofPtr->getTrialDisp();
        for (i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                (*U)(loc) = disp(i);
            }
        }

        const Vector &vel = dofPtr->getTrialVel();
        for (i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                (*Udot)(loc) = vel(i);
            }
        }

        const Vector &accel = dofPtr->getTrialAccel();
        for (i = 0; i < idSize; i++) {
            int loc = id(i);
            if (loc >= 0) {
                (*Udotdot)(loc) = accel(i);
            }
        }
    }

    return 0;
}

// AddingSensitivity:BEGIN //////////////////////////////
int PFEMIntegrator::revertToStart()
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

int 
PFEMIntegrator::formSensitivityRHS(int passedGradNumber)
{
    sensitivityFlag = 1;


    // Set a couple of data members
    gradNumber = passedGradNumber;

    // Get pointer to the SOE
    LinearSOE *theSOE = this->getLinearSOE();


    // Get the analysis model
    AnalysisModel *theModel = this->getAnalysisModel();



    // Randomness in external load (including randomness in time series)
    // Get domain
    Domain *theDomain = theModel->getDomainPtr();

    // Loop through nodes to zero the unbalaced load
    Node *nodePtr;
    NodeIter &theNodeIter = theDomain->getNodes();
    while ((nodePtr = theNodeIter()) != 0)
	nodePtr->zeroUnbalancedLoad();


    // Loop through load patterns to add external load sensitivity
    LoadPattern *loadPatternPtr;
    LoadPatternIter &thePatterns = theDomain->getLoadPatterns();
    double time;
    while((loadPatternPtr = thePatterns()) != 0) {
        time = theDomain->getCurrentTime();
        loadPatternPtr->applyLoadSensitivity(time);
    }


    // Randomness in element/material contributions
    // Loop through FE elements
    FE_Element *elePtr;
    FE_EleIter &theEles = theModel->getFEs();    
    while((elePtr = theEles()) != 0) {
        theSOE->addB(  elePtr->getResidual(this),  elePtr->getID()  );
    }


    // Loop through DOF groups (IT IS IMPORTANT THAT THIS IS DONE LAST!)
    DOF_Group *dofPtr;
    DOF_GrpIter &theDOFs = theModel->getDOFs();
    while((dofPtr = theDOFs()) != 0) {
        theSOE->addB(  dofPtr->getUnbalance(this),  dofPtr->getID()  );
    }


    // Reset the sensitivity flag
    sensitivityFlag = 0;

    return 0;
}

int 
PFEMIntegrator::formIndependentSensitivityRHS()
{
    return 0;
}

int
PFEMIntegrator::formIndependentSensitivityLHS(int statFlag)
{
    // set the sensitivity flag
    sensitivityFlag = 1;

    // call formTangent
    int res = this->formTangent(statFlag);

    // Reset the sensitivity flag
    sensitivityFlag = 0;
    
    return res;
}

int 
PFEMIntegrator::saveSensitivity(const Vector & dVNew,int gradNum,int numGrads)
{
    // Recover sensitivity results from previous step
    AnalysisModel *myModel = this->getAnalysisModel();

    // Compute new acceleration and velocity vectors:
    int vectorSize = dVNew.Size();
    Vector dANew(vectorSize);

    // dudotdot = 1/dt*dv{n+1} - 1/dt*dvn
    dANew.addVector(0.0, dVNew, c3);
    dANew.addVector(1.0, dVn, -c3);

    // du       = dun + dt*dv{n+1}
    dUn.addVector(1.0, dVNew, c1);

    // dv
    dVn = dVNew;

    // Now we can save vNew, vdotNew and vdotdotNew
    DOF_GrpIter &theDOFGrps = myModel->getDOFs();
    DOF_Group 	*dofPtr1;
    while ( (dofPtr1 = theDOFGrps() ) != 0)  {
        dofPtr1->saveSensitivity(dUn,dVn,dANew,gradNum,numGrads);
    }
	
    return 0;
}

int 
PFEMIntegrator::commitSensitivity(int gradNum, int numGrads)
{

    // Loop through the FE_Elements and set unconditional sensitivities
    AnalysisModel *theAnalysisModel = this->getAnalysisModel();
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != 0) {
        elePtr->commitSensitivity(gradNum, numGrads);
    }

    return 0;
}
// AddingSensitivity:END ////////////////////////////////

