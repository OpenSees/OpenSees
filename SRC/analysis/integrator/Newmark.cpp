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
#include <FE_EleIter.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>
#include <NodeIter.h>
#include <Domain.h>
#include <Node.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <elementAPI.h>
#include <fstream>
//#include<ReliabilityDomain.h>//Abbas
#include<Parameter.h>
#include<ParameterIter.h>//Abbas
static bool converged = false;
static int count = 0;

void *
OPS_Newmark(void)
{
  // Pointer to a uniaxial material that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want Newmark $gamma $beta <-form $typeUnknown>\n";
    return 0;
  }

  int dispFlag = 1;
  double dData[2];
  int numData = 2;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want Newmark $gamma $beta <-form $typeUnknown>\n";
    return 0;
  }
  
  if (argc == 2)
    theIntegrator = new Newmark(dData[0], dData[1]);
  else {
    //    char nextString[10];
    const char *nextString = OPS_GetString();
    //    OPS_GetString(nextString, 10);
    if (strcmp(nextString,"-form") == 0) {
      //      OPS_GetString(nextString, 10);
      nextString = OPS_GetString();
      if ((nextString[0] == 'D') || (nextString[0] == 'd')) 
	dispFlag = 1;
      else if ((nextString[0] == 'A') || (nextString[0] == 'a')) 
	dispFlag = 3;      
      else if ((nextString[0] == 'V') || (nextString[0] == 'v')) 
	dispFlag = 2;      
    }    
    theIntegrator = new Newmark(dData[0], dData[1], dispFlag);
  }

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating Newmark integrator\n";

  return theIntegrator;
}


Newmark::Newmark(int classTag)
    : TransientIntegrator(classTag),
      displ(true), gamma(0), beta(0), 
      c1(0.0), c2(0.0), c3(0.0), 
      Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
      determiningMass(false),
      sensitivityFlag(0), gradNumber(0), massMatrixMultiplicator(0),
      dampingMatrixMultiplicator(0), assemblyFlag(0), independentRHS(),
      dUn(), dVn(), dAn()
{
    
}


Newmark::Newmark(double _gamma, double _beta, bool dispFlag, bool aflag, int classTag_)
    : TransientIntegrator(classTag_),
      displ(dispFlag), gamma(_gamma), beta(_beta), 
      c1(0.0), c2(0.0), c3(0.0), 
      Ut(0), Utdot(0), Utdotdot(0), U(0), Udot(0), Udotdot(0),
      determiningMass(false),
      sensitivityFlag(0), gradNumber(0), massMatrixMultiplicator(0),
      dampingMatrixMultiplicator(0), assemblyFlag(aflag), independentRHS(),
      dUn(), dVn(), dAn()
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

    // clean up sensitivity
    if (massMatrixMultiplicator!=0)
	delete massMatrixMultiplicator;
    
    if (dampingMatrixMultiplicator!=0)
	delete dampingMatrixMultiplicator;
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
    if (displ == 1)  {
        c1 = 1.0;
        c2 = gamma/(beta*deltaT);
        c3 = 1.0/(beta*deltaT*deltaT);
    } else if (displ == 2) {
	c1 = deltaT*beta/gamma;
	c2 = 1.0;
	c3 = 1.0/(gamma*deltaT);
    } else if (displ == 3) {
        c1 = beta*deltaT*deltaT;
        c2 = gamma*deltaT;
        c3 = 1.0;
    }
    
    if (U == 0)  {
        opserr << "Newmark::newStep() - domainChange() failed or hasn't been called\n";
        return -3;	
    }
    
    // set response at t to be that at t+deltaT of previous step
    /*
    if (converged == true) {
      std::ofstream outfile;

      outfile.open("Newmark.dat",std::ofstream::out | std::ofstream::app);
      int size = U->Size();
      for (int i=0; i<size; i++)
	outfile << (*U)(i) << " ";
      outfile << "\n";
      outfile.close();
    }
    */

    converged = true;

    (*Ut) = *U;        
    (*Utdot) = *Udot;  
    (*Utdotdot) = *Udotdot;
    
    if (displ == 1 || displ == 2)  {    
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


const Vector &
Newmark::getVel(void)
{
  return *Udot;
}

int Newmark::revertToLastStep()
{
  // set response at t+deltaT to be that at t .. for next newStep
  converged = false;
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
    } else if (statusFlag == HALL_TANGENT)  {
        theEle->addKtToTang(c1*cFactor);
        theEle->addKiToTang(c1*iFactor);
        theEle->addCtoTang(c2);
        theEle->addMtoTang(c3);
    } else {
      opserr << "Newmark::formEleTangent - unknown FLAG\n";
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
	dUn.resize(size); dUn.Zero();
	dVn.resize(size); dVn.Zero();
	dAn.resize(size); dAn.Zero();
        
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
    if (displ == 1)  {
        (*U) += deltaU;

        Udot->addVector(1.0, deltaU, c2);
        Udotdot->addVector(1.0, deltaU, c3);
    } else if (displ == 2) {

	U->addVector(1.0, deltaU, c1);
	(*Udot) += deltaU;
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
    data(2) = displ;

    
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
    displ  = data(2);

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

int Newmark::formEleResidual(FE_Element* theEle)
{
    if(sensitivityFlag == 0) {  // no sensitivity
	this->TransientIntegrator::formEleResidual(theEle);
    } else {
	
	theEle->zeroResidual();

	// Compute the time-stepping parameters on the form
	// udotdot = a1*ui+1 + a2*ui + a3*udoti + a4*udotdoti
	// udot    = a5*ui+1 + a6*ui + a7*udoti + a8*udotdoti
	// (see p. 166 of Chopra)

	// The constants are:
	// a1 = 1.0/(beta*dt*dt)
	// a2 = -1.0/(beta*dt*dt)
	// a3 = -1.0/beta*dt
	// a4 = 1.0 - 1.0/(2.0*beta)
	// a5 = gamma/(beta*dt)
	// a6 = -gamma/(beta*dt)
	// a7 = 1.0 - gamma/beta
	// a8 = 1.0 - gamma/(2.0*beta)

	// We can make use of the data members c2 and c3 of this class. 
	// As long as disp==true, they are defined as:
	// c2 = gamma/(beta*dt)
	// c3 = 1.0/(beta*dt*dt)

	// So, the constants can be computed as follows:
	if (displ != 1) {
	    opserr << "ERROR: Newmark::formEleResidual() -- the implemented"
		   << " scheme only works if the displ variable is set to true." << endln;
	}
	double a2 = -c3;
	double a3 = -c2/gamma;
	double a4 = 1.0 - 1.0/(2.0*beta);
	double a6 = -c2;
	double a7 = 1.0 - gamma/beta;
	double dt = gamma/(beta*c2);
	double a8 = dt*(1.0 - gamma/(2.0*beta));

	// Pre-compute the vectors involving a2, a3, etc.
	//Vector tmp1 = V*a2 + Vdot*a3 + Vdotdot*a4;
	int vectorSize = U->Size();
	Vector dUn(vectorSize);
	Vector dVn(vectorSize);
	Vector dAn(vectorSize);
	int i, loc;

	AnalysisModel *myModel = this->getAnalysisModel();
	DOF_GrpIter &theDOFs = myModel->getDOFs();
	DOF_Group *dofPtr;
	while ((dofPtr = theDOFs()) != 0) {

		const ID &id = dofPtr->getID();
		int idSize = id.Size();
		const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
		for (i = 0; i < idSize; i++) {
			loc = id(i);
			if (loc >= 0) {
				dUn(loc) = dispSens(i);
			}
		}

		const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
		for (i = 0; i < idSize; i++) {
			loc = id(i);
			if (loc >= 0) {
				dVn(loc) = velSens(i);
			}
		}

		const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
		for (i = 0; i < idSize; i++) {
			loc = id(i);
			if (loc >= 0) {
				dAn(loc) = accelSens(i);
			}
		}
	}



	// Pre-compute the vectors involving a2, a3, etc.
	//Vector tmp1 = V*a2 + Vdot*a3 + Vdotdot*a4;
	Vector tmp1(vectorSize);
	tmp1.addVector(0.0, dUn, a2);
	tmp1.addVector(1.0, dVn, a3);
	tmp1.addVector(1.0, dAn, a4);
	//Vector tmp2 = V*a6 + Vdot*a7 + Vdotdot*a8;
	Vector tmp2(vectorSize);
	tmp2.addVector(0.0, dUn, a6);
	tmp2.addVector(1.0, dVn, a7);
	tmp2.addVector(1.0, dAn, a8);

	if (massMatrixMultiplicator == 0)
	    massMatrixMultiplicator = new Vector(tmp1.Size());
	if (dampingMatrixMultiplicator == 0)
	    dampingMatrixMultiplicator = new Vector(tmp2.Size());

	(*massMatrixMultiplicator) = tmp1;
	(*dampingMatrixMultiplicator) = tmp2;


	// Now we're ready to make calls to the FE Element:

	// The term -dPint/dh|u fixed
	theEle->addResistingForceSensitivity(gradNumber); 

	// The term -dM/dh*acc
	theEle->addM_ForceSensitivity(gradNumber, *Udotdot, -1.0);

	// The term -M*(a2*v + a3*vdot + a4*vdotdot)
	theEle->addM_Force(*massMatrixMultiplicator,-1.0);

	// The term -C*(a6*v + a7*vdot + a8*vdotdot)
	theEle->addD_Force(*dampingMatrixMultiplicator,-1.0);

	// The term -dC/dh*vel
	theEle->addD_ForceSensitivity(gradNumber, *Udot,-1.0);
		
    }

    return 0;
}

int
Newmark::formNodUnbalance(DOF_Group *theDof)
{

    if (sensitivityFlag == 0) {  // NO SENSITIVITY ANALYSIS

	this->TransientIntegrator::formNodUnbalance(theDof);

    }
    else {  // ASSEMBLE ALL TERMS

	theDof->zeroUnbalance();


	// The term -M*(a2*v + a3*vdot + a4*vdotdot)
	theDof->addM_Force(*massMatrixMultiplicator,-1.0);


	// The term -dM/dh*acc
	theDof->addM_ForceSensitivity(*Udotdot, -1.0);


	// The term -C*(a6*v + a7*vdot + a8*vdotdot)
	theDof->addD_Force(*dampingMatrixMultiplicator,-1.0);


	// The term -dC/dh*vel
	theDof->addD_ForceSensitivity(*Udot,-1.0);


	// In case of random loads (have already been formed by 'applyLoadSensitivity')
	theDof->addPtoUnbalance();

    }


    return 0;
}

int 
Newmark::formSensitivityRHS(int passedGradNumber)
{
    sensitivityFlag = 1;

    // Set a couple of data members
    gradNumber = passedGradNumber;

    // Get pointer to the SOE
    LinearSOE *theSOE = this->getLinearSOE();

    // Possibly set the independent part of the RHS
    if (assemblyFlag != 0) {
	theSOE->setB(independentRHS);
    }

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
Newmark::formIndependentSensitivityRHS()
{
    // For now; don't use this
/*
  sensitivityFlag = 2; // Tell subsequent methods what to be assembled

  // Get pointer to the SOE
  LinearSOE *theSOE = this->getLinearSOEPtr();


  // Get the analysis model
  AnalysisModel *theModel = this->getAnalysisModelPtr();

	
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


  // Set the data member of this class
  independentRHS = theSOE->getB();


  // Reset the sensitivity flag
  sensitivityFlag = 0;
*/

    return 0;
}

int 
Newmark::saveSensitivity(const Vector & vNew,int gradNum,int numGrads)
{

    // Compute Newmark parameters in general notation
    double a1 = c3;
    double a2 = -c3;
    double a3 = -c2/gamma;
    double a4 = 1.0 - 1.0/(2.0*beta);
    double a5 = c2;
    double a6 = -c2;
    double a7 = 1.0 - gamma/beta;
    double dt = gamma/(beta*c2);
    double a8 = dt*(1.0 - gamma/(2.0*beta));



	// Obtain sensitivity vectors from previous step modified by lei July 2018
	int vectorSize = U->Size();
	Vector dUn(vectorSize);
	Vector dVn(vectorSize);
	Vector dAn(vectorSize);
	int i, loc;

	AnalysisModel *myModel = this->getAnalysisModel();
	DOF_GrpIter &theDOFs = myModel->getDOFs();
	DOF_Group *dofPtr;
	while ((dofPtr = theDOFs()) != 0) {

		const ID &id = dofPtr->getID();
		int idSize = id.Size();
		const Vector &dispSens = dofPtr->getDispSensitivity(gradNumber);
		for (i = 0; i < idSize; i++) {
			loc = id(i);
			if (loc >= 0) {
				dUn(loc) = dispSens(i);
			}
		}

		const Vector &velSens = dofPtr->getVelSensitivity(gradNumber);
		for (i = 0; i < idSize; i++) {
			loc = id(i);
			if (loc >= 0) {
				dVn(loc) = velSens(i);
			}
		}

		const Vector &accelSens = dofPtr->getAccSensitivity(gradNumber);
		for (i = 0; i < idSize; i++) {
			loc = id(i);
			if (loc >= 0) {
				dAn(loc) = accelSens(i);
			}
		}
	}



    // Compute new acceleration and velocity vectors:
    Vector vdotNew(vectorSize);
    Vector vdotdotNew(vectorSize);
    //(*vdotdotNewPtr) = vNew*a1 + V*a2 + Vdot*a3 + Vdotdot*a4;
    vdotdotNew.addVector(0.0, vNew, a1);
    vdotdotNew.addVector(1.0, dUn, a2);
    vdotdotNew.addVector(1.0, dVn, a3);
    vdotdotNew.addVector(1.0, dAn, a4);
    
    //(*vdotNewPtr) = vNew*a5 + V*a6 + Vdot*a7 + Vdotdot*a8;
    vdotNew.addVector(0.0, vNew, a5);
    vdotNew.addVector(1.0, dUn, a6);
    vdotNew.addVector(1.0, dVn, a7);
    vdotNew.addVector(1.0, dAn, a8);

    // update
    dUn = vNew;
    dVn = vdotNew;
    dAn = vdotdotNew;

    // Now we can save vNew, vdotNew and vdotdotNew
    //AnalysisModel *myModel = this->getAnalysisModel();
    DOF_GrpIter &theDOFGrps = myModel->getDOFs();
    DOF_Group 	*dofPtr1;
    while ( (dofPtr1 = theDOFGrps() ) != 0)  {
	dofPtr1->saveSensitivity(vNew,vdotNew,vdotdotNew,gradNum,numGrads);
    }
	
    return 0;
}

int 
Newmark::commitSensitivity(int gradNum, int numGrads)
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

double
Newmark::getCFactor(void) {
  return c2;
}





int 
Newmark::computeSensitivities(void)
{
  //  opserr<<" computeSensitivity::start"<<endln; 
  LinearSOE *theSOE = this->getLinearSOE();
  
  /*
    if (theAlgorithm == 0) {
    opserr << "ERROR the FE algorithm must be defined before ";
    opserr << "the sensitivity algorithm\n";
    return -1;
    }
  */
  /*
  // Get pointer to the system of equations (SOE)
  LinearSOE *theSOE = theAlgorithm->getLinearSOEptr();
  if (theSOE == 0) {
  opserr << "ERROR the FE linearSOE must be defined before ";
  opserr << "the sensitivity algorithm\n";
  return -1;
  }
  
  // Get pointer to incremental integrator
  IncrementalIntegrator *theIncInt = theAlgorithm->getIncrementalIntegratorPtr();
  //	IncrementalIntegrator *theIncIntSens=theAlgorithm->getIncrementalIntegratorPtr();//Abbas
  if (theIncInt == 0 ) {
  opserr << "ERROR the FE integrator must be defined before ";
  opserr << "the sensitivity algorithm\n";
  return -1;
  }
  
  // Form current tangent at converged state
  // (would be nice with an if-statement here in case
  // the current tangent is already formed)
  if (this->formTangent(CURRENT_TANGENT) < 0){
  opserr << "WARNING SensitivityAlgorithm::computeGradients() -";
  opserr << "the Integrator failed in formTangent()\n";
  return -1;
  }
  */
  // Zero out the old right-hand side of the SOE
  theSOE->zeroB();
  
  if (this == 0) {
    opserr << "ERROR SensitivityAlgorithm::computeSensitivities() -";
    opserr << "the SensitivityIntegrator is NULL\n";
    return -1;
  }
  
  // Form the part of the RHS which are indepent of parameter
  this->formIndependentSensitivityRHS();
  AnalysisModel *theModel = this->getAnalysisModel();  //Abbas 
  Domain *theDomain=theModel->getDomainPtr();//Abbas
  ParameterIter &paramIter = theDomain->getParameters();
  //	opserr<<" get parameters "<<theDomain->getParameters()<<endln;//Abbas.......
  Parameter *theParam;
  // De-activate all parameters
  while ((theParam = paramIter()) != 0)
    theParam->activate(false);
  
  // Now, compute sensitivity wrt each parameter
  int numGrads = theDomain->getNumParameters();
  //opserr<<"the numGrads is "<<numGrads<<endln;//Abbas...............................
  paramIter = theDomain->getParameters();
  
  while ((theParam = paramIter()) != 0) {
    
    // Activate this parameter
    theParam->activate(true);
    
    // Zero the RHS vector
    theSOE->zeroB();
    
    // Get the grad index for this parameter
    int gradIndex = theParam->getGradIndex();
    //   opserr<<"gradNumber = "<<gradIndex<<endln;
    // Form the RHS
    this->formSensitivityRHS(gradIndex);
    
    // Solve for displacement sensitivity
    
    theSOE->solve();
    // Save sensitivity to nodes
    this->saveSensitivity( theSOE->getX(), gradIndex, numGrads );
    
    
    
    // Commit unconditional history variables (also for elastic problems; strain sens may be needed anyway)
    this->commitSensitivity(gradIndex, numGrads);
    
    // De-activate this parameter for next sensitivity calc
    theParam->activate(false);
    //  opserr<<"LoadControl::..........ComputeSensitivities. end"<<endln;
  }
  
  return 0;
}

