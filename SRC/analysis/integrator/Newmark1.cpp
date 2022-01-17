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
                                                                        
// $Revision: 1.9 $
// $Date: 2007-04-02 23:42:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/Newmark1.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/Newmark1.C
// 
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the implementation of the Newmark1 class.
//
// What: "@(#) Newmark1.C, revA"

#include <Newmark1.h>
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

void *
OPS_ADD_RUNTIME_VPV(OPS_Newmark1)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata != 2 && numdata != 6) {
	opserr << "WARNING integrator Newmark1 gamma beta <alphaM> <betaKcurrent> <betaKi> <betaKlastCommitted>\n";
	return 0;
    }

    double data[6] = {0,0,0,0,0,0};
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING integrator Newmark1 invalid double inputs\n";	  
	return 0;
    }

    double gamma = data[0];
    double beta = data[1];
    double alphaM = data[2], betaK=data[3], betaKi=data[4], betaKc=data[5];

    if (numdata == 2)
	return new Newmark1(gamma,beta);       
    else
	return new Newmark1(gamma,beta,alphaM,betaK,betaKi, betaKc);
}

Newmark1::Newmark1()
:TransientIntegrator(INTEGRATOR_TAGS_Newmark1),
 gamma(0), beta(0), 
 alphaM(0.0), betaK(0.0), betaKi(0.0),
 c1(0.0), c2(0.0), c3(0.0), c4(0.0),
 Up(0), Updot(0),  U(0), Udot(0), Udotdot(0)
{
    
}

Newmark1::Newmark1(double theGamma, double theBeta, bool dispFlag)
:TransientIntegrator(INTEGRATOR_TAGS_Newmark1),
 gamma(theGamma), beta(theBeta), 
 alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
 c1(0.0), c2(0.0), c3(0.0), c4(0.0),
 Up(0), Updot(0),  U(0), Udot(0), Udotdot(0)
{

}

Newmark1::Newmark1(double theGamma, double theBeta, 
		   double alpham, double betak, double betaki , double betakc)
:TransientIntegrator(INTEGRATOR_TAGS_Newmark1),
 gamma(theGamma), beta(theBeta), 
 alphaM(alpham), betaK(betak), betaKi(betaki), betaKc(betakc),
 c1(0.0), c2(0.0), c3(0.0), c4(0.0),
 Up(0), Updot(0),  U(0), Udot(0), Udotdot(0)
{

}

Newmark1::~Newmark1()
{
  // clean up the memory created
  if (Up != 0)
    delete Up;
  if (Updot != 0)
    delete Updot;
  if (U != 0)
    delete U;
  if (Udot != 0)
    delete Udot;
  if (Udotdot != 0)
    delete Udotdot;
}


int
Newmark1::initialize(void)
{
  return 0;
}

int
Newmark1::newStep(double deltaT)
{

  if (beta == 0 || gamma == 0 ) {
    opserr << "Newton::newStep() - error in variable\n";
    opserr << "gamma = " << gamma << " beta= " << beta << endln;
    return -1;
  }

  if (deltaT <= 0.0) {
    opserr << "Newmark1::newStep() - error in variable\n";
    opserr << "dT = " << deltaT << endln;
    return -2;	
  }

  // set the constants
  c1 = 1.0;
  c2 = gamma/(beta*deltaT);
  c3 = 1.0/(beta*deltaT*deltaT);
  c4 = gamma*deltaT;


  // set the new trial response quantities
  AnalysisModel *theModel = this->getAnalysisModel();

  if (U == 0) {
    opserr << "Newton::newStep() - domainChange() failed or hasn't been called\n";
    return -3;	
  }

  // determine predicted quantities at time t + delta t
  U->addVector(1.0, *Udot, deltaT);
  double a1 = deltaT * deltaT * (.5 - beta);
  U->addVector(1.0, *Udotdot, a1);

  double a2 = deltaT * (1 - gamma);
  Udot->addVector(1.0, *Udotdot,a2);

  Udotdot->Zero();

  *Up = *U;
  *Updot = *Udot;

  theModel->setResponse(*U,*Udot,*Udotdot);        

  // increment the time and apply the load
  double time = theModel->getCurrentDomainTime();
  time +=deltaT;

  if (theModel->updateDomain(time, deltaT) < 0) {
    opserr << "Newmark1::newStep() - failed to update the domain\n";
    return -4;
  }
  
  return 0;
}

int
Newmark1::revertToLastStep() {
  return this->domainChanged();
}

int
Newmark1::formEleTangent(FE_Element *theEle)
{
  theEle->zeroTangent();
  if (statusFlag == CURRENT_TANGENT) {
    theEle->addKtToTang(c1);
    theEle->addCtoTang(c2);
    theEle->addMtoTang(c3);
  } else if (statusFlag == INITIAL_TANGENT) {
    theEle->addKiToTang(c1);
    theEle->addCtoTang(c2);
    theEle->addMtoTang(c3);
  }

  return 0;
}    


int
Newmark1::formNodTangent(DOF_Group *theDof)
{
  theDof->zeroTangent();
  theDof->addMtoTang(c3);
  theDof->addCtoTang(c2);      

  return(0);
}    



int 
Newmark1::domainChanged()
{
  AnalysisModel *myModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();
  const Vector &x = theLinSOE->getX();
  int size = x.Size();

  // if damping factors exist set them in the ele & node of the domain
  if (alphaM != 0.0 || betaK != 0.0 || betaKi != 0.0 || betaKc != 0.0)
    myModel->setRayleighDampingFactors(alphaM, betaK, betaKi, betaKc);
  
  // create the new Vector objects
  if (U == 0 || U ->Size() != size) {

    // delete the old
    if (Up != 0)
      delete Up;
    if (Updot != 0)
      delete Updot;
    if (U != 0)
      delete U;
    if (Udot != 0)
      delete Udot;
    if (Udotdot != 0)
      delete Udotdot;
    
    // create the new
    Up = new Vector(size);
    Updot = new Vector(size);
    U = new Vector(size);
    Udot = new Vector(size);
    Udotdot = new Vector(size);

    // cheack we obtained the new
    if (Up == 0 || Up->Size() != size ||
	Updot == 0 || Updot->Size() != size ||
	U == 0 || U->Size() != size ||
	Udot == 0 || Udot->Size() != size ||
	Udotdot == 0 || Udotdot->Size() != size) {
      
      opserr << "Newmark1::domainChanged - ran out of memory\n";

      // delete the old
      if (Up != 0)
	delete Up;
      if (Updot != 0)
	delete Updot;
      if (U != 0)
	delete U;
      if (Udot != 0)
	delete Udot;
      if (Udotdot != 0)
	delete Udotdot;

      Up = 0; Updot = 0; 
      U = 0; Udot = 0; Udotdot = 0;
      return -1;
    }
  }        
    
  // now go through and populate U, Udot and Udotdot by iterating through
  // the DOF_Groups and getting the last committed velocity and accel

  DOF_GrpIter &theDOFs = myModel->getDOFs();
  DOF_Group *dofPtr;
    
  while ((dofPtr = theDOFs()) != 0) {
    const ID &id = dofPtr->getID();
    int idSize = id.Size();


	int i;
    const Vector &disp = dofPtr->getCommittedDisp();	
    for (i=0; i < idSize; i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*U)(loc) = disp(i);		
      }
    }

    const Vector &vel = dofPtr->getCommittedVel();
    for (i=0; i < idSize; i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*Udot)(loc) = vel(i);
      }
    }

    const Vector &accel = dofPtr->getCommittedAccel();	
    for (i=0; i < idSize; i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*Udotdot)(loc) = accel(i);
      }
    }

    /** NOTE WE CAN't DO TOGETHER BECAUSE DOF_GROUPS USING SINGLE VECTOR ******
    for (int i=0; i < id.Size(); i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*U)(loc) = disp(i);		
 	(*Udot)(loc) = vel(i);
 	(*Udotdot)(loc) = accel(i);
      }
    }
    *******************************************************************************/

  }    
  return 0;
}


int
Newmark1::update(const Vector &deltaU)
{
  AnalysisModel *theModel = this->getAnalysisModel();
  if (theModel == 0) {
    opserr << "WARNING Newmark1::update() - no AnalysisModel set\n";
    return -1;
  }	

  // check domainChanged() has been called, i.e. Ut will not be zero
  if (U == 0) {
    opserr << "WARNING Newmark1::update() - domainChange() failed or not called\n";
    return -2;
  }	

  // check deltaU is of correct size
  if (deltaU.Size() != U->Size()) {
    opserr << "WARNING Newmark1::update() - Vectors of incompatible size ";
    opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
    return -3;
  }
    
  //  determine the response at t+delta t
  (*U) += deltaU;

  (*Udotdot) = (*U);
  (*Udotdot) -= (*Up);
  (*Udotdot) *= c3;

  (*Udot) = (*Updot);
  Udot->addVector(1.0, *Udotdot,c4);
  
  // update the responses at the DOFs
  theModel->setResponse(*U,*Udot,*Udotdot);        
  if (theModel->updateDomain() < 0) {
    opserr << "Newmark1::newStep() - failed to update the domain\n";
    return -4;
  }
		    
  return 0;
}    

const Vector &
Newmark1::getVel()
{
  return *Udot;
}

int
Newmark1::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(7);
    data(0) = gamma;
    data(1) = beta;
    data(2) = 1.0;	
    data(3) = alphaM;
    data(4) = betaK;
    data(5) = betaKi;
    data(6) = betaKc;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
	opserr << "WARNING Newmark1::sendSelf() - could not send data\n";
	return -1;
    }	
    return 0;
}

int
Newmark1::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(7);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
	opserr << "WARNING Newmark1::recvSelf() - could not receive data\n";
	gamma = 0.5; beta = 0.25; 
	return -1;
    }
    
    gamma = data(0);
    beta = data(1);
    alphaM = data(3);
    betaK = data(4);
    betaKi = data(5);
    betaKc = data(6);
      
    return 0;
    
}

void
Newmark1::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t Newmark1 - currentTime: " << currentTime;
	s << "  gamma: " << gamma << "  beta: " << beta << endln;
	s << " c1: " << c1 << " c2: " << c2 << " c3: " << c3 << endln;
	s << "  Rayleigh Damping - alphaM: " << alphaM;
	s << "  betaK: " << betaK << "  betaKi: " << betaKi << endln;	    
    } else 
	s << "\t Newmark1 - no associated AnalysisModel\n";
}




