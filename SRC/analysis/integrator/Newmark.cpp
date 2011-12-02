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
                                                                        
// $Revision: 1.11 $
// $Date: 2003-03-06 20:32:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/Newmark.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/Newmark.C
// 
// Written: fmk 
// Created: 11/98
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

Newmark::Newmark()
:TransientIntegrator(INTEGRATOR_TAGS_Newmark),
 displ(true), gamma(0), beta(0), 
 alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 determiningMass(false)
{
    
}

Newmark::Newmark(double theGamma, double theBeta, bool dispFlag)
:TransientIntegrator(INTEGRATOR_TAGS_Newmark),
 displ(dispFlag),
 gamma(theGamma), beta(theBeta), 
 alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 determiningMass(false)
{

}

Newmark::Newmark(double theGamma, double theBeta, 
		 double alpham, double betak, 
		 double betaki, double betakc,
		 bool dispFlag)
:TransientIntegrator(INTEGRATOR_TAGS_Newmark),
 displ(dispFlag),
 gamma(theGamma), beta(theBeta), 
 alphaM(alpham), betaK(betak), betaKi(betaki), betaKc(betakc),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
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


int
Newmark::initialize(void)
{
  /***********************************************
  U->Zero();
  Udot->Zero();
  Udotdot->Zero();

  // determine the unbalance at initial time step
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  double time = theModel->getCurrentDomainTime();
  theModel->applyLoadDomain(time);
  this->formUnbalance();

  // set up the system wid mass matrix M
  LinearSOE *theLinSOE = this->getLinearSOEPtr();

  theLinSOE->zeroA();
  determiningMass = true;

  // loop through the DOF_Groups and add the unbalance
  DOF_GrpIter &theDOFs = theModel->getDOFs();
  DOF_Group *dofPtr;
    
  while ((dofPtr = theDOFs()) != 0) {
    dofPtr->zeroTangent();
    dofPtr->addMtoTang(1.0);
    if (theLinSOE->addA(dofPtr->getTangent(this),dofPtr->getID()) <0) {
      opserr << "TransientIntegrator::formTangent() - failed to addA:dof\n";
    }
  }    
    
  // loop through the FE_Elements getting them to add M to the tangent    
  FE_EleIter &theEles2 = theModel->getFEs();    
  FE_Element *elePtr;    
  while((elePtr = theEles2()) != 0)     {
    elePtr->zeroTangent();
    elePtr->addMtoTang(1.0);
    if (theLinSOE->addA(elePtr->getTangent(this),elePtr->getID()) < 0) {
      opserr << "TransientIntegrator::formTangent() - failed to addA:ele\n";
    }
  }

  determiningMass = false;

  // solve for the accelerations at initial time step
  theLinSOE->solve();
  (*Udotdot) = theLinSOE->getX();

  // set the new trial response quantities
  theModel->setResponse(*U,*Udot,*Udotdot);        
  theModel->updateDomain();
  ***************************************************/
  return 0;
}



int
Newmark::newStep(double deltaT)
{
  if (beta == 0 || gamma == 0 ) {
    opserr << "Newton::newStep() - error in variable\n";
    opserr << "gamma = " << gamma << " beta= " << beta << endln;
    return -1;
  }


  // get a pointer to the AnalysisModel
  AnalysisModel *theModel = this->getAnalysisModelPtr();

  if (displ == true) {
    if (deltaT <= 0.0) {
      opserr << "Newton::newStep() - error in variable\n";
      opserr << "dT = " << deltaT << endln;
      return -2;	
    }
    c1 = 1.0;
    c2 = gamma/(beta*deltaT);
    c3 = 1.0/(beta*deltaT*deltaT);
  } else {
    c1 = beta*deltaT*deltaT;
    c2 = gamma*deltaT;
    c3 = 1.0;
  }
    
  if (U == 0) {
    opserr << "Newton::newStep() - domainChange() failed or hasn't been called\n";
    return -3;	
  }

  // set response at t to be that at t+delta t of previous step
  (*Ut) = *U;        
  (*Utdot) = *Udot;  
  (*Utdotdot) = *Udotdot;  

  if (displ == true) {

    // set new velocity and accelerations at t + delta t
    double a1 = (1.0 - gamma/beta); 
    double a2 = (deltaT)*(1.0 - 0.5*gamma/beta);
    // (*Udot) *= a1;
    Udot->addVector(a1, *Utdotdot,a2);

    double a3 = -1.0/(beta*deltaT);
    double a4 = 1 - 0.5/beta;
    // (*Udotdot) *= a4;  
    Udotdot->addVector(a4, *Utdot,a3);
    
  } else {
    // set new displacement and velocity at t + delta t      
    double a1 = (deltaT*deltaT/2.0); 
    U->addVector(1.0, *Utdot,deltaT);
    U->addVector(1.0, *Utdotdot,a1);
    
    Udot->addVector(1.0, *Utdotdot,deltaT);
  }

  // set the new trial response quantities
  theModel->setResponse(*U,*Udot,*Udotdot);        

  // increment the time and apply the load
  double time = theModel->getCurrentDomainTime();
  time +=deltaT;
  if (theModel->updateDomain(time, deltaT) < 0) {
    opserr << "Newmark::newStep() - failed to update the domain\n";
    return -4;
  }
  
  return 0;
}



int
Newmark::revertToLastStep()
{
  // set response at t+delta t to be that at t .. for next newStep
  if (U != 0) {
    (*U) = *Ut;        
    (*Udot) = *Utdot;  
    (*Udotdot) = *Utdotdot;  
  }
  return 0;
}


int
Newmark::formEleTangent(FE_Element *theEle)
{
  if (determiningMass == true)
    return 0;

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
Newmark::formNodTangent(DOF_Group *theDof)
{
  if (determiningMass == true)
    return 0;

  theDof->zeroTangent();
  theDof->addMtoTang(c3);
  theDof->addCtoTang(c2);      

  return(0);
}    



int 
Newmark::domainChanged()
{
  AnalysisModel *myModel = this->getAnalysisModelPtr();
  LinearSOE *theLinSOE = this->getLinearSOEPtr();
  const Vector &x = theLinSOE->getX();
  int size = x.Size();

  // if damping factors exist set them in the ele & node of the domain
  if (alphaM != 0.0 || betaK != 0.0 || betaKi != 0.0 || betaKc != 0.0)
    myModel->setRayleighDampingFactors(alphaM, betaK, betaKi, betaKc);
  
  // create the new Vector objects
  if (Ut == 0 || Ut->Size() != size) {

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

    // cheack we obtained the new
    if (Ut == 0 || Ut->Size() != size ||
	Utdot == 0 || Utdot->Size() != size ||
	Utdotdot == 0 || Utdotdot->Size() != size ||
	U == 0 || U->Size() != size ||
	Udot == 0 || Udot->Size() != size ||
	Udotdot == 0 || Udotdot->Size() != size) {
      
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
Newmark::update(const Vector &deltaU)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    opserr << "WARNING Newmark::update() - no AnalysisModel set\n";
    return -1;
  }	

  // check domainChanged() has been called, i.e. Ut will not be zero
  if (Ut == 0) {
    opserr << "WARNING Newmark::update() - domainChange() failed or not called\n";
    return -2;
}	

  // check deltaU is of correct size
  if (deltaU.Size() != U->Size()) {
    opserr << "WARNING Newmark::update() - Vectors of incompatable size ";
    opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
    return -3;
  }
    
  //  determine the response at t+delta t
  if (displ == true) {
    (*U) += deltaU;
    Udot->addVector(1.0, deltaU,c2);
    Udotdot->addVector(1.0, deltaU,c3);
  } else {
    (*Udotdot) += deltaU;
    U->addVector(1.0, deltaU,c1);
    Udot->addVector(1.0, deltaU,c2);
  }

  // update the responses at the DOFs
  theModel->setResponse(*U,*Udot,*Udotdot);        
  if (theModel->updateDomain() < 0) {
    opserr << "Newmark::update() - failed to update the domain\n";
    return -4;
  }

  return 0;
}    

int
Newmark::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(8);
    data(0) = gamma;
    data(1) = beta;
    if (displ == true) 
      data(2) = 1.0;
    else
      data(2) = 0.0;

    data(3) = 1.0;	
    data(4) = alphaM;
    data(5) = betaK;
    data(6) = betaKi;
    data(7) = betaKc;
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
	opserr << "WARNING Newmark::sendSelf() - could not send data\n";
	return -1;
    }	
    return 0;
}

int
Newmark::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(8);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
	opserr << "WARNING Newmark::recvSelf() - could not receive data\n";
	gamma = 0.5; beta = 0.25; 
	return -1;
    }
    
    gamma = data(0);
    beta = data(1);
    if (data(2) == 1.0)
      displ = true;
    else
      displ = false;
    alphaM = data(4);
    betaK = data(5);
    betaKi = data(6);
    betaKc = data(7);
      
    return 0;
    
}

void
Newmark::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t Newmark - currentTime: " << currentTime;
	s << "  gamma: " << gamma << "  beta: " << beta << endln;
	s << " c1: " << c1 << " c2: " << c2 << " c3: " << c3 << endln;
	s << "  Rayleigh Damping - alphaM: " << alphaM;
	s << "  betaK: " << betaK << "   betaKi: " << betaKi << endln;	    
    } else 
	s << "\t Newmark - no associated AnalysisModel\n";
}

// AddingSensitivity:BEGIN //////////////////////////////
int
Newmark::revertToStart()
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

