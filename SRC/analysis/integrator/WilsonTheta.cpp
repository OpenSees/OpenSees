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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/WilsonTheta.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/WilsonTheta.C
// 
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the implementation of WilsonTheta.
//
// What: "@(#) WilsonTheta.C, revA"

#include <WilsonTheta.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

WilsonTheta::WilsonTheta()
:TransientIntegrator(INTEGRATOR_TAGS_WilsonTheta),
 theta(0), deltaT(0),
 rayleighDamping(false), alphaM(0), betaK(0), 
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0)
{
    
}

WilsonTheta::WilsonTheta(double _theta)
:TransientIntegrator(INTEGRATOR_TAGS_WilsonTheta),
 theta(_theta), deltaT(0.0), 
 rayleighDamping(false), alphaM(0), betaK(0), 
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0)
{
    
}

WilsonTheta::WilsonTheta(double _theta, double alpham, double betak)
:TransientIntegrator(INTEGRATOR_TAGS_WilsonTheta),
 theta(_theta), deltaT(0.0), 
 rayleighDamping(true), alphaM(alpham), betaK(betak), 
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0)
{
    
}

WilsonTheta::~WilsonTheta()
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
WilsonTheta::formEleResidual(FE_Element *theEle)
{
  theEle->zeroResidual();
  if (rayleighDamping == false) {
      theEle->addRIncInertiaToResidual();
  } else {
      theEle->addRIncInertiaToResidual();
      theEle->addKtForce(*Udot,-betaK);
      theEle->addM_Force(*Udot,-alphaM);
  }    
  return 0;
}    

int
WilsonTheta::formNodUnbalance(DOF_Group *theDof)
{
  theDof->zeroUnbalance();
  if (rayleighDamping == false) 
      theDof->addPIncInertiaToUnbalance();
  else {
      theDof->addPIncInertiaToUnbalance();
      theDof->addM_Force(*Udot,-alphaM);
  }
  return 0;
}    

int
WilsonTheta::newStep(double _deltaT)
{
  if (theta <= 0.0 ) {
    cerr << "Newton::newStep() - error in variable\n";
    cerr << "theta: " << theta << " <= 0.0\n";
    return -1;
  }

  
  if (_deltaT <= 0.0) {
    cerr << "Newton::newStep() - error in variable\n";
    cerr << "dT = " << deltaT << endl;
    return -2;	
  }

  deltaT = _deltaT;
  c1 = 1.0;
  c2 = 3.0/(theta*deltaT);
  c3 = 2*c2/(theta*deltaT);
    
  if (U == 0) {
    cerr << "Newton::newStep() - domainChange() failed or hasn't been called\n"; 
    return -3;	
  }

  // set response at t to be that at t+delta t of previous step
  (*Ut) = *U;        
  (*Utdot) = *Udot;  
  (*Utdotdot) = *Udotdot;  
    
  // set new velocity and accelerations at t + theta delta t
  //  (*Udot) *= -2.0;
  double a1 = -0.5*theta*deltaT; 
  Udot->addVector(-2.0,*Utdotdot,a1);

  //  (*Udotdot) *= -2.0;  
  double a2 = -6.0/theta/deltaT; 
  Udotdot->addVector(-2.0, *Utdot, a2);
    

  // set the new trial response quantities
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  theModel->setResponse(*U,*Udot,*Udotdot);        

  // increment the time and apply the load
  double time = theModel->getCurrentDomainTime();
  time += theta * deltaT;
  theModel->applyLoadDomain(time);

  theModel->updateDomain();
  
  return 0;
}

int
WilsonTheta::formEleTangent(FE_Element *theEle)
{
  theEle->zeroTangent();
  if (rayleighDamping == false) {
      theEle->addKtToTang(c1);
      theEle->addCtoTang(c2);
      theEle->addMtoTang(c3);
  } else {
      theEle->addKtToTang(c1 + c2*betaK);
      theEle->addMtoTang(c3 + c2*alphaM);
  }  

  return 0;
}    

int
WilsonTheta::formNodTangent(DOF_Group *theDof)
{
  theDof->zeroTangent();
  if (rayleighDamping == false) 
      theDof->addMtoTang(c3);
  else
      theDof->addMtoTang(c3 + c2*alphaM);        

  return(0);
}    

int 
WilsonTheta::domainChanged()
{
  AnalysisModel *myModel = this->getAnalysisModelPtr();
  LinearSOE *theLinSOE = this->getLinearSOEPtr();
  const Vector &x = theLinSOE->getX();
  int size = x.Size();
  
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
      
      cerr << "WilsonTheta::domainChanged - ran out of memory\n";

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
    const Vector &disp = dofPtr->getCommittedDisp();	
    const Vector &vel = dofPtr->getCommittedVel();
    const Vector &accel = dofPtr->getCommittedAccel();	
    const ID &id = dofPtr->getID();
    for (int i=0; i < id.Size(); i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*U)(loc) = disp(i);		
 	(*Udot)(loc) = vel(i);
 	(*Udotdot)(loc) = accel(i);
      }
    }
  }    
  
  return 0;
}


int
WilsonTheta::update(const Vector &deltaU)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    cerr << "WARNING WilsonTheta::update() - no AnalysisModel set\n";
    return -1;
  }	

  // check domainChanged() has been called, i.e. Ut will not be zero
  if (Ut == 0) {
    cerr << "WARNING WilsonTheta::update() - domainChange() failed or not called\n";
    return -2;
  }	
  
  // check deltaU is of correct size
  if (deltaU.Size() != U->Size()) {
    cerr << "WARNING WilsonTheta::update() - Vectors of incompatable size ";
    cerr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endl;
    return -3;
  }

  // determine the response at t + theta * deltaT
  (*U) += deltaU;
  Udot->addVector(1.0, deltaU,c2);
  Udotdot->addVector(1.0, deltaU,c3);
  
  // update the responses at the DOFs
  theModel->setResponse(*U,*Udot,*Udotdot);        
  theModel->updateDomain();
  
  return 0;
}    


int
WilsonTheta::commit(void)
{

  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    cerr << "WARNING WilsonTheta::commit() - no AnalysisModel set\n";
    return -1;
  }	  
  
  // determine the quantities at t + delta t 
  /* THIS IS WHAT IS IN BATHE'S BOOK - SEE BELOW FOR QUICKER SOLN
  double a1 = c3/theta;
  double a2 = c3*deltaT;	
  double a3 = (1 - 3.0/theta);
  (*Udotdot) = *U;
  (*Udotdot) -= *Ut;
  (*Udotdot) *= a1;
  Udotdot->addVector(*Utdot, -a2);
  Udotdot->addVector(*Utdotdot, a3);
  */
  
  /* REVISION BASED ON WHAT IS IN CHOPRA's BOOK - MAKES SENSE - LINEAR ACCEL */

  // determine response quantities at Ut+dt given Ut and Ut + theta dt
  double a1,a2;
  (*Udotdot) -= *Utdotdot;
  (*Udotdot) *= 1/theta;
  (*Udotdot) += *Utdotdot;
  
  (*Udot) = *Utdot;
  a1 = 0.5*deltaT;
  Udot->addVector(1.0, *Udotdot, a1);
  Udot->addVector(1.0, *Utdotdot, a1);


  (*U) = *Ut;
  U->addVector(1.0, *Utdot, deltaT);
  a2 = deltaT*deltaT/6.0;
  U->addVector(1.0, *Udotdot, a2);
  U->addVector(1.0, *Utdotdot, 2*a2);
  // update the t + delta t responses at the DOFs 
  theModel->setResponse(*U,*Udot,*Udotdot);        
  theModel->updateDomain();

  // set the time to be t+delta t
  double time = theModel->getCurrentDomainTime();
  time -= (theta -1) * deltaT;
  theModel->setCurrentDomainTime(time);

  // now invoke commit on the AnalysisModel
  return theModel->commitDomain();
}

int
WilsonTheta::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(4);
    data(0) = theta;
    if (rayleighDamping == true) {
	data(1) = 1.0;	
	data(2) = alphaM;
	data(3) = betaK;
    } else
	data(1) = 0.0;	    

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
	cerr << "WilsonTheta::sendSelf() - failed to send the data\n";
	return -1;
    }
    return 0;
}

int
WilsonTheta::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(1);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
	cerr << "WilsonTheta::recvSelf() - ";
	cerr << " failed to receive the Vector\n";
	theta = 1.0; rayleighDamping = false;
	return -1;
    }

    theta = data(0);
    if (data(1) == 1.0) {
	rayleighDamping = true;
	alphaM = data(2);
	betaK = data(3);
    } else
	rayleighDamping = false;	    

    return 0;
    
}

void
WilsonTheta::Print(ostream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t WilsonTheta - currentTime: " << currentTime;
	s << " theta: " << theta << endl;
	if (rayleighDamping == true) {
	    s << "  Rayleigh Damping - alphaM: " << alphaM;
	    s << "  betaK: " << betaK << endl;	    
	}	
    } else 
	s << "\t WilsonTheta - no associated AnalysisModel\n";
}

