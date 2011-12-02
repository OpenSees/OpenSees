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

Newmark1::Newmark1()
:TransientIntegrator(INTEGRATOR_TAGS_Newmark1),
 gamma(0), beta(0), 
 rayleighDamping(false), alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
 c1(0.0), c2(0.0), c3(0.0), c4(0.0),
 Up(0), Updot(0),  U(0), Udot(0), Udotdot(0)
{
    
}

Newmark1::Newmark1(double theGamma, double theBeta, bool dispFlag)
:TransientIntegrator(INTEGRATOR_TAGS_Newmark1),
 gamma(theGamma), beta(theBeta), 
 rayleighDamping(false), alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
 c1(0.0), c2(0.0), c3(0.0), c4(0.0),
 Up(0), Updot(0),  U(0), Udot(0), Udotdot(0)
{

}

Newmark1::Newmark1(double theGamma, double theBeta, 
		 double alpham, double betak, 
		   double betaki, double betakc)
:TransientIntegrator(INTEGRATOR_TAGS_Newmark1),
 gamma(theGamma), beta(theBeta), 
 rayleighDamping(true), 
 alphaM(alpham), betaK(betak), betaKi(betaki), betaKc(betakc),
 c1(0.0), c2(0.0), c3(0.0), c4(0.0),
 Up(0), Updot(0),  U(0), Udot(0), Udotdot(0)
{
    if (alpham == 0.0 && betak == 0.0 && betaki == 0.0 && betakc == 0.0)
	rayleighDamping = false;
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
Newmark1::formEleResidual(FE_Element *theEle)
{
  theEle->zeroResidual();
  if (rayleighDamping == false) {
      theEle->addRIncInertiaToResidual();
  } else {
      theEle->addRIncInertiaToResidual();
      theEle->addKtForce(*Udot, -betaK);
      theEle->addKcForce(*Udot, -betaKc);
      theEle->addKiForce(*Udot, -betaKi);
      theEle->addM_Force(*Udot, -alphaM);
  }    
  return 0;
}    

int
Newmark1::formNodUnbalance(DOF_Group *theDof)
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
Newmark1::initialize(void)
{

  // loop through the FE_Elements getting them to set Ki
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  FE_EleIter &theEles = theModel->getFEs();    
  FE_Element *elePtr;    
  while((elePtr = theEles()) != 0) 
    elePtr->setKi();
  return 0;
}

int
Newmark1::newStep(double deltaT)
{

  if (beta == 0 || gamma == 0 ) {
    cerr << "Newton::newStep() - error in variable\n";
    cerr << "gamma = " << gamma << " beta= " << beta << endl;
    return -1;
  }

  if (deltaT <= 0.0) {
    cerr << "Newmark1::newStep() - error in variable\n";
    cerr << "dT = " << deltaT << endl;
    return -2;	
  }

  // set the constants
  c1 = 1.0;
  c2 = gamma/(beta*deltaT);
  c3 = 1.0/(beta*deltaT*deltaT);
  c4 = gamma*deltaT;


  // set the new trial response quantities
  AnalysisModel *theModel = this->getAnalysisModelPtr();


  // if alphaKc is specified .. 
  //    loop over FE_Elements getting them to do a setKc()
  if (rayleighDamping == true && betaKc != 0.0) {
    FE_Element *elePtr;
    FE_EleIter &theEles = theModel->getFEs();    
    while((elePtr = theEles()) != 0)     
      elePtr->setKc();
  }

  if (U == 0) {
    cerr << "Newton::newStep() - domainChange() failed or hasn't been called\n";
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
  theModel->applyLoadDomain(time);

  theModel->updateDomain();
  
  return 0;
}

int
Newmark1::formEleTangent(FE_Element *theEle)
{
  theEle->zeroTangent();
  if (rayleighDamping == false) {
      theEle->addKtToTang(c1);
      theEle->addCtoTang(c2);
      theEle->addMtoTang(c3);
  } else {
      theEle->addKtToTang(c1 + c2*betaK);
      theEle->addMtoTang(c3 + c2*alphaM);
      theEle->addKiToTang(c2*betaKi);
      theEle->addKcToTang(c2*betaKc);
  }

  return 0;
}    


int
Newmark1::formNodTangent(DOF_Group *theDof)
{
  theDof->zeroTangent();
  if (rayleighDamping == false) 
      theDof->addMtoTang(c3);
  else
      theDof->addMtoTang(c3 + c2*alphaM);      

  return(0);
}    



int 
Newmark1::domainChanged()
{
  AnalysisModel *myModel = this->getAnalysisModelPtr();
  LinearSOE *theLinSOE = this->getLinearSOEPtr();
  const Vector &x = theLinSOE->getX();
  int size = x.Size();
  
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
      
      cerr << "Newmark1::domainChanged - ran out of memory\n";

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
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    cerr << "WARNING Newmark1::update() - no AnalysisModel set\n";
    return -1;
  }	

  // check domainChanged() has been called, i.e. Ut will not be zero
  if (U == 0) {
    cerr << "WARNING Newmark1::update() - domainChange() failed or not called\n";
    return -2;
  }	

  // check deltaU is of correct size
  if (deltaU.Size() != U->Size()) {
    cerr << "WARNING Newmark1::update() - Vectors of incompatable size ";
    cerr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endl;
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
  theModel->updateDomain();
		    
  return 0;
}    

int
Newmark1::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(7);
    data(0) = gamma;
    data(1) = beta;
    if (rayleighDamping == true) {
	data(2) = 1.0;	
	data(3) = alphaM;
	data(4) = betaK;
	data(5) = betaKi;
	data(6) = betaKc;
    } else
	data(2) = 0.0;	
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
	cerr << "WARNING Newmark1::sendSelf() - could not send data\n";
	return -1;
    }	
    return 0;
}

int
Newmark1::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(7);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
	cerr << "WARNING Newmark1::recvSelf() - could not receive data\n";
	gamma = 0.5; beta = 0.25; rayleighDamping = false;
	return -1;
    }
    
    gamma = data(0);
    beta = data(1);
    if (data(2) == 1.0) {
	rayleighDamping = true;
	alphaM = data(3);
	betaK = data(4);
	betaKi = data(5);
	betaKc = data(6);
    } else
	rayleighDamping = false;	
      
    return 0;
    
}

void
Newmark1::Print(ostream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t Newmark1 - currentTime: " << currentTime;
	s << "  gamma: " << gamma << "  beta: " << beta << endl;
	s << " c1: " << c1 << " c2: " << c2 << " c3: " << c3 << endl;
	if (rayleighDamping == true) {
	    s << "  Rayleigh Damping - alphaM: " << alphaM;
	    s << "  betaK: " << betaK << endl;	    
	}
    } else 
	s << "\t Newmark1 - no associated AnalysisModel\n";
}




