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
                                                                        
// $Revision: 1.7 $
// $Date: 2003-02-14 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/HHT1.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/HHT1.C
// 
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the implementation of the HHT1 class.
//
// What: "@(#) HHT1.C, revA"

#include <HHT1.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

HHT1::HHT1()
:TransientIntegrator(INTEGRATOR_TAGS_HHT1),
 alpha(0.5), gamma(1.0), beta(0), 
 alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{
    
}

HHT1::HHT1(double _alpha)
:TransientIntegrator(INTEGRATOR_TAGS_HHT1),
 alpha(_alpha), gamma(1.5-_alpha), beta((2-_alpha)*(2-_alpha)*0.25),
 alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{

}

HHT1::HHT1(double _alpha, double alpham, double betak, double betaki, double betakc)
:TransientIntegrator(INTEGRATOR_TAGS_HHT1),
 alpha(_alpha), gamma(1.5-_alpha), beta((2-_alpha)*(2-_alpha)*0.25), 
 alphaM(alpham), betaK(betak), betaKi(betaki), betaKc(betakc),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{

}

HHT1::~HHT1()
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
  if (Udotalpha != 0)
    delete Udotalpha;
}


int
HHT1::initialize(void)
{
  return 0;
}


int
HHT1::newStep(double deltaT)
{

  if (beta == 0 || gamma == 0 ) {
    opserr << "HHT1::newStep() - error in variable\n";
    opserr << "gamma = " << gamma << " beta= " << beta << endln;
    return -1;
  }
    
  if (deltaT <= 0.0) {
    opserr << "HHT1::newStep() - error in variable\n";
    opserr << "dT = " << deltaT << endln;
    return -2;	
  }
  c1 = 1.0;
  c2 = gamma/(beta*deltaT);
  c3 = 1.0/(beta*deltaT*deltaT);


  AnalysisModel *theModel = this->getAnalysisModelPtr();

  if (U == 0) {
    opserr << "HHT1::newStep() - domainChange() failed or hasn't been called\n";
    return -3;	
  }

  // set response at t to be that at t+delta t of previous step
  (*Ut) = *U;        
  (*Utdot) = *Udot;  
  (*Utdotdot) = *Udotdot;  
    
  // set new velocity and accelerations at t + delta t
  double a1 = (1.0 - gamma/beta); 
  double a2 = (deltaT)*(1.0 - 0.5*gamma/beta);

  //  (*Udot) *= a1;
  Udot->addVector(a1, *Utdotdot,a2);

  double a3 = -1.0/(beta*deltaT);
  double a4 = 1 - 0.5/beta;
  //  (*Udotdot) *= a4;  
   Udotdot->addVector(a4, *Utdot,a3);

  (*Ualpha) = *Ut;
  (*Udotalpha) = *Utdot;
  //  (*Udotalpha) *= (1 - alpha);
  Udotalpha->addVector((1-alpha), *Udot, alpha);

  // set the new trial response quantities

  theModel->setResponse(*Ualpha,*Udotalpha,*Udotdot);        

  // increment the time and apply the load
  double time = theModel->getCurrentDomainTime();
  time +=deltaT;
  if (theModel->updateDomain(time, deltaT) < 0) {
    opserr << "HHT::newStep() - failed to update the domain\n";
    return -4;
  }

  return 0;
}



int
HHT1::formEleTangent(FE_Element *theEle)
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
HHT1::formNodTangent(DOF_Group *theDof)
{
  theDof->zeroTangent();
  theDof->addMtoTang(c3);
  theDof->addCtoTang(c2);        
  
  return(0);
}    



int 
HHT1::domainChanged()
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
    if (Ualpha != 0)
      delete Ualpha;
    if (Udotalpha != 0)
      delete Udotalpha;
    
    // create the new
    Ut = new Vector(size);
    Utdot = new Vector(size);
    Utdotdot = new Vector(size);
    U = new Vector(size);
    Udot = new Vector(size);
    Udotdot = new Vector(size);
    Ualpha = new Vector(size);
    Udotalpha = new Vector(size);

    // check we obtained the new
    if (Ut == 0 || Ut->Size() != size ||
	Utdot == 0 || Utdot->Size() != size ||
	Utdotdot == 0 || Utdotdot->Size() != size ||
	U == 0 || U->Size() != size ||
	Udot == 0 || Udot->Size() != size ||
	Udotdot == 0 || Udotdot->Size() != size ||
	Ualpha == 0 || Ualpha->Size() != size ||
	Udotalpha == 0 || Udotalpha->Size() != size) {
  
      opserr << "HHT1::domainChanged - ran out of memory\n";

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
    if (Udotalpha != 0)
      delete Udotalpha;

      Ut = 0; Utdot = 0; Utdotdot = 0;
      U = 0; Udot = 0; Udotdot = 0; Udotalpha=0; Ualpha =0;
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
HHT1::update(const Vector &deltaU)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    opserr << "WARNING HHT1::update() - no AnalysisModel set\n";
    return -1;
  }	

  // check domainChanged() has been called, i.e. Ut will not be zero
  if (Ut == 0) {
    opserr << "WARNING HHT1::update() - domainChange() failed or not called\n";
    return -2;
  }	

  // check deltaU is of correct size
  if (deltaU.Size() != U->Size()) {
    opserr << "WARNING HHT1::update() - Vectors of incompatable size ";
    opserr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endln;
    return -3;
  }
    
  //  determine the response at t+delta t
  (*U) += deltaU;
  Udot->addVector(1.0, deltaU,c2);
  Udotdot->addVector(1.0, deltaU,c3);
  Ualpha->addVector(1.0, deltaU, alpha);
  Udotalpha->addVector(1.0, deltaU, c2*alpha);
  
  // update the responses at the DOFs
  theModel->setResponse(*Ualpha,*Udotalpha,*Udotdot);        
  if (theModel->updateDomain() < 0) {
    opserr << "HHT1::update() - failed to update the domain\n";
    return -4;
  }
		    
  return 0;
}    

int
HHT1::commit(void)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    opserr << "WARNING HHT1::commit() - no AnalysisModel set\n";
    return -1;
  }	  

  // update the responses at the DOFs
  theModel->setResponse(*U,*Udot,*Udotdot);        
  theModel->updateDomain();

  return theModel->commitDomain();
}

int
HHT1::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(8);
    data(0) = alpha;
    data(1) = beta;
    data(2) = gamma;
    data(3) = 1.0;	
    data(4) = alphaM;
    data(5) = betaK;
    data(6) = betaKi;
    data(7) = betaKc;

    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
	opserr << "WARNING HHT1::sendSelf() - could not send data\n";
	return -1;
    }	
    return 0;
}

int
HHT1::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(8);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
	opserr << "WARNING HHT1::recvSelf() - could not receive data\n";
	return -1;
    }

    alpha = data(0);
    beta = data(1);
    gamma = data(2);
    alpha = data(0);
    beta = data(1);
    gamma = data(2);
    
    alphaM = data(4);
    betaK = data(5);
    betaKi = data(6);
    betaKc = data(7);
    
    return 0;
}

void
HHT1::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t HHT1 - currentTime: " << currentTime << " alpha: ";
	s << alpha << " gamma: " << gamma << "  beta: " << beta << endln;
	s << "  Rayleigh Damping - alphaM: " << alphaM;
	s << "  betaK: " << betaK << " betaKi: " << betaKi << endln;	    
    } else 
	s << "\t HHT1 - no associated AnalysisModel\n";
}

