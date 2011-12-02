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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/HHT.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/HHT.C
// 
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the implementation of the HHT class.
//
// What: "@(#) HHT.C, revA"

#include <HHT.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

HHT::HHT()
:TransientIntegrator(INTEGRATOR_TAGS_HHT),
 alpha(0.5), gamma(1.0), beta(0), 
 rayleighDamping(false), alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{
    
}

HHT::HHT(double _alpha)
:TransientIntegrator(INTEGRATOR_TAGS_HHT),
 alpha(_alpha), gamma(1.5-_alpha), beta((2-_alpha)*(2-_alpha)*0.25),
 rayleighDamping(false), alphaM(0.0), betaK(0.0), betaKi(0.0), betaKc(0.0),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{

}

HHT::HHT(double _alpha, double alpham, double betak, double betaki, double betakc)
:TransientIntegrator(INTEGRATOR_TAGS_HHT),
 alpha(_alpha), gamma(1.5-_alpha), beta((2-_alpha)*(2-_alpha)*0.25),  
 rayleighDamping(true), 
 alphaM(alpham), betaK(betak), betaKi(betaki), betaKc(betakc),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utdot(0), Utdotdot(0),  U(0), Udot(0), Udotdot(0),
 Ualpha(0),Udotalpha(0)
{
    if (alpham == 0.0 && betak == 0.0 && betaki == 0.0 && betakc == 0.0)
	rayleighDamping = false;
}

HHT::~HHT()
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
HHT::formEleResidual(FE_Element *theEle)
{
  theEle->zeroResidual();
  if (rayleighDamping == false) {
      theEle->addRIncInertiaToResidual();
  } else {
      theEle->addRIncInertiaToResidual();
      theEle->addKtForce(*Udot,-betaK);
      theEle->addKcForce(*Udot, -betaKc);
      theEle->addKiForce(*Udot, -betaKi);
      theEle->addM_Force(*Udot,-alphaM);
  }    
  return 0;
}    

int
HHT::formNodUnbalance(DOF_Group *theDof)
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
HHT::initialize(void)
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
HHT::newStep(double deltaT)
{

  if (beta == 0 || gamma == 0 ) {
    cerr << "HHT::newStep() - error in variable\n";
    cerr << "gamma = " << gamma << " beta= " << beta << endl;
    return -1;
  }
    
  if (deltaT <= 0.0) {
    cerr << "HHT::newStep() - error in variable\n";
    cerr << "dT = " << deltaT << endl;
    return -2;	
  }
  c1 = 1.0;
  c2 = gamma/(beta*deltaT);
  c3 = 1.0/(beta*deltaT*deltaT);

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
    cerr << "HHT::newStep() - domainChange() failed or hasn't been called\n";
    return -3;	
  }

  // set response at t to be that at t+delta t of previous step
  (*Ut) = *U;        
  (*Utdot) = *Udot;  
  (*Utdotdot) = *Udotdot;  
    
  // set new velocity and accelerations at t + delta t
  double a1 = (1.0 - gamma/beta); 

  double a2 = (deltaT)*(1.0 - 0.5*gamma/beta);
  // (*Udot) *= a1;  
  Udot->addVector(a1,*Utdotdot,a2);

  double a3 = -1.0/(beta*deltaT);
  double a4 = 1 - 0.5/beta;
  //  (*Udotdot) *= a4;  
   Udotdot->addVector(a4,*Utdot,a3);

  (*Ualpha) = *Ut;
  (*Udotalpha) = *Utdot;
  //  (*Udotalpha) *= (1 - alpha);
  Udotalpha->addVector((1-alpha),*Udot, alpha);

  // set the new trial response quantities
  theModel->setResponse(*Ualpha,*Udotalpha,*Udotdot);        

  // increment the time and apply the load
  double time = theModel->getCurrentDomainTime();
  time +=deltaT;
  theModel->applyLoadDomain(time);

  theModel->updateDomain();

  return 0;
}



int
HHT::formEleTangent(FE_Element *theEle)
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
HHT::formNodTangent(DOF_Group *theDof)
{
  theDof->zeroTangent();
  if (rayleighDamping == false) 
      theDof->addMtoTang(c3);
  else
      theDof->addMtoTang(c3 + c2*alphaM);        
  
  return(0);
}    



int 
HHT::domainChanged()
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
  
      cerr << "HHT::domainChanged - ran out of memory\n";

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
HHT::update(const Vector &deltaU)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    cerr << "WARNING HHT::update() - no AnalysisModel set\n";
    return -1;
  }	

  // check domainChanged() has been called, i.e. Ut will not be zero
  if (Ut == 0) {
    cerr << "WARNING HHT::update() - domainChange() failed or not called\n";
    return -2;
  }	

  // check deltaU is of correct size
  if (deltaU.Size() != U->Size()) {
    cerr << "WARNING HHT::update() - Vectors of incompatable size ";
    cerr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endl;
    return -3;
  }
    
  //  determine the response at t+delta t
  (*U) += deltaU;
  Udot->addVector(1.0,deltaU,c2);
  Udotdot->addVector(1.0,deltaU,c3);
  Ualpha->addVector(1.0,deltaU, alpha);
  Udotalpha->addVector(1.0,deltaU, c2*alpha);
  
  // update the responses at the DOFs
  theModel->setResponse(*Ualpha,*Udotalpha,*Udotdot);        
  theModel->updateDomain();
		    
  return 0;
}    

int
HHT::commit(void)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    cerr << "WARNING HHT::commit() - no AnalysisModel set\n";
    return -1;
  }	  

  // update the responses at the DOFs
  theModel->setResponse(*U,*Udot,*Udotdot);        
  theModel->updateDomain();

  return theModel->commitDomain();
}

int
HHT::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(8);
    data(0) = alpha;
    data(1) = beta;
    data(2) = gamma;
    if (rayleighDamping == true) {
	data(3) = 1.0;	
	data(4) = alphaM;
	data(5) = betaK;
	data(6) = betaKi;
	data(7) = betaKc;
    } else
	data(3) = 0.0;	
    
    if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
	cerr << "WARNING HHT::sendSelf() - could not send data\n";
	return -1;
    }	
    return 0;
}

int
HHT::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(8);
    if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
	cerr << "WARNING HHT::recvSelf() - could not receive data\n";
	return -1;
    }

    alpha = data(0);
    beta = data(1);
    gamma = data(2);
    
    if (data(3) == 1.0) {
	rayleighDamping = true;
	alphaM = data(4);
	betaK = data(5);
	betaKi = data(6);
	betaKc = data(7);
    } else
	rayleighDamping = false;	
    
    return 0;
}

void
HHT::Print(ostream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t HHT - currentTime: " << currentTime << " alpha: ";
	s << alpha << " gamma: " << gamma << "  beta: " << beta << endl;
	if (rayleighDamping == true) {
	    s << "  Rayleigh Damping - alphaM: " << alphaM;
	    s << "  betaK: " << betaK << endl;	    
	}
    } else 
	s << "\t HHT - no associated AnalysisModel\n";
}

