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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/CentralDifference.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/CentralDifference.C
// 
// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the implementation of the CentralDifference class.
//
// What: "@(#) CentralDifference.C, revA"

#include <CentralDifference.h>
#include <FE_Element.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <AnalysisModel.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

CentralDifference::CentralDifference()
:TransientIntegrator(INTEGRATOR_TAGS_CentralDifference),
 c1(0.0), c2(0.0), c3(0.0), 
 Ut(0), Utm1(0),  U(0), Udot(0), Udotdot(0), deltaT(0)
{
    
}

CentralDifference::~CentralDifference()
{
  // clean up the memory created
  if (Ut != 0)
    delete Ut;
  if (Utm1 != 0)
    delete Utm1;
  if (U != 0)
    delete U;
  if (Udot != 0)
    delete Udot;
  if (Udotdot != 0)
    delete Udotdot;
}

int
CentralDifference::newStep(double _deltaT)
{
  
  deltaT = _deltaT;

  if (deltaT <= 0.0) {
    cerr << "Newton::newStep() - error in variable\n";
    cerr << "dT = " << deltaT << endl;
    return -2;	
  }

  c1 = 0.0;
  c2 = 0.5/deltaT;
  c3 = 1.0/(deltaT*deltaT);
    
  if (U == 0) {
    cerr << "Newton::newStep() - domainChange() failed or hasn't been called\n";
    return -3;	
  }

  // set disp at t and t-delta t
  (*Utm1) = *Ut;        
  (*Ut) = *U;  

  // set new velocity and accelerations at t 
  (*Udot) = *Ut;
  (*Udot) -= *Utm1;
  (*Udot) *= c2;

  (*Udotdot) = *Utm1;   
  (*Udotdot) -= *Ut;   
  (*Udotdot) *= c3;   
  Udotdot->addVector(1.0, *Ut, c3);
    

  // set the new trial response quantities
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  theModel->setResponse(*U,*Udot,*Udotdot);        

  // increment the time and apply the load
  double time = theModel->getCurrentDomainTime();
  theModel->applyLoadDomain(time);

  theModel->updateDomain();
  
  return 0;
}

int
CentralDifference::formEleTangent(FE_Element *theEle)
{
  theEle->zeroTangent();
  theEle->addCtoTang(c2);
  theEle->addMtoTang(c3);

  return 0;
}    

int
CentralDifference::formNodTangent(DOF_Group *theDof)
{
  theDof->zeroTangent();
  theDof->addMtoTang(c3);
  return(0);
}    

int 
CentralDifference::domainChanged()
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
    if (Utm1 != 0)
      delete Utm1;
    if (U != 0)
      delete U;
    if (Udot != 0)
      delete Udot;
    if (Udotdot != 0)
      delete Udotdot;
    
    // create the new
    Ut = new Vector(size);
    Utm1 = new Vector(size);
    U = new Vector(size);
    Udot = new Vector(size);
    Udotdot = new Vector(size);

    // cheack we obtained the new
    if (Ut == 0 || Ut->Size() != size ||
	Utm1 == 0 || Utm1->Size() != size ||
	U == 0 || U->Size() != size ||
	Udot == 0 || Udot->Size() != size ||
	Udotdot == 0 || Udotdot->Size() != size) {
      
      cerr << "CentralDifference::domainChanged - ran out of memory\n";

      // delete the old
      if (Ut != 0)
	delete Ut;
      if (Utm1 != 0)
	delete Utm1;
      if (U != 0)
	delete U;
      if (Udot != 0)
	delete Udot;
      if (Udotdot != 0)
	delete Udotdot;

      Ut = 0; Utm1 = 0; 
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
    const ID &id = dofPtr->getID();
    for (int i=0; i < id.Size(); i++) {
      int loc = id(i);
      if (loc >= 0) {
 	(*U)(loc) = disp(i);		
      }
    }
    // also need displacements at t-delta t
    cerr << "CentralDifference::update() - NEED disp at t-delta t\n";
  }    
  
  return 0;
}


int
CentralDifference::update(const Vector &deltaU)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    cerr << "WARNING CentralDifference::update() - no AnalysisModel set\n";
    return -1;
  }	
  
  // check domainChanged() has been called, i.e. Ut will not be zero
  if (Ut == 0) {
    cerr << "WARNING CentralDifference::update() - domainChange() failed or not called\n";
    return -2;
  }	

  // check deltaU is of correct size
  if (deltaU.Size() != U->Size()) {
    cerr << "WARNING CentralDifference::update() - Vectors of incompatable size ";
    cerr << " expecting " << U->Size() << " obtained " << deltaU.Size() << endl;
    return -3;
  }

  //  determine the displacement at t+delta t and others at time t
  (*U) = deltaU;
  Udot->addVector(1.0, deltaU,c2);
  Udotdot->addVector(1.0, deltaU,c3);
  
  // update the responses at the DOFs
  theModel->setResponse(*U,*Udot,*Udotdot);        
  theModel->updateDomain();
		    
  return 0;
}    

int
CentralDifference::commit(void)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    cerr << "WARNING CentralDifference::commit() - no AnalysisModel set\n";
    return -1;
  }	  
  
  // set response for time t, commit Domain and set time for next iteraration
  theModel->setResponse(*Ut, *Udot, *Udotdot);
  theModel->commitDomain();
  double time = theModel->getCurrentDomainTime() + deltaT;
  theModel->setCurrentDomainTime(time);

  return 0;
}

int
CentralDifference::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int
CentralDifference::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

void
CentralDifference::Print(ostream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t CentralDifference - currentTime: " << currentTime;
    } else 
	s << "\t CentralDifference - no associated AnalysisModel\n";
}

