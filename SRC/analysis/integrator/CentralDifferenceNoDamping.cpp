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
                                                                        
// $Revision$
// $Date$
// $URL$

// Written: fmk 
// Created: 11/98
//
// Description: This file contains the implementation of the CentralDifferenceNoDamping 
// class.

#include <CentralDifferenceNoDamping.h>
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
#include <elementAPI.h>
#define OPS_Export 


void *    OPS_CentralDifferenceNoDamping(void)
{
    // pointer to an integrator that will be returned
    TransientIntegrator *theIntegrator = 0;
    
    theIntegrator = new CentralDifferenceNoDamping();
    
    if (theIntegrator == 0)
        opserr << "WARNING - out of memory creating CentralDifferenceNoDamping integrator\n";
    
    return theIntegrator;
}


CentralDifferenceNoDamping::CentralDifferenceNoDamping()
:TransientIntegrator(INTEGRATOR_TAGS_CentralDifferenceNoDamping),
 updateCount(0), 
 U(0), Udot(0), Udotdot(0), deltaT(0)
{
    
}

CentralDifferenceNoDamping::~CentralDifferenceNoDamping()
{
  // clean up the memory created
  if (U != 0)
    delete U;
  if (Udot != 0)
    delete Udot;
  if (Udotdot != 0)
    delete Udotdot;
}

int
CentralDifferenceNoDamping::newStep(double _deltaT)
{
  updateCount = 0;
  
  deltaT = _deltaT;

  if (deltaT <= 0.0) {
    opserr << "CentralDifference::newStep() - error in variable\n";
    opserr << "dT = " << deltaT << endln;
    return -2;	
  }

  AnalysisModel *theModel = this->getAnalysisModel();
  double time = theModel->getCurrentDomainTime();
  theModel->applyLoadDomain(time);

  return 0;
}

int
CentralDifferenceNoDamping::formEleTangent(FE_Element *theEle)
{
  theEle->zeroTangent();
  theEle->addMtoTang();

  return 0;
}    

int
CentralDifferenceNoDamping::formNodTangent(DOF_Group *theDof)
{
  theDof->zeroTangent();
  theDof->addMtoTang();
  return(0);
}    

int
CentralDifferenceNoDamping::formEleResidual(FE_Element *theEle)
{
  theEle->zeroResidual();
  theEle->addRtoResidual();
  return 0;
}    

int
CentralDifferenceNoDamping::formNodUnbalance(DOF_Group *theDof)
{
  theDof->zeroUnbalance();
  theDof->addPtoUnbalance();
  return 0;
}    


int 
CentralDifferenceNoDamping::domainChanged()
{
  AnalysisModel *myModel = this->getAnalysisModel();
  LinearSOE *theLinSOE = this->getLinearSOE();
  const Vector &x = theLinSOE->getX();
  int size = x.Size();
  
  // create the new Vector objects
  if (U == 0 || U->Size() != size) {

    // delete the old
    if (U != 0)
      delete U;
    if (Udot != 0)
      delete Udot;
    if (Udotdot != 0)
      delete Udotdot;

    // create the new
    U = new Vector(size);
    Udot = new Vector(size);
    Udotdot = new Vector(size);

    // cheack we obtained the new
    if (U == 0 || U->Size() != size ||
	Udot == 0 || Udot->Size() != size ||
	Udotdot == 0 || Udotdot->Size() != size) {
      
      opserr << "CentralDifferenceNoDamping::domainChanged - ran out of memory\n";

      // delete the old
      if (U != 0)
	delete U;
      if (Udot != 0)
	delete U;
      if (Udotdot != 0)
	delete Udot;

      U = 0; Udot = 0; Udotdot = 0;
      return -1;
    }
  }        
    
  // now go through and populate U and Udot by iterating through
  // the DOF_Groups and getting the last committed velocity and accel

  DOF_GrpIter &theDOFs = myModel->getDOFs();
  DOF_Group *dofPtr;
    
  while ((dofPtr = theDOFs()) != 0) {
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
  }    
  
  return 0;
}


int
CentralDifferenceNoDamping::update(const Vector &X)
{
  updateCount++;
  if (updateCount > 1) {
    opserr << "ERROR CentralDifferenceNoDamping::update() - called more than once -";
    opserr << " Central Difference integraion schemes require a LINEAR solution algorithm\n";
    return -1;
  }
  
  AnalysisModel *theModel = this->getAnalysisModel();

  if (theModel == 0) {
    opserr << "ERROR CentralDifferenceNoDamping::update() - no AnalysisModel set\n";
    return -2;
  }	
  
  // check domainChanged() has been called, i.e. Ut will not be zero
  if (U == 0) {
    opserr << "WARNING CentralDifferenceNoDamping::update() - domainChange() failed or not called\n";
    return -2;
  }	

  // check deltaU is of correct size
  if (X.Size() != U->Size()) {
    opserr << "WARNING CentralDifferenceNoDamping::update() - Vectors of incompatible size ";
    opserr << " expecting " << U->Size() << " obtained " << X.Size() << endln;
    return -3;
  }

  //  determine the acceleration at time t 
  (*Udotdot) = X;

  //  determine the vel at t+ 0.5 * delta t 
  Udot->addVector(1.0, X, deltaT);
  
  //  determine the displacement at t+delta t 
  U->addVector(1.0, *Udot, deltaT);

  // update the disp & responses at the DOFs
  theModel->setDisp(*U);
  theModel->updateDomain();

  return 0;
}    

int
CentralDifferenceNoDamping::commit(void)
{
  AnalysisModel *theModel = this->getAnalysisModel();
  if (theModel == 0) {
    opserr << "WARNING CentralDifferenceNoDamping::commit() - no AnalysisModel set\n";
    return -1;
  }	  
  
  // update time in Domain to T + deltaT & commit the domain
  double time = theModel->getCurrentDomainTime() + deltaT;
  theModel->setCurrentDomainTime(time);

  return theModel->commitDomain();

  return 0;
}

int
CentralDifferenceNoDamping::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int
CentralDifferenceNoDamping::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

void
CentralDifferenceNoDamping::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModel();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t CentralDifferenceNoDamping - currentTime: " << currentTime;
    } else 
	s << "\t CentralDifferenceNoDamping - no associated AnalysisModel\n";
}

