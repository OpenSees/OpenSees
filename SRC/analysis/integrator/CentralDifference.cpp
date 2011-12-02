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
                                                                        
// $Revision: 1.4 $
// $Date: 2005-02-28 20:39:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/CentralDifference.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/98
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
 Utm1(0),  U(0), Udot(0), Udotdot(0), Y(0), Z(0), deltaT(0)
{
    
}

CentralDifference::~CentralDifference()
{
  // clean up the memory created
  if (Utm1 != 0)
    delete Utm1;
  if (U != 0)
    delete U;
  if (Udot != 0)
    delete Udot;
  if (Udotdot != 0)
    delete Udotdot;
  if (Y != 0)
    delete Y;
  if (Z != 0)
    delete Z;
}

int
CentralDifference::newStep(double _deltaT)
{
  updateCount = 0;
  
  deltaT = _deltaT;

  if (deltaT <= 0.0) {
    opserr << "CentralDifference::newStep() - error in variable\n";
    opserr << "dT = " << deltaT << endln;
    return -2;	
  }

  AnalysisModel *theModel = this->getAnalysisModelPtr();
  double time = theModel->getCurrentDomainTime();
  theModel->applyLoadDomain(time);

  c1 = 0.0;
  c2 = 0.5/deltaT;
  c3 = 1.0/(deltaT*deltaT);

  // determine the garbage quantities needed for vel & accel at nodes
  Y->addVector(0.0, *Utm1, -c2);
  
  Z->addVector(0.0, *U, -2.0);
  *Z += *Utm1;
  *Z *= c3;

  // update the responses at the DOFs; note that the set vel & accel are garbage
  theModel->setVel(*Y);
  theModel->setAccel(*Z);

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
  if (U == 0 || U->Size() != size) {

    if (Utm1 != 0)
      delete Utm1;
    if (U != 0)
      delete U;
    if (Udot != 0)
      delete Udot;
    if (Udotdot != 0)
      delete Udotdot;
    if (Y != 0)
      delete Y;
    if (Z != 0)
      delete Z;
    
    // create the new
    Utm1 = new Vector(size);
    U = new Vector(size);
    Udot = new Vector(size);
    Udotdot = new Vector(size);
    Y = new Vector(size);
    Z = new Vector(size);

    // cheack we obtained the new
    if (Utm1 == 0 || Utm1->Size() != size ||
	U == 0 || U->Size() != size ||
	Udot == 0 || Udot->Size() != size ||
	Udotdot == 0 || Udotdot->Size() != size ||
	Y == 0 || Y->Size() != size ||
	Z == 0 || Z->Size() != size) {
      
      opserr << "CentralDifference::domainChanged - ran out of memory\n";

      // delete the old
      if (Utm1 != 0)
	delete Utm1;
      if (U != 0)
	delete U;
      if (Udot != 0)
	delete Udot;
      if (Udotdot != 0)
	delete Udotdot;
      if (Y != 0)
	delete Y;
      if (Z != 0)
	delete Z;

      Utm1 = 0; 
      U = 0; 
      Y = 0;       Z = 0;
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
  }    

  opserr << "WARNING: CentralDifference::domainChanged() - assuming Ut-1 = 0\n";
  Utm1->Zero();

  return 0;
}


int
CentralDifference::update(const Vector &X)
{
  updateCount++;
  if (updateCount > 1) {
    opserr << "ERROR CentralDifference::update() - called more than once -";
    opserr << " Central Difference integration schemes require a LINEAR solution algorithm\n";
    return -1;
  }

  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    opserr << "WARNING CentralDifference::update() - no AnalysisModel set\n";
    return -1;
  }	
  
  // check domainChanged() has been called, i.e. Ut will not be zero
  if (U == 0) {
    opserr << "WARNING CentralDifference::update() - domainChange() failed or not called\n";
    return -2;
  }	


  // check deltaU is of correct size
  if (X.Size() != U->Size()) {
    opserr << "WARNING CentralDifference::update() - Vectors of incompatable size ";
    opserr << " expecting " << U->Size() << " obtained " << X.Size() << endln;
    return -3;
  }


  *Utm1 = *U;
  *U = X;

  *Udot = X;
  *Udot -= *Utm1;
  *Udot *= c2;

  Udotdot->addVector(0.0, *Utm1, 1.0);
  Udotdot->addVector(1.0, X, 1.0);
  Udotdot->addVector(1.0, *U, -2.0);
  *Udotdot *= c3;

  theModel->setResponse(*U, *Udot, *Udotdot);
  theModel->updateDomain();
		    
  return 0;
}    

int
CentralDifference::commit(void)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  if (theModel == 0) {
    opserr << "WARNING CentralDifference::commit() - no AnalysisModel set\n";
    return -1;
  }	  
  
  // set response for time t, commit Domain and set time for next iteraration
  double time = theModel->getCurrentDomainTime() + deltaT;
  theModel->setCurrentDomainTime(time);

  theModel->commitDomain();

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
CentralDifference::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentTime = theModel->getCurrentDomainTime();
	s << "\t CentralDifference - currentTime: " << currentTime;
    } else 
	s << "\t CentralDifference - no associated AnalysisModel\n";
}

