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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/LoadPath.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/LoadPath.h
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for LoadPath.
// LoadPath is an algorithmic class for perfroming a static analysis
// using a load control integration scheme.
//
// What: "@(#) LoadPath.h, revA"


#include <LoadPath.h>
#include <AnalysisModel.h>
#include <Vector.h>
#include <Channel.h>
#include <ID.h>
#include <stdlib.h>

LoadPath::LoadPath(Vector &theLoadPath)
:StaticIntegrator(INTEGRATOR_TAGS_LoadPath),
 loadPath(0), currentStep(0)
{
    loadPath = new Vector(theLoadPath);
    if (loadPath == 0 || loadPath->Size() == 0) {
	opserr << "LoadPath::LoadPath() - ran out of memory\n";
	exit(-1);
    }
}

LoadPath::LoadPath()
:StaticIntegrator(INTEGRATOR_TAGS_LoadPath),
 loadPath(0), currentStep(0)
{

}




LoadPath::~LoadPath()
{
    if (loadPath != 0)
	delete loadPath;
}

int 
LoadPath::newStep(void)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();    
    if (theModel == 0) {
	opserr << "LoadPath::newStep() - no associated AnalysisModel\n";
	return -1;
    }
    
    if (loadPath == 0) {
	opserr << "LoadPath::newStep() - no load path associated with object\n";
	return -2;
    }	
	

    double modelLambda = theModel->getCurrentDomainTime();

    double currentLambda;
    if (currentStep < loadPath->Size()) {
      
      if (currentStep > 0) {
	if (modelLambda == (*loadPath)(currentStep-1))
	  currentLambda = (*loadPath)(currentStep);  
        else
	  currentLambda = (*loadPath)(currentStep-1);  
      } else
	  currentLambda = (*loadPath)(currentStep);  
    }      
    else {
	currentLambda = 0.0;
	opserr << "LoadPath::newStep() - reached end of specified load path";
	opserr << " - setting lambda = 0.0 \n";
    }
    
    currentStep++;
    theModel->applyLoadDomain(currentLambda);
    
    return 0;
}
    
int
LoadPath::update(const Vector &deltaU)
{
    AnalysisModel *myModel = this->getAnalysisModelPtr();
    if (myModel == 0) {
	opserr << "WARNING LoadPath::update() ";
	opserr << "No AnalysisModel has been set\n";
	return -1;
    }

    myModel->incrDisp(deltaU);    
    myModel->updateDomain();
    return 0;
}


int
LoadPath::sendSelf(int cTag,
		   Channel &theChannel)
{
  ID data(2);
  data(0) = loadPath->Size();
  data(1) = currentStep;  
  if (theChannel.sendID(this->getDbTag(), cTag, data) < 0) {
      opserr << "LoadPath::sendSelf() - failed to send the ID\n";
      return -1;
  }
  
  if (theChannel.sendVector(this->getDbTag(), cTag, *loadPath) < 0) {
      opserr << "LoadPath::sendSelf() - failed to send the Vector\n";
      return -1;
  }  
  
  return 0;
}


int
LoadPath::recvSelf(int cTag,
		      Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  ID data(2);
  if (theChannel.recvID(this->getDbTag(), cTag, data) < 0) {
      opserr << "LoadPath::sendSelf() - failed to send the ID\n";
      return -1;
  }      
  int size = data(0);
  currentStep = data(1);
  
  loadPath = new Vector(size);
  if (loadPath == 0 || loadPath->Size() == 0) {
      opserr << "FATAL - LoadPath::recvSelf() - ran out of memory\n";
      exit(-1);
  }

  if (theChannel.recvVector(this->getDbTag(), cTag, *loadPath) < 0) {
      opserr << "LoadPath::sendSelf() - failed to send the Vector\n";
      return -1;
  }      
  
  return 0;
}



void
LoadPath::Print(OPS_Stream &s, int flag)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    if (theModel != 0) {
	double currentLambda = theModel->getCurrentDomainTime();
	s << "\t LoadPath - currentLambda: " << currentLambda << endln;
    } else 
	s << "\t LoadPath - no associated AnalysisModel\n";
    
}

