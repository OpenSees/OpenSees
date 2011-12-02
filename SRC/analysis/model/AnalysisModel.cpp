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
                                                                        
// $Revision: 1.8 $
// $Date: 2003-02-14 23:00:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/AnalysisModel.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/model/AnalysisModel.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Purpose: This file contains the class definition for AnalysisModel
// AnalysisModel is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints. These objects are all added to the AnalysisModel by a 
// ModelBuilder.
//
// What: "@(#) AnalysisModel.C, revA"

#include <stdlib.h>

#include <AnalysisModel.h>
#include <Domain.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <FE_EleIter.h>
#include <DOF_Graph.h>
#include <DOF_GroupGraph.h>
#include <Node.h>
#include <NodeIter.h>

//  AnalysisModel();
//	constructor

AnalysisModel::AnalysisModel(int theClassTag)
:MovableObject(theClassTag),
 myDomain(0), myDOFGraph(0), myGroupGraph(0),
 numFE_Ele(0), numDOF_Grp(0), numEqn(0),
 theFEs(0), theDOFs(0),
 theFEiter(*this), theDOFiter(*this)
{
    // for subclasses to use - they provide own container stuff
} 

AnalysisModel::AnalysisModel()
:MovableObject(AnaMODEL_TAGS_AnalysisModel),
 myDomain(0), myDOFGraph(0), myGroupGraph(0),
 numFE_Ele(0), numDOF_Grp(0), numEqn(0),
 theFEs(0), theDOFs(0),
 theFEiter(*this), theDOFiter(*this)
{
  sizeEle = 256;  // these arrays get enlarged as needed
  sizeDOF = 256;

  theFEs = new FE_Element *[sizeEle];
  if (theFEs == 0) {
      opserr << "FATAL:AnalysisModel::AnalysisModel()";
      opserr << " ran out of memory creating array for FE_Elements\n";
      exit(-1);
  }
  for (int i=0; i<sizeEle; i++) theFEs[i] = 0;

  theDOFs = new DOF_Group *[sizeDOF];
  if (theDOFs == 0) {
      opserr << "FATAL:AnalysisModel::AnalysisModel()";
      opserr << " ran out of memory creating array for DOF_Groups\n";
      exit(-1);
  }
  for (int j=0; j<sizeDOF; j++) theDOFs[j] = 0;
} 


// ~AnalysisModel();    
AnalysisModel::~AnalysisModel()
{
    if (theFEs != 0)
	delete [] theFEs;

    if (theDOFs != 0)
	delete [] theDOFs;

    if (myGroupGraph != 0) {
	delete myGroupGraph;    
    }	

    if (myDOFGraph != 0) {
	delete myDOFGraph;
    }
}    

void
AnalysisModel::setLinks(Domain &theDomain)
{
    myDomain = &theDomain;
}


// void addFE_Element(FE_Element *);
//	Method to add an element to the model.

bool
AnalysisModel::addFE_Element(FE_Element *theElement)
{

  // check we don't add a null pointer or this is a subclass
  // trying to use this method when it should'nt
  if (theElement == 0 || theFEs == 0)
      return false;
  
  // see if the current theFEs array is big enough, if not double
  // it's size by making a new one, copying the pointers and deleting old.
  if (numFE_Ele == sizeEle) { // we have to expand our componentArray array
      int newArraySize = 2 * sizeEle;
      FE_Element **newArray = new FE_Element *[newArraySize];
      if (newArray == 0) {
	  opserr << "AnalysisModel::addFE_Element -";
	  opserr << "could not allocate enough memory for new array\n";
	  exit(-1);
      }
      
      int i;
      for (i=0; i<newArraySize; i++)
	newArray[i] = 0;

      for (i=0; i<sizeEle; i++)
	newArray[i] = theFEs[i];

      delete  [] theFEs;
      theFEs = newArray;
      sizeEle = newArraySize;
  }

  // put the FE_Element into the list

  theFEs[numFE_Ele] = theElement;
  theElement->setAnalysisModel(*this);
  numFE_Ele++;

  return true;  // o.k.
}




// void addDOF_Group(DOF_Group *);
//	Method to add an element to the model.

bool
AnalysisModel::addDOF_Group(DOF_Group *theDOF_Group)
{

  // check we don't add a null pointer or this is a subclass trying
  // to use a method it should'nt be using
  if (theDOF_Group == 0 || theDOFs == 0)
      return false;
  
  // see if the current theFEs array is big enough, if not double
  // it's size by making a new one, copying the pointers and deleting old.

  if (numDOF_Grp == sizeDOF) { // we have to expand our componentArray array
      int newArraySize = 2 * sizeDOF;
      DOF_Group **newArray = new DOF_Group *[newArraySize];
      if (newArray == 0) {
	  opserr << "AnalysisModel::addDOF_Group -";
	  opserr << "could not allocate enough memory for new array\n";
	  exit(-1);
      }      
      int i;
      for (i=0; i<newArraySize; i++)
	newArray[i] = 0;

      for (i=0; i<sizeDOF; i++)
	newArray[i] = theDOFs[i];

      delete  [] theDOFs;
      theDOFs = newArray;
      sizeDOF = newArraySize;
  }

  // put the DOF_Group into the list

  theDOFs[numDOF_Grp] = theDOF_Group;
  numDOF_Grp++;
  return true;  // o.k.
}

void
AnalysisModel::clearAll(void) 
{
    // if the graphs have been constructed delete them
    if (myDOFGraph != 0)
	delete myDOFGraph;

    if (myGroupGraph != 0)
	delete myGroupGraph;    

    myDOFGraph = 0;
    myGroupGraph = 0;
    
    // remove the FE_Elements
    // they are not deleted as the ConstraintHandler may 
    // be able to reuse them
    for (int i=0; i<sizeEle; i++) {
	theFEs[i] = 0;
    }

    // remove the DOF_Groups 
    // they are not deleted as the ConstraintHandler may 
    // be able to reuse them.    
    for (int j=0; j<sizeDOF; j++) {
	theDOFs[j] =0;
    }

    numFE_Ele =0;
    numDOF_Grp = 0;
    numEqn = 0;    
}




int
AnalysisModel::getNumDOF_Groups(void) const
{
    return numDOF_Grp;
}


DOF_Group *
AnalysisModel::getDOF_GroupPtr(int tag)
{
    // check to see if in a simple position
    // this will be the case if DOF_Groups are created and added in order
    // which will probably be the typical situation (note true if DOF tags 
    // start from 0)
    if (tag < numDOF_Grp && tag > 0)
	if ((theDOFs[tag]->getTag()) == tag)
	    return theDOFs[tag];

    // else we have to search through and check till we find it
    for (int i=0; i<numDOF_Grp; i++)
	if ((theDOFs[i]->getTag()) == tag)
	    return theDOFs[i];
    
    return 0;
}


FE_EleIter &
AnalysisModel::getFEs()
{
    theFEiter.reset();
    return theFEiter;
}

DOF_GrpIter &
AnalysisModel::getDOFs()
{
    theDOFiter.reset();
    return theDOFiter;
}

void 
AnalysisModel::setNumEqn(int theNumEqn)
{
    numEqn = theNumEqn;
}

int 
AnalysisModel::getNumEqn(void) const
{
    return numEqn;
}


Graph &
AnalysisModel::getDOFGraph(void)
{
  if (myDOFGraph == 0) 
	myDOFGraph = new DOF_Graph(*this);
  return *myDOFGraph;
}


Graph &
AnalysisModel::getDOFGroupGraph(void)
{
    if (myGroupGraph == 0)
	myGroupGraph = new DOF_GroupGraph(*this);
    return *myGroupGraph;
}




void 
AnalysisModel::setResponse(const Vector &disp,
			   const Vector &vel, 
			   const Vector &accel)
{
    DOF_GrpIter &theDOFGrps = this->getDOFs();
    DOF_Group 	*dofPtr;
    
    while ((dofPtr = theDOFGrps()) != 0) {
	dofPtr->setNodeDisp(disp);
	dofPtr->setNodeVel(vel);
	dofPtr->setNodeAccel(accel);	
    }
}	
	
void 
AnalysisModel::setDisp(const Vector &disp)
{
    DOF_GrpIter &theDOFGrps = this->getDOFs();
    DOF_Group 	*dofPtr;

    while ((dofPtr = theDOFGrps()) != 0) 
	dofPtr->setNodeDisp(disp);
}	
	
void 
AnalysisModel::setVel(const Vector &vel)
{
        DOF_GrpIter &theDOFGrps = this->getDOFs();
    DOF_Group 	*dofPtr;
    
    while ((dofPtr = theDOFGrps()) != 0) 
	dofPtr->setNodeVel(vel);
}	
	

void 
AnalysisModel::setAccel(const Vector &accel)
{
    DOF_GrpIter &theDOFGrps = this->getDOFs();
    DOF_Group 	*dofPtr;
    
    while ((dofPtr = theDOFGrps()) != 0) 
	dofPtr->setNodeAccel(accel);	
}	

void 
AnalysisModel::incrDisp(const Vector &disp)
{
    DOF_GrpIter &theDOFGrps = this->getDOFs();
    DOF_Group 	*dofPtr;

    while ((dofPtr = theDOFGrps()) != 0) 
	dofPtr->incrNodeDisp(disp);
}	
	
void 
AnalysisModel::incrVel(const Vector &vel)
{
        DOF_GrpIter &theDOFGrps = this->getDOFs();
    DOF_Group 	*dofPtr;
    
    while ((dofPtr = theDOFGrps()) != 0) 
	dofPtr->incrNodeVel(vel);
}	
	
void 
AnalysisModel::incrAccel(const Vector &accel)
{
    DOF_GrpIter &theDOFGrps = this->getDOFs();
    DOF_Group 	*dofPtr;
    
    while ((dofPtr = theDOFGrps()) != 0) 
	dofPtr->incrNodeAccel(accel);	
}	


void 
AnalysisModel::setNumEigenvectors(int numEigenvectors)
{
    Node *theNode;
    NodeIter &theNodes = myDomain->getNodes();
    while ((theNode = theNodes()) != 0)
	theNode->setNumEigenvectors(numEigenvectors);
}	

void 
AnalysisModel::setEigenvalues(const Vector &eigenvalues)
{
    myDomain->setEigenvalues(eigenvalues);
}	

void 
AnalysisModel::setEigenvector(int mode, const Vector &eigenvalue)
{
    DOF_GrpIter &theDOFGrps = this->getDOFs();
    DOF_Group 	*dofPtr;
    
    while ((dofPtr = theDOFGrps()) != 0) 
	dofPtr->setEigenvector(mode, eigenvalue);	
}	

void 
AnalysisModel::applyLoadDomain(double pseudoTime)
{
    // check to see there is a Domain linked to the Model

    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::applyLoadDomain. No Domain linked.\n";
	return;
    }

    // invoke the method
    myDomain->applyLoad(pseudoTime);
}


int
AnalysisModel::updateDomain(void)
{
    // check to see there is a Domain linked to the Model

    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::updateDomain. No Domain linked.\n";
	return -1;
    }

    // invoke the method
    return myDomain->update();

}


int
AnalysisModel::updateDomain(double newTime, double dT)
{
    // check to see there is a Domain linked to the Model

    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::updateDomain. No Domain linked.\n";
	return -1;
    }

    // invoke the method
    return myDomain->update(newTime, dT);

}



int
AnalysisModel::commitDomain(void)
{
    // check to see there is a Domain linked to the Model
    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::commitDomain. No Domain linked.\n";
	return -1;
    }

    // invoke the method
    if (myDomain->commit() < 0) {
	opserr << "WARNING: AnalysisModel::commitDomain - Domain::commit() failed\n";
	return -2;
    }	

    return 0;
}

int
AnalysisModel::revertDomainToLastCommit(void)
{
    // check to see there is a Domain linked to the Model

    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::revertDomainToLastCommit.";
	opserr << " No Domain linked.\n";
	return -1;
    }

    // invoke the method
    if (myDomain->revertToLastCommit() < 0) {
	opserr << "WARNING: AnalysisModel::revertDomainToLastCommit.";
	opserr << " Domain::revertToLastCommit() failed.\n";
	return -2;
    }	
    return 0;
}

double
AnalysisModel::getCurrentDomainTime(void)
{
    // check to see there is a Domain linked to the Model

    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::getCurrentDomainTime.";
	opserr << " No Domain linked.\n";
	return 0.0;
    }

    // invoke the method
    return myDomain->getCurrentTime();
}


void
AnalysisModel::setCurrentDomainTime(double newTime)
{
    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::getCurrentDomainTime.";
	opserr << " No Domain linked.\n";
    }

    // invoke the method
    myDomain->setCurrentTime(newTime);
}



void
AnalysisModel::setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc)
{
    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::getCurrentDomainTime.";
	opserr << " No Domain linked.\n";
    }

    // invoke the method
    myDomain->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
}




Domain *
AnalysisModel::getDomainPtr(void) const
{
    return myDomain;
}


int
AnalysisModel::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}


int
AnalysisModel::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker) 
{
    return 0;
}

