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
                                                                        
// $Revision: 1.16 $
// $Date: 2009-08-26 00:00:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/model/AnalysisModel.cpp,v $
                                                                        
                                                                        
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

#include <ArrayOfTaggedObjects.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <FE_EleIter.h>
#include <Graph.h>
#include <Vertex.h>
#include <Node.h>
#include <NodeIter.h>
#include <ConstraintHandler.h>


#include <MapOfTaggedObjects.h>

#define START_EQN_NUM 0
#define START_VERTEX_NUM 0

//  AnalysisModel();
//	constructor

AnalysisModel::AnalysisModel(int theClassTag)
:MovableObject(theClassTag),
 myDomain(0), myHandler(0),
 myDOFGraph(0), myGroupGraph(0),
 numFE_Ele(0), numDOF_Grp(0), numEqn(0)
{
    theFEs     = new ArrayOfTaggedObjects(1024);
    theDOFs    =  new ArrayOfTaggedObjects(1024);
    theFEiter  = new FE_EleIter(theFEs);
    theDOFiter = new DOF_GrpIter(theDOFs);

    // for subclasses to use - they provide own container stuff
} 

AnalysisModel::AnalysisModel()
:MovableObject(AnaMODEL_TAGS_AnalysisModel),
 myDomain(0), myHandler(0),
 myDOFGraph(0), myGroupGraph(0),
 numFE_Ele(0), numDOF_Grp(0), numEqn(0)
{
  theFEs     = new ArrayOfTaggedObjects(256);
  theDOFs    = new ArrayOfTaggedObjects(256);
  theFEiter  = new FE_EleIter(theFEs);
  theDOFiter = new DOF_GrpIter(theDOFs);
} 


AnalysisModel::AnalysisModel(TaggedObjectStorage &theFes, TaggedObjectStorage &theDofs)
:MovableObject(AnaMODEL_TAGS_AnalysisModel),
 myDomain(0), myHandler(0),
 myDOFGraph(0), myGroupGraph(0),
 numFE_Ele(0), numDOF_Grp(0), numEqn(0)
{
  theFEs     = &theFes;
  theDOFs    = &theDofs;
  theFEiter  = new FE_EleIter(theFEs);
  theDOFiter = new DOF_GrpIter(theDOFs);
} 


// ~AnalysisModel();    
AnalysisModel::~AnalysisModel()
{
  if (theFEs != 0) {
    theFEs->clearAll();
    delete theFEs;
  }

  if (theDOFs != 0) {
    theDOFs->clearAll();
    delete theDOFs;
  }

  if (theFEiter != 0)
    delete theFEiter;

  if (theDOFiter != 0)
    delete theDOFiter;

  if (myGroupGraph != 0) {
    delete myGroupGraph;    
  }	
  
  if (myDOFGraph != 0) {
    delete myDOFGraph;
  }
}    

void
AnalysisModel::setLinks(Domain &theDomain, ConstraintHandler &theHandler)
{
    myDomain = &theDomain;
    myHandler = &theHandler;
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

  // check if an Element with a similar tag already exists in the Domain
  int tag = theElement->getTag();
  TaggedObject *other = theFEs->getComponentPtr(tag);
  if (other != 0) {
    opserr << "AnalysisModel::addFE_Element - element with tag " << tag << "already exists in model\n"; 
    return false;
  }

  // add the element to the container object for the elements
  bool result = theFEs->addComponent(theElement);
  if (result == true) {
    theElement->setAnalysisModel(*this);
    numFE_Ele++;
    return true;  // o.k.
  } else
    return false;


  return result;
}




// void addDOF_Group(DOF_Group *);
//	Method to add an element to the model.

bool
AnalysisModel::addDOF_Group(DOF_Group *theGroup)
{

  // check we don't add a null pointer or this is a subclass trying
  // to use a method it should'nt be using
  if (theGroup == 0 || theDOFs == 0)
      return false;
  

  // check if an Element with a similar tag already exists in the Domain
  int tag = theGroup->getTag();
  TaggedObject *other = theDOFs->getComponentPtr(tag);
  if (other != 0) {
    opserr << "AnalysisModel::addDOF_Group - group with tag " << tag << "already exists in model\n"; 
    return false;
  }

  // add the element to the container object for the elements
  bool result = theDOFs->addComponent(theGroup);
  if (result == true) {
    numDOF_Grp++;
    return true;  // o.k.
  } else
    return false;
}

void
AnalysisModel::clearAll(void) 
{
    // if the graphs have been constructed delete them
    if (myDOFGraph != 0)
	delete myDOFGraph;

    if (myGroupGraph != 0)
	delete myGroupGraph;    

    theFEs->clearAll();
    theDOFs->clearAll();

    myDOFGraph = 0;
    myGroupGraph = 0;
    
    numFE_Ele =0;
    numDOF_Grp = 0;
    numEqn = 0;    
}

void
AnalysisModel::clearDOFGraph(void) 
{
  if (myDOFGraph != 0)
    delete myDOFGraph;

    myDOFGraph = 0;
}

void
AnalysisModel::clearDOFGroupGraph(void) 
{
  if (myGroupGraph != 0)
    delete myGroupGraph;    
  
  myGroupGraph = 0;
}




int
AnalysisModel::getNumDOF_Groups(void) const
{
  return numDOF_Grp;
}


DOF_Group *
AnalysisModel::getDOF_GroupPtr(int tag)
{
  TaggedObject *other = theDOFs->getComponentPtr(tag);
  if (other == 0) {
    return 0;
  }
  DOF_Group *result = (DOF_Group *)other;
  return result;
}


FE_EleIter &
AnalysisModel::getFEs()
{
    theFEiter->reset();
    return *theFEiter;
}

DOF_GrpIter &
AnalysisModel::getDOFs()
{
    theDOFiter->reset();
    return *theDOFiter;
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
  if (myDOFGraph == 0) {
    int numVertex = this->getNumDOF_Groups();

    //    myDOFGraph = new Graph(numVertex);
    MapOfTaggedObjects *graphStorage = new MapOfTaggedObjects();
    myDOFGraph = new Graph(*graphStorage);

    //
    // create a vertex for each dof
    //
    
    DOF_Group *dofPtr =0;
    DOF_GrpIter &theDOFs = this->getDOFs();
    while ((dofPtr = theDOFs()) != 0) {
      const ID &id = dofPtr->getID();
      int size = id.Size();
      for (int i=0; i<size; i++) {
	int dofTag = id(i);
	if (dofTag >= START_EQN_NUM) {
	  Vertex *vertexPtr = myDOFGraph->getVertexPtr(dofTag);
	  if (vertexPtr == 0) {
	    Vertex *vertexPtr = new Vertex(dofTag, dofTag);      
	    if (vertexPtr == 0) {
	      opserr << "WARNING AnalysisModel::getDOFGraph";
	      opserr << " - Not Enough Memory to create " << i+1 << "th Vertex\n";
	      return *myDOFGraph;
	    }
	    if (myDOFGraph->addVertex(vertexPtr, false) == false) {
	      opserr << "WARNING AnalysisModel::getDOFGraph - error adding vertex\n";
	      return *myDOFGraph;
	    }
	  }
	}
      }
    }
    
    // now add the edges, by looping over the FE_elements, getting their
    // IDs and adding edges between DOFs for equation numbers >= START_EQN_NUM
    
    FE_Element *elePtr =0;
    FE_EleIter &eleIter = this->getFEs();
    int cnt = 0;
    
    while((elePtr = eleIter()) != 0) {
      const ID &id = elePtr->getID();
      cnt++;
      int size = id.Size();
      for (int i=0; i<size; i++) {
	int eqn1 = id(i);
	
	// if eqnNum of DOF is a valid eqn number add an edge
	// to all other DOFs with valid eqn numbers.
	
	if (eqn1 >=START_EQN_NUM) {
	  for (int j=i+1; j<size; j++) {
	    int eqn2 = id(j);
	    if (eqn2 >=START_EQN_NUM)
	      myDOFGraph->addEdge(eqn1-START_EQN_NUM+START_VERTEX_NUM,
				  eqn2-START_EQN_NUM+START_VERTEX_NUM);
	  }
	}
      }
    }
  }    

  return *myDOFGraph;
}


Graph &
AnalysisModel::getDOFGroupGraph(void)
{
  if (myGroupGraph == 0) {
    int numVertex = this->getNumDOF_Groups();

    if (numVertex == 0) {
	opserr << "WARNING AnalysisMode::getGroupGraph";
	opserr << "  - 0 vertices, has the Domain been populated?\n";
	exit(-1);
    }	

    //    myGroupGraph = new Graph(numVertex);
    MapOfTaggedObjects *graphStorage = new MapOfTaggedObjects();
    myGroupGraph = new Graph(*graphStorage);

    if (numVertex == 0) {
	opserr << "WARNING AnalysisMode::getGroupGraph";
	opserr << "  - out of memory\n";
	exit(-1);
    }	
	
    DOF_Group *dofPtr;

    // now create the vertices with a reference equal to the DOF_Group number.
    // and a tag which ranges from 0 through numVertex-1

    DOF_GrpIter &dofIter2 = this->getDOFs();
    int count = START_VERTEX_NUM;
    while ((dofPtr = dofIter2()) != 0) {
	int DOF_GroupTag = dofPtr->getTag();
	int DOF_GroupNodeTag = dofPtr->getNodeTag();
	int numDOF = dofPtr->getNumFreeDOF();
	Vertex *vertexPtr = new Vertex(DOF_GroupTag, DOF_GroupNodeTag, 0, numDOF);

	if (vertexPtr == 0) {
	    opserr << "WARNING DOF_GroupGraph::DOF_GroupGraph";
	    opserr << " - Not Enough Memory to create ";
	    opserr << count << "th Vertex\n";
	    return *myGroupGraph;
	}
	
	myGroupGraph->addVertex(vertexPtr);
    }

    // now add the edges, by looping over the Elements, getting their
    // IDs and adding edges between DOFs for equation numbers >= START_EQN_NUM
    
    FE_Element *elePtr;
    FE_EleIter &eleIter = this->getFEs();

    while((elePtr = eleIter()) != 0) {
	const ID &id = elePtr->getDOFtags();
	int size = id.Size();
	for (int i=0; i<size; i++) {
	    int dof1 = id(i);
	    for (int j=0; j<size; j++) 
		if (i != j) {
		    int dof2 = id(j);
		    myGroupGraph->addEdge(dof1,dof2);
		}
	}
    }
  }

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

const Vector &
AnalysisModel::getEigenvalues(void)
{
  return myDomain->getEigenvalues();
}	



const Vector *
AnalysisModel::getModalDampingFactors(void){
  return myDomain->getModalDampingFactors();
}

bool 
AnalysisModel::inclModalDampingMatrix(void)
{
  return myDomain->inclModalDampingMatrix();
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
    myHandler->applyLoad();
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
    int res = myDomain->update();
    if (res == 0)
      return myHandler->update();

    return res;
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

    int res = 0;
    myDomain->applyLoad(newTime);
    if (res == 0)
      res = myHandler->applyLoad();
    if (res == 0)
      res = myDomain->update();
    if (res == 0)
      res = myHandler->update();

    return res;
}


int
AnalysisModel::analysisStep(double dT)
{
    // check to see there is a Domain linked to the Model

    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::newStep. No Domain linked.\n";
	return -1;
    }

    // invoke the method
    return myDomain->analysisStep(dT);
}

int
AnalysisModel::eigenAnalysis(int numMode, bool generalized, bool findSmallest)
{
    // check to see there is a Domain linked to the Model

    if (myDomain == 0) {
	opserr << "WARNING: AnalysisModel::newStep. No Domain linked.\n";
	return -1;
    }

    // invoke the method
    return myDomain->eigenAnalysis(numMode, generalized, findSmallest);
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

