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
                                                                        
// $Revision: 1.17 $
// $Date: 2009-08-26 20:33:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/subdomain/Subdomain.cpp,v $
                                                                        
// Written: fmk 
// Revision: A
// Revision: B 03/98 - revised to allow parallel model generation
//
// Description: This file contains the class definition for Subdomain.
// Subdomain is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints that have been added to the subdomain.
//
// What: "@(#) Subdomain.C, revA"

#include <Subdomain.h>
#include <stdlib.h>

#include <DomainComponent.h>
#include <Element.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <MapOfTaggedObjects.h>
#include <DomainDecompositionAnalysis.h>
#include <FE_Element.h>
#include <SingleDomNodIter.h>
#include <classTags.h>
//#include <PartitionedModelBuilder.h>
#include <DOF_Group.h>
#include <ElementIter.h>

#include <EquiSolnAlgo.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>

#include <FE_Element.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>


Matrix Subdomain::badResult(1,1); // for returns from getStiff, getMass and getDamp


Subdomain::Subdomain(int tag)
:Element(tag,ELE_TAG_Subdomain),
 Domain(),
 mapBuilt(false),map(0),mappedVect(0),mappedMatrix(0),
 realCost(0.0),cpuCost(0),pageCost(0),
 theAnalysis(0), extNodes(0), theFEele(0) 
{

  //thePartitionedModelBuilder = 0;
    // init the arrays.
    internalNodes = new MapOfTaggedObjects();
    externalNodes = new MapOfTaggedObjects();
    //    realExternalNodes = new MapOfTaggedObjects();    
    
    internalNodeIter = new SingleDomNodIter(internalNodes);
    externalNodeIter = new SingleDomNodIter(externalNodes);    
    theNodIter = new SubdomainNodIter(*this);

    // check that space was available
    if (internalNodes == 0 || externalNodes == 0 ||
	internalNodeIter == 0 || externalNodeIter == 0 ||
	theNodIter == 0) {
      
      opserr << "Subdomain::Subdomain() - ran out of memory\n";
      exit(-1);
    }
}


Subdomain::Subdomain(int tag, 
		     TaggedObjectStorage &theInternalNodeStorage,
		     TaggedObjectStorage &theExternalNodeStorage,
		     TaggedObjectStorage &theElementsStorage,
		     TaggedObjectStorage &theLoadPatternsStorage,	      
		     TaggedObjectStorage &theMPsStorage,
		     TaggedObjectStorage &theSPsStorage)
  :Element(tag,ELE_TAG_Subdomain), 
   Domain(theExternalNodeStorage, theElementsStorage,
	  theLoadPatternsStorage, 
	  theMPsStorage,theSPsStorage),
   mapBuilt(false),map(0),mappedVect(0),mappedMatrix(0),
   internalNodes(&theInternalNodeStorage),
   externalNodes(&theExternalNodeStorage), 
   realCost(0.0),cpuCost(0),pageCost(0),
   theAnalysis(0), extNodes(0), theFEele(0)
{
  //thePartitionedModelBuilder = 0;
  //    realExternalNodes = new MapOfTaggedObjects(256);    
    
    internalNodeIter = new SingleDomNodIter(internalNodes);
    externalNodeIter = new SingleDomNodIter(externalNodes);    

    // check that space was available
    if (internalNodes == 0 || externalNodes == 0 ||
	internalNodeIter == 0 || externalNodeIter == 0 ||
	theNodIter == 0) {
	
	opserr << "Subdomain::Subdomain() - ran out of memory\n";
	exit(-1);
    }    

}    


Subdomain::~Subdomain()
{
  if (internalNodes != 0)
    delete internalNodes;

  if (externalNodes != 0)
    delete externalNodes;
  
  if (internalNodeIter != 0)
    delete internalNodeIter;
  
  if (externalNodeIter != 0)
    delete externalNodeIter;

  if (map != 0)
    delete map;
  if (mappedVect != 0)
    delete mappedVect;
  if (mappedMatrix != 0)
    delete mappedMatrix;
}


void
Subdomain::clearAll(void) 
{
  this->Domain::clearAll();

  if (internalNodes != 0)
    internalNodes->clearAll();

  if (externalNodes != 0)
    externalNodes->clearAll();
}
/*
int 
Subdomain::buildSubdomain(int numSubdomains, PartitionedModelBuilder &theBuilder)
{
  int result = theBuilder.buildSubdomain(this->getTag(), numSubdomains, *this);
  if (result == 0) {
    this->hasDomainChanged();
    this->invokeChangeOnAnalysis();    
  }
  return result;
}
*/


// void addNode(Node *);
//	Method to add a Node to the model.

bool
Subdomain::addNode(Node * node)
{
#ifdef _G3DEBUG  
//  int nodTag = node->getTag();
//  // check no other node exists with same tag
//  Node *nodePtr = this->getNodePtr(nodTag);
//  if (nodePtr != 0)
//    return false;
//      
//      // MISSING CODE
#endif
  
  bool result = internalNodes->addComponent(node);
  if (result == true) {
      node->setDomain(this);
      this->domainChange();    
  }

  return result;
}

bool 
Subdomain::addExternalNode(Node *thePtr)
{
#ifdef _G3DEBUG
  // check to see it has not alredy been added
	
  int nodTag = thePtr->getTag();
  TaggedObject *other = externalNodes->getComponentPtr(nodTag);
  if (other != 0)
    return false;
  
  other = internalNodes->getComponentPtr(nodTag);
  if (other != 0)
    return false;
	
#endif
    // create a dummy Node & try adding it to the external nodes
    Node *newDummy = new Node(*thePtr, false);
    if (newDummy == 0)
	return false;

    bool result = externalNodes->addComponent(newDummy);
    if (result == true) {
      //	result = realExternalNodes->addComponent(thePtr);
	newDummy->setDomain(this);
	this->domainChange();    
    }
    
    return result;
}



Node *
Subdomain::removeNode(int tag)
{
  TaggedObject *object = internalNodes->removeComponent(tag);
  if (object == 0) {
      object = externalNodes->removeComponent(tag);      
      if (object != 0) {
        //	  Node *dummy = (Node *)object;
	//	  object = realExternalNodes->removeComponent(tag);      	  
	Node *result = (Node *)object;
	this->domainChange();          
	//	  delete dummy;
	return result;	  
      }
  }
  else {
      this->domainChange();          
      Node *result = (Node *)object;
      return result;	  
  }
  
  return 0;  
}

NodeIter &
Subdomain::getNodes()
{
    theNodIter->reset();
    return *theNodIter;
}

Node **
Subdomain::getNodePtrs(void)
{
  opserr << "Subdomain::getNodePtrs() - should not be called\n";
  return 0;
}


Node *
Subdomain::getNode(int tag) 
{

  TaggedObject *object = internalNodes->getComponentPtr(tag);
  if (object == 0) {
      object = externalNodes->getComponentPtr(tag);
      if (object != 0) {
	  Node *result = (Node *)object;
	  return result;
      }
  }
  else {
      Node *result = (Node *)object;
      return result;
  }

  return 0;  
}

bool
Subdomain::hasNode(int tag) 
{
  if (this->getNode(tag) != 0)
    return true;
  else
    return false;
}  

bool
Subdomain::hasElement(int tag) 
{
  if (this->getElement(tag) != 0)
    return true;
  else
    return false;
}  


int 
Subdomain::getNumNodes(void) const
{
    return internalNodes->getNumComponents() +
	externalNodes->getNumComponents();
}

int
Subdomain::commit(void) 
{
    this->Domain::commit();
    
    NodeIter &theNodes = this->getNodes();
    Node *nodePtr;
    while ((nodePtr = theNodes()) != 0)
	nodePtr->commitState();

    return 0;
}

int
Subdomain::revertToLastCommit(void) 
{
    this->Domain::revertToLastCommit();
    
    NodeIter &theNodes = this->getNodes();
    Node *nodePtr;
    while ((nodePtr = theNodes()) != 0)
	nodePtr->revertToLastCommit();

    return 0;
}

int
Subdomain::revertToStart(void) 
{
    this->Domain::revertToLastCommit();

    NodeIter &theNodes = this->getNodes();
    Node *nodePtr;
    while ((nodePtr = theNodes()) != 0)
	nodePtr->revertToStart();

    return 0;
}

int
Subdomain::update(void)
{
  return this->Domain::update();
}

int
Subdomain::update(double newTime, double dT)
{
    return this->Domain::update(newTime, dT);
}

void
Subdomain::Print(OPS_Stream &s, int flag)
{
  s << "Current Subdomain Information for Subdomain: ";
  s << this->getTag() << "\n";

  s << "\nINTERNAL NODE DATA: NumNodes: ";
  s << internalNodes->getNumComponents() << "\n"; 
  internalNodes->Print(s);

  s << "\nEXTERNAL NODE DATA: NumNodes: ";
  s << externalNodes->getNumComponents() << "\n"; 
  externalNodes->Print(s);

  this->Domain::Print(s);
  s << "\nEnd Subdomain Information\n";
}


void Subdomain::Print(OPS_Stream &s, ID *nodeTags, ID *eleTags, int flag)
{
  if (nodeTags != 0) {
    int numNodes = nodeTags->Size();
    for (int i=0; i<numNodes; i++) {
      int nodeTag = (*nodeTags)(i);
      TaggedObject *theNode = internalNodes->getComponentPtr(nodeTag);
      if (theNode != 0)
	theNode->Print(s, flag);
      else {
	TaggedObject *theNode = externalNodes->getComponentPtr(nodeTag);
	if (theNode != 0)
	  theNode->Print(s, flag);
      }
    }
  }

  /*
  if (eleTags != 0) {
    int numEles = eleTags->Size();
    for (int i=0; i<numEles; i++) {
      int eleTag = (*eleTags)(i);
      Element *theEle = this->getElement(eleTag);
      if (theEle != 0)
	theEle->Print(s, flag);
    }
  }
  */

  this->Domain::Print(s, 0, eleTags, flag);
}



NodeIter &
Subdomain::getInternalNodeIter(void)
{
    internalNodeIter->reset();
    return *internalNodeIter;
}


NodeIter &
Subdomain::getExternalNodeIter(void)
{
    externalNodeIter->reset();
    return *externalNodeIter;
}



void
Subdomain::wipeAnalysis(void)
{
  if (theAnalysis != 0) {
    theAnalysis->clearAll();
    delete theAnalysis;
  }
  theAnalysis = 0;
}

void
Subdomain::setDomainDecompAnalysis(DomainDecompositionAnalysis &theNewAnalysis)
{
    theAnalysis = &theNewAnalysis;
//    this->Domain::setAnalysis(theNewAnalysis);
}


int 
Subdomain::setAnalysisAlgorithm(EquiSolnAlgo &theAlgorithm)
{
  if (theAnalysis != 0)
    return theAnalysis->setAlgorithm(theAlgorithm);

  return 0;
}

int 
Subdomain::setAnalysisIntegrator(IncrementalIntegrator &theIntegrator)
{
  if (theAnalysis != 0)
    return theAnalysis->setIntegrator(theIntegrator);
  return 0;
}

int 
Subdomain::setAnalysisLinearSOE(LinearSOE &theSOE)
{
  if (theAnalysis != 0)
    return theAnalysis->setLinearSOE(theSOE);

  return 0;
}


int 
Subdomain::setAnalysisEigenSOE(EigenSOE &theSOE)
{
  if (theAnalysis != 0)
    return theAnalysis->setEigenSOE(theSOE);

  return 0;
}

int 
Subdomain::setAnalysisConvergenceTest(ConvergenceTest &theTest)
{
  if (theAnalysis != 0)
    return theAnalysis->setConvergenceTest(theTest);
  return 0;
}

int
Subdomain::invokeChangeOnAnalysis(void)
{
    int result = 0;
    if (theAnalysis != 0)
	result = theAnalysis->domainChanged();
    
    mapBuilt = false;
    return result;
}


int 
Subdomain::getNumExternalNodes(void) const    
{
    return externalNodes->getNumComponents();
}

const ID &
Subdomain::getExternalNodes()
{
    // first we check that extNodes exists and is of correct size
    int numExt = externalNodes->getNumComponents();
    if (extNodes == 0) {
	extNodes = new ID(numExt);
	if (extNodes == 0 || extNodes->Size() != numExt) {
	    opserr << "Subdomain::getExternalNodes(): ";
	    opserr << " - ran out of memory for size " << numExt <<endln;
	    exit(-1);
	}
    }
    
    if (extNodes->Size() != numExt) {
	delete extNodes;
	extNodes = new ID(numExt);
	if (extNodes == 0 || extNodes->Size() != numExt) {
	    opserr << "Subdomain::getExternalNodes(): ";
	    opserr << " - ran out of memory for size " << numExt <<endln;
	    exit(-1);
	}
    }

    // we now set the values of extNodes to be the node tags of the 
    // external nodes

    NodeIter &theExtNodes = this->getExternalNodeIter();
    Node *nodPtr;
    int cnt = 0;
    
    while ((nodPtr = theExtNodes()) != 0) 
	(*extNodes)(cnt++) = nodPtr->getTag();

    // done
    ID &res = *extNodes;
    return res;
}



int 
Subdomain::getNumDOF(void)
{
  if (theAnalysis != 0)
    return theAnalysis->getNumExternalEqn();
  else  {
    //   opserr << "Subdomain::getNumDOF() - no StaticAnalysis has been set\n";
    return 0;
  }
}
    
int
Subdomain::commitState(void)    
{
    return this->commit();
}

const Matrix &
Subdomain::getTangentStiff(void)
{
    opserr << "Subdomain::getTangentStiff(void)";
    opserr << "DOES NOT DO ANYTHING";
    return badResult;
}

const Matrix &
Subdomain::getInitialStiff(void)
{
    opserr << "Subdomain::getSecantStiff(void)";
    opserr << "DOES NOT DO ANYTHING";
    return badResult;
}

const Matrix &
Subdomain::getDamp(void)
{
    opserr << "Subdomain::getDamp(void)";
    opserr << "DOES NOT DO ANYTHING";    
    return badResult;
}

const Matrix &
Subdomain::getMass(void)
{
    opserr << "Subdomain::getMass(void)";
    opserr << "DOES NOT DO ANYTHING";    
    return badResult;
}




void  
Subdomain::zeroLoad(void)
{
    opserr << "Subdomain::zeroLoad() - should not be called\n";
}


int	  
Subdomain::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr << "Subdomain::addLoad() - should not be called\n";
    return 0;
}

int	  
Subdomain::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}


const Vector &
Subdomain::getResistingForce(void)    
{
    if (theAnalysis == 0) {
	opserr << "Subdomain::getResistingForce() ";
	opserr << " - no StaticCondensationAnalysis has been set\n";
	exit(-1);
    }
    
    if (mapBuilt == false)
	this->buildMap();
      
    ID &theMap = *map;
    const Vector &anaResidual = theAnalysis->getResidual();
    int numDOF = this->getNumDOF();
    for (int i=0; i<numDOF; i++)
	(*mappedVect)(i) = anaResidual(theMap(i));
    //opserr << "Subdomain::getResidual() : " << *mappedVect;
    return *mappedVect;
}


const Vector &
Subdomain::getResistingForceIncInertia(void)    
{
    opserr << "Subdomain::getResistingForceWithInertia() ";
    opserr << " - should not be called\n";

    return this->getResistingForce();
}



bool 
Subdomain::isSubdomain(void)
{
    return true;
}


int 
Subdomain::setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc)
{
  return this->Domain::setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
}

int 
Subdomain::computeTang(void)
{   
  if (theAnalysis != 0) {
    //    theTimer.start();
    
    int res =0;
    res = theAnalysis->formTangent();
    
    return res;
    
  } else {
    opserr << "Subdomain::getcomputeTang() ";
    opserr << " - no StaticCondensationAnalysis has been set\n";
    return 0;
  }
}



int 
Subdomain::computeResidual(void)
{
  if (theAnalysis != 0) {
    //    theTimer.start();
    
    int res =0;
    res = theAnalysis->formResidual();
    
    //theTimer.pause();
    //    realCost += theTimer.getReal();
    //    cpuCost += theTimer.getCPU();
    //    pageCost += theTimer.getNumPageFaults();
    
    return res;
    
    } else {
      opserr << "Subdomain::computeResidual() ";
      opserr << " - no StaticCondensationAnalysis has been set\n";
      return 0;
    }
}
    

const Matrix &
Subdomain::getTang(void)    
{
    if (theAnalysis == 0) {
	opserr << "Subdomain::getTang() ";
	opserr << " - no StaticCondensationAnalysis has been set\n";
	exit(-1);
    }	

    if (mapBuilt == false)
	this->buildMap();

    ID &theMap = *map;
    const Matrix &anaTang = theAnalysis->getTangent();
    int numDOF = this->getNumDOF();
    for (int i=0; i<numDOF; i++)
	for (int j=0; j<numDOF; j++)
	    (*mappedMatrix)(i,j) = anaTang(theMap(i),theMap(j));

    return *mappedMatrix;
}


void
Subdomain::setFE_ElementPtr(FE_Element *theFE_Ele)
{
    theFEele = theFE_Ele;
}


FE_Element *
Subdomain::getFE_ElementPtr(void)
{
    return theFEele;
}



const Vector &
Subdomain::getLastExternalSysResponse(void)
{
    if (theFEele == 0) {
	opserr << "FATAL ERROR: Subdomain::getLastExternalSysResponse() :";
	opserr << " - no FE_Element *exists for a subdomain\n";
	opserr << " This is the responsibilty of the FE_ELement constructor\n";
	exit(0);
    }

    // get the response from the FE_ele for the nodal
    // quantities - WARNING this is expressed in global dof

    if (mapBuilt == false)
      this->buildMap();

    ID &theMap = *map;
    const Vector &localResponse = theFEele->getLastResponse();
    int numDOF = this->getNumDOF();
    for (int i=0; i<numDOF; i++)
      (*mappedVect)(theMap(i)) = localResponse(i);

    return *mappedVect;
}
    

int 
Subdomain::computeNodalResponse(void)
{
    int res =0;
    if (theAnalysis != 0) {
	res = theAnalysis->computeInternalResponse();
	return res;
    }
    else {
	opserr << "Subdomain::computeNodalResponse() ";
	opserr << "- no StaticAnalysis has been set\n"; 
	return 0;
    }
}


int 
Subdomain::analysisStep(double dT)
{
  if (theAnalysis != 0)
    return theAnalysis->analysisStep(dT);

  return 0;
}


int
Subdomain::eigenAnalysis(int numMode, bool generalized, bool findSmallest)
{
  if (theAnalysis != 0)
    return theAnalysis->eigenAnalysis(numMode, generalized, findSmallest);
  return 0;
}


bool
Subdomain::doesIndependentAnalysis(void)
{
  if (theAnalysis != 0)
    return theAnalysis->doesIndependentAnalysis();
  else
    return true;
}


int 
Subdomain::sendSelf(int cTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    if (theAnalysis != 0) {
	ID data(2);
	data(0) = theAnalysis->getClassTag();
	data(1) = 0;
	theChannel.sendID(dataTag, cTag, data);
        
	return theAnalysis->sendSelf(cTag, theChannel);
    }
    else {
	opserr << "Subdomain::sendSelf - no analysis set\n";
    
    }
    return -1;
}

int 
Subdomain::recvSelf(int cTag, Channel &theChannel, 
		    FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    ID data(2);
    theChannel.recvID(dataTag, cTag, data);
    if (data(1) == 0) {
      theAnalysis = theBroker.getNewDomainDecompAnalysis(data(0),*this);
      if (theAnalysis != 0)
	return theAnalysis->recvSelf(cTag, theChannel,theBroker);
    }
    return -1;
}

double    
Subdomain::getCost(void) 
{
    double lastRealCost = realCost;

    realCost = 0.0;
    cpuCost = 0.0;
    pageCost = 0;

    return lastRealCost;
}


int
Subdomain::buildMap(void)
{
  if (mapBuilt == false) {
	// determine the mapping between local dof and subdomain ana dof
        int numDOF = this->getNumDOF();
	if (map == 0) 
	  map = new ID(numDOF);
	if (map->Size() != numDOF) {
	  delete map;
	  map = new ID(numDOF);
	}

	//	int numExt = theAnalysis->getNumExternalEqn();
	int numInt = theAnalysis->getNumInternalEqn();

	const ID &theExtNodes = this->getExternalNodes();
	int numExtNodes = theExtNodes.Size();
	int locInMap =0;
	for (int i=0; i<numExtNodes; i++) {
	  Node *nodePtr = this->getNode(theExtNodes(i));
	  int numNodeDOF = nodePtr->getNumberDOF();
	  DOF_Group *theDOF = nodePtr->getDOF_GroupPtr();
	  const ID &theLocalID = theDOF->getID();
	  for (int j=0; j<numNodeDOF; j++){
	    int locInSubdomainExt = theLocalID(j)-numInt;
	    (*map)(locInMap)=locInSubdomainExt;
	    locInMap++;
	  }
	}
	mapBuilt = true;

	if (mappedVect == 0) 
	  mappedVect = new Vector(numDOF);
	if (mappedVect->Size() != numDOF) {
	  delete mappedVect;
	  mappedVect = new Vector(numDOF);
	}

	if (mappedMatrix == 0) 
	  mappedMatrix = new Matrix(numDOF,numDOF);
	if (mappedMatrix->noRows() != numDOF) {
	  delete mappedMatrix;
	  mappedMatrix = new Matrix(numDOF,numDOF);
	}
  }
  
  return 0;
}


DomainDecompositionAnalysis *
Subdomain::getDDAnalysis(void)
{
  return theAnalysis;
}

int 
Subdomain::addResistingForceToNodalReaction(bool inclInertia)
{
  return 0;
}

int  
Subdomain::updateParameter(int tag, int value){
  return this->Domain::updateParameter(tag, value);
}

int  
Subdomain::updateParameter(int tag, double value)
{
  return this->Domain::updateParameter(tag, value);
}
