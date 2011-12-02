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
                                                                        
// $Revision: 1.29 $
// $Date: 2005-03-30 20:23:15 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/Domain.cpp,v $
                                                                        
                                                                        
// FileP: ~/domain/domain/Domain.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Purpose: This file contains the class definition for Domain
// Domain is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints. These objects are all added to the Domain by a 
// ModelBuilder.
//
// What: "@(#) Domain.C, revA"

#include <stdlib.h>

#include <OPS_Globals.h>
#include <Domain.h>

#include <ElementIter.h>
#include <NodeIter.h>
#include <ElementalLoadIter.h>
#include <NodalLoadIter.h>
#include <Element.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>
#include <LoadPattern.h>

#include <ArrayOfTaggedObjects.h>
#include <ArrayOfTaggedObjectsIter.h>

#include <SingleDomEleIter.h>
#include <SingleDomNodIter.h>
#include <SingleDomSP_Iter.h>
#include <SingleDomMP_Iter.h>
#include <LoadPatternIter.h>
#include <SingleDomAllSP_Iter.h>

#include <Vertex.h>
#include <Graph.h>
#include <Recorder.h>
#include <MeshRegion.h>
#include <Analysis.h>
#include <FE_Datastore.h>
#include <FEM_ObjectBroker.h>

Domain::Domain()
:theRecorders(0), numRecorders(0),
 currentTime(0.0), committedTime(0.0), dT(0.0), currentGeoTag(0),
 hasDomainChangedFlag(false), theDbTag(0), lastGeoSendTag(-1),
 dbEle(0), dbNod(0), dbSPs(0), dbMPs(0), dbLPs(0),
 eleGraphBuiltFlag(false),  nodeGraphBuiltFlag(false), theNodeGraph(0), 
 theElementGraph(0), 
 theRegions(0), numRegions(0), commitTag(0),
 theBounds(6), theEigenvalues(0), theEigenvalueSetTime(0)
{
  
    // init the arrays for storing the domain components
    theElements = new ArrayOfTaggedObjects(4096);
    theNodes    = new ArrayOfTaggedObjects(4096);
    theSPs      = new ArrayOfTaggedObjects(256);
    theMPs      = new ArrayOfTaggedObjects(256);    
    theLoadPatterns = new ArrayOfTaggedObjects(32);

    // init the iters    
    theEleIter = new SingleDomEleIter(theElements);    
    theNodIter = new SingleDomNodIter(theNodes);
    theSP_Iter = new SingleDomSP_Iter(theSPs);
    theMP_Iter = new SingleDomMP_Iter(theMPs);
    theLoadPatternIter = new LoadPatternIter(theLoadPatterns);
    allSP_Iter = new SingleDomAllSP_Iter(*this);
    
    // check that there was space to create the data structures    
    if (theElements ==0 || theNodes == 0 || 
	theSPs == 0 || theMPs == 0 || 
	theEleIter == 0 || theNodIter == 0 || 
	theMP_Iter == 0 || theSP_Iter == 0 ||
	theLoadPatterns == 0 || theLoadPatternIter == 0) {	

      opserr << "Domain::Domain() - out of memory\n";
      exit(-1);
    }
    
    theBounds(0) = 0;
    theBounds(1) = 0;
    theBounds(2) = 0;
    theBounds(3) = 0;
    theBounds(4) = 0;    
    theBounds(5) = 0;            
}


Domain::Domain(int numNodes, int numElements, int numSPs, int numMPs,
	       int numLoadPatterns)
:theRecorders(0), numRecorders(0),
 currentTime(0.0), committedTime(0.0), dT(0.0), currentGeoTag(0),
 hasDomainChangedFlag(false), theDbTag(0), lastGeoSendTag(-1),
 dbEle(0), dbNod(0), dbSPs(0), dbMPs(0), dbLPs(0),
 eleGraphBuiltFlag(false), nodeGraphBuiltFlag(false), theNodeGraph(0), 
 theElementGraph(0),
 theRegions(0), numRegions(0), commitTag(0),
 theBounds(6), theEigenvalues(0), theEigenvalueSetTime(0)
{
    // init the arrays for storing the domain components
    theElements = new ArrayOfTaggedObjects(numElements);
    theNodes    = new ArrayOfTaggedObjects(numNodes);
    theSPs      = new ArrayOfTaggedObjects(numSPs);
    theMPs      = new ArrayOfTaggedObjects(numMPs);    
    theLoadPatterns = new ArrayOfTaggedObjects(numLoadPatterns);
    
    // init the iters
    theEleIter = new SingleDomEleIter(theElements);    
    theNodIter = new SingleDomNodIter(theNodes);
    theSP_Iter = new SingleDomSP_Iter(theSPs);
    theMP_Iter = new SingleDomMP_Iter(theMPs);
    theLoadPatternIter = new LoadPatternIter(theLoadPatterns);
    allSP_Iter = new SingleDomAllSP_Iter(*this);
    
    // check that there was space to create the data structures    
    if (theElements ==0 || theNodes == 0 || 
	theSPs == 0 || theMPs == 0 || 
	theEleIter == 0 || theNodIter == 0 || 
	theMP_Iter == 0 || theSP_Iter == 0 ||
	theLoadPatterns == 0 || theLoadPatternIter == 0) {	

	opserr << ("Domain::Domain(int, int, ...) - out of memory\n");
    }
    
    theBounds(0) = 0;
    theBounds(1) = 0;
    theBounds(2) = 0;
    theBounds(3) = 0;
    theBounds(4) = 0;    
    theBounds(5) = 0;            
}


Domain::Domain(TaggedObjectStorage &theNodesStorage,
	       TaggedObjectStorage &theElementsStorage,
	       TaggedObjectStorage &theMPsStorage,
	       TaggedObjectStorage &theSPsStorage,
	       TaggedObjectStorage &theLoadPatternsStorage)
:theRecorders(0), numRecorders(0),
 currentTime(0.0), committedTime(0.0), dT(0.0), currentGeoTag(0),
 hasDomainChangedFlag(false), theDbTag(0), lastGeoSendTag(-1),
 dbEle(0), dbNod(0), dbSPs(0), dbMPs(0), dbLPs(0),
 eleGraphBuiltFlag(false), nodeGraphBuiltFlag(false), theNodeGraph(0), 
 theElementGraph(0), 
 theElements(&theElementsStorage),
 theNodes(&theNodesStorage),
 theSPs(&theSPsStorage),
 theMPs(&theMPsStorage), 
 theLoadPatterns(&theLoadPatternsStorage),
 theRegions(0), numRegions(0), commitTag(0),
 theBounds(6), theEigenvalues(0), theEigenvalueSetTime(0)
{
    // init the iters    
    theEleIter = new SingleDomEleIter(theElements);    
    theNodIter = new SingleDomNodIter(theNodes);
    theSP_Iter = new SingleDomSP_Iter(theSPs);
    theMP_Iter = new SingleDomMP_Iter(theMPs);
    theLoadPatternIter = new LoadPatternIter(theLoadPatterns);
    allSP_Iter = new SingleDomAllSP_Iter(*this);

    // check that the containers are empty
    if (theElements->getNumComponents() != 0 ||
	theNodes->getNumComponents() != 0 ||
	theSPs->getNumComponents() != 0 ||
	theMPs->getNumComponents() != 0 ||
	theLoadPatterns->getNumComponents() != 0 ) {

	opserr << ("Domain::Domain(&, & ...) - out of memory\n");	
    }    	
	
    // check that there was space to create the data structures    
    if (theElements ==0 || theNodes == 0 || 
	theSPs == 0 || theMPs == 0 || 
	theEleIter == 0 || theNodIter == 0 ||
	theMP_Iter == 0 || theSP_Iter == 0 ||
	theLoadPatterns == 0 || theLoadPatternIter == 0) { 
    
	opserr << "FATAL Domain::Domain(TaggedObjectStorage, ...) - ";
	opserr << "Ran out of memory\n";
	exit(-1);
    }    
    
    theBounds(0) = 0;
    theBounds(1) = 0;
    theBounds(2) = 0;
    theBounds(3) = 0;
    theBounds(4) = 0;    
    theBounds(5) = 0;            
}    



Domain::Domain(TaggedObjectStorage &theStorage)
:theRecorders(0), numRecorders(0),
 currentTime(0.0), committedTime(0.0), dT(0.0), currentGeoTag(0),
 hasDomainChangedFlag(false), theDbTag(0), lastGeoSendTag(-1),
 dbEle(0), dbNod(0), dbSPs(0), dbMPs(0), dbLPs(0),
 eleGraphBuiltFlag(false), nodeGraphBuiltFlag(false), theNodeGraph(0), 
 theElementGraph(0), 
 theRegions(0), numRegions(0), commitTag(0),
 theBounds(6), theEigenvalues(0), theEigenvalueSetTime(0)
{
    // init the arrays for storing the domain components
    theStorage.clearAll(); // clear the storage just in case populated
    theElements = &theStorage;
    theNodes    = theStorage.getEmptyCopy();
    theSPs      = theStorage.getEmptyCopy();
    theMPs      = theStorage.getEmptyCopy();
    theLoadPatterns = theStorage.getEmptyCopy();    

    // init the iters    
    theEleIter = new SingleDomEleIter(theElements);    
    theNodIter = new SingleDomNodIter(theNodes);
    theSP_Iter = new SingleDomSP_Iter(theSPs);
    theMP_Iter = new SingleDomMP_Iter(theMPs);
    theLoadPatternIter = new LoadPatternIter(theLoadPatterns);
    allSP_Iter = new SingleDomAllSP_Iter(*this);

    // check that there was space to create the data structures    
    if (theElements ==0 || theNodes == 0 || 
	theSPs == 0 || theMPs == 0 || 
	theEleIter == 0 || theNodIter == 0 ||
	theMP_Iter == 0 || theSP_Iter == 0 ||
	theLoadPatterns == 0 || theLoadPatternIter == 0) { 
	
	opserr << ("Domain::Domain(ObjectStorage &) - out of memory\n");	
    }
    
    theBounds(0) = 0;
    theBounds(1) = 0;
    theBounds(2) = 0;
    theBounds(3) = 0;
    theBounds(4) = 0;    
    theBounds(5) = 0;            

    dbEle =0; dbNod =0; dbSPs =0; dbMPs =0; dbLPs = 0;
}





// ~Domain();    
//	destructor, this calls delete on all components of the model,
//	i.e. calls delete on all that is added to the model.
//	WARNING: if 3rd constructor, TaggedObjectStorage objects passed 
//      must have been created with new and nowhere else must the
//      destructor be called.

Domain::~Domain()
{

  // delete the objects in the domain
  this->clearAll();

    // delete all the storage objects
    // SEGMENT FAULT WILL OCCUR IF THESE OBJECTS WERE NOT CONSTRUCTED
    // USING NEW

    if (theElements != 0)
	delete theElements;    

    if (theNodes != 0)
	delete theNodes;

    if (theSPs != 0)
	delete theSPs;

    if (theMPs != 0)
	delete theMPs;

    if (theLoadPatterns != 0)
	delete theLoadPatterns;
    
    if (theEleIter != 0)
	delete theEleIter;
    
    if (theNodIter != 0)
	delete theNodIter;

    if (theSP_Iter != 0)
	delete theSP_Iter;
    
    if (theMP_Iter != 0)
	delete theMP_Iter;

    if (allSP_Iter != 0)
	delete allSP_Iter;

    if (theEigenvalues != 0)
      delete theEigenvalues;

	int i;
    for (i=0; i<numRecorders; i++)  
      delete theRecorders[i];
    
    if (theRecorders != 0) {
      delete [] theRecorders;
      theRecorders = 0;
    }

    for (i=0; i<numRegions; i++)  
      delete theRegions[i];
    
    if (theRegions != 0) {
      delete [] theRegions;
      theRegions = 0;
    }

    if (theNodeGraph != 0)
      delete theNodeGraph;

    if (theElementGraph != 0)
      delete theElementGraph;
  
    theRecorders = 0;
    numRecorders = 0;
}


// void addElement(Element *);
//	Method to add an element to the model.


bool
Domain::addElement(Element *element)
{
  int eleTag = element->getTag();

#ifdef _G3DEBUG
  // check all the elements nodes exist in the domain
  const ID &nodes = element->getExternalNodes();
  int numDOF = 0;
  for (int i=0; i<nodes.Size(); i++) {
      int nodeTag = nodes(i);
      Node *nodePtr = this->getNode(nodeTag);
      if (nodePtr == 0) {
	  opserr << "WARNING Domain::addElement - In element " << *element;
	  opserr << "\n no Node " << nodeTag << " exists in the domain\n";
	  return false;
      }
      numDOF += nodePtr->getNumberDOF();
  }   
#endif

  // check if an Element with a similar tag already exists in the Domain
  TaggedObject *other = theElements->getComponentPtr(eleTag);
  if (other != 0) {
    opserr << "Domain::addElement - element with tag " << eleTag << "already exists in model\n"; 
    return false;
  }

  // add the element to the container object for the elements
  bool result = theElements->addComponent(element);
  if (result == true) {
    element->setDomain(this);
    element->update();

    // finally check the ele has correct number of dof
#ifdef _G3DEBUG
    if (numDOF != element->getNumDOF()) { 
      
      opserr << "Domain::addElement - element " << eleTag << " - #DOF does not match with number at nodes\n";
      theElements->removeComponent(eleTag);
      return false;
    }
#endif      

    // mark the Domain as having been changed
    this->domainChange();
  } else 
    opserr << "Domain::addElement - element " << eleTag << "could not be added to container\n";      

  return result;
}



// void addNode(Node *);
//	Method to add a Node to the model.

bool
Domain::addNode(Node * node)
{
  int nodTag = node->getTag();

  TaggedObject *other = theNodes->getComponentPtr(nodTag);
  if (other != 0) {
    opserr << "Domain::addNode - node with tag " << nodTag << "already exists in model\n"; 
    return false;
  }
  
  bool result = theNodes->addComponent(node);
  if (result == true) {
      node->setDomain(this);
      this->domainChange();
      
      // see if the physical bounds are changed
      // note this assumes 0,0,0,0,0,0 as startup min,max values
      const Vector &crds = node->getCrds();
      int dim = crds.Size();
      if (dim >= 1) {
	  double x = crds(0);
	  if (x < theBounds(0)) theBounds(0) = x;
	  if (x > theBounds(3)) theBounds(3) = x;
      } 
      if (dim >= 2) {
	  double y = crds(1);
	  if (y < theBounds(1)) theBounds(1) = y;
	  if (y > theBounds(4)) theBounds(4) = y;	  
      } 
      if (dim == 3) {
	  double z = crds(2);
	  if (z < theBounds(2)) theBounds(2) = z;
	  if (z > theBounds(5)) theBounds(5) = z;	  
      }
      
  } else
    opserr << "Domain::addNode - node with tag " << nodTag << "could not be added to container\n";

  return result;
}


// void addSP_Constraint(SP_Constraint *);
//	Method to add a constraint to the model.
//

bool
Domain::addSP_Constraint(SP_Constraint *spConstraint)
{
#ifdef _G3DEBUG    
    // check the Node exists in the Domain
    int nodeTag = spConstraint->getNodeTag();
    Node *nodePtr = this->getNode(nodeTag);
    if (nodePtr == 0) {
      opserr << "Domain::addSP_Constraint - cannot add as node node with tag" <<
	nodeTag << "does not exist in model\n";       	
      return false;
    }

    // check that the DOF specified exists at the Node
    int numDOF = nodePtr->getNumberDOF();
    if (numDOF < spConstraint->getDOF_Number()) {
	opserr << "Domain::addSP_Constraint - cannot add as node node with tag" << 
	  nodeTag << "does not have associated constrained DOF\n"; 
	return false;
    }      
#endif

  // check that no other object with similar tag exists in model
  int tag = spConstraint->getTag();
  TaggedObject *other = theSPs->getComponentPtr(tag);
  if (other != 0) {
    opserr << "Domain::addSP_Constraint - cannot add as constraint with tag" << 
      tag << "already exists in model\n";             
    return false;
  }
  
  bool result = theSPs->addComponent(spConstraint);
  if (result == false) {
      opserr << "Domain::addSP_Constraint - cannot add constraint with tag" << 
	tag << "to the container\n";             
      return false;
  } 

  spConstraint->setDomain(this);
  this->domainChange();  

  return true;
}


// void addMP_Constraint(MP_Constraint *);
//	Method to add a constraint to the model.
//

bool
Domain::addMP_Constraint(MP_Constraint *mpConstraint)
{
#ifdef _G3DEBUG
    // perform the checks
    int nodeConstrained = mpConstraint->getNodeConstrained();
    Node *nodePtr = this->getNode(nodeConstrained);
    if (nodePtr == 0) {
      opserr << "Domain::addMP_Constraint -cannot add as constrained node with tag" <<
	nodeConstrained << "does not exist in model\n";       		
      return false;
    }
    
    int nodeRetained = mpConstraint->getNodeRetained();      
    nodePtr = this->getNode(nodeRetained);
    if (nodePtr == 0) {
      opserr << "Domain::addMP_Constraint - cannot add as retained node with tag" <<
	nodeRetained << "does not exist in model\n"; 	
      
      return false;
    }      
    // MISSING CODE
#endif

  // check that no other object with similar tag exists in model
  int tag = mpConstraint->getTag();
  TaggedObject *other = theMPs->getComponentPtr(tag);
  if (other != 0) {
    opserr << "Domain::addMP_Constraint - cannot add as constraint with tag" <<
      tag << "already exists in model";             
			      
    return false;
  }
  
  bool result = theMPs->addComponent(mpConstraint);
  if (result == true) {
      mpConstraint->setDomain(this);
      this->domainChange();
  } else
    opserr << "Domain::addMP_Constraint - cannot add constraint with tag" << 
      tag << "to the container\n";                   
			      
  return result;
}

bool 
Domain::addLoadPattern(LoadPattern *load)
{
    // first check if a load pattern with a similar tag exists in model
    int tag = load->getTag();
    TaggedObject *other = theLoadPatterns->getComponentPtr(tag);
    if (other != 0) {
      opserr << "Domain::addLoadPattern - cannot add as LoadPattern with tag" <<
	tag << "already exists in model\n";             
				
      return false;
    }    

    // now we add the load pattern to the container for load pattrens
    bool result = theLoadPatterns->addComponent(load);
    if (result == true) {
	load->setDomain(this);
	this->domainChange();
    }
    else 
      opserr << "Domain::addLoadPattern - cannot add LoadPattern with tag" <<
	tag << "to the container\n";                   	
			      
    return result;
}    


bool
Domain::addSP_Constraint(SP_Constraint *spConstraint, int pattern)
{
#ifdef _G3DEBUG
    // check the Node exists in the Domain
    int nodeTag = spConstraint->getNodeTag();
    Node *nodePtr = this->getNode(nodeTag);
    if (nodePtr == 0) {
      opserr << "Domain::addSP_Constraint - cannot add as node with tag" <<
	nodeTag << "does not exist in model\n";
	return false;
    }

    // check that the DOF specified exists at the Node
    int numDOF = nodePtr->getNumberDOF();
    if (numDOF < spConstraint->getDOF_Number()) {
      opserr << "Domain::addSP_Constraint - cannot add as node with tag" <<
	nodeTag << "does not have associated constrained DOF\n"; 

	return false;
    }      
#endif

  // now add it to the pattern
  TaggedObject *thePattern = theLoadPatterns->getComponentPtr(pattern);
  if (thePattern == 0) {
      opserr << "Domain::addSP_Constraint - cannot add as pattern with tag" <<
	pattern << "does not exist in domain\n"; 
			      
      return false;
  }
  LoadPattern *theLoadPattern = (LoadPattern *)thePattern;
  bool result = theLoadPattern->addSP_Constraint(spConstraint);
  if (result == false) {
    opserr << "Domain::addSP_Constraint - " << pattern << "pattern could not add the SP_Constraint\n"; 
			      
    return false;
  }

  spConstraint->setDomain(this);
  this->domainChange();  

  return true;
}

bool 
Domain::addNodalLoad(NodalLoad *load, int pattern)
{
    int nodTag = load->getNodeTag();
    Node *res = this->getNode(nodTag);
    if (res == 0) {
      opserr << "Domain::addNodalLoad() HI - no node with tag " << nodTag << 
	"exits in  the model, not adding the nodal load"  << *load << endln;
	return false;
    }

    // now add it to the pattern
    TaggedObject *thePattern = theLoadPatterns->getComponentPtr(pattern);
    if (thePattern == 0) {
      opserr << "Domain::addNodalLoad() - no pattern with tag" << 
	pattern << "in  the model, not adding the nodal load"  << *load << endln;
      
	return false;
    }
    LoadPattern *theLoadPattern = (LoadPattern *)thePattern;
    bool result = theLoadPattern->addNodalLoad(load);
    if (result == false) {
      opserr << "Domain::addNodalLoad() - pattern with tag" << 
	pattern << "could not add the load" << *load << endln;
				
      return false;
    }

    load->setDomain(this);    // done in LoadPattern::addNodalLoad()
    this->domainChange();

    return result;
}    


bool 
Domain::addElementalLoad(ElementalLoad *load, int pattern)
{
    // now add it to the pattern
    TaggedObject *thePattern = theLoadPatterns->getComponentPtr(pattern);
    if (thePattern == 0) {
      opserr << "Domain::addNodalLoad() - no pattern with tag " << pattern << 
	"exits in  the model, not adding the ele load " << *load << endln;

	return false;
    }
    LoadPattern *theLoadPattern = (LoadPattern *)thePattern;
    bool result = theLoadPattern->addElementalLoad(load);
    if (result == false) {
      opserr << "Domain::addNodalLoad() - no pattern with tag" << 
	pattern << "in  the model, not adding the ele load" << *load << endln;
      return false;
    }


    // load->setDomain(this); // done in LoadPattern::addElementalLoad()
    this->domainChange();
    return result;
}


/* GENERAL NOTE ON REMOVAL OF COMPONENTS:
**   downward casts (while bad) are o.k. as only the type
**   of components can be added to the storage objects, e.g.
**   only elements can be added to theElements therefore
**   casting a DomainComponent from theElements to an Element is o.k.
*/

void
Domain::clearAll(void) {
    // clear the loads and constraints from any load pattern
    LoadPatternIter &thePatterns = this->getLoadPatterns();
    LoadPattern *thePattern;
    while ((thePattern = thePatterns()) != 0)
	thePattern->clearAll();

    // clean out the containers
    theElements->clearAll();
    theNodes->clearAll();
    theSPs->clearAll();
    theMPs->clearAll();
    theLoadPatterns->clearAll();

    // remove the recorders
    int i;
    for (i=0; i<numRecorders; i++)
	delete theRecorders[i];
    numRecorders = 0; 
    
    if (theRecorders != 0) {
      delete [] theRecorders;
      theRecorders = 0;
    }

    for (i=0; i<numRegions; i++)
	delete theRegions[i];
    numRegions = 0;
    
    if (theRegions != 0) {
      delete [] theRegions;
      theRegions = 0;
    }

    // set the time back to 0.0
    currentTime = 0.0;
    committedTime = 0.0;
    dT = 0.0;

    // set the bounds around the origin
    theBounds(0) = 0;
    theBounds(1) = 0;
    theBounds(2) = 0;
    theBounds(3) = 0;
    theBounds(4) = 0;    
    theBounds(5) = 0;        

    currentGeoTag = 0;
    lastGeoSendTag = -1;

    // rest the flag to be as initial
    hasDomainChangedFlag = false;
    nodeGraphBuiltFlag = false;
    eleGraphBuiltFlag = false;

    dbEle =0; dbNod =0; dbSPs =0; dbMPs =0; dbLPs = 0;
}


Element *
Domain::removeElement(int tag)
{
  // remove the object from the container    
  TaggedObject *mc = theElements->removeComponent(tag);
  
  // if not there return 0
  if (mc == 0) 
      return 0;

  // otherwise mark the domain as having changed
  this->domainChange();
  
  // perform a downward cast to an Element (safe as only Element added to
  // this container, 0 the Elements DomainPtr and return the result of the cast  
  Element *result = (Element *)mc;
  //  result->setDomain(0);
  return result;
}

Node *
Domain::removeNode(int tag)
{
  // remove the object from the container
  TaggedObject *mc = theNodes->removeComponent(tag);
  
  // if not there return 0
  if (mc == 0) 
      return 0;  

  // mark the domain has having changed 
  this->domainChange();
  
  // perform a downward cast to a Node (safe as only Node added to
  // this container and return the result of the cast
  Node *result = (Node *)mc;
  // result->setDomain(0);
  return result;
}


SP_Constraint *
Domain::removeSP_Constraint(int tag)
{
    // remove the object from the container    
    TaggedObject *mc = theSPs->removeComponent(tag);
    
    // if not there return 0    
    if (mc == 0) 
	return 0;

    // mark the domain as having changed    
    this->domainChange();
    
    // perform a downward cast, set the objects domain pointer to 0
    // and return the result of the cast    
    SP_Constraint *result = (SP_Constraint *)mc;
    // result->setDomain(0);

    // should check that theLoad and result are the same    
    return result;
}

MP_Constraint *
Domain::removeMP_Constraint(int tag)
{
    // remove the object from the container        
    TaggedObject *mc = theMPs->removeComponent(tag);
    
    // if not there return 0    
    if (mc == 0) 
	return 0;

    // mark the domain as having changed    
    this->domainChange();
    
    // perform a downward cast, set the objects domain pointer to 0
    // and return the result of the cast        
    MP_Constraint *result = (MP_Constraint *)mc;
    // result->setDomain(0);
    return result;
}    

LoadPattern *
Domain::removeLoadPattern(int tag)
{
    // remove the object from the container            
    TaggedObject *obj = theLoadPatterns->removeComponent(tag);
    
    // if not there return 0    
    if (obj == 0)
	return 0;
    
    // perform a downward cast, set the objects domain pointer to 0
    // and return the result of the cast            
    LoadPattern *result = (LoadPattern *)obj;
    // result->setDomain(0);

    //
    // now set the Domain pointer for all loads and SP constraints 
    // in the loadPattern to be 0
    //
    
    NodalLoad *theNodalLoad;
    NodalLoadIter &theNodalLoads = result->getNodalLoads();
    while ((theNodalLoad = theNodalLoads()) != 0) {
      // theNodalLoad->setDomain(0);
    }

    ElementalLoad *theElementalLoad;
    ElementalLoadIter &theElementalLoads = result->getElementalLoads();
    while ((theElementalLoad = theElementalLoads()) != 0) {
      // theElementalLoad->setDomain(0);
    }

    int numSPs = 0;
    SP_Constraint *theSP_Constraint;
    SP_ConstraintIter &theSP_Constraints = result->getSPs();
    while ((theSP_Constraint = theSP_Constraints()) != 0) {
	numSPs++;
	// theSP_Constraint->setDomain(0);
    }

    // mark the domain has having changed if numSPs > 0
    // as the constraint handlers have to be redone
    if (numSPs > 0)
	this->domainChange();

    // finally return the load pattern
    return result;    
}    





NodalLoad *
Domain::removeNodalLoad(int tag, int loadPattern)
{
  // remove the object from the container            
  LoadPattern *theLoadPattern = this->getLoadPattern(loadPattern);
    
  // if not there return 0    
  if (theLoadPattern == 0)
    return 0;
    
  return theLoadPattern->removeNodalLoad(tag);
}    


ElementalLoad *
Domain::removeElementalLoad(int tag, int loadPattern)
{
  // remove the object from the container            
  LoadPattern *theLoadPattern = this->getLoadPattern(loadPattern);
    
  // if not there return 0    
  if (theLoadPattern == 0)
    return 0;
    
  return theLoadPattern->removeElementalLoad(tag);
}    


SP_Constraint *
Domain::removeSP_Constraint(int tag, int loadPattern)
{
  // remove the object from the container            
  LoadPattern *theLoadPattern = this->getLoadPattern(loadPattern);
    
  // if not there return 0    
  if (theLoadPattern == 0)
    return 0;
    
  SP_Constraint *theSP = theLoadPattern->removeSP_Constraint(tag);
  if (theSP != 0)
    this->domainChange();

  return theSP;
}    

ElementIter &
Domain::getElements()
{
    theEleIter->reset();    
    return *theEleIter;
}


NodeIter &
Domain::getNodes()
{
    theNodIter->reset();    
    return *theNodIter;
}

SP_ConstraintIter &
Domain::getSPs()
{
    theSP_Iter->reset();
    return *theSP_Iter;;
}

SP_ConstraintIter &
Domain::getDomainAndLoadPatternSPs()
{
    allSP_Iter->reset();
    return *allSP_Iter;;
}


MP_ConstraintIter &
Domain::getMPs()
{
    theMP_Iter->reset();
    return *theMP_Iter;;
}


LoadPatternIter &
Domain::getLoadPatterns()
{
    theLoadPatternIter->reset();
    return *theLoadPatternIter;;
}

/* GENERAL NOTE ON RETRIEVAL OF COMPONENT PTRs:
**   downward casts (while bad) are o.k. as only the type
**   of components can be added to the storage objects, e.g.
**   only elements can be added to theElements
*/

Element *
Domain::getElement(int tag) 
{
  TaggedObject *mc = theElements->getComponentPtr(tag);
  
  // if not there return 0 otherwise perform a cast and return that
  if (mc == 0) 
      return 0;
  Element *result = (Element *)mc;
  return result;
}


Node *
Domain::getNode(int tag) 
{
  TaggedObject *mc = theNodes->getComponentPtr(tag);

  // if not there return 0 otherwise perform a cast and return that  
  if (mc == 0) 
      return 0;  
  Node *result = (Node *)mc;
  return result;
}

SP_Constraint *
Domain::getSP_Constraint(int tag) 
{
  TaggedObject *mc = theSPs->getComponentPtr(tag);

  // if not there return 0 otherwise perform a cast and return that  
  if (mc == 0) 
      return 0;
  SP_Constraint *result = (SP_Constraint *)mc;
  return result;
}

MP_Constraint *
Domain::getMP_Constraint(int tag) 
{
  TaggedObject *mc = theMPs->getComponentPtr(tag);

  // if not there return 0 otherwise perform a cast and return that
  if (mc == 0) 
      return 0;
  MP_Constraint *result = (MP_Constraint *)mc;
  return result;
}

LoadPattern *
Domain::getLoadPattern(int tag) 
{
  TaggedObject *mc = theLoadPatterns->getComponentPtr(tag);
  // if not there return 0 otherwise perform a cast and return that  
  if (mc == 0) 
      return 0;
  LoadPattern *result = (LoadPattern *)mc;
  return result;
}


double
Domain::getCurrentTime(void) const
{
    return currentTime;
}

int
Domain::getCommitTag(void) const
{
  return commitTag;
}


int 
Domain::getNumElements(void) const
{
    return theElements->getNumComponents();
}
int 
Domain::getNumNodes(void) const
{
    return theNodes->getNumComponents();
}

int 
Domain::getNumSPs(void) const
{
    return theSPs->getNumComponents();
}


int 
Domain::getNumMPs(void) const
{
    return theMPs->getNumComponents();
}

int 
Domain::getNumLoadPatterns(void) const
{
    return theLoadPatterns->getNumComponents();
}


const Vector &
Domain::getPhysicalBounds(void)
{
    return theBounds;
}



Graph  &
Domain::getElementGraph(void)
{
    if (eleGraphBuiltFlag == false) {
	// if the current graph is out of date .. delete it so we can start again
	if (theElementGraph != 0) {
	    delete theElementGraph;
	    theElementGraph = 0;
	}	
	// create an empty graph 
        theElementGraph = new Graph(this->getNumElements()+START_VERTEX_NUM);

	if (theElementGraph == 0) {// if still 0 try a smaller one
	    theElementGraph = new Graph();
	    
	    if (theElementGraph == 0) { // if still 0 out of memory
		opserr << "Domain::getElementGraph() - out of memory\n";
		exit(-1);
	    }
	}

	// now build the graph
	if (this->buildEleGraph(theElementGraph) == 0)
	    eleGraphBuiltFlag = true;
	else
	    opserr << "Domain::getElementGraph() - failed to build the element graph\n";	    
    }
    
    // return the Graph
    return *theElementGraph;
}



Graph  &
Domain::getNodeGraph(void)
{
    if (nodeGraphBuiltFlag == false) {
	
	// if the current graph is out of date .. delete it so we can start again
	if (theNodeGraph != 0) {
	    delete theNodeGraph;
	    theNodeGraph = 0;
	}

	// try to get a graph as big as we should need
	theNodeGraph = new Graph(this->getNumNodes()+START_VERTEX_NUM);
	
	if (theNodeGraph == 0) { // if still 0 try a smaller one
	    theNodeGraph = new Graph();

	    if (theNodeGraph == 0) {// if still 0 out of memory
		opserr << "Domain::getNodeGraph() - out of memory\n";
		exit(-1);
	    }
	}

       // now build the graph
	if (this->buildNodeGraph(theNodeGraph) == 0)
	    nodeGraphBuiltFlag = true;
	else
	    opserr << "Domain::getNodeGraph() - failed to build the node graph\n";
    }

    // return the Graph
    return *theNodeGraph;
}


void
Domain::setCommitTag(int newTag)
{
    commitTag = newTag;
}

void
Domain::setCurrentTime(double newTime)
{
    currentTime = newTime;
    dT = currentTime - committedTime;
}

void
Domain::setCommittedTime(double newTime)
{
    committedTime = newTime;
    dT = currentTime - committedTime;
}


void
Domain::applyLoad(double timeStep)
{
    // set the current pseudo time in the domai to be newTime
    currentTime = timeStep;
    dT = currentTime - committedTime;

    //
    // first loop over nodes and elements getting them to first zero their loads
    //

    Node *nodePtr;
    NodeIter &theNodeIter = this->getNodes();
    while ((nodePtr = theNodeIter()) != 0)
	nodePtr->zeroUnbalancedLoad();

    Element *elePtr;
    ElementIter &theElemIter = this->getElements();    
    while ((elePtr = theElemIter()) != 0)
	if (elePtr->isSubdomain() == false)
	    elePtr->zeroLoad();    

    // now loop over load patterns, invoking applyLoad on them
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = this->getLoadPatterns();
    while((thePattern = thePatterns()) != 0)
      thePattern->applyLoad(timeStep);

    //
    // finally loop over the MP_Constraints and SP_Constraints of the domain
    //
    
    MP_ConstraintIter &theMPs = this->getMPs();
    MP_Constraint *theMP;
    while ((theMP = theMPs()) != 0)
	theMP->applyConstraint(timeStep);

    SP_ConstraintIter &theSPs = this->getSPs();
    SP_Constraint *theSP;
    while ((theSP = theSPs()) != 0)
	theSP->applyConstraint(timeStep);
}


void
Domain::setLoadConstant(void)
{
    // loop over all the load patterns that are currently added to the domain
    // getting them to set their loads as now constant
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = this->getLoadPatterns();
    while((thePattern = thePatterns()) != 0)
      thePattern->setLoadConstant();
}



int
Domain::initialize(void)
{
  Element *elePtr;
  ElementIter &theElemIter = this->getElements();    
  while ((elePtr = theElemIter()) != 0) 
    // lvalue needed here for M$ VC++ compiler -- MHS
    const Matrix &ret = elePtr->getInitialStiff();

  return 0;
}


int
Domain::setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc)
{
  int result = 0;
  Element *elePtr;
  ElementIter &theElemIter = this->getElements();    
  while ((elePtr = theElemIter()) != 0) 
    result += elePtr->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);

  Node *nodePtr;
  NodeIter &theNodeIter = this->getNodes();
  while ((nodePtr = theNodeIter()) != 0) {
    result += nodePtr->setRayleighDampingFactor(alphaM);
  }

  return result;
}



int
Domain::commit(void)
{
    // 
    // first invoke commit on all nodes and elements in the domain
    //
    Node *nodePtr;
    NodeIter &theNodeIter = this->getNodes();
    while ((nodePtr = theNodeIter()) != 0) {
      nodePtr->commitState();
    }

    Element *elePtr;
    ElementIter &theElemIter = this->getElements();    
    while ((elePtr = theElemIter()) != 0) {
      elePtr->commitState();
    }

    // set the new committed time in the domain
    committedTime = currentTime;
    dT = 0.0;

    // invoke record on all recorders
    for (int i=0; i<numRecorders; i++)
	theRecorders[i]->record(commitTag, currentTime);

    // update the commitTag
    commitTag++;
    return 0;
}

int
Domain::revertToLastCommit(void)
{
    // 
    // first invoke revertToLastCommit  on all nodes and elements in the domain
    //
    
    Node *nodePtr;
    NodeIter &theNodeIter = this->getNodes();
    while ((nodePtr = theNodeIter()) != 0)
	nodePtr->revertToLastCommit();
    
    Element *elePtr;
    ElementIter &theElemIter = this->getElements();    
    while ((elePtr = theElemIter()) != 0) {
	elePtr->revertToLastCommit();
    }

    // set the current time and load factor in the domain to last committed
    currentTime = committedTime;
    dT = 0.0;

    // apply load for the last committed time
    this->applyLoad(currentTime);

    return this->update();
}

int
Domain::revertToStart(void)
{
    // 
    // first invoke revertToLastCommit  on all nodes and 
    // elements in the domain
    //
    
    Node *nodePtr;
    NodeIter &theNodeIter = this->getNodes();
    while ((nodePtr = theNodeIter()) != 0) 
	nodePtr->revertToStart();

    Element *elePtr;
    ElementIter &theElements = this->getElements();    
    while ((elePtr = theElements()) != 0) {
	elePtr->revertToStart();
    }

// ADDED BY TERJE //////////////////////////////////
    // invoke 'restart' on all recorders
    for (int i=0; i<numRecorders; i++) 
		theRecorders[i]->restart();
/////////////////////////////////////////////////////

    // set the current time and load factor in the domain to last committed
    committedTime = 0;
    currentTime = 0;
    dT = 0.0;

    // apply load for the last committed time
    this->applyLoad(currentTime);

    return this->update();
}

int
Domain::update(void)
{
  // set the global constants
  ops_Dt = dT;
  ops_TheActiveDomain = this;

  int ok = 0;

  // invoke update on all the ele's
  ElementIter &theEles = this->getElements();
  Element *theEle;

  while ((theEle = theEles()) != 0) {
    ok += theEle->update();
  }

  if (ok != 0)
    opserr << "Domain::update - domain failed in update\n";

  return ok;
}


int
Domain::update(double newTime, double dT)
{
  this->applyLoad(newTime);
  this->update();
  return 0;
}



int
Domain::newStep(double dT)
{
  return 0;
}



int
Domain::setEigenvalues(const Vector &theValues)
{
  // make sure the eigen value vector is large enough
  if (theEigenvalues == 0 || theEigenvalues->Size() != theValues.Size()) {
    
    // if not zero delete the old and create a new one
    if (theEigenvalues != 0)
      delete theEigenvalues;

    // create the new vector
    theEigenvalues = new Vector(theValues);
  } else

    // otherwise just a straight assignment
    *theEigenvalues = theValues;


  // now set the time at which eigen values were determined to be current domain time
  theEigenvalueSetTime = this->getCurrentTime();

  return 0;
}


const Vector &
Domain::getEigenvalues(void) 
{
  // ensure the eigen values were set
  if (theEigenvalues == 0) {
    opserr << "Domain::getEigenvalues - Eigenvalues were never set\n";
    exit(-1);
  }

  return *theEigenvalues;
}  

double 
Domain::getTimeEigenvaluesSet(void) 
{
  return theEigenvalueSetTime;
}

void
Domain::setDomainChangeStamp(int newStamp)
{
    currentGeoTag = newStamp;
}


void
Domain::domainChange(void)
{
    hasDomainChangedFlag = true;
}


int
Domain::hasDomainChanged(void)
{	
    // if the flag indicating the domain has changed since the
    // last call to this method has changed, increment the integer
    // and reset the flag
    bool result = hasDomainChangedFlag;
    hasDomainChangedFlag = false;
    if (result == true) {
	currentGeoTag++;
	nodeGraphBuiltFlag = false;
	eleGraphBuiltFlag = false;
    }

    // return the integer so user can determine if domain has changed 
    // since their last call to this method
    return currentGeoTag;
}

void
Domain::Print(OPS_Stream &s, int flag) 
{

  s << "Current Domain Information\n";
  s << "\tCurrent Time: " << currentTime;
  s << "\ntCommitted Time: " << committedTime << endln;

  s << "\nNODE DATA: NumNodes: " << theNodes->getNumComponents() << "\n";
  theNodes->Print(s, flag);

  s << "\nELEMENT DATA: NumEle: " << theElements->getNumComponents() << "\n";
  theElements->Print(s, flag);
  
  s << "\nSP_Constraints: numConstraints: ";
  s << theSPs->getNumComponents() << "\n";
  theSPs->Print(s, flag);
  
  s << "\nMP_Constraints: numConstraints: ";
  s << theMPs->getNumComponents() << "\n";
  theMPs->Print(s, flag);

  s << "\nLOAD PATTERNS: num Patterns: ";
  s << theLoadPatterns->getNumComponents() << "\n\n";
  theLoadPatterns->Print(s, flag);
}

OPS_Stream &operator<<(OPS_Stream &s, Domain &M)
{
  M.Print(s);
  return s;
}


int
Domain::addRecorder(Recorder &theRecorder)
{
  if (theRecorder.setDomain(*this) != 0) {
    opserr << "Domain::addRecorder() - recorder could not be added\n";
    return -1;
  }

  Recorder **newRecorders = new Recorder *[numRecorders + 1]; 
  if (newRecorders == 0) {
    opserr << "Domain::addRecorder() - could not add ran out of memory\n";
    return -1;
  }
  
  for (int i=0; i<numRecorders; i++)
    newRecorders[i] = theRecorders[i];
  newRecorders[numRecorders] = &theRecorder;
  
  if (theRecorders != 0)
    delete [] theRecorders;
  
  theRecorders = newRecorders;
  numRecorders++;
  return 0;
}


int
Domain::removeRecorders(void)
{
    for (int i=0; i<numRecorders; i++)  
      delete theRecorders[i];
    
    if (theRecorders != 0) {
      delete [] theRecorders;
    }
  
    theRecorders = 0;
    numRecorders = 0;
    return 0;
}

int  
Domain::addRegion(MeshRegion &theRegion)
{
    MeshRegion **newRegions = new MeshRegion *[numRegions + 1]; 
    if (newRegions == 0) {
	opserr << "Domain::addRegion() - could not add ran out of memory\n";
	return -1;
    }
    
    for (int i=0; i<numRegions; i++)
	newRegions[i] = theRegions[i];
    newRegions[numRegions] = &theRegion;
    theRegion.setDomain(this);
    if (theRegions != 0)
      delete [] theRegions;
    
    theRegions = newRegions;
    numRegions++;
    return 0;
}

MeshRegion *
Domain::getRegion(int tag)
{
    for (int i=0; i<numRegions; i++)
      if (theRegions[i]->getTag() == tag)
	return theRegions[i];

    return 0;
}


int 
Domain::buildEleGraph(Graph *theEleGraph)
{
    int numVertex = this->getNumElements();

    // see if quick return

    if (numVertex == 0) 
	return 0;

    
    // create another vertices array which aids in adding edges
    
    int *theElementTagVertices = 0;
    int maxEleNum = 0;
    Element *elePtr;
    ElementIter &eleIter = this->getElements();
    while ((elePtr = eleIter()) != 0)
	if (elePtr->getTag() > maxEleNum)
	    maxEleNum = elePtr->getTag();

    theElementTagVertices = new int[maxEleNum+1];

    if (theElementTagVertices == 0) {
	opserr << "WARNING Domain::buildEleGraph ";
	opserr << " - Not Enough Memory for ElementTagVertices\n";
	return -1;
    }

    for (int j=0; j<=maxEleNum; j++) theElementTagVertices[j] = -1;
    opserr << "Domain::buildEleGraph numVertex maxEleNum " << numVertex << " " << maxEleNum << endln;
    // now create the vertices with a reference equal to the element number.
    // and a tag which ranges from 0 through numVertex-1

    ElementIter &eleIter2 = this->getElements();
    int count = START_VERTEX_NUM;
    while ((elePtr = eleIter2()) != 0) {
	int ElementTag = elePtr->getTag();
	Vertex *vertexPtr = new Vertex(count,ElementTag);

	if (vertexPtr == 0) {
	    opserr << "WARNING Domain::buildEleGraph";
	    opserr << " - Not Enough Memory to create ";
	    opserr << count << "th Vertex\n";
	    delete [] theElementTagVertices;
	    return -1;
	}

	theEleGraph->addVertex(vertexPtr);
	theElementTagVertices[ElementTag] = count++;
	
    }

    // We now need to determine which elements are asssociated with each node.
    // As this info is not in the Node interface we must build it; which we
    // do using vertices for each node, when we addVertex at thes nodes we
    // will not be adding vertices but element tags.

    Vertex **theNodeTagVertices = 0;
    int maxNodNum = 0;
    Node *nodPtr;
    NodeIter &nodeIter = this->getNodes();
    while ((nodPtr = nodeIter()) != 0)
	if (nodPtr->getTag() > maxNodNum)
	    maxNodNum = nodPtr->getTag();

    theNodeTagVertices = new Vertex *[maxNodNum+1];

    if (theNodeTagVertices == 0) {
	opserr << "WARNING Domain::buildEleGraph ";
	opserr << " - Not Enough Memory for NodeTagVertices\n";
	return -1;
    }
    
    for (int l=0; l<=maxNodNum; l++) theNodeTagVertices[l] = 0;

    // now create the vertices with a reference equal to the node number.
    // and a tag which ranges from 0 through numVertex-1 and placed in
    // theNodeTagVertices at a position equal to the node's tag.

    NodeIter &nodeIter2 = this->getNodes();
    count = START_VERTEX_NUM;
    while ((nodPtr = nodeIter2()) != 0) {
	int nodeTag = nodPtr->getTag();
	Vertex *vertexPtr = new Vertex(count++,nodeTag);
	theNodeTagVertices[nodeTag] = vertexPtr;

	if (vertexPtr == 0) {
	    opserr << "WARNING Domain::buildEleGraph";
	    opserr << " - Not Enough Memory to create ";
	    opserr << count << "th Node Vertex\n";
	    delete [] theNodeTagVertices;
	    return -1;
	}
    }

    // now add the the Elements to the nodes

    ElementIter &eleIter3 = this->getElements();

    while((elePtr = eleIter3()) != 0) {
	int eleTag = elePtr->getTag();
	const ID &id = elePtr->getExternalNodes();

	int size = id.Size();
	for (int i=0; i<size; i++) 
	    theNodeTagVertices[id(i)]->addEdge(eleTag);
    }



    // now add the edges to the vertices of our element graph;
    // this is done by looping over the Node vertices, getting their 
    // Adjacenecy and adding edges between elements with common nodes


    Vertex *vertexPtr;
    for (int k=0; k<=maxNodNum; k++)
	if ((vertexPtr = theNodeTagVertices[k]) != 0) {

	    const ID &id = vertexPtr->getAdjacency();

	    int size = id.Size();
	    for (int i=0; i<size; i++) {
		int Element1 = id(i);

		int vertexTag1 = theElementTagVertices[Element1];

		for (int j=0; j<size; j++) 
		    if (i != j) {

			int Element2 = id(j);
			int vertexTag2 = theElementTagVertices[Element2];

			// addEdge() adds for both vertices - do only once


			if (vertexTag1 > vertexTag2) 
			    theEleGraph->addEdge(vertexTag1,vertexTag2);
			    theEleGraph->addEdge(vertexTag2,vertexTag1);			
		    }
	    }
	}

    // done now delete theElementTagVertices, the node Vertices and
    // theNodeTagVertices
   
    delete [] theElementTagVertices;    
    
    for (int i=0; i<=maxNodNum; i++)
	if ((vertexPtr = theNodeTagVertices[i]) != 0) 
	    delete vertexPtr;
	    
    delete [] theNodeTagVertices;
    
    return 0;
    
}

int 
Domain::buildNodeGraph(Graph *theNodeGraph)
{
    int numVertex = this->getNumNodes();

    if (numVertex == 0) {
	return 0;
    }	
	
    // create another vertices array which aids in adding edges
    
    int *theNodeTagVertices = 0;
    int maxNodNum = 0;
    Node *nodPtr;
    NodeIter &nodeIter = this->getNodes();
    while ((nodPtr = nodeIter()) != 0)
	if (nodPtr->getTag() > maxNodNum)
	    maxNodNum = nodPtr->getTag();

    theNodeTagVertices = new int [maxNodNum+1];

    if (theNodeTagVertices == 0) {
	opserr << "WARNING Domain::buildNodeGraph ";
	opserr << " - Not Enough Memory for NodeTagVertices\n";
	return -1;
    }
    
    for (int j=0; j<=maxNodNum; j++) theNodeTagVertices[j] = -1;

    // now create the vertices with a reference equal to the node number.
    // and a tag which ranges from START_VERTEX_NUM through 
    // numNodes+START_VERTEX_NUM

    NodeIter &nodeIter2 = this->getNodes();
    int count = START_VERTEX_NUM;
    while ((nodPtr = nodeIter2()) != 0) {
	int nodeTag = nodPtr->getTag();
	Vertex *vertexPtr = new Vertex(count,nodeTag);

	if (vertexPtr == 0) {
	    opserr << "WARNING Domain::buildNodeGraph";
	    opserr << " - Not Enough Memory to create ";
	    opserr << count << "th Vertex\n";
	    delete [] theNodeTagVertices;
	    return -1;
	}

	// add the vertex to the graph
	theNodeGraph->addVertex(vertexPtr);
	theNodeTagVertices[nodeTag] = count++;
    }

    // now add the edges, by looping over the Elements, getting their
    // IDs and adding edges between all elements who share a node.
    
    Element *elePtr;
    ElementIter &eleIter = this->getElements();

    while((elePtr = eleIter()) != 0) {
	const ID &id = elePtr->getExternalNodes();

	int size = id.Size();
	for (int i=0; i<size; i++) {
	    int node1 = id(i);
	    int vertexTag1 = theNodeTagVertices[node1];
	    
	    for (int j=0; j<size; j++) 
		if (i != j) {

		    int node2 = id(j);
		    int vertexTag2 = theNodeTagVertices[node2];
		    
		    // addEdge() adds for both vertices - do only once
		    if (vertexTag1 > vertexTag2) 
			theNodeGraph->addEdge(vertexTag1,vertexTag2);
		}
	}
    }

    // done now delete theNodeTagVertices
    delete [] theNodeTagVertices;
    
    return 0;
}


int 
Domain::sendSelf(int cTag, Channel &theChannel)
{
  // update the commitTag and currentGeoTag
  commitTag = cTag;

  this->hasDomainChanged();

  // first we send info about the current domain flag and the number of
  // elements, nodes, constraints and load patterns currently in the domain
  int numEle, numNod, numSPs, numMPs, numLPs;
  numNod = theNodes->getNumComponents();
  numEle = theElements->getNumComponents();
  numSPs = theSPs->getNumComponents();
  numMPs = theMPs->getNumComponents();  
  numLPs = theLoadPatterns->getNumComponents();

  ID domainData(11);
  domainData(0) = currentGeoTag;

  domainData(1) = numNod;
  domainData(2) = numEle;
  domainData(3) = numSPs;
  domainData(4) = numMPs;
  domainData(5) = numLPs;

  // add the database tag for the ID's storing node, element, constraints
  // and loadpattern data into domainData
  // NOTE: if these still 0 get new ones from the channel
  if (dbNod == 0) {
    dbNod = theChannel.getDbTag();
    dbEle = theChannel.getDbTag();
    dbSPs = theChannel.getDbTag();
    dbMPs = theChannel.getDbTag();
    dbLPs = theChannel.getDbTag();
  } 

  domainData(6) = dbNod;
  domainData(7) = dbEle;
  domainData(8) = dbSPs;
  domainData(9) = dbMPs;
  domainData(10) = dbLPs;

  if (theChannel.sendID(theDbTag, commitTag, domainData) < 0) {
    opserr << "Domain::send - channel failed to send the initial ID\n";
    return -1;
  }    

  // send the time information
  Vector domainTime(1);
  domainTime(0) = committedTime;

  if (theChannel.sendVector(theDbTag, commitTag, domainTime) < 0) {
    opserr << "Domain::send - channel failed to send the time Vector\n";
    return -2;
  }    

  // now check if data defining the objects in the domain needs to be sent 
  // NOTE THIS APPROACH MAY NEED TO CHANGE FOR VERY LARGE PROBLEMS IF CHANNEL CANNOT
  // HANDLE VERY LARGE ID OBJECTS.
  if (lastGeoSendTag != currentGeoTag) {
    
    //
    // into an ID we are gonna place the class and db tags for each node so can rebuild
    // this ID we then send to the channel
    //

    // create the ID and get the node iter
    if (numNod != 0) {
      ID nodeData(numNod*2);
      Node *theNode;
      NodeIter &theNodes = this->getNodes();
      int loc =0;

      // loop over nodes in domain adding their classTag and dbTag to the ID
      while ((theNode = theNodes()) != 0) {
	nodeData(loc) = theNode->getClassTag();
	int dbTag = theNode->getDbTag();
	
	// if dbTag still 0 get one from Channel; 
	// if this tag != 0 set the dbTag in node
	if (dbTag == 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theNode->setDbTag(dbTag);
	}
      
	nodeData(loc+1) = dbTag;
	loc+=2;
      }    

      // now send the ID
      if (theChannel.sendID(dbNod, currentGeoTag, nodeData) < 0) {
	opserr << "Domain::send - channel failed to send the node ID\n";
	return -2;
      }
    }

    // we do the same for elements as we did for nodes above .. see comments
    // for nodes if you can't figure whats going on!

    if (numEle != 0) {
      ID elementData(numEle*2);
      Element *theEle;
      ElementIter &theElements = this->getElements();
      int loc = 0;
    
      while ((theEle = theElements()) != 0) {
	elementData(loc) = theEle->getClassTag();
	int dbTag = theEle->getDbTag();

	if (dbTag == 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theEle->setDbTag(dbTag);
	}
      
	elementData(loc+1) = dbTag;
	loc+=2;
      }

      // now send the ID
      if (theChannel.sendID(dbEle, currentGeoTag, elementData) < 0) {
	opserr << "Domain::send - channel failed to send the element ID\n";
	return -3;
      }
    }

    // we do the same for SP_Constraints as for Nodes above .. see comments
    // for nodes if you can't figure whats going on!    
    
    if (numSPs != 0) {
      ID spData(numSPs*2);
      SP_Constraint *theSP;
      SP_ConstraintIter &theSPs = this->getSPs();
      int loc = 0;
    
      while ((theSP = theSPs()) != 0) {
	spData(loc) = theSP->getClassTag();
	int dbTag = theSP->getDbTag();

	if (dbTag == 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theSP->setDbTag(dbTag);
	}
	
	spData(loc+1) = dbTag;
	loc+=2;
      }    

      if (theChannel.sendID(dbSPs, currentGeoTag, spData) < 0) {
	opserr << "Domain::send - channel failed to send the SP_Constraint ID\n";
	return -4;
      }
    }

    // we do the same for MP_Constraints as for Nodes above .. see comments
    // for nodes if you can't figure whats going on!    
    
    if (numMPs != 0) {
      ID mpData(numMPs*2);
      MP_Constraint *theMP;
      MP_ConstraintIter &theMPs = this->getMPs();
      int loc = 0;
    
      while ((theMP = theMPs()) != 0) {
	mpData(loc) = theMP->getClassTag();
	int dbTag = theMP->getDbTag();
	
	if (dbTag == 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theMP->setDbTag(dbTag);
	}
      
	mpData(loc+1) = dbTag;
	loc+=2;
      }    

      if (theChannel.sendID(dbMPs, currentGeoTag, mpData) < 0) {
	opserr << "Domain::send - channel failed to send the MP_Constraint ID\n";
	return -5;
      }
    }

    // we do the same for LoadPatterns as we did for Nodes above .. see comments
    // for nodes if you can't figure whats going on!    

    if (numLPs != 0) {
      ID lpData(numLPs*2);
      LoadPattern *theLP;
      LoadPatternIter &theLPs = this->getLoadPatterns();
      int loc = 0;
    
      while ((theLP = theLPs()) != 0) {
	lpData(loc) = theLP->getClassTag();
	int dbTag = theLP->getDbTag();

	if (dbTag == 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theLP->setDbTag(dbTag);
	}
      
	lpData(loc+1) = dbTag;
	loc+=2;
      }    

      if (theChannel.sendID(dbLPs, currentGeoTag, lpData) < 0) {
	opserr << "Domain::send - channel failed to send the LoadPattern ID\n";
	return -6;
      }    
    }


    // now so that we don't do this next time if nothing in the domain has changed
    lastGeoSendTag = currentGeoTag;
  }

  //
  // now we invoke sendSelf on each of the objects .. 
  // NOTE: don't have to set the dbTags of the objects as just done this above
  //

  // send the nodes
  Node *theNode;
  NodeIter &theNodes = this->getNodes();
  while ((theNode = theNodes()) != 0) {
    if (theNode->sendSelf(commitTag, theChannel) < 0) {
      opserr << "Domain::send - node with tag " << theNode->getTag() << " failed in sendSelf\n";
      return -7;
    }
  }

  // send the elements
  Element *theEle;
  ElementIter &theElements = this->getElements();
  while ((theEle = theElements()) != 0) {
    if (theEle->sendSelf(commitTag, theChannel) < 0) {
      opserr << "Domain::send - element with tag " << theEle->getTag() << " failed in sendSelf\n";
      return -8;
    }
  }

  // send the single point constraints
  SP_Constraint *theSP;
  SP_ConstraintIter &theSPs = this->getSPs();
  while ((theSP = theSPs()) != 0) {
    if (theSP->sendSelf(commitTag, theChannel) < 0) {
      opserr << "Domain::send - SP_Constraint with tag " << theSP->getTag() << " failed in sendSelf\n";
      return -9;
    }
  }    

  // send the multi point constraints
  MP_Constraint *theMP;
  MP_ConstraintIter &theMPs = this->getMPs();
  while ((theMP = theMPs()) != 0) {
    if (theMP->sendSelf(commitTag, theChannel) < 0) {
      opserr << "Domain::send - MP_Constraint with tag " << theMP->getTag() << " failed in sendSelf\n";
      return -10;
    }
  }    

  // send the load patterns
  LoadPattern *theLP;
  LoadPatternIter &theLPs = this->getLoadPatterns();
  while ((theLP = theLPs()) != 0) {
    if (theLP->sendSelf(commitTag, theChannel) < 0) {
      opserr << "Domain::send - LoadPattern with tag " << theLP->getTag() << " failed in sendSelf\n";
      return -11;
    }
  }  

  // if get here we were successfull
  return commitTag;
}


int 
Domain::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) 
{
  // set the commitTag in the domain to cTag & update the getTag if needed
  commitTag = cTag;
  this->hasDomainChanged();

  // first we get the data about the state of the domain for this commitTag
  ID domainData(11);
  if (theChannel.recvID(theDbTag, commitTag, domainData) < 0) {
    opserr << "Domain::recv - channel failed to recv the initial ID\n";
    return -1;
  }

  // recv the time information
  Vector domainTime(1);
  if (theChannel.recvVector(theDbTag, commitTag, domainTime) < 0) {
    opserr << "Domain::send - channel failed to recv thetime Vector\n";
    return -1;
  }    

  currentTime = domainTime(0);
  committedTime = currentTime;

  // 
  // now if the currentGeoTag does not agree with whats in the domain
  // we must wipe everything in the domain and recreate the domain based on the info from the channel
  //
  if (currentGeoTag == 0 || domainData(0) != currentGeoTag) {
    
    // set the currrentGeoTag
    int geoTag = domainData(0);

    int i, loc;
    int numEle, numNod, numSPs, numMPs, numLPs;

    // if receiving set lastGeoSendTag to be equal to currentGeoTag 
    // at time all the data was sent if not we must clear out the objects and rebuild
    lastGeoSendTag = domainData(0);

    // clear out the all the components in the current domain
    this->clearAll();

    currentTime = domainTime(0);
    committedTime = currentTime;

    // 
    // now we rebuild the nodes
    //
    
    // first get the information from the domainData about the nodes
    numNod = domainData(1);
    dbNod = domainData(6);
    
    if (numNod != 0) {
      ID nodeData(2*numNod);

      // now receive the ID about the nodes, class tag and dbTags
      if (theChannel.recvID(dbNod, geoTag, nodeData) < 0) {
	opserr << "Domain::recv - channel failed to recv the node ID\n";
	return -2;
      }

      // now for each node we 1) get a new node of the correct type from the ObjectBroker
      // 2) ensure the node exists and set it's dbTag, 3) we invoke recvSelf on this new 
      // blank node and 4) add this node to the domain
      loc = 0;
      for (i=0; i<numNod; i++) {
	int classTag = nodeData(loc);
	int dbTag = nodeData(loc+1);
      
	Node *theNode = theBroker.getNewNode(classTag);

	if (theNode == 0) {
	  opserr << "Domain::recv - cannot create node with classTag " << classTag << endln;
	  return -2;
	}			

	theNode->setDbTag(dbTag);
      
	if (theNode->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "Domain::recv - node with dbTag " << dbTag << " failed in recvSelf\n";
	  return -2;
	}			

	if (this->addNode(theNode) == false) {
	  opserr << "Domain::recv - could not add node with tag " << theNode->getTag() << " into domain\n!";
	  return -3;
	}			

	loc+=2;
      }
    }

    // 
    // now we rebuild the elements .. same as nodes above .. see above if can't understand!!
    //
    
    numEle = domainData(2);
    dbEle = domainData(7);

    if (numEle != 0) {
      ID eleData(2*numEle);

      if (theChannel.recvID(dbEle, geoTag, eleData) < 0) {
	opserr << "Domain::recv - channel failed to recv the Ele ID\n";
	return -2;
      }

      loc = 0;
      for (i=0; i<numEle; i++) {
	int classTag = eleData(loc);
	int dbTag = eleData(loc+1);
      
	Element *theEle = theBroker.getNewElement(classTag);
	if (theEle == 0) {
	  opserr << "Domain::recv - cannot create element with classTag " << classTag << endln;
	  return -2;
	}			
	theEle->setDbTag(dbTag);
      
	if (theEle->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "Domain::recv - Ele with dbTag " << dbTag << " failed in recvSelf()\n";
	  return -2;
	}			

	if (this->addElement(theEle) == false) {
	  opserr << "Domain::recv - could not add Ele with tag " << theEle->getTag() << " into domain!\n";
	  return -3;
	}			

	loc+=2;
      }
    }

    // 
    // now we rebuild the SP_Constraints .. same as nodes above .. see above if can't understand!!
    //
    
    numSPs = domainData(3);
    dbSPs = domainData(8);
    if (numSPs != 0) {
      ID spData(2*numSPs);

      if (theChannel.recvID(dbSPs, geoTag, spData) < 0) {
	opserr << "Domain::recv - channel failed to recv the SP_Constraints ID\n";
	return -2;
      }

      loc = 0;
      for (i=0; i<numSPs; i++) {
	int classTag = spData(loc);
	int dbTag = spData(loc+1);
      
	SP_Constraint *theSP = theBroker.getNewSP(classTag);
	if (theSP == 0) {
	  opserr << "Domain::recv - cannot create SP_Constraint with classTag " << classTag << endln;
	  return -2;
	}			
	theSP->setDbTag(dbTag);
      
	if (theSP->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "Domain::recv - SP_Constraint with dbTag " << dbTag << " failed in recvSelf\n";
	  return -2;
	}			

	if (this->addSP_Constraint(theSP) == false) {
	  opserr << "Domain::recv - could not add SP_Constraint with tag " << theSP->getTag() << " into domain!\n";
	  return -3;
	}			

	loc+=2;
      }
    }


    // 
    // now we rebuild the MP_Constraints .. same as nodes above .. see above if can't understand!!
    //
    
    numMPs = domainData(4);
    dbMPs = domainData(9);

    if (numMPs != 0) {
      ID mpData(2*numMPs);

      if (theChannel.recvID(dbMPs, geoTag, mpData) < 0) {
	opserr << "Domain::recv - channel failed to recv the MP_Constraints ID\n";
	return -2;
      }

      loc = 0;
      for (i=0; i<numMPs; i++) {
	int classTag = mpData(loc);
	int dbTag = mpData(loc+1);
      
	MP_Constraint *theMP = theBroker.getNewMP(classTag);
	if (theMP == 0) {
	  opserr << "Domain::recv - cannot create MP_Constraint with classTag " << classTag << endln;
	  return -2;
	}			
	theMP->setDbTag(dbTag);
      
	if (theMP->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "Domain::recv - MP_Constraint with dbTag " << dbTag << " failed in recvSelf\n";
	  return -2;
	}			

	if (this->addMP_Constraint(theMP) == false) {
	  opserr << "Domain::recv - could not add MP_Constraint with tag " << theMP->getTag() << " into domain!\n";
	  return -3;
	}			
	
	loc+=2;
      }
    }

    // 
    // now we rebuild the LoadPatterns .. same as nodes above .. see above if can't understand!!
    //
    
    numLPs = domainData(5);
    dbLPs = domainData(10);

    if (numLPs != 0) {
      ID lpData(2*numLPs);
      
      if (theChannel.recvID(dbLPs, geoTag, lpData) < 0) {
	opserr << "Domain::recv - channel failed to recv the MP_Constraints ID\n";
	return -2;
      }

      loc = 0;
      for (i=0; i<numLPs; i++) {
	int classTag = lpData(loc);
	int dbTag = lpData(loc+1);

	LoadPattern *theLP = theBroker.getNewLoadPattern(classTag);
	if (theLP == 0) {
	  opserr << "Domain::recv - cannot create MP_Constraint with classTag  " << classTag << endln;
	  return -2;
	}			
	theLP->setDbTag(dbTag);
      
	if (theLP->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "Domain::recv - LoadPattern with dbTag " << dbTag << " failed in recvSelf\n";
	  return -2;
	}			

	if (this->addLoadPattern(theLP) == false) {
	  opserr << "Domain::recv - could not add LoadPattern with tag " << theLP->getTag() <<  " into the Domain\n";
	  return -3;
	}			

	loc+=2;
      }
    }

    // set the currentGeoTag & mark domainChangeFlag as false
    // this way if restoring froma a database and domain has not changed for the analysis
    // the analysis will not have to to do a domainChanged() operation
    currentGeoTag = domainData(0);

    lastGeoSendTag = currentGeoTag;
    hasDomainChangedFlag = false;

  } else {

    // in this block .. we have the components they just have to recv themselves again
    
    Node *theNode;
    NodeIter &theNodes = this->getNodes();
    while ((theNode = theNodes()) != 0) {
      if (theNode->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "Domain::recv - node with tag " << theNode->getTag() << " failed in recvSelf\n";
	return -7;
      }
    }

    Element *theEle;
    ElementIter &theElements = this->getElements();
    while ((theEle = theElements()) != 0) {
      if (theEle->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "Domain::recv - element with tag " << theEle->getTag() <<  " failed in recvSelf\n";
	return -8;
      }
      theEle->update();
    }

    SP_Constraint *theSP;
    SP_ConstraintIter &theSPs = this->getSPs();
    while ((theSP = theSPs()) != 0) {
      if (theSP->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "Domain::recv - SP_Constraint with tag " << theSP->getTag() << " failed in recvSelf\n";
	return -9;
      }
    }    

    MP_Constraint *theMP;
    MP_ConstraintIter &theMPs = this->getMPs();
    while ((theMP = theMPs()) != 0) {
      if (theMP->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "Domain::recv - MP_Constraint with tag " << theMP->getTag() << " failed in recvSelf\n";
	return -10;
      }
    }    

    LoadPattern *theLP;
    LoadPatternIter &theLPs = this->getLoadPatterns();
    while ((theLP = theLPs()) != 0) {
      if (theLP->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "Domain::recv - LoadPattern with tag" << theLP->getTag() << " failed in recvSelf";
	return -11;
      }
    }  
  } 

  // now set the domains lastGeoSendTag and currentDomainChangedFlag
  lastGeoSendTag = currentGeoTag;  

  // if get here we were successfull
  return 0;
}


int
Domain::calculateNodalReactions(bool inclInertia)
{
  Node *theNode;
  Element *theElement;

  NodeIter &theNodes = this->getNodes();
  while ((theNode = theNodes()) != 0) {
    theNode->resetReactionForce(inclInertia);
  }

  ElementIter &theElements = this->getElements();
  while ((theElement = theElements()) != 0)
    theElement->addResistingForceToNodalReaction(inclInertia);
  return 0;
}
