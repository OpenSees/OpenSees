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
                                                                        
// $Revision: 1.61 $
// $Date: 2010-09-16 00:07:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/Domain.cpp,v $
                                                                        
// Written: fmk 
//
// Purpose: This file contains the class definition for Domain
// Domain is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints. These objects are all added to the Domain by a 
// ModelBuilder.
//
// What: "@(#) Domain.C, revA"

#include <stdlib.h>
#include <math.h>

#include <OPS_Globals.h>
#include <Domain.h>
#include <DummyStream.h>

#include <ElementIter.h>
#include <NodeIter.h>
#include <ElementalLoadIter.h>
#include <NodalLoadIter.h>
#include <Element.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <Pressure_Constraint.h>
#include <MP_Constraint.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <Response.h>

#include <MapOfTaggedObjects.h>
#include <MapOfTaggedObjectsIter.h>

#include <SingleDomEleIter.h>
#include <SingleDomNodIter.h>
#include <SingleDomSP_Iter.h>
#include <SingleDomPC_Iter.h>
#include <SingleDomMP_Iter.h>
#include <LoadPatternIter.h>
#include <SingleDomAllSP_Iter.h>
#include <SingleDomParamIter.h>

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>

#include <Vertex.h>
#include <Matrix.h>
#include <Graph.h>
#include <Recorder.h>
#include <MeshRegion.h>
#include <Analysis.h>
#include <FE_Datastore.h>
#include <FEM_ObjectBroker.h>

//
// global variables
//

Domain       *ops_TheActiveDomain = 0;
double        ops_Dt = 0.0;
bool          ops_InitialStateAnalysis = false;
int           ops_Creep = 0;

Domain::Domain()
:theRecorders(0), numRecorders(0),
 currentTime(0.0), committedTime(0.0), dT(0.0), currentGeoTag(0),
 hasDomainChangedFlag(false), theDbTag(0), lastGeoSendTag(-1),
 dbEle(0), dbNod(0), dbSPs(0), dbPCs(0), dbMPs(0), dbLPs(0), dbParam(0),
 eleGraphBuiltFlag(false),  nodeGraphBuiltFlag(false), theNodeGraph(0), 
 theElementGraph(0), 
 theRegions(0), numRegions(0), commitTag(0),
 theBounds(6), theEigenvalues(0), theEigenvalueSetTime(0), 
 theModalDampingFactors(0), inclModalMatrix(false),
 lastChannel(0),
 paramIndex(0), paramSize(0), numParameters(0)
{
  
    // init the arrays for storing the domain components
    theElements = new MapOfTaggedObjects();
    theNodes    = new MapOfTaggedObjects();
    theSPs      = new MapOfTaggedObjects();
    thePCs      = new MapOfTaggedObjects();
    theMPs      = new MapOfTaggedObjects();    
    theLoadPatterns = new MapOfTaggedObjects();
    theParameters   = new MapOfTaggedObjects();

    // init the iters    
    theEleIter = new SingleDomEleIter(theElements);    
    theNodIter = new SingleDomNodIter(theNodes);
    theSP_Iter = new SingleDomSP_Iter(theSPs);
    thePC_Iter = new SingleDomPC_Iter(thePCs);
    theMP_Iter = new SingleDomMP_Iter(theMPs);
    theLoadPatternIter = new LoadPatternIter(theLoadPatterns);
    allSP_Iter = new SingleDomAllSP_Iter(*this);
    theParamIter = new SingleDomParamIter(theParameters);
    
    // check that there was space to create the data structures    
    if (theElements ==0 || theNodes == 0 || 
	theSPs == 0 || theMPs == 0 || thePCs == 0 ||
	theEleIter == 0 || theNodIter == 0 || 
	theMP_Iter == 0 || theSP_Iter == 0 || thePC_Iter == 0 ||
	theLoadPatterns == 0 || theLoadPatternIter == 0 ||
	theParameters == 0) {	

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
 dbEle(0), dbNod(0), dbSPs(0), dbPCs(0), dbMPs(0), dbLPs(0), dbParam(0),
 eleGraphBuiltFlag(false), nodeGraphBuiltFlag(false), theNodeGraph(0), 
 theElementGraph(0),
 theRegions(0), numRegions(0), commitTag(0),
 theBounds(6), theEigenvalues(0), theEigenvalueSetTime(0), 
 theModalDampingFactors(0), inclModalMatrix(false),
 lastChannel(0), paramIndex(0), paramSize(0), numParameters(0)
{
    // init the arrays for storing the domain components
    theElements = new MapOfTaggedObjects();
    theNodes    = new MapOfTaggedObjects();
    theSPs      = new MapOfTaggedObjects();
    thePCs      = new MapOfTaggedObjects();
    theMPs      = new MapOfTaggedObjects();    
    theLoadPatterns = new MapOfTaggedObjects();
    theParameters   = new MapOfTaggedObjects();
    
    // init the iters
    theEleIter = new SingleDomEleIter(theElements);    
    theNodIter = new SingleDomNodIter(theNodes);
    theSP_Iter = new SingleDomSP_Iter(theSPs);
    thePC_Iter = new SingleDomPC_Iter(thePCs);
    theMP_Iter = new SingleDomMP_Iter(theMPs);
    theLoadPatternIter = new LoadPatternIter(theLoadPatterns);
    allSP_Iter = new SingleDomAllSP_Iter(*this);
    theParamIter = new SingleDomParamIter(theParameters);
    
    // check that there was space to create the data structures    
    if (theElements ==0 || theNodes == 0 || 
	theSPs == 0 || theMPs == 0 || thePCs == 0 ||
	theEleIter == 0 || theNodIter == 0 || 
	theMP_Iter == 0 || theSP_Iter == 0 || thePC_Iter == 0 ||
	theLoadPatterns == 0 || theLoadPatternIter == 0 ||
	theParameters == 0) {	

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
 dbEle(0), dbNod(0), dbSPs(0), dbPCs(0), dbMPs(0), dbLPs(0), dbParam(0),
 eleGraphBuiltFlag(false), nodeGraphBuiltFlag(false), theNodeGraph(0), 
 theElementGraph(0), 
 theElements(&theElementsStorage),
 theNodes(&theNodesStorage),
 theSPs(&theSPsStorage),
 theMPs(&theMPsStorage), 
 theLoadPatterns(&theLoadPatternsStorage),
 theRegions(0), numRegions(0), commitTag(0),
 theBounds(6), theEigenvalues(0), theEigenvalueSetTime(0), 
 theModalDampingFactors(0), inclModalMatrix(false),
 lastChannel(0),paramIndex(0), paramSize(0), numParameters(0)
{
    // init the arrays for storing the domain components
    thePCs      = new MapOfTaggedObjects();

    // init the iters    
    theEleIter = new SingleDomEleIter(theElements);    
    theNodIter = new SingleDomNodIter(theNodes);
    theSP_Iter = new SingleDomSP_Iter(theSPs);
    thePC_Iter = new SingleDomPC_Iter(thePCs);
    theMP_Iter = new SingleDomMP_Iter(theMPs);
    theLoadPatternIter = new LoadPatternIter(theLoadPatterns);
    allSP_Iter = new SingleDomAllSP_Iter(*this);
    theParameters   = new MapOfTaggedObjects();    
    theParamIter = new SingleDomParamIter(theParameters);

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
	theSPs == 0 || theMPs == 0 || thePCs == 0 ||
	theEleIter == 0 || theNodIter == 0 ||
	theMP_Iter == 0 || theSP_Iter == 0 || thePC_Iter == 0 ||
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
 dbEle(0), dbNod(0), dbSPs(0), dbPCs(0), dbMPs(0), dbLPs(0), dbParam(0),
 eleGraphBuiltFlag(false), nodeGraphBuiltFlag(false), theNodeGraph(0), 
 theElementGraph(0), 
 theRegions(0), numRegions(0), commitTag(0),
 theBounds(6), theEigenvalues(0), theEigenvalueSetTime(0), 
 theModalDampingFactors(0), inclModalMatrix(false),
 lastChannel(0),paramIndex(0), paramSize(0), numParameters(0)
{
    // init the arrays for storing the domain components
    theStorage.clearAll(); // clear the storage just in case populated
    theElements = &theStorage;
    theNodes    = theStorage.getEmptyCopy();
    theSPs      = theStorage.getEmptyCopy();
    thePCs      = theStorage.getEmptyCopy();
    theMPs      = theStorage.getEmptyCopy();
    theLoadPatterns = theStorage.getEmptyCopy();    
    theParameters   = theStorage.getEmptyCopy();    

    // init the iters    
    theEleIter = new SingleDomEleIter(theElements);    
    theNodIter = new SingleDomNodIter(theNodes);
    theSP_Iter = new SingleDomSP_Iter(theSPs);
    thePC_Iter = new SingleDomPC_Iter(thePCs);
    theMP_Iter = new SingleDomMP_Iter(theMPs);
    theLoadPatternIter = new LoadPatternIter(theLoadPatterns);
    allSP_Iter = new SingleDomAllSP_Iter(*this);
    theParamIter = new SingleDomParamIter(theParameters);

    // check that there was space to create the data structures    
    if (theElements ==0 || theNodes == 0 || 
	theSPs == 0 || theMPs == 0 || thePCs == 0 ||
	theEleIter == 0 || theNodIter == 0 ||
	theMP_Iter == 0 || theSP_Iter == 0 || thePC_Iter == 0 ||
	theLoadPatterns == 0 || theLoadPatternIter == 0 ||
	theParameters == 0) { 
	
	opserr << ("Domain::Domain(ObjectStorage &) - out of memory\n");	
    }
    
    theBounds(0) = 0;
    theBounds(1) = 0;
    theBounds(2) = 0;
    theBounds(3) = 0;
    theBounds(4) = 0;    
    theBounds(5) = 0;            

    dbEle =0; dbNod =0; dbSPs =0; dbPCs = 0; dbMPs =0; dbLPs = 0; dbParam = 0;
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
  this->Domain::clearAll();

  // delete all the storage objects
  // SEGMENT FAULT WILL OCCUR IF THESE OBJECTS WERE NOT CONSTRUCTED
  // USING NEW
  
  if (theElements != 0)
    delete theElements;    
  
  if (theNodes != 0)
    delete theNodes;
  
  if (theSPs != 0)
    delete theSPs;

  if (thePCs != 0)
    delete thePCs;
  
  if (theMPs != 0)
    delete theMPs;
  
  if (theLoadPatterns != 0)
    delete theLoadPatterns;

  if (theParameters != 0)
    delete theParameters;
  
  if (theEleIter != 0)
    delete theEleIter;
  
  if (theNodIter != 0)
    delete theNodIter;
  
  if (theSP_Iter != 0)
    delete theSP_Iter;

  if (thePC_Iter != 0)
    delete thePC_Iter;
  
  if (theMP_Iter != 0)
    delete theMP_Iter;
  
  if (allSP_Iter != 0)
    delete allSP_Iter;
  
  if (theParamIter != 0)
    delete theParamIter;

  if (theEigenvalues != 0)
    delete theEigenvalues;

  if (theLoadPatternIter != 0)
      delete theLoadPatternIter;

  if (theModalDampingFactors != 0)
    delete theModalDampingFactors;
  
  int i;
  for (i=0; i<numRecorders; i++) 
    if (theRecorders[i] != 0)
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
  
  theRecorders = 0;
  numRecorders = 0;
}


// void addElement(Element *);
//	Method to add an element to the model.


bool
Domain::addElement(Element *element)
{
  int eleTag = element->getTag();
  ops_TheActiveElement = element;

  // check all the elements nodes exist in the domain
  const ID &nodes = element->getExternalNodes();
  int numDOF = 0;
  for (int i=0; i<nodes.Size(); i++) {
      int nodeTag = nodes(i);
      Node *nodePtr = this->getNode(nodeTag);
      if (nodePtr == 0) {
	opserr << "WARNING Domain::addElement - In element " << eleTag;
	  opserr << "\n no Node " << nodeTag << " exists in the domain\n";
	  return false;
      }
      numDOF += nodePtr->getNumberDOF();
  }   

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
  //#ifdef _G3DEBUG    
    // check the Node exists in the Domain
    int nodeTag = spConstraint->getNodeTag();
    int dof = spConstraint->getDOF_Number();

    Node *nodePtr = this->getNode(nodeTag);
    if (nodePtr == 0) {
      opserr << "Domain::addSP_Constraint - cannot add as node node with tag" <<
	nodeTag << "does not exist in model\n";       	
      return false;
    }

    // check that the DOF specified exists at the Node
    int numDOF = nodePtr->getNumberDOF();
    if (numDOF < dof) {
	opserr << "Domain::addSP_Constraint - cannot add as node with tag" << 
	  nodeTag << "does not have associated constrained DOF\n"; 
	return false;
    }      
    // #endif

    // check if an existing SP_COnstraint exists for that dof at the node
    bool found = false;
    SP_ConstraintIter &theExistingSPs = this->getSPs();
    SP_Constraint *theExistingSP = 0;
    while ((found == false) && ((theExistingSP = theExistingSPs()) != 0)) {
      int spNodeTag = theExistingSP->getNodeTag();
      int spDof = theExistingSP->getDOF_Number();
      if (nodeTag == spNodeTag && spDof == dof) {
	found = true;
      }
    }
    
    if (found == true) {
	opserr << "Domain::addSP_Constraint - cannot add as node already constrained in that dof by existing SP_Constraint\n";
	spConstraint->Print(opserr);
	return false;
    }

  // check that no other object with similar tag exists in model
  int tag = spConstraint->getTag();
  TaggedObject *other = theSPs->getComponentPtr(tag);
  if (other != 0) {
    opserr << "Domain::addSP_Constraint - cannot add as constraint with tag " << 
      tag << "already exists in model\n";             
    spConstraint->Print(opserr);

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

// void addPressure_Constraint(Pressure_Constraint *);
//	Method to add a constraint to the model.
//

bool
Domain::addPressure_Constraint(Pressure_Constraint *pConstraint)
{
#ifdef _G3DEBUG    
    // check the Node exists in the Domain
    int nodeTag = pConstraint->getTag();
    Node *nodePtr = this->getNode(nodeTag);
    if (nodePtr == 0) {
        opserr << "Domain::addPressure_Constraint - cannot add as node with tag";
        opserr << nodeTag << "does not exist in model\n";
        return false;
    }
#endif

    // check that no other object with similar tag exists in model
    int tag = pConstraint->getTag();
    TaggedObject *other = thePCs->getComponentPtr(tag);
    if (other != 0) {
        opserr << "Domain::addPressure_Constraint - cannot add as constraint with tag";
        opserr << tag << "already exists in model\n";             
        return false;
    }
  
    bool result = thePCs->addComponent(pConstraint);
    if (result == false) {
        opserr << "Domain::addPressure_Constraint - cannot add constraint with tag";
        opserr << tag << "to the container\n";
        return false;
    } 

    pConstraint->setDomain(this);
    this->domainChange();  

    return true;
}


int
Domain::addSP_Constraint(int axisDirn, double axisValue, 
			 const ID &fixityConditions, double tol)
{
  if (axisDirn < 0)
    return -1;

  NodeIter &theNodes = this->getNodes();
  Node *theNode;
  int numAddedSPs = 0;
  // for each node in the domain
  while ((theNode = theNodes()) != 0) {
    const Vector &theCrds = theNode->getCrds();
    int sizeCrds = theCrds.Size();
    int numDOF = theNode->getNumberDOF();
    int nodeTag = theNode->getTag();

    // check it has crds in axis specified
    if (axisDirn < sizeCrds) {
      double nodeCrdDirn = theCrds(axisDirn);

      // check if coordinate is within tol of value given
      if (fabs(nodeCrdDirn-axisValue) <= tol) {

	// foreach dof to be constrained create 
	for (int i=0; i<fixityConditions.Size(); i++) {
	  if ((i < numDOF) && (fixityConditions(i) == 1)) {

	    // check if an existing SP_COnstraint exists for that dof at the node
	    bool found = false;
	    SP_ConstraintIter &theExistingSPs = this->getSPs();
	    SP_Constraint *theExistingSP = 0;
	    while ((found == false) && ((theExistingSP = theExistingSPs()) != 0)) {
	      int spNodeTag = theExistingSP->getNodeTag();
	      int dof = theExistingSP->getDOF_Number();
	      if (nodeTag == spNodeTag && i == dof) {
		found = true;
	      }
	    }
	    
	    // if no sp constraint, create one and ass it
	    if (found == false) {

	      SP_Constraint *theSP = new SP_Constraint(nodeTag, i, 0.0, true);
	      
	      if (this->addSP_Constraint(theSP) == false) {
		opserr << "WARNING could not add SP_Constraint to domain for node " << theNode->getTag();
		delete theSP;
	      } else {
		numAddedSPs++;
	      }
	    } else {
	      ; // opserr << "Domain::addSP(AXIS) - constraint exists at node\n";
	    }
	  }
	}
      }
    }
  }
  this->domainChange();
  return numAddedSPs;
}


// void addMP_Constraint(MP_Constraint *);
//	Method to add a constraint to the model.
//

bool
Domain::addMP_Constraint(MP_Constraint *mpConstraint)
{
//#ifdef _G3DEBUG
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
    //#endif

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
Domain::addParameter(Parameter *theParam)
{
  int paramTag = theParam->getTag();

  if (paramTag == 0) {
    // don't add it .. just invoke setDomain on the parameter
    theParam->setDomain(this);
    return true;
  }

  // check if a Parameter with a similar tag already exists in the Domain
  TaggedObject *other = theParameters->getComponentPtr(paramTag);
  if (other != 0) {
    opserr << "Domain::addParameter - parameter with tag " << paramTag << "already exists in model\n"; 
    return false;
  }

  // add the param to the container object for the parameters
  bool result = theParameters->addComponent(theParam);

  if (result == false) {
    opserr << "Domain::addParameter - parameter " << paramTag << "could not be added to container\n";
    theParam->setDomain(this);
    return result;
  }

  // mark the Domain as having been changed
  //    this->domainChange();
  
  // Array is full or empty
  if (numParameters == paramSize) {
    
    // Increase size and allocate new array
    paramSize += paramSize_grow;
    int *tmp_paramIndex = new int[paramSize];
    
    // Copy values from old array to new
    for (int i = 0; i < numParameters; i++)
      tmp_paramIndex[i] = paramIndex[i];
    
    // Get rid of old array
    delete [] paramIndex;
    
    // Set pointer to new array
    paramIndex = tmp_paramIndex;
  }

  // Add to index
  paramIndex[numParameters] = paramTag;
  theParam->setGradIndex(numParameters);
  numParameters++;    

  if (strcmp(theParam->getType(),"FEModel") != 0) {
    //theParam->setGradIndex(-1);
  }

  theParam->setDomain(this);
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
      opserr << "Domain::addElementalLoad() - no pattern with tag " << pattern << 
	"exits in  the model, not adding the ele load " << *load << endln;

	return false;
    }
    LoadPattern *theLoadPattern = (LoadPattern *)thePattern;
    bool result = theLoadPattern->addElementalLoad(load);
    if (result == false) {
      opserr << "Domain::addElementalLoad() - no pattern with tag" << 
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
  thePCs->clearAll();
  theMPs->clearAll();
  theLoadPatterns->clearAll();
  theParameters->clearAll();
  numParameters = 0;

  // remove the recorders
  int i;
  for (i=0; i<numRecorders; i++)
	  if (theRecorders[i] != 0)
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

  this->setModalDampingFactors(0);

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
  
  dbEle =0; dbNod =0; dbSPs =0; dbPCs = 0; dbMPs =0; dbLPs = 0; dbParam = 0;

  currentGeoTag = 0;
  lastGeoSendTag = -1;
  lastChannel = 0;

  // rest the flag to be as initial
  hasDomainChangedFlag = false;
  nodeGraphBuiltFlag = false;
  eleGraphBuiltFlag = false;

  if (theNodeGraph != 0)
    delete theNodeGraph;
  theNodeGraph = 0;

  if (theElementGraph != 0)
    delete theElementGraph;
  theElementGraph = 0;
  
  dbEle =0; dbNod =0; dbSPs =0; dbPCs = 0; dbMPs =0; dbLPs = 0; dbParam = 0;
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


int
Domain::removeSP_Constraint(int theNode, int theDOF, int loadPatternTag)
{
  SP_Constraint *theSP =0;
  bool found = false;
  int spTag = 0;

  if (loadPatternTag == -1) {
    SP_ConstraintIter &theSPs = this->getSPs();
    while ((found == false) && ((theSP = theSPs()) != 0)) {
      int nodeTag = theSP->getNodeTag();
      int dof = theSP->getDOF_Number();
      if (nodeTag == theNode && dof == theDOF) {
	spTag = theSP->getTag();
	found = true;
      }
    }
    
  } else {

    LoadPattern *thePattern = this->getLoadPattern(loadPatternTag);
    if (thePattern != 0) {
      SP_ConstraintIter &theSPs = thePattern->getSPs();
      while ((found == false) && ((theSP = theSPs()) != 0)) {
	int nodeTag = theSP->getNodeTag();
	int dof = theSP->getDOF_Number();
	if (nodeTag == theNode && dof == theDOF) {
	  spTag = theSP->getTag();
	  found = true;
	}
      }
    }
  }

  if (found == true)
    theSP = this->removeSP_Constraint(spTag);

  // mark the domain has having changed regardless if SP constrain was there or not
  this->domainChange();

  if (theSP != 0) {
    delete theSP;
    return 1;
  }
   
  return 0;
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

Pressure_Constraint *
Domain::removePressure_Constraint(int tag)
{
    // remove the object from the container    
    TaggedObject *mc = thePCs->removeComponent(tag);
    
    // if not there return 0    
    if (mc == 0) 
	return 0;

    // mark the domain as having changed    
    this->domainChange();
    
    // perform a downward cast, set the objects domain pointer to 0
    // and return the result of the cast    
    Pressure_Constraint *result = (Pressure_Constraint *)mc;
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


int
Domain::removeMP_Constraints(int nodeTag)
{
  ID tagsToRemove(0); int sizeTags = 0;
  MP_Constraint *theMP = 0;
  MP_ConstraintIter &theMPIter = this->getMPs();
  while ((theMP = theMPIter()) != 0) {
    int cNode = theMP->getNodeConstrained();
    if (cNode == nodeTag) {
      int mpTag = theMP->getTag();
      tagsToRemove[sizeTags] = mpTag;
      sizeTags++;
    }
  }

  if (sizeTags == 0)
    return 0;

  for (int i=0; i<sizeTags; i++) {
    int  tag = tagsToRemove(i);
    TaggedObject *mc = theMPs->removeComponent(tag);
    if (mc != 0)
      delete mc;
  }
    
  // mark the domain as having changed    
  this->domainChange();
    
  return sizeTags;
}    


Parameter *
Domain::removeParameter(int tag)
{
  Parameter *theParam = (Parameter*) theParameters->getComponentPtr(tag);

  if (theParam != 0) {

    // Find where RV is located
    int index;
    for (index = 0; index < numParameters; index++) {
      if (paramIndex[index] == tag)
	break;
    }

    // Shift indices down by one
    for (int i = index; i < numParameters-1; i++) {
      paramIndex[i] = paramIndex[i+1];
      Parameter *otherParam = this->getParameterFromIndex(i);
      otherParam->setGradIndex(i);
    }

    // Now remove the component
    theParameters->removeComponent(tag);

    numParameters--;
  }

  return 0;

  /*
  // remove the object from the container    
  TaggedObject *mc = theParameters->removeComponent(tag);
  
  // if not there return 0
  if (mc == 0) 
      return 0;

  // otherwise mark the domain as having changed
  this->domainChange();
  
  // perform a downward cast to an Element (safe as only Element added to
  // this container, 0 the Elements DomainPtr and return the result of the cast  
  Parameter *result = (Parameter *)mc;
  //  result->setDomain(0);
  return result;
  */
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
    return *theSP_Iter;
}

Pressure_ConstraintIter &
Domain::getPCs()
{
    thePC_Iter->reset();
    return *thePC_Iter;
}

SP_ConstraintIter &
Domain::getDomainAndLoadPatternSPs()
{
    allSP_Iter->reset();
    return *allSP_Iter;
}


MP_ConstraintIter &
Domain::getMPs()
{
    theMP_Iter->reset();
    return *theMP_Iter;
}


LoadPatternIter &
Domain::getLoadPatterns()
{
    theLoadPatternIter->reset();
    return *theLoadPatternIter;
}

ParameterIter &
Domain::getParameters()
{
    theParamIter->reset();
    return *theParamIter;
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

Pressure_Constraint *
Domain::getPressure_Constraint(int tag) 
{
    TaggedObject *mc = thePCs->getComponentPtr(tag);

    // if not there return 0 otherwise perform a cast and return that  
    if (mc == 0) 
        return 0;
    Pressure_Constraint *result = (Pressure_Constraint *)mc;
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

Parameter *
Domain::getParameter(int tag) 
{
  TaggedObject *mc = theParameters->getComponentPtr(tag);

  // if not there return 0 otherwise perform a cast and return that
  if (mc == 0) 
      return 0;
  Parameter *result = (Parameter *)mc;
  return result;
}

Parameter *
Domain::getParameterFromIndex(int index)
{
  if (index >= 0 && index < numParameters)
    return this->getParameter(paramIndex[index]);

  else {
    opserr << "Domain::getParameterFromIndex -- index " << index << " out of bounds 0 ... " << numParameters-1 << endln;
    return 0;
  }

}

int
Domain::getParameterIndex(int tag)
{
  int index;

  // Find index of RV with specified tag
  for (index = 0; index < numParameters; index++) {
    if (paramIndex[index] == tag)
      break;
  }

  if (index == numParameters) {
    opserr << "Domain::getParameterIndex -- parameter with tag " << tag << " not found" << endln;
    return -1;
  }

  return index;
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
Domain::getNumPCs(void) const
{
    return thePCs->getNumComponents();
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

int 
Domain::getNumParameters(void) const
{
  return theParameters->getNumComponents();
}

const Vector &
Domain::getPhysicalBounds(void)
{
    return theBounds;
}

const Vector *
Domain::getNodeResponse(int nodeTag, NodeResponseType responseType)
{
  Node *theNode = this->getNode(nodeTag);
  if (theNode == 0)
    return NULL;
  else 
    return theNode->getResponse(responseType);
}


static Vector responseData(0);

const Vector *
Domain::getElementResponse(int eleTag, const char **argv, int argc)
{
  Element *theEle = this->getElement(eleTag);
  if (theEle == 0)
    return NULL;
  else  {

    if (argc == 1) {
      if (strcmp(argv[0],"forces") == 0) {
	return &(theEle->getResistingForce());
      } else if (strcmp(argv[0],"nodeTags") == 0) {
	const ID&theNodes = theEle->getExternalNodes();
	int size = theNodes.Size();
	if (responseData.Size() != size) 
	  responseData.resize(size);
	for (int i=0; i<size; i++)
	  responseData(i) = theNodes(i);
	return &responseData;
      }
    }
	
    DummyStream dummy;
    Response *theResponse = theEle->setResponse(argv, argc, dummy);
    if (theResponse == 0) {
      return 0;	  
    }

    if (theResponse->getResponse() < 0) {
      delete theResponse;
      return 0;
    }

    Information &eleInfo = theResponse->getInformation();
    //const Vector *data = &(eleInfo.getData());
    responseData = eleInfo.getData();
    delete theResponse;
    return &responseData;
  }
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
Domain::clearElementGraph(void) {
  if (theElementGraph != 0)
    delete theElementGraph;
  theElementGraph = 0;

  eleGraphBuiltFlag = false;
}

void
Domain::clearNodeGraph(void) {
  if (theNodeGraph != 0)
    delete theNodeGraph;
  theNodeGraph = 0;

  nodeGraphBuiltFlag = false;
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
Domain::setCreep(int newCreep)
{
  ops_Creep = newCreep;
}

int
Domain::getCreep(void) const
{
  return ops_Creep;
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
    while ((theSP = theSPs()) != 0) {
      theSP->applyConstraint(timeStep);
    }

    ops_Dt = dT;
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


void
Domain::unsetLoadConstant(void)
{
    // loop over all the load patterns that are currently added to the domain
    // getting them to set their loads as now constant
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = this->getLoadPatterns();
    while((thePattern = thePatterns()) != 0)
      thePattern->unsetLoadConstant();
}


int
Domain::initialize(void)
{
  Element *elePtr;
  ElementIter &theElemIter = this->getElements();    
  while ((elePtr = theElemIter()) != 0) 
    // lvalue needed here for M$ VC++ compiler -- MHS
	// and either the  VS2011 or intel compiler does not like it!
#ifndef _VS2011
    Matrix initM(elePtr->getInitialStiff());

#else
	 elePtr->getInitialStiff();
#endif


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
Domain::record(bool fromAnalysis)
{
  int res = 0;

  // invoke record on all recorders
  for (int i=0; i<numRecorders; i++)
    if (theRecorders[i] != 0)
      res += theRecorders[i]->record(commitTag, currentTime);
  
  // update the commitTag
  commitTag++;

  return res;
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
      if (theRecorders[i] != 0)
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
      if (theRecorders[i] != 0)
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
    ops_TheActiveElement = theEle;
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
Domain::updateParameter(int tag, int value)
{
  // get the object from the container 
  TaggedObject *mc = theParameters->getComponentPtr(tag);
  
  // if not there return 0
  if (mc == 0) 
      return 0;

  // convert to a parameter & update
  Parameter *result = (Parameter *)mc;
  int res = result->update(value);

  return res;
}

int
Domain::updateParameter(int tag, double value)
{
  // remove the object from the container    
  TaggedObject *mc = theParameters->getComponentPtr(tag);
  
  // if not there return 0
  if (mc == 0) {
	  opserr << "Domain::updateParameter(int tag, double value) - parameter with tag not present\n";
      return 0;
  }

  Parameter *theParam = (Parameter *)mc;
  int res =  theParam->update(value);
  return res;
}


int
Domain::analysisStep(double dT)
{
  return 0;
}

int
Domain::eigenAnalysis(int nuMode, bool generalized, bool findSmallest)
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

int
Domain::setModalDampingFactors(Vector *theValues, bool inclMatrix)
{
  // if theValues == 0, turn off modal damping
  if (theValues == 0) {
    if (theModalDampingFactors != 0)
      delete theModalDampingFactors;
    theModalDampingFactors = 0;
    inclModalMatrix = inclMatrix;
    return 0;
  }

  // make sure the eigen value vector is large enough
  if (theModalDampingFactors == 0 || theModalDampingFactors->Size() != theValues->Size()) {
    if (theModalDampingFactors != 0)
      delete theModalDampingFactors;
    theModalDampingFactors = new Vector(*theValues);
  } else {
    *theModalDampingFactors = *theValues;
  }

  inclModalMatrix = inclMatrix;

  return 0;
}

const Vector *
Domain::getModalDampingFactors(void)
{
  return theModalDampingFactors;
}

bool
Domain::inclModalDampingMatrix(void)
{
  return inclModalMatrix;
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


bool 
Domain::getDomainChangeFlag(void)
{
  return hasDomainChangedFlag;
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
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {

    s << "\t\"properties\": {\n";

    OPS_printUniaxialMaterial(s, flag);
    s << ",\n";   
    OPS_printNDMaterial(s, flag);
    s << ",\n";
    OPS_printSectionForceDeformation(s, flag);
    s << ",\n";   
    OPS_printCrdTransf(s, flag);      

    s << "\n\t},\n";
	s << "\t\"geometry\": {\n";

    int numToPrint = theNodes->getNumComponents();
    NodeIter &theNodess = this->getNodes();
    Node *theNode;
    int numPrinted = 0;
    s << "\t\t\"nodes\": [\n";
    while ((theNode = theNodess()) != 0) {    
      theNode->Print(s, flag);
      numPrinted += 1;
      if (numPrinted < numToPrint)
	s << ",\n";
      else
	s << "\n\t\t],\n";
    }


    Element *theEle;
    ElementIter &theElementss = this->getElements();
    numToPrint = theElements->getNumComponents();
    numPrinted = 0;
    s << "\t\t\"elements\": [\n";
    while ((theEle = theElementss()) != 0) {
      theEle->Print(s, flag);
      numPrinted += 1;
      if (numPrinted < numToPrint)
	s << ",\n";
      else
	s << "\n\t\t]\n";
      }

	s << "\t}\n";
	s << "}\n";
    s << "}\n";

	return;
  }
      
  
  s << "Current Domain Information\n";
  s << "\tCurrent Time: " << currentTime;
  s << "\ntCommitted Time: " << committedTime << endln;    
  s << "NODE DATA: NumNodes: " << theNodes->getNumComponents() << "\n";  
  theNodes->Print(s, flag);
  
  s << "ELEMENT DATA: NumEle: " << theElements->getNumComponents() << "\n";
  theElements->Print(s, flag);
  
  s << "\nSP_Constraints: numConstraints: " << theSPs->getNumComponents() << "\n";
  theSPs->Print(s, flag);
  
  s << "\nPressure_Constraints: numConstraints: " << thePCs->getNumComponents() << "\n";
  thePCs->Print(s, flag);
  
  s << "\nMP_Constraints: numConstraints: " << theMPs->getNumComponents() << "\n";
  theMPs->Print(s, flag);
  
  s << "\nLOAD PATTERNS: numPatterns: " << theLoadPatterns->getNumComponents() << "\n\n";
  theLoadPatterns->Print(s, flag);
  
  s << "\nPARAMETERS: numParameters: " << theParameters->getNumComponents() << "\n\n";
  theParameters->Print(s, flag);
}


void Domain::Print(OPS_Stream &s, ID *nodeTags, ID *eleTags, int flag)
{
  if (nodeTags != 0) {
    int numNodes = nodeTags->Size();
    for (int i=0; i<numNodes; i++) {
      int nodeTag = (*nodeTags)(i);
      TaggedObject *theNode = theNodes->getComponentPtr(nodeTag);
      if (theNode != 0)
	theNode->Print(s, flag);
    }
  }

  if (eleTags != 0) {
    int numEle = eleTags->Size();
    for (int i=0; i<numEle; i++) {
      int eleTag = (*eleTags)(i);
      TaggedObject *theEle = theElements->getComponentPtr(eleTag);
      if (theEle != 0)
	theEle->Print(s, flag);
    }
  }
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

  for (int i=0; i<numRecorders; i++) {
    if (theRecorders[i] == 0) {
      theRecorders[i] = &theRecorder;
      return 0;
    }
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
      if (theRecorders[i] != 0)
	delete theRecorders[i];
    
    if (theRecorders != 0) {
      delete [] theRecorders;
    }
  
    theRecorders = 0;
    numRecorders = 0;
    return 0;
}

int
Domain::removeRecorder(int tag)
{
  for (int i=0; i<numRecorders; i++) {
    if (theRecorders[i] != 0) {
      if (theRecorders[i]->getTag() == tag) {
	delete theRecorders[i];
	theRecorders[i] = 0;
	return 0;
      }
    }    
  }
  
  return -1;
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

void
Domain::getRegionTags(ID& rtags) const
{
    rtags.resize(numRegions);
    for(int i=0; i<numRegions; i++) {
        rtags(i) = theRegions[i]->getTag();
    }

}

typedef map<int, int> MAP_INT;
typedef MAP_INT::value_type   MAP_INT_TYPE;
typedef MAP_INT::iterator     MAP_INT_ITERATOR;

typedef map<int, ID *> MAP_ID;
typedef MAP_ID::value_type   MAP_ID_TYPE;
typedef MAP_ID::iterator     MAP_ID_ITERATOR;


int 
Domain::buildEleGraph(Graph *theEleGraph)
{
   // see if quick return
    int numVertex = this->getNumElements();
    if (numVertex == 0)
        return 0;

    //
    // iterate over the lements of the domain
    //  create a vertex with a unique tag for each element
    //  also create a map to hold element tag - vertex tag mapping
    //

    MAP_INT theEleToVertexMap;
    MAP_INT_ITERATOR theEleToVertexMapEle;

    Element *theEle;
    ElementIter &theElements = this->getElements();
    int count = START_VERTEX_NUM;
    while ((theEle = theElements()) != 0) {
      int eleTag = theEle->getTag();
      Vertex *vertexPtr = new Vertex(count, eleTag);

      if (vertexPtr == 0) {
        opserr << "WARNING Domain::buildEleGraph - Not Enough Memory to create the " << count << " vertex\n";
        return -1;
      }

      theEleGraph->addVertex(vertexPtr);
      theEleToVertexMapEle = theEleToVertexMap.find(eleTag);
      if (theEleToVertexMapEle == theEleToVertexMap.end()) {
        theEleToVertexMap.insert(MAP_INT_TYPE(eleTag, count));

        // check if successfully added
        theEleToVertexMapEle = theEleToVertexMap.find(eleTag);
        if (theEleToVertexMapEle == theEleToVertexMap.end()) {
          opserr << "Domain::buildEleGraph - map STL failed to add object with tag : " << eleTag << endln;
          return false;
        }

        count++;
      }
    }

    //
    // We now need to determine which elements are associated with each node.
    // As this info is not in the Node interface we must build it;
    //
    // again we will use an stl map, index will be nodeTag, object will be Vertex
    // do using vertices for each node, when we addVertex at these nodes we
    // will not be adding vertices but element tags.
    //

    MAP_ID theNodeToVertexMap;
    MAP_ID_ITERATOR theNodeEle;

    Node *nodPtr;

    // now create the vertices with a reference equal to the node number.
    // and a tag which ranges from 0 through numVertex-1 and placed in
    // theNodeTagVertices at a position equal to the node's tag.

    NodeIter &theNodes = this->getNodes();
    while ((nodPtr = theNodes()) != 0) {
      int nodeTag = nodPtr->getTag();
      ID *eleTags = new ID(0, 4);

      if (eleTags == 0) {
        opserr << "WARNING Domain::buildEleGraph - Not Enough Memory to create the " << count << " vertex\n";
        return -1;
      }

      theNodeEle = theNodeToVertexMap.find(nodeTag);
      if (theNodeEle == theNodeToVertexMap.end()) {
        theNodeToVertexMap.insert(MAP_ID_TYPE(nodeTag, eleTags));

        // check if successfully added
        theNodeEle = theNodeToVertexMap.find(nodeTag);
        if (theNodeEle == theNodeToVertexMap.end()) {
          opserr << "Domain::buildEleGraph - map STL failed to add object with tag : " << nodeTag << endln;
          return false;
        }
      }
    }

    // now add the the Elements to the node vertices

    ElementIter &eleIter3 = this->getElements();

    while((theEle = eleIter3()) != 0) {
      int eleTag = theEle->getTag();
      const ID &id = theEle->getExternalNodes();

      int size = id.Size();
      for (int i=0; i<size; i++) {
        int nodeTag = id(i);

        MAP_ID_ITERATOR theNodeEle;
        theNodeEle = theNodeToVertexMap.find(nodeTag);
        if (theNodeEle == theNodeToVertexMap.end()) {
          return -1;
        } else {
          ID *theNodeEleTags = (*theNodeEle).second;
          theNodeEleTags->insert(eleTag);
        }
      }
    }

    //
    // now add the edges to the vertices of our element graph;
    // this is done by looping over the Node vertices, getting their
    // Adjacenecy and adding edges between elements with common nodes
    //

    MAP_ID_ITERATOR currentComponent;
    currentComponent = theNodeToVertexMap.begin();
    while (currentComponent != theNodeToVertexMap.end()) {
      ID *id = (*currentComponent).second;

      int size = id->Size();
      for (int i=0; i<size; i++) {
        int eleTag1 = (*id)(i);

        theEleToVertexMapEle = theEleToVertexMap.find(eleTag1);
        if (theEleToVertexMapEle != theEleToVertexMap.end()) {
          int vertexTag1 = (*theEleToVertexMapEle).second;

          for (int j=0; j<size; j++)
            if (i != j) {
	      int eleTag2 = (*id)(j);
              theEleToVertexMapEle = theEleToVertexMap.find(eleTag2);
              if (theEleToVertexMapEle != theEleToVertexMap.end()) {
                int vertexTag2 = (*theEleToVertexMapEle).second;

                // addEdge() adds for both vertices - do only once

                if (vertexTag1 > vertexTag2) {
                  theEleGraph->addEdge(vertexTag1,vertexTag2);
		  theEleGraph->addEdge(vertexTag2,vertexTag1);
		}
              }
            }
        }
      }
      currentComponent++;
    }

    // clean up - delete the ID's associated with the nodes
    currentComponent = theNodeToVertexMap.begin();
    while (currentComponent != theNodeToVertexMap.end()) {
      delete (*currentComponent).second;
      currentComponent++;
    }

    return 0;

}

int
Domain::buildNodeGraph(Graph *theNodeGraph)
{
    int numVertex = this->getNumNodes();

    if (numVertex == 0) {
	return 0;
    }	
	
    Node *nodPtr;
    MAP_INT theNodeTagVertices;

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
  int numEle, numNod, numSPs, numPCs, numMPs, numLPs, numParam;
  numNod = theNodes->getNumComponents();
  numEle = theElements->getNumComponents();
  numSPs = theSPs->getNumComponents();
  numPCs = thePCs->getNumComponents();
  numMPs = theMPs->getNumComponents();  
  numLPs = theLoadPatterns->getNumComponents();
  numParam = theParameters->getNumComponents();

  ID domainData(15);
  domainData(0) = currentGeoTag;
  domainData(1) = numNod;
  domainData(2) = numEle;
  domainData(3) = numSPs;
  domainData(13) = numPCs;
  domainData(4) = numMPs;
  domainData(5) = numLPs;
  domainData(11) = numParam;

  // add the database tag for the ID's storing node, element, constraints
  // and loadpattern data into domainData
  // NOTE: if these still 0 get new ones from the channel
  if (dbNod == 0) {
    dbNod = theChannel.getDbTag();
    dbEle = theChannel.getDbTag();
    dbSPs = theChannel.getDbTag();
    dbPCs = theChannel.getDbTag();
    dbMPs = theChannel.getDbTag();
    dbLPs = theChannel.getDbTag();
    dbParam = theChannel.getDbTag();
  } 

  domainData(6) = dbNod;
  domainData(7) = dbEle;
  domainData(8) = dbSPs;
  domainData(14) = dbPCs;
  domainData(9) = dbMPs;
  domainData(10) = dbLPs;
  domainData(12) = dbParam;

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

  /*
  if (theChannel.isDatastore() == 1) {
    static ID theLastSendTag(1);
    if (theChannel.recvID(0,0,theLastSendTag) == 0)
      lastGeoSendTag = theLastSendTag(0);
    else
      lastGeoSendTag = -1;
  }
  */

  if (lastChannel != theChannel.getTag() || lastGeoSendTag != currentGeoTag) {

    lastChannel = theChannel.getTag();
    
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
    // for nodes if you can't figure what's going on!

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
    // for nodes if you can't figure what's going on!    
    
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

    // we do the same for Pressure_Constraints as for Nodes above .. see comments
    // for nodes if you can't figure what's going on!    
    
    if (numPCs != 0) {
        ID pData(numPCs*2);
        Pressure_Constraint *thePC;
        Pressure_ConstraintIter &thePCs = this->getPCs();
        int loc = 0;
    
        while ((thePC = thePCs()) != 0) {
            pData(loc) = thePC->getClassTag();
            int dbTag = thePC->getDbTag();

            if (dbTag == 0) {// go get a new tag and setDbTag in ele if this not 0 
                dbTag = theChannel.getDbTag();
                if (dbTag != 0)
                    thePC->setDbTag(dbTag);
            }
	
            pData(loc+1) = dbTag;
            loc+=2;
        }    
        
        if (theChannel.sendID(dbPCs, currentGeoTag, pData) < 0) {
            opserr << "Domain::send - channel failed to send the Pressure_Constraint ID\n";
            return -4;
        }
    }

    // we do the same for MP_Constraints as for Nodes above .. see comments
    // for nodes if you can't figure what's going on!    
    
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
    // for nodes if you can't figure what's going on!    


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

    if (numParam != 0) {
      ID paramData(numParam*2);
      Parameter *theP;
      ParameterIter &theParameters = this->getParameters();
      int loc = 0;
    
      while ((theP = theParameters()) != 0) {
	paramData(loc) = theP->getClassTag();
	int dbTag = theP->getDbTag();

	if (dbTag == 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theP->setDbTag(dbTag);
	}
      
	paramData(loc+1) = dbTag;
	loc+=2;
      }    

      if (theChannel.sendID(dbLPs, currentGeoTag, paramData) < 0) {
	opserr << "Domain::send - channel failed to send the LoadPattern ID\n";
	return -6;
      }    
  }
    // now so that we don't do this next time if nothing in the domain has changed
    lastGeoSendTag = currentGeoTag;
    /*
    if (theChannel.isDatastore() == 1) {
      static ID theLastSendTag(1);
      theLastSendTag(0) = lastGeoSendTag;
      theChannel.sendID(0,0, theLastSendTag);
    }
    */
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

  // send the pressure constraints
  Pressure_Constraint *thePC;
  Pressure_ConstraintIter &thePCs = this->getPCs();
  while ((thePC = thePCs()) != 0) {
      if (thePC->sendSelf(commitTag, theChannel) < 0) {
          opserr << "Domain::send - Pressure_Constraint with tag " << thePC->getTag() << " failed in sendSelf\n";
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

  // send the parameters
  Parameter *theParam;
  ParameterIter &theParams = this->getParameters();
  while ((theParam = theParams()) != 0) {
    if (theParam->sendSelf(commitTag, theChannel) < 0) {
      opserr << "Domain::send - Parameter with tag " << theParam->getTag() << " failed in sendSelf\n";
      return -12;
    }
  }  

  // if get here we were successful
  return commitTag;
}


int 
Domain::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) 
{
  // set the commitTag in the domain to cTag & update the getTag if needed
  commitTag = cTag;
  this->hasDomainChanged();

  // first we get the data about the state of the domain for this commitTag
  ID domainData(15);
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
  // now if the currentGeoTag does not agree with what's in the domain
  // we must wipe everything in the domain and recreate the domain based on the info from the channel
  //

  /*
  if (theChannel.isDatastore() == 1) {
    static ID theLastSendTag(1);
    if (theChannel.recvID(0,0,theLastSendTag) == 0)
      lastGeoSendTag = theLastSendTag(0);
  }
  */

  if (currentGeoTag == 0 || lastChannel != theChannel.getTag() || domainData(0) != currentGeoTag) {

    lastChannel = theChannel.getTag();
    
    // set the currrentGeoTag
    int geoTag = domainData(0);

    int i, loc;
    int numEle, numNod, numSPs, numPCs, numMPs, numLPs;

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
    // now we rebuild the Pressure_Constraints .. same as nodes above .. see above if can't understand!!
    //
    
    numPCs = domainData(13);
    dbPCs = domainData(14);
    if (numPCs != 0) {
      ID pData(2*numPCs);

      if (theChannel.recvID(dbPCs, geoTag, pData) < 0) {
          opserr << "Domain::recv - channel failed to recv the Pressure_Constraints ID\n";
          return -2;
      }

      loc = 0;
      for (i=0; i<numPCs; i++) {
          int classTag = pData(loc);
          int dbTag = pData(loc+1);
      
          Pressure_Constraint *thePC = theBroker.getNewPC(classTag);
          if (thePC == 0) {
              opserr << "Domain::recv - cannot create Pressure_Constraint with classTag " << classTag << endln;
              return -2;
          }			
          thePC->setDbTag(dbTag);
      
          if (thePC->recvSelf(commitTag, theChannel, theBroker) < 0) {
              opserr << "Domain::recv - Pressure_Constraint with dbTag " << dbTag << " failed in recvSelf\n";
              return -2;
          }	

          if (this->addPressure_Constraint(thePC) == false) {
              opserr << "Domain::recv - could not add Pressure_Constraint with tag " << thePC->getTag() << " into domain!\n";
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


    int numParameters = domainData(11);
    int dbParameters = domainData(12);

    if (numParameters != 0) {
      ID paramData(2*numParameters);
      
      if (theChannel.recvID(dbParameters, geoTag, paramData) < 0) {
	opserr << "Domain::recv - channel failed to recv the MP_Constraints ID\n";
	return -2;
      }

      loc = 0;
      for (i=0; i<numParameters; i++) {
	int classTag = paramData(loc);
	int dbTag = paramData(loc+1);

	Parameter *theParameter = theBroker.getParameter(classTag);
	if (theParameter == 0) {
	  opserr << "Domain::recv - cannot create MP_Constraint with classTag  " << classTag << endln;
	  return -2;
	}			
	theParameter->setDbTag(dbTag);
      
	if (theParameter->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  opserr << "Domain::recv - Parameter with dbTag " << dbTag << " failed in recvSelf\n";
	  return -2;
	}			

	if (this->addParameter(theParameter) == false) {
	  opserr << "Domain::recv - could not add LoadPattern with tag " << theParameter->getTag() <<  " into the Domain\n";
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

    Pressure_Constraint *thePC;
    Pressure_ConstraintIter &thePCs = this->getPCs();
    while ((thePC = thePCs()) != 0) {
        if (thePC->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "Domain::recv - Pressure_Constraint with tag " << thePC->getTag() << " failed in recvSelf\n";
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

    Parameter *theParam;
    ParameterIter &theParams = this->getParameters();
    while ((theParam = theParams()) != 0) {
      if (theParam->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "Domain::recv - Parameter with tag" << theParam->getTag() << " failed in recvSelf";
	return -12;
      }
    }  
  } 

  // now set the domains lastGeoSendTag and currentDomainChangedFlag
  lastGeoSendTag = currentGeoTag;  

  // if get here we were successful
  return 0;
}


double
Domain::getNodeDisp(int nodeTag, int dof, int &errorFlag)
{
  double result = 0.0;
  errorFlag = 0;
  Node *theNode = this->getNode(nodeTag);
  if (theNode == 0) {
    errorFlag = -1;
    return 0.0;
  }
  const Vector &disp = theNode->getTrialDisp();
  if (dof < disp.Size() && dof >= 0) {
    result = disp(dof); 
  }  
  
  return result;
}

int 
Domain::setMass(const Matrix &mass, int nodeTag)
{
  Node *theNode = this->getNode(nodeTag);
  if (theNode == 0) {
    return -1;
  }
  return theNode->setMass(mass);  
}

int
Domain::calculateNodalReactions(int flag)
{

  // apply load again! (for case ele load removed and record before an analysis)
  this->applyLoad(committedTime);

  Node *theNode;
  Element *theElement;
  NodeIter &theNodes = this->getNodes();
  while ((theNode = theNodes()) != 0) {
    theNode->resetReactionForce(flag);
  }

  ElementIter &theElements = this->getElements();
  while ((theElement = theElements()) != 0)
    if (theElement->isSubdomain() == false)
      theElement->addResistingForceToNodalReaction(flag);
  return 0;
}

//added by SAJalali
Recorder*
Domain::getRecorder(int tag)
{
	Recorder* res = 0;

	// invoke record on all recorders
	for (int i = 0; i < numRecorders; i++)
	{
		if (theRecorders[i] == 0)
			break;
		if (theRecorders[i]->getTag() == tag)
		{
			res = theRecorders[i];
			break;
		}
	}

	return res;
}




int Domain::activateElements(const ID& elementList)
{
    ElementIter& iter = getElements();
    Element* theElement;
    for (int i = 0; i < elementList.Size(); ++i)
    {
        int eleTag = elementList(i);
        Element* theElement = this->getElement(eleTag);
        if (theElement != 0)
        {
            theElement->activate();
        }
    }
    return 0;
}



int Domain::deactivateElements(const ID& elementList)
{
    // ElementIter& iter = getElements();
    Element* theElement;
    for (int i = 0; i < elementList.Size(); ++i)
    {
        int eleTag = elementList(i);
        Element* theElement = this->getElement(eleTag);
        if (theElement != 0)
        {
            theElement->deactivate();
        }
    }
    return 0;
}
