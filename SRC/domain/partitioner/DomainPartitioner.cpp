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
                                                                        
// $Revision: 1.11 $
// $Date: 2009-10-12 23:51:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/partitioner/DomainPartitioner.cpp,v $
                                                                        
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for DomainPartitioner.
// A DomainPartitioner is an object used to partition a PartitionedDomain.
//
// What: "@(#) DomainPartitioner.C, revA"

#include <DomainPartitioner.h>
 
#include <stdlib.h>
#include <GraphPartitioner.h>
#include <PartitionedDomain.h>
#include <Subdomain.h>
#include <SubdomainIter.h>
#include <Node.h>
#include <Element.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <MP_ConstraintIter.h>
#include <SP_ConstraintIter.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <Graph.h>
#include <Vector.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>
#include <NodalLoadIter.h>
#include <ElementalLoadIter.h>
#include <LoadBalancer.h>
#include <LoadPatternIter.h>
#include <LoadPattern.h>
 
#include <Timer.h>

#include <MapOfTaggedObjects.h>

class NodeLocations: public TaggedObject
{
public:
  NodeLocations(int tag);
  void Print(OPS_Stream &s, int flag =0);  
  int addPartition(int partition);
  ID nodePartitions;
  int numPartitions;
};



NodeLocations::NodeLocations(int tag)
:TaggedObject(tag), 
 nodePartitions(0,1), 
 numPartitions(0)
{

}

void 
NodeLocations::Print(OPS_Stream &s, int flag)
{
  s << "NodeLocations tag: " << this->getTag() << " partitions: " << nodePartitions;
}

int
NodeLocations::addPartition(int partition)
{
  if (nodePartitions.insert(partition) != 1)
    numPartitions++;
  return 0;
}

DomainPartitioner::DomainPartitioner(GraphPartitioner &theGraphPartitioner)
:myDomain(0),thePartitioner(theGraphPartitioner),theBalancer(0),
 theElementGraph(0), theBoundaryElements(0), 
 theNodeLocations(0),elementPlace(0), numPartitions(0), partitionFlag(false), usingMainDomain(false)
{

}    

DomainPartitioner::DomainPartitioner(GraphPartitioner &theGraphPartitioner,
				     LoadBalancer &theLoadBalancer)
:myDomain(0),thePartitioner(theGraphPartitioner),theBalancer(&theLoadBalancer),
 theElementGraph(0), theBoundaryElements(0),
 theNodeLocations(0),elementPlace(0), numPartitions(0), partitionFlag(false), usingMainDomain(false)
{
    // set the links the loadBalancer needs
    theLoadBalancer.setLinks(*this);
}    


DomainPartitioner::~DomainPartitioner()
{
  if (theBoundaryElements != 0) {
    for (int i=0; i<numPartitions; i++)
      if (theBoundaryElements[i] != 0)
	delete theBoundaryElements[i];
    delete []theBoundaryElements;
  }
}

void 
DomainPartitioner::setPartitionedDomain(PartitionedDomain &theDomain)
{
    myDomain = &theDomain;
}

int
DomainPartitioner::partition(int numParts, bool usingMain, int mainPartitionTag, int specialElementTag)
{

  usingMainDomain = usingMain;
  mainPartition = mainPartitionTag;

  // first we ensure the partitioned domain has numpart subdomains
  // with tags 1 through numparts
  for (int i=1; i<=numParts; i++) {
    if (i != mainPartition) {
      Subdomain *subdomainPtr = myDomain->getSubdomainPtr(i);
      if (subdomainPtr == 0) {
	opserr << "DomainPartitioner::partition - No Subdomain: ";
	opserr << i << " exists\n";
	return -1;
      }
    }
  }

  // we get the ele graph from the domain and partition it
  //    Graph &theEleGraph = myDomain->getElementGraph();
  //    theElementGraph = new Graph(myDomain->getElementGraph());

  theElementGraph = &(myDomain->getElementGraph());

  int theError = thePartitioner.partition(*theElementGraph, numParts);

  if (theError < 0) {
    opserr << "DomainPartitioner::partition";
    opserr << " - the graph partioner failed to partition the ";
    opserr << "element graph\n";
    return -10+theError;
  }

  /* print graph */
  //  opserr << "DomainPartitioner::partition - eleGraph: \n";
  //  theElementGraph->Print(opserr, 4);
  
  VertexIter &theVertices1 = theElementGraph->getVertices();
  Vertex *vertexPtr = 0;
  bool moreThanOne = false;
  
  vertexPtr = theVertices1();
  int vertexOnePartition  = 0;
  if (vertexPtr != 0)
    vertexOnePartition  = vertexPtr->getColor();  
  while ((moreThanOne == false) && ((vertexPtr = theVertices1()) != 0)) {
    int partition = vertexPtr->getColor();
    if (partition != vertexOnePartition ) {
      moreThanOne = true;
    }
  }

  if (moreThanOne == false) {
    opserr <<"DomainPartitioner::partition - too few elements for model to be partitioned\n";
    return -1;
  }

  int specialElementColor = 1;
  if (specialElementTag != 0) {
    bool found = false;
    VertexIter &theVerticesSpecial = theElementGraph->getVertices();
    while ((found == false) && ((vertexPtr = theVerticesSpecial()) != 0)) {
      int eleTag = vertexPtr->getRef();
      if (eleTag == specialElementTag) {
	found = true;
	int vertexColor = vertexPtr->getColor();
	if (vertexColor != 1)
	  //	  specialElementColor = vertexColor;
	  vertexPtr->setColor(1);
      }
    }
  }
  
      
  // we create empty graphs for the numParts subdomains,
  // in the graphs we place the vertices for the elements on the boundaries
  
  // we do not invoke the destructor on the individual graphs as 
  // this would invoke the destructor on the individual vertices

  if (theBoundaryElements != 0)
    delete [] theBoundaryElements;
  
  theBoundaryElements = new Graph * [numParts];
  if (theBoundaryElements == 0) {
    opserr << "DomainPartitioner::partition(int numParts)";
    opserr << " - ran out of memory\n";
    numPartitions = 0;  
    return -1;
  }

  for (int l=0; l<numParts; l++) {
    theBoundaryElements[l] = new Graph(2048); // graphs can grow larger; just an estimate
    
    if (theBoundaryElements[l] == 0) {
      opserr << "DomainPartitioner::partition(int numParts)";
      opserr << " - ran out of memory\n";
      numPartitions = 0;
      return -1;
    }
  }
  
  numPartitions = numParts;

  //  opserr << "DomainPartitioner::partition() - nodes \n";  
  
  // we now create a MapOfTaggedObjectStorage to store the NodeLocations
  // and create a new NodeLocation for each node; adding it to the map object

  theNodeLocations = new MapOfTaggedObjects();
  if (theNodeLocations == 0) {
    opserr << "DomainPartitioner::partition(int numParts)";
    opserr << " - ran out of memory creating MapOfTaggedObjectStorage for node locations\n";
    numPartitions = 0;
    return -1;
  }

  NodeIter &theNodes = myDomain->getNodes();
  Node *nodePtr;
  while ((nodePtr = theNodes()) != 0) {
    NodeLocations *theNodeLocation = new NodeLocations(nodePtr->getTag());
    if (theNodeLocation == 0) {
      opserr << "DomainPartitioner::partition(int numParts)";
      opserr << " - ran out of memory creating NodeLocation for node: " << nodePtr->getTag() << endln;
      numPartitions = 0;
      return -1;
    }
    if (theNodeLocations->addComponent(theNodeLocation) == false) {
      opserr << "DomainPartitioner::partition(int numParts)";
      opserr << " - failed to add NodeLocation to Map for Node: " << nodePtr->getTag() << endln;
      numPartitions = 0;
      return -1;
    }
  }

  //
  // we now iterate through the vertices of the element graph
  // to see if the vertex is a boundary vertex or not - if it is
  // we add to the appropriate graph created above. We also set the
  // value the color variable of each of the external nodes connected 
  // to the element to a value which will indicate that that node will
  // have to be added to the subdomain.
  //
  
  VertexIter &theVertexIter = theElementGraph->getVertices();
  while ((vertexPtr = theVertexIter()) != 0) {
    int eleTag = vertexPtr->getRef();
    int vertexColor = vertexPtr->getColor();
    
    const ID &adjacency = vertexPtr->getAdjacency();
    int size = adjacency.Size();
    for (int i=0; i<size; i++) {
      Vertex *otherVertex = theElementGraph->getVertexPtr(adjacency(i));
      if (otherVertex->getColor() != vertexColor) {
	theBoundaryElements[vertexColor-1]->addVertex(vertexPtr,false);
	i = size;
      }
    }
    
    Element *elePtr = myDomain->getElement(eleTag);
    const ID &nodes = elePtr->getExternalNodes();
    size = nodes.Size();
    for (int j=0; j<size; j++) {
      int nodeTag = nodes(j);
      TaggedObject *theTaggedObject = theNodeLocations->getComponentPtr(nodeTag);
      if (theTaggedObject == 0) {
	opserr << "DomainPartitioner::partition(int numParts)";
	opserr << " - failed to find NodeLocation in Map for Node: " << nodePtr->getTag() << " -- A BUG!!\n";
	numPartitions = 0;
	return -1;	
      }
      NodeLocations *theNodeLocation = (NodeLocations *)theTaggedObject;
      theNodeLocation->addPartition(vertexColor);
    }
  }

  // now go through the MP_Constraints and ensure the retained node is in every 
  // partition the constrained node is in
  MP_ConstraintIter &theMPs = myDomain->getMPs();
  MP_Constraint *mpPtr;
  while ((mpPtr = theMPs()) != 0) {
    int retained = mpPtr->getNodeRetained();
    int constrained = mpPtr->getNodeConstrained();
    
    TaggedObject *theRetainedObject = theNodeLocations->getComponentPtr(retained);      
    TaggedObject *theConstrainedObject = theNodeLocations->getComponentPtr(constrained);
    
    if (theRetainedObject == 0 || theConstrainedObject == 0) {
      opserr << "DomainPartitioner::partition(int numParts)";
      if (theRetainedObject == 0)
	opserr << " - failed to find NodeLocation in Map for Node: " << retained << " -- A BUG!!\n";
      if (theConstrainedObject == 0)
	opserr << " - failed to find NodeLocation in Map for Node: " << constrained << " -- A BUG!!\n";
      numPartitions = 0;
      return -1;	
    }
    
    NodeLocations *theRetainedLocation = (NodeLocations *)theRetainedObject;
    NodeLocations *theConstrainedLocation = (NodeLocations *)theConstrainedObject;
    ID &theConstrainedNodesPartitions = theConstrainedLocation->nodePartitions;
    int numPartitions = theConstrainedNodesPartitions.Size();
    for (int i=0; i<numPartitions; i++) {
      theRetainedLocation->addPartition(theConstrainedNodesPartitions(i));
    }
  }

  // we now add the nodes, 
  TaggedObjectIter &theNodeLocationIter = theNodeLocations->getComponents();
  TaggedObject *theNodeObject;

  while ((theNodeObject = theNodeLocationIter()) != 0) {
    NodeLocations *theNodeLocation = (NodeLocations *)theNodeObject;

    int nodeTag = theNodeLocation->getTag();
    ID &nodePartitions = theNodeLocation->nodePartitions;
    int numPartitions = theNodeLocation->numPartitions;

    for (int i=0; i<numPartitions; i++) {
      int partition = nodePartitions(i);	  
      if (partition != mainPartition) {      
	Subdomain *theSubdomain = myDomain->getSubdomainPtr(partition); 
	if (numPartitions == 1) {
	  Node *nodePtr = myDomain->removeNode(nodeTag);
	  theSubdomain->addNode(nodePtr);
	} else {
	  Node *nodePtr = myDomain->getNode(nodeTag);
	  theSubdomain->addExternalNode(nodePtr);	  
	}
      }
    }
  }

  // we now move the elements 
  VertexIter &theVertices = theElementGraph->getVertices();
  while ((vertexPtr = theVertices()) != 0) {
    // move the element
    int partition = vertexPtr->getColor();
    if (partition != mainPartition) {          
      int eleTag = vertexPtr->getRef();

      //      opserr << "removing ele: " << eleTag << endln;
      
      Element *elePtr = myDomain->removeElement(eleTag);  
      //      opserr << *elePtr;

      if (elePtr != 0) {
	//	opserr << "adding ele - start\n";
	Subdomain *theSubdomain = myDomain->getSubdomainPtr(partition);  
	theSubdomain->addElement(elePtr);

	//	opserr << "adding ele - done\n";
      } else {
	opserr << "DomainPartitioner::partioner - element GONE! - eleTag " << eleTag << endln;
      }
    } 
  }

  // now we go through the load patterns and move NodalLoad
  // 1) make sure each subdomain has a copy of the partitioneddomains load patterns.
  // 2) move nodal loads
  // 3) move SP_Constraints
  
  LoadPatternIter &theLoadPatterns = myDomain->getLoadPatterns();
  LoadPattern *theLoadPattern;
  while ((theLoadPattern = theLoadPatterns()) != 0) {
    int loadPatternTag = theLoadPattern->getTag();

    
    // check that each subdomain has a loadPattern with a similar tag and class tag
    for (int i=1; i<=numParts; i++) {
      if (i != mainPartition) {
	Subdomain *theSubdomain = myDomain->getSubdomainPtr(i);
	LoadPattern *loadPatternCopy = theSubdomain->getLoadPattern(loadPatternTag);
	if (loadPatternCopy == 0) {
	  LoadPattern *newLoadPattern = theLoadPattern->getCopy();
	  if (newLoadPattern == 0) {
	    opserr << "DomaiPartitioner::partition - out of memory creating LoadPatterns\n";
 	    return -1;
	  }
	  theSubdomain->addLoadPattern(newLoadPattern);
	}
      }
    }

    // now remove any nodal loads that correspond to internal nodes in a subdomain
    // and add them to the appropriate loadpattern in the subdomain
    
    NodalLoadIter &theNodalLoads = theLoadPattern->getNodalLoads();
    NodalLoad *theNodalLoad;
    while ((theNodalLoad = theNodalLoads()) != 0) {
      int nodeTag = theNodalLoad->getNodeTag();

      TaggedObject *theTaggedObject = theNodeLocations->getComponentPtr(nodeTag);
      if (theTaggedObject == 0) {
	opserr << "DomainPartitioner::partition(int numParts)";
	opserr << " - failed to find NodeLocation in Map for Node: " << nodeTag << " -- A BUG!!\n";
	numPartitions = 0;
	return -1;	
      }
    
      NodeLocations *theNodeLocation = (NodeLocations *)theTaggedObject;
      ID &nodePartitions = theNodeLocation->nodePartitions;
      int numPartitions = theNodeLocation->numPartitions;
      for (int i=0; i<numPartitions; i++) {
	int partition = nodePartitions(i);	  
	if (partition != mainPartition) {      
	  if (numPartitions == 1) {
	    Subdomain *theSubdomain = myDomain->getSubdomainPtr(partition);
	    theLoadPattern->removeNodalLoad(theNodalLoad->getTag());
	    if ((theSubdomain->addNodalLoad(theNodalLoad, loadPatternTag)) != true)
	      opserr << "DomainPartitioner::partition() - failed to add Nodal Load\n";
	  }
	}
      }      
    }

  
    SP_ConstraintIter &theSPs = theLoadPattern->getSPs();
    SP_Constraint *spPtr;
    while ((spPtr = theSPs()) != 0) {
      int nodeTag = spPtr->getNodeTag();
      
      TaggedObject *theTaggedObject = theNodeLocations->getComponentPtr(nodeTag);
      if (theTaggedObject == 0) {
	opserr << "DomainPartitioner::partition(int numParts)";
	opserr << " - failed to find NodeLocation in Map for Node: " << nodeTag << " -- A BUG!!\n";
	numPartitions = 0;
	return -1;	
      }
      
      NodeLocations *theNodeLocation = (NodeLocations *)theTaggedObject;
      ID &nodePartitions = theNodeLocation->nodePartitions;
      int numPartitions = theNodeLocation->numPartitions;
      for (int i=0; i<numPartitions; i++) {
	int partition = nodePartitions(i);	  
	if (partition != mainPartition) {      
	  Subdomain *theSubdomain = myDomain->getSubdomainPtr(partition); 
	  if (numPartitions == 1) 
	    theLoadPattern->removeSP_Constraint(spPtr->getTag());
	  int res = theSubdomain->addSP_Constraint(spPtr, loadPatternTag);
	  if (res < 0)
	    opserr << "DomainPartitioner::partition() - failed to add SP Constraint\n";
	}
      }    
    }  

    ElementalLoadIter &theLoads = theLoadPattern->getElementalLoads();
    ElementalLoad *theLoad;
    while ((theLoad = theLoads()) != 0) {
      int loadEleTag = theLoad->getElementTag();

      SubdomainIter &theSubdomains = myDomain->getSubdomains();
      Subdomain *theSub;
      bool added = false;
      while (((theSub = theSubdomains()) != 0) && (added == false)) {
	bool res = theSub->hasElement(loadEleTag);
	if (res == true) {
	  theLoadPattern->removeElementalLoad(theLoad->getTag());
	  theSub->addElementalLoad(theLoad, loadPatternTag);
	  if (res < 0)
	    opserr << "DomainPartitioner::partition() - failed to add ElementalLoad\n";
	  added = true;
	}
      }   
    }
  }

  // add the single point constraints, 
  
  SP_ConstraintIter &theDomainSP = myDomain->getSPs();
  SP_Constraint *spPtr;
  while ((spPtr = theDomainSP()) != 0) {
    int nodeTag = spPtr->getNodeTag();

    TaggedObject *theTaggedObject = theNodeLocations->getComponentPtr(nodeTag);
    if (theTaggedObject == 0) {
      opserr << "DomainPartitioner::partition(int numParts)";
      opserr << " - failed to find NodeLocation in Map for Node: " << nodeTag << " -- A BUG!!\n";
      numPartitions = 0;
      return -1;	
    }
    
    NodeLocations *theNodeLocation = (NodeLocations *)theTaggedObject;
    ID &nodePartitions = theNodeLocation->nodePartitions;
    int numPartitions = theNodeLocation->numPartitions;
    for (int i=0; i<numPartitions; i++) {
      int partition = nodePartitions(i);	  

      if (partition != mainPartition) {      
	Subdomain *theSubdomain = myDomain->getSubdomainPtr(partition); 
	if (numPartitions == 1) {
	  myDomain->removeSP_Constraint(spPtr->getTag());
	}
	int res = theSubdomain->addSP_Constraint(spPtr);
	if (res < 0)
	  opserr << "DomainPartitioner::partition() - failed to add SP Constraint\n";
      }
    }    
  }  

  // move MP_Constraints - add an MP_Constraint to every partition a constrained node is in
  MP_ConstraintIter &moreMPs = myDomain->getMPs();
  while ((mpPtr = moreMPs()) != 0) {
    int constrained = mpPtr->getNodeConstrained();
    TaggedObject *theConstrainedObject = theNodeLocations->getComponentPtr(constrained);
    NodeLocations *theConstrainedLocation = (NodeLocations *)theConstrainedObject;
    ID &theConstrainedNodesPartitions = theConstrainedLocation->nodePartitions;
    int numPartitions = theConstrainedLocation->numPartitions;
    for (int i=0; i<numPartitions; i++) {
      int partition = theConstrainedNodesPartitions(i);
      if (partition != mainPartition) {
	Subdomain *theSubdomain = myDomain->getSubdomainPtr(partition);
	if (numPartitions == 1) 
	  myDomain->removeMP_Constraint(mpPtr->getTag());
	int res = theSubdomain->addMP_Constraint(mpPtr);
	if (res < 0)
	  opserr << "DomainPartitioner::partition() - failed to add MP Constraint\n";
      }
    }
  }

  // now we go through all the subdomains and tell them to update
  // their analysis for the new layouts
  
  SubdomainIter &theSubDomains = myDomain->getSubdomains();
  Subdomain *theSubDomain;
  while ((theSubDomain = theSubDomains()) != 0) 
    theSubDomain->domainChange();
  
  // we invoke change on the PartitionedDomain
  myDomain->domainChange();

  myDomain->clearElementGraph();
    
  // we are done
  partitionFlag = true;

  return 0;
}


int
DomainPartitioner::balance(Graph &theWeightedPGraph)
{
    int res = 0;

    // check that the object did the partitioning
    if (partitionFlag == false) {
  opserr << "DomainPartitioner::balance(const Vector &load)";
  opserr << " - not partitioned or DomainPartitioner did not partition\n";
  return -1;
    }

    if (theBalancer != 0) {

	// call on the LoadBalancer to partition		
	res = theBalancer->balance(theWeightedPGraph);
	    
	// now invoke domainChanged on Subdomains and PartitionedDomain
	SubdomainIter &theSubDomains = myDomain->getSubdomains();
	Subdomain *theSubDomain;

	while ((theSubDomain = theSubDomains()) != 0) 
	  theSubDomain->domainChange();
	
	// we invoke change on the PartitionedDomain
	myDomain->domainChange();
    }

    return res;
}



int 
DomainPartitioner::getNumPartitions(void) const
{
    return numPartitions;
}



Graph &
DomainPartitioner::getPartitionGraph(void)
{
    if (myDomain == 0) {
      opserr << "ERROR: DomainPartitioner::getPartitionGraph(void)";
      opserr << " - No domain has been set";
      exit(-1);
    }
    return myDomain->getSubdomainGraph();
}

Graph &
DomainPartitioner::getColoredGraph(void)
{
    if (myDomain == 0) {
      opserr << "ERROR: DomainPartitioner::getPartitionGraph(void)";
      opserr << " - No domain has been set";
      exit(0);
    }
    
    return myDomain->getElementGraph();
}


int 
DomainPartitioner::swapVertex(int from, int to, int vertexTag,
			      bool adjacentVertexNotInOther)
{
  opserr << "DomainPartitioner::swapVertex() " << from << " " << to << " " << vertexTag << endln;

  /*
  // check that the object did the partitioning
  if (partitionFlag == false) {
    opserr << "DomainPartitioner::balance(const Vector &load)";
    opserr << " - not partitioned or DomainPartitioner did not partition\n";
    return -1;
  }
  
  // check that the subdomain exist in partitioned domain
  Subdomain *fromSubdomain = myDomain->getSubdomainPtr(from);
  if (fromSubdomain == 0) {
    opserr << "DomainPartitioner::swapVertex - No from Subdomain: ";
    opserr << from << " exists\n";
    return -2;
  }
  Subdomain *toSubdomain = myDomain->getSubdomainPtr(to);
  if (toSubdomain == 0) {
    opserr << "DomainPartitioner::swapVertex - No to Subdomain: ";
    opserr << to << " exists\n";
    return -3;
  }    
  
  // remove the vertex from the boundary if there
  // if not there we have to get a pointer to it from the ele graph.
  
  //    Graph &theEleGraph = myDomain->getElementGraph();    
  Graph *fromBoundary = theBoundaryElements[from-1];
  Graph *toBoundary = theBoundaryElements[to-1];    
  Vertex *vertexPtr;
  
  // get a pointer to the vertex in the element graph
  if (adjacentVertexNotInOther  == false) {
    vertexPtr = fromBoundary->removeVertex(vertexTag,false);    
    if (vertexPtr == 0) 
      vertexPtr = theElementGraph->getVertexPtr(vertexTag);

    if (vertexPtr == 0)  // if still 0 no vertex given by tag exists
      return -4;
  } else { // get a pointer and check vertex not adjacent to another partition
    vertexPtr = fromBoundary->getVertexPtr(vertexTag);    
    if (vertexPtr == 0) 
      vertexPtr = theElementGraph->getVertexPtr(vertexTag);
    if (vertexPtr == 0)  // if still 0 no vertex given by tag exists
      return -4;	
    const ID &adjacent = vertexPtr->getAdjacency();
    bool inTo = false;
    bool inOther = false;
    int adjacentSize = adjacent.Size();
    for (int i=0; i<adjacentSize; i++) {
      Vertex *other = theElementGraph->getVertexPtr(adjacent(i));
      if (other->getColor() == to) 
	inTo = true;
      else if (other->getColor() != from) {
	inOther = true;
	i = adjacentSize;
      }
    }
    if (inTo != true || inOther == true) // we cannot remove the vertex
      return -5;
    // if can remove the vertex from the fromBoundary
    Vertex *theVertex = fromBoundary->removeVertex(vertexTag,false);
  }
  
  int eleTag = vertexPtr->getRef();
  
  //
  // in the FROM subdomain we:
  //  1. remove the element
  //  2. remove all nodes connected to the element divide color by prime
  //  3. see if have to add nodes back as external
  //  4. remove from the PartitionedDomain if node was external
  //  5. add new elements to the boundary vertices for from
  
  //  1. remove the element from the fromSubdomain
    //          this gives us a pointer to the element as well.

    Element *elePtr = fromSubdomain->removeElement(eleTag);
    if (elePtr == 0) // if ele not there we can just return but ERROR it should be
  return -6;

    //  2. remove all nodes connected to the element divide color by prime    
    //  3. add back if still in subdomain
    //  4. if node external remove from PartitionedDomains external nodes
    int primeFrom = primes(from);
    const ID  &nodes1 = elePtr->getExternalNodes(); // create a copy, need after send ele
    ID  nodes(nodes1);

    int nodesSize = nodes.Size();
    Node **nodesArray = new Node *[nodesSize];
    for (int i=0; i<nodesSize; i++) {
  int nodeTag = nodes(i);

  // remove the node
  Node *nodePtr = fromSubdomain->removeNode(nodeTag);

  (*nodePlace)(nodeTag) /= primeFrom;
  nodesArray[i] = nodePtr;
        // opserr << "DomainPartitioner::swapVertex -NODE " << *nodePtr;
  // add back as external if still needed
  if ((*nodePlace)(nodeTag)%primeFrom == 0)
      fromSubdomain->addExternalNode(nodePtr);

  // if it was external remove from PartitionedDomains external nodes
  int partition = (*nodePlace)(nodeTag);
  for (int j=1; j<=numPartitions; j++) 
      if (j != from) {
    int prime = primes(j);
    if (partition%prime == 0) {
        myDomain->removeExternalNode(nodeTag);
        j = numPartitions+1;
    }
      }
    }

    // 5. add new elements to boundary vertices of from if connected to
    //    vertex we just removed and of color of from and not already in boundary

    const ID &eleAdjacent = vertexPtr->getAdjacency();
    int eleAdjacentSize = eleAdjacent.Size();

    for (int a=0; a<eleAdjacentSize; a++) {
	int otherEleVertexTag = eleAdjacent(a);
	Vertex *other = fromBoundary->getVertexPtr(otherEleVertexTag);
	if (other == 0) {
	    other = theElementGraph->getVertexPtr(otherEleVertexTag);
	    if (other->getColor() == from)
		fromBoundary->addVertex(other,false);
	}
    }


    
    // in the TO subdomain we:
    //  1. remove all nodes connected to the element, may or may not be there
    //     if there don't divide by primeTo. 
    //  2. mark all connecting nodes of elements as belonging to To
    //  3. add the nodes, checking if external or internal
    //  3. add the element
    //  4. change the vertex color to be to and add vertex to boundary
    //     also see if we can reduce the size of boundary
    //     of boundary of to vertices

    //  1. remove all nodes connected to the element, may or may not be there
    //      if there divide by primeTo, do this as may have been added as extnl
    int primeTo = primes(to);
    for (int l=0; l<nodesSize; l++) {    
  int nodeTag = nodesArray[l]->getTag();
  int nodeColor = (*nodePlace)(nodeTag);
  if (nodeColor % primeTo == 0) {
      Node *nodePtr;
      nodePtr = toSubdomain->removeNode(nodes(l));
  }
    }


    for (int m=0; m<nodesSize; m++) {
  Node *nodePtr = nodesArray[m];
  int nodeTag = nodePtr->getTag();
  
  //  2. mark all connecting nodes of elements as belonging to To      
  (*nodePlace)(nodeTag) *= primeTo;

  //  3. add the nodes, checking if external or internal  
  int internal = 0;
  for (int n=1; n<=numPartitions; n++) 
      if (n != to) {
    int prime = primes(n);
    if ((*nodePlace)(nodeTag)%prime == 0) { // its external
        internal = 1;
        n = numPartitions+1;
    }
      }
  
  if (internal == 0)
      toSubdomain->addNode(nodePtr);
  else {
      toSubdomain->addExternalNode(nodePtr);
      myDomain->addNode(nodePtr);
  }
    }

    delete [] nodesArray; // DELETE nodesArray here should not use nodePtrs
                          // after here as node objects deleted if parallel    
    
    //  3. add the element
    
    toSubdomain->addElement(elePtr);
    
    //  4. change the vertex color to be to and add vertex to boundary
    //     also see if we can reduce the size of to's boundary vertices


    vertexPtr->setColor(to);
    toBoundary->addVertex(vertexPtr,false);

    for (int n=0; n<eleAdjacentSize; n++) {
	int otherEleVertexTag = eleAdjacent(n);
	Vertex *other = toBoundary->removeVertex(otherEleVertexTag,false);
	if (other != 0) {
	    const ID &othersAdjacency = other->getAdjacency();
	    int otherSize = othersAdjacency.Size();
	    for (int b=0; b<otherSize; b++) {
		int anotherEleVertexTag = othersAdjacency(b);
		Vertex *otherOther = theElementGraph->getVertexPtr(anotherEleVertexTag);
		if (otherOther->getColor() != to) {
		    toBoundary->addVertex(other,false);
		    b=otherSize;
		}
	    }
	}
    }


    // now remove any SP_Constraints that may have been external to 
    // PartitionedDomain and are now internal to toSubdomain
    SP_ConstraintIter &theSPs = myDomain->getSPs();
    SP_Constraint *spPtr;
    while ((spPtr = theSPs()) != 0) {
  int nodeTag = spPtr->getNodeTag();
  for (int i=0; i<nodesSize; i++) {
      if (nodeTag == nodes(i)) {
    int internal = 0;
    for (int j=1; j<=numPartitions; j++) 
        if (j != to) {
      int prime = primes(j);
      if ((*nodePlace)(nodeTag)%prime == 0) {
          internal = 1;
          j = numPartitions+1;
      }
        }
    if (internal == 0) { // add to toSubdomain if inteernal
        myDomain->removeSP_Constraint(spPtr->getTag());
        toSubdomain->addSP_Constraint(spPtr);
    }
      }
  }
    }

    
    
    // now remove any SP_Constraints that may have been internal to fromSubdomain
    // and are now external to PartitionedDomain or internal to toSubdomain 
    SP_ConstraintIter &theSPs2 = fromSubdomain->getSPs();
    while ((spPtr = theSPs2()) != 0) {
  int nodeTag = spPtr->getNodeTag();
  for (int i=0; i<nodesSize; i++) {
      if (nodeTag == nodes(i)) {
    fromSubdomain->removeSP_Constraint(spPtr->getTag());
    int internal = 0;
    for (int j=1; j<=numPartitions; j++) 
        if (j != to) {
      int prime = primes(j);
      if ((*nodePlace)(nodeTag)%prime == 0) {
          internal = 1;
          j = numPartitions+1;
      }
        }
    if (internal == 0) 
        toSubdomain->addSP_Constraint(spPtr);
    else
        myDomain->addSP_Constraint(spPtr);
      }
  }
    }



  */
    /*
    // now remove any NodalLoads that may have been external to 
    // PartitionedDomain and are now internal to toSubdomain
    NodalLoadIter &theNodalLoads = myDomain->getNodalLoads();
    NodalLoad *loadPtr;
    while ((loadPtr = theNodalLoads()) != 0) {
  int nodeTag = loadPtr->getNodeTag();
  for (int i=0; i<nodesSize; i++) {
      if (nodeTag == nodes(i)) {
    int internal = 0;
    for (int j=1; j<=numPartitions; j++) 
        if (j != to) {
      int prime = primes(j);
      if ((*nodePlace)(nodeTag)%prime == 0) {
          internal = 1;
          j = numPartitions+1;
      }
        }
    if (internal == 0) { // add to toSubdomain if inteernal
        myDomain->removeNodalLoad(loadPtr->getTag());
        toSubdomain->addNodalLoad(loadPtr);
    }
      }
  }
    }


    // now remove any NodalLoads that may have been internal to fromSubdomain
    // and are now external to PartitionedDomain or internal to toSubdomain 

    NodalLoadIter &theNodalLoads2 = myDomain->getNodalLoads();
    while ((loadPtr = theNodalLoads2()) != 0) {
  int nodeTag = loadPtr->getNodeTag();
  for (int i=0; i<nodesSize; i++) {
      if (nodeTag == nodes(i)) {
    fromSubdomain->removeNodalLoad(loadPtr->getTag());
    int internal = 0;
    for (int j=1; j<=numPartitions; j++) 
        if (j != to) {
      int prime = primes(j);
      if ((*nodePlace)(nodeTag)%prime == 0) {
          internal = 1;
          j = numPartitions+1;
      }
        }
    if (internal == 0) 
        toSubdomain->addNodalLoad(loadPtr);
    else
        myDomain->addNodalLoad(loadPtr);
      }  
  }
    }
    */

    return 0;

}



// method to move from from to to, all elements on the interface of 
// from that are adjacent with to.

int 
DomainPartitioner::swapBoundary(int from, int to, bool adjacentVertexNotInOther)
          
{
  opserr << "DomainPartitioner::swapBoundary: from "  << from << " to " << to << endln;
  /***********************************************************************************
    opserr << "DomainPartitioner::swapBoundary from " << from << "  to " << to << endln;
    Timer timer;
    timer.start();

    // check that the object did the partitioning
    if (partitionFlag == false) {
  opserr << "DomainPartitioner::balance(const Vector &load)";
  opserr << " - not partitioned or DomainPartitioner did not partition\n";
  return -1;
    }

    // check that the subdomains exist in partitioned domain
    Subdomain *fromSubdomain = myDomain->getSubdomainPtr(from);
    if (fromSubdomain == 0) {
  opserr << "DomainPartitioner::swapBoundary - No from Subdomain: ";
  opserr << from << " exists\n";
  return -2;
    }

    Subdomain *toSubdomain = myDomain->getSubdomainPtr(to);
    if (toSubdomain == 0) {
  opserr << "DomainPartitioner::swapBoundary - No to Subdomain: ";
  opserr << to << " exists\n";
  return -3;
    }    

    // get the element graph
    //    Graph &theEleGraph = myDomain->getElementGraph();    
    Graph *fromBoundary = theBoundaryElements[from-1];
    Graph *toBoundary = theBoundaryElements[to-1];

    // into a new graph place the vertices that are to be swapped
    Graph *swapVertices = new Graph(fromBoundary->getNumVertex());

    VertexIter &swappableVertices = fromBoundary->getVertices();
    Vertex *vertexPtr;

    while ((vertexPtr = swappableVertices()) != 0) {
	if (adjacentVertexNotInOther == false) {   
	    const ID &adjacency=vertexPtr->getAdjacency();
	    int size = adjacency.Size();
	    for (int i=0; i<size; i++) {
		int otherTag = adjacency(i);
		Vertex *otherVertex = toBoundary->getVertexPtr(otherTag);
		if (otherVertex != 0) {
		    swapVertices->addVertex(vertexPtr,false);
		    i = size;
		}
	    }
	}
	else {
	    const ID &adjacent = vertexPtr->getAdjacency();
	    bool inTo = false;
	    bool inOther = false;
	    int adjacentSize = adjacent.Size();
	    for (int i=0; i<adjacentSize; i++) {
		Vertex *other = theElementGraph->getVertexPtr(adjacent(i));
		if (other->getColor() == to) 
		    inTo = true;
		else if (other->getColor() != from) {
		    inOther = true;
		    i = adjacentSize;
		}
	    }
	    if (inTo == true && inOther == false) { // we remove the vertex
		swapVertices->addVertex(vertexPtr,false);
		vertexPtr->setColor(to);
	    }
	}
    }

    //
    // in the FROM subdomain we:
    //  1. remove the elements corresponding to the vertices
    //  2. remove all nodes connected to the elements divide color by prime
    //  3. see if have to add nodes back as external
    //  4. remove from the PartitionedDomain if node was external
    //  5. add new elements to the boundary vertices for from

    // 1. first remove all the elements from from subdomain
    VertexIter &verticesToSwap = swapVertices->getVertices();

    int numEleToSwap = swapVertices->getNumVertex();
    Element **elementsArray = new Element *[numEleToSwap];

    int count =0; 
    int numCannotRemove =0;
    for (int i=0; i<numEleToSwap; i++) {
      vertexPtr = verticesToSwap();
  if (vertexPtr == 0) 
      return -1;

  int vertexTag = vertexPtr->getTag();
  int eleTag = vertexPtr->getRef();

  // remove the vertex from fromBoundary, set color to to
  // and add to toBoundary
  bool inFrom = true;
  vertexPtr = fromBoundary->removeVertex(vertexTag,false);


  if (inFrom == true) {
      vertexPtr->setColor(to);
      toBoundary->addVertex(vertexPtr,false);
      Element *elePtr = fromSubdomain->removeElement(eleTag);
      if (elePtr == 0) // if ele not there we can just return IT SHOULD BE
    return -4;
      elementsArray[count] = elePtr;
      count++;
  }
  else  // we cannot remove the element
      numCannotRemove++;
    }
    numEleToSwap -= numCannotRemove;

    //  2. remove all nodes connected to the elements divide color by prime    
    //     for each element and add back if still in subdomain
    //  3. add back if still in subdomain
    //  4. if node external remove from PartitionedDomains external nodes

    // first determine which nodes to remove
    int primeFrom = primes(from);
    int primeTo = primes(to);
    int numNodToSwap = 0;
    ID nodesToRemove(0,128);
    for (int j=0; j<numEleToSwap; j++) {
  Element *elePtr = elementsArray[j];
  const ID  &nodes = elePtr->getExternalNodes(); 
  
  int nodesSize = nodes.Size();

  for (int k=0; k<nodesSize; k++) {
      int nodeTag = nodes(k);
      int loc = nodesToRemove.getLocation(nodeTag);
      if (loc < 0) {
    nodesToRemove[numNodToSwap] = nodeTag;
    numNodToSwap++;
      }
      (*nodePlace)(nodeTag) /= primeFrom;    
      (*nodePlace)(nodeTag) *= primeTo;        
  }
    }

    // remove the nodes
    Node **nodesArray = new Node *[numNodToSwap];  
    for (int k=0; k<numNodToSwap; k++) {
  int nodeTag = nodesToRemove(k); 
  Node *nodePtr = fromSubdomain->removeNode(nodeTag);
  nodesArray[k] = nodePtr;

  // add the node back to from if still external
  if ((*nodePlace)(nodeTag)%primeFrom == 0)
      fromSubdomain->addExternalNode(nodePtr);

  // if it was external remove from PartitionedDomains external nodes
  int partition = (*nodePlace)(nodeTag);
  for (int j=1; j<=numPartitions; j++) 
      if (j != from) {
    int prime = primes(j);
    if (partition%prime == 0) {
        myDomain->removeExternalNode(nodeTag);
        j = numPartitions+1;
    }
      }
    }

    // 5. add new vertices to boundary vertices of from if connected to
    //    vertex we just removed and of color of from and not already in boundary

    VertexIter &verticesToSwap2 = swapVertices->getVertices();
    while ((vertexPtr = verticesToSwap2()) != 0) {
  const ID &vertexAdjacent = vertexPtr->getAdjacency();
  int vertexAdjacentSize = vertexAdjacent.Size();

	for (int a=0; a<vertexAdjacentSize; a++) {
	    int otherEleVertexTag = vertexAdjacent(a);
	    Vertex *other = fromBoundary->getVertexPtr(otherEleVertexTag);
	    if (other == 0) {
		other = theElementGraph->getVertexPtr(otherEleVertexTag);
		if (other->getColor() == from)
		    fromBoundary->addVertex(other,false);
	    }
	}	
    }


    // in the TO subdomain we:
    //  1. remove all nodes connected to the element, may or may not be there
    //     if there don't divide by primeTo. 
    //  2. mark all connecting nodes of elements as belonging to To
    //  3. add the nodes, checking if external or internal
    //  3. add the element
    //  4. change the vertex color to be to and add vertex to boundary
    //     also see if we can reduce the size of boundary
    //     of boundary of to vertices

    //  1. remove all nodes connected to the element, may or may not be there
    //      if there divide by primeTo, do this as may have been added as extnl

    for (int l=0; l<numNodToSwap; l++) {    
  int nodeTag = nodesArray[l]->getTag();
  Node *nodePtr;
  nodePtr = toSubdomain->removeNode(nodesToRemove(l));
    }


    for (int m=0; m<numNodToSwap; m++) {
  Node *nodePtr = nodesArray[m];
  int nodeTag = nodePtr->getTag();
  
  //  3. add the nodes, checking if external or internal  
  int internal = 0;
  for (int n=1; n<=numPartitions; n++) 
      if (n != to) {
    int prime = primes(n);
    if ((*nodePlace)(nodeTag)%prime == 0) { // its external
        internal = 1;
        n = numPartitions+1;
    }
      }
  
  if (internal == 0)
      toSubdomain->addNode(nodePtr);
  else {
      toSubdomain->addExternalNode(nodePtr);
      myDomain->addNode(nodePtr);
  }
    }


    delete [] nodesArray; // DELETE nodesArray here should not use nodePtrs
                          // after here as node objects deleted if parallel    
    
    //  3. add the elements
    for (int a=0; a<numEleToSwap; a++)
  toSubdomain->addElement(elementsArray[a]);


    //  4. change the vertex color to be to and add vertex to boundary
    //     also see if we can reduce the size of to's boundary vertices

    VertexIter &verticesToSwap3 = swapVertices->getVertices();
    while ((vertexPtr = verticesToSwap3()) != 0) {
	const ID &vertexAdjacent = vertexPtr->getAdjacency();
	int vertexAdjacentSize = vertexAdjacent.Size();
	for (int n=0; n<vertexAdjacentSize; n++) {
	    int otherEleVertexTag = vertexAdjacent(n);
	    Vertex *other = toBoundary->removeVertex(otherEleVertexTag,false);
	    if (other != 0) {
		const ID &othersAdjacency = other->getAdjacency();
		int otherSize = othersAdjacency.Size();
		for (int b=0; b<otherSize; b++) {
		    int anotherEleVertexTag = othersAdjacency(b);
		    Vertex *otherOther 
			= theElementGraph->getVertexPtr(anotherEleVertexTag);
		    if (otherOther->getColor() != to) {
			toBoundary->addVertex(other,false);
			b=otherSize;
		    }
		}
	    }
	}
    }

    // now remove any SP_Constraints that may have been external to 
    // PartitionedDomain and are now internal to toSubdomain
    SP_ConstraintIter &theSPs = myDomain->getSPs();
    SP_Constraint *spPtr;
    while ((spPtr = theSPs()) != 0) {
  int nodeTag = spPtr->getNodeTag();
  int loc = nodesToRemove.getLocation(nodeTag);
  if (loc >= 0) {
      int internal = 0;
      for (int j=1; j<=numPartitions; j++) 
    if (j != to) {
        int prime = primes(j);
        if ((*nodePlace)(nodeTag)%prime == 0) {
      internal = 1;
      j = numPartitions+1;
        }
    }
      if (internal == 0) { // add to toSubdomain if inteernal
    myDomain->removeSP_Constraint(spPtr->getTag());
    toSubdomain->addSP_Constraint(spPtr);
      }
  }
    }

    
    // now remove any SP_Constraints that may have been internal to fromSubdomain
    // and are now external to PartitionedDomain or internal to toSubdomain 
    SP_ConstraintIter &theSPs2 = fromSubdomain->getSPs();
    while ((spPtr = theSPs2()) != 0) {
  int nodeTag = spPtr->getNodeTag();
  int loc = nodesToRemove.getLocation(nodeTag);
  if (loc >= 0) {
      fromSubdomain->removeSP_Constraint(spPtr->getTag());
      int internal = 0;
      for (int j=1; j<=numPartitions; j++) 
    if (j != to) {
        int prime = primes(j);
        if ((*nodePlace)(nodeTag)%prime == 0) {
      internal = 1;
      j = numPartitions+1;
        }
    }
      if (internal == 0) 
    toSubdomain->addSP_Constraint(spPtr);
      else
    myDomain->addSP_Constraint(spPtr);
  }  
    }

    // now remove any NodalLoads that may have been external to 
    // PartitionedDomain and are now internal to toSubdomain
    NodalLoadIter &theNodalLoads = myDomain->getNodalLoads();
    NodalLoad *loadPtr;
    while ((loadPtr = theNodalLoads()) != 0) {
  int nodeTag = loadPtr->getNodeTag();
  int loc = nodesToRemove.getLocation(nodeTag);
  if (loc >= 0) {

      int internal = 0;
      for (int j=1; j<=numPartitions; j++) 
    if (j != to) {
        int prime = primes(j);
        if ((*nodePlace)(nodeTag)%prime == 0) {
      internal = 1;
      j = numPartitions+1;
        }
    }
      if (internal == 0) { // add to toSubdomain if inteernal
    myDomain->removeNodalLoad(loadPtr->getTag());
    toSubdomain->addNodalLoad(loadPtr);
      }
  }
    }


    // now remove any NodalLoads that may have been internal to fromSubdomain
    // and are now external to PartitionedDomain or internal to toSubdomain 

    NodalLoadIter &theNodalLoads2 = myDomain->getNodalLoads();
    while ((loadPtr = theNodalLoads2()) != 0) {
  int nodeTag = loadPtr->getNodeTag();
  int loc = nodesToRemove.getLocation(nodeTag);
  if (loc >= 0) {
      fromSubdomain->removeNodalLoad(loadPtr->getTag());
      int internal = 0;
      for (int j=1; j<=numPartitions; j++) 
    if (j != to) {
        int prime = primes(j);
        if ((*nodePlace)(nodeTag)%prime == 0) {
      internal = 1;
      j = numPartitions+1;
        }
    }
      if (internal == 0) 
    toSubdomain->addNodalLoad(loadPtr);
      else
    myDomain->addNodalLoad(loadPtr);
  }  
    }


    delete swapVertices;

    timer.pause();
    opserr << "DomainPartitioner::swapBoundary DONE" << timer.getReal() << endln;

  ***************************************************************************/
    return 0;
}


int 
DomainPartitioner::releaseVertex(int from, 
				 int vertexTag, 
				 Graph &theWeightedPartitionGraph,
				 bool mustReleaseToLighter,
				 double factorGreater,
				 bool adjacentVertexNotInOther)
{
  // check that the object did the partitioning
  if (partitionFlag == false) {
    opserr << "DomainPartitioner::balance(const Vector &load)";
    opserr << " - not partitioned or DomainPartitioner did not partition\n";
    return -1;
  }
  
  // we first check the vertex is on the fromBoundary
  Subdomain *fromSubdomain = myDomain->getSubdomainPtr(from);
  if (fromSubdomain == 0) {
    opserr << "DomainPartitioner::swapVertex - No from Subdomain: ";
    opserr << from << " exists\n";
    return -1;
  }

  // get the vertex from the boundary vertices of from
  // Graph &theEleGraph = myDomain->getElementGraph();    
  Graph *fromBoundary = theBoundaryElements[from-1];
    
  Vertex *vertexPtr = fromBoundary->getVertexPtr(vertexTag);
  if (vertexPtr == 0) 
    vertexPtr = theElementGraph->getVertexPtr(vertexTag);
  
  if (vertexPtr == 0)  // if still 0 no vertex given by tag exists
    return -3;

  ID attraction(numPartitions+1);
  attraction.Zero();

  // determine the attraction to the other partitions
  const ID &adjacent = vertexPtr->getAdjacency();
  int numAdjacent = adjacent.Size();
  for (int i=0; i<numAdjacent; i++) {
    int otherTag = adjacent(i);
    Vertex *otherVertex = theElementGraph->getVertexPtr(otherTag);
    int otherPartition = otherVertex->getColor();
    if (otherPartition != from)
      attraction(otherPartition) += 1;
  }
  
  // determine the other partition the vertex is most attracted to
  int partition = 1;
  int maxAttraction = attraction(1);
  for (int j=2; j<=numPartitions; j++)
    if (attraction(j) > maxAttraction) {
      partition = j;
      maxAttraction = attraction(j);
    }

  // swap the vertex
  if (mustReleaseToLighter == false)
    return swapVertex(from, partition, vertexTag, adjacentVertexNotInOther);
  
  else { // check the other partition has a lighter load
    Vertex *fromVertex = theWeightedPartitionGraph.getVertexPtr(from);
    Vertex *toVertex = theWeightedPartitionGraph.getVertexPtr(partition);	    
    
    double fromWeight = fromVertex->getWeight();
    double toWeight  = toVertex->getWeight();
    
    if (fromWeight == toWeight)
      opserr << "DomainPartitioner::releaseVertex - TO CHANGE >= to >\n";

    if (fromWeight >= toWeight) {
      if (toWeight == 0.0) 
	return swapVertex(from,partition,vertexTag,adjacentVertexNotInOther);	    
      if (fromWeight/toWeight > factorGreater)        
	return swapVertex(from,partition,vertexTag,adjacentVertexNotInOther);	    
    }
  }
  
  return 0;
}



int 
DomainPartitioner::releaseBoundary(int from, 
           Graph &theWeightedPartitionGraph,  
           bool mustReleaseToLighter,
           double factorGreater,
           bool adjacentVertexNotInOther)
{
    // check that the object did the partitioning
    if (partitionFlag == false) {
      opserr << "DomainPartitioner::balance(const Vector &load)";
      opserr << " - not partitioned or DomainPartitioner did not partition\n";
      return -1;
    }

    // we first get a pointer to fromSubdomain & the fromBoundary
    Subdomain *fromSubdomain = myDomain->getSubdomainPtr(from);
    if (fromSubdomain == 0) {
      opserr << "DomainPartitioner::swapVertex - No from Subdomain: ";
      opserr << from << " exists\n";
      return -1;
    }

    //    Graph &theEleGraph = myDomain->getElementGraph();    
    Graph *fromBoundary = theBoundaryElements[from-1];

    // into a new graph place the vertices on the fromBoundary
    // we cannot use fromBoundary as this would empty all the nodes
    // as fromBoundary changes in called methods

    Graph *swapVertices = new Graph(fromBoundary->getNumVertex());

    VertexIter &swappableVertices = fromBoundary->getVertices();
    Vertex *vertexPtr;

    while ((vertexPtr = swappableVertices()) != 0) 
      swapVertices->addVertex(vertexPtr,false);

    // release all the vertices in the swapVertices
    VertexIter &verticesToSwap = swapVertices->getVertices();
    while ((vertexPtr = verticesToSwap()) != 0)
      releaseVertex(from,
		    vertexPtr->getTag(),
		    theWeightedPartitionGraph,
		    mustReleaseToLighter,
		    factorGreater,
		    adjacentVertexNotInOther);
    
    delete swapVertices;    

    return 0;
}

