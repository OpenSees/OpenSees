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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:00:59 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/partitioner/DomainPartitioner.cpp,v $
                                                                        
                                                                        
 // File: ~/domain/partitioner/DomainPartitioner.C
// 
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
 
DomainPartitioner::DomainPartitioner(GraphPartitioner &theGraphPartitioner)
:myDomain(0),thePartitioner(theGraphPartitioner),theBalancer(0),
 theBoundaryElements(0), nodePlace(0),elementPlace(0), 
 numPartitions(0), primes(10), partitionFlag(false)
{
    // set primes to handle up to 8 partitions, can enlarge if needed
    primes(1) = 2; primes(2) = 3; primes(3) = 5; primes(4) = 7; primes(5) = 11;
    primes(6) = 13; primes(7) = 17; primes(8) = 23; primes(9) = 31;

    opserr << "DomainPartitioner::DomainPartitioner - ";
    opserr << "does not deal with ele loads or mp_constraints yet\n";
}    

DomainPartitioner::DomainPartitioner(GraphPartitioner &theGraphPartitioner,
				     LoadBalancer &theLoadBalancer)
:myDomain(0),thePartitioner(theGraphPartitioner),theBalancer(&theLoadBalancer),
 theBoundaryElements(0),
 nodePlace(0),elementPlace(0), numPartitions(0),primes(10), partitionFlag(false)
{
    // set primes to handle up to 8 partitions, can enlarge if needed
    primes(1) = 2; primes(2) = 3; primes(3) = 5; primes(4) = 7; primes(5) = 11;
    primes(6) = 13; primes(7) = 17; primes(8) = 23; primes(9) = 31;

    
    // set the links the loadBalancer needs
    theLoadBalancer.setLinks(*this);
}    


DomainPartitioner::~DomainPartitioner()
{

}

void 
DomainPartitioner::setPartitionedDomain(PartitionedDomain &theDomain)
{
    myDomain = &theDomain;
}

int
DomainPartitioner::partition(int numParts)
{

    // first we ensure the partitioned domain has numpart subdomains
    // with tags 1 through numparts
    for (int i=1; i<=numParts; i++) {
	Subdomain *subdomainPtr = myDomain->getSubdomainPtr(i);
	if (subdomainPtr == 0) {
	    opserr << "DomainPartitioner::partition - No Subdomain: ";
	    opserr << i << " exists\n";
	    return -1;
	}
    }

    // Timer timer;
    // timer.start();

    // we get the ele graph from the domain and partition it
    Graph &theEleGraph = myDomain->getElementGraph();
    
    //    timer.pause(); opserr << "partition:get eleGraph " << timer.getReal() << timer.getCPU() << endln;
    
    int theError = thePartitioner.partition(theEleGraph,numParts);
    if (theError < 0) {
	opserr << "DomainPartitioner::partition";
	opserr << " - the graph partioner failed to partition the ";
	opserr << "element graph\n";
	return -10+theError;
    }

    //    timer.pause(); opserr << "partition:partition " << timer.getReal() << timer.getCPU() << endln;
    
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

    // we iterate through the nodes setting all their partition values to 1
    // the partition value will be used to determine what partition the node 
    // will be added to and whether or not the node is a boundary node.

    if (numParts+1 >= primes.Size()) { // colors 1 through numParts 
	opserr << "DomainPartitioner::partition(int numParts)";
	opserr << " - primes is  not yet big enough and function not yet written\n";
	return -1;
    }
    // we now iterate through the nodes to first determine 
    // min and max node numbers

    NodeIter &theNodes = myDomain->getNodes();
    Node *nodePtr;
    int maxNodeTag =0;
    int minNodeTag = 0;

    while ((nodePtr = theNodes()) != 0) {
	int tag = nodePtr->getTag();
	if (tag > maxNodeTag) maxNodeTag = tag;
	if (tag < minNodeTag) minNodeTag = tag;
    }

    if (minNodeTag < 0) {
	    opserr << "DomainPartitioner::partition(int numParts)";
	    opserr << " - minNodeTag < 0 \n";
	    numPartitions = 0;
	    return -1;	
    }

    // we now create an ID which will contain information
    // about which partition the node is located in
    
    nodePlace = new ID(maxNodeTag+1);
    if (nodePlace == 0 || nodePlace->Size() < maxNodeTag+1) {
	opserr << "DomainPartitioner::partition(int numParts)";
	opserr << " - ran out of memory\n";
	numPartitions = 0;
	return -1;
    }
    for (int m=minNodeTag; m<=maxNodeTag; m++)
	(*nodePlace)(m) = 1;

    //
    // we now iterate through the vertices of the element graph
    // to see if the vertex is a boundary vertex or not - if it is
    // we add to the appropriate graph created above. We also set the
    // value the color variable of each of the external nodes connected 
    // to the element to a value which will indicate that that node will
    // have to be added to the subdomain.
    //

    VertexIter &theVertexIter = theEleGraph.getVertices();
    Vertex *vertexPtr;
    while ((vertexPtr = theVertexIter()) != 0) {
	int eleTag = vertexPtr->getRef();
	int vertexColor = vertexPtr->getColor();

	const ID &adjacency = vertexPtr->getAdjacency();
	int size = adjacency.Size();
	for (int i=0; i<size; i++) {
	    Vertex *otherVertex = theEleGraph.getVertexPtr(adjacency(i));
	    if (otherVertex->getColor() != vertexColor) {
		theBoundaryElements[vertexColor-1]->addVertex(vertexPtr,false);
		i = size;
	    }
	}

	Element *elePtr = myDomain->getElement(eleTag);

	const ID &nodes = elePtr->getExternalNodes();
	size = nodes.Size();
	int primeVertexColor = primes(vertexColor);
	for (int j=0; j<size; j++) {
	    int node = nodes(j);
	    (*nodePlace)(node) *= primeVertexColor;
	}
    }
    
    // now go through the MP_Constraints and ensure:
    // 	1. if in different partitions both on boundary
    //  2. if constrained on boundary - retained also on boundary
    MP_ConstraintIter &theMPs = myDomain->getMPs();
    MP_Constraint *mpPtr;
    while ((mpPtr = theMPs()) != 0) {
	int retained = mpPtr->getNodeRetained();
	int constrained = mpPtr->getNodeConstrained();
	
	if ((*nodePlace)(retained) != (*nodePlace)(constrained)) { 
	    // we could have a problem 1. or 2. above
	    
	    // for now just put both nodes on boundary
	    (*nodePlace)(retained) *= primes(numParts+1);
	    (*nodePlace)(constrained) *= primes(numParts+1);
        }
    }

    //    timer.pause(); opserr << "partition:figure out nodes " << timer.getReal() << timer.getCPU() << endln;
    
    // we now add the nodes, by iterating through the nodes and for
    // each node determining if it is an internal or external node
    // to each subdomain
    NodeIter &theNodes2 = myDomain->getNodes();
    while ((nodePtr = theNodes2()) != 0) {
	int nodeTag = nodePtr->getTag();
	int partition = (*nodePlace)(nodeTag);
	
	for (int j=1; j<=numParts; j++) {

	    int prime = primes(j);
	    // add the nodes
	    
	    if (partition%prime == 0) { // it belongs to partition i
		
		// determine if internal or external and add accordingly
		int internal = 0; // assume internal
		for (int k=1; k<= numParts+1; k++) {  // numParts+1 for MP_Constraints
		    if (k != j) {
			prime = primes(k);
			if (partition%prime == 0) { // its external
			    internal = 1;
			    k = numParts+1;
			}
		    }
		}
		
		if (internal == 0) {  // its an internal node
		    Subdomain *theSubdomain = myDomain->getSubdomainPtr(j); 
		    nodePtr = myDomain->removeNode(nodePtr->getTag());
		    theSubdomain->addNode(nodePtr);
		    j = numParts+1;
		}
		else {  // its external, dont remove from partitioned domain
		    Subdomain *theSubdomain = myDomain->getSubdomainPtr(j); 
		    theSubdomain->addExternalNode(nodePtr);
		}
	    }
	}
    }


    // now we go through the load patterns and move NodalLoads and MP_Constraints
    // 1) make sure each subdomain has a copy of the partitioneddomains load patterns.
    // 2) move nodal loads
    // 3) move SP_Constraints
    
    LoadPatternIter &theLoadPatterns = myDomain->getLoadPatterns();
    LoadPattern *theLoadPattern;
    while ((theLoadPattern = theLoadPatterns()) != 0) {
      int loadPatternTag = theLoadPattern->getTag();

      // check that each subdomain has a loadPattern with a similar tag and class tag
      for (int i=1; i<=numParts; i++) {
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

      // now remove any nodal loads that correspond to internal nodes in a subdomain
      // and add them to the appropriate loadpattern in the subdomain

      NodalLoadIter &theNodalLoads = theLoadPattern->getNodalLoads();
      NodalLoad *theNodalLoad;
      while ((theNodalLoad = theNodalLoads()) != 0) {
	int nodeTag = theNodalLoad->getNodeTag();
	int partition = (*nodePlace)(nodeTag);

	for (int i=1; i<=numParts; i++) {
	  int prime = primes(i);
	  
	  if (partition % prime == 0) { // node in partition
	    // now test if internal, if so add
	    int internal = 0;
	    for (int j=1; j<=numParts; j++) 
	      if (i != j) {
		prime = primes(j);
		if (partition%prime == 0) { // its external
		  internal = 1;
		  j = numParts+1;
		}
	      }
	    if (internal == 0) {		
	      theLoadPattern->removeNodalLoad(theNodalLoad->getTag());
	      Subdomain *theSubdomain = myDomain->getSubdomainPtr(i);
	      theSubdomain->addNodalLoad(theNodalLoad, loadPatternTag);
	    }
	    else // its a boundary node
	      i = numParts+1;
	  }
	}
      }

      SP_ConstraintIter &theSPs = theLoadPattern->getSPs();
      SP_Constraint *spPtr;
      while ((spPtr = theSPs()) != 0) {
	int partition = (*nodePlace)(spPtr->getNodeTag());
	for (int i=1; i<=numParts; i++) {
	  int prime = primes(i);
	  if (partition%prime == 0) { // it belongs to partition
	    
	    // now test if internal, if so add
	    int internal = 0;
	    for (int j=1; j<=numParts; j++) 
	      if (i != j) {
		prime = primes(j);
		if (partition%prime == 0) { // its external
		  internal = 1;
		  j = numParts+1;
		}
	      }
	    if (internal == 0) {
	      Subdomain *theSubdomain = myDomain->getSubdomainPtr(i);	
	      theLoadPattern->removeSP_Constraint(spPtr->getTag());
	      theSubdomain->addSP_Constraint(spPtr, loadPatternTag);
	      i = numParts+1;
	    }
	  }
	}
      }

    }
      

    //    timer.pause(); opserr << "partition:add nodes " << timer.getReal() << timer.getCPU() << endln;
    
    // we now move the elements and any elemental Loads in the loadPatterns
    VertexIter &theVertices = theEleGraph.getVertices();
    while ((vertexPtr = theVertices()) != 0) {
      // move the element
      int partition = vertexPtr->getColor();
      int eleTag = vertexPtr->getRef();
      Element *elePtr = myDomain->removeElement(eleTag);	
      Subdomain *theSubdomain = myDomain->getSubdomainPtr(partition);	
      theSubdomain->addElement(elePtr);

      // if any corresponding elemental loads in the load patterns .. move the load as well
      LoadPatternIter &theLoadPatterns = myDomain->getLoadPatterns();
      LoadPattern *theLoadPattern;
      while ((theLoadPattern = theLoadPatterns()) != 0) {
	int loadPatternTag = theLoadPattern->getTag();
	ElementalLoadIter &theLoads = theLoadPattern->getElementalLoads();
	ElementalLoad *theLoad;
	while ((theLoad = theLoads()) != 0) {
	  opserr << "DomainPArtitioner::partition - REMOVE ELEMENTAL LOADS\n";
	  /*
	  if (theLoad->getElementTag() == eleTag)
	    theLoadPattern->removeElementalLoad(theLoad->getTag());
	  theSubdomain->addElementalLoad(theLoad, loadPatternTag);
	  */
	}
      }
    }
    
    // timer.pause(); opserr << "partition:add elements " << timer.getReal() << timer.getCPU() << endln;

    // add the single point constraints, only added if for an internal node in a subdomain

    SP_ConstraintIter &theSPs = myDomain->getSPs();
    SP_Constraint *spPtr;
    while ((spPtr = theSPs()) != 0) {
	int partition = (*nodePlace)(spPtr->getNodeTag());
	for (int i=1; i<=numParts; i++) {
	    int prime = primes(i);
	    if (partition%prime == 0) { // it belongs to partition
		
		// now test if internal, if so add
		int internal = 0;
		for (int j=1; j<=numParts; j++) 
		    if (i != j) {
			prime = primes(j);
			if (partition%prime == 0) { // its external
			    internal = 1;
			    j = numParts+1;
			}
		    }
		if (internal == 0) {
		    Subdomain *theSubdomain = myDomain->getSubdomainPtr(i);	
		    myDomain->removeSP_Constraint(spPtr->getTag());
		    theSubdomain->addSP_Constraint(spPtr);
		    i = numParts+1;
		}
	    }
	}
    }

 /******************************** all MP_Constraints are external 
    // add the MP_Constraints, only added if both internal to a subdomain
    MP_ConstraintIter &moreMPs = myDomain->getMPs();
    while ((mpPtr = moreMPs()) != 0) {
	int retained = mpPtr->getNodeRetained();
	int constrained = mpPtr->getNodeConstrained();
	int partRet = (*nodePlace)(retained);
	int partCon = (*nodePlace)(constrained);
	if (partRet == partCon) {
	    for (int i=1; i<=numParts; i++) {
		int prime = primes(i);
		if (partRet%prime == 0) { // it belongs to partition
		
		    // now test if internal, if so add
		    int internal = 0;
		    for (int j=1; j<=numParts; j++) 
			if (i != j) {
			    prime = primes(j);
			    if (partRet%prime == 0) { // its external
				internal = 1;
				j = numParts+1;
			    }
			}
		    if (internal == 0) {		
			Subdomain *theSubdomain = myDomain->getSubdomainPtr(i);
			myDomain->removeMP_Constraint(mpPtr->getTag());
			theSubdomain->addMP_Constraint(mpPtr);
			i = numParts+1;
		    }
		}
	    }	
	}
    }
    **************************************************************/
    
    // now we go through all the subdomains and tell them to update
    // their analysis for the new layouts

    SubdomainIter &theSubDomains = myDomain->getSubdomains();
    Subdomain *theSubDomain;
    while ((theSubDomain = theSubDomains()) != 0) 
	if (theSubDomain->hasDomainChanged() == true) 
	    theSubDomain->invokeChangeOnAnalysis();


    // we invoke change on the PartitionedDomain
    /***** it's up to the analysis object to determine if domain has changed
    if (myDomain->hasDomainChanged() == true) 
	myDomain->invokeChangeOnAnalysis();
    ****************************************************************/

    //    timer.pause(); opserr << "partition:domain change " << timer.getReal() << timer.getCPU() << endln;

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

	Timer theTimer;
	theTimer.start();
	
	// call on the LoadBalancer to partition		
	res = theBalancer->balance(theWeightedPGraph);
	theTimer.pause();
	opserr << "DomainPartitioner::balance 1 real: " << theTimer.getReal() << endln;
	    
	// now invoke domainChanged on Subdomains and PartitionedDomain
	SubdomainIter &theSubDomains = myDomain->getSubdomains();
	Subdomain *theSubDomain;
	while ((theSubDomain = theSubDomains()) != 0) 
	    if (theSubDomain->hasDomainChanged() == true) 
		theSubDomain->invokeChangeOnAnalysis();

	theTimer.pause();
	opserr << "DomainPartitioner::balance 1 real: " << theTimer.getReal() << endln;

	/************ up to analysis to determine if domain has changed
	// we invoke change on the PartitionedDomain
	if (myDomain->hasDomainChanged() == true) 
	    myDomain->invokeChangeOnAnalysis();	    
	    ******************/
	theTimer.pause();
	opserr << "DomainPartitioner::balance real: " << theTimer.getReal() << endln;
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
  /**********************************************************************
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

    Graph &theEleGraph = myDomain->getElementGraph();    
    Graph *fromBoundary = theBoundaryElements[from-1];
    Graph *toBoundary = theBoundaryElements[to-1];    
    Vertex *vertexPtr;

    // get a pointer to the vertex in the element graph
    if (adjacentVertexNotInOther  == false) {
	vertexPtr = fromBoundary->removeVertex(vertexTag,false);    
	if (vertexPtr == 0) 
	    vertexPtr = theEleGraph.getVertexPtr(vertexTag);
    
	if (vertexPtr == 0)  // if still 0 no vertex given by tag exists
	    return -4;
    } else { // get a pointer and check vertex not adjacent to another partition
	vertexPtr = fromBoundary->getVertexPtr(vertexTag);    
	if (vertexPtr == 0) 
	    vertexPtr = theEleGraph.getVertexPtr(vertexTag);
	if (vertexPtr == 0)  // if still 0 no vertex given by tag exists
	    return -4;	
	const ID &adjacent = vertexPtr->getAdjacency();
	bool inTo = false;
        bool inOther = false;
	int adjacentSize = adjacent.Size();
	for (int i=0; i<adjacentSize; i++) {
	    Vertex *other = theEleGraph.getVertexPtr(adjacent(i));
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
    //	1. remove the element
    //	2. remove all nodes connected to the element divide color by prime
    //  3. see if have to add nodes back as external
    //  4. remove from the PartitionedDomain if node was external
    //  5. add new elements to the boundary vertices for from
    
    //  1. remove the element from the fromSubdomain
    //          this gives us a pointer to the element as well.

    Element *elePtr = fromSubdomain->removeElement(eleTag);
    if (elePtr == 0) // if ele not there we can just return but ERROR it should be
	return -6;

    //	2. remove all nodes connected to the element divide color by prime    
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
	    other = theEleGraph.getVertexPtr(otherEleVertexTag);
	    if (other->getColor() == from)
		fromBoundary->addVertex(other,false);
	}
    }


    
    // in the TO subdomain we:
    //	1. remove all nodes connected to the element, may or may not be there
    //	   if there don't divide by primeTo. 
    //	2. mark all connecting nodes of elements as belonging to To
    //	3. add the nodes, checking if external or internal
    //  3. add the element
    //  4. change the vertex color to be to and add vertex to boundary
    //     also see if we can reduce the size of boundary
    //     of boundary of to vertices

    //	1. remove all nodes connected to the element, may or may not be there
    //	    if there divide by primeTo, do this as may have been added as extnl
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
		Vertex *otherOther = theEleGraph.getVertexPtr(anotherEleVertexTag);
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

  **********************************************************************/
    return 0;

}



// method to move from from to to, all elements on the interface of 
// from that are adjacent with to.

int 
DomainPartitioner::swapBoundary(int from, int to, bool adjacentVertexNotInOther)
			    
{
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
    Graph &theEleGraph = myDomain->getElementGraph();    
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
		Vertex *other = theEleGraph.getVertexPtr(adjacent(i));
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
    //	1. remove the elements corresponding to the vertices
    //	2. remove all nodes connected to the elements divide color by prime
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

    //	2. remove all nodes connected to the elements divide color by prime    
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
		other = theEleGraph.getVertexPtr(otherEleVertexTag);
		if (other->getColor() == from)
		    fromBoundary->addVertex(other,false);
	    }
	}	
    }


    // in the TO subdomain we:
    //	1. remove all nodes connected to the element, may or may not be there
    //	   if there don't divide by primeTo. 
    //	2. mark all connecting nodes of elements as belonging to To
    //	3. add the nodes, checking if external or internal
    //  3. add the element
    //  4. change the vertex color to be to and add vertex to boundary
    //     also see if we can reduce the size of boundary
    //     of boundary of to vertices

    //	1. remove all nodes connected to the element, may or may not be there
    //	    if there divide by primeTo, do this as may have been added as extnl

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
			= theEleGraph.getVertexPtr(anotherEleVertexTag);
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
    Graph &theEleGraph = myDomain->getElementGraph();    
    Graph *fromBoundary = theBoundaryElements[from-1];


    Vertex *vertexPtr = fromBoundary->getVertexPtr(vertexTag);
    if (vertexPtr == 0) 
	vertexPtr = theEleGraph.getVertexPtr(vertexTag);
    
    if (vertexPtr == 0)  // if still 0 no vertex given by tag exists
	return -3;
    
    ID attraction(numPartitions+1);
    attraction.Zero();


    // determine the attraction to the other partitions
    const ID &adjacent = vertexPtr->getAdjacency();
    int numAdjacent = adjacent.Size();
    for (int i=0; i<numAdjacent; i++) {
	int otherTag = adjacent(i);
	Vertex *otherVertex = theEleGraph.getVertexPtr(otherTag);
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
	return swapVertex(from,partition,vertexTag,adjacentVertexNotInOther);

    else { // check the other partition has a lighter load
	Vertex *fromVertex = theWeightedPartitionGraph.getVertexPtr(from);
	Vertex *toVertex = theWeightedPartitionGraph.getVertexPtr(partition);	    

        double fromWeight = fromVertex->getWeight();
        double toWeight  = toVertex->getWeight();

	if (fromWeight > toWeight) {
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
    opserr << "DomainPartitioner::releaseBoundary from " << from << endln;
    Timer timer;
    timer.start();

    
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

    Graph &theEleGraph = myDomain->getElementGraph();    
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

    timer.pause();
    opserr << "DomainPartitioner::releaseBoundary DONE" << timer.getReal() << endln;    

    return 0;
}

