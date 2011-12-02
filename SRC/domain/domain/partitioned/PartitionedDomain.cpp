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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/partitioned/PartitionedDomain.cpp,v $
                                                                        
                                                                        
// File: ~/domain/domain/partitioned/PartitionedDomain.C
// 
// Written: fmk 
// Created: Wed Sep 25 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for PartitionedDomain.
// PartitionedDomain is an abstract class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints just like a normal domain. In addition the domain provides
// a method to partition the domain into Subdomains.
//
// ModelBuilder. There are no partitions in a PartitionedDomain.
//
// What: "@(#) PartitionedDomain.C, revA"

#include <PartitionedDomain.h>
#include <stdlib.h>

#include <DomainPartitioner.h>
#include <Element.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <ArrayOfTaggedObjects.h>
#include <ArrayOfTaggedObjectsIter.h>
#include <Subdomain.h>
#include <DomainPartitioner.h>
#include <PartitionedDomain.h>
#include <PartitionedDomainEleIter.h>
#include <PartitionedDomainSubIter.h>
#include <SingleDomEleIter.h>
#include <Vertex.h>
#include <Graph.h>

PartitionedDomain::PartitionedDomain(DomainPartitioner &thePartitioner)
:Domain(),
 theSubdomains(0),theDomainPartitioner(&thePartitioner),
 theSubdomainIter(0), mySubdomainGraph(0)
{
    elements = new ArrayOfTaggedObjects(1024);    
    theSubdomains = new ArrayOfTaggedObjects(32);
    theSubdomainIter = new PartitionedDomainSubIter(theSubdomains);

    mainEleIter = new SingleDomEleIter(elements);    
    theEleIter = new PartitionedDomainEleIter(this);
    
    if (theSubdomains == 0 || elements == 0 ||
	theSubdomainIter == 0 || theDomainPartitioner == 0 ||
	theEleIter == 0 || mainEleIter == 0) {
	
	cerr << "FATAL: PartitionedDomain::PartitionedDomain ";
	cerr << "  - ran out of memory\n";
	exit(-1);
    }
}


PartitionedDomain::PartitionedDomain(int numNodes, int numElements,
				     int numSPs, int numMPs, int numLoadPatterns,
				     int numSubdomains,
				     DomainPartitioner &thePartitioner)

:Domain(numNodes,0,numSPs,numMPs,numLoadPatterns),
 theSubdomains(0),theDomainPartitioner(&thePartitioner),
 theSubdomainIter(0), mySubdomainGraph(0)
{
    elements = new ArrayOfTaggedObjects(numElements);    
    theSubdomains = new ArrayOfTaggedObjects(numSubdomains);
    theSubdomainIter = new PartitionedDomainSubIter(theSubdomains);

    mainEleIter = new SingleDomEleIter(elements);    
    theEleIter = new PartitionedDomainEleIter(this);
    
    if (theSubdomains == 0 || elements == 0 ||
	theSubdomainIter == 0 || theDomainPartitioner == 0 ||
	theEleIter == 0 || mainEleIter == 0) {
	
	cerr << "FATAL: PartitionedDomain::PartitionedDomain(int ..) ";
	cerr << "  - ran out of memory\n";
	exit(-1);
    }
}




PartitionedDomain::~PartitionedDomain()
{
    if (elements != 0)
	delete elements;

    if (theSubdomains != 0)
	delete theSubdomains;
    
    if (theSubdomainIter != 0)
	delete theSubdomainIter;
    
    if (theEleIter != 0)
	delete theEleIter;

}

void
PartitionedDomain::clearAll(void)
{
  this->Domain::clearAll();

  elements->clearAll();
  theSubdomains->clearAll();
}
    


bool 
PartitionedDomain::addElement(Element *elePtr)
{

  if (elePtr->isSubdomain() == true)
    return this->addSubdomain((Subdomain *)elePtr);

  int eleTag = elePtr->getTag();
#ifdef _DEBUG      
    if (check == true) {

	// check ele Tag >= 0
	if (eleTag < 0) {
	    cerr << "Domain::addElement - Element " << eleTag;
	    cerr << " tag must be >= 0\n";
	    return false;
	}      
	
	// check its not in this or any of the subdomains
	// MISSING CODE	
	
	// check all the elements nodes exist in the domain
	const ID &nodes = elePtr->getExternalNodes();
	for (int i=0; i<nodes.Size(); i++) {
	    int nodeTag = nodes(i);
	    Node *nodePtr = this->getNode(nodeTag);
	    if (nodePtr == 0) {
		cerr << "Domain::addElement - In element " << eleTag;
		cerr << " no node " << nodeTag << " exists in the domain\n";
		return false;
	    }      	
	}
	
    }
#endif
    
    TaggedObject *other = elements->getComponentPtr(eleTag);
    if (other != 0)
	return false;
    
    bool result = elements->addComponent(elePtr);
    if (result == true) {
	elePtr->setDomain(this);
	this->domainChange();
    }
    
    return result;
}    



bool 
PartitionedDomain::addNode(Node *nodePtr)
{
#ifdef _DEBUG    
    if (check == true) {
	// check its not in this or any of the subdomains

	// MISSING CODE	
    }
#endif
    return (this->Domain::addNode(nodePtr));    
}

Element *
PartitionedDomain::removeElement(int tag)
{
    // we first see if its in the original domain
    TaggedObject *res = elements->removeComponent(tag);
    Element *result = 0;
    if (res != 0) {
	result = (Element *)res;
	this->domainChange();
	return result;
    }

    // if not there we must check all the other subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->removeElement(tag);
	    if (result != 0) {
		return result;
	    }
	}
    }

    // its not there
    return 0;
}    


Node *
PartitionedDomain::removeNode(int tag)    
{
    // we first remove it form the original domain (in case on boundary)
    Node *result = this->Domain::removeNode(tag);

    // we must also try removing from the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    Node *res = theSub->removeNode(tag);
	    if (res != 0) 
		result = res;
	}
    }
    
    if (result != 0) 
	this->domainChange();
    
    return result;
}    

SP_Constraint *
PartitionedDomain::removeSP_Constraint(int tag)
{
    // we first see if its in the original domain
    SP_Constraint *result = this->Domain::removeSP_Constraint(tag);
    if (result != 0) {
	this->domainChange();
	return result;
    }
	

    // if not there we must check all the other subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->removeSP_Constraint(tag);
	    if (result != 0) {
		return result;
	    }
	}
    }

    // its not there
    return 0;
}


MP_Constraint *
PartitionedDomain::removeMP_Constraint(int tag)
{
    // we first see if its in the original domain
    MP_Constraint *result = this->Domain::removeMP_Constraint(tag);
    if (result != 0) {
	this->domainChange();
	return result;
    }

    // if not there we must check all the other subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);		
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->removeMP_Constraint(tag);
	    if (result != 0) {
		return result;
	    }
	}
    }    

    // its not there
    return 0;
}









// public member functions which have to be modified
ElementIter       &
PartitionedDomain::getElements()
{
    theEleIter->reset();
    return *theEleIter;
}    


Element  *
PartitionedDomain::getElement(int tag) 
{
    // we first see if its in the original domain
    TaggedObject *res = elements->getComponentPtr(tag);
    Element *result =0;
    if (res != 0) {
	result = (Element *)res;
	return result;
    }

    // go through the other subdomains until we find it or we run out of subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    result = theSub->getElement(tag);
	    if (result != 0) 
		return result;
	}
    }
    
    // its not there
    return 0;
}


int 		
PartitionedDomain::getNumElements(void) const
{
    int result = elements->getNumComponents();

    // add the number of subdomains
    result +=  theSubdomains->getNumComponents();
    return result;
}

void
PartitionedDomain::applyLoad(double timeStep)
{
    this->Domain::applyLoad(timeStep);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->applyLoad(timeStep);
	}
    }
}


void
PartitionedDomain::setCommitTag(int newTag)
{
    this->Domain::setCommitTag(newTag);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->setCommitTag(newTag);
	}
    }
}



void
PartitionedDomain::setCurrentTime(double newTime)
{
    this->Domain::setCurrentTime(newTime);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->setCurrentTime(newTime);
	}
    }
}


void
PartitionedDomain::setCommittedTime(double newTime)
{
    this->Domain::setCommittedTime(newTime);

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->setCommittedTime(newTime);
	}
    }
}


void
PartitionedDomain::setLoadConstant(void)
{
    this->Domain::setLoadConstant();

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->setLoadConstant();
	}
    }
}


int
PartitionedDomain::update(void)
{
    this->Domain::update();

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    theSub->computeNodalResponse();
	    theSub->update();
	}
    }
    return 0;
}

int
PartitionedDomain::commit(void)
{
    int result = this->Domain::commit();
    if (result < 0) {
	cerr << "PartitionedDomain::commit(void) - failed in Domain::commit()\n";
	return result;
    }

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    int res = theSub->commit();
	    if (res < 0) {
		cerr << "PartitionedDomain::commit(void)";
		cerr << " - failed in Subdomain::commit()\n";
		return res;
	    }	    
	}
    }

    // now we load balance if we have subdomains
    int numSubdomains = this->getNumSubdomains();
    if (numSubdomains != 0) 
	theDomainPartitioner->balance(this->getSubdomainGraph());
    
    return 0;
}


int
PartitionedDomain::revertToLastCommit(void)
{
    int result = this->Domain::revertToLastCommit();
    if (result < 0) {
	cerr << "PartitionedDomain::revertToLastCommit(void) - failed in Domain::revertToLastCommit()\n";
	return result;
    }

    // do the same for all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);	
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    Subdomain *theSub = (Subdomain *)theObject;	    
	    int res = theSub->revertToLastCommit();
	    if (res < 0) {
		cerr << "PartitionedDomain::revertToLastCommit(void)";
		cerr << " - failed in Subdomain::revertToLastCommit()\n";
		return res;
	    }	    
	}
    }

    return 0;
}



void 
PartitionedDomain::Print(ostream &s, int flag)
{
    this->Domain::Print(s, flag);

  s << "\nELEMENT DATA: NumEle: " << elements->getNumComponents() << "\n";
  elements->Print(s);

	
    // print all the subdomains
    if (theSubdomains != 0) {
	ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
	TaggedObject *theObject;
	while ((theObject = theSubsIter()) != 0) {
	    theObject->Print(s, flag);
	}
    }
}




int 
PartitionedDomain::partition(int numPartitions)
{
    // need to create element graph before create new subdomains
    // DO NOT REMOVE THIS LINE __ EVEN IF COMPILER WARNING ABOUT UNUSED VARIABLE
    Graph &theEleGraph = this->getElementGraph();

    // check to see that they have ones with the correct tags
    if (theSubdomains != 0) {
	for (int i=1; i<=numPartitions; i++) {
	    TaggedObject *theObject = theSubdomains->getComponentPtr(i);
	    if (theObject == 0) { // create a subdomain with appropriate tag
		cerr << "PartitionedDomain::partition - GET A COPY\n";
		exit(-1);
	    }
	}
    } 
    // now we call partition on the domainPartitioner which does the partitioning
    DomainPartitioner &thePartitioner = this->getPartitioner();
    thePartitioner.setPartitionedDomain(*this);

    return thePartitioner.partition(numPartitions);
}   

bool 
PartitionedDomain::addSubdomain(Subdomain *theSubdomain)
{
    int eleTag = theSubdomain->getTag();
    TaggedObject *other = theSubdomains->getComponentPtr(eleTag);
    if (other != 0)
	return false;
  
    bool result = theSubdomains->addComponent(theSubdomain);
    if (result == true) {
	theSubdomain->setDomain(this);
	this->domainChange();
    }
  
  return result;

}

int 
PartitionedDomain::getNumSubdomains(void)
{
    return theSubdomains->getNumComponents();
}
    
Subdomain *
PartitionedDomain::getSubdomainPtr(int tag)
{
    TaggedObject *mc = theSubdomains->getComponentPtr(tag);
    if (mc == 0) return 0;
    Subdomain *result = (Subdomain *)mc;
    return result;
}    

SubdomainIter &
PartitionedDomain::getSubdomains(void)
{
    theSubdomainIter->reset();
    return *theSubdomainIter;
}



DomainPartitioner &
PartitionedDomain::getPartitioner(void) const
{
    return *theDomainPartitioner;
}
	



int 
PartitionedDomain::buildEleGraph(Graph *theEleGraph)
{
    int numVertex = elements->getNumComponents();

    // see if quick return

    if (numVertex == 0) 
	return 0;
    
    // create another vertices array which aids in adding edges
    
    int *theElementTagVertices = 0;
    int maxEleNum = 0;
    
    TaggedObject *tagdObjPtr;
    TaggedObjectIter &theEles = elements->getComponents();
    while ((tagdObjPtr = theEles()) != 0)
	if (tagdObjPtr->getTag() > maxEleNum)
	    maxEleNum = tagdObjPtr->getTag();

    theElementTagVertices = new int[maxEleNum+1];

    if (theElementTagVertices == 0) {
	cerr << "WARNING Domain::buildEleGraph ";
	cerr << " - Not Enough Memory for ElementTagVertices\n";
	return -1;
    }

    for (int j=0; j<=maxEleNum; j++) theElementTagVertices[j] = -1;

    // now create the vertices with a reference equal to the element number.
    // and a tag which ranges from 0 through numVertex-1
    
    TaggedObjectIter &theEles2 = elements->getComponents();
    
    int count = START_VERTEX_NUM;
    while ((tagdObjPtr = theEles2()) != 0) {
	int ElementTag = tagdObjPtr->getTag();
	Vertex *vertexPtr = new Vertex(count,ElementTag);

	if (vertexPtr == 0) {
	    cerr << "WARNING Domain::buildEleGraph";
	    cerr << " - Not Enough Memory to create ";
	    cerr << count << "th Vertex\n";
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
	cerr << "WARNING Domain::buildEleGraph ";
	cerr << " - Not Enough Memory for NodeTagVertices\n";
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
	    cerr << "WARNING Domain::buildEleGraph";
	    cerr << " - Not Enough Memory to create ";
	    cerr << count << "th Node Vertex\n";
	    delete [] theNodeTagVertices;
	    return -1;
	}
    }

    // now add the the Elements to the nodes
    Element *elePtr;
    TaggedObjectIter &theEles3 = elements->getComponents();
    
    while((tagdObjPtr = theEles3()) != 0) {
	elePtr = (Element *)tagdObjPtr;
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



// a method which will only remove a node from the partitioned domain
// it does not touch the subdomains .. can be dangerous to use.
Node *
PartitionedDomain::removeExternalNode(int tag)
{
    return (this->Domain::removeNode(tag));        
}

Graph &
PartitionedDomain::getSubdomainGraph(void)
{

    // delete the old always - only object that will 
    // use this is a DomainBalancer & it is always looking for latest
    if (mySubdomainGraph != 0) {
	delete mySubdomainGraph;
	mySubdomainGraph = 0;
    }

    // create a new graph
    if (mySubdomainGraph == 0)
        mySubdomainGraph = new Graph(this->getNumSubdomains()+START_VERTEX_NUM);

    if (mySubdomainGraph == 0) // if still 0 try a smaller one
        mySubdomainGraph = new Graph();    
    
    int numVertex = theSubdomains->getNumComponents();

    // see if quick return

    if (numVertex == 0) 
	return *mySubdomainGraph;
    
    // create another vertices array which aids in adding edges
    
    int *theElementTagVertices = 0;
    int maxEleNum = 0;
    
    TaggedObject *tagdObjPtr;
    TaggedObjectIter &theEles = theSubdomains->getComponents();
    while ((tagdObjPtr = theEles()) != 0)
	if (tagdObjPtr->getTag() > maxEleNum)
	    maxEleNum = tagdObjPtr->getTag();

    theElementTagVertices = new int[maxEleNum+1];

    if (theElementTagVertices == 0) {
	cerr << "WARNING PartitionedDomain::buildEleGraph ";
	cerr << " - Not Enough Memory for ElementTagVertices\n";
	exit(-1);
    }

    for (int j=0; j<=maxEleNum; j++) theElementTagVertices[j] = -1;

    // now create the vertices with a reference equal to the subdomain number.
    // and a tag equal to the subdomain number and a weighed according to 
    // the subdomain cost 
    
    TaggedObjectIter &theEles2 = theSubdomains->getComponents();


    while ((tagdObjPtr = theEles2()) != 0) {
	Subdomain *theSub = (Subdomain *)tagdObjPtr; // upward cast ok as
	                                     // only subdomais can be added

	int ElementTag = tagdObjPtr->getTag();
	Vertex *vertexPtr = new Vertex(ElementTag,ElementTag,theSub->getCost()); 

	if (vertexPtr == 0) {
	    cerr << "WARNING Domain::buildEleGraph";
	    cerr << " - Not Enough Memory to create ";
	    cerr << ElementTag << "th Vertex\n";
	    delete [] theElementTagVertices;
	    exit(-1);
	}

	mySubdomainGraph->addVertex(vertexPtr);
	theElementTagVertices[ElementTag] = ElementTag;
	
    }

    // We now need to determine which theSubdomains are asssociated with each node.
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
	cerr << "WARNING Domain::buildEleGraph ";
	cerr << " - Not Enough Memory for NodeTagVertices\n";
	exit(-1);
    }

    for (int l=0; l<=maxNodNum; l++) theNodeTagVertices[l] = 0;

    // now create the vertices with a reference equal to the node number.
    // and a tag which ranges from 0 through numVertex-1 and placed in
    // theNodeTagVertices at a position equal to the node's tag.

    NodeIter &nodeIter2 = this->getNodes();
    int count = START_VERTEX_NUM;
    while ((nodPtr = nodeIter2()) != 0) {
	int nodeTag = nodPtr->getTag();
	Vertex *vertexPtr = new Vertex(count++,nodeTag);
	theNodeTagVertices[nodeTag] = vertexPtr;

	if (vertexPtr == 0) {
	    cerr << "WARNING Domain::buildEleGraph";
	    cerr << " - Not Enough Memory to create ";
	    cerr << count << "th Node Vertex\n";
	    delete [] theNodeTagVertices;
	    exit(-1);
	}
    }

    // now add the the TheSubdomains to the nodes
    Element *elePtr;
    TaggedObjectIter &theEles3 = theSubdomains->getComponents();
    
    while((tagdObjPtr = theEles3()) != 0) {
	elePtr = (Element *)tagdObjPtr;
	int eleTag = elePtr->getTag();
	const ID &id = elePtr->getExternalNodes();

	int size = id.Size();
	for (int i=0; i<size; i++) 
	    theNodeTagVertices[id(i)]->addEdge(eleTag);
    }

    // now add the edges to the vertices of our element graph;
    // this is done by looping over the Node vertices, getting their 
    // Adjacenecy and adding edges between theSubdomains with common nodes


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
			    mySubdomainGraph->addEdge(vertexTag1,vertexTag2);
			    mySubdomainGraph->addEdge(vertexTag2,vertexTag1);			
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
    
    return *mySubdomainGraph;
}

