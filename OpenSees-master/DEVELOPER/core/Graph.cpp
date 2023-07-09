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
// $Date: 2007-01-09 19:20:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/Graph.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class implementation for Graph.
//
// What: "@(#) Graph.C, revA"

#include <stdlib.h>

#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <MapOfTaggedObjects.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>

Graph::Graph()
  :myVertices(0), theVertexIter(0), numEdge(0), nextFreeTag(START_VERTEX_NUM)
{
    myVertices = new MapOfTaggedObjects();
    theVertexIter = new VertexIter(myVertices);
}


Graph::Graph(int numVertices)
  :myVertices(0), theVertexIter(0), numEdge(0), nextFreeTag(START_VERTEX_NUM)
{
    myVertices = new MapOfTaggedObjects();
    theVertexIter = new VertexIter(myVertices);
}


Graph::Graph(TaggedObjectStorage &theVerticesStorage)
  :myVertices(&theVerticesStorage), theVertexIter(0), numEdge(0), nextFreeTag(START_VERTEX_NUM)
{
  TaggedObject *theObject;
  TaggedObjectIter &theObjects = theVerticesStorage.getComponents();
  while ((theObject = theObjects()) != 0) 
    if (theObject->getTag() > nextFreeTag)
      nextFreeTag = theObject->getTag() + 1;

  theVerticesStorage.clearAll();
  theVertexIter = new VertexIter(myVertices);
}
    

Graph::Graph(Graph &other) 
  :myVertices(0), theVertexIter(0), numEdge(0), nextFreeTag(START_VERTEX_NUM)
{
  myVertices = new MapOfTaggedObjects();
  theVertexIter = new VertexIter(myVertices);

  VertexIter &otherVertices = other.getVertices();
  Vertex *vertexPtr;

  // loop through other creating vertices if tag not the same in this
  while ((vertexPtr = otherVertices()) != 0) {
    int vertexTag = vertexPtr->getTag();
    int vertexRef = vertexPtr->getRef();
    vertexPtr = new Vertex(vertexTag, vertexRef);
    if (vertexPtr == 0) {
      opserr << "Graph::Graph - out of memory\n";
      return;
    }
    this->addVertex(vertexPtr, false);
  }

  // loop through other adding all the edges that exist in other
  VertexIter &otherVertices2 = other.getVertices();
  while ((vertexPtr = otherVertices2()) != 0) {
    int vertexTag = vertexPtr->getTag();
    const ID &adjacency = vertexPtr->getAdjacency();
    for (int i=0; i<adjacency.Size(); i++) {
      if (this->addEdge(vertexTag, adjacency(i)) < 0) {
	opserr << "Graph::merge - could not add an edge!\n";
	return;
      }
    }
  }
}

Graph::~Graph()
{
    // invoke delete on the Vertices
    myVertices->clearAll();
    
    if (myVertices != 0)
	delete myVertices;
    
    if (theVertexIter != 0)
	delete theVertexIter;
}


// bool addVertex(int vertexTag, int vertexRef, 
//		 int vertexWeight=0, int vertexColor = 0) 
// Method to add a vertex to the graph. If the adjacency list
// of the vertex is not empty the graph will first check to see all
// vertices in the the the vertices adjacency list exist in the graph
// before the vertex is added. 

bool
Graph::addVertex(Vertex *vertexPtr, bool checkAdjacency)
{
    // check the vertex * and its adjacency list
    if (vertexPtr == 0) {
	opserr << "WARNING Graph::addVertex";
	opserr << " - attempting to add a NULL vertex*\n";
	return false;
    }

    if (checkAdjacency == true) {
	if (vertexPtr->getDegree() != 0) {
	    const ID &adjacency = vertexPtr->getAdjacency();
	    int size = adjacency.Size();
	    for (int i=0; i<size; i++) {
		Vertex *other = this->getVertexPtr(adjacency(i));
		if (other == 0) {
		    opserr << "WARNING Graph::addVertex";
		    opserr << " - vertex with adjacent vertex not in graph\n";
		    return false;
		}		
	    }
	}
    }


    bool result = myVertices->addComponent(vertexPtr);
    if (result == false) {
      opserr << *this;
      opserr << "BAD VERTEX\n: " << *vertexPtr;
	opserr << "WARNING Graph::addVertex";
	opserr << " - vertex could not be stored in TaggedObjectStorage object\n";
    }


    // check nextFreeTag
    if (vertexPtr->getTag() >= nextFreeTag)
      nextFreeTag = vertexPtr->getTag() + 1;

    return result;
}


// int addEdge(int vertexTag, int otherVertexTag);
// Causes the Graph to add an edge {\em (vertexTag,otherVertexTag)} to
// the Graph. A check is first made to see if vertices with tags given by
// {\em vertexTag} and {\em otherVertexTag} exist in the graph. If they
// do not exist a $-1$ is returned, otherwise the method invokes {\em
// addEdge()} on each of the corresponding vertices in the 
// graph. Returns $0$ if sucessfull, a negative number if not.

int 
Graph::addEdge(int vertexTag, int otherVertexTag)
{
    // get pointers to the vertices, if one does not exist return

    Vertex *vertex1 = this->getVertexPtr(vertexTag);
    Vertex *vertex2 = this->getVertexPtr(otherVertexTag);
    if ((vertex1 == 0) || (vertex2 == 0)) {
	opserr << "WARNING Graph::addEdge() - one or both of the vertices ";
	opserr << vertexTag << " " << otherVertexTag << " not in Graph\n";
	return -1;
    }

    // add an edge to each vertex
    int result = vertex1->addEdge(otherVertexTag);
	if (result == 1)
		return 0;  // already there
	else if (result == 0) {  // added to vertexTag now add to other
		if ((result = vertex2->addEdge(vertexTag)) == 0) {
			numEdge++;
		}
		else {
			opserr << " WARNING Graph::addEdge() - " << vertexTag;
			opserr << " added to " << otherVertexTag;
			opserr << " adjacency - but already there in otherVertexTag!.\n";
			opserr << *this; exit(0);
			return -2;
		}
	} else {
			opserr << " WARNING Graph::addEdge() - " << vertexTag;
			opserr << " added to " << otherVertexTag;
			opserr << " adjacency - but not vica versa!.\n";
			opserr << *this; exit(0);
			return -2;
	}
    return result;
}


Vertex *
Graph::getVertexPtr(int vertexTag)
{
    TaggedObject *res = myVertices->getComponentPtr(vertexTag);
    if (res == 0) return 0;
    Vertex *result = (Vertex *)res;
    return result;
}


VertexIter &
Graph::getVertices(void) 
{
    // reset the iter and then return it
    theVertexIter->reset();
    return *theVertexIter;
}


int 
Graph::getNumVertex(void) const
{
    return myVertices->getNumComponents();
}

int 
Graph::getNumEdge(void) const
{
    return numEdge;
}

int 
Graph::getFreeTag(void) 
{
  return nextFreeTag;
}

Vertex *
Graph::removeVertex(int tag, bool flag)
{
    TaggedObject *mc = myVertices->removeComponent(tag);
    if (mc == 0) return 0;
    Vertex *result = (Vertex *)mc;
    
    if (flag == true) { // remove all edges associated with the vertex
	opserr << "Graph::removeVertex(int tag, bool flag = true)";
	opserr << " - no code to remove edges yet\n";
	return 0;
    }
    return result;
}


int
Graph::merge(Graph &other) {

  int result =0;
  VertexIter &otherVertices = other.getVertices();
  Vertex *vertexPtrOther;

  // loop through other creating vertices if tag not the same in this
  while ((vertexPtrOther = otherVertices()) != 0) {
    int vertexTag = vertexPtrOther->getTag();
    Vertex *vertexPtr = this->getVertexPtr(vertexTag);
    if (vertexPtr == 0) {
      int vertexRef = vertexPtrOther->getRef();
      vertexPtr = new Vertex(vertexTag, vertexRef);
      if (vertexPtr == 0) {
	opserr << "Graph::merge - out of memory\n";
	return -1;
      }
      this->addVertex(vertexPtr, false);
    }
  }


  // loop through other adding all the edges that exist in other
  VertexIter &otherVertices2 = other.getVertices();
  while ((vertexPtrOther = otherVertices2()) != 0) {
    int vertexTag = vertexPtrOther->getTag();
    const ID &adjacency = vertexPtrOther->getAdjacency();
    for (int i=0; i<adjacency.Size(); i++) {
      if (this->addEdge(vertexTag, adjacency(i)) < 0) {
	opserr << "Graph::merge - could not add an edge!\n";
	return -2;	
      }
    }
  }
  
  return result;
}


void 
Graph::Print(OPS_Stream &s, int flag)
{
    myVertices->Print(s, flag);
}


OPS_Stream &operator<<(OPS_Stream &s, Graph &M)
{
  M.Print(s);
  return s;
}


int 
Graph::sendSelf(int commitTag, Channel &theChannel)
{
  // check not a datastore .. 
  if (theChannel.isDatastore() != 0) {
    opserr << "Graph::sendSelf() - does not at present send to a database\n";
    return -1;
  }

  int numVertex = this->getNumVertex();

  // send numEdge & the number of vertices
  static ID idData(2);
  idData(0) = numEdge;
  idData(1) = numVertex;

  if (theChannel.sendID(0, commitTag, idData) < 0) {
    opserr << "Graph::sendSelf() - failed to send the id\n";
    return -3;
  }

  if (numVertex != 0) {
    int *vertexData = new int[5 * numVertex + 2 * numEdge];
    Vector vertexWeights(numVertex);
    if (vertexData != 0) {
      VertexIter &theVertices = this->getVertices();
      Vertex *vertexPtr;
      int adjacencyLocation = 5 * numVertex;
      int vertexLocation = 0;
      int weightLoc = 0;
      while ((vertexPtr = theVertices()) != 0) {
	int tag = vertexPtr->getTag();
	int color = vertexPtr->getColor();
	int ref = vertexPtr->getRef();
	int tmp = vertexPtr->getTmp();
	const ID &adjacency = vertexPtr->getAdjacency();
	int adjSize = adjacency.Size();
	vertexData[vertexLocation++] = tag;
	vertexData[vertexLocation++] = ref;
	vertexData[vertexLocation++] = color;
	vertexData[vertexLocation++] = tmp;
	vertexData[vertexLocation++] = adjSize;
	for (int i=0; i<adjSize; i++)
	  vertexData[adjacencyLocation++] = adjacency(i);	  
	vertexWeights[weightLoc++] = vertexPtr->getWeight();

      }  

      ID verticesData(vertexData, 5*numVertex + 2*numEdge, true);
      if (theChannel.sendID(0, commitTag, verticesData) < 0) {
	opserr << "Graph::sendSelf() - failed to send the id\n";
	return -3;
      }

      if (theChannel.sendVector(0, commitTag, vertexWeights) < 0) {
	opserr << "Graph::sendSelf() - failed to send the id\n";
	return -3;
      }
    }
  }

  /* ************ OLD --- SENDING IND VERTICES *********************
  VertexIter &theVertices = this->getVertices();
  Vertex *vertexPtr;
  while ((vertexPtr = theVertices()) != 0) {
    if (vertexPtr->sendSelf(commitTag, theChannel) < 0) {
      opserr << "Graph::sendSelf() - failed to send a vertex: " << *vertexPtr;
      return -3;
    }
  }
  ********************************************************************/

  return 0;
}


int 
Graph::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) 
{
  // check not from a datastore
  if (theChannel.isDatastore() != 0) {
    opserr << "Graph::recvSelf() - at present does not receive from a database\n";
    return -1;
  }

  // check blank
  if (this->getNumVertex() != 0) {
    opserr << "Graph::recvSelf() - can only receive to an empty graph at present\n";

    numEdge = 0;
    myVertices->clearAll();
  }

  // recv numEdge & numVertices
  static ID idData(2);
  if (theChannel.recvID(0, commitTag, idData) < 0) {
    opserr << "Graph::recvSelf() - failed to receive the id\n";
    return -3;
  }

  numEdge = idData(0);
  int numVertex = idData(1);

  if (numVertex != 0) {
    
    int *vertexData = new int[5 * numVertex + 2 * numEdge];
    if (vertexData != 0) {
      ID verticesData(vertexData, 5*numVertex + 2 * numEdge, true);
      if (theChannel.recvID(0, commitTag, verticesData) < 0) {
	opserr << "Graph::recvSelf() - failed to receive the id\n";
	return -3;
      }
      Vector vertexWeights(numVertex);
      if (theChannel.recvVector(0, commitTag, vertexWeights) < 0) {
	opserr << "Graph::recvSelf() - failed to receive the weights\n";
	return -3;
      }
      
      int adjacencyLocation = 5 * numVertex;
      int vertexLocation = 0;
      for (int i=0; i<numVertex; i++) {
	int tag = vertexData[vertexLocation++];
	int ref = vertexData[vertexLocation++];
	int color = vertexData[vertexLocation++];
	int tmp = vertexData[vertexLocation++];
	int adjSize = vertexData[vertexLocation++];
	Vertex *theVertex = new Vertex(tag, ref);
	if (theVertex == 0) {
	  opserr << "Graph::recvSelf() - out of memory\n";
	  return -4;
	}
	theVertex->setColor(color);
	theVertex->setTmp(tmp);
	theVertex->setWeight(vertexWeights(i));
	for (int i=0; i<adjSize; i++) {
	  int edge = vertexData[adjacencyLocation++];
	  theVertex->addEdge(edge);
	}
	this->addVertex(theVertex, false);
      }
    } else {
      opserr << "Graph::recvSelf() - out of memory\n";
      return -5;
    }
  }
  /* ************ OLD --- RECEIVING IND VERTICES *********************

  // for each vertex to be received, create it, receive it and then add it to the graph
  for (int i=0; i<numVertex; i++) {
    Vertex *theVertex = new Vertex(0, 0);
    if (theVertex == 0) {
      opserr << "Graph::recvSelf() - out of memory\n";
      return -4;
    }
    if (theVertex->recvSelf(commitTag, theChannel, theBroker) < 0) {
      opserr << "Graph::recvSelf() - vertex failed to receive itself\n";      
      return -5;
    }
    this->addVertex(theVertex, false);
  }
  ****************************************************************** */

  return 0;
}
