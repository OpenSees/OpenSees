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
// $Date: 2003-02-14 23:01:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/Graph.cpp,v $
                                                                        
                                                                        
// File: ~/graph/graph/Graph.C
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class implementation for Graph.
//
// What: "@(#) Graph.C, revA"

#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <ArrayOfTaggedObjects.h>

Graph::Graph()
:myVertices(0), theVertexIter(0), numEdge(0)
{
    myVertices = new ArrayOfTaggedObjects(32);
    theVertexIter = new VertexIter(myVertices);
}


Graph::Graph(int numVertices)
:myVertices(0), theVertexIter(0), numEdge(0)
{
    myVertices = new ArrayOfTaggedObjects(numVertices);
    theVertexIter = new VertexIter(myVertices);
}


Graph::Graph(TaggedObjectStorage &theVerticesStorage)
:myVertices(&theVerticesStorage), theVertexIter(0), numEdge(0)
{
    theVerticesStorage.clearAll();
    theVertexIter = new VertexIter(myVertices);
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
	opserr << "WARNING ArrayGraph::addVertex";
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
		    opserr << "WARNING ArrayGraph::addVertex";
		    opserr << " - vertex with adjacent vertex not in graph\n";
		    return false;
		}		
	    }
	}
    }

    bool result = myVertices->addComponent(vertexPtr);
    if (result == false) {
	opserr << "WARNING ArrayGraph::addVertex";
	opserr << " - vertex could not be stored in TaggedObjectStorage object\n";
    }
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
    int result;
    if ((result = vertex1->addEdge(otherVertexTag)) == 0) {
	if ((result = vertex2->addEdge(vertexTag)) == 0) 
	    numEdge++;
	else {
	    opserr << "WARNING Graph::addEdge() - " << vertexTag;
	    opserr << " has not been added to " << otherVertexTag;
	    opserr << " adjacency - yet vica-versa ok.\n";
	    return -2;
	}
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

