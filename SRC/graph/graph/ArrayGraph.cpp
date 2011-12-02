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
// $Date: 2003-02-14 23:01:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/ArrayGraph.cpp,v $
                                                                        
                                                                        
// File: ~/graph/ArrayGraph.C
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for ArrayGraph.
// The vertices in an ArrayGraph are held in an array. This is more efficient
// than holding them in a List data structure, but problems can arise with
// large Graphs in getting enough contiguous memory for the array.
//
// What: "@(#) ArrayGraph.C, revA"

#include <ArrayGraph.h>
#include <Vertex.h>
#include <AnalysisModel.h>
#include <DOF_Group.h>
#include <FE_Element.h>

ArrayGraph::ArrayGraph(int arraySize)
:numVertex(0), numEdge(0), sizeVertices(arraySize), lastEmpty(0),
 theVertices(0), myIter(*this) 
{
    // we now try and get an array of size arraySize
    theVertices = new Vertex *[arraySize];
    if (theVertices == 0)  {
	opserr << "Warning ArrayGraph::ArrayGraph";
	opserr << " - no contiguous memory block big enough available\n";
	sizeVertices = 0;
    }
    
    // zero the pointers
    
    for (int i=0; i<arraySize; i++)
	theVertices[i] = 0;
}

ArrayGraph::~ArrayGraph()
{
    // delete all the vertices, then delete theVertices array
    
    if (theVertices != 0) {
	for (int i=0; i<numVertex; i++)
	    if (theVertices[i] != 0)
		delete theVertices[i];
	delete [] theVertices;
    }
}


// bool addVertex(int vertexTag, int vertexRef, 
//		 int vertexWeight=0, int vertexColor = 0) 
// Method to add a vertex to the graph. If the adjacency list
// of the vertex is not empty the graph will first check to see all
// vertices in the the the vertices adjacency list exist in the graph
// before the vertex is added. It then checks if it neeeds a new array
// and if so creates one, i.e. if the {\em arraySize} $=$ {\em
// numVertex} it creates a new array, whose size is double the original
// and copies the pointers to the vertices, before invoking {\em
// delete()} on the old array. It now tries to add the vertex in the
// array at location {\em vertexTag}. If this fails it adds at the first
// empty location it comes to. Returns a 0 if successfull addition, a
// $-1$ otherwise and a message to opserr explaining the problem. \\ 


bool
ArrayGraph::addVertex(Vertex *vertexPtr)
{
    // check the vertex * and its adjacency list
    if (vertexPtr == 0) {
	opserr << "WARNING ArrayGraph::addVertex";
	opserr << " - attempting to add a NULL vertex*\n";
	return false;
    }
    
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

    // check if we have room to place the vertex
    if (numVertex == sizeVertices) {

	int newSize = sizeVertices*2;
	Vertex **newVertices = new Vertex *[newSize];
	
	if (newVertices == 0) {
	    opserr << "WARNING ArrayGraph::addVertex";
	    opserr << " - out of contiguous memory could not create a new array";
	    delete vertexPtr;
	    return false;
	}

	// copy the old and 0 the extra, then delete the old
	for (int i=0; i<sizeVertices; i++)
	    newVertices[i] = theVertices[i];
	for (int j=sizeVertices; j<newSize; j++)
	    newVertices[j] = 0;

	delete [] theVertices;
	
	theVertices = newVertices;
        sizeVertices = newSize;
    }
    
    // now see if we can add the Vertex into the array in its vertexTag location
    int vertexTag = vertexPtr->getTag();
    if ((vertexTag >= 0) && (vertexTag < sizeVertices) && 
	(theVertices[vertexTag] == 0)) {

	theVertices[vertexTag]= vertexPtr;
	numVertex++;
	return 0;

    } else {
	
	// we have to serach through the array till we find an empty spot
	
	for (int i=0; i<sizeVertices; i++) 
	    if (theVertices[i] == 0) {
		
		lastEmpty = i+1; // stores the lastEmpty place
		theVertices[i] = vertexPtr;
		numVertex++;
		return true;
	    }
    }

    // we should never get here

    return false;
}

// Vertex *getVertexPtr(int vertexTag);} 
// A method which returns a pointer to the vertex whose tag is given by 
// vertexTag. The method first looks at location {\em vertexTag} for the
// vertex, otherwise it must search through the array until it finds the
// vertex it is looking for. If no such vertex exists in the graph $0$ is
// returned.\\ 

Vertex *
ArrayGraph::getVertexPtr(int vertexTag) 
{
    // check first to see if it's in a nice position
    if ((vertexTag >= 0) && (vertexTag < sizeVertices) &&
	(theVertices[vertexTag] != 0) &&
	(theVertices[vertexTag]->getTag() == vertexTag)) {

	return theVertices[vertexTag];	
    }
    // it's not nicely positioned, we have to search
    // through theVertices until we find it
    
    else 
	for (int i=0; i<sizeVertices; i++)
	    if ((theVertices[i] != 0) && 
		(theVertices[i]->getTag() == vertexTag)){

		return theVertices[i];
	    }

    // else the vertex is not there
    
    return 0;
}

// int addEdge(int vertexTag, int otherVertexTag);
// Causes the Graph to add an edge {\em (vertexTag,otherVertexTag)} to
// the Graph. A check is first made to see if vertices with tags given by
// {\em vertexTag} and {\em otherVertexTag} exist in the graph. If they
// do not exist a $-1$ is returned, otherwise the method invokes {\em
// addEdge()} on each of the corresponding vertices in the 
// graph. Returns $0$ if sucessfull, a negative number if not.

int 
ArrayGraph::addEdge(int vertexTag, int otherVertexTag)
{
    // get pointers to the vertices, if one does not exist return

    Vertex *vertex1 = this->getVertexPtr(vertexTag);
    Vertex *vertex2 = this->getVertexPtr(otherVertexTag);
    if ((vertex1 == 0) || (vertex2 == 0))
	return -1;

    // add an edge to each vertex
    int result;
    if ((result = vertex1->addEdge(otherVertexTag)) == 0)
	if ((result = vertex2->addEdge(vertexTag)) == 0)
	    numEdge++;

    return result;
}

// VertexIter \&getVertices(void);} 
// A method which first invokes {\em reset()} on the graphs ArrayVertexIter
// and then returns a reference to this iter.\\

VertexIter &
ArrayGraph::getVertices(void) 
{
    // reset the iter and then return it
    myIter.reset();
    return myIter;
}


int 
ArrayGraph::getNumVertex(void) const
{
    return numVertex;
}


int 
ArrayGraph::getNumEdge(void) const
{
    return numEdge;
}


int 
ArrayGraph::getArraySize(void) const
{
    return sizeVertices;
}


void 
ArrayGraph::Print(OPS_Stream &s) const
{
    s << numVertex << " " << numEdge << endln;
    
    Vertex *vertexPtr;
    
    // loop over the vertices and print each
    
    for (int i=0; i<sizeVertices; i++) {
	vertexPtr = theVertices[i];
	if (vertexPtr != 0)
	    vertexPtr->Print(s);
    }
}

OPS_Stream &operator<<(OPS_Stream &s, const ArrayGraph &M)
{
  M.Print(s);
  return s;
}

