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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/numberer/MyRCM.cpp,v $
                                                                        
                                                                        
// File: ~/graph/numberer/MyRCM.C
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for MyRCM.
// MyRCM is an object to perform the Reverse Cuthill-McKee numbering
// scheme on the vertices of a graph. This is done by invoking the
// number() method with the Graph.
//
// What: "@(#) MyRCM.C, revA"

#include <MyRCM.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

// Constructor
MyRCM::MyRCM(int startVertex, bool minDegreeFlag)
:GraphNumberer(GraphNUMBERER_TAG_MyRCM),
 numVertex(-1), theRefResult(0), startVertexTag(startVertex),
 minDegree(minDegreeFlag)
{

}

// Destructor
MyRCM::~MyRCM()
{
    if (theRefResult != 0)
	delete theRefResult;
}

void 
MyRCM::setStartVertex(int startVertex)
{
    startVertexTag = startVertex;
}

void 
MyRCM::setMinDegreeFlag(bool flag)
{
    minDegree = flag;
}


// const ID &number(Graph &theGraph,int startVertexTag = -1,
//                  bool minDegree = false)
//    Method to perform the Reverse Cuthill-mcKenn numbering scheme. The
// user can supply a starting vertex, if none is provided the first vertex
// returned by the iter is used. If minDegree flag is set to true, at each 
// level set the adjacent vertices not yet added from a vertex in the previos
// level set are added in descending degree. The result of the numbering scheme
// is returned in an ID which contains the references for the vertices. If not
// enough memory to allocate a new ID an ID of size 0 is returned.
//
// side effects: this routine changes the color of the vertices.

const ID &
MyRCM::number(Graph &theGraph, int startVertex)
{
    // first check our size, if not same make new
    
    if (numVertex != theGraph.getNumVertex()) {

	// delete the old
	if (theRefResult != 0)
	    delete theRefResult;
	
	numVertex = theGraph.getNumVertex();
	theRefResult = new ID(numVertex);

	if (theRefResult == 0) {
	    opserr << "ERROR:  MyRCM::number - Out of Memory\n";
	    theRefResult = new ID(0);
	    numVertex = 0;
	    return *theRefResult;
	}
    }

    // see if we can do quick return
    
    if (numVertex == 0) 
	return *theRefResult;
	    

    // we first set the Tmp of all vertices to -1, indicating
    // they have not yet been added.
    
    Vertex *vertexPtr;
    VertexIter &vertexIter = theGraph.getVertices();
    
    while ((vertexPtr = vertexIter()) != 0)
	vertexPtr->setTmp(-1);

    // we now set up; setting our markers and getting first vertex
    if (startVertex != -1)
	startVertexTag = startVertex;
    
    if (startVertexTag != -1) {
	vertexPtr = theGraph.getVertexPtr(startVertexTag);
	if (vertexPtr == 0) {
	    opserr << "WARNING:  MyRCM::number - No vertex with tag ";
	    opserr << startVertexTag << "Exists - using first come from iter\n";
	    startVertexTag = -1;
	}
    }	

    // if no starting vertex use the first one we get from the VertexIter
    
    VertexIter &vertexIter2 = theGraph.getVertices();    
    if (startVertexTag == -1) 
	vertexPtr = vertexIter2();

    int currentMark = numVertex-1;  // marks current vertex visiting.
    int nextMark = currentMark -1;  // indiactes where to put next Tag in ID.
    (*theRefResult)(currentMark) = vertexPtr->getTag();
    vertexPtr->setTmp(currentMark);

    // we continue till the ID is full

    while (nextMark >= 0) {

	// get the current vertex and its adjacency
	
	vertexPtr = theGraph.getVertexPtr((*theRefResult)(currentMark));
	const ID &adjacency = vertexPtr->getAdjacency();

	// go through the vertices adjacency and add vertices which
	// have not yet been Tmp'ed to the (*theRefResult)

	int size = adjacency.Size();
	for (int i=0; i<size; i++) {
	    
	    int vertexTag = adjacency(i);
	    vertexPtr = theGraph.getVertexPtr(vertexTag);
	    if ((vertexPtr->getTmp()) == -1) {
		vertexPtr->setTmp(nextMark);
		(*theRefResult)(nextMark--) = vertexTag;
	    }
	}

	// go to the next vertex
	//  we decrement because we are doing reverse Cuthill-McKee
	
	currentMark--;

	// check to see if graph is disconneted
	
	if ((currentMark == nextMark) && (currentMark >= 0)) {
	    opserr << "WARNING:  MyRCM::number - Disconnected graph\n";
	    
	    // loop over iter till we get a vertex not yet Tmped
	    
	    while (((vertexPtr = vertexIter2()) != 0) && 
		   (vertexPtr->getTmp() != -1)) 
		;
	    
	    nextMark--;
	    vertexPtr->setTmp(currentMark);
	    (*theRefResult)(currentMark) = vertexPtr->getTag();
	}

    }
    
    // now set the vertex references instead of the vertex tags
    // in the result, we change the Tmp to indicate number and return
    
    for (int i=0; i<numVertex; i++) {
	int vertexTag = (*theRefResult)(i);
	vertexPtr = theGraph.getVertexPtr(vertexTag);
	vertexPtr->setTmp(i+1); // 1 through numVertex
	(*theRefResult)(i) = vertexPtr->getRef();
    }

    return *theRefResult;
}



int
MyRCM::sendSelf(int tag, Channel &theChannel)
{
    return 0;
}

int
MyRCM::recvSelf(int tag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}







const ID &
MyRCM::number(Graph &theGraph, const ID &startVertices)
{
    // first check our size, if not same make new
    
    if (numVertex != theGraph.getNumVertex()) {

	// delete the old
	if (theRefResult != 0)
	    delete theRefResult;
	
	numVertex = theGraph.getNumVertex();
	theRefResult = new ID(numVertex);

	if (theRefResult == 0) {
	    opserr << "ERROR:  MyRCM::number - Out of Memory\n";
	    theRefResult = new ID(0);
	    numVertex = 0;
	    return *theRefResult;
	}
    }

    // see if we can do quick return
    
    if (numVertex == 0) 
	return *theRefResult;

    ID copyStart(startVertices);

    // we determine which node to start with
    int minStartVertexTag =0;
    int minAvgProfile = 0;
    int startVerticesSize = startVertices.Size();
    
    for (int j=0; j<startVerticesSize; j++) {
	// we first set the Tmp of all vertices to -1, indicating
	// they have not yet been added.
    
	Vertex *vertexPtr;
	VertexIter &vertexIter = theGraph.getVertices();
    
	while ((vertexPtr = vertexIter()) != 0)
	    vertexPtr->setTmp(-1);

	// we now set up; setting our markers and set first vertices
	VertexIter &vertexIter2 = theGraph.getVertices();    
	int currentMark = numVertex-1;  // marks current vertex visiting.
	int nextMark = currentMark-1;

	for (int k=0; k<startVerticesSize; k++)
	    if (k != j)
		copyStart(k) = 0;
	    else	
		copyStart(k) = 1;
	
	vertexPtr = theGraph.getVertexPtr(startVertices(j));	
	(*theRefResult)(currentMark) = vertexPtr->getTag();
	vertexPtr->setTmp(currentMark);

	int numFromStart = 1;
	int avgProfile = 1;	
	while (numFromStart < startVerticesSize) {
	    // get the current vertex and its adjacency

	    vertexPtr = theGraph.getVertexPtr((*theRefResult)(currentMark));
	    const ID &adjacency = vertexPtr->getAdjacency();

	    // go through the vertices adjacency and add vertices which
	    // have not yet been Tmp'ed to the (*theRefResult)

	    int size = adjacency.Size();
	    for (int i=0; i<size; i++) {
		int vertexTag = adjacency(i);
		int loc =startVertices.getLocation(vertexTag);
		if (loc >= 0) {
		    vertexPtr = theGraph.getVertexPtr(vertexTag);
		    if ((vertexPtr->getTmp()) == -1) {
			vertexPtr->setTmp(nextMark);
			copyStart(loc) = 1;
			numFromStart++;
			avgProfile += currentMark - nextMark;			
			(*theRefResult)(nextMark--) = vertexTag;
		    }
		}
	    }

	    // go to the next vertex
	    //  we decrement because we are doing reverse Cuthill-McKee
	
	    currentMark--;

	    // check to see if graph is disconneted
	
	    if (currentMark == nextMark && numFromStart < startVerticesSize) {
		// loop over iter till we get a vertex not yet included
		
		for (int l=0; l<startVerticesSize; l++)
		    if (copyStart(l) == 0) {
			int vertexTag = startVertices(l);			
			vertexPtr = theGraph.getVertexPtr(vertexTag);			
			nextMark--;
			copyStart(l) = 1;
			vertexPtr->setTmp(currentMark);
			numFromStart++;
			(*theRefResult)(currentMark) = vertexPtr->getTag();
			l =startVerticesSize;
		    }
	    }
	}
	
	currentMark = numVertex-1; // set current to the first again
	nextMark =  numVertex - startVerticesSize -1;

	// we continue till the ID is full

	while (nextMark >= 0) {
	    // get the current vertex and its adjacency
	
	    vertexPtr = theGraph.getVertexPtr((*theRefResult)(currentMark));
	    const ID &adjacency = vertexPtr->getAdjacency();

	    // go through the vertices adjacency and add vertices which
	    // have not yet been Tmp'ed to the (*theRefResult)

	    int size = adjacency.Size();
	    for (int i=0; i<size; i++) {
	    
		int vertexTag = adjacency(i);
		vertexPtr = theGraph.getVertexPtr(vertexTag);
		if ((vertexPtr->getTmp()) == -1) {
		    vertexPtr->setTmp(nextMark);
		    avgProfile += currentMark - nextMark;

		    (*theRefResult)(nextMark--) = vertexTag;
		}
	    }
	    // go to the next vertex
	    //  we decrement because we are doing reverse Cuthill-McKee
	
	    currentMark--;

	    // check to see if graph is disconneted
	
	    if ((currentMark == nextMark) && (currentMark >= 0)) {
	    
		// loop over iter till we get a vertex not yet Tmped
		
		while (((vertexPtr = vertexIter2()) != 0) && 
		       (vertexPtr->getTmp() != -1)) 
		    ;
		
		nextMark--;
		vertexPtr->setTmp(currentMark);
		(*theRefResult)(currentMark) = vertexPtr->getTag();
	    }
	}
	
	if (j == 0 || minAvgProfile > avgProfile) {
	    minStartVertexTag = startVertices(j);
	    minAvgProfile = avgProfile;
	}

    }

    
    // now we numebr based on minStartVErtexTag

    // we first set the Tmp of all vertices to -1, indicating
    // they have not yet been added.
    
    Vertex *vertexPtr;
    VertexIter &vertexIter = theGraph.getVertices();
    
    while ((vertexPtr = vertexIter()) != 0)
	vertexPtr->setTmp(-1);

    // we now set up; setting our markers and set first vertices
    VertexIter &vertexIter2 = theGraph.getVertices();    
    int currentMark = numVertex-1;  // marks current vertex visiting.
    int nextMark = currentMark-1;
    
    vertexPtr = theGraph.getVertexPtr(minStartVertexTag);	
    (*theRefResult)(currentMark) = vertexPtr->getTag();
    vertexPtr->setTmp(currentMark);
    currentMark--;	
    
    int loc = startVertices.getLocation(minStartVertexTag);
    for (int k=0; k<startVerticesSize; k++)
	if (k != loc)
	    copyStart(k) = 0;
     

    int numFromStart = 1;
    while (numFromStart < startVerticesSize) {
	// get the current vertex and its adjacency

	vertexPtr = theGraph.getVertexPtr((*theRefResult)(currentMark));
	const ID &adjacency = vertexPtr->getAdjacency();

	// go through the vertices adjacency and add vertices which
	// have not yet been Tmp'ed to the (*theRefResult)

	int size = adjacency.Size();
	for (int i=0; i<size; i++) {
	    int vertexTag = adjacency(i);
	    int loc =startVertices.getLocation(vertexTag);
	    if (loc >= 0) {
		vertexPtr = theGraph.getVertexPtr(vertexTag);
		if ((vertexPtr->getTmp()) == -1) {
		    vertexPtr->setTmp(nextMark);
		    copyStart(loc) = 1;
		    numFromStart++;
		    (*theRefResult)(nextMark--) = vertexTag;
		}
	    }
	}

	// go to the next vertex
	//  we decrement because we are doing reverse Cuthill-McKee
	
	currentMark--;

	// check to see if graph is disconneted
	
	if (currentMark == nextMark && numFromStart < startVerticesSize) {
	    // loop over iter till we get a vertex not yet included
		
	    for (int l=0; l<startVerticesSize; l++)
		if (copyStart(l) == 0) {
		    int vertexTag = startVertices(l);			
		    vertexPtr = theGraph.getVertexPtr(vertexTag);			
		    nextMark--;
		    copyStart(l) = 1;
		    vertexPtr->setTmp(currentMark);
		    numFromStart++;
		    (*theRefResult)(currentMark) = vertexPtr->getTag();
		    l =startVerticesSize;
		}
	}	
    }
	
    currentMark = numVertex-1; // set current to the first again
    nextMark =  numVertex - startVerticesSize -1;


    currentMark = numVertex-1; // set current to the first again
    
    // we continue till the ID is full
    while (nextMark >= 0) {

	// get the current vertex and its adjacency
	
	vertexPtr = theGraph.getVertexPtr((*theRefResult)(currentMark));
	const ID &adjacency = vertexPtr->getAdjacency();
	
	// go through the vertices adjacency and add vertices which
	// have not yet been Tmp'ed to the (*theRefResult)
	
	int size = adjacency.Size();
	for (int i=0; i<size; i++) {
	    
	    int vertexTag = adjacency(i);
	    vertexPtr = theGraph.getVertexPtr(vertexTag);
	    if ((vertexPtr->getTmp()) == -1) {
		vertexPtr->setTmp(nextMark);
		(*theRefResult)(nextMark--) = vertexTag;
	    }
	}

	// go to the next vertex
	//  we decrement because we are doing reverse Cuthill-McKee
	
	currentMark--;
	
	// check to see if graph is disconneted
	
	if ((currentMark == nextMark) && (currentMark >= 0)) {
	    opserr << "WARNING:  MyRCM::number - Disconnected graph ";
	    
	    // loop over iter till we get a vertex not yet Tmped
	    
	    while (((vertexPtr = vertexIter2()) != 0) && 
		   (vertexPtr->getTmp() != -1)) 
		;
		
	    nextMark--;
	    vertexPtr->setTmp(currentMark);
	    (*theRefResult)(currentMark) = vertexPtr->getTag();
	}
	
    }

    // now set the vertex references instead of the vertex tags
    // in the result, we change the Tmp to indicate number and return
    
    for (int m=0; m<numVertex; m++) {
	int vertexTag = (*theRefResult)(m);
	vertexPtr = theGraph.getVertexPtr(vertexTag);
	vertexPtr->setTmp(m+1); // 1 through numVertex
	(*theRefResult)(m) = vertexPtr->getRef();
    }

    return *theRefResult;
}



