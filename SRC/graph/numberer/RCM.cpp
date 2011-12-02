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
// $Date: 2003-10-06 20:54:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/numberer/RCM.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for RCM.
// RCM is an object to perform the Reverse Cuthill-McKee numbering
// scheme on the vertices of a graph. This is done by invoking the
// number() method with the Graph.
//
// What: "@(#) RCM.C, revA"

#include <RCM.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

// Constructor
RCM::RCM(bool gps)
:GraphNumberer(GraphNUMBERER_TAG_RCM),
 numVertex(-1), theRefResult(0), GPS(gps)
{

}

// Destructor
RCM::~RCM()
{
    if (theRefResult != 0)
	delete theRefResult;
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
RCM::number(Graph &theGraph, int startVertex)
{
    // first check our size, if not same make new
    if (numVertex != theGraph.getNumVertex()) {

	// delete the old
	if (theRefResult != 0)
	    delete theRefResult;
	
	numVertex = theGraph.getNumVertex();
	theRefResult = new ID(numVertex);

	if (theRefResult == 0) {
	    opserr << "ERROR:  RCM::number - Out of Memory\n";
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
    int startVertexTag;
    startVertexTag = startVertex;
    
    if (startVertexTag != -1) {
	vertexPtr = theGraph.getVertexPtr(startVertexTag);
	if (vertexPtr == 0) {
	    opserr << "WARNING:  RCM::number - No vertex with tag ";
	    opserr << startVertexTag << "Exists - using first come from iter\n";
	    startVertexTag = -1;
	}
    }	

    // if no starting vertex use the first one we get from the VertexIter
    VertexIter &vertexIter2 = theGraph.getVertices();    
    Vertex *start;

    if (startVertexTag == -1) {
	vertexPtr = vertexIter2();
	start = vertexPtr;

	// if GPS true use gibbs-poole-stodlmyer determine the last 
	// level set assuming a starting vertex and then use one of the 
	// nodes in this set to base the numbering on	
	if (GPS == true) {	
	    
	    int currentMark = numVertex-1;  // marks current vertex visiting.
	    int nextMark = currentMark -1;  // marks where to put next Tag
	    int startLastLevelSet = nextMark;
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
		
		if (startLastLevelSet == currentMark)
		    startLastLevelSet = nextMark;

		// check to see if graph is disconneted
	
		if ((currentMark == nextMark) && (currentMark >= 0)) {

		    // loop over iter till we get a vertex not yet Tmped
		    while (((vertexPtr = vertexIter2()) != 0) && 
			   (vertexPtr->getTmp() != -1)) 
			;
	    
		    nextMark--;
		    startLastLevelSet = nextMark;
		    vertexPtr->setTmp(currentMark);
		    (*theRefResult)(currentMark) = vertexPtr->getTag();

		}
	    
	    }

	    // create an id of the last level set
	    if (startLastLevelSet > 0) {
		ID lastLevelSet(startLastLevelSet);
		for (int i=0; i<startLastLevelSet; i++)
		    lastLevelSet(i) = (*theRefResult)(i);
		
		return this->number(theGraph,lastLevelSet);
	    }

	}
	else // we start with just the first vertex we get
	    vertexPtr = start;
    }


    VertexIter &vertexIter3 = theGraph.getVertices();
    Vertex *otherPtr;

    // set to -1 all the vertices     
    while ((otherPtr = vertexIter3()) != 0)
	otherPtr->setTmp(-1);


    VertexIter &vertexIter4 = theGraph.getVertices();

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
	  // opserr << "WARNING:  RCM::number - Disconnected graph -2 \n ";
	    
	  // loop over iter till we get a vertex not yet Tmped
	    
	  while (((vertexPtr = vertexIter4()) != 0) && 
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
	(*theRefResult)(i) = vertexPtr->getTag();
    }

    return *theRefResult;
}



int
RCM::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
RCM::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}


const ID &
RCM::number(Graph &theGraph, const ID &startVertices)
{

    // first check our size, if not same make new
    if (numVertex != theGraph.getNumVertex()) {

	// delete the old
	if (theRefResult != 0)
	    delete theRefResult;
	
	numVertex = theGraph.getNumVertex();
	theRefResult = new ID(numVertex);

	if (theRefResult == 0) {
	    opserr << "ERROR:  RCM::number - Out of Memory\n";
	    theRefResult = new ID(0);
	    numVertex = 0;
	    return *theRefResult;
	}
    }

    // see if we can do quick return
    
    if (numVertex == 0) 
	return *theRefResult;

    // determine one that gives the min avg profile	    
    int minStartVertexTag =0;
    int minAvgProfile = 0;
    int startVerticesSize = startVertices.Size();
    Vertex *vertexPtr;
    
    int startVertexTag =0;
    
    for (int i=0; i<startVerticesSize; i++) {
	// we first set the Tmp of all vertices to -1, indicating
	// they have not yet been added.
	VertexIter &vertexIter = theGraph.getVertices();
    
	while ((vertexPtr = vertexIter()) != 0)
	    vertexPtr->setTmp(-1);

        startVertexTag = startVertices(i);

	if (startVertexTag != -1) {
	    vertexPtr = theGraph.getVertexPtr(startVertexTag);
	    if (vertexPtr == 0) {
		opserr << "WARNING:  RCM::number - No vertex with tag ";
		opserr << startVertexTag << "Exists - using first come from iter\n";
		startVertexTag = -1;
	    }
	}	

	int currentMark = numVertex-1;  // marks current vertex visiting.
	int nextMark = currentMark -1;  // indiactes where to put next Tag in ID.
	(*theRefResult)(currentMark) = vertexPtr->getTag();
	vertexPtr->setTmp(currentMark);
	int avgProfile = 0;
	VertexIter &vertexIter2 = theGraph.getVertices();    

	// we continue till the ID is full

	while (nextMark >= 0) {

	    // get the current vertex and its adjacency
	    vertexPtr = theGraph.getVertexPtr((*theRefResult)(currentMark));
	    const ID &adjacency = vertexPtr->getAdjacency();

	    // go through the vertices adjacency and add vertices which
	    // have not yet been Tmp'ed to the (*theRefResult)

	    int size = adjacency.Size();
	    for (int j=0; j<size; j++) {
		
		int vertexTag = adjacency(j);
		vertexPtr = theGraph.getVertexPtr(vertexTag);
		if ((vertexPtr->getTmp()) == -1) {
		    vertexPtr->setTmp(nextMark);
		    avgProfile += (currentMark - nextMark);
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

	if (i == 0 || minAvgProfile > avgProfile) {
	    minStartVertexTag = startVertexTag;
	    minAvgProfile = avgProfile;
	}
	
    }

	
    // we number based on minStartVertexTag
    minAvgProfile = 0;

    if (minStartVertexTag != startVertexTag) {

	// we could just call the other numbering routine - 
	// we will    just write it out again for now - CHANGE LATER
	startVertexTag = minStartVertexTag;

	// set all unmarked
	VertexIter &vertexIter = theGraph.getVertices();
	while ((vertexPtr = vertexIter()) != 0)
	    vertexPtr->setTmp(-1);
	
	vertexPtr = theGraph.getVertexPtr(startVertexTag);

	int currentMark = numVertex-1;  // marks current vertex visiting.
	int nextMark = currentMark -1;  // indiactes where to put next Tag in ID.
	(*theRefResult)(currentMark) = vertexPtr->getTag();
	vertexPtr->setTmp(currentMark);
	VertexIter &vertexIter2 = theGraph.getVertices();    

	// we continue till the ID is full

	while (nextMark >= 0) {

	    // get the current vertex and its adjacency
	
	    vertexPtr = theGraph.getVertexPtr((*theRefResult)(currentMark));
	    const ID &adjacency = vertexPtr->getAdjacency();

	    // go through the vertices adjacency and add vertices which
	    // have not yet been Tmp'ed to the (*theRefResult)

	    int size = adjacency.Size();
	    for (int j=0; j<size; j++) {
		
		int vertexTag = adjacency(j);
		vertexPtr = theGraph.getVertexPtr(vertexTag);
		if ((vertexPtr->getTmp()) == -1) {
		    vertexPtr->setTmp(nextMark);
		    minAvgProfile += (currentMark-nextMark);
		    (*theRefResult)(nextMark--) = vertexTag;
		}
	    }

	    // go to the next vertex
	    //  we decrement because we are doing reverse Cuthill-McKee
	
	    currentMark--;

	    // loop over iter till we get a vertex not yet Tmped
	    if ((currentMark == nextMark) && (currentMark >= 0)) {
	      // opserr << "WARNING:  RCM::number - Disconnected graph\n";
		
	      while (((vertexPtr = vertexIter2()) != 0) && 
		     (vertexPtr->getTmp() != -1)) 
		;
	      
	      nextMark--;
	      vertexPtr->setTmp(currentMark);
	      (*theRefResult)(currentMark) = vertexPtr->getTag();
	    }
	}	    
    }	


    // now set the vertex references instead of the vertex tags
    // in the result, we change the Tmp to indicate number and return
    
    for (int j=0; j<numVertex; j++) {
	int vertexTag = (*theRefResult)(j);
	vertexPtr = theGraph.getVertexPtr(vertexTag);
	vertexPtr->setTmp(j+1); // 1 through numVertex
	(*theRefResult)(j) = vertexPtr->getTag();
    }	

    return *theRefResult;
}









