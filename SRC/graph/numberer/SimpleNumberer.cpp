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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/numberer/SimpleNumberer.cpp,v $
                                                                        
                                                                        
// File: ~/graph/numberer/SimpleNumberer.C
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for SimpleNumberer.
// SimpleNumberer is an object to perform a simple numbering of the vertices.
// It does this by using the graphs VertexIter and assigning the numbers as
// it comes across the vertices.
//
// What: "@(#) SimpleNumberer.C, revA"

#include <SimpleNumberer.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

// Constructor
SimpleNumberer::SimpleNumberer()
:GraphNumberer(GraphNUMBERER_TAG_SimpleNumberer),
 numVertex(-1), theRefResult(0)
{

}

// Destructor
SimpleNumberer::~SimpleNumberer()
{
    if (theRefResult != 0)
	delete theRefResult;
}

// const ID &number(Graph &theGraph)
//                  


const ID &
SimpleNumberer::number(Graph &theGraph, int lastVertex)
{
    // first check our size, if not same make new
    
    if (numVertex != theGraph.getNumVertex()) {

	if (theRefResult != 0)
	    delete theRefResult;
	
	numVertex = theGraph.getNumVertex();
	theRefResult = new ID(numVertex);

	if (theRefResult == 0) {
	    opserr << "ERROR:  SimpleNumberer::number - Out of Memory\n";
	    theRefResult = new ID(0);
	    numVertex = 0;
	    return *theRefResult;
	}
    }

    // see if we can do quick return
    
    if (numVertex == 0) 
	return *theRefResult;
	    

    // Now we go through the iter and assign the numbers

    if (lastVertex != -1) {
	opserr << "WARNING:  SimpleNumberer::number -";
	opserr << " - does not deal with lastVertex";
    }
    
    Vertex *vertexPtr;
    VertexIter &vertexIter = theGraph.getVertices();
    int count = 0;
    
    while ((vertexPtr = vertexIter()) != 0) {

	(*theRefResult)(count++) = vertexPtr->getRef();
	vertexPtr->setColor(count);
    }
    
    return *theRefResult;
}


int
SimpleNumberer::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
SimpleNumberer::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

const ID &
SimpleNumberer::number(Graph &theGraph, const ID &startVertices)
{
    return this->number(theGraph);
}
