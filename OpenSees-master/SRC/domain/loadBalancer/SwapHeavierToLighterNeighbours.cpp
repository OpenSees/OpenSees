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
// $Date: 2003-02-14 23:00:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/loadBalancer/SwapHeavierToLighterNeighbours.cpp,v $
                                                                        
                                                                        
 // File: ~/domain/loadBalancer/HeavierToSmallerNeighbours.C
// 
// Written: fmk 
// Created: Fri Aug 29 17:43:25 1997
// Revision: A
//
// Description: This file contains the class definition for 
// HeavierToSmallerNeighbours. A HeavierToSmallerNeighbours is an object 
// is used to partition a PartitionedDomain.
//
// What: "@(#) HeavierToSmallerNeighbours.C, revA"

#include <SwapHeavierToLighterNeighbours.h>
#include <Graph.h>
#include <VertexIter.h>
#include <Vertex.h>
#include <ID.h> 


SwapHeavierToLighterNeighbours::SwapHeavierToLighterNeighbours()
 :numReleases(1), factorGreater(1.0)
{

}

SwapHeavierToLighterNeighbours::
SwapHeavierToLighterNeighbours(double factGreater, 
			       int releases)
 :numReleases(releases),factorGreater(factGreater)
{
    if (releases < 1)
	numReleases = 1;
}


SwapHeavierToLighterNeighbours::~SwapHeavierToLighterNeighbours()
{
    
}

int
SwapHeavierToLighterNeighbours::balance(Graph &theWeightedGraph)
{
    // check to see a domain partitioner has been set
    DomainPartitioner *thePartitioner = this->getDomainPartitioner();
    if (thePartitioner == 0) {
	opserr << "SwapHeavierToLighterNeighbours::balance";
	opserr << "- No DomainPartitioner has been set\n"; 
	return -1;
    }

    int res = 0;

    for (int ii=0; ii<numReleases; ii++) {
	VertexIter &theVertices = theWeightedGraph.getVertices();
	Vertex *vertexPtr;
	while ((vertexPtr = theVertices()) != 0) {
	    int vertexTag = vertexPtr->getTag();
	    double vertexLoad = vertexPtr->getWeight();
	    const ID &adjacency = vertexPtr->getAdjacency();
	    int size = adjacency.Size();
	    for (int j=0; j<size; j++) {
		int otherVertexTag = adjacency(j);
		Vertex *otherVertexPtr 
		    = theWeightedGraph.getVertexPtr(otherVertexTag);
		double otherVertexLoad = otherVertexPtr->getWeight();
		
		if (vertexLoad > otherVertexLoad && otherVertexLoad != 0) 
		    if (vertexLoad/otherVertexLoad > factorGreater) {
			res = thePartitioner->
			    swapBoundary(vertexTag,otherVertexTag);
			if (res < 0) {
			    opserr << "WARNING SwapHeavierToLighterNeighbours";
			    opserr << "::balance - DomainPartitioner returned ";
			    opserr << res << endln;
			    return res;
			}
		    }
		
		if (vertexLoad != 0 && otherVertexLoad == 0)  {
		    res = thePartitioner->
			swapBoundary(vertexTag,otherVertexTag); 
		    if (res < 0) {
			opserr << "WARNING SwapHeavierToLighterNeighbours";
			opserr << "::balance - DomainPartitioner returned ";
			opserr << res << endln;
			return res;
		    }
		}
	    }		
	}
    }


    return res;

}
