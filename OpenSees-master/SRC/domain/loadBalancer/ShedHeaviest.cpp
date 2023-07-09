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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/loadBalancer/ShedHeaviest.cpp,v $
                                                                        
                                                                        
 // File: ~/domain/loadBalancer/ShedHeaviest.C
// 
// Written: fmk 
// Created: Fri Aug 29 17:43:25 1997
// Revision: A
//
// Description: This file contains the class definition for ShedHeaviest.
// A ShedHeaviest is an object used to partition a PartitionedDomain.
//
// What: "@(#) ShedHeaviest.C, revA"

#include <ShedHeaviest.h>
#include <Graph.h>
#include <VertexIter.h>
#include <Vertex.h>
 
ShedHeaviest::ShedHeaviest()
 :numReleases(1),factorGreater(1.0),disallowDisconnectedGraphs(true)
{
    
}

ShedHeaviest::ShedHeaviest(double fact, int releases, bool disallowDisconnected)
 :numReleases(releases),
  factorGreater(fact),
  disallowDisconnectedGraphs(disallowDisconnected)
{
    if (releases < 0)
	numReleases = 0;
}

ShedHeaviest::~ShedHeaviest()
{
    
}

int
ShedHeaviest::balance(Graph &theWeightedGraph)
{
    // check to see a domain partitioner has been set
    DomainPartitioner *thePartitioner = this->getDomainPartitioner();
    if (thePartitioner == 0) {
	opserr << "ShedHeaviest::balance - No DomainPartitioner has been set\n";
	return -1;
    }

    // determine the max loaded partition
    VertexIter &theVertices = theWeightedGraph.getVertices();
    Vertex *vertexPtr = theVertices();
    int maxPartition = vertexPtr->getTag();
    double maxLoad = vertexPtr->getWeight();
    while ((vertexPtr = theVertices()) != 0)
	if (vertexPtr->getWeight() > maxLoad) {
	    maxLoad = vertexPtr->getWeight();
	    maxPartition = vertexPtr->getTag();
	}
    
    // release the boundary numReleases times
    int res = 0;
    for (int j=0; j<numReleases; j++) {
	res = thePartitioner->
	    releaseBoundary(maxPartition,theWeightedGraph,
			    true,factorGreater);  
	
	if (res < 0) {
	    opserr << "WARNING ShedHeaviest::balance() ";
	    opserr << " - DomainPartitioner::releaseBoundary returned ";
	    opserr << res << endln;
	    j = numReleases;
	}
    }
    
    return res;
}
