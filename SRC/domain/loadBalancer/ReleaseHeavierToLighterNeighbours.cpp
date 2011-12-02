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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/loadBalancer/ReleaseHeavierToLighterNeighbours.cpp,v $
                                                                        
                                                                        
 // File: ~/domain/loadBalancer/ReleaseHeavierToLighterNeighbours.C
// 
// Written: fmk 
// Created: Fri Aug 29 17:43:25 1997
// Revision: A
//
// Description: This file contains the class definition for ReleaseHeavierToLighterNeighbours.
// A ReleaseHeavierToLighterNeighbours is an object used to partition a PartitionedDomain.
//
// What: "@(#) ReleaseHeavierToLighterNeighbours.C, revA"

#include <ReleaseHeavierToLighterNeighbours.h>
#include <Graph.h>
#include <VertexIter.h>
#include <Vertex.h>
 
ReleaseHeavierToLighterNeighbours::ReleaseHeavierToLighterNeighbours()
 :numReleases(1), factorGreater(1.0), disallowDisconnectedGraphs(true)
{
    
    
}

ReleaseHeavierToLighterNeighbours::
ReleaseHeavierToLighterNeighbours(double factGreater, 
				  int releases,
				  bool disallowDisconnected)
 :numReleases(releases),factorGreater(factGreater),
  disallowDisconnectedGraphs(disallowDisconnected)  
{
    if (releases < 0)
	numReleases = 0;
}

ReleaseHeavierToLighterNeighbours::~ReleaseHeavierToLighterNeighbours()
{
    
}

int
ReleaseHeavierToLighterNeighbours::balance(Graph &theWeightedGraph)
{
    // check to see a domain partitioner has been set
    DomainPartitioner *thePartitioner = this->getDomainPartitioner();    
    if (thePartitioner == 0) {
	opserr << "ReleaseHeavierToLighterNeighbours::balance ";
	opserr << "- No DomainPartitioner has been set\n";
	return -1;
    }

    int res = 0;
    int numPartitions = theWeightedGraph.getNumVertex();
    for (int i=1; i<=numPartitions; i++) {
	res = thePartitioner->
	    releaseBoundary(i,theWeightedGraph,
			    true,
			    factorGreater);  
	if (res < 0) {
	    opserr << "WARNING ReleaseHeavierToLighterNeighbours";
	    opserr << "::balance - DomainPartitioner returned ";
	    opserr << res << endln;
	    return res;    
	}
    }
    
    return res;
}
