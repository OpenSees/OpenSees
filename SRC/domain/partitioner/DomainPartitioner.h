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
// $Date: 2006-01-10 00:31:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/partitioner/DomainPartitioner.h,v $
                                                                        
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for DomainPartitioner.
// A DomainPartitioner is an object used to partition a PartitionedDomain.
//
// What: "@(#) DomainPartitioner.h, revA"

#ifndef DomainPartitioner_h
#define DomainPartitioner_h

#ifndef _bool_h
#include <bool.h>
#endif

#include <ID.h>

class GraphPartitioner;
class LoadBalancer;
class PartitionedDomain;
class Vector;
class Graph;
class TaggedObjectStorage;

class DomainPartitioner
{
  public:
    
    DomainPartitioner(GraphPartitioner &theGraphPartitioner,
    		      LoadBalancer &theLoadBalancer);
    
    DomainPartitioner(GraphPartitioner &theGraphPartitioner);
    virtual  ~DomainPartitioner();    

    virtual void setPartitionedDomain(PartitionedDomain &theDomain);
    virtual int partition(int numParts, bool useMainDomain = false, int mainPartition = 0, int specialElementTag = 0);

    virtual int balance(Graph &theWeightedSubdomainGraph);

    // public member functions needed by the load balancer
    virtual int getNumPartitions(void) const;
    virtual Graph &getPartitionGraph(void);
    virtual Graph &getColoredGraph(void);
    
    virtual int  swapVertex(int from, 
			    int to, 
			    int vertexTag,
			    bool adjacentVertexNotInOther = true);

    virtual int	 swapBoundary(int from, 
			      int to,
			      bool adjacentVertexNotInOther = true);

    virtual int  releaseVertex(int from, 
			       int vertexTag, 
			       Graph &theWeightedPartitionGraph, 
			       bool mustReleaseToLighter = true,
			       double factorGreater = 1.0,
			       bool adjacentVertexNotInOther = true);

    virtual int releaseBoundary(int from, 
			       Graph &theWeightedPartitionGraph, 
			       bool mustReleaseToLighter = true,
			       double factorGreater = 1.0,
			       bool adjacentVertexNotInOther = true);
				 

    virtual GraphPartitioner* getGraphPartitioner();
  protected:    
    
  private:
    PartitionedDomain *myDomain; 
    GraphPartitioner  &thePartitioner;
    LoadBalancer      *theBalancer;    

    Graph *theElementGraph;
    Graph **theBoundaryElements; 
    
    TaggedObjectStorage *theNodeLocations;
    ID *elementPlace;
    int numPartitions;
    ID primes;
    bool partitionFlag;
    
    bool usingMainDomain;
    int mainPartition;
};

#endif


