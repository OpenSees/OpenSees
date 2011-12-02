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
                                                                        
// $Revision: 1.1 $
// $Date: 2005-11-29 21:55:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/numberer/ParallelNumberer.h,v $
                                                                        
// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for ParallelNumberer.
// ParallelNumberer is a subclass of DOF_Numberer. The ParallelNumberer numbers
// the dof of a partitioned domain, where the partitions are on different processors
// and each processor has a ParallelNumberer. The ParalellNumberer sitting on P0, 
// collects each partition graph from P1 through Pn-1, merges them into 1 large graph, 
// & then numbers this graph. The ParallelNumberers sitting on P1 through Pn-1 then 
// receive the mapping info for the dof tag and dof numbering from P0.
//
// What: "@(#) ParallelNumberer.h, revA"

#ifndef ParallelNumberer_h
#define ParallelNumberer_h

class Graph;

#include <DOF_Numberer.h>

class ParallelNumberer: public DOF_Numberer
{
  public:
    ParallelNumberer(int domainTag, int numSubdomains, Channel **theChannel);
    ParallelNumberer(GraphNumberer &theGraphNumberer);
    ParallelNumberer();

    ~ParallelNumberer();

    int numberDOF(int lastDOF = -1);
    int numberDOF(ID &lastDOFs);    

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    virtual int setProcessID(int domainTag);
    virtual int setChannels(int numChannels, Channel **theChannels);

  protected:
    int mergeSubGraph(Graph &theGraph, Graph &theSubGraph, ID &vertexTags, ID &vertexRefs, ID &theSubdomainMap);

  private:
    GraphNumberer *theNumberer;

    int processID;
    int numChannels;
    Channel **theChannels;
};

#endif

