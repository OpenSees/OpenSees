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
// $Date: 2009-12-10 00:40:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/numberer/AMDNumberer.cpp,v $
                                                                        

// Description: This file contains the class definition for AMD.

// What: "@(#) AMD.C, revA"

#include <amd.h>

#include <AMDNumberer.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

// Constructor
AMD::AMD()
:GraphNumberer(GraphNUMBERER_TAG_AMD)
{
  opserr << "AMD::AMD()\n";
}

// Destructor
AMD::~AMD()
{

}


// const ID &number(Graph &theGraph,int startVertexTag = -1,
//                  bool minDegree = false)
const ID &
AMD::number(Graph &theGraph, int startVertex)
{
  int numVertex = theGraph.getNumVertex();

  if (numVertex == 0) 
    return theResult;

  theResult.resize(numVertex);

  int nnz = 0;
  Vertex *vertexPtr;
  VertexIter &vertexIter = theGraph.getVertices();
  
  while ((vertexPtr = vertexIter()) != 0) {
    const ID &adjacency = vertexPtr->getAdjacency();
    nnz += adjacency.Size();
  }

  int *P = new int[numVertex];
  int *Ap = new int[numVertex+1];
  int *Ai = new int[nnz];
  double Control[AMD_CONTROL];
  double Info[AMD_INFO];

  VertexIter &vertexIter2 = theGraph.getVertices();

  nnz = 0;
  int count = 1;
  Ap[0] = 0;

  while ((vertexPtr = vertexIter2()) != 0) {
    const ID &adjacency = vertexPtr->getAdjacency();
    for (int i=0; i<adjacency.Size(); i++)
      Ai[nnz++] = adjacency(i);
    Ap[count++] = nnz;
  }  

  amd_order(numVertex, Ap, Ai, P, (double *)NULL, (double *)NULL);
  
  for (int i=0; i<numVertex; i++)
    theResult[i] = P[i];

  delete [] P;
  delete [] Ap;
  delete [] Ai;

  return theResult;
}



int
AMD::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
AMD::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}


const ID &
AMD::number(Graph &theGraph, const ID &startVertices)
{
  opserr << "WARNING:  AMD::number - Not implemented with startVertices";
  return theResult;
}






