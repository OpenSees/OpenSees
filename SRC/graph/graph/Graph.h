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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/Graph.h,v $
                                                                        
                                                                        
#ifndef Graph_h
#define Graph_h

// File: ~/graph/graph/Graph.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for Graph.
// The Graph class provides the abstraction of a graph, a collection of
// vertices and edges. The Graph class is a container class which stores
// and provides access to Vertex objects. The Vertices contain information 
// about the edges in this design.
//
// What: "@(#) Graph.h, revA"

#ifndef _bool_h
#include "bool.h"
#endif

#include <OPS_Stream.h>

class Vertex;
class VertexIter;
class TaggedObjectStorage;


class Graph
{
  public:
    Graph();
    Graph(int numVertices);    
    Graph(TaggedObjectStorage &theVerticesStorage);
    virtual ~Graph();

    virtual bool addVertex(Vertex *vertexPtr, bool checkAdjacency = true);
    virtual int addEdge(int vertexTag, int otherVertexTag);
    
    virtual Vertex *getVertexPtr(int vertexTag);
    virtual VertexIter &getVertices(void);
    virtual int getNumVertex(void) const;
    virtual int getNumEdge(void) const;
    virtual Vertex *removeVertex(int tag, bool removeEdgeFlag = true);
    
    virtual void Print(OPS_Stream &s, int flag =0);
    friend OPS_Stream &operator<<(OPS_Stream &s, Graph &M);    
    
  protected:
    
  private:
    TaggedObjectStorage *myVertices;
    VertexIter *theVertexIter;
    int numEdge;
};

#endif

