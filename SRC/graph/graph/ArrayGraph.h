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
// $Date: 2003-02-14 23:01:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/ArrayGraph.h,v $
                                                                        
                                                                        
// File: ~/graph/graph/ArrayGraph.h
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for ArrayGraph.
// The vertices in an ArrayGraph are held in an array. This is more efficient
// than holding them in a List data structure, but problems can arise with
// large Graphs in getting enough contiguous memory for the array.
//
// What: "@(#) ArrayGraph.h, revA"

#ifndef ArrayGraph_h
#define ArrayGraph_h

#include <Graph.h>
#include <ArrayVertexIter.h>

class ArrayGraph: public Graph
{
  public:
    ArrayGraph(int arraySize);
    virtual ~ArrayGraph();

    virtual bool addVertex(Vertex *vertexPtr);
    virtual int addEdge(int vertexTag, int otherVertexTag);
    
    virtual Vertex *getVertexPtr(int vertexTag);
    virtual VertexIter &getVertices(void);
    int getNumVertex(void) const;
    int getNumEdge(void) const;

    virtual void Print(OPS_Stream &s) const;
    friend OPS_Stream &operator<<(OPS_Stream &s, const ArrayGraph &M);    
    
    friend class ArrayVertexIter;    
    
  protected:
    int getArraySize(void) const;
    
  private:
    int numVertex;
    int numEdge;
    int sizeVertices;
    int lastEmpty;
    Vertex **theVertices;
    ArrayVertexIter myIter;

};

#endif

