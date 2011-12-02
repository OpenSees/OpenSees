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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/ArrayVertexIter.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/model/simple/ArrayVertexIter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// ArrayVertexIter. ArrayVertexIter is a class for iterating through the 
// vertices of an ArrayGraph. 


#include <ArrayGraph.h>
#include <ArrayVertexIter.h>
#include <Vertex.h>

// ArrayVertexIter(ArrayGraph &theGraph):
//	constructor that takes the graph, just the basic iter

ArrayVertexIter::ArrayVertexIter(ArrayGraph &theGraph)
  :myGraph(theGraph), currIndex(0), numDone(0)
{
}


ArrayVertexIter::~ArrayVertexIter()
{
}    

void
ArrayVertexIter::reset(void)
{
    currIndex = 0;
    numDone = 0;
}

Vertex *
ArrayVertexIter::operator()(void)
{
  // check if we still have vertices in the Graph
  // if not return 0, indicating we are done
  if (numDone >= myGraph.numVertex)
    return 0;

  // search through domains ele list till we find the next element
  while ((currIndex < myGraph.sizeVertices) 
	 && (myGraph.theVertices[currIndex] == 0))
      currIndex++;

  // if not at the end of the list return the element
  if (currIndex < myGraph.sizeVertices) {
      Vertex *result = myGraph.theVertices[currIndex];
      numDone++; currIndex++;
      return(result);
  }
  return (0);
}

    
    




