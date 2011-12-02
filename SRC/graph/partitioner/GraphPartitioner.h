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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/partitioner/GraphPartitioner.h,v $
                                                                        
                                                                        
// File: ~/graph/partitioner/GraphPartitioner.h
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for GraphPartitioner.
// GraphPartitioner is an abstract base class. Its subtypes are responsible for
// partioning the vertices of a graph. The partitioning is done in the method
// partition which sets the colors of the vertices of the graph to colors 1
// through numParrtitions.
//
// What: "@(#) GraphPartitioner.h, revA"

#ifndef GraphPartitioner_h
#define GraphPartitioner_h

class ID;
class Graph;

class GraphPartitioner
{
  public:
    GraphPartitioner() {};
    virtual ~GraphPartitioner() {};
    
    virtual int partition(Graph &theGraph, int numPart) =0;

  protected:
    
  private:
    
};

#endif

