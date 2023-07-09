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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/VertexIter.h,v $
                                                                        
                                                                        
// File: ~/analysis/model/VertexIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for VertexIter.
// VertexIter is an abstract base class. A VertexIter is an iter 
// for returning the vertices of a graph.
// VertexIters must be written for each subclass of Graph.


#ifndef VertexIter_h
#define VertexIter_h

class Vertex;
class TaggedObjectStorage;
class TaggedObjectIter;

class VertexIter 
{
  public:
    VertexIter();    
    VertexIter(TaggedObjectStorage *);
    virtual ~VertexIter();

    virtual void reset(void);
    virtual Vertex *operator()(void);

  protected:
    
  private:
    TaggedObjectIter &myIter;
};

#endif

