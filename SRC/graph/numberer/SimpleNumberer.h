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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/numberer/SimpleNumberer.h,v $
                                                                        
                                                                        
// File: ~/graph/numberer/SimpleNumberer.h
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for SimpleNumberer.
// SimpleNumberer is an object to perform a simple numbering of the vertices.
// It does this by using the graphs VertexIter and assigning the numbers as
// it comes across the vertices.
//
// What: "@(#) SimpleNumberer.h, revA"

#ifndef SimpleNumberer_h
#define SimpleNumberer_h

#include <GraphNumberer.h>

class SimpleNumberer: public GraphNumberer
{
  public:
    SimpleNumberer(); 
    ~SimpleNumberer();
    
    const ID &number(Graph &theGraph, int lastVertex = -1);
    const ID &number(Graph &theGraph, const ID &startVertices);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
  protected:
    
  private:
    int numVertex;
    ID *theRefResult;
};

#endif

