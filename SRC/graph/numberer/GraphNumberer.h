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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/numberer/GraphNumberer.h,v $
                                                                        
                                                                        
// File: ~/graph/numberer/GraphNumberer.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for GraphNumberer.
// GraphNumberer is an abstract base class. Its subtypes are responsible for
// numbering the vertices of a graph.
//
// What: "@(#) GraphNumberer.h, revA"

#ifndef GraphNumberer_h
#define GraphNumberer_h
#include <MovableObject.h>

class ID;
class Graph;
class Channel;
class ObjectBroker;

class GraphNumberer : public MovableObject
{
  public:
    GraphNumberer(int classTag);
    virtual ~GraphNumberer();
    
    virtual const ID &number(Graph &theGraph, int lastVertex = -1) =0;
    virtual const ID &number(Graph &theGraph, const ID &lastVertices) =0;
    
  protected:
    
  private:
    
};

#endif

