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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/graph/FE_VertexIter.h,v $
                                                                        
                                                                        
// File: ~/analysis/model/simple/FE_VertexIter.h
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for FE_VertexIter.
// FE_VertexIter is an iter for returning the vertices of an object of class
// FE_Graph. FE_VertexIter must be written for each subclass of FE_Graph:
// wherin the vertices are stored differently to that in FE_Graph.

#ifndef FE_VertexIter_h
#define FE_VertexIter_h

#include <VertexIter.h>

class FE_Graph;

class FE_VertexIter: public VertexIter
{
  public:
    FE_VertexIter(FE_Graph &theFE_Graph);
    virtual ~FE_VertexIter();
    
    virtual void reset(void);
    virtual Vertex *operator()(void);
    
  private:
    FE_Graph &myGraph;
    int currIndex;
    int numDone;
};

#endif





