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
// $Source: /usr/local/cvs/OpenSees/SRC/graph/numberer/MyRCM.h,v $
                                                                        
                                                                        
// File: ~/graph/numberer/MyRCM.h
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the class definition for MyRCM.
// MyRCM is an object to perform the Reverse Cuthill-McKee numbering
// scheme on the vertices of a graph. This is done by invoking the
// number() method with the Graph to be numbered.
//
// Side effects: numberer() changes the Tmp values of the vertices to
// the number assigned to that vertex.
//
// What: "@(#) MyRCM.h, revA"

#ifndef MyRCM_h
#define MyRCM_h

#include <GraphNumberer.h>

#ifndef _bool_h
#include <bool.h>
#endif

class MyRCM: public GraphNumberer
{
  public:
    MyRCM(int startVertex = -1, bool minDegreeFlag = false); 
    ~MyRCM();

    void setStartVertex(int startVertex);
    void setMinDegreeFlag(bool flag);    

    const ID &number(Graph &theGraph, int lastVertex = -1);
    const ID &number(Graph &theGraph, const ID &lastVertices);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
  protected:
    
  private:
    int numVertex;
    ID *theRefResult;
    int startVertexTag;
    bool minDegree;
};

#endif

