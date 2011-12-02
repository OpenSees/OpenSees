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
// $Date: 2006-05-26 18:26:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/numberer/RCM.h,v $
                                                                        
                                                                        
// File: ~/graph/numberer/RCM.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for RCM.
// RCM is an object to perform the Reverse Cuthill-McKee numbering
// scheme on the vertices of a graph. This is done by invoking the
// number() method with the Graph to be numbered.
//
// Side effects: numberer() changes the Tmp values of the vertices to
// the number assigned to that vertex.
//
// What: "@(#) RCM.h, revA"

#ifndef RCM_h
#define RCM_h

#include <GraphNumberer.h>

#ifndef _bool_h
#include <bool.h>
#endif

class RCM: public GraphNumberer
{
  public:
    RCM(bool GPS = false); 
    ~RCM();

    const ID &number(Graph &theGraph, int lastVertex = -1);
    const ID &number(Graph &theGraph, const ID &lastVertices);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:
    
  private:
    
    int numVertex;
    ID *theRefResult;
    bool GPS; // flag for gibbs-poole-stodlymer
};

#endif

