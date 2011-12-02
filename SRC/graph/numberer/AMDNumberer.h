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
                                                                        
// $Revision: 1.1 $
// $Date: 2009-12-10 00:40:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/numberer/AMDNumberer.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 12/09
// Revision: A
//
// Description: This file contains the class definition for AMD.
// AMD is an object to perform the Approx. Min Degree Ordering
//
// What: "@(#) AMD.h, revA"

#ifndef AMD_h
#define AMD_h

#include <GraphNumberer.h>

#ifndef _bool_h
#include <bool.h>
#endif

#include <ID.h>

class AMD: public GraphNumberer
{
  public:
    AMD(); 
    ~AMD();

    const ID &number(Graph &theGraph, int lastVertex = -1);
    const ID &number(Graph &theGraph, const ID &lastVertices);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:
    
  private:
    ID theResult;
};

#endif

