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
// $Date: 2000-12-12 07:26:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/ImposedMotionSP1.cpp,v $
                                                                        
                                                                        
// File: ~/domain/constraints/ImposedMotionSP1.C
//
// Written: fmk 
// Created: 11/00
// Revision: A
//
// Purpose: This file contains the implementation of class ImposedMotionSP1.

#include <ImposedMotionSP1.h>
#include <classTags.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <GroundMotion.h>
#include <Node.h>
#include <Domain.h>

// constructor for FEM_ObjectBroker
ImposedMotionSP1::ImposedMotionSP1()
:SP_Constraint(CNSTRNT_TAG_ImposedMotionSP1),
 theGroundMotion(0), theNode(0), theGroundMotionResponse(3), destroyMotion(0)
{
    // does nothing else
}

// constructor for a subclass to use
ImposedMotionSP1::ImposedMotionSP1(int tag, int node, int ndof, 
				 GroundMotion &theMotion, bool killMotion)
:SP_Constraint(tag, node, ndof, CNSTRNT_TAG_ImposedMotionSP1),
 theNode(0), theGroundMotionResponse(3), destroyMotion(0)
{
  theGroundMotion = &theMotion;
  
  if (killMotion == true)
    destroyMotion = 1;
}


ImposedMotionSP1::~ImposedMotionSP1()
{
  if (destroyMotion == 1)
    delete theGroundMotion;
}



double
ImposedMotionSP1::getValue(void)
{
  // always return 0.0 - applyConstraint() sets the values at Node 
    return theGroundMotionResponse(0);
}


int
ImposedMotionSP1::applyConstraint(double time)
{
    // on first 
    if (theNode == 0) {
	Domain *theDomain = this->getDomain();

	theNode = theDomain->getNode(nodeTag);
	if (theNode == 0) {
	    
	    return -1;
	}
    }

    // now get the response from the ground motion
    theGroundMotionResponse = theGroundMotion->getDispVelAccel(time);

    return 0;
}


bool
ImposedMotionSP1::isHomogeneous(void) const
{
  return false;
}


int 
ImposedMotionSP1::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
ImposedMotionSP1::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
ImposedMotionSP1::Print(ostream &s, int flag) 
{
    s << "ImposedMotionSP1: " << this->getTag();
}








