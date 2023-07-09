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
                                                                        
// $Revision: 1.5 $
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/ImposedMotionSP1.cpp,v $
                                                                        
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
#include <LoadPattern.h>
#include <ID.h>

// constructor for FEM_ObjectBroker
ImposedMotionSP1::ImposedMotionSP1()
:SP_Constraint(CNSTRNT_TAG_ImposedMotionSP1),
 theGroundMotion(0), theNode(0), theGroundMotionResponse(3), destroyMotion(0)
{
    // does nothing else
}

// constructor for a subclass to use
ImposedMotionSP1::ImposedMotionSP1(int node, int ndof, int pattern, int motion)
:SP_Constraint(node, ndof, CNSTRNT_TAG_ImposedMotionSP1),
 groundMotionTag(motion), patternTag(pattern),
 theGroundMotion(0), theNode(0), theGroundMotionResponse(3)
{

}


ImposedMotionSP1::~ImposedMotionSP1()
{

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
  if (theGroundMotion == 0 || theNode == 0) {
    Domain *theDomain = this->getDomain();
    
    theNode = theDomain->getNode(nodeTag);
    if (theNode == 0) {
      
      return -1;
    }
    LoadPattern *theLoadPattern = theDomain->getLoadPattern(patternTag);
    if (theLoadPattern == 0)
      return -3;
    
    theGroundMotion = theLoadPattern->getMotion(groundMotionTag);
    if (theGroundMotion == 0)
      return -4;
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
  int dbTag = this->getDbTag();
  int result = 0;
  result = this->SP_Constraint::sendSelf(cTag, theChannel);
  if (result < 0) {
    opserr << "ImposedMotionSP1::sendSelf() - base SP_Constraint class failed\n";
    return -1;
  }
  
  static ID myExtraData(2);
  myExtraData(0) = groundMotionTag;
  myExtraData(1) = patternTag;
  if (theChannel.sendID(dbTag, cTag, myExtraData) < 0) {
    opserr << "ImposedMotionSP1::sendSelf() - failed to send extra data\n";
    return -1;
  }

  return 0;
}

int 
ImposedMotionSP1::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  int result = 0;
  result = this->SP_Constraint::recvSelf(cTag, theChannel, theBroker);
  if (result < 0) {
    opserr << "ImposedMotionSP1::recvSelf() - base SP_Constraint class failed\n";
    return -1;
  }
  
  static ID myExtraData(2);
  if (theChannel.recvID(dbTag, cTag, myExtraData) < 0) {
    opserr << "ImposedMotionSP::sendSelf() - failed to send extra data\n";
    return -1;
  }
  groundMotionTag = myExtraData(0);
  patternTag = myExtraData(1);

  return 0;
}

void
ImposedMotionSP1::Print(OPS_Stream &s, int flag) 
{
    s << "ImposedMotionSP1: " << this->getTag();
    s << "\t Node: " << this->getNodeTag();
    s << " DOF: " << this->getDOF_Number() << endln;    
}








