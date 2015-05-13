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
                                                                        
// $Revision: 1.9 $
// $Date: 2010-04-23 22:50:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/constraints/ImposedMotionSP.cpp,v $
                                                                        
// Written: fmk 
// Created: 11/00
// Revision: A
//
// Purpose: This file contains the implementation of class ImposedMotionSP.

#include <ImposedMotionSP.h>
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
ImposedMotionSP::ImposedMotionSP()
:SP_Constraint(CNSTRNT_TAG_ImposedMotionSP),
 groundMotionTag(0), patternTag(0),
 theGroundMotion(0), theNode(0), theNodeResponse(0), theGroundMotionResponse(3)

{
    // does nothing else
}

// constructor for a subclass to use
ImposedMotionSP::ImposedMotionSP(int node, int ndof, int pattern, int motion)
:SP_Constraint(node, ndof, CNSTRNT_TAG_ImposedMotionSP),
 groundMotionTag(motion), patternTag(pattern),
 theGroundMotion(0), theNode(0), theNodeResponse(0), theGroundMotionResponse(3)
{

}


ImposedMotionSP::~ImposedMotionSP()
{
  if (theNodeResponse != 0)
    delete theNodeResponse;
}



double
ImposedMotionSP::getValue(void)
{
  // no longer return 0 for TransformationConstraints also set response at nodes
    return theGroundMotionResponse(0);    
}


int
ImposedMotionSP::applyConstraint(double time)
{
  // on first 
  if (theGroundMotion == 0 || theNode == 0 || theNodeResponse == 0) {
    Domain *theDomain = this->getDomain();

    theNode = theDomain->getNode(nodeTag);
    if (theNode == 0) {
      opserr << "ImposedMotionSP::applyConstraint() - node " << nodeTag << " does not exist\n";
      return -1;
    }

    int numNodeDOF = theNode->getNumberDOF();

    if (dofNumber < 0 || numNodeDOF <= dofNumber) {
      opserr << "ImposedMotionSP::applyConstraint() - dof number " << dofNumber++ << " at node " << nodeTag << " not valid\n";
      return -2;
    }

    theNodeResponse = new Vector(numNodeDOF);
    if (theNodeResponse == 0) {
      opserr << "ImposedMotionSP::applyConstraint() - out of memory\n";
      return -2;
    }

    LoadPattern *theLoadPattern = theDomain->getLoadPattern(patternTag);
    if (theLoadPattern == 0) {
      opserr << "ImposedMotionSP::applyConstraint() - no load pattern\n";
      return -3;
    }

    theGroundMotion = theLoadPattern->getMotion(groundMotionTag);
    if (theGroundMotion == 0) {
      opserr << "ImposedMotionSP::applyConstraint() - no ground motion\n";
      return -4;
    }
  }

  if (theNodeResponse == 0) 
    return -1;
  
  // now get the response from the ground motion
  theGroundMotionResponse = theGroundMotion->getDispVelAccel(time);
  
  //
  // now set the responses at the node
  //
  
  /* ***********************************************************
   * disp response the responsibility of constraint handler
   
   *theNodeResponse = theNode->getTrialDisp();
   (*theNodeResponse)(dofNumber) = theGroundMotionResponse(0);
   theNode->setTrialDisp(*theNodeResponse);
  *************************************************************/
  
  *theNodeResponse = theNode->getTrialVel();
  (*theNodeResponse)(dofNumber) = theGroundMotionResponse(1);
  theNode->setTrialVel(*theNodeResponse);    
  
  *theNodeResponse = theNode->getTrialAccel();
  (*theNodeResponse)(dofNumber) = theGroundMotionResponse(2);
  theNode->setTrialAccel(*theNodeResponse);        

  return 0;
}


bool
ImposedMotionSP::isHomogeneous(void) const
{
  return false;
}


int 
ImposedMotionSP::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  int result = 0;
  result = this->SP_Constraint::sendSelf(cTag, theChannel);
  if (result < 0) {
    opserr << "ImposedMotionSP::sendSelf() - base SP_Constraint class failed\n";
    return -1;
  }
  
  static ID myExtraData(2);
  myExtraData(0) = groundMotionTag;
  myExtraData(1) = patternTag;
  if (theChannel.sendID(dbTag, cTag, myExtraData) < 0) {
    opserr << "ImposedMotionSP::sendSelf() - failed to send extra data\n";
    return -1;
  }

  return 0;
}

int 
ImposedMotionSP::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();
  int result = 0;
  result = this->SP_Constraint::recvSelf(cTag, theChannel, theBroker);
  if (result < 0) {
    opserr << "ImposedMotionSP::recvSelf() - base SP_Constraint class failed\n";
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
ImposedMotionSP::Print(OPS_Stream &s, int flag) 
{
    s << "ImposedMotionSP: " << this->getTag();
    s << "\t Node: " << this->getNodeTag();
    s << " DOF: " << this->getDOF_Number() << endln;
}








