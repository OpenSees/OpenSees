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
                                                                        
// $Revision: 1.7 $
// $Date: 2005-11-07 21:36:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/UniformExcitation.cpp,v $
                                                                        
                                                                        
// File: ~/domain/load/UniformExcitation.h
//
// Written: fmk 11/98
// Revised:
//
// Purpose: This file contains the class definition for UniformExcitation.
// UniformExcitation is an abstract class.

#include <UniformExcitation.h>
#include <GroundMotion.h>
#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

UniformExcitation::UniformExcitation()
:EarthquakePattern(0, PATTERN_TAG_UniformExcitation), 
  theMotion(0), theDof(0), vel0(0.0)
{

}


UniformExcitation::UniformExcitation(GroundMotion &_theMotion, 
				   int dof, int tag, double velZero)
:EarthquakePattern(tag, PATTERN_TAG_UniformExcitation), 
  theMotion(&_theMotion), theDof(dof), vel0(velZero)
{
  // add the motion to the list of ground motions
  this->addMotion(*theMotion);
}


UniformExcitation::~UniformExcitation()
{

}


void
UniformExcitation::setDomain(Domain *theDomain) 
{
  this->LoadPattern::setDomain(theDomain);

  // now we go through and set all the node velocities to be vel0
  if (vel0 != 0.0) {
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    Vector newVel(1);
    int currentSize = 1;
    while ((theNode = theNodes()) != 0) {
      int numDOF = theNode->getNumberDOF();
      if (numDOF != currentSize) 
	newVel.resize(numDOF);
      
      newVel = theNode->getVel();
      newVel(theDof) = vel0;
      theNode->setTrialVel(newVel);
      theNode->commitState();
    }
  }
}

void
UniformExcitation::applyLoad(double time)
{
  Domain *theDomain = this->getDomain();
  if (theDomain == 0)
    return;

//  if (numNodes != theDomain->getNumNodes()) {
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      theNode->setNumColR(1);
      theNode->setR(theDof, 0, 1.0);
    }
//  }

  this->EarthquakePattern::applyLoad(time);

  return;
}


void
UniformExcitation::applyLoadSensitivity(double time)
{
  Domain *theDomain = this->getDomain();
  if (theDomain == 0)
    return;

//  if (numNodes != theDomain->getNumNodes()) {
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      theNode->setNumColR(1);
      theNode->setR(theDof, 0, 1.0);
    }
//  }

  this->EarthquakePattern::applyLoadSensitivity(time);

  return;
}



int 
UniformExcitation::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector data(5);
  data(0) = this->getTag();
  data(1) = theDof;
  data(2) = vel0;
  data(3) = theMotion->getClassTag();
  
  int motionDbTag = theMotion->getDbTag();
  if (motionDbTag == 0) {
    motionDbTag = theChannel.getDbTag();
    theMotion->setDbTag(motionDbTag);
  }
  data(4) = motionDbTag;

  int res = theChannel.sendVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "UniformExcitation::sendSelf() - channel failed to send data\n";
    return res;
  }
      
  res = theMotion->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "UniformExcitation::sendSelf() - ground motion to send self\n";
    return res;
  }

  return 0;
}


int 
UniformExcitation::recvSelf(int commitTag, Channel &theChannel, 
			   FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector data(5);
  int res = theChannel.recvVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "UniformExcitation::recvSelf() - channel failed to recv data\n";
    return res;
  }

  this->setTag(data(0));
  theDof = data(1);
  vel0 = data(2);
  int motionClassTag = data(3);
  int motionDbTag = data(4);

  if (theMotion == 0 || theMotion->getClassTag() != motionClassTag) {
    if (theMotion != 0)
      delete theMotion;
    theMotion = theBroker.getNewGroundMotion(motionClassTag);
    if (theMotion == 0) {
      opserr << "UniformExcitation::recvSelf() - could not create a grond motion\n";
      return -3;
    }

    // have to set the motion in EarthquakePattern base class
    if (numMotions == 0) 
      this->addMotion(*theMotion);
    else
      theMotions[0] = theMotion;
  }

  theMotion->setDbTag(motionDbTag);
  res = theMotion->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
      opserr << "UniformExcitation::recvSelf() - motion could not receive itself \n";
      return res;
  }

  return 0;
}


void 
UniformExcitation::Print(OPS_Stream &s, int flag)
{
  s << "UniformExcitation  " << this->getTag() << " - Not Printing the GroundMotion\n";
}

LoadPattern *
UniformExcitation::getCopy(void)
{
  LoadPattern *theCopy = new UniformExcitation(*theMotion, theDof, this->getTag());
   return theCopy;
}
//  LocalWords:  OpenSees
