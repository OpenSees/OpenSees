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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-02-29 20:47:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/UniformRotationalExcitation.cpp,v $
                                                                        
                                                                        
// File: ~/domain/load/UniformRotationalExcitation.h
//
// Written: fmk 11/98
// Revised:
//
// Purpose: This file contains the class definition for UniformRotationalExcitation.
// UniformRotationalExcitation is an abstract class.

#include <UniformRotationalExcitation.h>
#include <GroundMotion.h>
#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <SP_ConstraintIter.h>
#include <SP_Constraint.h>

UniformRotationalExcitation::UniformRotationalExcitation()
:EarthquakePattern(0, PATTERN_TAG_UniformRotationalExcitation), 
 theMotion(0), theDof(0), zCRD(0), vel0(0.0), fact(0.0)
{

}


UniformRotationalExcitation::UniformRotationalExcitation(GroundMotion &_theMotion, 
							 int dof, 
							 int tag, 
							 double velZero,
							 double theFactor)
:EarthquakePattern(tag, PATTERN_TAG_UniformRotationalExcitation), 
 theMotion(&_theMotion), theDof(dof), vel0(velZero), fact(theFactor)
{
  // add the motion to the list of ground motions
  this->addMotion(*theMotion);
}


UniformRotationalExcitation::~UniformRotationalExcitation()
{

}


const GroundMotion *
UniformRotationalExcitation::getGroundMotion(void)
{
  return theMotion;
}

int
UniformRotationalExcitation::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMotion->setParameter(argv, argc, param);
}

/*
int
UniformRotationalExcitation::updateParameter(int parameterID, Information &info)
{
  return theMotion->updateParameter(parameterID, info);
}

int
UniformRotationalExcitation::activateParameter(int pparameterID)
{
  return theMotion->activateParameter(pparameterID);
}
*/

void
UniformRotationalExcitation::setDomain(Domain *theDomain) 
{
  this->LoadPattern::setDomain(theDomain);

  // now we go through and set all the node velocities to be vel0 
  // for those nodes not fixed in the dirn!
  if (vel0 != 0.0) {

    SP_ConstraintIter &theSPs = theDomain->getSPs();
    SP_Constraint *theSP;
    ID constrainedNodes(0);
    int count = 0;
    while ((theSP=theSPs()) != 0) {
      if (theSP->getDOF_Number() == theDof) {
	constrainedNodes[count] = theSP->getNodeTag();
	count++;
      }
    }


    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    Vector newVel(1);
    int currentSize = 1;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      if (constrainedNodes.getLocation(tag) < 0) {
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
}

void
UniformRotationalExcitation::applyLoad(double time)
{
  Domain *theDomain = this->getDomain();
  if (theDomain == 0)
    return;

//  if (numNodes != theDomain->getNumNodes()) {
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      theNode->setNumColR(1);

      const Vector &crds=theNode->getCrds();
      int ndm = crds.Size();
      
      if (ndm == 2) {
	double xCrd = crds(0);
	double zCrd = crds(1);	  
	
	theNode->setR(0, 0, -fact*zCrd);	
	theNode->setR(1, 0, fact*xCrd);	
	theNode->setR(2, 0, fact);	
	
      } else if (ndm == 3 && theDof < 2) {
	
	if (theDof == 0) {
	  double xCrd = crds(0);
	  double zCrd = crds(2);	  
	  
	  theNode->setR(0, 0, fact*zCrd);	
	  theNode->setR(2, 0, -fact*xCrd);
	  theNode->setR(4, 0, fact);
	  
	} else if (theDof == 1) {
	  
	  double yCrd = crds(1);
	  double zCrd = crds(2);	  
	  
	  theNode->setR(1, 0, -fact*zCrd);	
	  theNode->setR(2, 0, fact*yCrd);
	  theNode->setR(3, 0, fact);
	}
      }
    }

  this->EarthquakePattern::applyLoad(time);

  return;
}


void
UniformRotationalExcitation::applyLoadSensitivity(double time)
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
UniformRotationalExcitation::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector data(6);
  data(0) = this->getTag();
  data(1) = theDof;
  data(2) = vel0;
  data(5) = fact;
  data(3) = theMotion->getClassTag();
  
  int motionDbTag = theMotion->getDbTag();
  if (motionDbTag == 0) {
    motionDbTag = theChannel.getDbTag();
    theMotion->setDbTag(motionDbTag);
  }
  data(4) = motionDbTag;

  int res = theChannel.sendVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "UniformRotationalExcitation::sendSelf() - channel failed to send data\n";
    return res;
  }
      
  res = theMotion->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "UniformRotationalExcitation::sendSelf() - ground motion to send self\n";
    return res;
  }

  return 0;
}


int 
UniformRotationalExcitation::recvSelf(int commitTag, Channel &theChannel, 
			   FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector data(6);
  int res = theChannel.recvVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "UniformRotationalExcitation::recvSelf() - channel failed to recv data\n";
    return res;
  }

  this->setTag(data(0));
  theDof = data(1);
  vel0 = data(2);
  fact = data(5);
  int motionClassTag = data(3);
  int motionDbTag = data(4);

  if (theMotion == 0 || theMotion->getClassTag() != motionClassTag) {
    if (theMotion != 0)
      delete theMotion;
    theMotion = theBroker.getNewGroundMotion(motionClassTag);
    if (theMotion == 0) {
      opserr << "UniformRotationalExcitation::recvSelf() - could not create a grond motion\n";
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
      opserr << "UniformRotationalExcitation::recvSelf() - motion could not receive itself \n";
      return res;
  }

  return 0;
}


void 
UniformRotationalExcitation::Print(OPS_Stream &s, int flag)
{
  s << "UniformRotationalExcitation  " << this->getTag() << " - Not Printing the GroundMotion\n";
}

LoadPattern *
UniformRotationalExcitation::getCopy(void)
{
  LoadPattern *theCopy = new UniformRotationalExcitation(*theMotion, theDof, this->getTag(), vel0, fact);
  return theCopy;
}
//  LocalWords:  OpenSees
