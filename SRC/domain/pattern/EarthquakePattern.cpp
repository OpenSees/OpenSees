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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/EarthquakePattern.cpp,v $
                                                                        
                                                                        
// File: ~/domain/pattern/EarthquakePattern.cpp
//
// Written: fmk 11/98
// Revised:
//
// Purpose: This file contains the class definition for EarthquakePattern.
// EarthquakePattern is an abstract class.

#include <EarthquakePattern.h>
#include <GroundMotion.h>

#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <stdlib.h>
#include <Channel.h>
#include <ErrorHandler.h>

EarthquakePattern::EarthquakePattern(int tag, int _classTag)
:LoadPattern(tag, _classTag), theMotions(0), numMotions(0), uDotG(0), uDotDotG(0)
{

}


EarthquakePattern::~EarthquakePattern()
{
  // invoke the destructor on all ground motions supplied
  for (int i=0; i<numMotions; i++)
      delete theMotions[i];
  
  if (theMotions != 0)
    //    free ((void *)theMotions);
    delete [] theMotions;
}



void 
EarthquakePattern::applyLoad(double time)
{
  // see if quick return, i.e. no Ground Motions or domain set
  if (numMotions == 0)
    return;

  Domain *theDomain = this->getDomain();
  if (theDomain == 0)
    return;


  // set the vel and accel vector
  for (int i=0; i<numMotions; i++) {
    (*uDotG)(i) = theMotions[i]->getVel(time);
    (*uDotDotG)(i) = theMotions[i]->getAccel(time);
  }

  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  while ((theNode = theNodes()) != 0) 
    theNode->addInertiaLoadToUnbalance(*uDotDotG, 1.0);
  

  ElementIter &theElements = theDomain->getElements();
  Element *theElement;
  while ((theElement = theElements()) != 0) 
    theElement->addInertiaLoadToUnbalance(*uDotDotG);
}
    
int
EarthquakePattern::addMotion(GroundMotion &theMotion)
{
  // make space for new
  GroundMotion **newMotions = new GroundMotion *[numMotions+1];
  //  GroundMotion **newMotions = (GroundMotion **)malloc(sizeof(GroundMotion *)*(numMotions+1));
  if (newMotions == 0) {
    cerr << "EarthquakePattern::addMotion - could not add new, out of mem\n";
    return -1;
  }
  
  // copy old
  for (int i=0; i<numMotions; i++)
    newMotions[i] = theMotions[i];

  // add the new motion to new
  newMotions[numMotions] = &theMotion;

  // delete the old
  if (theMotions != 0)
    delete [] theMotions;

  // reset
  theMotions = newMotions;
  numMotions++;

  // delete old vectors and create new ones
  if (uDotG != 0)
    delete uDotG;
  uDotG = new Vector(numMotions);

  if (uDotDotG != 0)
    delete uDotDotG;
  uDotDotG = new Vector(numMotions);

  if (uDotDotG == 0 || uDotDotG->Size() == 0 || uDotG == 0 || uDotG->Size() == 0) {
    cerr << "EarthquakePattern::addMotion - ran out of memory creating vectors\n";
    numMotions = 0;
    return -2;
  }
  return 0;
}


bool
EarthquakePattern::addSP_Constraint(SP_Constraint *)
{
  cerr << "EarthquakePattern::addSP_Constraint() - cannot add SP_Constraint to EQ pattern\n";
  return false;
}

bool
EarthquakePattern::addNodalLoad(NodalLoad *)
{
  cerr << "EarthquakePattern::addNodalLoad() - cannot add NodalLoad to EQ pattern\n";  
  return false;
}

bool
EarthquakePattern::addElementalLoad(ElementalLoad *)
{
  cerr << "EarthquakePattern::addElementalLoad() - cannot add ElementalLoad to EQ pattern\n";    
  return false;
}


/* **********************************************************************************************
int 
EarthquakePattern::sendSelf(int commitTag, Channel &theChannel)
{
  // first send the tag and info about the number of ground motions
  int myDbTag = this->getDbTag();
  ID theData(2);
  theData(0) = this->getTag();
  theData(1) = numMotions;

  if (theChannel.sendID(myDbTag, commitTag, theData) < 0) {
    g3ErrorHandler->warning("EarthquakePattern::sendSelf - channel failed to send the initial ID");
    return -1;
  }    

  // now for each motion we send it's classsss tag and dbtag
  ID theMotionsData(2*numMotions);
  for (int i=0; i<numMotions; i++) {
    theMotionsData[i] = theMotions[i]->getClassTag();
    int motionsDbTag = theMotions[i]->getDbTag();
    if (motionsDbTag == 0) {
      motionsDbTag = theChannel.getDbTag();
      theMotions[i]->setDbTag(motionsDbTag);
    }
    theMotionsData[i+numMotions] = motionsDbTag;
  }

  if (theChannel.sendID(myDbTag, commitTag, theMotionsData) < 0) {
    g3ErrorHandler->warning("EarthquakePattern::sendSelf - channel failed to send the motions ID");
    return -1;
  }    

  // now we send each motion
  for (int j=0; j<numMotions; j++)
    if (theMotions[j]->sendSelf(commitTag, theChannel) < 0) {
      g3ErrorHandler->warning("EarthquakePattern::sendSelf - motion no: %d failed in sendSelf", j);
      return -1;
    }

  // if get here successfull
  return 0;
}

int 
EarthquakePattern::recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  // first get the tag and info about the number of ground motions from the Channel
  int myDbTag = this->getDbTag();
  ID theData(2);
  if (theChannel.recvID(myDbTag, commitTag, theData) < 0) {
    g3ErrorHandler->warning("EarthquakePattern::recvSelf - channel failed to recv the initial ID");
    return -1;
  }    
  
  // set current tag
  this->setTag(theData(0));

  // now get info about each channel
  ID theMotionsData (2*theData(1));
  if (theChannel.recvID(myDbTag, commitTag, theMotionsData) < 0) {
    g3ErrorHandler->warning("EarthquakePattern::recvSelf - channel failed to recv the motions ID");
    return -1;
  }    


  if (numMotions != theData(1)) {

    //
    // we must delete the old motions and create new ones and then invoke recvSelf on these new ones
    //

    if (numMotions != 0) {
      for (int i=0; i<numMotions; i++)
	delete theMotions[i];
      delete [] theMotions;
    }
    numMotions = theData[1];
    theMotions = new (GroundMotion *)[numMotions];
    if (theMotions == 0) {
      g3ErrorHandler->warning("EarthquakePattern::recvSelf - out of memory creating motion array of size %d\n", numMotions);
      numMotions = 0;
      return -1;
    }    

    for (int i=0; i<numMotions; i++) {
      theMotions[i] = theBroker.getNewGroundMotion(theMotionsData[i]);
      if (theMotions[i] == 0) {
	g3ErrorHandler->warning("EarthquakePattern::recvSelf - out of memory creating motion array of size %d\n", numMotions);
	numMotions = 0;
	return -1;
      }    
      theMotions[i]->setDbTag(theMotionsData[i+numMotions]);
      if (theMotions[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	g3ErrorHandler->warning("EarthquakePattern::recvSelf - motion no: %d failed in recvSelf", i);
	numMotions = 0;
	return -1;
      }      
    }

  } else {

    // 
    // we invoke rrecvSelf on the motions, note if a motion not of correct class
    // we must invoke the destructor on the motion and create a new one of correct type
    //

    for (int i=0; i<numMotions; i++) {
      if (theMotions[i]->getClassTag() == theMotionsData[i]) {
	if (theMotions[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  g3ErrorHandler->warning("EarthquakePattern::recvSelf - motion no: %d failed in recvSelf", i);
	  return -1;
	}      
      } else {
	// motion not of correct type
	delete theMotions[i];
	theMotions[i] = theBroker.getNewGroundMotion(theMotionsData[i]);
	if (theMotions[i] == 0) {
	  g3ErrorHandler->warning("EarthquakePattern::recvSelf - out of memory creating motion array of size %d\n", numMotions);
	  numMotions = 0;
	  return -1;
	}    
	theMotions[i]->setDbTag(theMotionsData[i+numMotions]);
	if (theMotions[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	  g3ErrorHandler->warning("EarthquakePattern::recvSelf - motion no: %d failed in recvSelf", i);
	  numMotions = 0;
	  return -1;
	}      
      }
    }    
  }

  // if get here successfull
  return 0;
}


***************************************************************************************** */

void 
EarthquakePattern::Print(ostream &s, int flag)
{
  cerr << "EarthquakePattern::Print() - not yet implemented\n";    
}

/* ****************************************************************************************
// method to obtain a blank copy of the LoadPattern
LoadPattern *
EarthquakePattern::getCopy(void)
{
  EarthquakePattern *theCopy = new EarthquakePattern(0, 0);
  theCopy->setTag(this->getTag);
  theCopy->numMotions = numMotions;
  theCopy->theMotions = new (GroundMotion *)[numMotions];
  for (int i=0; i<numMotions; i++)
    theCopy->theMotions[i] = theMotions[i];

  return 0;
}
***************************************************************************************** */
