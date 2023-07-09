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
                                                                        
// $Revision: 1.4 $
// $Date: 2006-01-04 21:59:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/MultiSupportPattern.cpp,v $
                                                                        
// Written: fmk 11/00
// Revised:
//
// Purpose: This file contains the class definition for MultiSupportPattern.
// MultiSupportPattern is an abstract class.

#include <MultiSupportPattern.h>
#include <GroundMotion.h>

#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <stdlib.h>
#include <Channel.h>
#include <ErrorHandler.h>
#include <elementAPI.h>

void* OPS_MultiSupportPattern()
{
    int tag;

    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"insufficient number of args\n";
	return 0;
    }
    
    // get tag 
    int numData = 1;
    if(OPS_GetIntInput(&numData,&tag) < 0) {
	opserr << "WARNING failed to get multi support pattern tag\n";
	return 0;
    }

    return new MultiSupportPattern(tag);
}


MultiSupportPattern::MultiSupportPattern(int tag, int _classTag)
  :LoadPattern(tag, _classTag), 
   theMotions(0), theMotionTags(0,16), numMotions(0), dbMotions(0)
{

}


MultiSupportPattern::MultiSupportPattern(int tag)
  :LoadPattern(tag, PATTERN_TAG_MultiSupportPattern), 
   theMotions(0), theMotionTags(0,16), numMotions(0), dbMotions(0)
{

}


MultiSupportPattern::MultiSupportPattern()
  :LoadPattern(0, PATTERN_TAG_MultiSupportPattern), 
   theMotions(0), theMotionTags(0,16), numMotions(0), dbMotions(0)
{

}

MultiSupportPattern::~MultiSupportPattern()
{
  // invoke the destructor on all ground motions supplied
  for (int i=0; i<numMotions; i++)
      delete theMotions[i];
  
  if (theMotions != 0)
    //    free ((void *)theMotions);
    delete [] theMotions;
}



void 
MultiSupportPattern::applyLoad(double time)
{
  SP_Constraint *sp;
  SP_ConstraintIter &theIter = this->getSPs();
  while ((sp = theIter()) != 0) {
    sp->applyConstraint(time);
  }
}
    

int
MultiSupportPattern::addMotion(GroundMotion &theMotion, int tag)
{
  // ensure no motion with given tag already added
  if (theMotionTags.getLocation(tag) >= 0) {
    opserr << "MultiSupportPattern::addMotion - could not add new, motion with same tag exists\n";
    return -1;
  }    

  // make space for new
  GroundMotion **newMotions = new GroundMotion *[numMotions+1];
  //  GroundMotion **newMotions = (GroundMotion **)malloc(sizeof(GroundMotion *)*(numMotions+1));
  if (newMotions == 0) {
    opserr << "MultiSupportPattern::addMotion - could not add new, out of mem\n";
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
  theMotionTags[numMotions] = tag;
  numMotions++;

  return 0;
}



GroundMotion *
MultiSupportPattern::getMotion(int tag)
{
  int loc = theMotionTags.getLocation(tag);
  if (loc < 0)
    return 0;
  else
    return theMotions[loc];

}


bool
MultiSupportPattern::addNodalLoad(NodalLoad *)
{
  opserr << "MultiSupportPattern::addNodalLoad() - cannot add NodalLoad to EQ pattern\n";  
  return false;
}

bool
MultiSupportPattern::addElementalLoad(ElementalLoad *)
{
  opserr << "MultiSupportPattern::addElementalLoad() - cannot add ElementalLoad to EQ pattern\n";    
  return false;
}



int 
MultiSupportPattern::sendSelf(int commitTag, Channel &theChannel)
{
  // get my current database tag

  // NOTE - dbTag equals 0 if not sending to a database OR has not yet been sent
  int myDbTag = this->getDbTag();

  if (this->LoadPattern::sendSelf(commitTag, theChannel) < 0) {
    opserr << "MultiSupportPattern::sendSelf() - LoadPattern class failed in sendSelf()";
    return -1;
  }

  static ID myData(3);
  myData(0) = numMotions;
  if (dbMotions == 0)
    dbMotions = theChannel.getDbTag();
  myData(1) = dbMotions;
  
  if (theChannel.sendID(myDbTag, commitTag, myData) < 0) {
    opserr << "MultiSupportPattern::sendSelf - channel failed to send the initial ID\n";
    return -1;
  }    

  // create the ID and get the node iter
  if (numMotions != 0) {
    ID motionData(numMotions*3);
    GroundMotion *theMotion;
    for (int i=0; i<numMotions; i++) {
      theMotion = theMotions[i];
      motionData(i*3) = theMotion->getClassTag();
      int dbTag = theMotion->getDbTag();
	
      // if dbTag still 0 get one from Channel; 
      // if this tag != 0 set the dbTag in node
      if (dbTag == 0 && myDbTag != 0) {// go get a new tag and setDbTag in ele if this not 0 
	dbTag = theChannel.getDbTag();
	if (dbTag != 0)
	  theMotion->setDbTag(dbTag);
      }
      motionData(i*3+1) = dbTag;
      motionData(i*3+2) = theMotionTags(i);
    }    

    // now send the ID
    if (theChannel.sendID(dbMotions, commitTag, motionData) < 0) {
      opserr << "MultiSupportPattern::sendSelf - channel failed to send the NodalLoads ID\n";
      return -4;
    }

    
    for (int j=0; j<numMotions; j++) {
      theMotion = theMotions[j];
      if (theMotion->sendSelf(commitTag, theChannel) < 0) {
	opserr << "MultiSupportPattern::sendSelf - ground motion  failed in sendSelf\n";
	return -7;
      }
    }
  }

  // if we get here we are successful
  return 0;
}

int 
MultiSupportPattern::recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  // get my current database tag
  // NOTE - dbTag equals 0 if not sending to a database OR has not yet been sent
  int myDbTag = this->getDbTag();
  
  if (this->LoadPattern::recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "MultiSupportPattern::recvSelf() - LoadPattern class failed in sendSelf()";
    return -1;
  }
  
  // clear out the all the components in the current load pattern
  if (theMotions != 0) {
    for (int i=0; i<numMotions; i++)
      if (theMotions[i] != 0)
	delete theMotions[i];
    delete [] theMotions;
    numMotions = 0;
  }

  // 
  // now we rebuild the motions
  //

  static ID myData(3);
  if (theChannel.recvID(myDbTag, commitTag, myData) < 0) {
    opserr << "MultiSupportPattern::sendSelf - channel failed to send the initial ID\n";
    return -1;
  }    

  numMotions = myData(0);
  dbMotions = myData(1);

  if (numMotions != 0) {
    ID motionData(numMotions*3);

    // now send the ID
    if (theChannel.recvID(dbMotions, commitTag, motionData) < 0) {
      opserr << "MultiSupportPattern::sendSelf - channel failed to send the NodalLoads ID\n";
      return -4;
    }

    theMotions = new GroundMotion *[numMotions];
    for (int i=0; i<numMotions; i++)
      theMotions[i] = 0;

    if (theMotions == 0) {
      opserr << "MultiSupportPattern::recvSelf() - out of memory\n";
      return -1;
    }

    GroundMotion *theMotion;
    for (int i=0; i<numMotions; i++) {
      int motionClassTag = motionData(i*3);

      theMotion = theBroker.getNewGroundMotion(motionClassTag);
      if (theMotion == 0) {
	opserr << "MultiSupportPattern::sendSelf - failed to get a ground motion from object broker\n";
	return -1;
      }

      theMotion->setDbTag(motionData(i*3+1));
      int tag = motionData(i*3+2);

      if (theMotion->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "MultiSupportPattern::sendSelf - ground motion failed in recvSelf\n";
	return -7;
      }

      theMotionTags[i] = tag;
      theMotions[i] = theMotion;

      //      this->addMotion(*theMotion, tag);
    }    
  }

  return 0;
}

void 
MultiSupportPattern::Print(OPS_Stream &s, int flag)
{
  s << "MultiSupportPattern  tag: " << this->getTag() <<
    "   numMotions: " << numMotions << endln;

  SP_Constraint *sp;
  SP_ConstraintIter &theIter = this->getSPs();
  while ((sp = theIter()) != 0)
    sp->Print(s, flag);
  
}

LoadPattern *
MultiSupportPattern::getCopy(void)
{
  LoadPattern *theCopy = new MultiSupportPattern(this->getTag());
  for (int i=0; i<numMotions; i++) {
    theCopy->addMotion(*theMotions[i], theMotionTags[i]);
  }  
  return theCopy;
}
