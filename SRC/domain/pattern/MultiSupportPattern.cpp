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
// $Date: 2000-12-12 07:44:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/MultiSupportPattern.cpp,v $
                                                                        
// Written: fmk 11/00
// Revised:
//
// Purpose: This file contains the class definition for MultiSupportPattern.
// MultiSupportPattern is an abstract class.

#include <MultiSupportPattern.h>
#include <GroundMotion.h>

#include <Domain.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <stdlib.h>
#include <Channel.h>
#include <ErrorHandler.h>

MultiSupportPattern::MultiSupportPattern(int tag, int _classTag)
  :LoadPattern(tag, _classTag), 
  theMotions(0), theMotionTags(0,16), numMotions(0)
{

}


MultiSupportPattern::MultiSupportPattern(int tag)
  :LoadPattern(tag, PATTERN_TAG_MultiSupportPattern), 
  theMotions(0), theMotionTags(0,16), numMotions(0)
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
  while ((sp = theIter()) != 0)
    sp->applyConstraint(time);
}
    

int
MultiSupportPattern::addMotion(GroundMotion &theMotion, int tag)
{
  // ensure no motion with given tag already added
  if (theMotionTags.getLocation(tag) >= 0) {
    cerr << "MultiSupportPattern::addMotion - could not add new, motion wih same tag exists\n";
    return -1;
  }    

  // make space for new
  GroundMotion **newMotions = new GroundMotion *[numMotions+1];
  //  GroundMotion **newMotions = (GroundMotion **)malloc(sizeof(GroundMotion *)*(numMotions+1));
  if (newMotions == 0) {
    cerr << "MultiSupportPattern::addMotion - could not add new, out of mem\n";
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
  cerr << "MultiSupportPattern::addNodalLoad() - cannot add NodalLoad to EQ pattern\n";  
  return false;
}

bool
MultiSupportPattern::addElementalLoad(ElementalLoad *)
{
  cerr << "MultiSupportPattern::addElementalLoad() - cannot add ElementalLoad to EQ pattern\n";    
  return false;
}



int 
MultiSupportPattern::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
MultiSupportPattern::recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
MultiSupportPattern::Print(ostream &s, int flag)
{
  s << "MultiSupportPattern  tag: " << this->getTag() << endl;
  SP_Constraint *sp;
  SP_ConstraintIter &theIter = this->getSPs();
  while ((sp = theIter()) != 0)
    sp->Print(s, flag);
}

LoadPattern *
MultiSupportPattern::getCopy(void)
{
  LoadPattern *theCopy = new MultiSupportPattern(this->getTag());
  return theCopy;
}
