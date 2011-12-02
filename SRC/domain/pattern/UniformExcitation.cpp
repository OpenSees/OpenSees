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
// $Date: 2000-12-12 07:47:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/UniformExcitation.cpp,v $
                                                                        
                                                                        
// File: ~/domain/load/UniformExcitation.h
//
// Written: fmk 11/98
// Revised:
//
// Purpose: This file contains the class definition for UniformExcitation.
// UniformExcitation is an abstract class.

#include <UniformExcitation.h>
#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>

UniformExcitation::UniformExcitation(GroundMotion &_theMotion, 
				   int dof, int tag, double velZero)
:EarthquakePattern(tag, LOAD_TAG_UniformExcitation), 
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



int 
UniformExcitation::sendSelf(int commitTag, Channel &theChannel)
{
  cerr << "UniformExcitation::sendSelf() - not yet implemented\n";
  return 0;
}


int 
UniformExcitation::recvSelf(int commitTag, Channel &theChannel, 
			   FEM_ObjectBroker &theBroker)
{
  cerr << "UniformExcitation::recvSelf() - not yet implemented\n";
  return 0;
}


void 
UniformExcitation::Print(ostream &s, int flag)
{
  s << "UniformExcitation - Need to Print the GroundMotion\n";
}

LoadPattern *
UniformExcitation::getCopy(void)
{
  LoadPattern *theCopy = new UniformExcitation(*theMotion, theDof, this->getTag());
   return theCopy;
}
