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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: August 2000
//
// Description: This file contains the interface for UnloadingRule,
// which models hysteretic unloading stiffness degradation.

#include <UnloadingRule.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>
#include <api/runtimeAPI.h>

static MapOfTaggedObjects theUnloadingRuleObjects;

bool OPS_addUnloadingRule(UnloadingRule *newComponent)
{
  return theUnloadingRuleObjects.addComponent(newComponent);
}

UnloadingRule *OPS_getUnloadingRule(int tag)
{
  TaggedObject *theResult = theUnloadingRuleObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "UnloadingRule *getUnloadingRule(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  UnloadingRule *theMat = (UnloadingRule *)theResult;

  return theMat;  
}

void
OPS_ADD_RUNTIME_VXV(OPS_clearAllUnloadingRule)
{
  theUnloadingRuleObjects.clearAll();
}

UnloadingRule::UnloadingRule(int tag, int classTag)
:MaterialState(tag,classTag)
{

}

UnloadingRule::~UnloadingRule()
{

}

UnloadingRule*
UnloadingRule::getCopy(UniaxialMaterial *u)
{
	return this->getCopy();
}
