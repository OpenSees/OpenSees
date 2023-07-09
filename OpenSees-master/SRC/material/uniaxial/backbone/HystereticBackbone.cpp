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
// $Date: 2008-11-09 06:05:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/HystereticBackbone.cpp,v $

// Written: MHS
// Created: Aug 2000
//
// Description: This file contains the interface for HystereticBackbone,
// which represents a backbone curve for hysteretic models.

#include <HystereticBackbone.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>

static MapOfTaggedObjects theHystereticBackboneObjects;

bool OPS_addHystereticBackbone(HystereticBackbone *newComponent) {
  return theHystereticBackboneObjects.addComponent(newComponent);
}

HystereticBackbone *OPS_getHystereticBackbone(int tag) {

  TaggedObject *theResult = theHystereticBackboneObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "HystereticBackbone *getHystereticBackbone(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  HystereticBackbone *theMat = (HystereticBackbone *)theResult;

  return theMat;
}

void OPS_clearAllHystereticBackbone(void) {
  theHystereticBackboneObjects.clearAll();
}


HystereticBackbone::HystereticBackbone (int tag, int classTag):
  TaggedObject(tag), MovableObject(classTag)
{
  
}

HystereticBackbone::~HystereticBackbone()
{
  
}

int 
HystereticBackbone::setVariable (char *argv)
{
  return -1;
}

int
HystereticBackbone::getVariable (int varID, double &theValue)
{
  return -1;
}

int
HystereticBackbone::setParameter(char **argv, int argc, Information &eleInformation)
{
  return -1;
}

int
HystereticBackbone::updateParameter(int responseID, Information &eleInformation)
{
  return -1;
}
