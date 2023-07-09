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

// $Revision: 1.0 $
// $Date: 2019-01-28 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/damping/Damping.cpp,v $

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 02/2020
// Revision: A
//
// Description: This file contains the implementation for the Damping class.
// Damping provides the abstraction of an elemental damping imposition.
// It is an abstract base class and thus no objects of it's type can be instatiated
// It has pure virtual functions which  must be implemented in it's derived classes.
// Reference:
// Yuan Tian, Yuli Huang, Zhe Qu, Yifan Fei, Xinzheng Lu,
// High-performance uniform damping model for response history analysis in OpenSees,
// Journal of Earthquake Engineering,
// 2022,
// https://doi.org/10.1080/13632469.2022.2124557
//

#include <Damping.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>

static MapOfTaggedObjects theDampingObjects;

bool OPS_addDamping(Damping *newComponent)
{
  return theDampingObjects.addComponent(newComponent);
}

bool OPS_removeDamping(int tag)
{
  TaggedObject* obj = theDampingObjects.removeComponent(tag);
  if (obj != 0)
  {
    delete obj;
    return true;
  }
  return false;
}

Damping *OPS_getDamping(int tag)
{
  TaggedObject *theResult = theDampingObjects.getComponentPtr(tag);
  if (theResult == 0)
  {
    opserr << "Damping *getDamping(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  Damping *theSeries = (Damping *)theResult;

  return theSeries;
}

void OPS_clearAllDamping(void)
{
  theDampingObjects.clearAll();
}

void OPS_printDamping(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON)
  {
    s << "\t\t\"Dampings\": [\n";        
    MapOfTaggedObjectsIter theObjects = theDampingObjects.getIter();
    theObjects.reset();
    TaggedObject *theObject;
    int count = 0;
    int numComponents = theDampingObjects.getNumComponents();    
    while ((theObject = theObjects()) != 0) {
      Damping *theDamping = (Damping *)theObject;
      theDamping->Print(s, flag);
      if (count < numComponents-1) s << ",\n";
      count++;      
    }
    s << "\n\t\t]";
  }
}


// constructor:
Damping::Damping(int tag, int classTag):TaggedObject(tag), MovableObject(classTag)
{
}

// destructor:
Damping::~Damping()
{
}

