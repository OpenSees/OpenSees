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

// $Revision: 1.5 $
// $Date: 2009-08-19 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CrdTransf.cpp,v $

// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the implementation for the CrdTransf class.
// CrdTransf provides the abstraction of a frame 
// coordinate transformation. It is an abstract base class and 
// thus no objects of  it's type can be instatiated. 

#include <CrdTransf.h>
#include <Vector.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>

static MapOfTaggedObjects theCrdTransfObjects;

bool 
OPS_addCrdTransf(CrdTransf *newComponent) {
  return theCrdTransfObjects.addComponent(newComponent);
}

CrdTransf *
OPS_getCrdTransf(int tag) {

  TaggedObject *theResult = theCrdTransfObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "CrdTransf *getCrdTransf(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  CrdTransf *theSeries = (CrdTransf *)theResult;

  return theSeries;
}

void 
OPS_clearAllCrdTransf(void) {
  theCrdTransfObjects.clearAll();
}


// constructor:
CrdTransf::CrdTransf(int tag, int classTag):TaggedObject(tag), MovableObject(classTag)
{
}

// destructor:
CrdTransf::~CrdTransf()
{
}

const Vector &
CrdTransf::getBasicDisplSensitivity(int gradNumber)
{
    opserr << "WARNING CrdTransf::getBasicDisplSensitivity() - this method "
        << " should not be called." << endln;
    
    static Vector dummy(1);
    return dummy;
}

const Vector &
CrdTransf::getGlobalResistingForceShapeSensitivity(const Vector &pb,
						   const Vector &p0,
						   int gradNumber)
{
    opserr << "ERROR CrdTransf::getGlobalResistingForceSensitivity() - has not been"
        << " implemented yet for the chosen transformation." << endln;
    
    static Vector dummy(1);
    return dummy;
}

const Vector &
CrdTransf::getGlobalResistingForceShapeSensitivity(const Vector &pb,
						   const Vector &p0)
{
    opserr << "ERROR CrdTransf::getGlobalResistingForceSensitivity() - has not been"
        << " implemented yet for the chosen transformation." << endln;
    
    static Vector dummy(1);
    return dummy;
}

const Vector &
CrdTransf::getBasicTrialDispShapeSensitivity(void)
{
    opserr << "ERROR CrdTransf::getBasicTrialDispShapeSensitivity() - has not been"
        << " implemented yet for the chosen transformation." << endln;
    
    static Vector dummy(1);
    return dummy;
}



// --Quan
const Vector &
CrdTransf::getBasicDisplSensitivity(int gradNumber, int)
{
    opserr << "WARNING CrdTransf::getBasicDisplSensitivity() - this method "
        << " should not be called." << endln;
    
    static Vector dummy(1);
    return dummy;
}
