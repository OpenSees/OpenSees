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
// $Date: 2006-09-05 22:32:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/limitCurve/LimitCurve.cpp,v $

// Written: kje 
// Created: 08/01
// Revision: A
//
// Description: This file contains the class implementation for 
// LimitCurve.
//
// What: "@(#) LimitCurve.C, revA"

#include <LimitCurve.h>
#include <string.h>
#include <Information.h>


#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>
#include <api/runtimeAPI.h>

static MapOfTaggedObjects theLimitCurveObjects;


bool OPS_addLimitCurve(LimitCurve *newComponent) {
  return theLimitCurveObjects.addComponent(newComponent);
}

LimitCurve *OPS_getLimitCurve(int tag) {

  TaggedObject *theResult = theLimitCurveObjects.getComponentPtr(tag);
  if (theResult == 0) {
    opserr << "LimitCurve *getLimitCurve(int tag) - none found with tag: " << tag << endln;
    return 0;
  }
  LimitCurve *theMat = (LimitCurve *)theResult;

  return theMat;
}

void
OPS_ADD_RUNTIME_VXV(OPS_clearAllLimitCurve)
{
  theLimitCurveObjects.clearAll();
}





LimitCurve::LimitCurve(int tag, int clasTag)
:TaggedObject(tag), MovableObject(clasTag)
{

}

LimitCurve::~LimitCurve()
{
  // does nothing
}
