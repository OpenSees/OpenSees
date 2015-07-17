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
// $Date: 2007-02-02 01:30:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/NewElement.cpp,v $
                                                                        
// Written: fmk 
// Created: 08/01
//
// Description: This file contains the implementation for the NewElement class.
//
// What: "@(#) NewElement.cpp, revA"

#include "NewElement.h"
#include <Information.h>
#include <ElementResponse.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>

static int NUM_NODE =2;
static int NUM_DOF  =4;

// constructors:
NewElement::NewElement(int tag)
 :Element(tag,ELE_TAG_NewElement),     
  connectedExternalNodes(2), theMatrix(NUM_DOF, NUM_DOF), theVector(NUM_DOF)
{

}

NewElement::NewElement()
 :Element(0,ELE_TAG_NewElement),     
  connectedExternalNodes(2), theMatrix(NUM_DOF, NUM_DOF), theVector(NUM_DOF)
{

}

//  destructor:
NewElement::~NewElement()
{

}


int
NewElement::getNumExternalNodes(void) const
{
    return NUM_NODE;
}

const ID &
NewElement::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
NewElement::getNodePtrs(void)
{
  return 0;
}

int
NewElement::getNumDOF(void) 
{
    return NUM_DOF;
}


void
NewElement::setDomain(Domain *theDomain)
{
    // call the base class method
    this->DomainComponent::setDomain(theDomain);

}   	 


int
NewElement::commitState()
{
  int retVal = 0;

  // call the base class method
  retVal = this->Element::commitState();
  if (retVal < 0) {
    opserr << "NewElement::commitState() - failed in base class\n";
    return retVal;
  }

  return retVal;
}

int
NewElement::revertToLastCommit()
{
  return 0;
}

int
NewElement::revertToStart()
{
  return 0;
}

int
NewElement::update(void)
{
  return 0;
}


const Matrix &
NewElement::getTangentStiff(void)
{
  return theMatrix;
}

const Matrix &
NewElement::getInitialStiff(void)
{
  return theMatrix;
}
    
void 
NewElement::zeroLoad(void)
{
  return;
}

int 
NewElement::addLoad(const Vector &addP)
{
  return 0;
}

int 
NewElement::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
NewElement::getResistingForce()
{	
  return theVector;
}


const Vector &
NewElement::getResistingForceIncInertia()
{	
  return theVector;
}


int
NewElement::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
NewElement::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}


int
NewElement::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  return 0;
}


void
NewElement::Print(OPS_Stream &s, int flag)
{
  return;
}


Response*
NewElement::setResponse(const char **argv, int argc, OPS_Stream &S)
{
	return this->Element::setResponse(argv, argc, S);
}


int 
NewElement::getResponse(int responseID, Information &eleInfo)
{
	return this->Element::getResponse(responseID, eleInfo);
}


int
NewElement::setParameter(const char **argv, int argc, Parameter &param)
{
  return 0;
}
    

int
NewElement::updateParameter(int parameterID, Information &info)
{
  return -1;
}
