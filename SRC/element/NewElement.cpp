/********************************************************************************
(C) Copyright 2001-2022, The Regents of the University of California    
All Rights Reserved.                                               

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
OF THE POSSIBILITY OF SUCH DAMAGE.

******************************************************************************* */
                                                                        
// Written: fmk 
// Created: 08/01

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
  theNodeTags(2), theMatrix(NUM_DOF, NUM_DOF), theVector(NUM_DOF)
{

}

NewElement::NewElement()
 :Element(0,ELE_TAG_NewElement),     
  theNodeTags(2), theMatrix(NUM_DOF, NUM_DOF), theVector(NUM_DOF)
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
    return theNodeTags;
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
    

const Vector &
NewElement::getResistingForce()
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


