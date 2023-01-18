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

#include "PlanarTruss.h"
#include <Information.h>
#include <ElementResponse.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <OPS_Globals.h>
#include <OPS_Stream.h>

#include <math.h>

static int NUM_NODE =2;
static int NUM_DOF  =4;

// constructors:
PlanarTruss::PlanarTruss(int tag, int iNode, int jNode,
			 UniaxialMaterial *theMat, double a) 
 :Element(tag,ELE_TAG_PlanarTruss),     
  connectedExternalNodes(2), theMatrix(NUM_DOF, NUM_DOF),
  theVector(NUM_DOF), A(a), L(0.), trans(1,4)
{
  connectedExternalNodes(0) = iNode;
  connectedExternalNodes(1) = jNode;
  theMaterial = theMat->getCopy();
  theNodes[0] = 0;
  theNodes[1] = 0;	
}

PlanarTruss::PlanarTruss()
 :Element(0,ELE_TAG_PlanarTruss),     
  connectedExternalNodes(2), theMatrix(NUM_DOF, NUM_DOF), theVector(NUM_DOF),
  A(0.), L(0.), trans(1,4), theMaterial(0)  
{
  theNodes[0] = 0;
  theNodes[1] = 0;
}

//  destructor:
PlanarTruss::~PlanarTruss()
{
  if (theMaterial != 0)
    delete theMaterial;
}


int
PlanarTruss::getNumExternalNodes(void) const
{
    return NUM_NODE;
}

const ID &
PlanarTruss::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
PlanarTruss::getNodePtrs(void)
{
  return theNodes;
}

int
PlanarTruss::getNumDOF(void) 
{
    return NUM_DOF;
}


void
PlanarTruss::setDomain(Domain *theDomain)
{
    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    //  determine L & transormation matrix
    int node1 = connectedExternalNodes(0);
    int node2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(node1);
    theNodes[1] = theDomain->getNode(node2);
    
    const Vector &crd1 = theNodes[0]->getCrds();
    const Vector &crd2 = theNodes[1]->getCrds();
    double dx = crd2(0)-crd1(0); 
    double dy = crd2(1)-crd1(1);
    L = sqrt(dx*dx + dy * dy); // requires #include <math.h>
    double cs = dx/L; double sn = dy/L;
    trans(0,0)=-cs; trans(0,1)=-sn;
    trans(0,2) = cs; trans(0,3)=sn;
    this->update();  
}   	 


int
PlanarTruss::commitState()
{
  return theMaterial->commitState();
}

int
PlanarTruss::revertToLastCommit()
{
  return theMaterial->revertToLastCommit();;
}

int
PlanarTruss::revertToStart()
{
  return theMaterial->revertToStart();
}

int
PlanarTruss::update(void)
{
  const Vector &disp1 = theNodes[0]->getTrialDisp();
  const Vector &disp2 = theNodes[1]->getTrialDisp();
  double dLength = 0.0;
  for (int i=0; i<2; i++) 
    dLength -= (disp2(i)-disp1(i)) * trans(0,i);
  double strain = dLength / L;
  
  return theMaterial->setTrialStrain(strain);
}


const Matrix &
PlanarTruss::getTangentStiff(void)
{
  double E = theMaterial->getTangent();
  theMatrix = trans ^ trans;
  theMatrix *= A*E/L;
  return theMatrix;
}

const Matrix &
PlanarTruss::getInitialStiff(void)
{
  double E = theMaterial->getInitialTangent();
  theMatrix = trans ^ trans;
  theMatrix *= A*E/L;
  return theMatrix;
}

const Vector &
PlanarTruss::getResistingForce()
{
  double force = A*theMaterial->getStress();
  for (int i=0; i<4; i++)
    theVector(i) = trans(0,i)*force;
  return theVector;
}

int
PlanarTruss::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
PlanarTruss::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}


void
PlanarTruss::Print(OPS_Stream &s, int flag)
{
  s << "Element: PlanarTruss tag: " << this->getTag() << " (nodes: " <<
    connectedExternalNodes[0] << " " << connectedExternalNodes(1) << ")\n";
  s << this->getResistingForce();
  return;
}


Response*
PlanarTruss::setResponse(const char **argv, int argc, OPS_Stream &S)
{
  return this->Element::setResponse(argv, argc, S);
}


int 
PlanarTruss::getResponse(int responseID, Information &eleInfo)
{
  return this->Element::getResponse(responseID, eleInfo);
}
