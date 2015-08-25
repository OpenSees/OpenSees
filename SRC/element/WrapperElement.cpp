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
                                                                        
// $Revision: 1.4 $
// $Date: 2008-12-19 17:28:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/WrapperElement.cpp,v $
                                                                        
// Written: fmk 
#include "WrapperElement.h"
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

extern modelState theModelState;

Vector WrapperElement::Rvector;    // vector to hold the applied load P
Matrix WrapperElement::Kmatrix;
Matrix WrapperElement::Mmatrix; // AGGIUNTA
ID WrapperElement::connectedNodes;

WrapperElement::WrapperElement(const char *name, eleObject *ele)
  :Element(ele->tag,ELEMENT_TAGS_WrapperElement),
 funcName(0), theEle(ele),
 theNodes(0), u(0), R(0), K(0), M(0)
{

  funcName = new char[strlen(name)+1];
  if (funcName != 0)
    strcpy(funcName, name);
}

WrapperElement::WrapperElement()
:Element(0,ELEMENT_TAGS_WrapperElement), 
 funcName(0), theEle(0),
 theNodes(0), u(0), R(0), K(0), M(0)
{

}

WrapperElement::~WrapperElement()
{
  if (funcName != 0)
    delete [] funcName;

  if (u != 0)
    delete [] u;
  if (R != 0)
    delete [] R;
  if (K != 0)
    delete [] K;
  if (M != 0)
    delete [] M;
}

int
WrapperElement::getNumExternalNodes(void) const
{
  if (theEle != 0)
    return theEle->nNode;
  else
    return -1;
}

const ID &
WrapperElement::getExternalNodes(void)
{
  connectedNodes.setData(theEle->node, theEle->nNode);
  return connectedNodes;
}

Node **
WrapperElement::getNodePtrs(void)
{
  return theNodes;
}

int
WrapperElement::getNumDOF(void)
{
  return theEle->nDOF;
}

void
WrapperElement::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    if (theNodes != 0) {
      delete [] theNodes;
      theNodes = 0;
    }
    return;
  }    

  int numNodes = theEle->nNode;

  theNodes = new Node *[numNodes];
  for (int i=0; i<numNodes; i++) {
    Node *theNode = theDomain->getNode(theEle->node[i]);
    if (theNode == 0) {
      opserr << "WARNING WrapperElement::setDomain(Domain *theDomain) - node: ";
      opserr << theEle->node[i] << " does not exist in domain for ele " << *this;
      return;
    }
    // set up the pointer to the node
    theNodes[i] = theNode;
  }

  // call the DomainComponent class method THIS IS VERY IMPORTANT
  this->DomainComponent::setDomain(theDomain);
  
  u = new double [theEle->nDOF];
  R = new double [theEle->nDOF];
  K = new double [theEle->nDOF*theEle->nDOF];
  M = new double [theEle->nDOF*theEle->nDOF];
  if (u == 0 || R == 0 || K == 0 || M == 0) {
      opserr << "WARNING WrapperElement::setDomain(Domain *theDomain) - out of memory\n";
      return;
  }    
}

int
WrapperElement::commitState()
{
  // get the current load factor
  Domain *theDomain=this->getDomain();
  double time = theDomain->getCurrentTime();
  double dt = theDomain->getCurrentTime() - time;
  
  theModelState.time = time;
  theModelState.dt = dt;

  // invoke the element routine
  int isw = ISW_COMMIT;
  int error = 0;
  theEle->eleFunctPtr(theEle, &theModelState, K , R, &isw, &error);
  return error;
}

int
WrapperElement::revertToLastCommit()
{
  // get the current load factor
  Domain *theDomain=this->getDomain();
  double time = theDomain->getCurrentTime();
  double dt = theDomain->getCurrentTime() - time;
  
  theModelState.time = time;
  theModelState.dt = dt;

  // invoke the element routine
  int isw = ISW_REVERT;
  int error = 0;
  theEle->eleFunctPtr(theEle, &theModelState, K , R, &isw, &error);
  return error;
}

int
WrapperElement::revertToStart()
{
  // get the current load factor
  Domain *theDomain=this->getDomain();
  double time = theDomain->getCurrentTime();
  double dt = theDomain->getCurrentTime() - time;
  
  theModelState.time = time;
  theModelState.dt = dt;

  // invoke the element routine
  int isw = ISW_REVERT_TO_START;
  int error = 0;
  theEle->eleFunctPtr(theEle, &theModelState, K , R, &isw, &error);
  return error;
}

const Matrix &
WrapperElement::getTangentStiff(void)
{
    //// set the vector
    //Kmatrix.setData(K, theEle->nDOF, theEle->nDOF);

    //return  Kmatrix;

	//this->update();
    // set the vector
    Kmatrix.setData(K, theEle->nDOF, theEle->nDOF);

    return  Kmatrix;
}

const Matrix &
WrapperElement::getMass(void)
{
    // get the current load factor
    Domain *theDomain=this->getDomain();
    double time = theDomain->getCurrentTime();
    double dt = theDomain->getCurrentTime() - time;

    theModelState.time = time;
    theModelState.dt = dt;

    // invoke the element routine
    int isw = ISW_FORM_MASS;
    int error = 0;
    theEle->eleFunctPtr(theEle, &theModelState, M , R, &isw, &error);

	//this->update();
    // set the matrix
    Mmatrix.setData(M,theEle->nDOF, theEle->nDOF);

    return  Mmatrix;
}

const Vector &
WrapperElement::getResistingForce()
{		
    // set the vector
    Rvector.setData(R,theEle->nDOF);

    return  Rvector;
}

int
WrapperElement::sendSelf(int commitTag, Channel &theChannel)
{
  opserr << "WrapperElement::sendSelf() - not yet implemented\n";
  return -1;
}

int
WrapperElement::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  opserr << "WrapperElement::recvSelf() - not yet implemented\n";
  return -1;
}

int
WrapperElement::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    return 0;
}


void
WrapperElement::Print(OPS_Stream &s, int flag)
{
  s << "WrapperElement tag: " << this->getTag() << " Nodes: " << this->getExternalNodes();
}

int
WrapperElement::update()
{

    // get the current load factor
    Domain *theDomain=this->getDomain();
    double time = theDomain->getCurrentTime();
    double dt = theDomain->getCurrentTime() - time;

    theModelState.time = time;
    theModelState.dt = dt;
    // zero the matrix
    int error = 0;
    int isw = ISW_FORM_TANG_AND_RESID;
    theEle->eleFunctPtr(theEle, &theModelState, K , R, &isw, &error);
    
    Rvector.setData(R, theEle->nDOF);

    return error;
}

const Matrix &
WrapperElement::getInitialStiff(void)
{
  if (Ki == 0)
    Ki = new Matrix(this->getTangentStiff());

  if (Ki == 0) {
    opserr << "FATAL WrapperElement::getInitialStiff() -";
    opserr << "ran out of memory\n";
    exit(-1);
  }  
    
  return *Ki;
}




