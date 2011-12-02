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
                                                                        
// $Revision: 1.25 $
// $Date: 2008-12-18 22:45:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Element.cpp,v $
                                                                        
                                                                        
// Written: fmk 11/95
//
// Purpose: This file contains the class definition for Element.
// Element is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// The interface:
//

#include <stdlib.h>

#include "Element.h"
#include "ElementResponse.h"
//#include <Renderer.h>
#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Domain.h>

Matrix **Element::theMatrices; 
Vector **Element::theVectors1; 
Vector **Element::theVectors2; 
int  Element::numMatrices(0);

// Element(int tag, int noExtNodes);
// 	constructor that takes the element's unique tag and the number
//	of external nodes for the element.

Element::Element(int tag, int cTag) 
  :DomainComponent(tag, cTag), alphaM(0.0), 
  betaK(0.0), betaK0(0.0), betaKc(0.0), 
  Kc(0), index(-1), nodeIndex(-1)
{
    // does nothing
}


Element::~Element() 
{
  if (Kc != 0)
    delete Kc;
}

int
Element::commitState(void)
{
  if (Kc != 0)
    *Kc = this->getTangentStiff();
  
  return 0;
}

int
Element::update(void)
{
    return 0;
}

int
Element::revertToStart(void)
{
  return 0;
}


void 
Element::zeroLoad(void)
{

}

int 
Element::addInertiaLoadToUnbalance(const Vector &accel)
{
  return -1;
}


int
Element::setRayleighDampingFactors(double alpham, double betak, double betak0, double betakc)
{
  alphaM = alpham;
  betaK  = betak;
  betaK0 = betak0;
  betaKc = betakc;

  // check that memory has been allocated to store compute/return
  // damping matrix & residual force calculations
  if (index == -1) {
    int numDOF = this->getNumDOF();

    for (int i=0; i<numMatrices; i++) {
      Matrix *aMatrix = theMatrices[i];
      if (aMatrix->noRows() == numDOF) {
	index = i;
	i = numMatrices;
      }
    }
    if (index == -1) {
      Matrix **nextMatrices = new Matrix *[numMatrices+1];
      if (nextMatrices == 0) {
	opserr << "Element::getTheMatrix - out of memory\n";
      }
	  int j;
      for (j=0; j<numMatrices; j++)
	nextMatrices[j] = theMatrices[j];
      Matrix *theMatrix = new Matrix(numDOF, numDOF);
      if (theMatrix == 0) {
	opserr << "Element::getTheMatrix - out of memory\n";
	exit(-1);
      }
      nextMatrices[numMatrices] = theMatrix;

      Vector **nextVectors1 = new Vector *[numMatrices+1];
      Vector **nextVectors2 = new Vector *[numMatrices+1];
      if (nextVectors1 == 0 || nextVectors2 == 0) {
	opserr << "Element::getTheVector - out of memory\n";
	exit(-1);
      }

      for (j=0; j<numMatrices; j++) {
	nextVectors1[j] = theVectors1[j];
	nextVectors2[j] = theVectors2[j];
      }
	
      Vector *theVector1 = new Vector(numDOF);
      Vector *theVector2 = new Vector(numDOF);
      if (theVector1 == 0 || theVector2 == 0) {
	opserr << "Element::getTheVector - out of memory\n";
	exit(-1);
      }

      nextVectors1[numMatrices] = theVector1;
      nextVectors2[numMatrices] = theVector2;

      if (numMatrices != 0) {
	delete [] theMatrices;
	delete [] theVectors1;
	delete [] theVectors2;
      }
      index = numMatrices;
      numMatrices++;
      theMatrices = nextMatrices;
      theVectors1 = nextVectors1;
      theVectors2 = nextVectors2;
    }
  }

  // if need storage for Kc go get it
  if (betaKc != 0.0) {  
    if (Kc == 0) 
      Kc = new Matrix(this->getTangentStiff());
    if (Kc == 0) {
      opserr << "WARNING - ELEMENT::setRayleighDampingFactors - out of memory\n";
      betaKc = 0.0;
    }

    // if don't need storage for Kc & have allocated some for it, free the memory
  } else if (Kc != 0) { 
      delete Kc;
      Kc = 0;
  }

  return 0;
}

const Matrix &
Element::getDamp(void) 
{
  if (index  == -1) {
    this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  }

  // now compute the damping matrix
  Matrix *theMatrix = theMatrices[index]; 
  theMatrix->Zero();
  if (alphaM != 0.0)
    theMatrix->addMatrix(0.0, this->getMass(), alphaM);
  if (betaK != 0.0)
    theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);      
  if (betaK0 != 0.0)
    theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);      
  if (betaKc != 0.0)
    theMatrix->addMatrix(1.0, *Kc, betaKc);      

  // return the computed matrix
  return *theMatrix;
}



const Matrix &
Element::getMass(void)
{
  if (index  == -1) {
    this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  }

  // zero the matrix & return it
  Matrix *theMatrix = theMatrices[index]; 
  theMatrix->Zero();
  return *theMatrix;
}

const Vector &
Element::getResistingForceIncInertia(void) 
{
  if (index == -1) {
    this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  }

  Matrix *theMatrix = theMatrices[index]; 
  Vector *theVector = theVectors2[index];
  Vector *theVector2 = theVectors1[index];

  //
  // perform: R = P(U) - Pext(t);
  //

  (*theVector) = this->getResistingForce();

  //
  // perform: R = R - M * a
  //

  int loc = 0;
  Node **theNodes = this->getNodePtrs();
  int numNodes = this->getNumExternalNodes();

  int i;
  for (i=0; i<numNodes; i++) {
    const Vector &acc = theNodes[i]->getAccel();
    for (int i=0; i<acc.Size(); i++) {
      (*theVector2)(loc++) = acc(i);
    }
  }
  theVector->addMatrixVector(1.0, this->getMass(), *theVector2, +1.0);

  //
  // perform: R = R + (alphaM * M + betaK0 * K0 + betaK * K) * v
  //            = R + D * v
  //

  // determine the vel vector from ele nodes
  loc = 0;
  for (i=0; i<numNodes; i++) {
    const Vector &vel = theNodes[i]->getTrialVel();
    for (int i=0; i<vel.Size(); i++) {
      (*theVector2)(loc++) = vel[i];
    }
  }

  // now compute the damping matrix
  theMatrix->Zero();
  if (alphaM != 0.0)
    theMatrix->addMatrix(0.0, this->getMass(), alphaM);
  if (betaK != 0.0)
    theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);      
  if (betaK0 != 0.0)
    theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);      
  if (betaKc != 0.0)
    theMatrix->addMatrix(1.0, *Kc, betaKc);      

  // finally the D * v
  theVector->addMatrixVector(1.0, *theMatrix, *theVector2, 1.0);

  return *theVector;
}


const Vector &
Element::getRayleighDampingForces(void) 
{

  if (index == -1) {
    this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  }

  Matrix *theMatrix = theMatrices[index]; 
  Vector *theVector = theVectors2[index];
  Vector *theVector2 = theVectors1[index];

  //
  // perform: R = (alphaM * M + betaK0 * K0 + betaK * K) * v
  //            = D * v
  //

  // determine the vel vector from ele nodes
  Node **theNodes = this->getNodePtrs();
  int numNodes = this->getNumExternalNodes();
  int loc = 0;
  for (int i=0; i<numNodes; i++) {
    const Vector &vel = theNodes[i]->getTrialVel();
    for (int i=0; i<vel.Size(); i++) {
      (*theVector2)(loc++) = vel[i];
    }
  }

  // now compute the damping matrix
  theMatrix->Zero();
  if (alphaM != 0.0)
    theMatrix->addMatrix(0.0, this->getMass(), alphaM);
  if (betaK != 0.0)
    theMatrix->addMatrix(1.0, this->getTangentStiff(), betaK);      
  if (betaK0 != 0.0)
    theMatrix->addMatrix(1.0, this->getInitialStiff(), betaK0);      
  if (betaKc != 0.0)
    theMatrix->addMatrix(1.0, *Kc, betaKc);      

  // finally the D * v
  theVector->addMatrixVector(0.0, *theMatrix, *theVector2, 1.0);

  return *theVector;
}

int 
Element::addLoad(ElementalLoad *theLoad, double loadFactor) {
  return 0;
}

/*
int 
Element::addInertiaLoadToUnbalance(const Vector &accel)
{
  // some vectors to hold the load increment and RV
  int ndof = this->getNumDOF();
  Vector load(ndof);
  Vector RV(ndof);

  // 
  // for each node we will add it's R*accel contribution to RV
  //

  const ID &theNodes = this->getExternalNodes();
  int numNodes = theNodes.Size();
  int loc = 0;
  Domain *theDomain = this->getDomain();
  for (int i=0; i<numNodes; i++) {
    Node *theNode = theDomain->getNode(theNodes(i));
    if (theNode == 0)
      return -1;
    else {
      int numNodeDOF = theNode->getNumberDOF();
      const Vector &nodeRV = theNode->getRV(accel);
      for (int j=0; j<numNodeDOF; j++)
#ifdef _G3DEBUG
 	if (loc<ndof)
#endif
	  RV(loc++) = nodeRV(j);
    }
  }

  //
  // now we determine - M * R * accel
  //
  const Matrix &mass = this->getMass();
  load = mass * RV;
  load *= -1.0;

  return this->addLoad(load);
}
*/

bool
Element::isSubdomain(void)
{
    return false;
}

Response*
Element::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType",this->getClassType());
  output.attr("eleTag",this->getTag());
  int numNodes = this->getNumExternalNodes();
  const ID &nodes = this->getExternalNodes();
  static char nodeData[32];

  for (int i=0; i<numNodes; i++) {
    sprintf(nodeData,"node%d",i+1);
    output.attr(nodeData,nodes(i));
  }

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {
    const Vector &force = this->getResistingForce();
    int size = force.Size();
    for (int i=0; i<size; i++) {
      sprintf(nodeData,"P%d",i+1);
      output.tag("ResponseType",nodeData);
    }
    theResponse = new ElementResponse(this, 1, this->getResistingForce());
  }

  else if (strcmp(argv[0],"dampingForce") == 0 || strcmp(argv[0],"dampingForces") == 0) {
    const Vector &force = this->getResistingForce();
    int size = force.Size();
    for (int i=0; i<size; i++) {
      sprintf(nodeData,"P%d",i+1);
      output.tag("ResponseType",nodeData);
    }
    theResponse = new ElementResponse(this, 2, this->getResistingForce());
  }
  output.endTag();
  return theResponse;
}

int
Element::getResponse(int responseID, Information &eleInfo)
{
  int loc;

  switch (responseID) {
  case 1: // global forces
    return eleInfo.setVector(this->getResistingForce());
  case 2:
    return eleInfo.setVector(this->getRayleighDampingForces());
  default:
    return -1;
  }
}

int
Element::getResponseSensitivity(int responseID, int gradIndex,
				Information &eleInfo)
{
  return -1;
}

// AddingSensitivity:BEGIN //////////////////////////////////////////
const Vector &
Element::getResistingForceSensitivity(int gradIndex)
{
	static Vector dummy(1);
	return dummy;
}

const Matrix &
Element::getInitialStiffSensitivity(int gradIndex)
{
	static Matrix dummy(1,1);
	return dummy;
}

const Matrix &
Element::getMassSensitivity(int gradIndex)
{
	static Matrix dummy(1,1);
	return dummy;
}

int
Element::commitSensitivity(int gradIndex, int numGrads)
{
	return -1;
}

int
Element::addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool tag)
{
	return -1;
}



// AddingSensitivity:END ///////////////////////////////////////////


const Matrix &
Element::getDampSensitivity(int gradIndex) 
{
  if (index  == -1) {
    this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  }

  // now compute the damping matrix
  Matrix *theMatrix = theMatrices[index]; 
  theMatrix->Zero();
  if (alphaM != 0.0) {
    theMatrix->addMatrix(0.0, this->getMassSensitivity(gradIndex), alphaM);
  }
  if (betaK != 0.0) {
	theMatrix->addMatrix(1.0, this->getTangentStiff(), 0.0); // Don't use this and DDM      
	opserr << "Rayleigh damping with non-zero betaCurrentTangent is not compatible with DDM sensitivity analysis" << endln;
  }
  if (betaK0 != 0.0) {
    theMatrix->addMatrix(1.0, this->getInitialStiffSensitivity(gradIndex), betaK0);      
  }
  if (betaKc != 0.0) {
    theMatrix->addMatrix(1.0, *Kc, 0.0);      // Don't use this and DDM   
    opserr << "Rayleigh damping with non-zero betaCommittedTangent is not compatible with DDM sensitivity analysis" << endln;
  }

  // return the computed matrix
  return *theMatrix;
}


int 
Element::addResistingForceToNodalReaction(int flag) 
{
  int result = 0;
  int numNodes = this->getNumExternalNodes();
  Node **theNodes = this->getNodePtrs();

  //
  // first we create the nodes in static arrays as presume
  // we are going to do this many times & save new/delete
  //

  if (nodeIndex == -1) {
    int numLastDOF = -1;
    for (int i=0; i<numNodes; i++) {
      int numNodalDOF = theNodes[i]->getNumberDOF();

      if (numNodalDOF != 0 && numNodalDOF != numLastDOF) {
	// see if an existing vector will do
	int j;
	for (j=0; j<numMatrices; j++) {
	  if (theVectors1[j]->Size() == numNodalDOF)
	    nodeIndex = j;
	    j = numMatrices+2;
	}

	// if not we need to enlarge the bool of temp vectors, matrices
	if (j != (numMatrices+2)) {
	    Matrix **nextMatrices = new Matrix *[numMatrices+1];
	    if (nextMatrices == 0) {
	      opserr << "Element::addResistingForceToNodalReaction - out of memory\n";
	    }
	    int j;
	    for (j=0; j<numMatrices; j++)
	      nextMatrices[j] = theMatrices[j];
	    Matrix *theMatrix = new Matrix(numNodalDOF, numNodalDOF);
	    if (theMatrix == 0) {
	      opserr << "Element::addResistingForceToNodalReaction - out of memory\n";
	      exit(-1);
	    }
	    nextMatrices[numMatrices] = theMatrix;
	    
	    Vector **nextVectors1 = new Vector *[numMatrices+1];
	    Vector **nextVectors2 = new Vector *[numMatrices+1];
	    if (nextVectors1 == 0 || nextVectors2 == 0) {
	      opserr << "Element::addResistingForceToNodalReaction - out of memory\n";
	      exit(-1);
	    }
	    
	    for (j=0; j<numMatrices; j++) {
	      nextVectors1[j] = theVectors1[j];
	      nextVectors2[j] = theVectors2[j];
	    }
	    
	    Vector *theVector1 = new Vector(numNodalDOF);
	    Vector *theVector2 = new Vector(numNodalDOF);
	    if (theVector1 == 0 || theVector2 == 0) {
	      opserr << "Element::addResistingForceToNodalReaction - out of memory\n";
	      exit(-1);
	    }
	    
	    nextVectors1[numMatrices] = theVector1;
	    nextVectors2[numMatrices] = theVector2;
	    
	    if (numMatrices != 0) {
	      delete [] theMatrices;
	      delete [] theVectors1;
	      delete [] theVectors2;
	    }
	    nodeIndex = numMatrices;
	    numMatrices++;
	    theMatrices = nextMatrices;
	    theVectors1 = nextVectors1;
	    theVectors2 = nextVectors2;
	}
      }
      numLastDOF = numNodalDOF;
    }
  }

  //
  // now determine the resisting force
  //

  const Vector *theResistingForce;
  if (flag == 0)
    theResistingForce = &(this->getResistingForce());
  else if (flag == 1)
    theResistingForce = &(this->getResistingForceIncInertia());
  else if (flag == 2)
    theResistingForce = &(this->getRayleighDampingForces());

  if (nodeIndex == -1) {
    opserr << "LOGIC ERROR Element::addResistingForceToNodalReaction() -HUH!\n";
    return -1;
  }
  
  //
  // iterate over the elements nodes; determine nodes contribution & add it
  //

  int nodalDOFCount = 0;

  for (int i=0; i<numNodes; i++) {
    Node *theNode = theNodes[i];

    int numNodalDOF = theNode->getNumberDOF();
    Vector *theVector = theVectors1[nodeIndex];
    if (theVector->Size() != numNodalDOF) {
      for (int j=0; j<numMatrices; j++)
	if (theVectors1[j]->Size() == numNodalDOF) {
	  j = numMatrices;
	  theVector = theVectors1[j];
	}
    }
    for (int j=0; j<numNodalDOF; j++) {
      (*theVector)(j) = (*theResistingForce)(nodalDOFCount);
      nodalDOFCount++;
    }
    result +=theNode->addReactionForce(*theVector, 1.0);
  }

  return result;
}
