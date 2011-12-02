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
// $Date: 2006-08-04 19:13:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/truss/Truss.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the implementation for the Truss class.
//
// What: "@(#) Truss.C, revA"

#include "Truss.h"
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

//#include <fstream>

// initialise the class wide variables
Matrix Truss::trussM2(2,2);
Matrix Truss::trussM4(4,4);
Matrix Truss::trussM6(6,6);
Matrix Truss::trussM12(12,12);
Vector Truss::trussV2(2);
Vector Truss::trussV4(4);
Vector Truss::trussV6(6);
Vector Truss::trussV12(12);

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the truss end nodes.
Truss::Truss(int tag, 
	     int dim,
	     int Nd1, int Nd2, 
	     UniaxialMaterial &theMat,
	     double a, double r)
 :Element(tag,ELE_TAG_Truss),     
  theMaterial(0), connectedExternalNodes(2),
  dimension(dim), numDOF(0), theLoad(0),
  theMatrix(0), theVector(0),
  L(0.0), A(a), rho(r)
{
    // get a copy of the material and check we obtained a valid copy
    theMaterial = theMat.getCopy();
    if (theMaterial == 0) {
      opserr << "FATAL Truss::Truss - " << tag <<
	"failed to get a copy of material with tag " << theMat.getTag() << endln;
      exit(-1);
    }
    
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL Truss::Truss - " <<  tag << "failed to create an ID of size 2\n";
      exit(-1);
    }

    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;        

    // set node pointers to NULL
    for (int i=0; i<2; i++)
      theNodes[i] = 0;

    cosX[0] = 0.0;
    cosX[1] = 0.0;
    cosX[2] = 0.0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	theLoadSens = 0;
// AddingSensitivity:END //////////////////////////////////////
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
Truss::Truss()
:Element(0,ELE_TAG_Truss),     
 theMaterial(0),connectedExternalNodes(2),
 dimension(0), numDOF(0), theLoad(0),
 theMatrix(0), theVector(0),
  L(0.0), A(0.0), rho(0.0)
{
    // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL Truss::Truss - failed to create an ID of size 2\n";
      exit(-1);
  }

  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  cosX[0] = 0.0;
  cosX[1] = 0.0;
  cosX[2] = 0.0;

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	theLoadSens = 0;
// AddingSensitivity:END //////////////////////////////////////
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
Truss::~Truss()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    if (theMaterial != 0)
	delete theMaterial;
    if (theLoad != 0)
	delete theLoad;
    if (theLoadSens != 0)
	delete theLoadSens;
}


int
Truss::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
Truss::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
Truss::getNodePtrs(void) 
{
  return theNodes;
}

int
Truss::getNumDOF(void) 
{
    return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the truss element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
Truss::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	theNodes[0] = 0;
	theNodes[1] = 0;
	L = 0;
	return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	
    
    // if can't find both - send a warning message
    if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
      if (theNodes[0] == 0)
	opserr <<"Truss::setDomain() - truss" << this->getTag() << " node " << Nd1 <<
	  "does not exist in the model\n";
      else
	opserr <<"Truss::setDomain() - truss" << this->getTag() << " node " << Nd2 <<
	  "does not exist in the model\n";

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	

      return;
    }

    // now determine the number of dof and the dimesnion    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2) {
      opserr <<"WARNING Truss::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
	"have differing dof at ends for truss " << this->getTag() << endln;

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
	
      return;
    }	

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    // now set the number of dof for element and set matrix and vector pointer
    if (dimension == 1 && dofNd1 == 1) {
	numDOF = 2;    
	theMatrix = &trussM2;
	theVector = &trussV2;
    }
    else if (dimension == 2 && dofNd1 == 2) {
	numDOF = 4;
	theMatrix = &trussM4;
	theVector = &trussV4;	
    }
    else if (dimension == 2 && dofNd1 == 3) {
	numDOF = 6;	
	theMatrix = &trussM6;
	theVector = &trussV6;		
    }
    else if (dimension == 3 && dofNd1 == 3) {
	numDOF = 6;	
	theMatrix = &trussM6;
	theVector = &trussV6;			
    }
    else if (dimension == 3 && dofNd1 == 6) {
	numDOF = 12;	    
	theMatrix = &trussM12;
	theVector = &trussV12;			
    }
    else {
      opserr <<"WARNING Truss::setDomain cannot handle " << dimension << " dofs at nodes in " << 
	dofNd1  << " problem\n";

      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
      return;
    }

    if (theLoad == 0)
      theLoad = new Vector(numDOF);
    else if (theLoad->Size() != numDOF) {
      delete theLoad;
      theLoad = new Vector(numDOF);
    }

    if (theLoad == 0) {
      opserr << "Truss::setDomain - truss " << this->getTag() << 
	"out of memory creating vector of size" << numDOF << endln;
      exit(-1);
      return;
    }          
    
    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    if (dimension == 1) {
      double dx = end2Crd(0)-end1Crd(0);	

      L = sqrt(dx*dx);
      
      if (L == 0.0) {
	opserr <<"WARNING Truss::setDomain() - truss " << this->getTag() << " has zero length\n";
	return;
      }

      cosX[0] = 1.0;
    }
    else if (dimension == 2) {
	double dx = end2Crd(0)-end1Crd(0);
	double dy = end2Crd(1)-end1Crd(1);	
    
	L = sqrt(dx*dx + dy*dy);
    
	if (L == 0.0) {
	  opserr <<"WARNING Truss::setDomain() - truss " << this->getTag() << " has zero length\n";
	  return;
	}
	
	cosX[0] = dx/L;
	cosX[1] = dy/L;
    }
    else {
	double dx = end2Crd(0)-end1Crd(0);
	double dy = end2Crd(1)-end1Crd(1);	
	double dz = end2Crd(2)-end1Crd(2);		
    
	L = sqrt(dx*dx + dy*dy + dz*dz);
    
	if (L == 0.0) {
	  opserr <<"WARNING Truss::setDomain() - truss " << this->getTag() << " has zero length\n";
	  return;
	}
	
	cosX[0] = dx/L;
	cosX[1] = dy/L;
	cosX[2] = dz/L;
    }
}   	 


int
Truss::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "Truss::commitState () - failed in base class";
  }    
  retVal = theMaterial->commitState();
  return retVal;
}

int
Truss::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}

int
Truss::revertToStart()
{
    return theMaterial->revertToStart();
}

int
Truss::update(void)
{
    // determine the current strain given trial displacements at nodes
    double strain = this->computeCurrentStrain();
    double rate = this->computeCurrentStrainRate();
    return theMaterial->setTrialStrain(strain, rate);
}


const Matrix &
Truss::getTangentStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    double E = theMaterial->getTangent();

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;

    int numDOF2 = numDOF/2;
    double temp;
    double EAoverL = E*A/L;
    for (int i = 0; i < dimension; i++) {
      for (int j = 0; j < dimension; j++) {
	temp = cosX[i]*cosX[j]*EAoverL;
	stiff(i,j) = temp;
	stiff(i+numDOF2,j) = -temp;
	stiff(i,j+numDOF2) = -temp;
	stiff(i+numDOF2,j+numDOF2) = temp;
      }
    }

    return stiff;
}


const Matrix &
Truss::getInitialStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    double E = theMaterial->getInitialTangent();

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;

    int numDOF2 = numDOF/2;
    double temp;
    double EAoverL = E*A/L;
    for (int i = 0; i < dimension; i++) {
      for (int j = 0; j < dimension; j++) {
	temp = cosX[i]*cosX[j]*EAoverL;
	stiff(i,j) = temp;
	stiff(i+numDOF2,j) = -temp;
	stiff(i,j+numDOF2) = -temp;
	stiff(i+numDOF2,j+numDOF2) = temp;
      }
    }

    return *theMatrix;
}

const Matrix &
Truss::getDamp(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    double eta = theMaterial->getDampTangent();

    // come back later and redo this if too slow
    Matrix &damp = *theMatrix;

    int numDOF2 = numDOF/2;
    double temp;
    double etaAoverL = eta*A/L;
    for (int i = 0; i < dimension; i++) {
      for (int j = 0; j < dimension; j++) {
	temp = cosX[i]*cosX[j]*etaAoverL;
	damp(i,j) = temp;
	damp(i+numDOF2,j) = -temp;
	damp(i,j+numDOF2) = -temp;
	damp(i+numDOF2,j+numDOF2) = temp;
      }
    }

    return damp;
}


const Matrix &
Truss::getMass(void)
{   
  // zero the matrix
  Matrix &mass = *theMatrix;
  mass.Zero();    
  
  // check for quick return
  if (L == 0.0 || rho == 0.0) { // - problem in setDomain() no further warnings
    return mass;
  }    
  
  double M = 0.5*rho*L;
  int numDOF2 = numDOF/2;
  for (int i = 0; i < dimension; i++) {
    mass(i,i) = M;
    mass(i+numDOF2,i+numDOF2) = M;
  }

  return mass;
}

void 
Truss::zeroLoad(void)
{
  theLoad->Zero();
}

int 
Truss::addLoad(ElementalLoad *theLoad, double loadFactor)

{  
  opserr <<"Truss::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  
  return -1;
}

int 
Truss::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for a quick return
    if (L == 0.0 || rho == 0.0) 
	return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);    

  int nodalDOF = numDOF/2;
    
#ifdef _G3DEBUG    
  if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size()) {
    opserr <<"Truss::addInertiaLoadToUnbalance " <<
      "matrix and vector sizes are incompatable\n";
    return -1;
  }
#endif
    
  double M = 0.5*rho*L;
    // want to add ( - fact * M R * accel ) to unbalance
    for (int i=0; i<dimension; i++) {
	double val1 = Raccel1(i);
	double val2 = Raccel2(i);	
	
	// perform - fact * M*(R * accel) // remember M a diagonal matrix
	val1 *= -M;
	val2 *= -M;
	
	(*theLoad)(i) += val1;
	(*theLoad)(i+nodalDOF) += val2;
    }	

    return 0;
}


int 
Truss::addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool somethingRandomInMotions)
{

  if (theLoadSens == 0) {
    theLoadSens = new Vector(numDOF);
  }
  else {
    theLoadSens->Zero();
  }
  
  
  if (somethingRandomInMotions) {
    
    
    // check for a quick return
    if (L == 0.0 || rho == 0.0) 
      return 0;
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);    
    
    int nodalDOF = numDOF/2;
    
#ifdef _G3DEBUG    
    if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size()) {
      opserr << "Truss::addInertiaLoadToUnbalance " <<
	"matrix and vector sizes are incompatable\n";
      return -1;
    }
#endif
    
	double M  = 0.5*rho*L;
    // want to add ( - fact * M R * accel ) to unbalance
    for (int i=0; i<dimension; i++) {
      double val1 = Raccel1(i);
      double val2 = Raccel2(i);	
      
      // perform - fact * M*(R * accel) // remember M a diagonal matrix
      val1 *= M;
      val2 *= M;
      
      (*theLoadSens)(i) = val1;
      (*theLoadSens)(i+nodalDOF) = val2;
    }	
  }
  else {
    
    // check for a quick return
    if (L == 0.0 || rho == 0.0) 
      return 0;
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);    
    
    int nodalDOF = numDOF/2;
    
#ifdef _G3DEBUG    
    if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size()) {
      opserr << "Truss::addInertiaLoadToUnbalance " <<
	"matrix and vector sizes are incompatable\n";
      return -1;
    }
#endif
    
    double massDerivative = 0.0;
    if (parameterID == 2) {
      massDerivative = 0.5*L;
    }
      
    // want to add ( - fact * M R * accel ) to unbalance
    for (int i=0; i<dimension; i++) {
      double val1 = Raccel1(i);
      double val2 = Raccel2(i);	
      
      // perform - fact * M*(R * accel) // remember M a diagonal matrix
      
      val1 *= massDerivative;
      val2 *= massDerivative;
      
      (*theLoadSens)(i) = val1;
      (*theLoadSens)(i+nodalDOF) = val2;
    }	
  }
  return 0;
}

const Vector &
Truss::getResistingForce()
{	
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theVector->Zero();
	return *theVector;
    }
    
    // R = Ku - Pext
    // Ku = F * transformation
    double force = A*theMaterial->getStress();
    int numDOF2 = numDOF/2;
    double temp;
    for (int i = 0; i < dimension; i++) {
      temp = cosX[i]*force;
      (*theVector)(i) = -temp;
      (*theVector)(i+numDOF2) = temp;
    }

    // subtract external load:  Ku - P
    (*theVector) -= *theLoad;

    return *theVector;
}


const Vector &
Truss::getResistingForceIncInertia()
{	
    this->getResistingForce();
    
    // now include the mass portion
    if (L != 0.0 && rho != 0.0) {
	
	const Vector &accel1 = theNodes[0]->getTrialAccel();
	const Vector &accel2 = theNodes[1]->getTrialAccel();	
	
	int numDOF2 = numDOF/2;
	double M = 0.5*rho*L;
	for (int i = 0; i < dimension; i++) {
	    (*theVector)(i) += M*accel1(i);
	    (*theVector)(i+numDOF2) += M*accel2(i);
	}

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	  (*theVector) += this->getRayleighDampingForces();
    }  else {

      // add the damping forces if rayleigh damping
      if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	(*theVector) += this->getRayleighDampingForces();
    }
    
    return *theVector;
}

int
Truss::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // truss packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(7);
  data(0) = this->getTag();
  data(1) = dimension;
  data(2) = numDOF;
  data(3) = A;
  data(6) = rho;
  
  data(4) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();

  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(5) = matDbTag;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING Truss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }	      

  // truss then sends the tags of it's two end nodes

  res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr <<"WARNING Truss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  // finally truss asks it's material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr <<"WARNING Truss::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }

  return 0;
}

int
Truss::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();

  // truss creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(7);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr <<"WARNING Truss::recvSelf() - failed to receive Vector\n";
    return -1;
  }	      

  this->setTag((int)data(0));
  dimension = (int)data(1);
  numDOF = (int)data(2);
  A = data(3);
  rho = data(6);
  
  // truss now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr <<"WARNING Truss::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally truss creates a material object of the correct type,
  // sets its database tag and asks this new object to recveive itself.

  int matClass = (int)data(4);
  int matDb = (int)data(5);

  // check if we have a material object already & if we do if of right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

    // if old one .. delete it
    if (theMaterial != 0)
      delete theMaterial;

    // create a new material object
    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      opserr <<"WARNING Truss::recvSelf() - " << this->getTag() 
	<< " failed to get a blank Material of type " << matClass << endln;
      return -3;
    }
  }

  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr <<"WARNING Truss::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
    return -3;    
  }

  return 0;
}


int
Truss::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // ensure setDomain() worked
    if (L == 0.0)
       return 0;

    // first determine the two end points of the truss based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    static Vector v1(3);
    static Vector v2(3);

    if (displayMode == 1 || displayMode == 2) {
      const Vector &end1Disp = theNodes[0]->getDisp();
      const Vector &end2Disp = theNodes[1]->getDisp();    

	for (int i=0; i<dimension; i++) {
	    v1(i) = end1Crd(i)+end1Disp(i)*fact;
	    v2(i) = end2Crd(i)+end2Disp(i)*fact;    
	}
	
	// compute the strain and axial force in the member
	double strain, force;
	if (L == 0.0) {
	    strain = 0.0;
	    force = 0.0;
	} else {
	    strain = this->computeCurrentStrain();
	    theMaterial->setTrialStrain(strain);
	    force = A*theMaterial->getStress();    
	}
    
	if (displayMode == 2) // use the strain as the drawing measure
	  return theViewer.drawLine(v1, v2, strain, strain);	
	else { // otherwise use the axial force as measure
	  return theViewer.drawLine(v1,v2, force, force);
	}
    } else if (displayMode < 0) {
      int mode = displayMode  *  -1;
      const Matrix &eigen1 = theNodes[0]->getEigenvectors();
      const Matrix &eigen2 = theNodes[1]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < dimension; i++) {
	  v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	  v2(i) = end2Crd(i) + eigen2(i,mode-1)*fact;    
	}    
      } else {
	for (int i = 0; i < dimension; i++) {
	  v1(i) = end1Crd(i);
	  v2(i) = end2Crd(i);
	}    
      }
    }
    return 0;
}


void
Truss::Print(OPS_Stream &s, int flag)
{
    // compute the strain and axial force in the member
    double strain, force;
    strain = theMaterial->getStrain();
    force = A * theMaterial->getStress();
    
    if (flag == 0) { // print everything
	s << "Element: " << this->getTag(); 
	s << " type: Truss  iNode: " << connectedExternalNodes(0);
	s << " jNode: " << connectedExternalNodes(1);
	s << " Area: " << A << " Mass/Length: " << rho;
	
	s << " \n\t strain: " << strain;
	s << " axial load: " <<  force;
	if (L != 0.0) {
	  int numDOF2 = numDOF/2;
	  double temp;
	  for (int i = 0; i < dimension; i++) {
	    temp = cosX[i]*force;
	    (*theVector)(i) = -temp;
	    (*theVector)(i+numDOF2) = temp;
	  }
	  s << " \n\t unbalanced load: " << *theVector;	
	}

	s << " \t Material: " << *theMaterial;
	s << endln;
    } else if (flag == 1) {
	s << this->getTag() << "  " << strain << "  ";
	s << force << endln;
    }
}

double
Truss::computeCurrentStrain(void) const
{
    // NOTE method will not be called if L == 0

    // determine the strain
    const Vector &disp1 = theNodes[0]->getTrialDisp();
    const Vector &disp2 = theNodes[1]->getTrialDisp();	

    double dLength = 0.0;
    for (int i = 0; i < dimension; i++)
      dLength += (disp2(i)-disp1(i))*cosX[i];
  
    // this method should never be called with L == 0
    return dLength/L;
}

double
Truss::computeCurrentStrainRate(void) const
{
    // NOTE method will not be called if L == 0

    // determine the strain
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();	

    double dLength = 0.0;
    for (int i = 0; i < dimension; i++){
	dLength += (vel2(i)-vel1(i))*cosX[i];
    }

    // this method should never be called with L == 0
    return dLength/L;
}

Response*
Truss::setResponse(const char **argv, int argc, Information &eleInfo, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","Truss");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

  //
  // we compare argv[0] for known response types for the Truss
  //

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"axialForce") == 0) {
    opserr << "HELLO\n";
    output.tag("ResponseType", "N");
    theResponse =  new ElementResponse(this, 1, 0.0);

  } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	     strcmp(argv[0],"deformation") == 0) {

    output.tag("ResponseType", "eps");
    theResponse = new ElementResponse(this, 2, 0.0);

  // a material quantity    
  } else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"-material") == 0) {

    theResponse =  theMaterial->setResponse(&argv[1], argc-1, eleInfo, output);
  }

  output.endTag();
  return theResponse;
}

int 
Truss::getResponse(int responseID, Information &eleInfo)
{
  switch (responseID) {
    case 1:
      return eleInfo.setDouble(A * theMaterial->getStress());
      
    case 2:
      return eleInfo.setDouble(L * theMaterial->getStrain());
      
    case 3:
      return eleInfo.setMatrix(this->getTangentStiff());

    default:
      return 0;
  }
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
Truss::setParameter (const char **argv, int argc, Information &info)
{
    if (argc < 1)
        return -1;

    // Cross sectional area of the truss itself
    if (strcmp(argv[0],"A") == 0) {
        info.theType = DoubleType;
        return 1;
    }

    // Mass densitity (per unit volume) of the truss itself
    if (strcmp(argv[0],"rho") == 0) {
        info.theType = DoubleType;
        return 2;
    }

    // a material parameter
    if (strcmp(argv[0],"-material") == 0 || strcmp(argv[0],"material") == 0) {
      int ok = theMaterial->setParameter(&argv[1], argc-1, info);
      if (ok < 0)
	return -1;
      else
	return ok + 100;
    } 
    
    // otherwise parameter is unknown for the Truss class
    else
      return -1;
}

int
Truss::updateParameter (int parameterID, Information &info)
{
  switch (parameterID) {
    case -1:
      return -1;
      
    case 1:
        this->A = info.theDouble;
        return 0;

    case 2:
        this->rho = info.theDouble;
        return 0;

    default:
      if (parameterID >= 100)
	  return theMaterial->updateParameter(parameterID-100, info);
      else
	  return -1;
  }
}
int
Truss::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	// The identifier needs to be passed "downwards" also when it's zero
	if (passedParameterID == 0 || passedParameterID == 1 || passedParameterID == 2) {
		theMaterial->activateParameter(0);
	}

	// If the identifier is non-zero and the parameter belongs to the material
	else if ( passedParameterID > 100) {
		theMaterial->activateParameter(passedParameterID-100);
	}

	return 0;
}


const Matrix &
Truss::getKiSensitivity(int gradNumber)
{
  Matrix &stiff = *theMatrix;
  stiff.Zero();
    
  if (parameterID == 0) {
  }
  else if (parameterID == 1) {
    // If cross sectional area is random
    double E = theMaterial->getInitialTangent();
    
    int numDOF2 = numDOF/2;
    double temp;
    double EoverL = E/L;
    for (int i = 0; i < dimension; i++) {
      for (int j = 0; j < dimension; j++) {
	temp = cosX[i]*cosX[j]*EoverL;
	stiff(i,j) = temp;
	stiff(i+numDOF2,j) = -temp;
	stiff(i,j+numDOF2) = -temp;
	stiff(i+numDOF2,j+numDOF2) = temp;
      }
    }
  }
  else if (parameterID == 2) {
    // Nothing here when 'rho' is random
  }
  else {
    double Esens = theMaterial->getInitialTangentSensitivity(gradNumber);

    int numDOF2 = numDOF/2;
    double temp;
    double EAoverL = Esens*A/L;
    for (int i = 0; i < dimension; i++) {
      for (int j = 0; j < dimension; j++) {
	temp = cosX[i]*cosX[j]*EAoverL;
	stiff(i,j) = temp;
	stiff(i+numDOF2,j) = -temp;
	stiff(i,j+numDOF2) = -temp;
	stiff(i+numDOF2,j+numDOF2) = temp;
      }
    }
  }
  
  return stiff;
}

const Matrix &
Truss::getMassSensitivity(int gradNumber)
{
  Matrix &mass = *theMatrix;
  mass.Zero();

  if (parameterID == 2) {
    double massDerivative = 0.5*L;

    int numDOF2 = numDOF/2;
    for (int i = 0; i < dimension; i++) {
      mass(i,i) = massDerivative;
      mass(i+numDOF2,i+numDOF2) = massDerivative;
    }
  }

  return mass;
}

const Vector &
Truss::getResistingForceSensitivity(int gradNumber)
{
	theVector->Zero();

	// Initial declarations
	int i;
	double stressSensitivity, temp1, temp2;

	// Make sure the material is up to date
	double strain = this->computeCurrentStrain();
	double rate = this->computeCurrentStrainRate();
	theMaterial->setTrialStrain(strain,rate);

	// Contribution from material
	stressSensitivity = theMaterial->getStressSensitivity(gradNumber,true);

	// Check if a nodal coordinate is random
	double dcosXdh[3];
	dcosXdh[0] = 0.0;
	dcosXdh[1] = 0.0;
	dcosXdh[2] = 0.0;
	
	int nodeParameterID0 = theNodes[0]->getCrdsSensitivity();
	int nodeParameterID1 = theNodes[1]->getCrdsSensitivity();
	if (nodeParameterID0 != 0 || nodeParameterID1 != 0) {
	
	  double dx = L*cosX[0];
	  double dy = L*cosX[1];
	  //double dz = L*cosX[2];

		// Compute derivative of transformation matrix (assume 4 dofs)
		if (nodeParameterID0 == 1) { // here x1 is random
			temp1 = (-L+dx*dx/L)/(L*L);
			temp2 = dx*dy/(L*L*L);
			//dtdh(0) = -temp1;
			//dtdh(1) = -temp2;
			//dtdh(2) = temp1;
			//dtdh(3) = temp2;
			dcosXdh[0] = temp1;
			dcosXdh[1] = temp2;
			dcosXdh[2] = 0.0;
		}
		if (nodeParameterID0 == 2) { // here y1 is random
			temp1 = (-L+dy*dy/L)/(L*L);
			temp2 = dx*dy/(L*L*L);
			//dtdh(0) = -temp2;
			//dtdh(1) = -temp1;
			//dtdh(2) = temp2;
			//dtdh(3) = temp1;
			dcosXdh[0] = temp2;
			dcosXdh[1] = temp1;
			dcosXdh[2] = 0.0;
		}
		if (nodeParameterID1 == 1) { // here x2 is random
			temp1 = (L-dx*dx/L)/(L*L);
			temp2 = -dx*dy/(L*L*L);
			//dtdh(0) = -temp1;
			//dtdh(1) = -temp2;
			//dtdh(2) = temp1;
			//dtdh(3) = temp2;
			dcosXdh[0] = temp1;
			dcosXdh[1] = temp2;
			dcosXdh[2] = 0.0;
		}
		if (nodeParameterID1 == 2) { // here y2 is random
			temp1 = (L-dy*dy/L)/(L*L);
			temp2 = -dx*dy/(L*L*L);
			//dtdh(0) = -temp2;
			//dtdh(1) = -temp1;
			//dtdh(2) = temp2;
			//dtdh(3) = temp1;
			dcosXdh[0] = temp2;
			dcosXdh[1] = temp1;
			dcosXdh[2] = 0.0;
		}

		const Vector &disp1 = theNodes[0]->getTrialDisp();
		const Vector &disp2 = theNodes[1]->getTrialDisp();
		double dLengthDerivative = 0.0;
		for (i = 0; i < dimension; i++) {
			dLengthDerivative += (disp2(i)-disp1(i))*dcosXdh[i];
		}

		double materialTangent = theMaterial->getTangent();
		double strainSensitivity;

		if (nodeParameterID0 == 1) {		// here x1 is random
			strainSensitivity = (dLengthDerivative*L+strain*dx)/(L*L);
		}
		if (nodeParameterID0 == 2) {	// here y1 is random
			strainSensitivity = (dLengthDerivative*L+strain*dy)/(L*L);
		}
		if (nodeParameterID1 == 1) {		// here x2 is random
			strainSensitivity = (dLengthDerivative*L-strain*dx)/(L*L);
		}
		if (nodeParameterID1 == 2) {	// here y2 is random
			strainSensitivity = (dLengthDerivative*L-strain*dy)/(L*L);
		}
		stressSensitivity += materialTangent * strainSensitivity;
	}


	// Compute sensitivity depending on 'parameter'
	double stress = theMaterial->getStress();
	int numDOF2 = numDOF/2;
	double temp;
	if (parameterID == 1) {			// Cross-sectional area
	  for (i = 0; i < dimension; i++) {
	    temp = (stress + A*stressSensitivity)*cosX[i];
	    (*theVector)(i) = -temp;
	    (*theVector)(i+numDOF2) = temp;
	  }
	}
	else {		// Density, material parameter or nodal coordinate
	  for (i = 0; i < dimension; i++) {
	    temp = A*(stressSensitivity*cosX[i] + stress*dcosXdh[i]);
	    (*theVector)(i) = -temp;
	    (*theVector)(i+numDOF2) = temp;
	  }
	}

	// subtract external load sensitivity
	if (theLoadSens == 0) {
		theLoadSens = new Vector(numDOF);
	}
	(*theVector) -= *theLoadSens;

	return *theVector;
}

int
Truss::commitSensitivity(int gradNumber, int numGrads)
{
	// Initial declarations
	int i; 
	double strainSensitivity, temp1, temp2;

	// Displacement difference between the two ends
	double strain = this->computeCurrentStrain();
	double dLength = strain*L;

	// Displacement sensitivity difference between the two ends
	double sens1;
	double sens2;
	double dSensitivity = 0.0;
	for (i=0; i<dimension; i++){
	  sens1 = theNodes[0]->getDispSensitivity(i+1, gradNumber);
	  sens2 = theNodes[1]->getDispSensitivity(i+1, gradNumber);
	  dSensitivity += (sens2-sens1)*cosX[i];
	}

	strainSensitivity = dSensitivity/L;

	// Check if a nodal coordinate is random
	int nodeParameterID0 = theNodes[0]->getCrdsSensitivity();
	int nodeParameterID1 = theNodes[1]->getCrdsSensitivity();
	if (nodeParameterID0 != 0 || nodeParameterID1 != 0) {

	  double dx = L*cosX[0];
	  double dy = L*cosX[1];
	  //double dz = L*cosX[2];

		// Compute derivative of transformation matrix (assume 4 dofs)
		double dcosXdh[3];

		if (nodeParameterID0 == 1) { // here x1 is random
			temp1 = (-L+dx*dx/L)/(L*L);
			temp2 = dx*dy/(L*L*L);
			//dtdh(0) = -temp1;
			//dtdh(1) = -temp2;
			//dtdh(2) = temp1;
			//dtdh(3) = temp2;
			dcosXdh[0] = temp1;
			dcosXdh[1] = temp2;
			dcosXdh[2] = 0.0;
		}
		if (nodeParameterID0 == 2) { // here y1 is random
			temp1 = (-L+dy*dy/L)/(L*L);
			temp2 = dx*dy/(L*L*L);
			//dtdh(0) = -temp2;
			//dtdh(1) = -temp1;
			//dtdh(2) = temp2;
			//dtdh(3) = temp1;
			dcosXdh[0] = temp2;
			dcosXdh[1] = temp1;
			dcosXdh[2] = 0.0;
		}

		if (nodeParameterID1 == 1) { // here x2 is random
			temp1 = (L-dx*dx/L)/(L*L);
			temp2 = -dx*dy/(L*L*L);
			//dtdh(0) = -temp1;
			//dtdh(1) = -temp2;
			//dtdh(2) = temp1;
			//dtdh(3) = temp2;
			dcosXdh[0] = temp1;
			dcosXdh[1] = temp2;
			dcosXdh[2] = 0.0;
		}
		if (nodeParameterID1 == 2) { // here y2 is random
			temp1 = (L-dy*dy/L)/(L*L);
			temp2 = -dx*dy/(L*L*L);
			//dtdh(0) = -temp2;
			//dtdh(1) = -temp1;
			//dtdh(2) = temp2;
			//dtdh(3) = temp1;
			dcosXdh[0] = temp2;
			dcosXdh[1] = temp1;
			dcosXdh[2] = 0.0;
		}

		const Vector &disp1 = theNodes[0]->getTrialDisp();
		const Vector &disp2 = theNodes[1]->getTrialDisp();
		double dLengthDerivative = 0.0;
		for (i = 0; i < dimension; i++){
			dLengthDerivative += (disp2(i)-disp1(i))*dcosXdh[i];
		}

		strainSensitivity += dLengthDerivative/L;

		if (nodeParameterID0 == 1) {		// here x1 is random
			strainSensitivity += dLength/(L*L*L)*dx;
		}
		if (nodeParameterID0 == 2) {	// here y1 is random
			strainSensitivity += dLength/(L*L*L)*dy;
		}
		if (nodeParameterID1 == 1) {		// here x2 is random
			strainSensitivity -= dLength/(L*L*L)*dx;
		}
		if (nodeParameterID1 == 2) {	// here y2 is random
			strainSensitivity -= dLength/(L*L*L)*dy;
		}
	}
	
	// Pass it down to the material
	theMaterial->commitSensitivity(strainSensitivity, gradNumber, numGrads);

	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////
