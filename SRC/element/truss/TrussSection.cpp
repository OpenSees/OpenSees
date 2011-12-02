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
// $Date: 2000-12-18 10:40:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/truss/TrussSection.cpp,v $
                                                                        
                                                                        
// File: ~/element/truss/TrussSection.C
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the implementation for the TrussSection class.
//
// What: "@(#) TrussSection.C, revA"

#include <TrussSection.h>
#include <Information.h>

#include <string.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <SectionForceDeformation.h>
#include <ElementResponse.h>

#include <G3Globals.h>

#include <Renderer.h>

#include <math.h>
#include <stdlib.h>

// initialise the class wide variables
Matrix TrussSection::trussM2(2,2);
Matrix TrussSection::trussM3(3,3);
Matrix TrussSection::trussM4(4,4);
Matrix TrussSection::trussM6(6,6);
Matrix TrussSection::trussM12(12,12);
Vector TrussSection::trussV2(2);
Vector TrussSection::trussV3(3);
Vector TrussSection::trussV4(4);
Vector TrussSection::trussV6(6);
Vector TrussSection::trussV12(12);

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the truss end nodes.
/*TrussSection::TrussSection(int tag, 
			   int dim,
			   int Nd1, int Nd2, 
			   SectionForceDeformation &theSect,
			   double rho)
:Element(tag,ELE_TAG_TrussSection),     
 connectedExternalNodes(2),
  dimension(dim), numDOF(0), theLoad(0),
 theMatrix(0), theVector(0), t(0),
 L(0.0), M(rho), end1Ptr(0), end2Ptr(0), theSection(0), sectDim(2)
{
    // get a copy of the section and check we obtained a valid copy
    theSection = theSect.getCopy();

    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;        
}
*/
TrussSection::TrussSection(int tag, 
			   int dim,
			   int Nd1, int Nd2, 
			   SectionForceDeformation &theSect,
			   double rho)
:Element(tag,ELE_TAG_TrussSection),     
 connectedExternalNodes(2),
  dimension(dim), numDOF(0), theLoad(0), 
 theMatrix(0), theVector(0), t(0),
 L(0.0), M(rho), end1Ptr(0), end2Ptr(0), theSection(0)
{
    // get a copy of the material and check we obtained a valid copy
    theSection = theSect.getCopy();
    if (theSection == 0) 
      g3ErrorHandler->fatal("FATAL TrussSection::TrussSection - failed to get a copy of material %d\n",
			    theSect.getTag());

	int order = theSection->getOrder();
	const ID &code = theSection->getType();
	
	int i;
	for (i = 0; i < order; i++)
		if (code(i) == SECTION_RESPONSE_P)
			break;

	if (i == order)
		g3ErrorHandler->warning("TrussSection::TrussSection - section does not provide axial response");

    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)
      g3ErrorHandler->fatal("FATAL TrussSection::TrussSection - failed to create an ID of correct size\n");
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;        
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
TrussSection::TrussSection()
:Element(0,ELE_TAG_TrussSection),     
 connectedExternalNodes(2),
 dimension(0), numDOF(0),
 theMatrix(0), theVector(0), t(0), 
 L(0.0), M(0.0), end1Ptr(0), end2Ptr(0), theSection(0)
{
    // ensure the connectedExternalNode ID is of correct size 
    if (connectedExternalNodes.Size() != 2)
      g3ErrorHandler->fatal("FATAL TrussSection::TrussSection - failed to create an ID of correct size\n");
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
TrussSection::~TrussSection()
{
    if (theSection != 0)
	delete theSection;
    if (t != 0)
	delete t;
}


int
TrussSection::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
TrussSection::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

int
TrussSection::getNumDOF(void) 
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
TrussSection::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	end1Ptr = 0;
	end2Ptr = 0;
	L = 0;
	return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);	

    // if nodes not in domain, warning message & set default numDOF as 2
    if ((end1Ptr == 0) || (end2Ptr == 0)){
      if (end1Ptr == 0)
        g3ErrorHandler->warning("TrussSection::setDomain() - Nd1: %d does not exist in ",Nd1);
      else
        g3ErrorHandler->warning("TrussSection::setDomain() - Nd2: %d does not exist in ",Nd2);

      g3ErrorHandler->warning("model for truss with id %d\n",this->getTag());

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	

      return;
    }

    // now determine the number of dof and the dimesnion    
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();	

    if (dofNd1 != dofNd2) {
      g3ErrorHandler->warning("WARNING TrussSection::setDomain(): nodes %d and %d %s %d\n",Nd1, Nd2,
			      "have differing dof at ends for truss",this->getTag());	

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
	
      return;
    }	

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    // now set the number of dof for element and set matrix and vector pointers
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
      g3ErrorHandler->warning("WARNING TrussSection::setDomain cannot handle %d dofs at nodes in %d d problem\n",
			      dimension,dofNd1);

      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
      return;
    }

    // create a transformation matrix for the element
    t = new Matrix(1,numDOF);
    if (t == 0 || (t->noCols() != numDOF)) {
      g3ErrorHandler->fatal("FATAL TrussSection::setDomain out of memory creating T matrix (1 x %d)\n",numDOF);
      return;
    }      
    
    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	

    if (dimension == 1) {
	Matrix &trans = *t;      
	trans(0,0) = -1;
	trans(0,1) = 1;

	double dx = end2Crd(0)-end1Crd(0);	
	L = sqrt(dx*dx);
	
	if (L == 0.0) {
	  g3ErrorHandler->warning("WARNING TrussSection::setDomain() - truss %d has zero length\n",
				  this->getTag());
	  return;
	}	

    } else if (dimension == 2) {
	double dx = end2Crd(0)-end1Crd(0);
	double dy = end2Crd(1)-end1Crd(1);	
    
	L = sqrt(dx*dx + dy*dy);
    
	if (L == 0.0) {
	  g3ErrorHandler->warning("WARNING TrussSection::setDomain() - truss %d has zero length\n",
				  this->getTag());
	  return;
	}
	
	double cs = dx/L;
	double sn = dy/L;

	Matrix &trans = *t;
	if (numDOF == 4) {
	    trans(0,0) = -cs;
	    trans(0,1) = -sn;
	    trans(0,2) = cs;
	    trans(0,3) = sn;
	} else { // it must be 6
	    trans(0,0) = -cs;
	    trans(0,1) = -sn;
	    trans(0,2) = 0.0;
	    trans(0,3) = cs;
	    trans(0,4) = sn;	
	    trans(0,5) = 0.0;
	}     

    } else {
	double dx = end2Crd(0)-end1Crd(0);
	double dy = end2Crd(1)-end1Crd(1);	
	double dz = end2Crd(2)-end1Crd(2);		
    
	L = sqrt(dx*dx + dy*dy + dz*dz);
    
	if (L == 0.0) {
	  g3ErrorHandler->warning("WARNING TrussSection::setDomain() - truss %d has zero length\n",
				  this->getTag());
	  return;
	}
	
	double cx = dx/L;
	double cy = dy/L;
	double cz = dz/L;	

	Matrix &trans = *t;	
	if (numDOF == 6) {
	    trans(0,0) = -cx;
	    trans(0,1) = -cy;
	    trans(0,2) = -cz;
	    trans(0,3) = cx;
	    trans(0,4) = cy;
	    trans(0,5) = cz;	    
	} else { // it must be 12
	    trans(0,0) = -cx;
	    trans(0,1) = -cy;
	    trans(0,2) = -cz;
	    trans(0,3) = 0;
	    trans(0,4) = 0;
	    trans(0,5) = 0;	    
	    trans(0,6) = cx;
	    trans(0,7) = cy;
	    trans(0,8) = cz;
	    trans(0,9) = 0;
	    trans(0,10) = 0;
	    trans(0,11) = 0;	    	    
	}     
    }


    // create the load vector
    theLoad = new Vector(numDOF);
    if (theLoad == 0) {
	g3ErrorHandler->fatal("TrussSection::setDomain - truss %d %s %d\n",
			      this->getTag(), 
			      "out of memory creating vector of size",
			      numDOF);	
      return;
    }          
    
    // determine the nodal mass for lumped mass approach
    //    double A = theSection->getA();
    M = M * L/2;
}   	 


int
TrussSection::commitState()
{
    return theSection->commitState();
}

int
TrussSection::revertToLastCommit()
{
    return theSection->revertToLastCommit();
}

int
TrussSection::revertToStart()
{
    return theSection->revertToStart();
}


const Matrix &
TrussSection::getTangentStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    // determine the current strain given trial displacements at nodes
    double strain = this->computeCurrentStrain();

	int order = theSection->getOrder();
	const ID &code = theSection->getType();
	
	Vector e (order);
	
	int i;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			e(i) = strain;
	}
	
	theSection->setTrialSectionDeformation(e);
    
	const Matrix &k = theSection->getSectionTangent();
	double AE = 0.0;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			AE += k(i,i);
	}

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;
    Matrix &trans = *t;

    stiff = trans^trans;

    stiff *= AE/L;  

    return *theMatrix;
}

const Matrix &
TrussSection::getSecantStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    // determine the current strain given trial displacements at nodes
    double strain = this->computeCurrentStrain();

	int order = theSection->getOrder();
	const ID &code = theSection->getType();
	
	Vector e (order);
	
	int i;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			e(i) = strain;
	}
	
	theSection->setTrialSectionDeformation(e);
    
	const Matrix &k = theSection->getSectionTangent();
	double AE = 0.0;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			AE += k(i,i);
	}
    
    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;
    Matrix &trans = *t;
    stiff = trans^trans;
    stiff *= AE/L;  
    
    return *theMatrix;
}
    
const Matrix &
TrussSection::getDamp(void)
{
  theMatrix->Zero();
  return *theMatrix; 
}


const Matrix &
TrussSection::getMass(void)
{   // zero the matrix
    theMatrix->Zero();    
  
    // check for quick return
    if (L == 0.0 || M == 0.0) { // - problem in setDomain() no further warnings
	return *theMatrix;
    }    

    Matrix &mass = *theMatrix;
    if (dimension == 1 && numDOF == 2) {
	mass(0,0) = M; 
	mass(1,1) = M;
    }
    else if (dimension == 2 && numDOF == 4) {
	mass(0,0) = M; 
	mass(1,1) = M;
	mass(2,2) = M; 
	mass(3,3) = M;	
    }
    else if (dimension == 2 && numDOF == 6) {
	mass(0,0) = M; 
	mass(1,1) = M;
	mass(2,2) = 0; 
	mass(3,3) = M;
	mass(4,4) = M; 
	mass(5,5) = 0;		
    }
    else if (dimension == 3 && numDOF == 6) {
	mass(0,0) = M; 
	mass(1,1) = M;
	mass(2,2) = M; 
	mass(3,3) = M;
	mass(4,4) = M; 
	mass(5,5) = M;		
    }
    else if (dimension == 3 && numDOF == 12) {
	mass(0,0) = M; 
	mass(1,1) = M;
	mass(2,2) = M; 
	mass(3,3) = 0;
	mass(4,4) = 0; 
	mass(5,5) = 0;		
	mass(6,6) = M; 
	mass(7,7) = M;
	mass(8,8) = M; 
	mass(9,9) = 0;
	mass(10,10) = 0; 
	mass(11,11) = 0;			
    }
    
    return *theMatrix; // so it will compile
}



void 
TrussSection::zeroLoad(void)
{
  theLoad->Zero();
}


int 
TrussSection::addLoad(const Vector &addP)
{
#ifdef _G3DEBUG    
    if (dimension != addP.Size()) {
	g3ErrorHandler->warning("TrussSection::addPtoUnbalance %s\n",
				"matrix and vector sizes are incompatable");
	return -1;
    }
#endif
    (*theLoad) += addP;
	return 0;
}

int 
TrussSection::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for a quick return
    if (L == 0.0 || M == 0.0) 
	return 0;

    // get R * accel from the nodes
    const Vector &Raccel1 = end1Ptr->getRV(accel);
    const Vector &Raccel2 = end2Ptr->getRV(accel);    

    int nodalDOF = numDOF/2;
    
#ifdef _G3DEBUG    
    if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size()) {
	g3ErrorHandler->warning("TrussSection::addInertiaLoadToUnbalance %s\n",
				"matrix and vector sizes are incompatable");
	return -1;
    }
#endif
    
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


const Vector &
TrussSection::getResistingForce()
{	
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theVector->Zero();
	return *theVector;
    }
    
    // determine the current strain
    double strain = this->computeCurrentStrain();

	int order = theSection->getOrder();
	const ID &code = theSection->getType();
	
	Vector e (order);
	
	int i;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			e(i) = strain;
	}
	
	theSection->setTrialSectionDeformation(e);
    
	const Vector &s = theSection->getStressResultant();
	double force = 0.0;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			force += s(i);
	}

    for (i=0; i<numDOF; i++)
	(*theVector)(i) = (*t)(0,i)*force;

    // add P
    (*theVector) -= *theLoad;

    return *theVector;
}


const Vector &
TrussSection::getResistingForceIncInertia()
{	
    this->getResistingForce();
    
    // now include the mass portion
    if (L != 0.0 && M != 0.0) {
	
	// remember we set M = M*A*L/2 in setDoamin()
	const Vector &accel1 = end1Ptr->getTrialAccel();
	const Vector &accel2 = end2Ptr->getTrialAccel();	
	
	int dof = dimension;
	int start = numDOF/2;
	for (int i=0; i<dof; i++) {
	    (*theVector)(i) = (*theVector)(i) - M*accel1(i);
	    (*theVector)(i+start) = (*theVector)(i+start) - M*accel2(i);
	}
    }    
    
    return *theVector;
}


int
TrussSection::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // truss packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(6);
  data(0) = this->getTag();
  data(1) = dimension;
  data(2) = numDOF;
  if (L != 0)
      data(3) = M * 2 / (L);
  else
      data(3) = M;
  data(4) = theSection->getClassTag();
  int matDbTag = theSection->getDbTag();

  // NOTE: we do have to ensure that the Section has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theSection->setDbTag(matDbTag);
  }
  data(5) = matDbTag;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING TrussSection::sendSelf() - %d failed to send Vector\n",this->getTag());
    return -1;
  }	      

  // truss then sends the tags of it's two end nodes

  res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING TrussSection::sendSelf() - %d failed to send ID\n",this->getTag());
    return -2;
  }

  // finally truss asks it's Section object to send itself

  res = theSection->sendSelf(commitTag, theChannel);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING TrussSection::sendSelf() - %d failed to send its Section\n",this->getTag());
    return -3;
  }

  return 0;
}

int
TrussSection::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{

  int res;
  int dataTag = this->getDbTag();

  // truss creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(6);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING TrussSection::recvSelf() - failed to receive Vector\n");
    return -1;
  }	      

  this->setTag((int)data(0));
  dimension = data(1);
  numDOF = data(2);
  M = data(3);

  // truss now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING TrussSection::recvSelf() - %d failed to receive ID\n", this->getTag());
    return -2;
  }

  // finally truss creates a new section object of the correct type,
  // sets its database tag and asks this new object to recveive itself.
  int sectClass = data(4);
  int sectDb = data(5);

  // Get new section if null
  if (theSection == 0)
	  theSection = theBroker.getNewSection(sectClass);

  // Check that section is of right type
  else if (theSection->getClassTag() != sectClass) {
	  delete theSection;
	  theSection = theBroker.getNewSection(sectClass);
  }
  
  // Check if either allocation failed
  if (theSection == 0) {
    g3ErrorHandler->warning("WARNING TrussSection::recvSelf() - %d failed to get a blank Section of type %d\n", 
			    this->getTag(), sectClass);
    return -3;
  }

  theSection->setDbTag(sectDb); // note: we set the dbTag before we receive the Section
  res = theSection->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING TrussSection::recvSelf() - %d failed to receive its Section\n", this->getTag());
    return -3;
  }

  return 0;
}


int
TrussSection::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // ensure setDomain() worked
    if (L == 0.0)
       return 0;

    // first determine the two end points of the truss based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	
    const Vector &end1Disp = end1Ptr->getDisp();
    const Vector &end2Disp = end2Ptr->getDisp();    

    if (displayMode == 1 || displayMode == 2) {
	static Vector v1(3);
	static Vector v2(3);
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

		int order = theSection->getOrder();
		const ID &code = theSection->getType();

		Vector e (order);
	
		int i;
		for (i = 0; i < order; i++) {
			if (code(i) == SECTION_RESPONSE_P)
				e(i) = strain;
		}
	
		theSection->setTrialSectionDeformation(e);
    
		const Vector &s = theSection->getStressResultant();
		for (i = 0; i < order; i++) {
			if (code(i) == SECTION_RESPONSE_P)
				force += s(i);
		}
	}
    
	if (displayMode == 2) // use the strain as the drawing measure
	    return theViewer.drawLine(v1, v2, strain, strain);	
	else { // otherwise use the axial force as measure
	    return theViewer.drawLine(v1,v2, force, force);
	}
    }
    return 0;
}


void
TrussSection::Print(ostream &s, int flag)
{
    // compute the strain and axial force in the member
    double strain, force;
    if (L == 0.0) {
	strain = 0;
	force = 0.0;
    } else {
		strain = this->computeCurrentStrain();	

		int order = theSection->getOrder();
		const ID &code = theSection->getType();

		Vector e (order);
	
		int i;
		for (i = 0; i < order; i++) {
			if (code(i) == SECTION_RESPONSE_P)
				e(i) = strain;
		}
	
		theSection->setTrialSectionDeformation(e);
    
		const Vector &s = theSection->getStressResultant();
		for (i = 0; i < order; i++) {
			if (code(i) == SECTION_RESPONSE_P)
				force += s(i);
		}
    }
    
    for (int i=0; i<numDOF; i++)
	(*theVector)(i) = (*t)(0,i)*force;
    
    if (flag == 0) { // print everything
	s << "Element: " << this->getTag(); 
	s << " type: TrussSection  iNode: " << connectedExternalNodes(0);
	s << " jNode: " << connectedExternalNodes(1);
	s << " Total Mass: " << M*2;
	
	s << " \n\t strain: " << strain;
	s << " axial load: " << force;
	if (theVector != 0) 
	    s << " \n\t unbalanced load: " << *theVector;	
	s << " \t Section: " << *theSection;
	s << endl;
    } else if (flag == 1) {
	s << this->getTag() << "  " << strain << "  ";
	s << force << endl;
    }
}

double
TrussSection::computeCurrentStrain(void) const
{
    // NOTE method will not be called if L == 0

    // determine the strain
    const Vector &disp1 = end1Ptr->getTrialDisp();
    const Vector &disp2 = end2Ptr->getTrialDisp();	

    double dLength = 0.0;
    for (int i=0; i<dimension; i++){
	dLength -= (disp2(i)-disp1(i))* (*t)(0,i);
    }

    // this method should never be called with L == 0
    return dLength/L;
}

Response*
TrussSection::setResponse(char **argv, int argc, Information &eleInformation)
{
    //
    // we compare argv[0] for known response types for the Truss
    //

    // axial force
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 || 
		strcmp(argv[0],"axialForce") == 0) {
		eleInformation.theType = DoubleType;
		return new ElementResponse(this, 1);
    } 

    else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
		strcmp(argv[0],"deformation") == 0) {
		eleInformation.theType = DoubleType;
		return new ElementResponse(this, 2);
    }     
    
    // tangent stiffness matrix
    else if (strcmp(argv[0],"stiff") ==0) {
		Matrix *newMatrix = new Matrix(*theMatrix);
		if (newMatrix == 0) {
			g3ErrorHandler->warning("WARNING TrussSection::setResponse() - %d out of memory creating matrix\n",
				this->getTag());
		return 0;
		}

		eleInformation.theMatrix = newMatrix;
		eleInformation.theType = MatrixType;
		return new ElementResponse(this, 3);
    } 

    // a section quantity    
    else if (strcmp(argv[0],"section") ==0)
		return theSection->setResponse(&argv[1], argc-1, eleInformation);
    
    // otherwise response quantity is unknown for the Truss class
    else
		return 0;
}

int 
TrussSection::getResponse(int responseID, Information &eleInformation)
{
 double strain, force;
 
  switch (responseID) {
    case 1:
      if (L == 0.0) {
	  strain = 0;
	  force = 0.0;
      } else {
	  strain = this->computeCurrentStrain();	
		int order = theSection->getOrder();
		const ID &code = theSection->getType();

		Vector e (order);
	
		int i;
		for (i = 0; i < order; i++) {
			if (code(i) == SECTION_RESPONSE_P)
				e(i) = strain;
		}
	
		theSection->setTrialSectionDeformation(e);
    
		const Vector &s = theSection->getStressResultant();
		for (i = 0; i < order; i++) {
			if (code(i) == SECTION_RESPONSE_P)
				force += s(i);
		}

      }      
      eleInformation.theDouble = force;    
      return 0;

    case 2:
      if (L == 0.0) {
	  strain = 0;
      } else {
	  strain = this->computeCurrentStrain();	
      }
      eleInformation.theDouble = strain*L;    
      return 0;
      
    case 3:
      if (L == 0.0) { // - problem in setDomain() no further warnings
	  theMatrix->Zero();
	  *(eleInformation.theMatrix) = *theMatrix;
	  return 0;
      } else {
    
	  // determine the current strain given trial displacements at nodes
	  strain = this->computeCurrentStrain();

	int order = theSection->getOrder();
	const ID &code = theSection->getType();
	
	Vector e (order);
	
	int i;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			e(i) = strain;
	}
	
	theSection->setTrialSectionDeformation(e);
    
	const Matrix &k = theSection->getSectionTangent();
	double AE = 0.0;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			AE += k(i,i);
	}
      
	  // come back later and redo this if too slow
	  Matrix &stiff = *theMatrix;
	  Matrix &trans = *t;
	  stiff = trans^trans;
	  stiff *= AE/L;  
	  *(eleInformation.theMatrix) = *theMatrix;      
	  return 0;
      }
      
    default:
      if (responseID >= 100)
	  return theSection->getResponse(responseID-100, eleInformation);
      else
	  return -1;
  }
}

int
TrussSection::setParameter (char **argv, int argc, Information &info)
{
    // a material parameter
    if (strcmp(argv[0],"section") == 0) {
		int ok = theSection->setParameter(&argv[1], argc-1, info);
		if (ok < 0)
			return -1;
		else
			return ok + 100;
    } 
    
    // otherwise parameter is unknown for the TrussSection class
    else
		return -1;

}
    
int
TrussSection::updateParameter (int parameterID, Information &info)
{
  switch (parameterID) {
    case -1:
      return -1;
      
    default:
      if (parameterID >= 100)
	  return theSection->updateParameter(parameterID-100, info);
      else
	  return -1;
  }
}
