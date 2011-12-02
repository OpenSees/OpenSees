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
                                                                        
// $Revision: 1.3 $
// $Date: 2000-12-18 10:40:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/truss/Truss.cpp,v $
                                                                        
                                                                        
// File: ~/element/truss/Truss.C
// 
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

#include <G3Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

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
	     double a, double rho)
 :Element(tag,ELE_TAG_Truss),     
  theMaterial(0), connectedExternalNodes(2),
  dimension(dim), numDOF(0), theLoad(0),
  theMatrix(0), theVector(0), t(0),
  L(0.0), A(a), M(rho), end1Ptr(0), end2Ptr(0)
{
    // get a copy of the material and check we obtained a valid copy
    theMaterial = theMat.getCopy();
    if (theMaterial == 0) 
      g3ErrorHandler->fatal("FATAL Truss::Truss - %d %s %d\n", tag,
			    "failed to get a copy of material with tag ",
			    theMat.getTag());

    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2)
      g3ErrorHandler->fatal("FATAL Truss::Truss - %d %s\n", tag,
			    "failed to create an ID of size 2");
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;        
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
Truss::Truss()
:Element(0,ELE_TAG_Truss),     
 theMaterial(0),connectedExternalNodes(2),
 dimension(0), numDOF(0),
 theMatrix(0), theVector(0), t(0), 
 L(0.0), A(0.0), M(0.0), end1Ptr(0), end2Ptr(0)
{
    // ensure the connectedExternalNode ID is of correct size 
    if (connectedExternalNodes.Size() != 2)
      g3ErrorHandler->fatal("FATAL Truss::Truss - %s\n",
			    "failed to create an ID of size 2");
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
    if (t != 0)
	delete t;
    if (theLoad != 0)
	delete theLoad;
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

    // if can't find both - send a warning message
    if ((end1Ptr == 0) || (end2Ptr == 0)) {
	g3ErrorHandler->warning("Truss::setDomain() - truss %d node %d %s\n",
				this->getTag(), Nd1,
				"does not exist in the model");

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	

      return;
    }

    // now determine the number of dof and the dimesnion    
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2) {
      g3ErrorHandler->warning("WARNING Truss::setDomain(): nodes %d and %d %s %d\n",Nd1, Nd2,
			      "have differing dof at ends for truss",this->getTag());

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
      g3ErrorHandler->warning("WARNING Truss::setDomain cannot handle %d dofs at nodes in %d d problem\n",
			      dimension,dofNd1);

      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
      return;
    }


    // create a transformation matrix for the element
    t = new Matrix(1,numDOF);
    theLoad = new Vector(numDOF);
    
    if (t == 0 || (t->noCols() != numDOF)) {
      g3ErrorHandler->fatal("Truss::setDomain - truss %d %s %d\n",
			    this->getTag(), 
			    "out of memory creating T matrix of size 1 x",
			    numDOF);	
      
      return;
    }      
    
    if (theLoad == 0) {
	g3ErrorHandler->fatal("Truss::setDomain - truss %d %s %d\n",
			      this->getTag(), 
			      "out of memory creating vector of size",
			      numDOF);	
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
	  g3ErrorHandler->warning("WARNING Truss::setDomain() - truss %d has zero length\n",
				  this->getTag());
	  return;
	}	

    } else if (dimension == 2) {
	double dx = end2Crd(0)-end1Crd(0);
	double dy = end2Crd(1)-end1Crd(1);	
    
	L = sqrt(dx*dx + dy*dy);
    
	if (L == 0.0) {
	  g3ErrorHandler->warning("WARNING Truss::setDomain() - truss %d has zero length\n",
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
	  g3ErrorHandler->warning("WARNING Truss::setDomain() - truss %d has zero length\n",
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
    
    // determine the nodal mass for lumped mass approach
    M = M * A * L/2;
}   	 


int
Truss::commitState()
{
    return theMaterial->commitState();
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


const Matrix &
Truss::getTangentStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    // determine the current strain given trial displacements at nodes
    double strain = this->computeCurrentStrain();

    // get the current E from the material for this strain
    theMaterial->setTrialStrain(strain);
    double E = theMaterial->getTangent();

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;
    Matrix &trans = *t;

    stiff = trans^trans;

    stiff *= A*E/L;  

    return *theMatrix;
}

const Matrix &
Truss::getSecantStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    // determine the current strain given trial displacements at nodes
    double strain = this->computeCurrentStrain();

    // get the current E from the material for this strain
    theMaterial->setTrialStrain(strain);
    double stress = theMaterial->getStress();    
    double E = stress/strain;

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;
    Matrix &trans = *t;
    stiff = trans^trans;
    stiff *= A*E/L;  
    
    return *theMatrix;
}
    
const Matrix &
Truss::getDamp(void)
{
    theMatrix->Zero();
    return *theMatrix; 
}


const Matrix &
Truss::getMass(void)
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
	mass(3,3) = M;
	mass(4,4) = M; 
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
	mass(6,6) = M; 
	mass(7,7) = M;
	mass(8,8) = M; 
    }
    
    return *theMatrix; // so it will compile
}



void 
Truss::zeroLoad(void)
{
    theLoad->Zero();
}

int 
Truss::addLoad(const Vector &addP)
{
#ifdef _G3DEBUG    
    if (dimension != addP.Size()) {
	g3ErrorHandler->warning("Truss::addPtoUnbalance %s\n",
				"matrix and vector sizes are incompatable");
	return -1;
    }
#endif
    (*theLoad) += addP;
	return 0;
}

int 
Truss::addInertiaLoadToUnbalance(const Vector &accel)
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
	g3ErrorHandler->warning("Truss::addInertiaLoadToUnbalance %s\n",
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
Truss::getResistingForce()
{	
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theVector->Zero();
	return *theVector;
    }
    
    // determine the current strain
    double strain = this->computeCurrentStrain();

    // get the current E from the material for this strain
    theMaterial->setTrialStrain(strain);

    // R = Ku - Pext
    // Ku = F * transformation
    double force = A*theMaterial->getStress();
    for (int i=0; i<numDOF; i++)
	(*theVector)(i) = (*t)(0,i)*force;

    // subtract external load:  Ku - P
    (*theVector) -= *theLoad;
    
    return *theVector;
}


const Vector &
Truss::getResistingForceIncInertia()
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
	    (*theVector)(i) = (*theVector)(i) + M*accel1(i);
	    (*theVector)(i+start) = (*theVector)(i+start) + M*accel2(i);
	}
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
  if (L != 0)
      data(6) = M * 2 / (L*A);
  else
      data(6) = M;
  
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
    g3ErrorHandler->warning("WARNING Truss::sendSelf() - %d failed to send Vector\n",this->getTag());
    return -1;
  }	      

  // truss then sends the tags of it's two end nodes

  res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING Truss::sendSelf() - %d failed to send Vector\n",this->getTag());
    return -2;
  }

  // finally truss asks it's material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING Truss::sendSelf() - %d failed to send its Material\n",this->getTag());
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
    g3ErrorHandler->warning("WARNING Truss::recvSelf() - failed to receive Vector\n");
    return -1;
  }	      

  this->setTag((int)data(0));
  dimension = data(1);
  numDOF = data(2);
  A = data(3);
  M = data(6);
  
  // truss now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING Truss::recvSelf() - %d failed to receive ID\n", this->getTag());
    return -2;
  }

  // finally truss creates a material object of the correct type,
  // sets its database tag and asks this new object to recveive itself.

  int matClass = data(4);
  int matDb = data(5);

  // check if we have a material object already & if we do if of right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

    // if old one .. delete it
    if (theMaterial != 0)
      delete theMaterial;

    // create a new material object
    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      g3ErrorHandler->warning("WARNING Truss::recvSelf() - %d failed to get a blank Material of type %d\n", 
			      this->getTag(), matClass);
      return -3;
    }
  }

  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    g3ErrorHandler->warning("WARNING Truss::recvSelf() - %d failed to receive its Material\n", this->getTag());
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
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	
    const Vector &end1Disp = end1Ptr->getDisp();
    const Vector &end2Disp = end2Ptr->getDisp();    

    if (displayMode == 1 || displayMode == 2) {
	Vector v1(3);
	Vector v2(3);
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
    }
    return 0;
}


void
Truss::Print(ostream &s, int flag)
{
    // compute the strain and axial force in the member
    double strain, force;
    if (L == 0.0) {
	strain = 0;
	force = 0.0;
    } else {
	strain = this->computeCurrentStrain();
	theMaterial->setTrialStrain(strain);
	force = A*theMaterial->getStress();    
    }
    
    for (int i=0; i<numDOF; i++)
	(*theVector)(i) = (*t)(0,i)*force;
    
    if (flag == 0) { // print everything
	s << "Element: " << this->getTag(); 
	s << " type: Truss  iNode: " << connectedExternalNodes(0);
	s << " jNode: " << connectedExternalNodes(1);
	s << " Area: " << A << " Total Mass: " << M*2;
	
	s << " \n\t strain: " << strain;
	s << " axial load: " << A*theMaterial->getStress();
	if (theVector != 0) 
	    s << " \n\t unbalanced load: " << *theVector;	
	s << " \t Material: " << *theMaterial;
	s << endl;
    } else if (flag == 1) {
	s << this->getTag() << "  " << strain << "  ";
	s << A*theMaterial->getStress() << endl;
    }
}

double
Truss::computeCurrentStrain(void) const
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
Truss::setResponse(char **argv, int argc, Information &eleInfo)
{
    //
    // we compare argv[0] for known response types for the Truss
    //

    // force (axialForce)
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"axialForce") == 0)
		return new ElementResponse(this, 1, 0.0);

    else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
		strcmp(argv[0],"deformation") == 0)
		return new ElementResponse(this, 2, 0.0);

    // tangent stiffness matrix
    else if (strcmp(argv[0],"stiff") == 0)
		return new ElementResponse(this, 3, *theMatrix);

    // a material quantity    
    else if (strcmp(argv[0],"material") == 0)
		return theMaterial->setResponse(&argv[1], argc-1, eleInfo);
    
	else
		return 0;
}

int 
Truss::getResponse(int responseID, Information &eleInfo)
{
 double strain;
 
  switch (responseID) {
    case 1:
      //strain = this->computeCurrentStrain();
      //theMaterial->setTrialStrain(strain);
      //eleInformation.theDouble = A*theMaterial->getStress();    
      //return 0;
	  return eleInfo.setDouble(A * theMaterial->getStress());
      
    case 2:
      //strain = this->computeCurrentStrain();
      //eleInformation.theDouble = strain*L;    
      //return 0;      
	  return eleInfo.setDouble(L * theMaterial->getStrain());
      
    case 3:
      //if (eleInformation.theMatrix != 0)
	  //*(eleInformation.theMatrix) = this->getTangentStiff();
      //return 0;      
	  return eleInfo.setMatrix(this->getTangentStiff());

    default:
	  return 0;
  }
}


int
Truss::setParameter (char **argv, int argc, Information &info)
{
    // a material parameter
    if (strcmp(argv[0],"material") == 0) {
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
      
    default:
      if (parameterID >= 100)
	  return theMaterial->updateParameter(parameterID-100, info);
      else
	  return -1;
  }
}
