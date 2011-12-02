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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:09 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/TclPlaneTruss/MyTruss.cpp,v $                                                                        
                                                                        
// File: ~/example/tcl/MyTruss.C
// 
// Written: fmk 
// Created: 02/99
// Revision: A
//
// Description: This file contains the implementation for the MyTruss class.
//
//
// What: "@(#) MyTruss.C, revA"


// we specify what header files we need
#include "MyTruss.h"

#include <G3Globals.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// initialise the class wide variables
Matrix MyTruss::trussK(4,4);
Matrix MyTruss::trussM(4,4);
Matrix MyTruss::trussD(4,4);
Vector MyTruss::trussR(4);

// typical constructor
MyTruss::MyTruss(int tag, 
                 int Nd1, int Nd2, 
                 UniaxialMaterial &theMat, 
                 double a,
		 double m)
:Element(tag,ELE_TAG_MyTruss),     
 externalNodes(2),
 t(0), L(0.0), A(a), M(m), end1Ptr(0), end2Ptr(0), load(4)
{	
    // get a copy of the material object for our own use
    theMaterial = theMat.getCopy();
    if (theMaterial == 0) 
      g3ErrorHandler->fatal("FATAL MyTruss::MyTruss() - out of memory, could not get a copy of the Material\n");
    
    // fill in the ID containing external node info with node id's    
    if (externalNodes.Size() != 2)
      g3ErrorHandler->fatal("FATAL MyTruss::MyTruss() - out of memory, could not create an ID of size 2\n");

    externalNodes(0) = Nd1;
    externalNodes(1) = Nd2;        
}

// constructor which should be invoked by an FE_ObjectBroker only
MyTruss::MyTruss()
:Element(0,ELE_TAG_MyTruss),     
 theMaterial(0),
 externalNodes(2),
 t(0), L(0.0), A(0.0), M(0.0), end1Ptr(0), end2Ptr(0)
{
  if (externalNodes.Size() != 2)
    g3ErrorHandler->fatal("FATAL MyTruss::MyTruss() - out of memory, could not create an ID of size 2\n");
}

//  destructor - provided to clean up any memory
MyTruss::~MyTruss()
{
    // clean up the memory associated with the element, this is
    // memory the MyTruss objects allocates in it's constructor, e.g.
    // the t matrix, and memory allocated by other objects that the 
    // MyTruss object is responsible for cleaning up, i.e. the MaterialObject.

    if (t != 0)
	delete t;
    
    if (theMaterial != 0)
	delete theMaterial;    
}

int
MyTruss::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
MyTruss::getExternalNodes(void) 
{
    return externalNodes;
}

int
MyTruss::getNumDOF(void) {
    return 4;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void
MyTruss::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	end1Ptr = 0;
	end2Ptr = 0;
	L = 0.0;
	return;
    }
    
    // first ensure nodes exist in Domain and set the node pointers
    int Nd1 = externalNodes(0);
    int Nd2 = externalNodes(1);
    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);	
    if (end1Ptr == 0) {
        g3ErrorHandler->warning("WARNING Truss::setDomain() - at truss %d node %d does not exist in domain\n",
				this->getTag(), Nd1);
	L = 0.0;
	return;  // don't go any further - otherwise segemntation fault
    }
    if (end2Ptr == 0) {        
        g3ErrorHandler->warning("WARNING Truss::setDomain() - at truss %d node %d does not exist in domain\n",
				this->getTag(), Nd1);
	L = 0.0;
	return;  // don't go any further - otherwise segemntation fault
    }	
    
    // call the DomainComponent class method THIS IS VERY IMPORTANT
    this->DomainComponent::setDomain(theDomain);

    // ensure connected nodes have correct number of dof's
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();	
    if ((dofNd1 != 2) || (dofNd2 != 2)) {
      g3ErrorHandler->warning("MyTruss::setDomain(): 2 dof required at nodes, %d and %d provided\n",
			      dofNd1, dofNd2);
    }	

    // create the transformation matrix and ensure we obtained it
    t = new Matrix(1,4);

    if ((t == 0) || ( t->noCols() != 4)) {
      g3ErrorHandler->warning("WARNING MyTruss::setDomain() - out of memory creating transformation matrix\n");
      L = 0.0;
      return;  // don't go any further - otherwise segemntation fault
    }
    
    // now determine the length & transformation matrix
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	

    double dx = end2Crd(0)-end1Crd(0);
    double dy = end2Crd(1)-end1Crd(1);	
    
    L = sqrt(dx*dx + dy*dy);
    
    if (L == 0.0) {
      g3ErrorHandler->warning("WARNING MyTruss::setDomain() - MyTruss %d has zero length\n", this->getTag());
      return;  // don't go any further - otherwise divide by 0 error
    }
	
    double cs = dx/L;
    double sn = dy/L;

    Matrix &trans = *t;
    trans(0,0) = cs;
    trans(0,1) = sn;    
    trans(0,2) = -cs;
    trans(0,3) = -sn;
}   	 


int
MyTruss::commitState()
{
    return theMaterial->commitState();
}

int
MyTruss::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}

int
MyTruss::revertToStart()
{
    return theMaterial->revertToStart();
}


const Matrix &
MyTruss::getTangentStiff(void)
{
    if (L == 0.0) { // length = zero - problem in setDomain() warning message already printed
	trussK.Zero();
	return trussK;
    }

    // determine the current strain given trial displacements at nodes
    double strain = this->computeCurrentStrain();

    // get the current E from the material for this strain
    theMaterial->setTrialStrain(strain);
    double E = theMaterial->getTangent();

    // form the tangent stiffness matrix
    Matrix &trans = *t;

    trussK = trans^trans;
    trussK *= A*E/L;  

    // return the matrix
    return trussK;
}

const Matrix &
MyTruss::getSecantStiff(void)
{
    if (L == 0.0) { // length = zero - problem in setDomain()
	trussK.Zero();
	return trussK;
    }
    
    // determine the current strain given trial displacements at nodes
    double strain = this->computeCurrentStrain();

    // get the current E from the material for this strain
    theMaterial->setTrialStrain(strain);
    double stress = theMaterial->getStress();    
    double E = stress/strain;

    // form the tangent stiffness matrix
    Matrix &trans = *t;

    trussK = trans^trans;
    trussK *= A*E/L;  

    // return the matrix
    return trussK;
}
    
const Matrix &
MyTruss::getDamp(void)
{
  return trussD;
}


const Matrix &
MyTruss::getMass(void)
{ 
    if (L == 0.0 || M == 0.0) { // length = zero - problem in setDomain()
	trussM.Zero();
	return trussM;
    }
    
    // determine mass matrix assuming lumped mass
    double nodeMass = M * A * L/2;
    for (int i=0; i<4; i++) trussM(i,i) = nodeMass;
    
    return trussM;
}

void 
MyTruss::zeroLoad(void)
{
    load.Zero();
}

int
MyTruss::addLoad(const Vector &addLoad)
{
    load += addLoad;
    return 0;
}

const Vector &
MyTruss::getResistingForce()
{	
    if (L == 0.0) { // if length == 0, problem in setDomain()
	trussR.Zero();
	return trussR;
    }
    
    // determine the current strain
    double strain = this->computeCurrentStrain();

    // get the current E from the material for this strain
    theMaterial->setTrialStrain(strain);

    // force = F * transformation 
    double force = A*theMaterial->getStress();
    for (int i=0; i<4; i++)
	trussR(i) = (*t)(0,i)*force;

    return trussR;
}



const Vector &
MyTruss::getResistingForceIncInertia()
{	
    // detrmine the resisting force sans mass
    this->getResistingForce();
    
    // now include the mass portion
    if (L != 0.0 && M != 0.0) {
	double nodeMass = M * A * L/2;
	
	const Vector &accel1 = end1Ptr->getTrialAccel();
	const Vector &accel2 = end2Ptr->getTrialAccel();	
	
	for (int i=0; i<2; i++) {
	    trussR(i) = trussR(i) - nodeMass*accel1(i);
	    trussR(i+2) = trussR(i+2) - nodeMass*accel2(i);	    
	}
    }

    return trussR;
}



int
MyTruss::sendSelf(int commitTag, Channel &theChannel)
{
    int res;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // MyTruss packs it's data into a Vector and sends this to theChannel
    // along with it's dbTag and the commitTag passed in the arguments

    Vector data(5);
    data(0) = this->getTag();
    data(1) = A;
    data(4) = M;
    data(2) = theMaterial->getClassTag();
    int matDbTag = theMaterial->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
	matDbTag = theChannel.getDbTag();
	if (matDbTag != 0)
	    theMaterial->setDbTag(matDbTag);
    }
    data(3) = matDbTag;

    res = theChannel.sendVector(dataTag, commitTag, data);
    if (res < 0) {
      g3ErrorHandler->warning("WARNING MyTruss::sendSelf() - failed to send Vector\n");
      return -1;
    }	      

    // MyTruss then sends the tags of it's two end nodes
    res = theChannel.sendID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      g3ErrorHandler->warning("WARNING MyTruss::sendSelf() - failed to send ID\n");
      return -2;
    }

    // finally MyTruss asks it's material object to send itself
    res = theMaterial->sendSelf(commitTag, theChannel);
    if (res < 0) {
      g3ErrorHandler->warning("WARNING MyTruss::sendSelf() - failed to send the Material\n");
      return -3;
    }

    return 0;
}

int
MyTruss::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res;
    int dataTag = this->getDbTag();

    // MyTruss creates a Vector, receives the Vector and then sets the 
    // internal data with the data in the Vector

    Vector data(5);
    res = theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
      g3ErrorHandler->warning("WARNING MyTruss::recvSelf() - failed to receive Vector\n");
      return -1;
    }	      

    this->setTag((int)data(0));
    A = data(1);
    M = data(4);
    
    // MyTruss now receives the tags of it's two external nodes
    res = theChannel.recvID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      g3ErrorHandler->warning("WARNING MyTruss::recvSelf() - failed to receive ID\n");
      return -2;
    }

    // we create a material object of the correct type,
    // sets its database tag and asks this new object to recveive itself.
    int matClass = data(2);
    int matDb = data(3);

    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      g3ErrorHandler->warning("WARNING MyTruss::recvSelf() - failed to create a Material\n");
      return -3;
    }

    // we set the dbTag before we receive the material  - this is important
    theMaterial->setDbTag(matDb); 
    res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      g3ErrorHandler->warning("WARNING MyTruss::recvSelf() - failed to receive the Material\n");
      return -3;
    }

    return 0;
}


int
MyTruss::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // check setDomain() was successfull
    if (L == 0.0)
       return 0;

    // first determine the two end points of the truss based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	
    const Vector &end1Disp = end1Ptr->getDisp();
    const Vector &end2Disp = end2Ptr->getDisp();    

    Vector v1(3);
    Vector v2(3);
    for (int i=0; i<2; i++) {
      v1(i) = end1Crd(i)+end1Disp(i)*fact;
      v2(i) = end2Crd(i)+end2Disp(i)*fact;    
    }

    // compute the strain and axial force in the member
    double strain, force;
    strain = this->computeCurrentStrain();
    theMaterial->setTrialStrain(strain);
    force = A*theMaterial->getStress();    
    
    if (displayMode == 2) // use the strain as the drawing measure
      return theViewer.drawLine(v1, v2, strain, strain);	
    else { // otherwise use the axial force as measure
      return theViewer.drawLine(v1,v2, force, force);
    }
}


void
MyTruss::Print(ostream &s, int flag)
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

    for (int i=0; i<4; i++)
      trussR(i) = (*t)(0,i)*force;

    if (flag == 0) { // print everything
      s << "Element: " << this->getTag(); 
      s << " type: MyTruss  iNode: " << externalNodes(0);
      s << " jNode: " << externalNodes(1);
      s << " Area: " << A;
      if (M != 0) s << " Mass (PerUnitVolume): " << M;	
	
      s << " \n\t strain: " << strain;
      s << " axial load: " <<  force;
      s << " \n\t unbalanced load: " << trussR;
      s << " \t Material: " << *theMaterial;
      s << endl;
    } else if (flag == 1) { // just print ele id, strain and force
      s << this->getTag() << "  " << strain << "  " << force << endl;
    }
}



int 
MyTruss::setResponse(char **argv, int argc, Information &eleInformation)
{
    //
    // we compare argv[0] for known response types for the Truss
    //

    // axial force
    if (strcmp(argv[0],"axialForce") ==0) {
	eleInformation.theType = DoubleType;
	return 1;
    } 

    // tangent stiffness matrix
    else if (strcmp(argv[0],"stiff") ==0) {
	Matrix *newMatrix = new Matrix(trussK);
	if (newMatrix == 0) {
	  g3ErrorHandler->warning("WARNING Truss::setResponse() - out of memory creating matrix\n");
	  return -1;
	}
	eleInformation.theMatrix = newMatrix;
	eleInformation.theType = MatrixType;
	return 2;
    } 

    // a material quantity    
    else if (strcmp(argv[0],"material") ==0) {
      int ok = theMaterial->setResponse(&argv[1], argc-1, eleInformation);
      if (ok < 0)
	return -1;
      else
	return ok + 100; // note - we add 100 to valid material response so that can idetify as meterial response
    } 
    
    // otherwise response quantity is unknown for the Truss class
    else
	return -1;
}



int 
MyTruss::getResponse(int responseID, Information &eleInformation)
{
 double strain;
 
  switch (responseID) {
    case -1:
      return -1;
      
    case 1:
      strain = this->computeCurrentStrain();
      theMaterial->setTrialStrain(strain);
      eleInformation.theDouble = A*theMaterial->getStress();    
      return 0;
      
    case 2:
      if (eleInformation.theMatrix != 0)
	  *(eleInformation.theMatrix) = this->getTangentStiff();
      return 0;      

    default:
      if (responseID >= 100) // we added 100 to valid material response
	  return theMaterial->getResponse(responseID-100, eleInformation);
      else
	  return -1;
  }
}

double
MyTruss::computeCurrentStrain(void) const
{
    // NOTE this method will never be called with L == 0.0

    // determine the strain
    const Vector &disp1 = end1Ptr->getTrialDisp();
    const Vector &disp2 = end2Ptr->getTrialDisp();	

    double dLength = 0.0;
    for (int i=0; i<2; i++){
      dLength += (disp2(i)-disp1(i))* (*t)(0,i);
    }
    
    double strain = dLength/L;

    return strain;
}


