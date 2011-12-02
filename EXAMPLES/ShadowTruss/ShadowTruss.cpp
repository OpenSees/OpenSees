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
// $Date: 2003-10-15 00:38:07 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/ShadowTruss/ShadowTruss.cpp,v $
                                                                        
// Written: fmk 
// Created: 08/03
//
// Description: This file contains the implementation for the ShadowTruss class.
//
// What: "@(#) ShadowTruss.C, revA"


// we specify what header files we need
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <OPS_Globals.h>

#include "ShadowTruss.h"
#include "ShadowActorTruss.h"

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Response.h>


// initialise the class wide variables
Matrix ShadowTruss::trussK(4,4);
Matrix ShadowTruss::trussM(4,4);
Vector ShadowTruss::trussR(4);
Node * ShadowTruss::theNodes[2];

// typical constructor
ShadowTruss::ShadowTruss(int tag, 
			 int Nd1, int Nd2, 
			 UniaxialMaterial &theMat,
			 double a, 
			 double m,
			 Channel &theChannel, 
			 FEM_ObjectBroker &theObjectBroker)
:Element(tag, ELE_TAG_ShadowTruss),     
 Shadow(theChannel, theObjectBroker),
 msgData(2),
 externalNodes(2),
 trans(1, 4), L(0.0), A(a), M(m), end1Ptr(0), end2Ptr(0), theLoad(0)
{	
  // fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 2)
    opserr << "FATAL ShadowTruss::ShadowTruss() - out of memory, could not create an ID of size 2\n";

  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;        

  // send area & material to remote
  msgData(0) = ShadowActorTruss_setMaterial;
  msgData(1) = theMat.getClassTag();
  this->sendID(msgData);
  static Vector data(1);
  data(0) = A;
  this->sendVector(data);
  this->sendObject(theMat);

}

//  destructor - provided to clean up any memory
ShadowTruss::~ShadowTruss()
{
  msgData(0) = ShadowActorTruss_DIE;
  this->sendID(msgData);

  if (theLoad != 0)
    delete theLoad;
}

int
ShadowTruss::getNumExternalNodes(void) const
{
  return 2;
}

const ID &
ShadowTruss::getExternalNodes(void) 
{
  return externalNodes;
}


Node **
ShadowTruss::getNodePtrs(void) 
{
  theNodes[0] = end1Ptr;
  theNodes[1] = end2Ptr;
  
  return theNodes;
}

int
ShadowTruss::getNumDOF(void) {
    return 4;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void
ShadowTruss::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	return;
    }
    
    // first ensure nodes exist in Domain and set the node pointers
    int Nd1 = externalNodes(0);
    int Nd2 = externalNodes(1);
    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);	
    if (end1Ptr == 0) {
      opserr << "WARNING Truss::setDomain() - at truss " << this->getTag() << " node " <<
	Nd1 << " does not exist in domain\n";
				
      return;  // don't go any further - otherwise segemntation fault
    }
    if (end2Ptr == 0) {        
      opserr << "WARNING Truss::setDomain() - at truss " << this->getTag() << " node " <<
	Nd2 << " does not exist in domain\n";

      return;
    }	
    
    // call the DomainComponent class method THIS IS VERY IMPORTANT
    this->DomainComponent::setDomain(theDomain);

    // ensure connected nodes have correct number of dof's
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();	
    if ((dofNd1 != 2) || (dofNd2 != 2)) {
      opserr << "ShadowTruss::setDomain(): 2 dof required at nodes, " << dofNd1 << " and "
	     <<  dofNd2 << " provided\n";
      
    }	


    // now determine the length & transformation matrix
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	

    double dx = end2Crd(0)-end1Crd(0);
    double dy = end2Crd(1)-end1Crd(1);	
    
    L = sqrt(dx*dx + dy*dy);
    
    if (L == 0.0) {
      opserr << "WARNING ShadowTruss::setDomain() - ShadowTruss " << this->getTag()<< " has zero length\n";
      return;  // don't go any further - otherwise divide by 0 error
    }
	
    double cs = dx/L;
    double sn = dy/L;

    trans(0,0) = -cs;
    trans(0,1) = -sn;    
    trans(0,2) = cs;
    trans(0,3) = sn;

    // determine the nodal mass for lumped mass approach
    M = M * A * L/2;

    // create a vector to hop applied loads
    theLoad = new Vector(4);

    // send coord data to remote
    static Vector coords(4);
    coords(0) = end1Crd(0);
    coords(1) = end1Crd(1);
    coords(2) = end2Crd(0);
    coords(3) = end2Crd(1);
    msgData(0) = ShadowActorTruss_setDomain;
    this->sendID(msgData);
    this->sendVector(coords);
}   	 


int
ShadowTruss::commitState()
{
  msgData(0) = ShadowActorTruss_commitState;
  return this->sendID(msgData);
}

int
ShadowTruss::revertToLastCommit()
{
  msgData(0) = ShadowActorTruss_revertToLastCommit;
  return this->sendID(msgData);
}

int
ShadowTruss::revertToStart()
{
  msgData(0) = ShadowActorTruss_revertToStart;
  return this->sendID(msgData);
}

int
ShadowTruss::update()
{
  // determine the current strain given trial displacements at nodes
  double strain = this->computeCurrentStrain();

  // send strain to remote
  msgData(0) = ShadowActorTruss_update;
  this->sendID(msgData);
  static Vector data(1);
  data(0) = strain;
  this->sendVector(data);

  return 0;
}


const Matrix &
ShadowTruss::getTangentStiff(void)
{
  msgData(0) = ShadowActorTruss_getTangentStiff;
  this->sendID(msgData);
  this->recvMatrix(trussK);

  // return the matrix
  return trussK;
}

const Matrix &
ShadowTruss::getInitialStiff(void)

{  
  msgData(0) = ShadowActorTruss_getInitialStiff;
  this->sendID(msgData);
  this->recvMatrix(trussK);

  // return the matrix
  return trussK;
}

    
const Matrix &
ShadowTruss::getMass(void)
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
ShadowTruss::zeroLoad(void)
{
  // does nothing - no elemental loads
}

int
ShadowTruss::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // does nothing - no elemental loads
  return 0;
}


int 
ShadowTruss::addInertiaLoadToUnbalance(const Vector &accel)
{
  // check for a quick return
  if (L == 0.0 || M == 0.0) 
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = end1Ptr->getRV(accel);
  const Vector &Raccel2 = end2Ptr->getRV(accel);    

  int nodalDOF = 2;
    
#ifdef _G3DEBUG    
  if (nodalDOF != Raccel1.Size() || nodalDOF != Raccel2.Size()) {
    opserr <<"Truss::addInertiaLoadToUnbalance " <<
      "matrix and vector sizes are incompatable\n";
    return -1;
  }
#endif
    
  // want to add ( - fact * M R * accel ) to unbalance
  for (int i=0; i<2; i++) {
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
ShadowTruss::getResistingForce()
{	
  msgData(0) = ShadowActorTruss_getResistingForce;
  this->sendID(msgData);
  this->recvVector(trussR);

  return trussR;
}



int
ShadowTruss::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
ShadowTruss::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}


int
ShadowTruss::displaySelf(Renderer &theViewer, int displayMode, float fact)
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

    
    return theViewer.drawLine(v1,v2, 1.0, 1.0);
}


void
ShadowTruss::Print(OPS_Stream &s, int flag)
{
    if (flag == 0) { // print everything
      s << "Element: " << this->getTag(); 
      s << " type: ShadowTruss  iNode: " << externalNodes(0);
      s << " jNode: " << externalNodes(1);
      s << " Area: " << A;
      if (M != 0) s << " Mass (PerUnitVolume): " << M;	
      s << endln;
    }
}



Response *
ShadowTruss::setResponse(const char **argv, int argc, Information &eleInformation)
{
    //
    // we compare argv[0] for known response types for the Truss
    //

    // axial force
    if (strcmp(argv[0],"strain") ==0) 
      return new ElementResponse(this, 1, 0.0);

    else
      return 0;
}



int 
ShadowTruss::getResponse(int responseID, Information &eleInfo)
{
 double strain;
 
  switch (responseID) {
    case -1:
      return -1;
      
    case 1:
      return eleInfo.setDouble(strain);

    default:
      return 0;
  }
}

double
ShadowTruss::computeCurrentStrain(void) const
{
    // NOTE this method will never be called with L == 0.0

    // determine the strain
    const Vector &disp1 = end1Ptr->getTrialDisp();
    const Vector &disp2 = end2Ptr->getTrialDisp();	

    double dLength = 0.0;
    for (int i=0; i<2; i++){
      dLength -= (disp2(i)-disp1(i)) * trans(0,i);
    }
    
    double strain = dLength/L;

    return strain;
}


