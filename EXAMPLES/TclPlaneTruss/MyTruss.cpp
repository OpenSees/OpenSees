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
// $Date: 2003-02-19 21:48:31 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/TclPlaneTruss/MyTruss.cpp,v $                                                                        
                                                                        
// Written: fmk 
// Created: 02/99
//
// Description: This file contains the implementation for the MyTruss class.
//
// What: "@(#) MyTruss.C, revA"


// we specify what header files we need
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <OPS_Globals.h>

#include "MyTruss.h"

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>



// initialise the class wide variables
Matrix MyTruss::trussK(4,4);
Matrix MyTruss::trussM(4,4);
Vector MyTruss::trussR(4);
Node * MyTruss::theNodes[2];

// typical constructor
MyTruss::MyTruss(int tag, 
                 int Nd1, int Nd2, 
                 UniaxialMaterial &theMat, 
                 double a,
		 double m)

:Element(tag,ELE_TAG_MyTruss),     
 externalNodes(2),
 trans(1,4), L(0.0), A(a), M(m), end1Ptr(0), end2Ptr(0), theLoad(0)
{	
  // get a copy of the material object for our own use
  theMaterial = theMat.getCopy();
  if (theMaterial == 0) 
    opserr << "FATAL MyTruss::MyTruss() - out of memory, could not get a copy of the Material\n";
    
  // fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 2)
    opserr << "FATAL MyTruss::MyTruss() - out of memory, could not create an ID of size 2\n";

  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;        
}

// constructor which should be invoked by an FE_ObjectBroker only
MyTruss::MyTruss()
:Element(0,ELE_TAG_MyTruss),     
 theMaterial(0),
 externalNodes(2),
 trans(1,4), L(0.0), A(0.0), M(0.0), end1Ptr(0), end2Ptr(0), theLoad(0)
{
  if (externalNodes.Size() != 2)
    opserr << "FATAL MyTruss::MyTruss() - out of memory, could not create an ID of size 2\n";
}

//  destructor - provided to clean up any memory
MyTruss::~MyTruss()
{
    // clean up the memory associated with the element, this is
    // memory the MyTruss objects allocates and memory allocated 
    // by other objects that the MyTruss object is responsible for 
    // cleaning up, i.e. the MaterialObject.

    if (theMaterial != 0)
	delete theMaterial;    
    if (theLoad != 0)
      delete theLoad;
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


Node **
MyTruss::getNodePtrs(void) 
{
  theNodes[0] = end1Ptr;
  theNodes[1] = end2Ptr;

  return theNodes;
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
      opserr << "MyTruss::setDomain(): 2 dof required at nodes, " << dofNd1 << " and "
	     <<  dofNd2 << " provided\n";
      
    }	

    // now determine the length & transformation matrix
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	

    double dx = end2Crd(0)-end1Crd(0);
    double dy = end2Crd(1)-end1Crd(1);	
    
    L = sqrt(dx*dx + dy*dy);
    
    if (L == 0.0) {
      opserr << "WARNING MyTruss::setDomain() - MyTruss " << this->getTag()<< " has zero length\n";
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

int
MyTruss::update()
{
  // determine the current strain given trial displacements at nodes
  double strain = this->computeCurrentStrain();

  // set the strain in the materials
  theMaterial->setTrialStrain(strain);

  return 0;
}


const Matrix &
MyTruss::getTangentStiff(void)
{
    if (L == 0.0) { // length = zero - problem in setDomain() warning message already printed
	trussK.Zero();
	return trussK;
    }

    // get the current E from the material for the last updated strain
    double E = theMaterial->getTangent();

    // form the tangent stiffness matrix
    trussK = trans^trans;
    trussK *= A*E/L;  

    // return the matrix
    return trussK;
}

const Matrix &
MyTruss::getInitialStiff(void)
{
    if (L == 0.0) { // length = zero - problem in setDomain() warning message already printed
	trussK.Zero();
	return trussK;
    }

    // get the current E from the material for the last updated strain
    double E = theMaterial->getInitialTangent();

    // form the tangent stiffness matrix
    trussK = trans^trans;
    trussK *= A*E/L;  

    // return the matrix
    return trussK;
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
  // does nothing - no elemental loads
}

int
MyTruss::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // does nothing - no elemental loads
  return 0;
}


int 
MyTruss::addInertiaLoadToUnbalance(const Vector &accel)
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
      "matrix and vector sizes are incompatible\n";
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
MyTruss::getResistingForce()
{	
    if (L == 0.0) { // if length == 0, problem in setDomain()
	trussR.Zero();
	return trussR;
    }

    // want: R = Ku - Pext

    // force = F * transformation 
    double force = A*theMaterial->getStress();
    for (int i=0; i<4; i++)
	trussR(i) = trans(0,i)*force;

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
      opserr << "WARNING MyTruss::sendSelf() - failed to send Vector\n";
      return -1;
    }	      

    // MyTruss then sends the tags of it's two end nodes
    res = theChannel.sendID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING MyTruss::sendSelf() - failed to send ID\n";
      return -2;
    }

    // finally MyTruss asks it's material object to send itself
    res = theMaterial->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING MyTruss::sendSelf() - failed to send the Material\n";
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
      opserr << "WARNING MyTruss::recvSelf() - failed to receive Vector\n";
      return -1;
    }	      

    this->setTag((int)data(0));
    A = data(1);
    M = data(4);
    
    // MyTruss now receives the tags of it's two external nodes
    res = theChannel.recvID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING MyTruss::recvSelf() - failed to receive ID\n";
      return -2;
    }

    // we create a material object of the correct type,
    // sets its database tag and asks this new object to recveive itself.
    int matClass = data(2);
    int matDb = data(3);

    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      opserr << "WARNING MyTruss::recvSelf() - failed to create a Material\n";
      return -3;
    }

    // we set the dbTag before we receive the material  - this is important
    theMaterial->setDbTag(matDb); 
    res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING MyTruss::recvSelf() - failed to receive the Material\n";
      return -3;
    }

    return 0;
}


int
MyTruss::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // check setDomain() was successful
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

    if (displayMode == 3) { // use the strain as the drawing measure
	double strain = theMaterial->getStrain();
        return theViewer.drawLine(v1, v2, strain, strain);	
    } else if (displayMode == 2) { // otherwise use the material stress
        double stress = A*theMaterial->getStress();
        return theViewer.drawLine(v1,v2, stress, stress);
    } else { // use the axial force
        double force = A * theMaterial->getStress();
        return theViewer.drawLine(v1,v2, force, force);
    }
}


void
MyTruss::Print(OPS_Stream &s, int flag)
{
    // compute the strain and axial force in the member
    double strain, force;
    if (L == 0.0) {
      strain = 0;
      force = 0.0;
    } else {
      strain = theMaterial->getStrain();
      force = A * theMaterial->getStress();    
    }

    for (int i=0; i<4; i++)
      trussR(i) = trans(0,i)*force;

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
      s << endln;
    } else if (flag == 1) { // just print ele id, strain and force
      s << this->getTag() << "  " << strain << "  " << force << endln;
    }
}



Response *
MyTruss::setResponse(char **argv, int argc, Information &eleInformation)
{
    //
    // we compare argv[0] for known response types for the Truss
    //

    // axial force
    if (strcmp(argv[0],"axialForce") ==0) 
      return new ElementResponse(this, 1, 0.0);

    // a material quantity    
    else if (strcmp(argv[0],"material") == 0)
      return theMaterial->setResponse(&argv[1], argc-1, eleInformation);

    else
	return 0;
}



int 
MyTruss::getResponse(int responseID, Information &eleInfo)
{
 double strain;
 
  switch (responseID) {
    case -1:
      return -1;
      
    case 1:
      theMaterial->setTrialStrain(strain);
      return eleInfo.setDouble(A * theMaterial->getStress());

    default:
      return 0;
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
      dLength -= (disp2(i)-disp1(i)) * trans(0,i);
    }
    
    double strain = dLength/L;

    return strain;
}


