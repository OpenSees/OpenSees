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
                                                                        
// Written: fmk 
//
// What: "@(#) Truss2D.C, revA"


// we specify what header files we need
#include "Truss2D.h"
#include <elementAPI.h>
#include <G3Globals.h>

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// initialise the class wide variables
Matrix Truss2D::trussK(4,4);
Vector Truss2D::trussR(4);

#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numMyTruss = 0;

OPS_Export void *
OPS_Truss2D()
{
  // print out a message about who wrote this element & any copyright info wanted
  if (numMyTruss == 0) {
    opserr << "Truss2D element - Written by fmk UC Berkeley Copyright 2008 - Use at your Own Peril\n";
    numMyTruss++;
  }

  Element *theTruss = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  if (numRemainingArgs == 0) { // parallel processing
    theTruss = new Truss2D();
    return theTruss;
  }

  if (numRemainingArgs != 5) {
    opserr << "ERROR - Truss2D not enough args provided, want: element Truss2D tag? iNode? jNode? Area? matTag?\n";
    numMyTruss++;
  }

  // get the id and end nodes 
  int iData[4];
  double dData[1];
  int numData;

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element area for element" << eleTag << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr << "WARNING error reading element material tag for element " << eleTag << endln;
    return 0;
  }

  int matID = iData[3];
  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matID);
  
  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matID << "not found for element " << eleTag << endln;
    return 0;
  }

  // now create the truss and add it to the Domain

  theTruss = new Truss2D(eleTag, iData[1], iData[2], *theMaterial, dData[0]);

  if (theTruss == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    delete theMaterial;
    return 0;
  }

  return theTruss;
}


// typical constructor
Truss2D::Truss2D(int tag, 
		   int Nd1, int Nd2, 
		   UniaxialMaterial &theMat, 
		   double a)
:Element(tag, 0),     
 externalNodes(2),
 trans(1,4), L(0.0), A(a)
{	
  // get a copy of the material object for our own use
  theMaterial = theMat.getCopy();
  if (theMaterial == 0) {
   opserr << "FATAL Truss2D::Truss2D() - out of memory, could not get a copy of the Material\n";
   exit(-1);
  }
    
  // fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 2) {
    opserr << "FATAL Truss2D::Truss2D() - out of memory, could not create an ID of size 2\n";
    exit(-1);
  }

  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;        

  theNodes[0] = 0; 
  theNodes[1] = 0;
}

// constructor which should be invoked by an FE_ObjectBroker only
Truss2D::Truss2D()
:Element(0, 0),     
 theMaterial(0),
 externalNodes(2),
 trans(1,4), L(0.0), A(0.0)
{
  theNodes[0] = 0; 
  theNodes[1] = 0;
}

//  destructor - provided to clean up any memory
Truss2D::~Truss2D()
{
    // clean up the memory associated with the element, this is
    // memory the Truss2D objects allocates and memory allocated 
    // by other objects that the Truss2D object is responsible for 
    // cleaning up, i.e. the MaterialObject.

    if (theMaterial != 0)
	delete theMaterial;    
}

int
Truss2D::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
Truss2D::getExternalNodes(void) 
{
  return externalNodes;
}

Node **
Truss2D::getNodePtrs(void) 
{
  return theNodes;
}

int
Truss2D::getNumDOF(void) {
    return 4;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void
Truss2D::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	return;
    }
    
    // first ensure nodes exist in Domain and set the node pointers
    Node *end1Ptr, *end2Ptr;
    int Nd1 = externalNodes(0);
    int Nd2 = externalNodes(1);
    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);	
    if (end1Ptr == 0) {
      opserr << "WARNING Truss2D::setDomain() - at truss " << this->getTag() << " node " <<
	Nd1 << "  does not exist in domain\n";
				
	return;  // don't go any further - otherwise segemntation fault
    }
    if (end2Ptr == 0) {        
      opserr << "WARNING Truss2D::setDomain() - at truss " << this->getTag() << " node " <<
	Nd2 << "  does not exist in domain\n";

	return;  // don't go any further - otherwise segemntation fault
    }	
    theNodes[0] = end1Ptr;
    theNodes[1] = end2Ptr;
    // call the DomainComponent class method THIS IS VERY IMPORTANT
    this->DomainComponent::setDomain(theDomain);

    // ensure connected nodes have correct number of dof's
    int dofNd1 = end1Ptr->getNumberDOF();
    int dofNd2 = end2Ptr->getNumberDOF();	
    if ((dofNd1 != 2) || (dofNd2 != 2)) {
      opserr << "Truss2D::setDomain(): 2 dof required at nodes\n";
      return;
    }	

    // now determine the length & transformation matrix
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	

    double dx = end2Crd(0)-end1Crd(0);
    double dy = end2Crd(1)-end1Crd(1);	
    
    L = sqrt(dx*dx + dy*dy);
    
    if (L == 0.0) {
      opserr << "WARNING Truss2D::setDomain() - Truss2D " << this->getTag() << 
	" has zero length\n";
      return;  // don't go any further - otherwise divide by 0 error
    }
	
    double cs = dx/L;
    double sn = dy/L;

    trans(0,0) = -cs;
    trans(0,1) = -sn;    
    trans(0,2) = cs;
    trans(0,3) = sn;
}   	 


int
Truss2D::commitState()
{
    return theMaterial->commitState();
}

int
Truss2D::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}

int
Truss2D::revertToStart()
{
    return theMaterial->revertToStart();
}

int
Truss2D::update()
{
  // determine the current strain given trial displacements at nodes
  double strain = this->computeCurrentStrain();

  // set the strain in the materials
  theMaterial->setTrialStrain(strain);

  return 0;
}


const Matrix &
Truss2D::getTangentStiff(void)
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
Truss2D::getInitialStiff(void)
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

const Vector &
Truss2D::getResistingForce()
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
Truss2D::sendSelf(int commitTag, Channel &theChannel)
{
    int res;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // Truss2D packs it's data into a Vector and sends this to theChannel
    // along with it's dbTag and the commitTag passed in the arguments

    Vector data(5);
    data(0) = this->getTag();
    data(1) = A;
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
      opserr << "WARNING Truss2D::sendSelf() - failed to send Vector\n";
      return -1;
    }	      

    // Truss2D then sends the tags of it's two end nodes
    res = theChannel.sendID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING Truss2D::sendSelf() - failed to send ID\n";
      return -2;
    }

    // finally Truss2D asks it's material object to send itself
    res = theMaterial->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING Truss2D::sendSelf() - failed to send the Material\n";
      return -3;
    }

    return 0;
}

int
Truss2D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res;
    int dataTag = this->getDbTag();

    // Truss2D creates a Vector, receives the Vector and then sets the 
    // internal data with the data in the Vector

    Vector data(5);
    res = theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
      opserr << "WARNING Truss2D::recvSelf() - failed to receive Vector\n";
      return -1;
    }	      

    this->setTag((int)data(0));
    A = data(1);
    
    // Truss2D now receives the tags of it's two external nodes
    res = theChannel.recvID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING Truss2D::recvSelf() - failed to receive ID\n";
      return -2;
    }

    // we create a material object of the correct type,
    // sets its database tag and asks this new object to recveive itself.
    int matClass = data(2);
    int matDb = data(3);

    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      opserr << "WARNING Truss2D::recvSelf() - failed to create a Material\n";
      return -3;
    }

    // we set the dbTag before we receive the material  - this is important
    theMaterial->setDbTag(matDb); 
    res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING Truss2D::recvSelf() - failed to receive the Material\n";
      return -3;
    }

    return 0;
}

void
Truss2D::Print(OPS_Stream &s, int flag)
{
  s << "Element: " << this->getTag(); 
  s << " type: Truss2D  iNode: " << externalNodes(0);
  s << " jNode: " << externalNodes(1);
  s << " Area: " << A;
  s << " \t Material: " << *theMaterial;
}



Response *
Truss2D::setResponse(const char **argv, int argc, OPS_Stream &output)
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
  } else if (strcmp(argv[0],"axialForce") ==0) 
      return new ElementResponse(this, 3, 0.0);

  output.endTag();
  return theResponse;
}



int 
Truss2D::getResponse(int responseID, Information &eleInfo)
{
 double strain;
 
  switch (responseID) {
  case -1:
    return -1;
  case 1: // global forces
    return eleInfo.setVector(this->getResistingForce());
  case 2:
    return eleInfo.setVector(this->getRayleighDampingForces());
  case 3:
    theMaterial->setTrialStrain(strain);
    return eleInfo.setDouble(A * theMaterial->getStress());
  default:
    return 0;
  }
}

double
Truss2D::computeCurrentStrain(void) const
{
    // NOTE this method will never be called with L == 0.0

    // determine the strain
    const Vector &disp1 = theNodes[0]->getTrialDisp();
    const Vector &disp2 = theNodes[1]->getTrialDisp();	

    double dLength = 0.0;
    for (int i=0; i<2; i++){
      dLength -= (disp2(i)-disp1(i)) * trans(0,i);
    }
    
    double strain = dLength/L;

    return strain;
}


