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
                                                                        
// $Revision: 1.8 $
// $Date: 2009-07-23 23:05:21 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewElement/cpp/TrussCPP.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the implementation for the TrussCPP class.
//
// What: "@(#) TrussCPP.C, revA"


// we specify what header files we need
#include "TrussCPP.h"
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
Matrix TrussCPP::trussK(4,4);
Matrix TrussCPP::trussM(4,4);
Matrix TrussCPP::trussD(4,4);
Vector TrussCPP::trussR(4);

#ifdef _USRDLL
#include <windows.h>
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif


OPS_Export void
localInit() 
{
  OPS_Error("TrussCPP element - Written by fmk UC Berkeley Copyright 2008 - Use at your Own Peril\n", 1);
}

OPS_Export void *
OPS_TrussCPP()
{
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
  Element *theTruss = new TrussCPP(eleTag, iData[1], iData[2], *theMaterial, dData[0]);


  if (theTruss == 0) {
    opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
    return 0;
  }

  return theTruss;
}


// typical constructor
TrussCPP::TrussCPP(int tag, 
                 int Nd1, int Nd2, 
                 UniaxialMaterial &theMat, 
                 double a)
:Element(tag,ELE_TAG_TrussCPP),     
 externalNodes(2),
 trans(1,4), L(0.0), A(a)
{	
  // get a copy of the material object for our own use
  theMaterial = theMat.getCopy();
  if (theMaterial == 0) {
   opserr << "FATAL TrussCPP::TrussCPP() - out of memory, could not get a copy of the Material\n";
   exit(-1);
  }
    
  // fill in the ID containing external node info with node id's    
  if (externalNodes.Size() != 2) {
    opserr << "FATAL TrussCPP::TrussCPP() - out of memory, could not create an ID of size 2\n";
    exit(-1);
  }

  externalNodes(0) = Nd1;
  externalNodes(1) = Nd2;        

  theNodes[0] = 0; 
  theNodes[1] = 0;
}

// constructor which should be invoked by an FE_ObjectBroker only
TrussCPP::TrussCPP()
:Element(0,ELE_TAG_TrussCPP),     
 theMaterial(0),
 externalNodes(2),
 trans(1,4), L(0.0), A(0.0)
{
  theNodes[0] = 0; 
  theNodes[1] = 0;
}

//  destructor - provided to clean up any memory
TrussCPP::~TrussCPP()
{
    // clean up the memory associated with the element, this is
    // memory the TrussCPP objects allocates and memory allocated 
    // by other objects that the TrussCPP object is responsible for 
    // cleaning up, i.e. the MaterialObject.

    if (theMaterial != 0)
	delete theMaterial;    
}

int
TrussCPP::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
TrussCPP::getExternalNodes(void) 
{
  return externalNodes;
}

Node **
TrussCPP::getNodePtrs(void) 
{
  return theNodes;
}

int
TrussCPP::getNumDOF(void) {
    return 4;
}

// method: setDomain()
//    to set a link to the enclosing Domain, ensure nodes exist in Domain
//    and set pointers to these nodes, also determines the length and 
//    transformation Matrix.
void
TrussCPP::setDomain(Domain *theDomain)
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
      opserr << "WARNING Truss::setDomain() - at truss " << this->getTag() << " node " <<
	Nd1 << "  does not exist in domain\n";
				
	return;  // don't go any further - otherwise segemntation fault
    }
    if (end2Ptr == 0) {        
      opserr << "WARNING Truss::setDomain() - at truss " << this->getTag() << " node " <<
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
      opserr << "TrussCPP::setDomain(): 2 dof required at nodes\n";
      return;
    }	

    // now determine the length & transformation matrix
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	

    double dx = end2Crd(0)-end1Crd(0);
    double dy = end2Crd(1)-end1Crd(1);	
    
    L = sqrt(dx*dx + dy*dy);
    
    if (L == 0.0) {
      opserr << "WARNING TrussCPP::setDomain() - TrussCPP " << this->getTag() << 
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
TrussCPP::commitState()
{
    return theMaterial->commitState();
}

int
TrussCPP::revertToLastCommit()
{
    return theMaterial->revertToLastCommit();
}

int
TrussCPP::revertToStart()
{
    return theMaterial->revertToStart();
}

int
TrussCPP::update()
{
  // determine the current strain given trial displacements at nodes
  double strain = this->computeCurrentStrain();

  // set the strain in the materials
  theMaterial->setTrialStrain(strain);

  return 0;
}


const Matrix &
TrussCPP::getTangentStiff(void)
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
TrussCPP::getInitialStiff(void)
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
TrussCPP::getDamp(void)
{
  return trussD;
}


const Matrix &
TrussCPP::getMass(void)
{ 
  trussM.Zero();
  return trussM;
}

const Vector &
TrussCPP::getResistingForce()
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
TrussCPP::sendSelf(int commitTag, Channel &theChannel)
{
    int res;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // TrussCPP packs it's data into a Vector and sends this to theChannel
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
      opserr << "WARNING TrussCPP::sendSelf() - failed to send Vector\n";
      return -1;
    }	      

    // TrussCPP then sends the tags of it's two end nodes
    res = theChannel.sendID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING TrussCPP::sendSelf() - failed to send ID\n";
      return -2;
    }

    // finally TrussCPP asks it's material object to send itself
    res = theMaterial->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING TrussCPP::sendSelf() - failed to send the Material\n";
      return -3;
    }

    return 0;
}

int
TrussCPP::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res;
    int dataTag = this->getDbTag();

    // TrussCPP creates a Vector, receives the Vector and then sets the 
    // internal data with the data in the Vector

    Vector data(5);
    res = theChannel.recvVector(dataTag, commitTag, data);
    if (res < 0) {
      opserr << "WARNING TrussCPP::recvSelf() - failed to receive Vector\n";
      return -1;
    }	      

    this->setTag((int)data(0));
    A = data(1);
    
    // TrussCPP now receives the tags of it's two external nodes
    res = theChannel.recvID(dataTag, commitTag, externalNodes);
    if (res < 0) {
      opserr << "WARNING TrussCPP::recvSelf() - failed to receive ID\n";
      return -2;
    }

    // we create a material object of the correct type,
    // sets its database tag and asks this new object to recveive itself.
    int matClass = data(2);
    int matDb = data(3);

    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      opserr << "WARNING TrussCPP::recvSelf() - failed to create a Material\n";
      return -3;
    }

    // we set the dbTag before we receive the material  - this is important
    theMaterial->setDbTag(matDb); 
    res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "WARNING TrussCPP::recvSelf() - failed to receive the Material\n";
      return -3;
    }

    return 0;
}

void
TrussCPP::Print(OPS_Stream &s, int flag)
{
  s << "Element: " << this->getTag(); 
  s << " type: TrussCPP  iNode: " << externalNodes(0);
  s << " jNode: " << externalNodes(1);
  s << " Area: " << A;
  s << " \t Material: " << *theMaterial;
}



Response *
TrussCPP::setResponse(const char **argv, int argc, OPS_Stream &s)
{
    //
    // we compare argv[0] for known response types for the Truss
    //

    // axial force
    if (strcmp(argv[0],"axialForce") ==0) 
      return new ElementResponse(this, 1, 0.0);

    // a material quantity    
    /*
    else if (strcmp(argv[0],"material") == 0)
      return theMaterial->setResponse(&argv[1], argc-1, eleInformation,s);
    */

    else
	return 0;
}



int 
TrussCPP::getResponse(int responseID, Information &eleInfo)
{
 double strain;
 
  switch (responseID) {
    case -1:
      return -1;
      
    case 1:
      return eleInfo.setDouble(A * theMaterial->getStress());

    default:
      return 0;
  }
}

double
TrussCPP::computeCurrentStrain(void) const
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


