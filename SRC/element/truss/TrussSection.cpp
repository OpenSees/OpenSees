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
                                                                        
// $Revision$
// $Date$
// $URL$
                                                                        
                                                                        
// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the implementation for the TrussSection class.
//
// What: "@(#) TrussSection.C, revA"

#include <TrussSection.h>
#include <Information.h>
#include <Parameter.h>

#include <string.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <SectionForceDeformation.h>
#include <ElementResponse.h>

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


#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_TrussSectionElement()
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 4) {
    opserr << "Invalid Args want: element TrussSection $tag $iNode $jNode $sectTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;	
  }

  int    iData[4];
  double rho = 0.0;
  int ndm = OPS_GetNDM();
  int doRayleigh = 0; // by default rayleigh not done
  int cMass = 0; // by default use lumped mass matrix

  int numData = 4;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode, sectTag) in element TrussSection " << endln;
    return 0;
  }

  SectionForceDeformation *theSection = OPS_getSectionForceDeformation(iData[3]);
    
  if (theSection == 0) {
    opserr << "WARNING: Invalid section not found element TrussSection " << iData[0] << " $iNode $jNode " << 
      iData[3] << " <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
    return 0;
  }
  
  numRemainingArgs -= 4;
  while (numRemainingArgs > 1) {
    const char *argvS = OPS_GetString();
  
    if (strcmp(argvS,"-rho") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &rho) != 0) {
	opserr << "WARNING Invalid rho in element TrussSection " << iData[0] << 
	  " $iNode $jNode $secTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	return 0;
      }
    } else if (strcmp(argvS,"-cMass") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &cMass) != 0) {
	opserr << "WARNING: Invalid cMass in element TrussSection " << iData[0] << 
	  " $iNode $jNode $sectTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	return 0;
      }
    } else if (strcmp(argvS,"-doRayleigh") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &doRayleigh) != 0) {
	opserr << "WARNING: Invalid doRayleigh in element TrussSection " << iData[0] << 
	  " $iNode $jNode $sectTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
	return 0;
      }
    } else {
      opserr << "WARNING: Invalid option " << argvS << "  in: element TrussSection " << iData[0] << 
	" $iNode $jNode $secTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
      return 0;
    }
    numRemainingArgs -= 2;
  }

  // now create the TrussSection
  theElement = new TrussSection(iData[0], ndm, iData[1], iData[2], *theSection, rho, doRayleigh, cMass);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element TrussSection " << iData[0] << 
      " $iNode $jNode $secTag <-rho $rho> <-cMass $flag> <-doRayleigh $flag>\n";
  }

  return theElement;
}


TrussSection::TrussSection(int tag, int dim,
			   int Nd1, int Nd2,
			   SectionForceDeformation &theSect,
			   double r, int damp, int cm)
:Element(tag,ELE_TAG_TrussSection),     
  connectedExternalNodes(2),
  dimension(dim), numDOF(0),
  theLoad(0), theMatrix(0), theVector(0),
  L(0.0), rho(r), doRayleighDamping(damp),
  cMass(cm), theSection(0), initialDisp(0)
{
    // get a copy of the material and check we obtained a valid copy
    theSection = theSect.getCopy();
    if (theSection == 0) {
      opserr << "FATAL TrussSection::TrussSection - failed to get a copy of material " << 
	theSect.getTag() << endln;
      exit(-1);
    }
    int order = theSection->getOrder();
    const ID &code = theSection->getType();
    
    int i;
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_P)
	break;
    
    if (i == order)
      opserr << "TrussSection::TrussSection - section does not provide axial response\n";

    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL TrussSection::TrussSection - failed to create an ID of correct size\n";
      exit(-1);
    }

    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;        

    // set node pointers to NULL
    for (i=0; i<2; i++)
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
TrussSection::TrussSection()
:Element(0,ELE_TAG_TrussSection),     
 connectedExternalNodes(2),
  dimension(0), numDOF(0),
  theLoad(0), theMatrix(0), theVector(0),
  L(0.0), rho(0.0), doRayleighDamping(0),
  cMass(0), theSection(0), initialDisp(0)
{
    // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL TrussSection::TrussSection - failed to create an ID of correct size\n";
      exit(-1);
  }

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

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
TrussSection::~TrussSection()
{
  if (theSection != 0)
    delete theSection;
  if (theLoad != 0)
    delete theLoad;
  if (theLoadSens != 0)
    delete theLoadSens;  
  if (initialDisp != 0)
    delete [] initialDisp;
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

Node **
TrussSection::getNodePtrs(void) 
{
  return theNodes;
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

    // if nodes not in domain, warning message & set default numDOF as 2
    if ((theNodes[0] == 0) || (theNodes[1] == 0)){
      if (theNodes[0] == 0)
        opserr << "TrussSection::setDomain() - Nd1: " << Nd1 << " does not exist in Domain\n";
      else
        opserr << "TrussSection::setDomain() - Nd1: " << Nd2 << " does not exist in Domain\n";

      opserr << " for truss with id " << this->getTag() << endln;

      // fill this in so don't segment fault later
      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	

      return;
    }

    // now determine the number of dof and the dimesnion    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    if (dofNd1 != dofNd2) {
      opserr << "WARNING TrussSection::setDomain(): nodes " << Nd1 << " and " <<
	Nd2 << "have differing dof at ends for truss " << this->getTag() << endln;	

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
      opserr << "WARNING TrussSection::setDomain cannot handle " << dimension << 
	" dofs at nodes in " << dofNd1 << " d problem\n"; 

      numDOF = 2;    
      theMatrix = &trussM2;
      theVector = &trussV2;	
      return;
    }

    // now determine the length, cosines and fill in the transformation
    // NOTE t = -t(every one else uses for residual calc)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();

    if (dimension == 1) {
      double dx = end2Crd(0)-end1Crd(0);	
      if (initialDisp == 0) {
	double iDisp = end2Disp(0)-end1Disp(0);
	
	if (iDisp != 0) {
	  initialDisp = new double[1];
	  initialDisp[0] = iDisp;
	  dx += iDisp;
	}
      }
      
      L = sqrt(dx*dx);
      
      if (L == 0.0) {
	opserr << "WARNING TrussSection::setDomain() - truss " << this->getTag() << " has zero length\n";
	return;
      }	

	cosX[0] = 1.0;

    } else if (dimension == 2) {
      double dx = end2Crd(0)-end1Crd(0);
      double dy = end2Crd(1)-end1Crd(1);	

      if (initialDisp == 0) {
	double iDispX = end2Disp(0)-end1Disp(0);
	double iDispY = end2Disp(1)-end1Disp(1);
	if (iDispX != 0 || iDispY != 0) {
	  initialDisp = new double[2];
	  initialDisp[0] = iDispX;
	  initialDisp[1] = iDispY;
	  dx += iDispX;
	  dy += iDispY;
	}
      }
    
      L = sqrt(dx*dx + dy*dy);
      
      if (L == 0.0) {
	opserr << "WARNING TrussSection::setDomain() - truss " << this->getTag() << " has zero length\n";
	return;
      }
	
      cosX[0] = dx/L;
      cosX[1] = dy/L;

    } else {
      double dx = end2Crd(0)-end1Crd(0);
      double dy = end2Crd(1)-end1Crd(1);	
      double dz = end2Crd(2)-end1Crd(2);		
      
      if (initialDisp == 0) {
	double iDispX = end2Disp(0)-end1Disp(0);
	double iDispY = end2Disp(1)-end1Disp(1);      
	double iDispZ = end2Disp(2)-end1Disp(2);      
	if (iDispX != 0 || iDispY != 0 || iDispZ != 0) {
	  initialDisp = new double[3];
	  initialDisp[0] = iDispX;
	  initialDisp[1] = iDispY;
	  initialDisp[2] = iDispZ;
	  dx += iDispX;
	  dy += iDispY;
	  dz += iDispZ;
	}
      }

      L = sqrt(dx*dx + dy*dy + dz*dz);
      
      if (L == 0.0) {
	opserr << "WARNING TrussSection::setDomain() - truss " << this->getTag() << " has zero length\n";
	return;
      }
	
	cosX[0] = dx/L;
	cosX[1] = dy/L;
	cosX[2] = dz/L;	
    }

    // create the load vector
    if (theLoad == 0)
      theLoad = new Vector(numDOF);
    else if (theLoad->Size() != numDOF) {
      delete theLoad;
      theLoad = new Vector(numDOF);
    }

    if (theLoad == 0) {
      opserr << "TrussSection::setDomain - truss " << this->getTag() << 
	"out of memory creating vector of size" << numDOF << endln;
      exit(-1);
      return;
    }          
    
    this->update();
}   	 


int
TrussSection::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "TrussSection::commitState () - failed in base class";
  }    
  retVal = theSection->commitState();
  return retVal;
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




int
TrussSection::update()
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	return -1;
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
    
    return theSection->setTrialSectionDeformation(e);
}


const Matrix &
TrussSection::getTangentStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    int order = theSection->getOrder();
    const ID &code = theSection->getType();
	
    const Matrix &k = theSection->getSectionTangent();
    double AE = 0.0;
    int i;
    for (i = 0; i < order; i++) {
      if (code(i) == SECTION_RESPONSE_P)
	AE += k(i,i);
    }

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;

    int numDOF2 = numDOF/2;
    double temp;
    AE /= L;
    for (i = 0; i < dimension; i++) {
      for (int j = 0; j < dimension; j++) {
	temp = cosX[i]*cosX[j]*AE;
	stiff(i,j) = temp;
	stiff(i+numDOF2,j) = -temp;
	stiff(i,j+numDOF2) = -temp;
	stiff(i+numDOF2,j+numDOF2) = temp;
      }
    }

    return *theMatrix;
}

const Matrix &
TrussSection::getInitialStiff(void)
{
    if (L == 0.0) { // - problem in setDomain() no further warnings
	theMatrix->Zero();
	return *theMatrix;
    }
    
    int order = theSection->getOrder();
    const ID &code = theSection->getType();
	
    const Matrix &k = theSection->getInitialTangent();
    double AE = 0.0;
    int i;
    for (i = 0; i < order; i++) {
      if (code(i) == SECTION_RESPONSE_P)
	AE += k(i,i);
    }

    // come back later and redo this if too slow
    Matrix &stiff = *theMatrix;

    int numDOF2 = numDOF/2;
    double temp;
    AE /= L;
    for (i = 0; i < dimension; i++) {
      for (int j = 0; j < dimension; j++) {
	temp = cosX[i]*cosX[j]*AE;
	stiff(i,j) = temp;
	stiff(i+numDOF2,j) = -temp;
	stiff(i,j+numDOF2) = -temp;
	stiff(i+numDOF2,j+numDOF2) = temp;
      }
    }

    return *theMatrix;
}
    
const Matrix &
TrussSection::getDamp(void)
{   
  if (doRayleighDamping == 1)
    return this->Element::getDamp();
  
  theMatrix->Zero();
  return *theMatrix;
}


const Matrix &
TrussSection::getMass(void)
{
  // zero the matrix
  Matrix &mass = *theMatrix;
  mass.Zero();    
  
  // check for quick return
  if (L == 0.0 || rho == 0.0) { // - problem in setDomain() no further warnings
    return mass;
  }
  
  if (cMass == 0)  {
    // lumped mass matrix
    double m = 0.5*rho*L;
    int numDOF2 = numDOF/2;
    for (int i = 0; i < dimension; i++) {
      mass(i,i) = m;
      mass(i+numDOF2,i+numDOF2) = m;
    }
  } else  {
    // consistent mass matrix
    double m = rho*L/6.0;
    int numDOF2 = numDOF/2;
    for (int i = 0; i < dimension; i++) {
      mass(i,i) = 2.0*m;
      mass(i,i+numDOF2) = m;
      mass(i+numDOF2,i) = m;
      mass(i+numDOF2,i+numDOF2) = 2.0*m;
    }
  }
  
  return mass;
}



void 
TrussSection::zeroLoad(void)
{
  theLoad->Zero();
}


int 
TrussSection::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "TrussSection::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  return -1;
}


int 
TrussSection::addInertiaLoadToUnbalance(const Vector &accel)
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
    opserr <<"TrussSection::addInertiaLoadToUnbalance " <<
      "matrix and vector sizes are incompatible\n";
    return -1;
  }
#endif
  
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    double m = 0.5*rho*L;
    for (int i=0; i<dimension; i++) {
      (*theLoad)(i) -= m*Raccel1(i);
      (*theLoad)(i+nodalDOF) -= m*Raccel2(i);
    }
  } else  {
    double m = rho*L/6.0;
    for (int i=0; i<dimension; i++) {
      (*theLoad)(i) -= 2.0*m*Raccel1(i) + m*Raccel2(i);
      (*theLoad)(i+nodalDOF) -= m*Raccel1(i) + 2.0*m*Raccel2(i);
    }
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
    
    int order = theSection->getOrder();
    const ID &code = theSection->getType();
	
    const Vector &s = theSection->getStressResultant();
    double force = 0.0;
    int i;
    for (i = 0; i < order; i++) {
      if (code(i) == SECTION_RESPONSE_P)
	force += s(i);
    }

    int numDOF2 = numDOF/2;
    double temp;
    for (i = 0; i < dimension; i++) {
      temp = cosX[i]*force;
      (*theVector)(i) = -temp;
      (*theVector)(i+numDOF2) = temp;
    }

    // subtract external load
    (*theVector) -= *theLoad;
  
    return *theVector;
}



const Vector &
TrussSection::getResistingForceIncInertia()
{	
  this->getResistingForce();
  
  // now include the mass portion
  if (L != 0.0 && rho != 0.0) {
    
    // add inertia forces from element mass
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();	
    
    int numDOF2 = numDOF/2;
    
    if (cMass == 0)  {
      // lumped mass matrix
      double m = 0.5*rho*L;
      for (int i = 0; i < dimension; i++) {
        (*theVector)(i) += m*accel1(i);
        (*theVector)(i+numDOF2) += m*accel2(i);
      }
    } else  {
      // consistent mass matrix
      double m = rho*L/6.0;
      for (int i=0; i<dimension; i++) {
        (*theVector)(i) += 2.0*m*accel1(i) + m*accel2(i);
        (*theVector)(i+numDOF2) += m*accel1(i) + 2.0*m*accel2(i);
      }
    }
    
    // add the damping forces if rayleigh damping
    if (doRayleighDamping == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
      theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
  } else {
    
    // add the damping forces if rayleigh damping
    if (doRayleighDamping == 1 && (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
      theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
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

  static Vector data(11);
  data(0) = this->getTag();
  data(1) = dimension;
  data(2) = numDOF;
  data(5) = rho;
  data(6) = doRayleighDamping;
  data(7) = cMass;

  data(3) = theSection->getClassTag();
  int matDbTag = theSection->getDbTag();

  // NOTE: we do have to ensure that the Section has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theSection->setDbTag(matDbTag);
  }
  data(4) = matDbTag;

  if (initialDisp != 0) {
    for (int i=0; i<dimension; i++) {
      data[8+i] = initialDisp[i];
    }
  }

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING TrussSection::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }	      

  // truss then sends the tags of it's two end nodes
  res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr << "WARNING TrussSection::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return -2;
  }

  // finally truss asks it's Section object to send itself
  res = theSection->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "WARNING TrussSection::sendSelf() - " << this->getTag() << " failed to send its Section\n";
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

  static Vector data(11);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING TrussSection::recvSelf() - failed to receive Vector\n";
    return -1;
  }	      

  this->setTag((int)data(0));
  dimension = (int)data(1);
  numDOF = (int)data(2);
  rho = data(5);
  doRayleighDamping = (int)data(6);
  cMass = (int)data(7);

  initialDisp = new double[dimension];
  for (int i=0; i<dimension; i++)
    initialDisp[i] = 0.0;
  
  int initial = 0;
  for (int i=0; i<dimension; i++) {
    if (data(8+i) != 0.0) {
      initial = 1;
    }
  }
  
  if (initial != 0) {
    for (int i=0; i<dimension; i++) {
      initialDisp[i] = data(8+i);
    }    
  }

  // truss now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr << "WARNING TrussSection::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally truss creates a new section object of the correct type,
  // sets its database tag and asks this new object to recveive itself.

  int sectClass = (int)data(3);
  int sectDb = (int)data(4);

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
    opserr << "WARNING TrussSection::recvSelf() - " << this->getTag() << 
      " failed to get a blank Section of type " << sectClass << endln;
    return -3;
  }

  theSection->setDbTag(sectDb); // note: we set the dbTag before we receive the Section
  res = theSection->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "WARNING TrussSection::recvSelf() - " << this->getTag() << " failed to receive its Section\n";
    return -3;
  }

  return 0;
}


int
TrussSection::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **mode, int numMode)
{
    // ensure setDomain() worked
    if (L == 0.0)
       return 0;

    // first determine the two end points of the truss based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    static Vector v1(3);
    static Vector v2(3);
    float d1 = 0.0;
    float d2 = 0.0;

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    if (displayMode == 1 || displayMode == 2) {

        // compute the strain and axial force in the member
        double strain, force;
        if (L == 0.0) {
            strain = 0.0;
            force = 0.0;
        } else {
	        strain = this->computeCurrentStrain();
            force = 0.0;
	
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
	        return theViewer.drawLine(v1, v2, (float)strain, (float)strain, this->getTag());	
        else // otherwise use the axial force as measure
	        return theViewer.drawLine(v1,v2, (float)force, (float)force, this->getTag());

    }
    else {
        return theViewer.drawLine(v1, v2, 0.0, 0.0, this->getTag());
    }
}


void
TrussSection::Print(OPS_Stream &s, int flag)
{
    // compute the strain and axial force in the member
    double strain, force;
    if (L == 0.0) {
	strain = 0;
	force = 0.0;
    } else {
		strain = this->computeCurrentStrain();
        force = 0.0;

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
    
    double temp;
    int numDOF2 = numDOF/2;
	if (theVector != 0) {
		for (int i=0; i<dimension; i++) {
			temp = force*cosX[i];
			(*theVector)(i) = -force;
			(*theVector)(i+numDOF2) = force;
		}
	}
     
    if (flag == OPS_PRINT_CURRENTSTATE) { // print everything
        s << "Element: " << this->getTag();
        s << " type: TrussSection  iNode: " << connectedExternalNodes(0);
        s << " jNode: " << connectedExternalNodes(1);
        s << " Mass density/length: " << rho;
        s << " cMass: " << cMass;
        
        s << " \n\t strain: " << strain;
        s << " axial load: " << force;
        if (theVector != 0)
            s << " \n\t unbalanced load: " << *theVector;
        s << " \t Section: " << *theSection;
        s << endln;
    }
    
    if (flag == 1) {
        s << this->getTag() << "  " << strain << "  ";
        s << force << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"TrussSection\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"massperlength\": " << rho << ", ";
        s << "\"section\": \"" << theSection->getTag() << "\"}";
    }
}

double
TrussSection::computeCurrentStrain(void) const
{
    // NOTE method will not be called if L == 0

    // determine the strain
    const Vector &disp1 = theNodes[0]->getTrialDisp();
    const Vector &disp2 = theNodes[1]->getTrialDisp();	

    double dLength = 0.0;
    if (initialDisp == 0)
      for (int i = 0; i < dimension; i++)
	dLength += (disp2(i)-disp1(i))*cosX[i];
    else
      for (int i = 0; i < dimension; i++)
	dLength += (disp2(i)-disp1(i)-initialDisp[i])*cosX[i];

    // this method should never be called with L == 0
    return dLength/L;
}

Response*
TrussSection::setResponse(const char **argv, int argc, OPS_Stream &output)
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

    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForce") == 0) || (strcmp(argv[0],"globalForces") == 0)){
            char outputData[10];
            int numDOFperNode = numDOF/2;
            for (int i=0; i<numDOFperNode; i++) {
                sprintf(outputData,"P1_%d", i+1);
                output.tag("ResponseType", outputData);
            }
            for (int j=0; j<numDOFperNode; j++) {
                sprintf(outputData,"P2_%d", j+1);
                output.tag("ResponseType", outputData);
            }
            theResponse =  new ElementResponse(this, 1, Vector(numDOF));

    } else if ((strcmp(argv[0],"localForce") == 0) || (strcmp(argv[0],"localForces") == 0) ) {
            theResponse =  new ElementResponse(this, 11, Vector(numDOF));
	    
    } else if ((strcmp(argv[0],"axialForce") == 0) || 
	       (strcmp(argv[0],"basicForce") == 0) || 
	       (strcmp(argv[0],"basicForces") == 0)) {
      output.tag("ResponseType", "N");
      theResponse =  new ElementResponse(this, 2, Vector(1));

    } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformation") == 0 ||
        strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"basicDefo") == 0 ||
        strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0) {

            output.tag("ResponseType", "U");
            theResponse = new ElementResponse(this, 3, Vector(1));

    } else if (strcmp(argv[0],"basicStiffness") == 0) {

            output.tag("ResponseType", "K");
            theResponse = new ElementResponse(this, 4, Matrix(1,1));

    // a section quantity
    } else if (strcmp(argv[0], "section") == 0) {
        output.tag("GaussPointOutput");
        output.attr("number", 1);
        output.attr("eta", 0.0);
        
        if (argc > 1) {
            // we need at least one more argument otherwise 
            // there is no need to forward this call to the section
            if (argc > 2) {
                // if we have 2 or more extra arguments, the first one 
                // could be an integer. In this case we check to see if it is the section id
                // (only 1 in this case)
                int sectionNum = atoi(argv[1]);
                if (sectionNum == 0) {
                    // if it is not a number we forward the call to the section as usual
                    theResponse = theSection->setResponse(&argv[1], argc - 1, output);
                }
                else {
                    // it is a number. Now we have to make sure it is within the allowed range
                    // for this element (in this case it can only be 1)
                    // If it is > 1, then we MUST return NULL, because the MPCO recorder iteratively
                    // uses this call to understand how many fibers we have in a section
                    if (sectionNum == 1) {
                        theResponse = theSection->setResponse(&argv[2], argc - 2, output);
                    }
                }
            }
            else {
                // otherwise forward it as usual
                theResponse = theSection->setResponse(&argv[1], argc - 1, output);
            }
        }
        output.endTag();
    }

    output.endTag();
    return theResponse;

}

int 
TrussSection::getResponse(int responseID, Information &eleInfo)
{
    double strain, force;
    static Vector sVec(1);
    static Vector fVec(1);
    static Matrix kVec(1,1);

    switch (responseID) {
  case 1:
      return eleInfo.setVector(this->getResistingForce());

    case 11: {
      Vector P(numDOF);
      int order = theSection->getOrder();
      const ID &code = theSection->getType();

      const Vector &s = theSection->getStressResultant();
      force = 0.0;
      int i;
      for (i = 0; i < order; i++) {
	if (code(i) == SECTION_RESPONSE_P)
	  force += s(i);
      }

      P(numDOF/2) = force;
      P(0) = -P(numDOF/2);
      return eleInfo.setVector(P);
    }
  case 2:
      if (L == 0.0) {
          strain = 0.0;
          force = 0.0;
      } else {

          int order = theSection->getOrder();
          const ID &code = theSection->getType();

          const Vector &s = theSection->getStressResultant();
          force = 0.0;
          int i;
          for (i = 0; i < order; i++) {
              if (code(i) == SECTION_RESPONSE_P)
                  force += s(i);
          }

      }      
      fVec(0) = force;
      return eleInfo.setVector(fVec);    

  case 3:
      if (L == 0.0) {
          strain = 0.0;
      } else {
          strain = this->computeCurrentStrain();
      }
      sVec(0) = L*strain;
      return eleInfo.setVector(sVec);

    case 4:
      if (L == 0.0) {
          force = 0.0;
      } else {

          int order = theSection->getOrder();
          const ID &code = theSection->getType();

          const Matrix &ks = theSection->getSectionTangent();
          force = 0.0;
          int i;
          for (i = 0; i < order; i++) {
	    if (code(i) == SECTION_RESPONSE_P)
	      force += ks(i,i);
          }

      }      
      kVec(0,0) = force/L;
      return eleInfo.setMatrix(kVec);

  default:
      return -1;
    }
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
TrussSection::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;
  
  // Mass densitity of the truss
  if (strcmp(argv[0],"rho") == 0)
    return param.addObject(2, this);
  
  // Explicit specification of a material parameter
  if (strstr(argv[0],"material") != 0 || strstr(argv[0],"section") != 0) {
    
    if (argc < 2)
      return -1;
    
    else
      return theSection->setParameter(&argv[1], argc-1, param);
  } 
  
  // Otherwise, send it to the material
  else
    return theSection->setParameter(argv, argc, param);
}

int
TrussSection::updateParameter (int parameterID, Information &info)
{
  switch (parameterID) {
  case 2:
    rho = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
TrussSection::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}


const Matrix &
TrussSection::getKiSensitivity(int gradIndex)
{
  Matrix &stiff = *theMatrix;
  stiff.Zero();
    
  if (parameterID == 0) {
  }
  else if (parameterID == 2) {
    // Nothing here when 'rho' is random
  }
  else {
    int order = theSection->getOrder();
    const ID &code = theSection->getType();
	
    const Matrix &k = theSection->getInitialTangentSensitivity(gradIndex);
    double AdEdh = 0.0;
    int i;
    for (i = 0; i < order; i++) {
      if (code(i) == SECTION_RESPONSE_P)
	AdEdh += k(i,i);
    }

    int numDOF2 = numDOF/2;
    double temp;
    double EAoverL = AdEdh/L;
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
TrussSection::getMassSensitivity(int gradNumber)
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
TrussSection::getResistingForceSensitivity(int gradIndex)
{
	theVector->Zero();

	// Initial declarations
	int i;
	double stressSensitivity, temp1, temp2;

	// Make sure the material is up to date
	double strain = this->computeCurrentStrain();
	//double rate = this->computeCurrentStrainRate();
	//theMaterial->setTrialStrain(strain,rate);

	// Contribution from material
	//stressSensitivity = theMaterial->getStressSensitivity(gradIndex,true);

	int order = theSection->getOrder();
	const ID &code = theSection->getType();

	const Vector &dsdh = theSection->getStressResultantSensitivity(gradIndex, true);
	double dNdh = 0.0;
	for (i = 0; i < order; i++) {
	  if (code(i) == SECTION_RESPONSE_P)
	    dNdh += dsdh(i);
	}

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

		//double materialTangent = theMaterial->getTangent();
		const Matrix &ks = theSection->getSectionTangent();
		double EA = 0.0;
		for (i = 0; i < order; i++) {
		  if (code(i) == SECTION_RESPONSE_P)
		    EA += ks(i,i);
		}

		double strainSensitivity = 0.0;

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
		//stressSensitivity += materialTangent * strainSensitivity;
		stressSensitivity += EA * strainSensitivity;
	}


	// Compute sensitivity depending on 'parameter'
	//double stress = theMaterial->getStress();
	double N = 0.0;
	const Vector &s = theSection->getStressResultant();
	for (i = 0; i < order; i++) {
	  if (code(i) == SECTION_RESPONSE_P)
	    N += s(i);
	}

	int numDOF2 = numDOF/2;
	double temp;
	if (parameterID == 1) {			// Cross-sectional area

	}
	else {		// Density, material parameter or nodal coordinate
	  for (i = 0; i < dimension; i++) {
	    temp = dNdh*cosX[i] + N*dcosXdh[i];
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
TrussSection::commitSensitivity(int gradIndex, int numGrads)
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
	  sens1 = theNodes[0]->getDispSensitivity(i+1, gradIndex);
	  sens2 = theNodes[1]->getDispSensitivity(i+1, gradIndex);
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
	int order = theSection->getOrder();
	const ID &code = theSection->getType();

	Vector dedh(order);
	for (int i = 0; i < order; i++)
	  if (code(i) == SECTION_RESPONSE_P)
	    dedh(i) = strainSensitivity;

	return theSection->commitSensitivity(dedh, gradIndex, numGrads);
}

int 
TrussSection::addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool somethingRandomInMotions)
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
	"matrix and vector sizes are incompatible\n";
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
	"matrix and vector sizes are incompatible\n";
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

// AddingSensitivity:END /////////////////////////////////////////////
