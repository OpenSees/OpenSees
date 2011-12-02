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
                                                                        
// $Revision: 1.16 $
// $Date: 2006-03-21 22:19:12 $
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

TrussSection::TrussSection(int tag, 
			   int dim,
			   int Nd1, int Nd2, 
			   SectionForceDeformation &theSect,
			   double r)
:Element(tag,ELE_TAG_TrussSection),     
  connectedExternalNodes(2),
  dimension(dim), numDOF(0), theLoad(0), 
 theMatrix(0), theVector(0),
  L(0.0), rho(r), theSection(0)
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
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
TrussSection::TrussSection()
:Element(0,ELE_TAG_TrussSection),     
 connectedExternalNodes(2),
  dimension(0), numDOF(0), theLoad(0),
 theMatrix(0), theVector(0),
 L(0.0), rho(0.0), theSection(0)
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
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
TrussSection::~TrussSection()
{
    if (theSection != 0)
	delete theSection;
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

    if (dimension == 1) {
	double dx = end2Crd(0)-end1Crd(0);	
	L = sqrt(dx*dx);
	
	if (L == 0.0) {
	  opserr << "WARNING TrussSection::setDomain() - truss " << this->getTag() << " has zero length\n";
	  return;
	}	

	cosX[0] = 1.0;

    } else if (dimension == 2) {
	double dx = end2Crd(0)-end1Crd(0);
	double dy = end2Crd(1)-end1Crd(1);	
    
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
TrussSection::getMass(void)
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
      opserr << "TrussSection::addInertiaLoadToUnbalance " <<
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

    // add P
    (*theVector) -= *theLoad;

    return *theVector;
}



const Vector &
TrussSection::getResistingForceIncInertia()
{	
    this->getResistingForce();
    
    // now include the mass portion
    if (L != 0.0 && rho != 0.0) {
	
	const Vector &accel1 = theNodes[0]->getTrialAccel();
	const Vector &accel2 = theNodes[1]->getTrialAccel();	
	
	double M = 0.5*rho*L;
	int dof = dimension;
	int start = numDOF/2;
	for (int i=0; i<dof; i++) {
	    (*theVector)(i) += M*accel1(i);
	    (*theVector)(i+start) += M*accel2(i);
	}
    }    

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      *theVector += this->getRayleighDampingForces();

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
  data(3) = rho;
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

  static Vector data(6);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING TrussSection::recvSelf() - failed to receive Vector\n";
    return -1;
  }	      

  this->setTag((int)data(0));
  dimension = (int)data(1);
  numDOF = (int)data(2);
  rho = data(3);

  // truss now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr << "WARNING TrussSection::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally truss creates a new section object of the correct type,
  // sets its database tag and asks this new object to recveive itself.
  int sectClass = (int)data(4);
  int sectDb = (int)data(5);

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
TrussSection::displaySelf(Renderer &theViewer, int displayMode, float fact)
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
TrussSection::Print(OPS_Stream &s, int flag)
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
    
    double temp;
    int numDOF2 = numDOF/2;
    for (int i=0; i<dimension; i++) {
      temp = force*cosX[i];
      (*theVector)(i) = -force;
      (*theVector)(i+numDOF2) = force;
    }
    
    if (flag == 0) { // print everything
	s << "Element: " << this->getTag(); 
	s << " type: TrussSection  iNode: " << connectedExternalNodes(0);
	s << " jNode: " << connectedExternalNodes(1);
	s << " Mass density/length: " << rho;
	
	s << " \n\t strain: " << strain;
	s << " axial load: " << force;
	if (theVector != 0) 
	    s << " \n\t unbalanced load: " << *theVector;	
	s << " \t Section: " << *theSection;
	s << endln;
    } else if (flag == 1) {
	s << this->getTag() << "  " << strain << "  ";
	s << force << endln;
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
    for (int i=0; i<dimension; i++){
	dLength += (disp2(i)-disp1(i))*cosX[i];
    }

    // this method should never be called with L == 0
    return dLength/L;
}

Response*
TrussSection::setResponse(const char **argv, int argc, Information &eleInformation)
{
  //
  // we compare argv[0] for known response types for the Truss
  //

  // axial force
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 || 
      strcmp(argv[0],"axialForce") == 0) 
    return new ElementResponse(this, 1, 0);
  
  else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	   strcmp(argv[0],"deformation") == 0) 
    return new ElementResponse(this, 2, 0);
  
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
      
    default:
      if (responseID >= 100)
	  return theSection->getResponse(responseID-100, eleInformation);
      else
	  return -1;
  }
}

int
TrussSection::setParameter (const char **argv, int argc, Information &info)
{
    // a material parameter
    if (strcmp(argv[0],"section") == 0 || strcmp(argv[0],"-section") == 0) {
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
