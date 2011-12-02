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
                                                                        
// $Revision: 1.9 $
// $Date: 2003-02-25 23:33:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthSection.cpp,v $
                                                                        
// Written: MHS
// Created: Sept 2000
//
// Description: This file contains the implementation for the 
// ZeroLengthSection class.

#include <ZeroLengthSection.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <SectionForceDeformation.h>
#include <Renderer.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

Matrix ZeroLengthSection::K6(6,6);
Matrix ZeroLengthSection::K12(12,12);

Vector ZeroLengthSection::P6(6);
Vector ZeroLengthSection::P12(12);

//  Constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the ZeroLengthSection end nodes.

ZeroLengthSection::ZeroLengthSection(int tag, int dim, int Nd1, int Nd2, 
	       const Vector& x, const Vector& yprime, 
		   SectionForceDeformation& sec) : 
Element(tag, ELE_TAG_ZeroLengthSection),
connectedExternalNodes(2),
dimension(dim), numDOF(0), 
transformation(3,3), A(0), v(0), K(0), P(0),
theSection(0), order(0)
{
	// Obtain copy of section model
	theSection = sec.getCopy();
	
	if (theSection == 0) {
	  opserr << "ZeroLengthSection::ZeroLengthSection -- failed to get copy of section\n";
	  exit(-1);
	}

	// Get the section order
	order = theSection->getOrder();

	// Set up the transformation matrix of direction cosines
	this->setUp(Nd1, Nd2, x, yprime);
}

ZeroLengthSection::ZeroLengthSection() : 
Element(0, ELE_TAG_ZeroLengthSection),
connectedExternalNodes(2),
dimension(0), numDOF(0), 
transformation(3,3), A(0), v(0), K(0), P(0),
theSection(0), order(0)
{

}

ZeroLengthSection::~ZeroLengthSection()
{
    // invoke the destructor on any objects created by the object
    // that the object still holds a pointer to
    
	if (theSection != 0)
		delete theSection;
	if (A != 0)
		delete A;
	if (v != 0)
		delete v;
}

int
ZeroLengthSection::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ZeroLengthSection::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
ZeroLengthSection::getNodePtrs(void) 
{
  return theNodes;
}

int
ZeroLengthSection::getNumDOF(void) 
{
    return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the ZeroLengthSection element, we set matrix and vector pointers,
//    allocate space for t matrix and define it as the basic deformation-
//    displacement transformation matrix.
void
ZeroLengthSection::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	

    // if can't find both - send a warning message
    if (theNodes[0] == 0 || theNodes[1] == 0) {
      if (theNodes[0] == 0) 
	opserr << "ZeroLengthSection::setDomain() -- Nd2: " << Nd2 << " does not exist in ";
      else
	opserr << "ZeroLengthSection::setDomain() -- Nd2: " << Nd2 << " does not exist in ";
		
      opserr << "model for ZeroLengthSection with id " << this->getTag() << endln;
		
      return;
    }

    // now determine the number of dof and the dimension    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2) {
      opserr << "ZeroLengthSection::setDomain() -- nodes " << Nd1 << " and " << Nd2 << 
	"have differing dof at ends for ZeroLengthSection " << this->getTag() << endln;
      return;
    }	

    numDOF = 2*dofNd1;

    if (numDOF != 6 && numDOF != 12)
      opserr << "ZeroLengthSection::setDomain() -- element only works for 3 (2d) or 6 (3d) dof per node\n";
    
    // Set pointers to class wide objects
    if (numDOF == 6) {
      P = &P6;
      K = &K6;
    }
    else {
      P = &P12;
      K = &K12;
    }
    
    // Check that length is zero within tolerance
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector     diff = end1Crd - end2Crd;
    double L  = diff.Norm();
    double v1 = end1Crd.Norm();
    double v2 = end2Crd.Norm();
    double vm;
    
    vm = (v1<v2) ? v2 : v1;
    
    if (L > LENTOL*vm)
      opserr << "ZeroLengthSection::setDomain() -- Element " << this->getTag() << 
	"has L= " << L << ", which is greater than the tolerance\n";
	

// call the base class method
this->DomainComponent::setDomain(theDomain);

// Set up the A matrix
	this->setTransformation();
}   	 

int
ZeroLengthSection::commitState()
{
    int retVal=0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "ZeroLength::commitState () - failed in base class\n";
    }    
  // Commit the section
  retVal += theSection->commitState();
  return retVal;
}

int
ZeroLengthSection::revertToLastCommit()
{
  // Revert the section
  return theSection->revertToLastCommit();
}

int
ZeroLengthSection::revertToStart()
{
  // Revert the section to start
  return theSection->revertToStart();
}

const Matrix &
ZeroLengthSection::getTangentStiff(void)
{
	// Compute section deformation vector
	this->computeSectionDefs();

	// Set trial section deformation
	theSection->setTrialSectionDeformation(*v);

	// Get section tangent stiffness, the element basic stiffness
	const Matrix &kb = theSection->getSectionTangent();

	// Compute element stiffness ... K = A^*kb*A
	K->addMatrixTripleProduct(0.0, *A, kb, 1.0);

	return *K;
}

const Matrix &
ZeroLengthSection::getInitialStiff(void)
{
  // Get section tangent stiffness, the element basic stiffness
  const Matrix &kb = theSection->getInitialTangent();
  
  // Compute element stiffness ... K = A^*kb*A
  K->addMatrixTripleProduct(0.0, *A, kb, 1.0);
	
  return *K;
}

void 
ZeroLengthSection::zeroLoad(void)
{
	// does nothing now
}

int 
ZeroLengthSection::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "ZeroLengthSection::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  return -1;
}


int 
ZeroLengthSection::addInertiaLoadToUnbalance(const Vector &accel)
{
  // does nothing as element has no mass yet!
  return 0;
}

const Vector &
ZeroLengthSection::getResistingForce()
{
	// Compute section deformation vector
	this->computeSectionDefs();

	// Set trial section deformation
	theSection->setTrialSectionDeformation(*v);

	// Get section stress resultants, the element basic forces
	const Vector &q = theSection->getStressResultant();

	// Compute element resisting force ... P = A^*q
	P->addMatrixTransposeVector(0.0, *A, q, 1.0);

	return *P;
}


const Vector &
ZeroLengthSection::getResistingForceIncInertia()
{	
    this->getResistingForce();
    
    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      *P += this->getRayleighDampingForces();

    return *P;
}


int
ZeroLengthSection::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// ZeroLengthSection packs its data into an ID and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
	static ID idData(8);

	idData(0) = this->getTag();
	idData(1) = dimension;
	idData(2) = numDOF;
	idData(3) = order;
	idData(4) = connectedExternalNodes(0);
	idData(5) = connectedExternalNodes(1);
	idData(6) = theSection->getClassTag();
	
	int secDbTag = theSection->getDbTag();
	if (secDbTag == 0) {
		secDbTag = theChannel.getDbTag();
		if (secDbTag != 0)
			theSection->setDbTag(secDbTag);
	}
	idData(7) = secDbTag;

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
	  opserr << "ZeroLengthSection::sendSelf -- failed to send ID data\n";
			
		return res;
	}

	// Send the 3x3 direction cosine matrix, have to send it since it is only set
	// in the constructor and not setDomain()
	res += theChannel.sendMatrix(dataTag, commitTag, transformation);
	if (res < 0) {
	  opserr << "ZeroLengthSection::sendSelf -- failed to send transformation Matrix\n";
	  return res;
	}

	// Send the section
	res += theSection->sendSelf(commitTag, theChannel);
	if (res < 0) {
	  opserr << "ZeroLengthSection::sendSelf -- failed to send Section\n";
	  return res;
	}

	return res;
}

int
ZeroLengthSection::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
  
	int dataTag = this->getDbTag();

	// ZeroLengthSection creates an ID, receives the ID and then sets the 
	// internal data with the data in the ID

	static ID idData(8);

	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
	  opserr << "ZeroLengthSection::recvSelf -- failed to receive ID data\n";
	  return res;
	}

	res += theChannel.recvMatrix(dataTag, commitTag, transformation);
	if (res < 0) {
	  opserr << "ZeroLengthSection::recvSelf -- failed to receive transformation Matrix\n";
	  return res;
	}

	this->setTag(idData(0));
	dimension = idData(1);
	numDOF = idData(2);
	connectedExternalNodes(0) = idData(4);
	connectedExternalNodes(1) = idData(5);

	// Check that there is correct number of materials, reallocate if needed
	if (order != idData(3)) {

		order = idData(3);

		// Allocate transformation matrix
		if (A != 0)
			delete A;

		A = new Matrix(order, numDOF);

		if (A == 0) {
		  opserr << "ZeroLengthSection::recvSelf -- failed to allocate transformation Matrix\n";
		  exit(-1);
		}

		// Allocate section deformation vector
		if (v != 0)
			delete v;

		v = new Vector(order);

		if (v == 0) {
		  opserr << "ZeroLengthSection::recvSelf -- failed to allocate deformation Vector\n";
		  exit(-1);
		}

		if (numDOF == 6) {
			P = &P6;
			K = &K6;
		}
		else {
			P = &P12;
			K = &K12;
		}
	}

	int secClassTag = idData(6);

	// If null, get a new one from the broker
	if (theSection == 0)
		theSection = theBroker.getNewSection(secClassTag);

	// If wrong type, get a new one from the broker
	if (theSection->getClassTag() != secClassTag) {
		delete theSection;
		theSection = theBroker.getNewSection(secClassTag);
	}

	// Check if either allocation failed from broker
	if (theSection == 0) {
	  opserr << "ZeroLengthSection::recvSelf -- failed to allocate new Section\n";
	  return -1;
	}

	// Receive the section
	theSection->setDbTag(idData(7));
	res += theSection->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
	  opserr << "ZeroLengthSection::recvSelf -- failed to receive Section\n";
	  return res;
	}

	return res;
}

int
ZeroLengthSection::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // ensure setDomain() worked
    if (theNodes[0] == 0 || theNodes[1] == 0)
		return 0;

    // first determine the two end points of the ZeroLengthSection based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();    

    if (displayMode == 1 || displayMode == 2) {
		static Vector v1(3);
		static Vector v2(3);
		for (int i = 0; i < dimension; i++) {
			v1(i) = end1Crd(i)+end1Disp(i)*fact;
			v2(i) = end2Crd(i)+end2Disp(i)*fact;    
		}
		
		return theViewer.drawLine(v1, v2, 0.0, 0.0);
    }

    return 0;
}

void
ZeroLengthSection::Print(OPS_Stream &s, int flag)
{
	s << "ZeroLengthSection, tag: " << this->getTag() << endln;
	s << "\tConnected Nodes: " << connectedExternalNodes << endln;
	s << "\tSection, tag: " << theSection->getTag() << endln;
}

Response*
ZeroLengthSection::setResponse(const char **argv, int argc, Information &eleInformation)
{
	// element forces
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
		return new ElementResponse(this, 1, *P);

    // element stiffness matrix
    else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
		return new ElementResponse(this, 2, *K);

    else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
		strcmp(argv[0],"deformation") == 0) {
		return new ElementResponse(this, 3, Vector(order));
	}     

	else if (strcmp(argv[0],"section") == 0)
		return theSection->setResponse(&argv[1], argc-1, eleInformation);

	else 
		return 0;
}

int 
ZeroLengthSection::getResponse(int responseID, Information &eleInfo)
{
  switch (responseID) {
	case 1:
		return eleInfo.setVector(this->getResistingForce());
		
    case 2:
		return eleInfo.setMatrix(this->getTangentStiff());

	case 3:
		this->computeSectionDefs();
		return eleInfo.setVector(*v);

    default:
		return -1;
  }
}

// Private methods


// Establish the external nodes and set up the transformation matrix
// for orientation
void
ZeroLengthSection::setUp(int Nd1, int Nd2, const Vector &x, const Vector &yp)
{ 
    // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2) {
    opserr << "ZeroLengthSection::setUp -- failed to create an ID of correct size\n";
    exit(-1);
  }
    
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
	int i;
    for (i=0; i<2; i++)
      theNodes[i] = 0;

    // check that vectors for orientation are correct size
    if ( x.Size() != 3 || yp.Size() != 3 )
      opserr << "ZeroLengthSection::setUp -- incorrect dimension of orientation vectors\n";
			

    // establish orientation of element for the tranformation matrix
    // z = x cross yp
    static Vector z(3);
    z(0) = x(1)*yp(2) - x(2)*yp(1);
    z(1) = x(2)*yp(0) - x(0)*yp(2);
    z(2) = x(0)*yp(1) - x(1)*yp(0);

    // y = z cross x
    static Vector y(3);
    y(0) = z(1)*x(2) - z(2)*x(1);
    y(1) = z(2)*x(0) - z(0)*x(2);
    y(2) = z(0)*x(1) - z(1)*x(0);

    // compute length(norm) of vectors
    double xn = x.Norm();
    double yn = y.Norm();
    double zn = z.Norm();

    // check valid x and y vectors, i.e. not parallel and of zero length
    if (xn == 0 || yn == 0 || zn == 0)
      opserr << "ZeroLengthSection::setUp -- invalid vectors to constructor\n";
    
    // create transformation matrix of direction cosines
    for (i = 0; i < 3; i++) {
		transformation(0,i) = x(i)/xn;
		transformation(1,i) = y(i)/yn;
		transformation(2,i) = z(i)/zn;
	}
}

// Set basic deformation-displacement transformation matrix for section
void 
ZeroLengthSection::setTransformation(void)
{
	// Allocate transformation matrix
	if (A != 0)
		delete A;

	A = new Matrix(order, numDOF);

	if (A == 0)
	  opserr << "ZeroLengthSection::setTransformation -- failed to allocate transformation Matrix\n";
			

	// Allocate section deformation vector
	if (v != 0)
		delete v;

	v = new Vector(order);

	// Get the section code
	const ID &code = theSection->getType();
		
	// Set a reference to make the syntax nicer
	Matrix &tran = *A;
	
	tran.Zero();

	// Loop over the section code
	for (int i = 0; i < order; i++) {

		// Fill in row i of A based on section code
		switch(code(i)) {

		// The in-plane transformations
		case SECTION_RESPONSE_MZ:
			if (numDOF == 6) {
				tran(i,3) = 0.0;
				tran(i,4) = 0.0;
				tran(i,5) = transformation(2,2);
			}
			else if (numDOF == 12) {
				tran(i,9) = transformation(2,0);
				tran(i,10) = transformation(2,1);
				tran(i,11) = transformation(2,2);
			}
			break;
		case SECTION_RESPONSE_P:
			if (numDOF == 6) {
				tran(i,3) = transformation(0,0);
				tran(i,4) = transformation(0,1);
				tran(i,5) = 0.0;
			}
			else if (numDOF == 12) {
				tran(i,6) = transformation(0,0);
				tran(i,7) = transformation(0,1);
				tran(i,8) = transformation(0,2);
			}
			break;
		case SECTION_RESPONSE_VY:
			if (numDOF == 6) {
				tran(i,3) = transformation(1,0);
				tran(i,4) = transformation(1,1);
				tran(i,5) = 0.0;
			}
			else if (numDOF == 12) {
				tran(i,6) = transformation(1,0);
				tran(i,7) = transformation(1,1);
				tran(i,8) = transformation(1,2);
			}
			break;

		// The out-of-plane transformations
		case SECTION_RESPONSE_MY:
			if (numDOF == 12) {
				tran(i,9) = transformation(1,0);
				tran(i,10) = transformation(1,1);
				tran(i,11) = transformation(1,2);
			}
			break;
		case SECTION_RESPONSE_VZ:
			if (numDOF == 12) {
				tran(i,6) = transformation(2,0);
				tran(i,7) = transformation(2,1);
				tran(i,8) = transformation(2,2);
			}
			break;
		case SECTION_RESPONSE_T:
			if (numDOF == 12) {
				tran(i,9) = transformation(0,0);
				tran(i,10) = transformation(0,1);
				tran(i,11) = transformation(0,2);
			}
			break;
		default:
			break;
		}

		// Fill in first half of transformation matrix with negative sign
		for (int j = 0; j < numDOF/2; j++ )
			tran(i,j) = -tran(i,j+numDOF/2);
	}
}
		     
void
ZeroLengthSection::computeSectionDefs(void)
{
	// Get nodal displacements
	const Vector &u1 = theNodes[0]->getTrialDisp();
	const Vector &u2 = theNodes[1]->getTrialDisp();

	// Compute differential displacements
	const Vector diff = u2 - u1;

	// Set some references to make the syntax nicer
	Vector &def = *v;
	const Matrix &tran = *A;

	def.Zero();

	// Compute element basic deformations ... v = A*(u2-u1)
	for (int i = 0; i < order; i++)
		for (int j = 0; j < numDOF/2; j++)
			def(i) += -diff(j)*tran(i,j);
}
