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
// $Date: 2010-02-04 01:17:47 $
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
#include <elementAPI.h>

Matrix ZeroLengthSection::K6(6,6);
Matrix ZeroLengthSection::K12(12,12);

Vector ZeroLengthSection::P6(6);
Vector ZeroLengthSection::P12(12);

void * OPS_ADD_RUNTIME_VPV(OPS_ZeroLengthSection)
{
    int ndm = OPS_GetNDM();
    
    if(OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"insufficient arguments for ZeroLengthSection\n";
	return 0;
    }

    // get eleTag,iNode,jNode,secTag
    int iData[4];
    int numData = 4;
    if(OPS_GetIntInput(&numData,&iData[0]) < 0) {
	opserr<<"WARNING: invalid integer inputs\n";
	return 0;
    }

    // options
    Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
    Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
    double *x_ptr=&x(0), *y_ptr=&y(0);
    int doRayleighDamping = 1;
    while(OPS_GetNumRemainingInputArgs() > 1) {
	const char* type = OPS_GetString();
	if(strcmp(type, "-orient") == 0) {
	    if(OPS_GetNumRemainingInputArgs() > 5) {
		numData = 3;
		if(OPS_GetDoubleInput(&numData,x_ptr) < 0) {
		    opserr<<"WARNING: invalid double inputs\n";
		    return 0;
		}
		if(OPS_GetDoubleInput(&numData,y_ptr) < 0) {
		    opserr<<"WARNING: invalid double inputs\n";
		    return 0;
		}
	    }
	} else if(strcmp(type, "-doRayleigh") == 0) {
	    numData = 1;
	    if(OPS_GetIntInput(&numData,&doRayleighDamping) < 0) {
		opserr<<"WARNING: invalid integer inputs\n";
		return 0;
	    }
	}
    }

    // get section
    SectionForceDeformation* theSection = OPS_getSectionForceDeformation(iData[3]);
    if(theSection == 0) {
	opserr << "zeroLengthSection -- no section with tag " << iData[0] << " exists in Domain\n";
	return 0;
    }

    return new ZeroLengthSection(iData[0],ndm,iData[1],iData[2],x,y,*theSection,doRayleighDamping);
}

//  Constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the ZeroLengthSection end nodes.

ZeroLengthSection::ZeroLengthSection(int tag, int dim, int Nd1, int Nd2, 
				     const Vector& x, const Vector& yprime, 
				     SectionForceDeformation& sec,
				     int doRayleigh) 
 :Element(tag, ELE_TAG_ZeroLengthSection),
  connectedExternalNodes(2),
  dimension(dim), numDOF(0), 
  transformation(3,3), useRayleighDamping(doRayleigh),A(0), v(0), K(0), P(0),
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
    Vector diff = end1Crd - end2Crd;
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
ZeroLengthSection::update()		// MSN: added to allow error identification in setTrialSectionDeformation()
{
	// Compute section deformation vector
	this->computeSectionDefs();

	// Set trial section deformation
	if (theSection->setTrialSectionDeformation(*v) < 0) {
		opserr << "WARNING! ZeroLengthSection::update() - element: " << this->getTag() << " failed in setTrialSectionDeformation\n";
		return -1;
	}

	return 0;
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
	// this->computeSectionDefs();	// MSN: commented out beause the method "update()" was added to the class

	// Set trial section deformation
	// theSection->setTrialSectionDeformation(*v);	// MSN: commented out beause the method "update()" was added to the class

	// Get section tangent stiffness, the element basic stiffness
	const Matrix &kb = theSection->getSectionTangent();

	// Compute element stiffness ... K = A^*kb*A
	K->addMatrixTripleProduct(0.0, *A, kb, 1.0);

	return *K;
}

const Matrix &
ZeroLengthSection::getDamp()
{	
  if (useRayleighDamping == 1)
    return this->Element::getDamp();

  K->Zero();
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
	// this->computeSectionDefs();	// MSN: commented out beause the method "update()" was added to the class

	// Set trial section deformation
	// theSection->setTrialSectionDeformation(*v);	// MSN: commented out beause the method "update()" was added to the class

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

    if (useRayleighDamping == 1)
      if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
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
	static ID idData(9);

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
	idData(8) = useRayleighDamping;


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

	static ID idData(9);

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
	useRayleighDamping =idData(8);

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
ZeroLengthSection::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // ensure setDomain() worked
    if (theNodes[0] == 0 || theNodes[1] == 0)
		return 0;

	// get the end point display coords    
	static Vector v1(3);
	static Vector v2(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);

	// draw the line
	return theViewer.drawLine(v1, v2, 0.0, 0.0, this->getTag());
}

void
ZeroLengthSection::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "ZeroLengthSection, tag: " << this->getTag() << endln;
        s << "\tConnected Nodes: " << connectedExternalNodes << endln;
        s << "\tSection, tag: " << theSection->getTag() << endln;
        theSection->Print(s, flag);
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ZeroLengthSection\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"section\": \"" << theSection->getTag() << "\", ";
        s << "\"transMatrix\": [[";
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                if (j < 2)
                    s << transformation(i, j) << ", ";
                else if (j == 2 && i < 2)
                    s << transformation(i, j) << "], [";
                else if (j == 2 && i == 2)
                    s << transformation(i, j) << "]]}";
            }
        }
    }
}

Response*
ZeroLengthSection::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","ZeroLengthSection");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    char outputData[5];
    // element forces
    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0)
        || (strcmp(argv[0],"globalForces") == 0) || (strcmp(argv[0],"globalforces") == 0)) {

            for (int i=0; i<P->Size(); i++) {
                sprintf(outputData,"P%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 1, *P);

    } else if ((strcmp(argv[0],"basicForce") == 0) || (strcmp(argv[0],"basicForces") == 0) ||
	       (strcmp(argv[0],"localForce") == 0) || (strcmp(argv[0],"localForces") == 0)) {

        for (int i=0; i<order; i++) {
            sprintf(outputData,"P%d",i+1);
            output.tag("ResponseType",outputData);
        }
        theResponse = new ElementResponse(this, 2, Vector(order));

    } else if (strcmp(argv[0],"basicStiffness") == 0) {

      theResponse = new ElementResponse(this, 13, Matrix(order,order));


    } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
        strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"basicDeformation") == 0) {

            for (int i=0; i<order; i++) {
                sprintf(outputData,"e%d",i+1);
                output.tag("ResponseType",outputData);
            }
            theResponse = new ElementResponse(this, 3, Vector(order));

    // a section quantity
    } else if (strcmp(argv[0],"section") == 0) {
        theResponse = theSection->setResponse(&argv[1], argc-1, output);
    }

    output.endTag();
    return theResponse;

}

int 
ZeroLengthSection::getResponse(int responseID, Information &eleInfo)
{
    Vector q(order);
    Matrix kb(order,order);

    switch (responseID) {
    case 1:
        return eleInfo.setVector(this->getResistingForce());

    case 2:
        // Get section stress resultants, the element basic forces
        q = theSection->getStressResultant();
        return eleInfo.setVector(q);

    case 3:
        this->computeSectionDefs();
        return eleInfo.setVector(*v);

    case 13:
      kb = theSection->getSectionTangent();
      return eleInfo.setMatrix(kb);


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
			

    // establish orientation of element for the transformation matrix
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
	Vector diff = u2 - u1;

	// Set some references to make the syntax nicer
	Vector &def = *v;
	const Matrix &tran = *A;

	def.Zero();

	// Compute element basic deformations ... v = A*(u2-u1)
	for (int i = 0; i < order; i++)
		for (int j = 0; j < numDOF/2; j++)
			def(i) += -diff(j)*tran(i,j);
}

int
ZeroLengthSection::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  return theSection->setParameter(argv, argc, param);
}

const Vector &
ZeroLengthSection::getResistingForceSensitivity(int gradIndex)
{
  // Compute section deformation vector
  // this->computeSectionDefs();	// MSN: commented out beause the method "update()" was added to the class

  // Set trial section deformation
  // theSection->setTrialSectionDeformation(*v);	// MSN: commented out beause the method "update()" was added to the class

  // Get section stress resultants, the element basic forces
  const Vector &dqdh = theSection->getStressResultantSensitivity(gradIndex, true);

  // Compute element resisting force ... P = A^*q
  P->addMatrixTransposeVector(0.0, *A, dqdh, 1.0);

  return *P;
}
    
int
ZeroLengthSection::commitSensitivity(int gradIndex, int numGrads)
{
  // Get nodal displacement sensitivity
  Vector diff(numDOF/2);
  for (int i = 0; i < numDOF/2; i++) {
    diff(i) = theNodes[1]->getDispSensitivity(i+1,gradIndex) - theNodes[0]->getDispSensitivity(i+1,gradIndex);
  }

  // Set some references to make the syntax nicer
  Vector &dedh = *v;
  const Matrix &tran = *A;

  dedh.Zero();

  // Compute element basic deformations ... v = A*(u2-u1)
  for (int i = 0; i < order; i++)
    for (int j = 0; j < numDOF/2; j++)
      dedh(i) += -diff(j)*tran(i,j);

  return theSection->commitSensitivity(dedh, gradIndex, numGrads);
}
