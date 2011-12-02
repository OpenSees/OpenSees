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
                                                                        
// $Revision: 1.7 $
// $Date: 2003-03-12 19:20:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/truss/CorotTrussSection.cpp,v $
                                                                        
// Written: MHS 
// Created: May 2001
//
// Description: This file contains the class implementation for CorotTrussSection.

#include <CorotTrussSection.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <SectionForceDeformation.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

Matrix CorotTrussSection::M2(2,2);
Matrix CorotTrussSection::M4(4,4);
Matrix CorotTrussSection::M6(6,6);
Matrix CorotTrussSection::M12(12,12);

Vector CorotTrussSection::V2(2);
Vector CorotTrussSection::V4(4);
Vector CorotTrussSection::V6(6);
Vector CorotTrussSection::V12(12);

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the CorotTrussSection end nodes.
CorotTrussSection::CorotTrussSection(int tag, int dim,
			   int Nd1, int Nd2, 
			   SectionForceDeformation &theSec,
			   double r)
  :Element(tag,ELE_TAG_CorotTrussSection),     
  theSection(0), connectedExternalNodes(2),
  numDOF(0), numDIM(dim),
  Lo(0.0), Ln(0.0), rho(r), R(3,3),
  theMatrix(0), theVector(0)
{
  // get a copy of the material and check we obtained a valid copy
  theSection = theSec.getCopy();
  if (theSection == 0) {
    opserr << "FATAL CorotTrussSection::CorotTrussSection - " << tag << 
      "failed to get a copy of material with tag " << theSec.getTag() << endln;
    exit(-1);
  }
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2) {
    opserr << "FATAL CorotTrussSection::CorotTrussSection - " <<  tag <<
      "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;        

  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
CorotTrussSection::CorotTrussSection()
  :Element(0,ELE_TAG_CorotTrussSection),     
  theSection(0),connectedExternalNodes(2),
  numDOF(0), numDIM(0),
  Lo(0.0), Ln(0.0), rho(0.0), R(3,3),
  theMatrix(0), theVector(0)
{
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2) {
    opserr << "FATAL CorotTrussSection::CorotTrussSection - failed to create an ID of size 2\n";
    exit(-1);
  }			  
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
CorotTrussSection::~CorotTrussSection()
{
  // invoke the destructor on any objects created by the object
  // that the object still holds a pointer to
  if (theSection != 0)
    delete theSection;
}

int
CorotTrussSection::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
CorotTrussSection::getExternalNodes(void) 
{
	return connectedExternalNodes;
}


Node **
CorotTrussSection::getNodePtrs(void) 
{
  return theNodes;
}

int
CorotTrussSection::getNumDOF(void) 
{
	return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the CorotTrussSection element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
CorotTrussSection::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    Lo = 0.0;
    Ln = 0.0;
    return;
  }
  
  // first set the node pointers
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);	
  
  // if can't find both - send a warning message
  if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
    opserr << "CorotTrussSection::setDomain() - CorotTrussSection " << 			  
      this->getTag() << " node doe not exist in the model\n";
    
    // fill this in so don't segment fault later
    numDOF = 6;    
    
    return;
  }
  
  // now determine the number of dof and the dimesnion    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != dofNd2) {
    opserr << "WARNING CorotTrussSection::setDomain(): nodes have differing dof at ends for CorotTrussSection" <<
      this->getTag() << endln;
    
    // fill this in so don't segment fault later
    numDOF = 6;    
    
    return;
  }	
  
  if (numDIM == 1 && dofNd1 == 1) {
    numDOF = 2;
    theMatrix = &M2;
    theVector = &V2;
  }
  else if (numDIM == 2 && dofNd1 == 2) {
    numDOF = 4;
    theMatrix = &M4;
    theVector = &V4;
  }
  else if (numDIM == 2 && dofNd1 == 3) {
    numDOF = 6;
    theMatrix = &M6;
    theVector = &V6;
  }
  else if (numDIM == 3 && dofNd1 == 3) {
    numDOF = 6;
    theMatrix = &M6;
    theVector = &V6;
  }
  else if (numDIM == 3 && dofNd1 == 6) {
    numDOF = 12;
    theMatrix = &M12;
    theVector = &V12;
  }
  else {
    opserr << "CorotTrussSection::setDomain -- nodal DOF not compatible with element " << this->getTag() << endln;
    
    // fill this in so don't segment fault later
    numDOF = 6;    
    
    return;
  }

	// call the base class method
	this->DomainComponent::setDomain(theDomain);

	// now determine the length, cosines and fill in the transformation
	// NOTE t = -t(every one else uses for residual calc)
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();

	// Determine global offsets
    double cosX[3];
    cosX[0] = 0.0;  cosX[1] = 0.0;  cosX[2] = 0.0;
    int i;
    for (i = 0; i < numDIM; i++) {
        cosX[i] += end2Crd(i)-end1Crd(i);
    }

	// Set undeformed and initial length
	Lo = cosX[0]*cosX[0] + cosX[1]*cosX[1] + cosX[2]*cosX[2];
	Lo = sqrt(Lo);
	Ln = Lo;

    // Initial offsets
   	d21[0] = Lo;
	d21[1] = 0.0;
	d21[2] = 0.0;

	// Set global orientation
	cosX[0] /= Lo;
	cosX[1] /= Lo;
	cosX[2] /= Lo;

	R(0,0) = cosX[0];
	R(0,1) = cosX[1];
	R(0,2) = cosX[2];

	// Element lies outside the YZ plane
	if (fabs(cosX[0]) > 0.0) {
		R(1,0) = -cosX[1];
		R(1,1) =  cosX[0];
		R(1,2) =  0.0;

		R(2,0) = -cosX[0]*cosX[2];
		R(2,1) = -cosX[1]*cosX[2];
		R(2,2) =  cosX[0]*cosX[0] + cosX[1]*cosX[1];
	}
	// Element is in the YZ plane
	else {
		R(1,0) =  0.0;
		R(1,1) = -cosX[2];
		R(1,2) =  cosX[1];

		R(2,0) =  1.0;
		R(2,1) =  0.0;
		R(2,2) =  0.0;
	}

	// Orthonormalize last two rows of R
	double norm;
	for (i = 1; i < 3; i++) {
		norm = sqrt(R(i,0)*R(i,0) + R(i,1)*R(i,1) + R(i,2)*R(i,2));
		R(i,0) /= norm;
		R(i,1) /= norm;
		R(i,2) /= norm;
	}
}

int
CorotTrussSection::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "CorotTrussSection::commitState () - failed in base class";
  }    
  retVal = theSection->commitState();
  return retVal;
}

int
CorotTrussSection::revertToLastCommit()
{
	// Revert the material
	return theSection->revertToLastCommit();
}

int
CorotTrussSection::revertToStart()
{
	// Revert the material to start
	return theSection->revertToStart();
}

int
CorotTrussSection::update(void)
{
	// Nodal displacements
	const Vector &end1Disp = theNodes[0]->getTrialDisp();
	const Vector &end2Disp = theNodes[1]->getTrialDisp();    

    // Initial offsets
	d21[0] = Lo;
	d21[1] = 0.0;
	d21[2] = 0.0;

   	// Update offsets in basic system due to nodal displacements
	int i;
    for (i = 0; i < numDIM; i++) {
        d21[0] += R(0,i)*(end2Disp(i)-end1Disp(i));
        d21[1] += R(1,i)*(end2Disp(i)-end1Disp(i));
        d21[2] += R(2,i)*(end2Disp(i)-end1Disp(i));
    }

	// Compute new length
	Ln = d21[0]*d21[0] + d21[1]*d21[1] + d21[2]*d21[2];
	Ln = sqrt(Ln);

	// Compute engineering strain
	double strain = (Ln-Lo)/Lo;

	int order = theSection->getOrder();
	const ID &code = theSection->getType();

	static double data[10];
	Vector e(data, order);
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			e(i) = strain;
		else
			e(i) = 0.0;
	}

	// Set material trial strain
	return theSection->setTrialSectionDeformation(e);
}

const Matrix &
CorotTrussSection::getTangentStiff(void)
{
    static Matrix kl(3,3);

    // Material stiffness
    //
    // Get material tangent
	int order = theSection->getOrder();
	const ID &code = theSection->getType();

	const Matrix &ks = theSection->getSectionTangent();
	const Vector &s  = theSection->getStressResultant();
	
	double EA = 0.0;
	double q = 0.0;

	int i,j;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P) {
			EA += ks(i,i);
			q  += s(i);
		}
	}

	EA /= (Ln*Ln*Lo);

    for (i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            kl(i,j) = EA*d21[i]*d21[j];

	// Geometric stiffness
	//
	// Get material stress
	double SA = q/(Ln*Ln*Ln);
	double SL = q/Ln;

    for (i = 0; i < 3; i++) {
        kl(i,i) += SL;
        for (j = 0; j < 3; j++)
            kl(i,j) -= SA*d21[i]*d21[j];
    }
    
    // Compute R'*kl*R
    static Matrix kg(3,3);
    kg.addMatrixTripleProduct(0.0, R, kl, 1.0);

    Matrix &K = *theMatrix;
    K.Zero();

    // Copy stiffness into appropriate blocks in element stiffness
    int numDOF2 = numDOF/2;
    for (i = 0; i < numDIM; i++) {
        for (j = 0; j < numDIM; j++) {
            K(i,j)                 =  kg(i,j);
            K(i,j+numDOF2)         = -kg(i,j);
            K(i+numDOF2,j)         = -kg(i,j);
            K(i+numDOF2,j+numDOF2) =  kg(i,j);
        }
    }

    return *theMatrix;
}

const Matrix &
CorotTrussSection::getInitialStiff(void)
{
    static Matrix kl(3,3);

    // Material stiffness
    //
    // Get material tangent
    int order = theSection->getOrder();
    const ID &code = theSection->getType();
    
    const Matrix &ks = theSection->getInitialTangent();
    
    double EA = 0.0;

    int i,j;
    for (i = 0; i < order; i++) {
      if (code(i) == SECTION_RESPONSE_P) {
	EA += ks(i,i);
      }
    }

    kl(0,0) = EA / Lo;

    // Compute R'*kl*R
    static Matrix kg(3,3);
    kg.addMatrixTripleProduct(0.0, R, kl, 1.0);
    
    Matrix &K = *theMatrix;
    K.Zero();
    
    // Copy stiffness into appropriate blocks in element stiffness
    int numDOF2 = numDOF/2;
    for (i = 0; i < numDIM; i++) {
      for (j = 0; j < numDIM; j++) {
	K(i,j)                 =  kg(i,j);
	K(i,j+numDOF2)         = -kg(i,j);
	K(i+numDOF2,j)         = -kg(i,j);
	K(i+numDOF2,j+numDOF2) =  kg(i,j);
      }
    }

    return *theMatrix;
}

const Matrix &
CorotTrussSection::getMass(void)
{
    Matrix &Mass = *theMatrix;
    Mass.Zero();

    // check for quick return
    if (Lo == 0.0 || rho == 0.0)
	return Mass;

    double M = 0.5*rho*Lo;
    int numDOF2 = numDOF/2;
    for (int i = 0; i < numDIM; i++) {
        Mass(i,i)                 = M;
        Mass(i+numDOF2,i+numDOF2) = M;
    }

    return *theMatrix;
}

void 
CorotTrussSection::zeroLoad(void)
{
	return;
}

int 
CorotTrussSection::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "CorotTrussSection::addLoad - load type unknown for truss with tag: " <<  this->getTag() << endln;
  return -1;
}

int 
CorotTrussSection::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

const Vector &
CorotTrussSection::getResistingForce()
{
	int order = theSection->getOrder();
	const ID &code = theSection->getType();

	const Vector &s  = theSection->getStressResultant();

	double SA = 0.0;
	
	int i;
	for (i = 0; i < order; i++) {
		if (code(i) == SECTION_RESPONSE_P)
			SA += s(i);
	}

	SA /= Ln;

    static Vector ql(3);

	ql(0) = d21[0]*SA;
	ql(1) = d21[1]*SA;
	ql(2) = d21[2]*SA;

    static Vector qg(3);
    qg.addMatrixTransposeVector(0.0, R, ql, 1.0);

    Vector &P = *theVector;
    P.Zero();

    // Copy forces into appropriate places
    int numDOF2 = numDOF/2;
    for (i = 0; i < numDIM; i++) {
        P(i)         = -qg(i);
        P(i+numDOF2) =  qg(i);
    }

    return *theVector;
}

const Vector &
CorotTrussSection::getResistingForceIncInertia()
{	
    Vector &P = *theVector;
    P = this->getResistingForce();
    
    if (rho != 0.0) {
	
      const Vector &accel1 = theNodes[0]->getTrialAccel();
      const Vector &accel2 = theNodes[1]->getTrialAccel();	
      
      double M = 0.5*rho*Lo;
      int numDOF2 = numDOF/2;
      for (int i = 0; i < numDIM; i++) {
	P(i)         += M*accel1(i);
	P(i+numDOF2) += M*accel2(i);
      }
    }

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      *theVector += this->getRayleighDampingForces();

    return *theVector;
}

int
CorotTrussSection::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int
CorotTrussSection::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	return -1;
}

int
CorotTrussSection::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	// ensure setDomain() worked
	if (Ln == 0.0)
		return 0;

	// first determine the two end points of the CorotTrussSection based on
	// the display factor (a measure of the distorted image)
	// store this information in 2 3d vectors v1 and v2
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();	
	const Vector &end1Disp = theNodes[0]->getDisp();
	const Vector &end2Disp = theNodes[1]->getDisp();    

	static Vector v1(3);
	static Vector v2(3);
	for (int i = 0; i < numDIM; i++) {
		v1(i) = end1Crd(i)+end1Disp(i)*fact;
		v2(i) = end2Crd(i)+end2Disp(i)*fact;    
	}

	return theViewer.drawLine(v1, v2, 1.0, 1.0);
}

void
CorotTrussSection::Print(OPS_Stream &s, int flag)
{
	s << "\nCorotTrussSection, tag: " << this->getTag() << endln;
	s << "\tConnected Nodes: " << connectedExternalNodes;
	s << "\tUndeformed Length: " << Lo << endln;
	s << "\tCurrent Length: " << Ln << endln;
	s << "\tMass Density/Length: " << rho << endln;
	s << "\tRotation matrix: " << endln;

	if (theSection) {
		s << "\tSection, tag: " << theSection->getTag() << endln;
		theSection->Print(s,flag);
	}
}

Response*
CorotTrussSection::setResponse(const char **argv, int argc, Information &eleInfo)
{
    // a material quantity    
	if (strcmp(argv[0],"section") == 0)
		return theSection->setResponse(&argv[1], argc-1, eleInfo);
    
	else
		return 0;
}

int 
CorotTrussSection::getResponse(int responseID, Information &eleInfo)
{
	return 0;
}
