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
                                                                        
// $Revision: 1.13 $
// $Date: 2010-02-06 19:08:26 $
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


#include <elementAPI.h>
#define OPS_Export 

OPS_Export void *
OPS_NewCorotTrussSectionElement()
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 4) {
    opserr << "Invalid Args want: element CorotTrussSection $tag $iNode $jNode $sectTag <-rho $rho> \n";
    return 0;	
  }

  int    iData[4];
  double rho = 0.0;
  int ndm = OPS_GetNDM();

  int numData = 4;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode, sectTag) in element CorotTrussSection " << endln;
    return 0;
  }

  SectionForceDeformation *theSection = OPS_GetSectionForceDeformation(iData[3]);
    
  if (theSection == 0) {
    opserr << "WARNING: Invalid section not found element CorotTrussSection " << iData[0] << " $iNode $jNode " << 
      iData[3] << " <-rho $rho> \n";
    return 0;
  }
  
  numRemainingArgs -= 4;
  while (numRemainingArgs > 1) {
    char argvS[10];
    if (OPS_GetString(argvS, 10) != 0) {
      opserr << "WARNING: Invalid optional string element CorotTrussSection " << iData[0] << 
	" $iNode $jNode $sectTag <-rho $rho>\n";
      return 0;
    } 
  
    if (strcmp(argvS,"-rho") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &rho) != 0) {
	opserr << "WARNING Invalid rho in element CorotTrussSection " << iData[0] << 
	  " $iNode $jNode $secTag <-rho $rho>\n";
	return 0;
      }
    } else {
      opserr << "WARNING: Invalid option " << argvS << "  in: element CorotTrussSection " << iData[0] << 
	" $iNode $jNode $secTag <-rho $rho>\n";
      return 0;
    }      
    numRemainingArgs -= 2;
  }

  //now create the ReinforcedConcretePlaneStress
  theElement = new CorotTrussSection(iData[0], ndm, iData[1], iData[2], *theSection, rho);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element CorotTrussSection " << iData[0] << 
      " $iNode $jNode $secTag <-rho $rho>\n";
  }

  return theElement;
}


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
    if (Lo == 0.0) { // - problem in setDomain() no further warnings
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

double
CorotTrussSection::computeCurrentStrain(void)
{
    // NOTE method will not be called if Lo == 0

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

    // this method should never be called with Lo == 0
    return (Ln-Lo)/Lo;
}

Response*
CorotTrussSection::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","Truss");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types for the CorotTruss
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

    } else if ((strcmp(argv[0],"axialForce") == 0) || (strcmp(argv[0],"basicForce") == 0) || 
        (strcmp(argv[0],"basicForces") == 0)) {
            output.tag("ResponseType", "N");
            theResponse =  new ElementResponse(this, 2, 0.0);

    } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformation") == 0 ||
        strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"basicDefo") == 0 ||
        strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0) {

            output.tag("ResponseType", "U");
            theResponse = new ElementResponse(this, 3, 0.0);

    // a section quantity
    }  else if (strcmp(argv[0],"section") ==0) {
        theResponse = theSection->setResponse(&argv[1], argc-1, output);

    }  

    output.endTag();
    return theResponse;
}

int 
CorotTrussSection::getResponse(int responseID, Information &eleInfo)
{
    double strain, force;

    switch (responseID) {
    case 1:
        return eleInfo.setVector(this->getResistingForce());

    case 2:
        if (Lo == 0.0) {
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
        return eleInfo.setDouble(force);    

    case 3:
        if (Lo == 0.0) {
            strain = 0.0;
        } else {
            strain = this->computeCurrentStrain();
        }
        return eleInfo.setDouble(Lo * strain);

    default:
        return -1;
    }
}
