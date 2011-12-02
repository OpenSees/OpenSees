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
                                                                        
// $Revision: 1.12 $
// $Date: 2006-08-04 19:13:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/truss/CorotTruss.cpp,v $
                                                                        
// Written: MHS 
// Created: May 2001
//
// Description: This file contains the class implementation for CorotTruss.

#include <CorotTruss.h>
#include <Information.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

Matrix CorotTruss::M2(2,2);
Matrix CorotTruss::M4(4,4);
Matrix CorotTruss::M6(6,6);
Matrix CorotTruss::M12(12,12);

Vector CorotTruss::V2(2);
Vector CorotTruss::V4(4);
Vector CorotTruss::V6(6);
Vector CorotTruss::V12(12);

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the CorotTruss end nodes.
CorotTruss::CorotTruss(int tag, int dim,
			   int Nd1, int Nd2, 
			   UniaxialMaterial &theMat,
			   double a, double r)
  :Element(tag,ELE_TAG_CorotTruss),     
  theMaterial(0), connectedExternalNodes(2),
  numDOF(0), numDIM(dim),
  Lo(0.0), Ln(0.0), 
  A(a), rho(r), R(3,3),
  theMatrix(0), theVector(0)
{
  // get a copy of the material and check we obtained a valid copy
  theMaterial = theMat.getCopy();
  if (theMaterial == 0) {
    opserr << "FATAL CorotTruss::CorotTruss - " <<  tag <<
      "failed to get a copy of material with tag " << theMat.getTag() << endln;
    exit(-1);
  }
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2) {
    opserr << "FATAL CorotTruss::CorotTruss - " <<  tag <<
      "failed to create an ID of size 2\n";
    exit(-1);
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;        

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
CorotTruss::CorotTruss()
  :Element(0,ELE_TAG_CorotTruss),     
  theMaterial(0),connectedExternalNodes(2),
  numDOF(0), numDIM(0),
  Lo(0.0), Ln(0.0), 
  A(0.0), rho(0.0), R(3,3),
  theMatrix(0), theVector(0)
{
  // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2) {
    opserr << "FATAL CorotTruss::CorotTruss - failed to create an ID of size 2\n";
    exit(-1);
  }

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
CorotTruss::~CorotTruss()
{
  // invoke the destructor on any objects created by the object
  // that the object still holds a pointer to
  if (theMaterial != 0)
    delete theMaterial;
}

int
CorotTruss::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
CorotTruss::getExternalNodes(void) 
{
	return connectedExternalNodes;
}

Node **
CorotTruss::getNodePtrs(void) 
{
  return theNodes;
}

int
CorotTruss::getNumDOF(void) 
{
	return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the CorotTruss element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
CorotTruss::setDomain(Domain *theDomain)
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
    opserr << "CorotTruss::setDomain() - CorotTruss " << this->getTag() << " node " <<
      Nd1 << "does not exist in the model \n";
    
    // fill this in so don't segment fault later
    numDOF = 6;    
    
    return;
  }
  
  // now determine the number of dof and the dimesnion    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != dofNd2) {
    opserr << "WARNING CorotTruss::setDomain(): nodes " << Nd1 <<
      " and " << Nd2 << "have differing dof at ends for CorotTruss " << this->getTag() << endln;
    
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
    opserr << " CorotTruss::setDomain -- nodal DOF " << dofNd1 << " not compatible with element\n";
    
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
CorotTruss::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "CorotTruss::commitState () - failed in base class\n";
  }    
  retVal = theMaterial->commitState();
  return retVal;
}

int
CorotTruss::revertToLastCommit()
{
	// Revert the material
	return theMaterial->revertToLastCommit();
}

int
CorotTruss::revertToStart()
{
	// Revert the material to start
	return theMaterial->revertToStart();
}

int
CorotTruss::update(void)
{
  // Nodal displacements
  const Vector &end1Disp = theNodes[0]->getTrialDisp();
  const Vector &end2Disp = theNodes[1]->getTrialDisp();    
  
  // Initial offsets
  d21[0] = Lo;
  d21[1] = 0.0;
  d21[2] = 0.0;
  
  // Update offsets in basic system due to nodal displacements
  for (int i = 0; i < numDIM; i++) {
    d21[0] += R(0,i)*(end2Disp(i)-end1Disp(i));
    d21[1] += R(1,i)*(end2Disp(i)-end1Disp(i));
    d21[2] += R(2,i)*(end2Disp(i)-end1Disp(i));
  }
  
  // Compute new length
  Ln = d21[0]*d21[0] + d21[1]*d21[1] + d21[2]*d21[2];
  Ln = sqrt(Ln);
  
  // Compute engineering strain
  double strain = (Ln-Lo)/Lo;
  
  // Set material trial strain
  return theMaterial->setTrialStrain(strain);
}

const Matrix &
CorotTruss::getTangentStiff(void)
{
    static Matrix kl(3,3);

    // Material stiffness
    //
    // Get material tangent
    double EA = A*theMaterial->getTangent();
    EA /= (Ln * Ln * Lo);

    int i,j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            kl(i,j) = EA*d21[i]*d21[j];

    // Geometric stiffness
    //
    // Get material stress
    double q = A*theMaterial->getStress();
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
CorotTruss::getInitialStiff(void)
{
    static Matrix kl(3,3);

    // Material stiffness
    kl.Zero();
    kl(0,0) = A * theMaterial->getInitialTangent() / Lo;

    // Compute R'*kl*R
    static Matrix kg(3,3);
    kg.addMatrixTripleProduct(0.0, R, kl, 1.0);

    Matrix &K = *theMatrix;
    K.Zero();

    // Copy stiffness into appropriate blocks in element stiffness
    int numDOF2 = numDOF/2;
    for (int i = 0; i < numDIM; i++) {
        for (int j = 0; j < numDIM; j++) {
            K(i,j)                 =  kg(i,j);
            K(i,j+numDOF2)         = -kg(i,j);
            K(i+numDOF2,j)         = -kg(i,j);
            K(i+numDOF2,j+numDOF2) =  kg(i,j);
        }
    }

    return *theMatrix;
}


const Matrix &
CorotTruss::getMass(void)
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
CorotTruss::zeroLoad(void)
{
	return;
}

int 
CorotTruss::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "CorotTruss::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  
  return -1;
}



int 
CorotTruss::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

const Vector &
CorotTruss::getResistingForce()
{
	// Get material stress
	double SA = A*theMaterial->getStress();
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
    for (int i = 0; i < numDIM; i++) {
        P(i)         = -qg(i);
        P(i+numDOF2) =  qg(i);
    }

    return *theVector;
}



const Vector &
CorotTruss::getResistingForceIncInertia()
{	
    Vector &P = *theVector;
    P = this->getResistingForce();
    
    if (rho != 0.0) {
	
      const Vector &accel1 = theNodes[0]->getTrialAccel();
      const Vector &accel2 = theNodes[1]->getTrialAccel();	
      
      double M = 0.5*rho*Lo;
      int numDOF2 = numDOF/2;
      for (int i = 0; i < numDIM; i++) {
	P(i)        += M*accel1(i);
	P(i+numDOF2) += M*accel2(i);
      }
    }

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      *theVector += this->getRayleighDampingForces();

    return *theVector;
}

int
CorotTruss::sendSelf(int commitTag, Channel &theChannel)
{
  int res;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // truss packs it's data into a Vector and sends this to theChannel
  // along with it's dbTag and the commitTag passed in the arguments

  static Vector data(7);
  data(0) = this->getTag();
  data(1) = numDIM;
  data(2) = numDOF;
  data(3) = A;
  data(6) = rho;
  
  data(4) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();

  // NOTE: we do have to ensure that the material has a database
  // tag if we are sending to a database channel.
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    if (matDbTag != 0)
      theMaterial->setDbTag(matDbTag);
  }
  data(5) = matDbTag;

  res = theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING Truss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -1;
  }	      

  // truss then sends the tags of it's two end nodes

  res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr << "WARNING Truss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return -2;
  }

  // finally truss asks it's material object to send itself
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "WARNING Truss::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    return -3;
  }

  return 0;
}

int
CorotTruss::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res;
  int dataTag = this->getDbTag();

  // truss creates a Vector, receives the Vector and then sets the 
  // internal data with the data in the Vector

  static Vector data(7);
  res = theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING Truss::recvSelf() - failed to receive Vector\n";
    return -1;
  }	      

  this->setTag((int)data(0));
  numDIM = (int)data(1);
  numDOF = (int)data(2);
  A = data(3);
  rho = data(6);
  
  // truss now receives the tags of it's two external nodes
  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
  if (res < 0) {
    opserr << "WARNING Truss::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return -2;
  }

  // finally truss creates a material object of the correct type,
  // sets its database tag and asks this new object to recveive itself.

  int matClass = (int)data(4);
  int matDb = (int)data(5);

  // check if we have a material object already & if we do if of right type
  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

    // if old one .. delete it
    if (theMaterial != 0)
      delete theMaterial;

    // create a new material object
    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
    if (theMaterial == 0) {
      opserr << "WARNING Truss::recvSelf() - " << this->getTag() << 
	"failed to get a blank Material of type: " << matClass << endln;
      return -3;
    }
  }

  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "WARNING Truss::recvSelf() - " << this->getTag() << " failed to receive its Material\n";
    return -3;    
  }

  return 0;
}

int
CorotTruss::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	// ensure setDomain() worked
	if (Ln == 0.0)
		return 0;

	// first determine the two end points of the CorotTruss based on
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
CorotTruss::Print(OPS_Stream &s, int flag)
{
	s << "\nCorotTruss, tag: " << this->getTag() << endln;
	s << "\tConnected Nodes: " << connectedExternalNodes;
	s << "\tSection Area: " << A << endln;
	s << "\tUndeformed Length: " << Lo << endln;
	s << "\tCurrent Length: " << Ln << endln;
	s << "\tMass Density/Length: " << rho << endln;
	s << "\tRotation matrix: " << endln;

	if (theMaterial) {
		s << "\tAxial Force: " << A*theMaterial->getStress() << endln;
		s << "\tUniaxialMaterial, tag: " << theMaterial->getTag() << endln;
		theMaterial->Print(s,flag);
	}
}

Response*
CorotTruss::setResponse(const char **argv, int argc, Information &eleInfo, OPS_Stream &output)
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

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"axialForce") == 0) {

    output.tag("ResponseType", "N");
    theResponse =  new ElementResponse(this, 1, 0.0);

  } else if (strcmp(argv[0],"defo") == 0 || strcmp(argv[0],"deformations") == 0 ||
	     strcmp(argv[0],"deformation") == 0) {

    output.tag("ResponseType", "eps");
    theResponse = new ElementResponse(this, 2, 0.0);

  // a material quantity    
  } else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"-material") == 0) {

    theResponse =  theMaterial->setResponse(&argv[1], argc-1, eleInfo, output);
  }

  output.endTag();
  return theResponse;
}

int 
CorotTruss::getResponse(int responseID, Information &eleInfo)
{
  switch (responseID) {
    case 1:
	  return eleInfo.setDouble(A * theMaterial->getStress());
      
    case 2:
	  return eleInfo.setDouble(Lo * theMaterial->getStrain());
      
    default:
	  return 0;
  }
}
