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
                                                                        
// Written: Yuan Lu and M. Panagiotou 2013
// Code Based on minor modifications to CorotTruss written MHS 2001
//
// Description: This file contains the class implementation for CorotTruss2.

#include <CorotTruss2.h>
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

Matrix CorotTruss2::M2(2,2);
Matrix CorotTruss2::M4(4,4);
Matrix CorotTruss2::M6(6,6);
Matrix CorotTruss2::M12(12,12);

Vector CorotTruss2::V2(2);
Vector CorotTruss2::V4(4);
Vector CorotTruss2::V6(6);
Vector CorotTruss2::V12(12);


#include <elementAPI.h>
#define OPS_Export 

void *
OPS_CorotTruss2(void)
{
  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs < 7) {
    opserr << "Invalid Args want: element CorotTruss2 $tag $iNode $jNode $auxN1 $auxN2 $A $matTag <-rho $rho>\n";
    return 0;	
  }

  int    iData[5];
  double A = 0.0;
  double rho = 0.0;
  int matTag = 0;
  int ndm = OPS_GetNDM();


  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer (tag, iNode, jNode, auxN1, auxN2) in element CorotTruss2 " << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDouble(&numData, &A) != 0) {
    opserr << "WARNING: Invalid A: element CorotTruss2 " << iData[0] << 
      " $iNode $jNode $A $matTag <-rho $rho> <-rayleig $flagh>\n";
    return 0;	
  }

  numData = 1;
  if (OPS_GetInt(&numData, &matTag) != 0) {
    opserr << "WARNING: Invalid matTag: element CorotTruss2 " << iData[0] << 
      " $iNode $jNode $A $matTag <-rho $rho> <-rayleig $flagh>\n";
    return 0;
  }

  UniaxialMaterial *theUniaxialMaterial = OPS_GetUniaxialMaterial(matTag);
    
  if (theUniaxialMaterial == 0) {
    opserr << "WARNING: Invalid material not found element CorotTruss2 " << iData[0] << " $iNode $jNode $auxN1 $auxN2 $A " << 
      matTag << " <-rho $rho> <-rayleigh $flagh>\n";
    return 0;
  }
  
  numRemainingArgs -= 7;
  while (numRemainingArgs > 1) {
    const char *argvS = OPS_GetString();
  
    if (strcmp(argvS,"-rho") == 0) {
      numData = 1;
      if (OPS_GetDouble(&numData, &rho) != 0) {
	opserr << "WARNING Invalid rho in element CorotTruss2 " << iData[0] << 
	  " $iNode $jNode $auxN1 $auxN2 $A $matTag <-rho $rho> <-rayleigh $flagh>\n";
	return 0;
      }
    } else {
      opserr << "WARNING: Invalid option " << argvS << "  in: element CorotTruss2 " << iData[0] << 
	" $iNode $jNode $auxN1 $auxN2 $A $matTag <-rho $rho> <-rayleigh $flagh>\n";
      return 0;
    }      
    numRemainingArgs -= 2;
  }

  //now create the ReinforcedConcretePlaneStress
  theElement = new CorotTruss2(iData[0], ndm, iData[1], iData[2], iData[3], iData[4], *theUniaxialMaterial, A, rho);

  if (theElement == 0) {
    opserr << "WARNING: out of memory: element CorotTruss2 " << iData[0] << 
      " $iNode $jNode $A $matTag <-rho $rho> \n";
  }

  return theElement;
}


// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the CorotTruss2 end nodes.
CorotTruss2::CorotTruss2(int tag, int dim,
			   int Nd1, int Nd2, int oNd1, int oNd2, 
			   UniaxialMaterial &theMat,
			   double a, double r)
  :Element(tag,ELE_TAG_CorotTruss2),     
  theMaterial(0), theBetaMaterial(0), 
  connectedExternalNodes(2), connectedExternalOtherNodes(2),
  numDOF(0), numDIM(dim),
  Lo(0.0), Ln(0.0), otherLength(0.0),
  A(a), rho(r), R(3,3),
  theMatrix(0), theVector(0)
{
  // get a copy of the material and check we obtained a valid copy
  theMaterial = theMat.getCopy();
  if (theMaterial == 0) {
    opserr << "FATAL CorotTruss2::CorotTruss2 - " <<  tag <<
      "failed to get a copy of material with tag " << theMat.getTag() << endln;
    exit(-1);
  } else if (theMaterial->getClassTag() == MAT_TAG_ConcretewBeta) {
	theBetaMaterial = (ConcretewBeta *) theMaterial;
  }
  
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2 || connectedExternalOtherNodes.Size() != 2) {
    opserr << "FATAL CorotTruss2::CorotTruss2 - " <<  tag <<
      " failed to create an ID of size 2\n";
    exit(-1);
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;  
  /// some work to save the other nodes:
  connectedExternalOtherNodes(0) = oNd1;
  connectedExternalOtherNodes(1) = oNd2;

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
  theOtherNodes[0] = 0;
  theOtherNodes[1] = 0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
CorotTruss2::CorotTruss2()
  :Element(0,ELE_TAG_CorotTruss2),     
  theMaterial(0), theBetaMaterial(0), 
  connectedExternalNodes(2), connectedExternalOtherNodes(2),
  numDOF(0), numDIM(0),
  Lo(0.0), Ln(0.0), otherLength(0.0),
  A(0.0), rho(0.0), R(3,3),
  theMatrix(0), theVector(0)
{
  // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2 || connectedExternalOtherNodes.Size() != 2) {
    opserr << "FATAL CorotTruss2::CorotTruss2 - failed to create an ID of size 2\n";
    exit(-1);
  }

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
  theOtherNodes[0] = 0;
  theOtherNodes[1] = 0;
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
CorotTruss2::~CorotTruss2()
{
  // invoke the destructor on any objects created by the object
  // that the object still holds a pointer to
  if (theMaterial != 0)
    delete theMaterial;
}

int
CorotTruss2::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
CorotTruss2::getExternalNodes(void) 
{
	return connectedExternalNodes;
}

Node **
CorotTruss2::getNodePtrs(void) 
{
  return theNodes;
}

int
CorotTruss2::getNumDOF(void) 
{
	return numDOF;
}

// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the CorotTruss2 element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
CorotTruss2::setDomain(Domain *theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    Lo = 0.0;
    Ln = 0.0;

	theOtherNodes[0] = 0;
	theOtherNodes[1] = 0;
	otherLength = 0.0;
    return;
  }
  
  // first set the node pointers
  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);
  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);	
  // and the other node pointers (will be reordered later)
  int oNd1 = connectedExternalOtherNodes(0);
  int oNd2 = connectedExternalOtherNodes(1);
  theOtherNodes[0] = theDomain->getNode(oNd1);
  theOtherNodes[1] = theDomain->getNode(oNd2);	
  
  // if can't find both - send a warning message
  if ((theNodes[0] == 0) || (theNodes[1] == 0) || theOtherNodes[0] == 0 || theOtherNodes[1] == 0) {
      if (theNodes[0] == 0)
	opserr <<"Truss2::setDomain() - truss" << this->getTag() << " node " << Nd1 <<
	  " does not exist in the model\n";
      else if (theNodes[1] == 0)
	opserr <<"Truss2::setDomain() - truss" << this->getTag() << " node " << Nd2 <<
	  " does not exist in the model\n";
	  else if (theOtherNodes[0] == 0)
	opserr <<"Truss2::setDomain() - truss" << this->getTag() << " node " << oNd1 <<
	  " does not exist in the model\n";
	  else 
	opserr <<"Truss2::setDomain() - truss" << this->getTag() << " node " << oNd2 <<
	  " does not exist in the model\n";

      // fill this in so don't segment fault later
      numDOF = 6;
	  theMatrix = &M6;
      theVector = &V6;
      return;
    }

  
  // now determine the number of dof and the dimesnion    
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();	
  
  // if differing dof at the ends - print a warning message
  if (dofNd1 != dofNd2) {
    opserr << "WARNING CorotTruss2::setDomain(): nodes " << Nd1 <<
      " and " << Nd2 << "have differing dof at ends for CorotTruss2 " << this->getTag() << endln;
    
    // fill this in so don't segment fault later
    numDOF = 6;    
    theMatrix = &M6;
    theVector = &V6;
    return;
  }	
  
  if (numDIM == 1 && dofNd1 == 1) {
    numDOF = 2;
    theMatrix = &M2;
    theVector = &V2;
  }
  else 
  if (numDIM == 2 && dofNd1 == 2) {
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
    opserr << " CorotTruss2::setDomain -- nodal DOF " << dofNd1 << " not compatible with element\n";
    
    // fill this in so don't segment fault later
    numDOF = 6;    
    theMatrix = &M6;
    theVector = &V6;	
    return;
  }

	// call the base class method
	this->DomainComponent::setDomain(theDomain);

	// now determine the length, cosines and fill in the transformation
	// NOTE t = -t(every one else uses for residual calc)
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();
	double dx = end2Crd(0)-end1Crd(0);
	double dy = end2Crd(1)-end1Crd(1);	
	double dz = end2Crd(2)-end1Crd(2);		

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

	// set up the normal strain nodes.
	 // determine the strain
    const Vector &Oend1Crd = theOtherNodes[0]->getCrds();
    const Vector &Oend2Crd = theOtherNodes[1]->getCrds();	
	double dx2 = Oend2Crd(0)-Oend1Crd(0);
	double dy2 = Oend2Crd(1)-Oend1Crd(1);	
	double dz2 = Oend2Crd(2)-Oend1Crd(2);	

	od21[0] = 0.0;
	od21[1] = 0.0;
	od21[2] = 0.0;

	// Update offsets in basic system due to nodal displacements
	for (int i = 0; i < numDIM; i++) {
		double deltaDisp = Oend1Crd(i) - Oend2Crd(i);
		od21[0] += deltaDisp*R(0,i);
		od21[1] += deltaDisp*R(1,i);
		od21[2] += deltaDisp*R(2,i);
	}

	otherLength = sqrt(od21[0]*od21[0] + od21[1]*od21[1] + od21[2]*od21[2]);
	theta = acos((dx*dx2+dy*dy2+dz*dz2)/(Lo*otherLength)); // cos*theta = a*b/(|a|*|b|)
}

int
CorotTruss2::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "CorotTruss2::commitState () - failed in base class\n";
  }    
  retVal = theMaterial->commitState();
  return retVal;
}

int
CorotTruss2::revertToLastCommit()
{
	// Revert the material
	return theMaterial->revertToLastCommit();
}

int
CorotTruss2::revertToStart()
{
	// Revert the material to start
	return theMaterial->revertToStart();
}

int
CorotTruss2::update(void)
{
  // Nodal displacements
  const Vector &end1Disp  = theNodes[0]->getTrialDisp();
  const Vector &end2Disp  = theNodes[1]->getTrialDisp();    
  const Vector &end1Vel   = theNodes[0]->getTrialVel();
  const Vector &end2Vel   = theNodes[1]->getTrialVel();	
  
  // Initial offsets
  d21[0] = Lo; d21[1] = d21[2] = 0.0;
  v21[0] = v21[1] = v21[2] = 0.0;
  
  // Update offsets in basic system due to nodal displacements
  for (int i = 0; i < numDIM; i++) {
    double deltaDisp = end2Disp(i) - end1Disp(i);
    d21[0] += deltaDisp*R(0,i);
    d21[1] += deltaDisp*R(1,i);
    d21[2] += deltaDisp*R(2,i);
    double deltaVel = end2Vel(i) - end1Vel(i);
    v21[0] += deltaVel*R(0,i);
    v21[1] += deltaVel*R(1,i);
    v21[2] += deltaVel*R(2,i);
  }
  
  // Compute new length
  Ln = sqrt(d21[0]*d21[0] + d21[1]*d21[1] + d21[2]*d21[2]);
  
  // Compute engineering strain and strain rate
  double strain = (Ln - Lo)/Lo;
  double rate = (d21[0]*v21[0] + d21[1]*v21[1] + d21[2]*v21[2])/Ln/Lo;
  
  // Set material trial strain
  if (theBetaMaterial && theta != 0.0) {
  	double strain_t = this->computeCurrentNormalStrain();
	strain_t = (strain_t-strain*fabs(cos(theta)))/(fabs(sin(theta)));
	return theBetaMaterial->setTrialStrainwBeta(strain, strain_t, rate);
  } else {
	return theMaterial->setTrialStrain(strain, rate);
  }
}

double
CorotTruss2::computeCurrentNormalStrain(void)
{
    // normal vector = (-cosX[1], cosX[0])
	if (otherLength == 0) {
		return 0;
	}

    // determine the strain
    const Vector &disp1 = theOtherNodes[0]->getTrialDisp();
    const Vector &disp2 = theOtherNodes[1]->getTrialDisp();	

      // Initial offsets
	double temp[3];
	temp[0] = od21[0];
	temp[1] = od21[1];
	temp[2] = od21[2];
  
	// Update offsets in basic system due to nodal displacements
	for (int i = 0; i < numDIM; i++) {
		double deltaDisp = disp1(i) - disp2(i);
		temp[0] += deltaDisp*R(0,i);
		temp[1] += deltaDisp*R(1,i);
		temp[2] += deltaDisp*R(2,i);
	}
  
    // this method should never be called with L == 0
	otherLength_new = sqrt(temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2]);
	double strain_t = (otherLength_new - otherLength)/otherLength;

    return strain_t;
}

const Matrix &
CorotTruss2::getTangentStiff(void)
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
CorotTruss2::getInitialStiff(void)
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
CorotTruss2::getMass(void)
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
CorotTruss2::zeroLoad(void)
{
	return;
}

int 
CorotTruss2::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "CorotTruss2::addLoad - load type unknown for truss with tag: " << this->getTag() << endln;
  
  return -1;
}



int 
CorotTruss2::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

const Vector &
CorotTruss2::getResistingForce()
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
CorotTruss2::getResistingForceIncInertia()
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
CorotTruss2::sendSelf(int commitTag, Channel &theChannel)
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
  res = theChannel.sendID(dataTag, commitTag, connectedExternalOtherNodes);
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
CorotTruss2::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
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
  res = theChannel.recvID(dataTag, commitTag, connectedExternalOtherNodes);
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
    } else if (theMaterial->getClassTag() == MAT_TAG_ConcretewBeta) {
		theBetaMaterial = (ConcretewBeta *) theMaterial;
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
CorotTruss2::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// ensure setDomain() worked
	if (Ln == 0.0)
		return 0;

    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

void
CorotTruss2::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << "\nCorotTruss2, tag: " << this->getTag() << endln;
        s << "\tConnected Nodes: " << connectedExternalNodes;
        s << "\tSection Area: " << A << endln;
        s << "\tUndeformed Length: " << Lo << endln;
        s << "\tCurrent Length: " << Ln << endln;
        s << "\tMass Density/Length: " << rho << endln;
        s << "\tRotation matrix: " << endln;
        
        if (theMaterial) {
            s << "\tAxial Force: " << A*theMaterial->getStress() << endln;
            s << "\tUniaxialMaterial, tag: " << theMaterial->getTag() << endln;
            theMaterial->Print(s, flag);
        }
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"CorotTruss2\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"A\": " << A << ", ";
        s << "\"massperlength\": " << rho << ", ";
        s << "\"material\": \"" << theMaterial->getTag() << "\"}";
    }
}

Response*
CorotTruss2::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","Truss");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    //
    // we compare argv[0] for known response types for the CorotTruss2
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

            // a material quantity
    }
    else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "-material") == 0) {
        output.tag("GaussPointOutput");
        output.attr("number", 1);
        output.attr("eta", 0.0);

        if (argc > 1) {
            // we need at least one more argument otherwise 
            // there is no need to forward this call to the material
            if (argc > 2) {
                // if we have 2 or more extra arguments, the first one 
                // could be an integer. In this case we check to see if it is the section id
                // (only 1 in this case)
                int sectionNum = atoi(argv[1]);
                if (sectionNum == 0) {
                    // if it is not a number we forward the call to the section as usual
                    theResponse = theMaterial->setResponse(&argv[1], argc - 1, output);
                }
                else {
                    // it is a number. Now we have to make sure it is within the allowed range
                    // for this element (in this case it can only be 1)
                    // If it is > 1, then we MUST return NULL, because the MPCO recorder iteratively
                    // uses this call to understand how many fibers we have in a section
                    if (sectionNum == 1) {
                        theResponse = theMaterial->setResponse(&argv[2], argc - 2, output);
                    }
                }
            }
            else {
                // otherwise forward it as usual
                theResponse = theMaterial->setResponse(&argv[1], argc - 1, output);
            }
        }
        output.endTag();
    }

    output.endTag();
    return theResponse;
}

int 
CorotTruss2::getResponse(int responseID, Information &eleInfo)
{
    double strain;

    switch (responseID) {
    case 1:
        return eleInfo.setVector(this->getResistingForce());

    case 2:
        return eleInfo.setDouble(A * theMaterial->getStress());

    case 3:
        if (Lo == 0.0) {
            strain = 0.0;
        } else {
            strain = theMaterial->getStrain();
        }
        return eleInfo.setDouble(Lo * strain);

    default:
        return 0;
    }
}
