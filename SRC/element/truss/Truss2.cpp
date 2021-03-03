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


// Written Y.Lu and M.Panagiotou 2013
//  minor mod to Truss written by fmk 1997

//
// Description: This file contains the implementation for the Truss2 class.
//

#include <Truss2.h>
#include <Information.h>
#include <Parameter.h>

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

//#include <fstream>

// initialise the class wide variables
Matrix Truss2::trussM2(2,2);
Matrix Truss2::trussM4(4,4);
Matrix Truss2::trussM6(6,6);
Matrix Truss2::trussM12(12,12);
Vector Truss2::trussV2(2);
Vector Truss2::trussV4(4);
Vector Truss2::trussV6(6);
Vector Truss2::trussV12(12);

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the truss end nodes.

#include <elementAPI.h>
#define OPS_Export 

void *
OPS_Truss2(void)
{
	Element *theElement = 0;

	int numRemainingArgs = OPS_GetNumRemainingInputArgs();

	if (numRemainingArgs < 7) {
		opserr << "Invalid Args want: element Truss2 $tag $iNode $jNode $auxN1 $auxN2 $A $matTag <-rho $rho> <-rayleigh $flag>\n";
		return 0;	
	}

	int    iData[5];
	double A = 0.0;
	double rho = 0.0;
	int matTag = 0;
	int doRayleigh = 0;
	int ndm = OPS_GetNDM();

	int numData = 5;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid integer (tag, iNode, jNode, auxN1, auxN2) in element Truss2 " << endln;
		return 0;
	}

	numData = 1;
	if (OPS_GetDouble(&numData, &A) != 0) {
		opserr << "WARNING: Invalid A: element Truss2 " << iData[0] << 
			" $iNode $jNode $auxN1 $auxN2 $A $matTag <-rho $rho> <-rayleig $flagh>\n";
		return 0;	
	}

	numData = 1;
	if (OPS_GetInt(&numData, &matTag) != 0) {
		opserr << "WARNING: Invalid matTag: element Truss2 " << iData[0] << 
			" $iNode $jNode $auxN1 $auxN2 $A $matTag <-rho $rho> <-rayleig $flagh>\n";
		return 0;
	}

	UniaxialMaterial *theUniaxialMaterial = OPS_GetUniaxialMaterial(matTag);

	if (theUniaxialMaterial == 0) {
		opserr << "WARNING: Invalid material not found element Truss2 " << iData[0] << " $iNode $jNode $auxN1 $auxN2 $A " << 
			matTag << " <-rho $rho> <-rayleig $flagh>\n";
		return 0;
	}

	numRemainingArgs -= 7;
	while (numRemainingArgs > 1) {
	  const char *argvS = OPS_GetString();

		if (strcmp(argvS,"-rho") == 0) {
			numData = 1;
			if (OPS_GetDouble(&numData, &rho) != 0) {
				opserr << "WARNING Invalid rho in element Truss " << iData[0] << 
					" $iNode $jNode $A $matTag <-rho $rho> <-doRayleigh $flagh>\n";
				return 0;
			}
		} else if (strcmp(argvS,"-doRayleigh") == 0) {
			numData = 1;
			if (OPS_GetInt(&numData, &doRayleigh) != 0) {
				opserr << "WARNING: Invalid doRayleigh in element Truss " << iData[0] << 
					" $iNode $jNode $A $matTag <-rho $rho> <-doRayleigh $flagh>\n";
				return 0;
			}
		} else {
			opserr << "WARNING: Invalid option " << argvS << "  in: element Truss " << iData[0] << 
				" $iNode $jNode $A $matTag <-rho $rho> <-doRayleigh $flagh>\n";
			return 0;
		}      
		numRemainingArgs -= 2;
	}

	//now create the ReinforcedConcretePlaneStress
	theElement = new Truss2(iData[0], ndm, iData[1], iData[2], iData[3], iData[4], *theUniaxialMaterial, A, rho, doRayleigh);

	if (theElement == 0) {
		opserr << "WARNING: out of memory: element Truss2 " << iData[0] << 
			" $iNode $jNode $auxN1 $auxN2 $A $matTag <-rho $rho>\n";
	}

	return theElement;
}

Truss2::Truss2(int tag, 
	int dim,
	int Nd1, int Nd2, int oNd1, int oNd2, 
	UniaxialMaterial &theMat,
	double a, double r, int damp)
	:Element(tag,ELE_TAG_Truss2),     
	theMaterial(0), theBetaMaterial(0), connectedExternalNodes(2), connectedExternalOtherNodes(2),
	dimension(dim), numDOF(0), theLoad(0),
	theMatrix(0), theVector(0),
	L(0.0), A(a), rho(r), doRayleighDamping(damp)
{
	// get a copy of the material and check we obtained a valid copy
	theMaterial = theMat.getCopy();
	if (theMaterial == 0) {
	  opserr << "FATAL Truss2::Truss2 - " << tag <<
	    "failed to get a copy of material with tag " << theMat.getTag() << endln;
	  exit(-1);
	} else if (theMaterial->getClassTag() == MAT_TAG_ConcretewBeta) {
	  theBetaMaterial = (ConcretewBeta *) theMaterial;
	}

	// ensure the connectedExternalNode ID is of correct size & set values
	if (connectedExternalNodes.Size() != 2 || connectedExternalOtherNodes.Size() != 2) {
		opserr << "FATAL Truss2::Truss2 - " <<  tag << "failed to create an ID of size 2\n";
		exit(-1);
	}

	connectedExternalNodes(0) = Nd1;
	connectedExternalNodes(1) = Nd2; 

	/// some work to save the other nodes:
	connectedExternalOtherNodes(0) = oNd1;
	connectedExternalOtherNodes(1) = oNd2;

	// set node pointers to NULL
	for (int i=0; i<2; i++) {
		theNodes[i] = 0;
		theOtherNodes[i] = 0;
	}

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
Truss2::Truss2()
	:Element(0,ELE_TAG_Truss2),     
	theMaterial(0),connectedExternalNodes(2), connectedExternalOtherNodes(2),
	dimension(0), numDOF(0), theLoad(0),
	theMatrix(0), theVector(0),
	L(0.0), A(0.0), rho(0.0)
{

	// ensure the connectedExternalNode ID is of correct size 
	if (connectedExternalNodes.Size() != 2 || connectedExternalOtherNodes.Size() != 2) {
		opserr << "FATAL Truss2::Trus2s - failed to create an ID of size 2\n";
		exit(-1);
	}

	for (int i=0; i<2; i++) {
		theNodes[i] = 0;
		theOtherNodes[i] = 0;
	}

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
Truss2::~Truss2()
{
	// invoke the destructor on any objects created by the object
	// that the object still holds a pointer to
	if (theMaterial != 0)
		delete theMaterial;
	if (theLoad != 0)
		delete theLoad;
	if (theLoadSens != 0)
		delete theLoadSens;
}


int
	Truss2::getNumExternalNodes(void) const
{
	return 2;
}

const ID &
	Truss2::getExternalNodes(void) 
{
	return connectedExternalNodes;
}

Node **
	Truss2::getNodePtrs(void) 
{
	return theNodes;
}

int
	Truss2::getNumDOF(void) 
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
	Truss2::setDomain(Domain *theDomain)
{
	// check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		L = 0;

		theOtherNodes[0] = 0;
		theOtherNodes[1] = 0;
		otherLength = 0;
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
		numDOF = 2;    
		theMatrix = &trussM2;
		theVector = &trussV2;	

		return;
	}

	// now determine the number of dof and the dimesnion    
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();	

	// if differing dof at the ends - print a warning message
	if (dofNd1 != dofNd2) {
		opserr <<"WARNING Truss2::setDomain(): nodes " << Nd1 << " and " << Nd2 <<
			"have differing dof at ends for truss " << this->getTag() << endln;

		// fill this in so don't segment fault later
		numDOF = 2;    
		theMatrix = &trussM2;
		theVector = &trussV2;	

		return;
	}	

	// call the base class method
	this->DomainComponent::setDomain(theDomain);

	// now set the number of dof for element and set matrix and vector pointer
	if (dimension == 1 && dofNd1 == 1) {
		numDOF = 2;    
		theMatrix = &trussM2;
		theVector = &trussV2;
	}
	else 
		if (dimension == 2 && dofNd1 == 2) {
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
			opserr <<"WARNING Truss2::setDomain cannot handle " << dimension << " dofs at nodes in " << 
				dofNd1  << " problem\n";

			numDOF = 2;    
			theMatrix = &trussM2;
			theVector = &trussV2;	
			return;
		}

		if (theLoad == 0)
			theLoad = new Vector(numDOF);
		else if (theLoad->Size() != numDOF) {
			delete theLoad;
			theLoad = new Vector(numDOF);
		}

		if (theLoad == 0) {
			opserr << "Truss2::setDomain - truss " << this->getTag() << 
				"out of memory creating vector of size" << numDOF << endln;
			exit(-1);
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
				opserr <<"WARNING Truss2::setDomain() - truss " << this->getTag() << " has zero length\n";
				return;
			}

			cosX[0] = 1.0;

			const Vector &end1Crd2 = theOtherNodes[0]->getCrds();
			const Vector &end2Crd2 = theOtherNodes[1]->getCrds();	
			dx = end2Crd2(0)-end1Crd2(0);
			otherLength = sqrt(dx*dx);
			otherCosX[0] = 0;

		} else if (dimension == 2) {
			double dx = end2Crd(0)-end1Crd(0);
			double dy = end2Crd(1)-end1Crd(1);	

			L = sqrt(dx*dx + dy*dy);

			if (L == 0.0) {
				opserr <<"WARNING Truss2::setDomain() - truss " << this->getTag() << " has zero length\n";
				return;
			}

			cosX[0] = dx/L;
			cosX[1] = dy/L;

			const Vector &end1Crd2 = theOtherNodes[0]->getCrds();
			const Vector &end2Crd2 = theOtherNodes[1]->getCrds();	
			double dx2 = end2Crd2(0)-end1Crd2(0);
			double dy2 = end2Crd2(1)-end1Crd2(1);	

			otherLength = sqrt(dx2*dx2 + dy2*dy2);

			if (otherLength == 0.0) {
				opserr <<"WARNING Truss2::setDomain() - truss " << this->getTag() << " has auxiliary nodes that are the same point\n";
				otherCosX[0] = 0;
				otherCosX[1] = 0;
				return;
			} 

			otherCosX[0] = dx2/otherLength;
			otherCosX[1] = dy2/otherLength;

			// project on normal of truss direction.
			theta = acos((dx*dx2+dy*dy2)/(L*otherLength)); // cos*theta = a*b/(|a|*|b|)
			if (theta == 0.0) {
				opserr << "WARNING Truss2::setDomain() - truss2 " << this->getTag() << " has theta = 0, disabling biaxial effects\n";
			}
		}

		else {
			double dx = end2Crd(0)-end1Crd(0);
			double dy = end2Crd(1)-end1Crd(1);	
			double dz = end2Crd(2)-end1Crd(2);		

			L = sqrt(dx*dx + dy*dy + dz*dz);

			if (L == 0.0) {
				opserr <<"WARNING Truss2::setDomain() - truss " << this->getTag() << " has zero length\n";
				return;
			}

			cosX[0] = dx/L;
			cosX[1] = dy/L;
			cosX[2] = dz/L;

			const Vector &end1Crd2 = theOtherNodes[0]->getCrds();
			const Vector &end2Crd2 = theOtherNodes[1]->getCrds();	
			double dx2 = end2Crd2(0)-end1Crd2(0);
			double dy2 = end2Crd2(1)-end1Crd2(1);	
			double dz2 = end2Crd2(2)-end1Crd2(2);	

			otherLength = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);

			if (otherLength == 0.0) {
				opserr <<"WARNING Truss2::setDomain() - truss " << this->getTag() << " has auxiliary nodes that are the same point\n";
				otherCosX[0] = 0;
				otherCosX[1] = 0;
				otherCosX[2] = 0;
				return;
			} 

			otherCosX[0] = dx2/otherLength;
			otherCosX[1] = dy2/otherLength;
			otherCosX[2] = dz2/otherLength;

			// project on normal of truss direction.
			theta = acos((dx*dx2+dy*dy2+dz*dz2)/(L*otherLength)); // cos*theta = a*b/(|a|*|b|)
			if (theta == 0.0) {
				opserr << "WARNING Truss2::setDomain() - truss2 " << this->getTag() << " has theta = 0, disabling biaxial effects\n";
			}
		}
}   	 


int
	Truss2::commitState()
{
	int retVal = 0;
	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "Truss2::commitState () - failed in base class";
	}    
	retVal = theMaterial->commitState();
	return retVal;
}

int
	Truss2::revertToLastCommit()
{
	return theMaterial->revertToLastCommit();
}

int
	Truss2::revertToStart()
{
	return theMaterial->revertToStart();
}

int
	Truss2::update(void)
{
	// determine the current strain given trial displacements at nodes
	double strain = this->computeCurrentStrain();
	double rate = this->computeCurrentStrainRate();
	if (theBetaMaterial && theta != 0.0) {
		double strain_t = this->computeCurrentNormalStrain();
		strain_t = (strain_t-strain*fabs(cos(theta)))/(fabs(sin(theta)));
		return theBetaMaterial->setTrialStrainwBeta(strain, strain_t, rate);
	} else {
		return theMaterial->setTrialStrain(strain, rate);
	}
}


const Matrix &
	Truss2::getTangentStiff(void)
{
	if (L == 0.0) { // - problem in setDomain() no further warnings
		theMatrix->Zero();
		return *theMatrix;
	}

	double E = theMaterial->getTangent();

	// come back later and redo this if too slow
	Matrix &stiff = *theMatrix;

	int numDOF2 = numDOF/2;
	double temp;
	double EAoverL = E*A/L;
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			temp = cosX[i]*cosX[j]*EAoverL;
			stiff(i,j) = temp;
			stiff(i+numDOF2,j) = -temp;
			stiff(i,j+numDOF2) = -temp;
			stiff(i+numDOF2,j+numDOF2) = temp;
		}
	}

	return stiff;
}


const Matrix &
	Truss2::getInitialStiff(void)
{
	if (L == 0.0) { // - problem in setDomain() no further warnings
		theMatrix->Zero();
		return *theMatrix;
	}

	double E = theMaterial->getInitialTangent();

	// come back later and redo this if too slow
	Matrix &stiff = *theMatrix;

	int numDOF2 = numDOF/2;
	double temp;
	double EAoverL = E*A/L;
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			temp = cosX[i]*cosX[j]*EAoverL;
			stiff(i,j) = temp;
			stiff(i+numDOF2,j) = -temp;
			stiff(i,j+numDOF2) = -temp;
			stiff(i+numDOF2,j+numDOF2) = temp;
		}
	}

	return *theMatrix;
}

const Matrix &
	Truss2::getDamp(void)
{
	if (L == 0.0) { // - problem in setDomain() no further warnings
		theMatrix->Zero();
		return *theMatrix;
	}

	theMatrix->Zero();

	if (doRayleighDamping == 1)
		*theMatrix = this->Element::getDamp();

	double eta = theMaterial->getDampTangent();

	// come back later and redo this if too slow
	Matrix &damp = *theMatrix;

	int numDOF2 = numDOF/2;
	double temp;
	double etaAoverL = eta*A/L;
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			temp = cosX[i]*cosX[j]*etaAoverL;
			damp(i,j) += temp;
			damp(i+numDOF2,j) += -temp;
			damp(i,j+numDOF2) += -temp;
			damp(i+numDOF2,j+numDOF2) += temp;
		}
	}

	return damp;
}


const Matrix &
	Truss2::getMass(void)
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
	Truss2::zeroLoad(void)
{
	theLoad->Zero();
}

int 
	Truss2::addLoad(ElementalLoad *theLoad, double loadFactor)

{  
	opserr <<"Truss2::addLoad - load type unknown for truss with tag: " << this->getTag() << endln; 
	return -1;
}

int 
	Truss2::addInertiaLoadToUnbalance(const Vector &accel)
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
		opserr <<"Truss2::addInertiaLoadToUnbalance " <<
			"matrix and vector sizes are incompatible\n";
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


int 
	Truss2::addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool somethingRandomInMotions)
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
			opserr << "Truss2::addInertiaLoadToUnbalance " <<
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
			opserr << "Truss2::addInertiaLoadToUnbalance " <<
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

const Vector &
	Truss2::getResistingForce()
{	
	if (L == 0.0) { // - problem in setDomain() no further warnings
		theVector->Zero();
		return *theVector;
	}

	// R = Ku - Pext
	// Ku = F * transformation
	double force = A*theMaterial->getStress();
	int numDOF2 = numDOF/2;
	double temp;
	for (int i = 0; i < dimension; i++) {
		temp = cosX[i]*force;
		(*theVector)(i) = -temp;
		(*theVector)(i+numDOF2) = temp;
	}
    
	return *theVector;
}


const Vector &
	Truss2::getResistingForceIncInertia()
{	
	this->getResistingForce();
    
	// subtract external load
	(*theVector) -= *theLoad;
    
	// now include the mass portion
	if (L != 0.0 && rho != 0.0) {

		const Vector &accel1 = theNodes[0]->getTrialAccel();
		const Vector &accel2 = theNodes[1]->getTrialAccel();	

		int numDOF2 = numDOF/2;
		double M = 0.5*rho*L;
		for (int i = 0; i < dimension; i++) {
			(*theVector)(i) += M*accel1(i);
			(*theVector)(i+numDOF2) += M*accel2(i);
		}

		// add the damping forces if rayleigh damping
		if (doRayleighDamping == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
			(*theVector) += this->getRayleighDampingForces();
	}  else {

		// add the damping forces if rayleigh damping
		if (doRayleighDamping == 1 && (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
			(*theVector) += this->getRayleighDampingForces();
	}

	return *theVector;
}

int Truss2::sendSelf(int commitTag, Channel &theChannel)
{
	// opserr << "Truss2::sendSelf - start\n";

	int res;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// truss packs it's data into a Vector and sends this to theChannel
	// along with it's dbTag and the commitTag passed in the arguments

	static Vector data(8);
	data(0) = this->getTag();
	data(1) = dimension;
	data(2) = numDOF;
	data(3) = A;
	data(6) = rho;
	if (doRayleighDamping == 0)
		data(7) = 0;
	else
		data(7) = 1;

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
		opserr <<"WARNING Truss2::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return -1;
	}	      

	// truss then sends the tags of it's two end nodes

	res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr <<"WARNING Truss2::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return -2;
	}
	res = theChannel.sendID(dataTag, commitTag, connectedExternalOtherNodes);
	if (res < 0) {
		opserr <<"WARNING Truss2::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return -2;
	}

	// finally truss asks it's material object to send itself
	res = theMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr <<"WARNING Truss2::sendSelf() - " << this->getTag() << " failed to send its Material\n";
		return -3;
	}
	
	//static ID endData(1);
	//theChannel.recvID(dataTag, commitTag, endData);
	//opserr << "Truss2::sendSelf - done\n";
	return 0;
}

int
	Truss2::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	//opserr << "Truss2::recvSelf - start\n";
	int res;
	int dataTag = this->getDbTag();

	// truss creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector

	static Vector data(8);
	res = theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr <<"WARNING Truss2::recvSelf() - failed to receive Vector\n";
		return -1;
	}	      

	this->setTag((int)data(0));
	dimension = (int)data(1);
	numDOF = (int)data(2);
	A = data(3);
	rho = data(6);
	if (data(7) == 0)
		doRayleighDamping = 0;
	else
		doRayleighDamping = 1;

	// truss now receives the tags of it's two external nodes
	res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr <<"WARNING Truss2::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return -2;
	}
	// truss now receives the tags of it's two external nodes
	res = theChannel.recvID(dataTag, commitTag, connectedExternalOtherNodes);
	if (res < 0) {
		opserr <<"WARNING Truss2::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return -2;
	}
//	opserr << connectedExternalNodes << connectedExternalOtherNodes;

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
			opserr <<"WARNING Truss2::recvSelf() - " << this->getTag() 
				<< " failed to get a blank Material of type " << matClass << endln;
			return -3;
		} else if (theMaterial->getClassTag() == MAT_TAG_ConcretewBeta) {
			theBetaMaterial = (ConcretewBeta *) theMaterial;
		}
	}

	theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
	res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		opserr <<"WARNING Truss2::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
		return -3;    
	}
//	opserr << "Truss2::recvSelf - done\n";
//	static ID endData(1);
//	theChannel.sendID(dataTag, commitTag, endData);

	return 0;
}


int
Truss2::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// ensure setDomain() worked
	if (L == 0.0)
		return 0;

	// get display coordinates
	static Vector v1(3);
	static Vector v2(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);

	// determine color and draw line
	if (displayMode == 1 || displayMode == 2) {
		// compute the strain and axial force in the member
		double strain, force;
		if (L == 0.0) {
			strain = 0.0;
			force = 0.0;
		}
		else {
			strain = this->computeCurrentStrain();
			theMaterial->setTrialStrain(strain);
			force = A * theMaterial->getStress();
		}
		if (displayMode == 2) {// use the strain as the drawing measure
			return theViewer.drawLine(v1, v2, (float)strain, (float)strain);
		}
		else { // otherwise use the axial force as measure
			return theViewer.drawLine(v1, v2, (float)force, (float)force);
		}
	}
	else {
		return theViewer.drawLine(v1, v2, 1.0, 1.0);
	}
}


void
	Truss2::Print(OPS_Stream &s, int flag)
{
	// compute the strain and axial force in the member
	double strain, force;
	strain = theMaterial->getStrain();
	force = A * theMaterial->getStress();

    if (flag == OPS_PRINT_CURRENTSTATE) { // print everything
        s << "Element: " << this->getTag();
        s << " type: Truss2  iNode: " << connectedExternalNodes(0);
        s << " jNode: " << connectedExternalNodes(1);
        s << " Area: " << A << " Mass/Length: " << rho;
        
        s << " \n\t strain: " << strain;
        s << " axial load: " << force;
        if (L != 0.0) {
            int numDOF2 = numDOF / 2;
            double temp;
            for (int i = 0; i < dimension; i++) {
                temp = cosX[i] * force;
                (*theVector)(i) = -temp;
                (*theVector)(i + numDOF2) = temp;
            }
            s << " \n\t unbalanced load: " << *theVector;
        }
        
        s << " \t Material: " << *theMaterial;
        s << endln;
    }
    
    if (flag == 1) {
        s << this->getTag() << "  " << strain << "  ";
        s << force << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"Truss2\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << ", " << connectedExternalOtherNodes(0) << ", " << connectedExternalOtherNodes(1) << "], ";
        s << "\"A\": " << A << ", ";
        s << "\"massperlength\": " << rho << ", ";
        s << "\"material\": \"" << theMaterial->getTag() << "\"}";
    }
}

double
	Truss2::computeCurrentStrain(void) const
{
	// NOTE method will not be called if L == 0

	// determine the strain
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();	

	double dLength = 0.0;
	for (int i = 0; i < dimension; i++)
		dLength += (disp2(i)-disp1(i))*cosX[i];

	// this method should never be called with L == 0
	return dLength/L;
}

double
	Truss2::computeCurrentNormalStrain(void) const
{
	// normal vector = (-cosX[1], cosX[0])
	if (otherLength == 0) {
		return 0;
	}

	// determine the strain
	const Vector &disp1 = theOtherNodes[0]->getTrialDisp();
	const Vector &disp2 = theOtherNodes[1]->getTrialDisp();	

	double dLength = 0.0;
	for (int i = 0; i < dimension; i++)
		dLength += (disp2(i)-disp1(i))*otherCosX[i];

	// this method should never be called with L == 0
	// double et = projFactor*dLength/otherLength;
	double et = dLength/otherLength;

	return et;
}

double
	Truss2::computeCurrentStrainRate(void) const
{
	// NOTE method will not be called if L == 0

	// determine the strain
	const Vector &vel1 = theNodes[0]->getTrialVel();
	const Vector &vel2 = theNodes[1]->getTrialVel();	

	double dLength = 0.0;
	for (int i = 0; i < dimension; i++){
		dLength += (vel2(i)-vel1(i))*cosX[i];
	}

	// this method should never be called with L == 0
	return dLength/L;
}

Response*
	Truss2::setResponse(const char **argv, int argc, OPS_Stream &output)
{

	Response *theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType","Truss2");
	output.attr("eleTag",this->getTag());
	output.attr("node1",connectedExternalNodes[0]);
	output.attr("node2",connectedExternalNodes[1]);

	//
	// we compare argv[0] for known response types for the Truss2
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

	} else if ((strcmp(argv[0],"axialForce") == 0) || 
		(strcmp(argv[0],"basicForce") == 0) || 
		(strcmp(argv[0],"localForce") == 0) || 
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
	Truss2::getResponse(int responseID, Information &eleInfo)
{
	double strain;

	switch (responseID) {
	case 1:
		return eleInfo.setVector(this->getResistingForce());

	case 2:
		return eleInfo.setDouble(A * theMaterial->getStress());

	case 3:
		if (L == 0.0) {
			strain = 0.0;
		} else {
			strain = theMaterial->getStrain();
		}
		return eleInfo.setDouble(L * strain);

	default:
		return 0;
	}
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
	Truss2::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 1)
		return -1;

	// Cross sectional area of the truss
	if (strcmp(argv[0],"A") == 0)
		return param.addObject(1, this);

	// Mass densitity of the truss
	if (strcmp(argv[0],"rho") == 0)
		return param.addObject(2, this);

	// Explicit specification of a material parameter
	if (strstr(argv[0],"material") != 0) {

		if (argc < 2)
			return -1;

		else
			return theMaterial->setParameter(&argv[1], argc-1, param);
	} 

	// Otherwise, send it to the material
	else
		return theMaterial->setParameter(argv, argc, param);
}

int
	Truss2::updateParameter (int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		A = info.theDouble;
		return 0;
	case 2:
		rho = info.theDouble;
		return 0;
	default:
		return -1;
	}
}

int
	Truss2::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}


const Matrix &
	Truss2::getKiSensitivity(int gradNumber)
{
	Matrix &stiff = *theMatrix;
	stiff.Zero();

	if (parameterID == 0) {
	}
	else if (parameterID == 1) {
		// If cross sectional area is random
		double E = theMaterial->getInitialTangent();

		int numDOF2 = numDOF/2;
		double temp;
		double EoverL = E/L;
		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				temp = cosX[i]*cosX[j]*EoverL;
				stiff(i,j) = temp;
				stiff(i+numDOF2,j) = -temp;
				stiff(i,j+numDOF2) = -temp;
				stiff(i+numDOF2,j+numDOF2) = temp;
			}
		}
	}
	else if (parameterID == 2) {
		// Nothing here when 'rho' is random
	}
	else {
		double Esens = theMaterial->getInitialTangentSensitivity(gradNumber);

		int numDOF2 = numDOF/2;
		double temp;
		double EAoverL = Esens*A/L;
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
	Truss2::getMassSensitivity(int gradNumber)
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
	Truss2::getResistingForceSensitivity(int gradNumber)
{
	theVector->Zero();

	// Initial declarations
	int i;
	double stressSensitivity, temp1, temp2;

	// Make sure the material is up to date
	double strain = this->computeCurrentStrain();
	double rate = this->computeCurrentStrainRate();
	theMaterial->setTrialStrain(strain,rate);

	// Contribution from material
	stressSensitivity = theMaterial->getStressSensitivity(gradNumber,true);

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

		double materialTangent = theMaterial->getTangent();
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
		stressSensitivity += materialTangent * strainSensitivity;
	}


	// Compute sensitivity depending on 'parameter'
	double stress = theMaterial->getStress();
	int numDOF2 = numDOF/2;
	double temp;
	if (parameterID == 1) {			// Cross-sectional area
		for (i = 0; i < dimension; i++) {
			temp = (stress + A*stressSensitivity)*cosX[i];
			(*theVector)(i) = -temp;
			(*theVector)(i+numDOF2) = temp;
		}
	}
	else {		// Density, material parameter or nodal coordinate
		for (i = 0; i < dimension; i++) {
			temp = A*(stressSensitivity*cosX[i] + stress*dcosXdh[i]);
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
	Truss2::commitSensitivity(int gradNumber, int numGrads)
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
		sens1 = theNodes[0]->getDispSensitivity(i+1, gradNumber);
		sens2 = theNodes[1]->getDispSensitivity(i+1, gradNumber);
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
	theMaterial->commitSensitivity(strainSensitivity, gradNumber, numGrads);

	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////
