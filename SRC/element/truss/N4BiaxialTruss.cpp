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

#include <N4BiaxialTruss.h>
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
#include <CompositeResponse.h>
#include <ConcretewBeta.h>

//#include <fstream>

// initialise the class wide variables
Matrix N4BiaxialTruss::trussM2(2,2);
Matrix N4BiaxialTruss::trussM8(8,8);
Matrix N4BiaxialTruss::trussM12(12,12);
Matrix N4BiaxialTruss::trussM24(24,24);
Vector N4BiaxialTruss::trussV2(2);
Vector N4BiaxialTruss::trussV4(4);
Vector N4BiaxialTruss::trussV6(6);
Vector N4BiaxialTruss::trussV8(8);
Vector N4BiaxialTruss::trussV12(12);
Vector N4BiaxialTruss::trussV24(24);

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the truss end nodes.

#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_N4BiaxialTruss)
{
	Element *theElement = 0;

	int numRemainingArgs = OPS_GetNumRemainingInputArgs();

	if (numRemainingArgs < 7) {
		opserr << "Invalid Args want: element N4BiaxialTruss $tag $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> <-doRayleigh $flag>\n";
		return 0;	
	}

	int    iData[5]; //tag, iNode, jNode, iGNode, jGNode
	double A = 0.0;
	double rho = 0.0;
	int matTag1 = 0;
	int matTag2 = 0;
	int doRayleigh = 0;
	int ndm = OPS_GetNDM();

	int numData = 5;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid integer (tag, iNode, jNode, iGNode, jGNode) in element N4BiaxialTruss " << endln;
		return 0;
	}

	numData = 1;
	if (OPS_GetDouble(&numData, &A) != 0) {
		opserr << "WARNING: Invalid A: element N4BiaxialTruss " << iData[0] << 
		" $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> <-doRayleigh $flag>\n";
		return 0;	
	}

	numData = 1;
	if (OPS_GetInt(&numData, &matTag1) != 0) {
		opserr << "WARNING: Invalid matTag1: element N4BiaxialTruss " << iData[0] << 
		" $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> <-doRayleigh $flag>\n";
		return 0;
	}

	UniaxialMaterial *theUniaxialMaterial_1 = OPS_GetUniaxialMaterial(matTag1);
	if (theUniaxialMaterial_1 == 0) {
		opserr << "WARNING: Invalid material not found element N4BiaxialTruss " << iData[0] << " $mattag1: " << matTag1 << " \n";
		return 0;
	}

	numRemainingArgs -= 6;
	while (numRemainingArgs > 1) {
	  const char *argvS = OPS_GetString();

		if (strcmp(argvS,"-rho") == 0) {
			numData = 1;
			if (OPS_GetDouble(&numData, &rho) != 0) {
				opserr << "WARNING Invalid rho in element N4BiaxialTruss " << iData[0] << 
				" $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> <-doRayleigh $flag>\n";
				return 0;
			}
		} else if (strcmp(argvS,"-doRayleigh") == 0) {
			numData = 1;
			if (OPS_GetInt(&numData, &doRayleigh) != 0) {
				opserr << "WARNING: Invalid doRayleigh in element N4BiaxialTruss " << iData[0] << 
				" $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> <-doRayleigh $flag>\n";
				return 0;
			}
		} else {
			opserr << "WARNING: Invalid option " << argvS << "  in: element N4BiaxialTruss " << iData[0] << 
			" $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> <-doRayleigh $flag>\n";
			return 0;
		}      
		numRemainingArgs -= 2;
	}

	//now create the ReinforcedConcretePlaneStress
	theElement = new N4BiaxialTruss(iData[0], ndm, iData[1], iData[2], iData[3], iData[4], *theUniaxialMaterial_1, A, rho, doRayleigh);

	if (theElement == 0) {
		opserr << "WARNING: out of memory: element N4BiaxialTruss " << iData[0] << 
		" $i1Node $j1Node $iG2Node $j2Node $A $matTag1 <-rho $rho> <-doRayleigh $flag>\n";
	}

	return theElement;
}

N4BiaxialTruss::N4BiaxialTruss(int tag, 
int dim,
int Nd1, int Nd2, 
int GNd1, int GNd2, 
UniaxialMaterial &theMat,
double a, double r, int damp)
:Element(tag,ELE_TAG_N4BiaxialTruss),     
theMaterial_1(0), theBetaMaterial_1(0), 
theMaterial_2(0), theBetaMaterial_2(0), 
connectedExternalNodes(4),
dimension(dim), numDOF(0), theLoad(0),
theMatrix(0), theVector(0), theVector2(0),
L(0.0), A(a), rho(r), doRayleighDamping(damp)
{
	// get a copy of the material and check we obtained a valid copy
	theMaterial_1 = theMat.getCopy();
	theMaterial_2 = theMat.getCopy();
	if ((theMaterial_1 == 0) || (theMaterial_2 == 0)) {
		opserr << "FATAL N4BiaxialTruss::N4BiaxialTruss - " << tag <<
		"failed to get a copy of material with tag " << theMat.getTag() << endln;
		exit(-1);
	} else if (theMat.getClassTag() == MAT_TAG_ConcretewBeta) {
	  theBetaMaterial_1 = (ConcretewBeta *) theMaterial_1;
	  theBetaMaterial_2 = (ConcretewBeta *) theMaterial_2;
	}
	
	// ensure the connectedExternalNode ID is of correct size & set values
	if (connectedExternalNodes.Size() != 4) {
		opserr << "FATAL N4BiaxialTruss::N4BiaxialTruss - " <<  tag << "failed to create an node ID array of size 4\n";
		exit(-1);
	}

	connectedExternalNodes(0) = Nd1;
	connectedExternalNodes(1) = Nd2;
	connectedExternalNodes(2) = GNd1;
	connectedExternalNodes(3) = GNd2;

	// set node pointers to NULL
	for (int i=0; i<4; i++)
	theNodes[i] = 0;

	cosX[0] = 0.0;
	cosX[1] = 0.0;
	cosX[2] = 0.0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
N4BiaxialTruss::N4BiaxialTruss()
:Element(0,ELE_TAG_N4BiaxialTruss),     
theMaterial_1(0), theBetaMaterial_1(0), 
theMaterial_2(0), theBetaMaterial_2(0), 
connectedExternalNodes(2),
dimension(0), numDOF(0), theLoad(0),
theMatrix(0), theVector(0), theVector2(0),
L(0.0), A(0.0), rho(0.0)
{
	// ensure the connectedExternalNode ID is of correct size 
	if (connectedExternalNodes.Size() != 4) {
		opserr << "FATAL N4BiaxialTruss::N4BiaxialTruss - failed to create an ID of size 2\n";
		exit(-1);
	}

	for (int i=0; i<4; i++)
	theNodes[i] = 0;

	cosX[0] = 0.0;
	cosX[1] = 0.0;
	cosX[2] = 0.0;
}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
N4BiaxialTruss::~N4BiaxialTruss()
{
	// invoke the destructor on any objects created by the object
	// that the object still holds a pointer to
	if (theMaterial_1 != 0)
	delete theMaterial_1;
	if (theMaterial_2 != 0)
	delete theMaterial_2;
	if (theLoad != 0)
	delete theLoad;
}


int
N4BiaxialTruss::getNumExternalNodes(void) const
{
	return 4;
}

const ID &
N4BiaxialTruss::getExternalNodes(void) 
{
	return connectedExternalNodes;
}

Node **
N4BiaxialTruss::getNodePtrs(void) 
{
	return theNodes;
}

int
N4BiaxialTruss::getNumDOF(void) 
{
	return numDOF;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the N4BiaxialTruss element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
N4BiaxialTruss::setDomain(Domain *theDomain)
{
	// check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		L = 0;
		return;
	}

	// first set the node pointers
	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);
	int GNd1 = connectedExternalNodes(2);
	int GNd2 = connectedExternalNodes(3);
	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);	
	theNodes[2] = theDomain->getNode(GNd1);	
	theNodes[3] = theDomain->getNode(GNd2);	
	
	// if can't find the nodes - send a warning message
	if ((theNodes[0] == 0)) {
		opserr <<"N4BiaxialTruss::setDomain() - N4BiaxialTruss" << this->getTag() << " node " << Nd1 <<
		"does not exist in the model\n";
		
		// fill this in so don't segment fault later
		numDOF = 2;    
		theMatrix = &trussM2;
		theVector = &trussV2;	
		theVector2 = &trussV2;	
		return;
	} else if ((theNodes[1] == 0)) {
		opserr <<"N4BiaxialTruss::setDomain() - N4BiaxialTruss" << this->getTag() << " node " << Nd2 <<
		"does not exist in the model\n";
		
		// fill this in so don't segment fault later
		numDOF = 2;    
		theMatrix = &trussM2;
		theVector = &trussV2;	
		theVector2 = &trussV2;	
		return;
	} else if ((theNodes[2] == 0)) {
		opserr <<"N4BiaxialTruss::setDomain() - N4BiaxialTruss" << this->getTag() << " node " << GNd1 <<
		"does not exist in the model\n";
		
		// fill this in so don't segment fault later
		numDOF = 2;    
		theMatrix = &trussM2;
		theVector = &trussV2;	
		theVector2 = &trussV2;	
		return;
	} else if ((theNodes[3] == 0)) {
		opserr <<"N4BiaxialTruss::setDomain() - N4BiaxialTruss" << this->getTag() << " node " << GNd2 <<
		"does not exist in the model\n";
		
		// fill this in so don't segment fault later
		numDOF = 2;    
		theMatrix = &trussM2;
		theVector = &trussV2;	
		theVector2 = &trussV2;	
		return;
	}

	// now determine the number of dof and the dimesnion    
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();	
	int dofNd3 = theNodes[2]->getNumberDOF();	
	int dofNd4 = theNodes[3]->getNumberDOF();	

	// if differing dof at the ends - print a warning message
	if (dofNd1 != dofNd2 || dofNd2 != dofNd3 || dofNd3 != dofNd4 || dofNd4 != dofNd1) {
		opserr <<"WARNING N4BiaxialTruss::setDomain(): nodes have differing dof at ends for N4BiaxialTruss " << this->getTag() << endln;

		// fill this in so don't segment fault later
		numDOF = 2;    
		theMatrix = &trussM2;
		theVector = &trussV2;	
		theVector2 = &trussV2;	
		
		return;
	}	

	// call the base class method
	this->DomainComponent::setDomain(theDomain);

	// now set the number of dof for element and set matrix and vector pointer
	if (dimension == 2 && dofNd1 == 2) {
		numDOF = 8;
		theMatrix = &trussM8;
		theVector = &trussV8;	
		theVector2 = &trussV4;	
	}
	else if (dimension == 2 && dofNd1 == 3) {
		numDOF = 12;	
		theMatrix = &trussM12;
		theVector = &trussV12;		
		theVector2 = &trussV6;
	}
	else if (dimension == 3 && dofNd1 == 3) {
		numDOF = 12;	
		theMatrix = &trussM12;
		theVector = &trussV12;			
		theVector2 = &trussV6;
	}
	else if (dimension == 3 && dofNd1 == 6) {
		numDOF = 24;	    
		theMatrix = &trussM24;
		theVector = &trussV24;			
		theVector2 = &trussV12;	
	}
	else {
		opserr <<"WARNING N4BiaxialTruss::setDomain cannot handle " << dimension << " dofs at nodes in " << 
		dofNd1  << " problem\n";

		numDOF = 2;    
		theMatrix = &trussM2;
		theVector = &trussV2;	
		theVector2 = &trussV2;	
		return;
	}

	if (theLoad == 0)
	theLoad = new Vector(numDOF);
	else if (theLoad->Size() != numDOF) {
		delete theLoad;
		theLoad = new Vector(numDOF);
	}

	if (theLoad == 0) {
		opserr << "N4BiaxialTruss::setDomain - N4BiaxialTruss " << this->getTag() << 
		"out of memory creating vector of size" << numDOF << endln;
		exit(-1);
		return;
	}          
	
	// now determine the length, directions, and fill in transformation matrix
	const Vector &end1Crd = theNodes[0]->getCrds();
	const Vector &end2Crd = theNodes[1]->getCrds();	
	const Vector &GN1Crd = theNodes[2]->getCrds();
	const Vector &GN2Crd = theNodes[3]->getCrds();	
	
	if (dimension == 2) {
		double dx = end2Crd(0) - end1Crd(0);
		double dy = end2Crd(1) - end1Crd(1);
		
		double dx2 = GN2Crd(0) - GN1Crd(0);
		double dy2 = GN2Crd(1) - GN1Crd(1);
		L = sqrt(dx*dx + dy*dy);
		L2 = sqrt(dx2*dx2 + dy2*dy2);
		
		if (L == 0.0 || L2 == 0.0) {
			opserr <<"WARNING N4BiaxialTruss::setDomain() - N4BiaxialTruss " << this->getTag() << " has zero length\n";
			return;
		} else if (fabs(L - L2) > 1.e-6) {
			opserr <<"WARNING N4BiaxialTruss::setDomain() - N4BiaxialTruss " << this->getTag() << " is not a rectangular element, error in diagonal length comparison\n";
			return;
		}
		cosX[0] = dx/L;
		cosX[1] = dy/L;
		cosX[2] = 0;
		cosX2[0] = dx2/L2;
		cosX2[1] = dy2/L2;
		cosX2[2] = 0;
		
		// forming local coordinates
		vectorX[0] = GN1Crd(0)-end1Crd(0);
		vectorX[1] = GN1Crd(1)-end1Crd(1);
		vectorX[2] = 0;
		
		vectorY[0] = GN2Crd(0)-end1Crd(0);
		vectorY[1] = GN2Crd(1)-end1Crd(1);
		vectorY[2] = 0;
	
		lengthX = sqrt(vectorX[0]*vectorX[0] + vectorX[1]*vectorX[1]);
		lengthY = sqrt(vectorY[0]*vectorY[0] + vectorY[1]*vectorY[1]);
		
		if (lengthX == 0.0 || lengthY == 0.0) {
			opserr <<"WARNING N4BiaxialTruss::setDomain() - N4BiaxialTruss " << this->getTag() << " is not rectangular, error in forming local coordinates\n";
			return;
		} 
		
		// normalizing direction vectors
		vectorX[0] = vectorX[0]/lengthX;
		vectorX[1] = vectorX[1]/lengthX;
		
		vectorY[0] = vectorY[0]/lengthY;
		vectorY[1] = vectorY[1]/lengthY;
		
		double perpCheck = vectorX[0]*vectorY[0] + vectorX[1]*vectorY[1];
		if (fabs(perpCheck) > 1.e-6) {
			opserr <<"WARNING N4BiaxialTruss::setDomain() - N4BiaxialTruss " << this->getTag() << " is not a rectangular element\n";
			opserr << "perpCheck returns " << perpCheck << "\n";
			return;
		}
	}
	else {
		double dx = end2Crd(0) - end1Crd(0);
		double dy = end2Crd(1) - end1Crd(1);
		double dz = end2Crd(2) - end1Crd(2);
		
		double dx2 = GN2Crd(0) - GN1Crd(0);
		double dy2 = GN2Crd(1) - GN1Crd(1);
		double dz2 = GN2Crd(2) - GN1Crd(2);
		L = sqrt(dx*dx + dy*dy + dz*dz);
		L2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
		
		if (L == 0.0 || L2 == 0.0) {
			opserr <<"WARNING N4BiaxialTruss::setDomain() - N4BiaxialTruss " << this->getTag() << " has zero length\n";
			return;
		} else if (fabs(L - L2) > 1.e-6) {
			opserr <<"WARNING N4BiaxialTruss::setDomain() - N4BiaxialTruss " << this->getTag() << " is not a rectangular element, error in diagonal length comparison\n";
			return;
		}
		
		cosX[0] = dx/L;
		cosX[1] = dy/L;
		cosX[2] = dz/L;
		cosX2[0] = dx2/L2;
		cosX2[1] = dy2/L2;
		cosX2[2] = dz2/L2;
		
		// forming local coordinates
		vectorX[0] = GN1Crd(0)-end1Crd(0);
		vectorX[1] = GN1Crd(1)-end1Crd(1);
		vectorX[2] = GN1Crd(2)-end1Crd(2);
		
		vectorY[0] = GN2Crd(0)-end1Crd(0);
		vectorY[1] = GN2Crd(1)-end1Crd(1);
		vectorY[2] = GN2Crd(2)-end1Crd(2);

		lengthX = sqrt(vectorX[0]*vectorX[0] + vectorX[1]*vectorX[1] + vectorX[2]*vectorX[2]);
		lengthY = sqrt(vectorY[0]*vectorY[0] + vectorY[1]*vectorY[1] + vectorY[2]*vectorY[2]);
		
		if (lengthX == 0.0 || lengthY == 0.0) {
			opserr <<"WARNING N4BiaxialTruss::setDomain() - N4BiaxialTruss " << this->getTag() << " is not rectangular, error in forming local coordinates\n";
			return;
		} 
				
		// normalizing direction vectors
		vectorX[0] = vectorX[0]/lengthX;
		vectorX[1] = vectorX[1]/lengthX;
		vectorX[2] = vectorX[2]/lengthX;
		
		vectorY[0] = vectorY[0]/lengthY;
		vectorY[1] = vectorY[1]/lengthY;
		vectorY[2] = vectorY[2]/lengthY;
		
		// checking the rectangularity of the elements
		double dx0 = GN2Crd(0)-end2Crd(0);
		double dx1 = GN2Crd(1)-end2Crd(1);
		dx2 = GN2Crd(2)-end2Crd(2);
		
		double dy0 = GN1Crd(0)-end2Crd(0);
		double dy1 = GN1Crd(1)-end2Crd(1);
		dy2 = GN1Crd(2)-end2Crd(2);
		
		double Lx2 = sqrt(dx0*dx0 + dx1*dx1 +dx2*dx2);
		double Ly2 = sqrt(dy0*dy0 + dy1*dy1 +dy2*dy2);
		
		double perpCheck = vectorX[0]*vectorY[0] + vectorX[1]*vectorY[1] + vectorX[2]*vectorY[2];
		if (fabs(perpCheck) > 1.e-6 || fabs(Lx2 - lengthX) > 1.e-6 || fabs(Ly2 - lengthY) > 1.e-6) {
			opserr <<"WARNING N4BiaxialTruss::setDomain() - N4BiaxialTruss " << this->getTag() << " is not a rectangular element\n";
			return;
		}
	}
	
	// fill in computation variables:
	oneOverL = 1.0/L;
	LxoverL = lengthX/L;
	LyoverL = lengthY/L;
	oneOver2Lx = 1.0/(2.0*lengthX);
	oneOver2Ly = 1.0/(2.0*lengthY);
}

int
N4BiaxialTruss::commitState()
{
	int retVal = 0;
	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "N4BiaxialTruss::commitState () - failed in base class";
	}    
	retVal = theMaterial_1->commitState();
	retVal += theMaterial_2->commitState();
	return retVal;
}

int
N4BiaxialTruss::revertToLastCommit()
{
	int retVal = theMaterial_1->revertToLastCommit();
	retVal += theMaterial_2->revertToLastCommit();
	return retVal;
}

int
N4BiaxialTruss::revertToStart()
{
	int retVal = theMaterial_1->revertToStart();
	retVal += theMaterial_2->revertToStart();
	return retVal;
}

int
N4BiaxialTruss::update(void)
{
	// determine the current strain given trial displacements at nodes
	int retVal = this->computeCurrentStrainBiaxial();
	retVal += this->computeCurrentStrainRate();
	
	if (theBetaMaterial_1) {
		retVal += theBetaMaterial_1->setTrialStrainwBeta(strain_1, normalStrain_1, strainRate_1);
	} else {
		retVal += theMaterial_1->setTrialStrain(strain_1, strainRate_1);
	}
	
	if (theBetaMaterial_2) {
		return theBetaMaterial_2->setTrialStrainwBeta(strain_2, normalStrain_2, strainRate_2);
	} else {
		return theMaterial_2->setTrialStrain(strain_2, strainRate_2);
	}
}


const Matrix &
N4BiaxialTruss::getTangentStiff(void)
{
	if (L == 0.0) { // - problem in setDomain() no further warnings
		theMatrix->Zero();
		return *theMatrix;
	}
	
	double E1 = theMaterial_1->getTangent();
	double E2 = theMaterial_2->getTangent();

	// come back later and redo this if too slow
	Matrix &stiff = *theMatrix;
	stiff.Zero();

	int numDOF2 = numDOF/4;
	double temp, temp2;
	double EAoverL = E1*A*oneOverL;
	double EAoverL2 = E2*A*oneOverL;
	
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			temp = cosX[i]*cosX[j]*EAoverL;
			temp2 = cosX2[i]*cosX2[j]*EAoverL2;
			stiff(i,		j)		   = temp;
			stiff(i+numDOF2,j)	       = -temp;
			stiff(i,		j+numDOF2) = -temp;
			stiff(i+numDOF2,j+numDOF2) = temp;
			stiff(i+2*numDOF2,j+2*numDOF2) = temp2;
			stiff(i+3*numDOF2,j+2*numDOF2) = -temp2;
			stiff(i+2*numDOF2,j+3*numDOF2) = -temp2;
			stiff(i+3*numDOF2,j+3*numDOF2) = temp2;
		}
	}
	
	return stiff;
}


const Matrix &
N4BiaxialTruss::getInitialStiff(void)
{
	if (L == 0.0) { // - problem in setDomain() no further warnings
		return *theMatrix;
		theMatrix->Zero();
	}
	
	double E1 = theMaterial_1->getInitialTangent();
	double E2 = theMaterial_2->getInitialTangent();

	// come back later and redo this if too slow
	Matrix &stiff = *theMatrix;
	stiff.Zero();    

	int numDOF2 = numDOF/4;
	double temp, temp2;
	double EAoverL = E1*A*oneOverL;
	double EAoverL2 = E2*A*oneOverL;
	
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			temp = cosX[i]*cosX[j]*EAoverL;
			temp2 = cosX2[i]*cosX2[j]*EAoverL2;
			stiff(i,		j)		   = temp;
			stiff(i+numDOF2,j)	       = -temp;
			stiff(i,		j+numDOF2) = -temp;
			stiff(i+numDOF2,j+numDOF2) = temp;
			stiff(i+2*numDOF2,j+2*numDOF2) = temp2;
			stiff(i+3*numDOF2,j+2*numDOF2) = -temp2;
			stiff(i+2*numDOF2,j+3*numDOF2) = -temp2;
			stiff(i+3*numDOF2,j+3*numDOF2) = temp2;
		}
	}
	
	return stiff;
}

const Matrix &
N4BiaxialTruss::getDamp(void)
{
	if (L == 0.0) { // - problem in setDomain() no further warnings
		theMatrix->Zero();
		return *theMatrix;
	}
	theMatrix->Zero();

	if (doRayleighDamping == 1) {
		*theMatrix = this->Element::getDamp();
	}

	double eta = theMaterial_1->getDampTangent();
	double eta2 = theMaterial_2->getDampTangent();

	// come back later and redo this if too slow
	Matrix &damp = *theMatrix;

	int numDOF2 = numDOF/4;
	double temp, temp2;
	double etaAoverL = eta*A*oneOverL;
	double etaAoverL2 = eta2*A*oneOverL;
	
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			temp = cosX[i]*cosX[j]*etaAoverL;
			temp2 = cosX[i]*cosX[j]*etaAoverL2;
			damp(i,j) += temp;
			damp(i+numDOF2,j) += -temp;
			damp(i,j+numDOF2) += -temp;
			damp(i+numDOF2,j+numDOF2) += temp;
			damp(i+2*numDOF2,j+2*numDOF2) += temp2;
			damp(i+3*numDOF2,j+2*numDOF2) += -temp2;
			damp(i+2*numDOF2,j+3*numDOF2) += -temp2;
			damp(i+3*numDOF2,j+3*numDOF2) += temp2;
		}
	}
	
	return damp;
}


const Matrix &
N4BiaxialTruss::getMass(void)
{   
	// zero the matrix
	Matrix &mass = *theMatrix;
	mass.Zero();    

	// check for quick return
	if (L == 0.0 || rho == 0.0) { // - problem in setDomain() no further warnings
		return mass;
	}

	double M = 0.5*rho*L;
	int numDOF2 = numDOF/4;
	for (int i = 0; i < dimension; i++) {
		mass(i,i) = M;
		mass(i+numDOF2,i+numDOF2) = M;
		mass(i+2*numDOF2,i+2*numDOF2) = M;
		mass(i+3*numDOF2,i+3*numDOF2) = M;
	}

	return mass;
}

void 
N4BiaxialTruss::zeroLoad(void)
{
	theLoad->Zero();
}

int 
N4BiaxialTruss::addLoad(ElementalLoad *theLoad, double loadFactor)
{  
	opserr <<"N4BiaxialTruss::addLoad - load type unknown for N4BiaxialTruss with tag: " << this->getTag() << endln; 
	return -1;
}

int 
N4BiaxialTruss::addInertiaLoadToUnbalance(const Vector &accel)
{
	// check for a quick return
	if (L == 0.0 || rho == 0.0) 
	return 0;

	// get R * accel from the nodes
	const Vector &Raccel1 = theNodes[0]->getRV(accel);
	const Vector &Raccel2 = theNodes[1]->getRV(accel);    
	const Vector &Raccel3 = theNodes[2]->getRV(accel);   
	const Vector &Raccel4 = theNodes[3]->getRV(accel);   

	int nodalDOF = numDOF/4;
	
	double M = 0.25*rho*L;
	// want to add ( - fact * M R * accel ) to unbalance
	for (int i=0; i<dimension; i++) {
		double val1 = Raccel1(i);
		double val2 = Raccel2(i);
		double val3 = Raccel3(i);	
		double val4 = Raccel4(i);		

		// perform - fact * M*(R * accel) // remember M a diagonal matrix
		val1 *= -M;
		val2 *= -M;
		val3 *= -M;
		val4 *= -M;
		
		(*theLoad)(i) += val1;
		(*theLoad)(i+  nodalDOF) += val2;
		(*theLoad)(i+2*nodalDOF) += val3;
		(*theLoad)(i+3*nodalDOF) += val4;
	}	

	return 0;
}

const Vector &
N4BiaxialTruss::getResistingForce()
{	
	if (L == 0.0) { // - problem in setDomain() no further warnings
		theVector->Zero();
		return *theVector;
	}
	
	// R = Ku - Pext
	// Ku = F * transformation
	double force1 = A*theMaterial_1->getStress();
	double force2 = A*theMaterial_2->getStress();
	int numDOF2 = numDOF/4;
	double temp;
	for (int i = 0; i < dimension; i++) {
		temp = cosX[i]*force1;
		(*theVector)(i) = -temp;
		(*theVector)(i+numDOF2) = temp;
		temp = cosX2[i]*force2;
		(*theVector)(i+2*numDOF2) = -temp;
		(*theVector)(i+3*numDOF2) = temp;
	}

	// subtract external load:  Ku - P
	(*theVector) -= *theLoad;

	return *theVector;
}


const Vector &
N4BiaxialTruss::getResistingForceIncInertia()
{	
	this->getResistingForce();

	// now include the mass portion
	if (L != 0.0 && rho != 0.0) {
		
		const Vector &accel1 = theNodes[0]->getTrialAccel();
		const Vector &accel2 = theNodes[1]->getTrialAccel();	
		const Vector &accel3 = theNodes[2]->getTrialAccel();	
		const Vector &accel4 = theNodes[3]->getTrialAccel();	
		
		int numDOF2 = numDOF/4;
		double M = 0.5*rho*L;
		for (int i = 0; i < dimension; i++) {
			(*theVector)(i) += M*accel1(i);
			(*theVector)(i+numDOF2) += M*accel2(i);
			(*theVector)(i+2*numDOF2) += M*accel3(i);
			(*theVector)(i+3*numDOF2) += M*accel4(i);
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

int
N4BiaxialTruss::sendSelf(int commitTag, Channel &theChannel)
{
	int res;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// N4BiaxialTruss packs it's data into a Vector and sends this to theChannel
	// along with it's dbTag and the commitTag passed in the arguments

	static Vector data(11);
	data(0) = this->getTag();
	data(1) = dimension;
	data(2) = numDOF;
	data(3) = A;
	data(4) = theMaterial_1->getClassTag();
	data(5) = theMaterial_2->getClassTag();
	data(6) = rho;

	if (doRayleighDamping == 0)
	data(7) = 0;
	else
	data(7) = 1;
	
	int matDbTag1 = theMaterial_1->getDbTag();
	int matDbTag2 = theMaterial_2->getDbTag();

	// NOTE: we do have to ensure that the material has a database
	// tag if we are sending to a database channel.
	if (matDbTag1 == 0) {
		matDbTag1 = theChannel.getDbTag();
		if (matDbTag1 != 0)
		theMaterial_1->setDbTag(matDbTag1);
	}	
	if (matDbTag2 == 0) {
		matDbTag2 = theChannel.getDbTag();
		if (matDbTag2 != 0)
		theMaterial_1->setDbTag(matDbTag2);
	}
	data(8) = matDbTag1;
	data(9) = matDbTag2;

	res = theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr <<"WARNING N4BiaxialTruss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return -1;
	}

	// N4BiaxialTruss then sends the tags of it's nodes
	res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr <<"WARNING N4BiaxialTruss::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return -2;
	}

	// finally N4BiaxialTruss asks it's material object to send itself
	res = theMaterial_1->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr <<"WARNING N4BiaxialTruss::sendSelf() - " << this->getTag() << " failed to send its Material_1\n";
		return -3;
	}
	res = theMaterial_2->sendSelf(commitTag, theChannel);
	if (res < 0) {
		opserr <<"WARNING N4BiaxialTruss::sendSelf() - " << this->getTag() << " failed to send its Material_2\n";
		return -3;
	}

	return 0;
}

int
N4BiaxialTruss::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res;
	int dataTag = this->getDbTag();

	// N4BiaxialTruss creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(11);
	res = theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr <<"WARNING N4BiaxialTruss::recvSelf() - failed to receive Vector\n";
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

	// N4BiaxialTruss now receives the tags of it's nodes
	res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr <<"WARNING N4BiaxialTruss::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return -2;
	}

	// finally N4BiaxialTruss creates a material object of the correct type,
	// sets its database tag and asks this new object to recveive itself.
	int matClass1 = (int)data(4);
	int matDb1 = (int)data(8);
	int matClass2 = (int)data(5);
	int matDb2 = (int)data(9);

	// check if we have a material object already & if we do if of right type
	if ((theMaterial_1 == 0) || (theMaterial_1->getClassTag() != matClass1)) {

		// if old one .. delete it
		if (theMaterial_1 != 0)
		delete theMaterial_1;

		// create a new material object
		theMaterial_1 = theBroker.getNewUniaxialMaterial(matClass1);
		if (theMaterial_1 == 0) {
			opserr <<"WARNING N4BiaxialTruss::recvSelf() - " << this->getTag() 
			<< " failed to get a blank Material of type " << matClass1 << endln;
			return -3;
		}
		
		if (theMaterial_1->getClassTag() == MAT_TAG_ConcretewBeta) {
			theBetaMaterial_1 = (ConcretewBeta *) theMaterial_1;
		}
		
	}
	theMaterial_1->setDbTag(matDb1); // note: we set the dbTag before we receive the material
	res = theMaterial_1->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		opserr <<"WARNING N4BiaxialTruss::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
		return -3;    
	}
	
	// Again with Material_2
	if ((theMaterial_2 == 0) || (theMaterial_2->getClassTag() != matClass2)) {

		// if old one .. delete it
		if (theMaterial_2 != 0)
		delete theMaterial_2;

		// create a new material object
		theMaterial_2 = theBroker.getNewUniaxialMaterial(matClass2);
		if (theMaterial_2 == 0) {
			opserr <<"WARNING N4BiaxialTruss::recvSelf() - " << this->getTag() 
			<< " failed to get a blank Material of type " << matClass2 << endln;
			return -3;
		}
		
		if (theMaterial_2->getClassTag() == MAT_TAG_ConcretewBeta) {
			theBetaMaterial_2 = (ConcretewBeta *) theMaterial_2;
		}
		
	}
	theMaterial_2->setDbTag(matDb2); // note: we set the dbTag before we receive the material
	res = theMaterial_2->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		opserr <<"WARNING N4BiaxialTruss::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
		return -3;    
	}

	return 0;
}


int
N4BiaxialTruss::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// ensure setDomain() worked
	if (L == 0.0)
	return 0;

	// get display coordinates
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);
	theNodes[2]->getDisplayCrds(v3, fact, displayMode);
	theNodes[3]->getDisplayCrds(v4, fact, displayMode);

	// determine color and draw lines
	int retVal = 0;
	if (displayMode == 1 || displayMode == 2) {
		// compute the strain and axial force in the member
		double force1, force2;
		if (L == 0.0) {
			strain_1 = 0.0;
			strain_2 = 0.0;
			force1 = 0.0;
			force2 = 0.0;
		}
		else {
			this->computeCurrentStrainBiaxial();
			theMaterial_1->setTrialStrain(strain_1);
			theMaterial_2->setTrialStrain(strain_2);
			force1 = A * theMaterial_1->getStress();
			force2 = A * theMaterial_2->getStress();
		}
		if (displayMode == 2) {// use the strain as the drawing measure
			retVal += theViewer.drawLine(v1, v2, (float)strain_1, (float)strain_1);
			retVal += theViewer.drawLine(v3, v4, (float)strain_2, (float)strain_2);
		}
		else { // otherwise use the axial force as measure
			retVal += theViewer.drawLine(v1, v2, (float)force1, (float)force1);
			retVal += theViewer.drawLine(v3, v4, (float)force2, (float)force2);
		}

	}
	else {
		retVal += theViewer.drawLine(v1, v2, 1.0, 1.0);
		retVal += theViewer.drawLine(v3, v4, 1.0, 1.0);
	}
	return retVal;
}

void
N4BiaxialTruss::Print(OPS_Stream &s, int flag)
{
	// compute the strain and axial force in the member
	double strain1, force1, strain2, force2;
	strain1 = theMaterial_1->getStrain();
	force1 = A * theMaterial_1->getStress();
	strain2 = theMaterial_2->getStrain();
	force2 = A * theMaterial_2->getStress();
	
	if (flag == 0) { // print everything
		s << endln;
		s << "Element: " << this->getTag(); 
		s << " type: Truss2  iNode: " << connectedExternalNodes(0);
		s << " jNode: " << connectedExternalNodes(1);
		s << " Area: " << A << " Mass/Length: " << rho;
		
		s << " \n\t strain: " << strain1;
		s << " axial load: " <<  force1;
		if (L != 0.0) {
			int numDOF2 = numDOF/4;
			double temp;
			for (int i = 0; i < dimension; i++) {
				temp = cosX[i]*force1;
				(*theVector2)(i) = -temp;
				(*theVector2)(i+numDOF2) = temp;
			}
			s << " \n\t unbalanced load: " << *theVector2;	
		}

		s << " \t Material: " << *theMaterial_1;
		s << endln;
		s << endln;
		
		s << "Element: " << this->getTag()+1; 
		s << " type: Truss2  iNode: " << connectedExternalNodes(2);
		s << " jNode: " << connectedExternalNodes(3);
		s << " Area: " << A << " Mass/Length: " << rho;
		
		s << " \n\t strain: " << strain2;
		s << " axial load: " <<  force2;
		if (L != 0.0) {
			int numDOF2 = numDOF/4;
			double temp;
			for (int i = 0; i < dimension; i++) {
				temp = cosX[i]*force1;
				(*theVector2)(i) = -temp;
				(*theVector2)(i+numDOF2) = temp;
			}
			s << " \n\t unbalanced load: " << *theVector2;	
		}

		s << " \t Material: " << *theMaterial_2;
		s << endln;
		s << endln;
	} else if (flag == 1) {
		s << this->getTag() << "  " << strain1 << "  ";
		s << force1 << endln;
		s << endln;
		s << this->getTag()+1 << "  " << strain2 << "  ";
		s << force2 << endln;
	} else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"N4BiaxialTruss\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << ", " << connectedExternalNodes(2) << ", " << connectedExternalNodes(3) << "], ";
		s << "\"A\": " << A << ", ";
		s << "\"massperlength\": " << rho << ", ";
		s << "\"material\": \"" << theMaterial_1->getTag() << "\"}";
	}
}

int
N4BiaxialTruss::computeCurrentStrainBiaxial(void)
{
	// NOTE method will not be called if Lx, Ly == 0

	// extract displacements
	const Vector &dispN1 = theNodes[0]->getTrialDisp();
	const Vector &dispN2 = theNodes[1]->getTrialDisp();	
	const Vector &dispGN1 = theNodes[2]->getTrialDisp();	
	const Vector &dispGN2 = theNodes[3]->getTrialDisp();	

	double dispXN1, dispXN2, dispXGN1, dispXGN2, dispYN1, dispYN2, dispYGN1, dispYGN2;

	// convert to local X, Y directions
	if (dimension == 2) {
		strain_1 = oneOverL*(( dispN2(0)- dispN1(0))* cosX[0] + ( dispN2(1)- dispN1(1))* cosX[1]);
		strain_2 = oneOverL*((dispGN2(0)-dispGN1(0))*cosX2[0] + (dispGN2(1)-dispGN1(1))*cosX2[1]);
	
		dispXN1  = vectorX[0]*dispN1(0)  + vectorX[1]*dispN1(1);
		dispXN2  = vectorX[0]*dispN2(0)  + vectorX[1]*dispN2(1);
		dispXGN1 = vectorX[0]*dispGN1(0) + vectorX[1]*dispGN1(1);
		dispXGN2 = vectorX[0]*dispGN2(0) + vectorX[1]*dispGN2(1);
		
		dispYN1  = vectorY[0]*dispN1(0)  + vectorY[1]*dispN1(1);
		dispYN2  = vectorY[0]*dispN2(0)  + vectorY[1]*dispN2(1);
		dispYGN1 = vectorY[0]*dispGN1(0) + vectorY[1]*dispGN1(1);
		dispYGN2 = vectorY[0]*dispGN2(0) + vectorY[1]*dispGN2(1);
		
	} else {
		strain_1 = oneOverL*(( dispN2(0)- dispN1(0))* cosX[0] + ( dispN2(1)- dispN1(1))* cosX[1] + ( dispN2(2)- dispN1(2))* cosX[2]); 
		strain_2 = oneOverL*((dispGN2(0)-dispGN1(0))*cosX2[0] + (dispGN2(1)-dispGN1(1))*cosX2[1] + (dispGN2(2)-dispGN1(2))*cosX2[2]);
		
		dispXN1  = vectorX[0]*dispN1(0)  + vectorX[1]*dispN1(1)  + vectorX[2]*dispN1(2);
		dispXN2  = vectorX[0]*dispN2(0)  + vectorX[1]*dispN2(1)  + vectorX[2]*dispN2(2);
		dispXGN1 = vectorX[0]*dispGN1(0) + vectorX[1]*dispGN1(1) + vectorX[2]*dispGN1(2);
		dispXGN2 = vectorX[0]*dispGN2(0) + vectorX[1]*dispGN2(1) + vectorX[2]*dispGN2(2);
		
		dispYN1  = vectorY[0]*dispN1(0)  + vectorY[1]*dispN1(1)  + vectorY[2]*dispN1(2);
		dispYN2  = vectorY[0]*dispN2(0)  + vectorY[1]*dispN2(1)  + vectorY[2]*dispN2(2);
		dispYGN1 = vectorY[0]*dispGN1(0) + vectorY[1]*dispGN1(1) + vectorY[2]*dispGN1(2);
		dispYGN2 = vectorY[0]*dispGN2(0) + vectorY[1]*dispGN2(1) + vectorY[2]*dispGN2(2);
	}
	
	// computing strains at mid-point
	double epsXX = -oneOver2Lx*dispXN1 + oneOver2Lx*dispXGN1 + oneOver2Lx*dispXN2 - oneOver2Lx*dispXGN2;
	
	double epsYY = -oneOver2Ly*dispYN1 - oneOver2Ly*dispYGN1 + oneOver2Lx*dispYN2 + oneOver2Ly*dispYGN2;
	
	double gammaXY = -oneOver2Ly*dispXN1 - oneOver2Lx*dispYN1 - oneOver2Ly*dispXGN1 + oneOver2Lx*dispYGN1 
	                + oneOver2Ly*dispXN2 + oneOver2Lx*dispYN2 + oneOver2Ly*dispXGN2 - oneOver2Lx*dispYGN2;
	
	normalStrain_1 = LyoverL*LyoverL*epsXX - LxoverL*LyoverL*gammaXY + LxoverL*LxoverL*epsYY;
	normalStrain_2 = LyoverL*LyoverL*epsXX + LxoverL*LyoverL*gammaXY + LxoverL*LxoverL*epsYY;

	return 0;
}

int
N4BiaxialTruss::computeCurrentStrainRate(void)
{
	// NOTE method will not be called if L == 0

	// determine the velocity
	const Vector &velN1 = theNodes[0]->getTrialVel();
	const Vector &velN2 = theNodes[1]->getTrialVel();	
	const Vector &velGN1 = theNodes[2]->getTrialVel();	
	const Vector &velGN2 = theNodes[3]->getTrialVel();	

	if (dimension == 2) {
		strainRate_1 = oneOverL*(( velN2(0)- velN1(0))* cosX[0] + ( velN2(1)- velN1(1))* cosX[1]);
		strainRate_2 = oneOverL*((velGN2(0)-velGN1(0))*cosX2[0] + (velGN2(1)-velGN1(1))*cosX2[1]);
	
	} else {
		strainRate_1 = oneOverL*(( velN2(0)- velN1(0))* cosX[0] + ( velN2(1)- velN1(1))* cosX[1] + ( velN2(2)- velN1(2))* cosX[2]); 
		strainRate_2 = oneOverL*((velGN2(0)-velGN1(0))*cosX2[0] + (velGN2(1)-velGN1(1))*cosX2[1] + (velGN2(2)-velGN1(2))*cosX2[2]);
	}
	return 0;
}

Response*
N4BiaxialTruss::setResponse(const char **argv, int argc, OPS_Stream &output)
{

	Response *theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType","N4BiaxialTruss");
	output.attr("eleTag",this->getTag());
	output.attr("truss1_node1",connectedExternalNodes[0]);
	output.attr("truss1_node2",connectedExternalNodes[1]);
	output.attr("truss2_node1",connectedExternalNodes[2]);
	output.attr("truss2_node2",connectedExternalNodes[3]);

	//
	// we compare argv[0] for known response types for the N4BiaxialTruss
	//

	if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
			|| (strcmp(argv[0],"globalForce") == 0) || (strcmp(argv[0],"globalForces") == 0)){
		char outputData[10];
		int numDOFperNode = numDOF/4;
		for (int i=0; i<numDOFperNode; i++) {
			sprintf(outputData,"T1_P1_%d", i+1);
			output.tag("ResponseType", outputData);
		}
		for (int j=0; j<numDOFperNode; j++) {
			sprintf(outputData,"T1_P2_%d", j+1);
			output.tag("ResponseType", outputData);
		}
		
		for (int i=0; i<numDOFperNode; i++) {
			sprintf(outputData,"T2_P1_%d", i+1);
			output.tag("ResponseType", outputData);
		}
		for (int j=0; j<numDOFperNode; j++) {
			sprintf(outputData,"T2_P2_%d", j+1);
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
	} else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"-material") == 0) {
		if (argc > 1) {
			// we need at least one more argument otherwise 
			// there is no need to forward this call to the material
			// by default assume the old call style for backward compatibility "material result"
			int offset = 1;
			bool is_valid = true;
			// in case the user specifies the gauss point id... "material 1 result"
			if (argc > 2) {
				int sectionNum = atoi(argv[1]);
				if (sectionNum == 1) {
					// this is the only supported gauss id
					offset = 2;
				}
				else if (sectionNum > 1) {
					// this is a number, but not within the valid range
					is_valid = false;
				}
				// if it is 0, then it is not a number, forward it as usual...
			}
			if (is_valid) {
				output.tag("GaussPointOutput");
				output.attr("number", 1);
				output.attr("eta", 0.0);
				CompositeResponse* theCResponse = new CompositeResponse();
				Response* theResponse1 = theMaterial_1->setResponse(&argv[offset], argc - offset, output);
				Response* theResponse2 = theMaterial_2->setResponse(&argv[offset], argc - offset, output);
				theCResponse->addResponse(theResponse1);
				theCResponse->addResponse(theResponse2);
				theResponse = theCResponse;
				output.endTag();
			}
		}
	}

	output.endTag();
	return theResponse;
}

int 
N4BiaxialTruss::getResponse(int responseID, Information &eleInfo)
{
	double strain;

	switch (responseID) {
	case 1:
		return eleInfo.setVector(this->getResistingForce());

	default:
		return 0;
	}
}
