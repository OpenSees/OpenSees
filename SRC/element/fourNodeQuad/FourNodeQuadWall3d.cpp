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

// $Revision: 1.35 $
// $Date: 2009-10-13 21:14:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/FourNodeQuadWall3d.cpp,v $

// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
// Description: This file contains the class definition for FourNodeQuadWall3d.

#include <FourNodeQuadWall3d.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <DummyStream.h>

#include <elementAPI.h>

#include <TclModelBuilder.h>

double FourNodeQuadWall3d::shp[3][4];
double FourNodeQuadWall3d::pts[4][2];
double FourNodeQuadWall3d::wts[4];
int FourNodeQuadWall3d::dirns[3];

/* TO DO:
1) FIGURE OUT LOCAL ORIENTATION AND RELATIONSHIP WITH FSAM INPUT PARAMETERS
2) ADD INPUT PARAMETER FOR EFFECTIVE THICKNESS
3) BODY FORCES AND PRESSURE LOAD (maybe keep this if it will work OK)
4) MASS, DAMPING, INERTIAL LOAD
*/

//////////////////////////////////////////////////////////////////////////////////////////
 // INPUT READ FOR PARALLEL COMPUTING - FIGURE OUT HOW TO ORGANIZE NODES
void* OPS_FourNodeQuadWall3d()
{

	Element* theEle = 0;

	int numRemainingArgs = OPS_GetNumRemainingInputArgs();
	if (numRemainingArgs == 0) { // parallel processing
		theEle = new FourNodeQuadWall3d();
		return theEle;
	}

	if (numRemainingArgs < 7) {
		opserr << "ERROR - FourNodeQuadWall3d not enough args provided, want: element FourNodeQuadWall3d tag? iNode? jNode? kNode? lNode? thickness? matID? <p? rho? b1? b2?>\n";
	}

	// get the id and end nodes 
	int iData[6];
	double dData[5];
	dData[1] = 0.0;
	dData[2] = 0.0;
	dData[3] = 0.0;
	dData[4] = 0.0;

	int numData;
	int matTag = 0;
	int eleTag = 0;
	//char *pType;

	//  int iNode, jNode, kNode, lNode, Node1, Node2, Node3, Node4;

	numData = 5;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING element FourNodeQuadWall3d : invalid element data\n";
		return 0;
	}
	eleTag = iData[0];

	numData = 1;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "WARNING element FourNodeQuadWall3d : invalid thickness for element: " << eleTag << "\n";
		return 0;
	}

	numData = 1;
	if (OPS_GetIntInput(&numData, &matTag) != 0) {
		opserr << "WARNING element FourNodeQuadWall3d : invalid matTag for element: " << eleTag << "\n";
		//delete [] pType;
		return 0;
	}


	NDMaterial* theMaterial = OPS_GetNDMaterial(matTag);

	if (theMaterial == 0) {
		opserr << "WARNING material with tag " << matTag << "not found for element " << eleTag << endln;
		return 0;
	}

	if (numRemainingArgs == 11) {
		numData = 4;
		if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
			opserr << "WARNING element FourNodeQuadWall3d : invalid optional args for element: " << eleTag << "\n";
			return 0;
		}
	}

	// now create the truss and add it to the Domain
	theEle = new FourNodeQuadWall3d(eleTag, iData[1], iData[2], iData[3], iData[4],
		*theMaterial,
		dData[0], dData[1], dData[2], dData[3], dData[4]);

	if (theEle == 0) {
		opserr << "WARNING ran out of memory creating element with tag " << eleTag << endln;
		delete theMaterial;
		return 0;
	}

	return theEle;
}
// 
///////////////////////////////////////////////////////////////////////////////////

// Typical Constructor  - ADD DRILL AS OPTIONAL INPUT
FourNodeQuadWall3d::FourNodeQuadWall3d(int tag, int nd1, int nd2, int nd3, int nd4,
	NDMaterial& m,
	double t, double STIFFMODOUT, double POISSON,
	double p, double r, double b1, double b2)
	:Element(tag, ELE_TAG_FourNodeQuadWall3d),
	theMaterial(0), connectedExternalNodes(4),
	Q(24), pressureLoad(24), thickness(t), stiffModOut(STIFFMODOUT), Poisson(POISSON),
	applyLoad(0), pressure(p), rho(r),
	Kel(24, 24), Kgl(24, 24), KPlBend(12, 12), Mel(24, 24), Mgl(24, 24), Pel(24), Pgl(24), T(24, 24), Tt(3, 3), Ttstrain(6, 3),
	nd1Crds(3), nd2Crds(3), nd3Crds(3), nd4Crds(3) {

	Kel.Zero();
	Kgl.Zero();
	KPlBend.Zero();
	Mel.Zero();
	Mgl.Zero();
	Pel.Zero();
	Pgl.Zero();

	// Node coordinates in local cs
	nd1Crds.Zero();
	nd2Crds.Zero();
	nd3Crds.Zero();
	nd4Crds.Zero();

	pts[0][0] = -0.5773502691896258;
	pts[0][1] = -0.5773502691896258;
	pts[1][0] = 0.5773502691896258;
	pts[1][1] = -0.5773502691896258;
	pts[2][0] = 0.5773502691896258;
	pts[2][1] = 0.5773502691896258;
	pts[3][0] = -0.5773502691896258;
	pts[3][1] = 0.5773502691896258;

	wts[0] = 1.0;
	wts[1] = 1.0;
	wts[2] = 1.0;
	wts[3] = 1.0;

	dirns[0] = 1;
	dirns[1] = 1;
	dirns[2] = 1;

	// Body forces
	b[0] = b1;
	b[1] = b2;

	// Allocate arrays of pointers to NDMaterials
	theMaterial = new NDMaterial * [4];

	if (theMaterial == 0) {
		opserr << "FourNodeQuadWall3d::FourNodeQuadWall3d - failed allocate material model pointer\n";
		exit(-1);
	}

	int i;
	for (i = 0; i < 4; i++) {

		// Get copies of the material model for each integration point
		theMaterial[i] = m.getCopy();

		// Check allocation
		if (theMaterial[i] == 0) {
			opserr << "FourNodeQuadWall3d::FourNodeQuadWall3d -- failed to get a copy of material model\n";
			exit(-1);
		}
	}

	// Set connected external node IDs
	connectedExternalNodes(0) = nd1;
	connectedExternalNodes(1) = nd2;
	connectedExternalNodes(2) = nd3;
	connectedExternalNodes(3) = nd4;

	for (i = 0; i < 4; i++) {
		theNodes[i] = 0;
	}

}

// Blank Contructor
FourNodeQuadWall3d::FourNodeQuadWall3d()
	:Element(0, ELE_TAG_FourNodeQuadWall3d),
	theMaterial(0), connectedExternalNodes(4),
	Q(24), pressureLoad(24), thickness(0.0), stiffModOut(0.0), Poisson(0.0),
	applyLoad(0), pressure(0.0), rho(0.0),
	Kel(24, 24), Kgl(24, 24), KPlBend(12, 12), Mel(24, 24), Mgl(24, 24), Pel(24), Pgl(24), T(24, 24), Tt(3, 3), Ttstrain(6, 3),
	nd1Crds(3), nd2Crds(3), nd3Crds(3), nd4Crds(3) {

	Kel.Zero();
	Kgl.Zero();
	KPlBend.Zero();
	Mel.Zero();
	Mgl.Zero();
	Pel.Zero();
	Pgl.Zero();

	// Node coordinates in local cs
	nd1Crds.Zero();
	nd2Crds.Zero();
	nd3Crds.Zero();
	nd4Crds.Zero();

	pts[0][0] = -0.577350269189626;
	pts[0][1] = -0.577350269189626;
	pts[1][0] = 0.577350269189626;
	pts[1][1] = -0.577350269189626;
	pts[2][0] = 0.577350269189626;
	pts[2][1] = 0.577350269189626;
	pts[3][0] = -0.577350269189626;
	pts[3][1] = 0.577350269189626;

	wts[0] = 1.0;
	wts[1] = 1.0;
	wts[2] = 1.0;
	wts[3] = 1.0;

	dirns[0] = 1;
	dirns[1] = 1;
	dirns[2] = 1;

	for (int i = 0; i < 4; i++)
		theNodes[i] = 0;
}


FourNodeQuadWall3d::~FourNodeQuadWall3d() {
	for (int i = 0; i < 4; i++) {
		if (theMaterial[i])
			delete theMaterial[i];
	}

	// Delete the array of pointers to NDMaterial pointer arrays
	if (theMaterial)
		delete[] theMaterial;
}

// Standard OpenSees routine
int FourNodeQuadWall3d::getNumExternalNodes() const {
	return 4;
}

// Standard OpenSees routine
const ID& FourNodeQuadWall3d::getExternalNodes() {
	return connectedExternalNodes;
}

// Standard OpenSees routine
Node** FourNodeQuadWall3d::getNodePtrs(void) {
	return theNodes;
}

// Standard OpenSees routine
int FourNodeQuadWall3d::getNumDOF() {
	return 24;
}

// Set Domain
void FourNodeQuadWall3d::setDomain(Domain* theDomain) {

	// Check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		return;
	}

	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);
	int Nd3 = connectedExternalNodes(2);
	int Nd4 = connectedExternalNodes(3);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);
	theNodes[2] = theDomain->getNode(Nd3);
	theNodes[3] = theDomain->getNode(Nd4);

	if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
		opserr << "FATAL ERROR FourNodeQuadWall3d (tag: " << this->getTag() << " ) a node does not exist\n";
		exit(-1);
	}

	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();
	int dofNd3 = theNodes[2]->getNumberDOF();
	int dofNd4 = theNodes[3]->getNumberDOF();

	if (dofNd1 != 6 || dofNd2 != 6 || dofNd3 != 6 || dofNd4 != 6) {
		opserr << "FATAL ERROR FourNodeQuadWall3d (tag: " << this->getTag() << " ) needs ndf = 3\n";
		exit(-1);
	}

	this->DomainComponent::setDomain(theDomain);

	// Compute consistent nodal loads due to pressure
	this->setPressureLoadAtNodes();

	const Vector& crds1 = theNodes[0]->getCrds();
	const Vector& crds2 = theNodes[1]->getCrds();
	const Vector& crds3 = theNodes[2]->getCrds();
	const Vector& crds4 = theNodes[3]->getCrds();

	if (crds1.Size() != 3 || crds2.Size() != 3 || crds3.Size() != 3 || crds4.Size() != 3) {
		opserr << "FATAL ERROR FourNodeQuadWall3d (tag: " << this->getTag() << " ) needs ndm = 3\n";
		exit(-1);
	}

	// Check element orientation: Currently the element can be in x-y, y-z, or x-z plane ONLY
	for (int i = 0; i < 3; i++)
		dirns[i] = 1;

	if ((crds1(0) == crds2(0)) && (crds2(0) == crds3(0)) && (crds3(0) == crds4(0)))
		dirns[0] = 0;
	if ((crds1(1) == crds2(1)) && (crds2(1) == crds3(1)) && (crds3(1) == crds4(1)))
		dirns[1] = 0;
	if ((crds1(2) == crds2(2)) && (crds2(2) == crds3(2)) && (crds3(2) == crds4(2)))
		dirns[2] = 0;

	int sum = 0;
	for (int i = 0; i < 3; i++) {
		if (dirns[i] != 0 && sum < 2)
			dirn[sum] = i;
		sum += dirns[i];
	}

	if (sum != 2) {
		opserr << "DIRNS: " << dirns[0] << " " << dirns[1] << " " << dirns[2];
		theNodes[0]->Print(opserr);
		theNodes[1]->Print(opserr);
		theNodes[2]->Print(opserr);
		theNodes[3]->Print(opserr);
		opserr << "FATAL ERROR FourNodeQuadWall3d (tag: " << this->getTag() << " ) needs four nodes to be in x-y, y-z, or x-z plane\n";
		exit(-1);
	}

	// Calculate Transformation Matrix
	setTransformationMatrix();

	// Compute plate bending stiffness matrix ...............................................
	// Get Concrete Young's Modulus
	theResponses = new Response * [1];
	if (theResponses == 0) {
		opserr << " quadWall3D::quadWall3D - failed allocate responses array\n";
		exit(-1);
	}

	OPS_Stream* theDummyStream = new DummyStream();
	const char** argv = new const char* [1];

	argv[0] = "getInputParameters"; // to get input parameters from concrete material
	theResponses[0] = theMaterial[0]->setResponse(argv, 1, *theDummyStream);

	if (theResponses[0] == 0) {
		opserr << " FSAM::FSAM - failed to set appropriate materials tag: " << this->getTag() << "\n";
		exit(-1);
	}

	// Get ConcreteCM material input variables
	theResponses[0]->getResponse();
	Information& theInfoInput = theResponses[0]->getInformation();
	const Vector InputNDMat = theInfoInput.getData();

	Vector InputNDMaterial(InputNDMat.Size());

	for (int i = 0; i < InputNDMat.Size(); i++)
		InputNDMaterial[i] = InputNDMat[i];

	double E = InputNDMaterial[9];

	// Compute plate bending stiffness matrix
	computePlateTangent(w, h, E, stiffModOut, Poisson, thickness);

	return;
}


int FourNodeQuadWall3d::commitState() {

	int retVal = 0;

	// call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "FourNodeQuadWall3d::commitState () - failed in base class";
	}

	// Loop over the integration points and commit the material states
	for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->commitState();

	return retVal;
}


int FourNodeQuadWall3d::revertToLastCommit() {

	int retVal = 0;

	// Loop over the integration points and revert to last committed state
	for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToLastCommit();

	return retVal;
}


int FourNodeQuadWall3d::revertToStart() {

	int retVal = 0;

	// Loop over the integration points and revert states to start
	for (int i = 0; i < 4; i++)
		retVal += theMaterial[i]->revertToStart();

	return retVal;
}

// Update displacement field - OK by KK
int FourNodeQuadWall3d::update() {

	// Node Trial Displacements
	const Vector& dispGL1 = theNodes[0]->getTrialDisp();
	const Vector& dispGL2 = theNodes[1]->getTrialDisp();
	const Vector& dispGL3 = theNodes[2]->getTrialDisp();
	const Vector& dispGL4 = theNodes[3]->getTrialDisp();

	// Pack all displacements into one vector
	Vector DispGL(24); // Vector of displacements in global cs
	DispGL.Zero();

	for (int i = 0; i < 6; i++) {
		DispGL(i) = dispGL1(i);
		DispGL(i + 6) = dispGL2(i);
		DispGL(i + 12) = dispGL3(i);
		DispGL(i + 18) = dispGL4(i);
	}

	/*// DEBUG
	for (int i = 0; i < 24; i++) {
		opserr << "Dgl" << i+1 << " = " << DispGL(i) << "\n";
	}
	opserr <<"\n"; // */

	// Compute vector of displacements in local cs
	Vector DispLOC(24);
	DispLOC.Zero();

	// Coordinate transformation - from global to local
	// Dloc = T Dgl
	DispLOC.addMatrixVector(1.0, T, DispGL, 1.0);

	static double u[2][4];

	// Store in-plane element deformations (no drill for now) - ADD DRILL DOFs (maybe not necessary here)
	u[0][0] = DispLOC(0);  // DOF 1
	u[1][0] = DispLOC(1);  // DOF 2
	u[0][1] = DispLOC(6);  // DOF 7
	u[1][1] = DispLOC(7);  // DOF 8
	u[0][2] = DispLOC(12); // DOF 13
	u[1][2] = DispLOC(13); // DOF 14
	u[0][3] = DispLOC(18); // DOF 19
	u[1][3] = DispLOC(19); // DOF 20

	/*// DEBUG
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 2; j++) {
			opserr << "Dm" << i+1 << j+1 << " = " << u[j][i]  << "\n";
		}
	}
	opserr <<"\n"; // */

	// Obtain strains for in-plane action
	static Vector eps(3);

	int ret = 0;

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		this->shapeFunction(pts[i][0], pts[i][1]);

		// Interpolate strains
		// eps = B*u;
		// eps.addMatrixVector(0.0, B, u, 1.0);
		eps.Zero();

		for (int beta = 0; beta < 4; beta++) {
			eps(0) += shp[0][beta] * u[0][beta];
			eps(1) += shp[1][beta] * u[1][beta];
			eps(2) += shp[0][beta] * u[1][beta] + shp[1][beta] * u[0][beta];
		}

		// Set the material strain
		ret += theMaterial[i]->setTrialStrain(eps);
	}

	return ret;
}


// Compute Plate Tangent - OK by KK
void FourNodeQuadWall3d::computePlateTangent(double w, double h, double E, double KmodOut, double nu, double t) {

	Matrix Kp(13, 13); // temporary 
	Kp.Zero();

	KPlBend.Zero();

	double factor = KmodOut * E * pow(t, 3.0) / (12.0 * (1.0 - pow(nu, 2.0)));

	// Closed-form solution for stiffness matrix - see Matlab derivation
	Kp(1, 1) = (4.0 * h) / pow(w, 3.0) - (pow(h, 2.0) * ((4.0 * nu * pow(w, 4)) / 5.0 - (14.0 * pow(w, 4.0)) / 5.0) - 4.0 * pow(w, 6.0)) / (pow(h, 3.0) * pow(w, 5.0));
	Kp(1, 2) = ((4.0 * nu) / 5.0 + 1.0 / 5.0) / w + (2.0 * w) / pow(h, 2.0);
	Kp(1, 3) = -((4.0 * nu) / 5.0 + 1.0 / 5.0) / h - (2.0 * h) / pow(w, 2.0);
	Kp(1, 4) = (pow(h, 2.0) * ((4.0 * nu * pow(w, 4.0)) / 5.0 - (14.0 * pow(w, 4.0)) / 5.0) + 2.0 * pow(w, 6.0)) / (pow(h, 3.0) * pow(w, 5.0)) - (4.0 * h) / pow(w, 3.0);
	Kp(1, 5) = w / pow(h, 2.0) - ((4.0 * nu) / 5.0 + 1.0 / 5.0) / w;
	Kp(1, 6) = (nu / 5.0 - 1.0 / 5.0) / h - (2.0 * h) / pow(w, 2.0);
	Kp(1, 7) = -(2.0 * h) / pow(w, 3.0) - (pow(h, 2.0) * ((4.0 * nu * pow(w, 4.0)) / 5.0 - (14.0 * pow(w, 4)) / 5.0) + 2.0 * pow(w, 6.0)) / (pow(h, 3.0) * pow(w, 5.0));
	Kp(1, 8) = (nu / 5.0 - 1.0 / 5.0) / w + w / pow(h, 2.0);
	Kp(1, 9) = -(nu / 5.0 - 1.0 / 5.0) / h - h / pow(w, 2.0);
	Kp(1, 10) = (2.0 * h) / pow(w, 3.0) + (pow(h, 2.0) * ((4.0 * nu * pow(w, 4.0)) / 5.0 - (14 * pow(w, 4.0)) / 5.0) - 4.0 * pow(w, 6.0)) / (pow(h, 3.0) * pow(w, 5.0));
	Kp(1, 11) = (2.0 * w) / pow(h, 2.0) - (nu / 5.0 - 1.0 / 5.0) / w;
	Kp(1, 12) = ((4.0 * nu) / 5.0 + 1.0 / 5.0) / h - h / pow(w, 2.0);

	Kp(2, 1) = Kp(1, 2);
	Kp(2, 2) = (4.0 * w) / (3.0 * h) - (h * ((4.0 * nu) / 15.0 - 4.0 / 15.0)) / w;
	Kp(2, 3) = -nu;
	Kp(2, 4) = w / pow(h, 2.0) - ((4.0 * nu) / 5.0 + 1.0 / 5.0) / w;
	Kp(2, 5) = (2.0 * w) / (3.0 * h) + (h * ((4.0 * nu) / 15.0 - 4.0 / 15.0)) / w;
	Kp(2, 6) = 0.0;
	Kp(2, 7) = -(nu / 5.0 - 1.0 / 5.0) / w - w / pow(h, 2.0);
	Kp(2, 8) = w / (3.0 * h) - (h * (nu / 15.0 - 1.0 / 15.0)) / w;
	Kp(2, 9) = 0.0;
	Kp(2, 10) = (nu / 5.0 - 1.0 / 5.0) / w - (2.0 * w) / pow(h, 2.0);
	Kp(2, 11) = (2.0 * w) / (3.0 * h) + (h * (nu / 15.0 - 1.0 / 15.0)) / w;
	Kp(2, 12) = 0.0;

	Kp(3, 1) = Kp(1, 3);
	Kp(3, 2) = Kp(2, 3);
	Kp(3, 3) = (4.0 * h) / (3.0 * w) - (8.0 * w * (nu / 2.0 - 1.0 / 2.0)) / (15.0 * h);
	Kp(3, 4) = (2.0 * h) / pow(w, 2.0) - (nu / 5.0 - 1.0 / 5.0) / h;
	Kp(3, 5) = 0.0;
	Kp(3, 6) = (2.0 * h) / (3.0 * w) + (w * (nu / 15.0 - 1.0 / 15.0)) / h;
	Kp(3, 7) = (nu / 5.0 - 1.0 / 5.0) / h + h / pow(w, 2.0);
	Kp(3, 8) = 0.0;
	Kp(3, 9) = h / (3.0 * w) - (w * (nu / 15.0 - 1.0 / 15.0)) / h;
	Kp(3, 10) = ((4.0 * nu) / 5.0 + 1.0 / 5.0) / h - h / pow(w, 2.0);
	Kp(3, 11) = 0.0;
	Kp(3, 12) = (2.0 * h) / (3.0 * w) + (8.0 * w * (nu / 2.0 - 1.0 / 2.0)) / (15.0 * h);

	Kp(4, 1) = Kp(1, 4);
	Kp(4, 2) = Kp(2, 4);
	Kp(4, 3) = Kp(3, 4);
	Kp(4, 4) = (4.0 * h) / pow(w, 3.0) - (pow(h, 2.0) * ((4.0 * nu * pow(w, 4.0)) / 5.0 - (14.0 * pow(w, 4.0)) / 5.0) - 4.0 * pow(w, 6.0)) / (pow(h, 3.0) * pow(w, 5.0));
	Kp(4, 5) = ((4.0 * nu) / 5.0 + 1.0 / 5.0) / w + (2.0 * w) / pow(h, 2.0);
	Kp(4, 6) = ((4.0 * nu) / 5.0 + 1.0 / 5.0) / h + (2.0 * h) / pow(w, 2.0);
	Kp(4, 7) = (2.0 * h) / pow(w, 3.0) + (pow(h, 2.0) * ((4.0 * nu * pow(w, 4.0)) / 5.0 - (14.0 * pow(w, 4.0)) / 5.0) - 4.0 * pow(w, 6.0)) / (pow(h, 3.0) * pow(w, 5.0));
	Kp(4, 8) = (2.0 * w) / pow(h, 2.0) - (nu / 5.0 - 1.0 / 5.0) / w;
	Kp(4, 9) = h / pow(w, 2.0) - ((4.0 * nu) / 5.0 + 1.0 / 5.0) / h;
	Kp(4, 10) = -(2.0 * h) / pow(w, 3.0) - (pow(h, 2.0) * ((4.0 * nu * pow(w, 4.0)) / 5.0 - (14.0 * pow(w, 4.0)) / 5.0) + 2.0 * pow(w, 6.0)) / (pow(h, 3.0) * pow(w, 5.0));
	Kp(4, 11) = (nu / 5.0 - 1.0 / 5.0) / w + w / pow(h, 2.0);
	Kp(4, 12) = (nu / 5.0 - 1.0 / 5.0) / h + h / pow(w, 2.0);

	Kp(5, 1) = Kp(1, 5);
	Kp(5, 2) = Kp(2, 5);
	Kp(5, 3) = Kp(3, 5);
	Kp(5, 4) = Kp(4, 5);
	Kp(5, 5) = (4.0 * w) / (3.0 * h) - (h * ((4.0 * nu) / 15.0 - 4.0 / 15.0)) / w;
	Kp(5, 6) = nu;
	Kp(5, 7) = (nu / 5.0 - 1.0 / 5.0) / w - (2.0 * w) / pow(h, 2.0);
	Kp(5, 8) = (2.0 * w) / (3.0 * h) + (h * (nu / 15.0 - 1.0 / 15.0)) / w;
	Kp(5, 9) = 0.0;
	Kp(5, 10) = -(nu / 5.0 - 1.0 / 5.0) / w - w / pow(h, 2.0);
	Kp(5, 11) = w / (3.0 * h) - (h * (nu / 15.0 - 1.0 / 15.0)) / w;
	Kp(5, 12) = 0.0;

	Kp(6, 1) = Kp(1, 6);
	Kp(6, 2) = Kp(2, 6);
	Kp(6, 3) = Kp(3, 6);
	Kp(6, 4) = Kp(4, 6);
	Kp(6, 5) = Kp(5, 6);
	Kp(6, 6) = (4.0 * h) / (3.0 * w) - (w * ((4.0 * nu) / 15.0 - 4.0 / 15.0)) / h;
	Kp(6, 7) = h / pow(w, 2.0) - ((4.0 * nu) / 5.0 + 1.0 / 5.0) / h;
	Kp(6, 8) = 0.0;
	Kp(6, 9) = (2.0 * h) / (3.0 * w) + (8.0 * w * (nu / 2.0 - 1.0 / 2.0)) / (15.0 * h);
	Kp(6, 10) = -(nu / 5.0 - 1.0 / 5.0) / h - h / pow(w, 2.0);
	Kp(6, 11) = 0.0;
	Kp(6, 12) = h / (3.0 * w) - (w * (nu / 15.0 - 1.0 / 15.0)) / h;

	Kp(7, 1) = Kp(1, 7);
	Kp(7, 2) = Kp(2, 7);
	Kp(7, 3) = Kp(3, 7);
	Kp(7, 4) = Kp(4, 7);
	Kp(7, 5) = Kp(5, 7);
	Kp(7, 6) = Kp(6, 7);
	Kp(7, 7) = (4.0 * h) / pow(w, 3.0) - (pow(h, 2.0) * ((4.0 * nu * pow(w, 4.0)) / 5.0 - (14.0 * pow(w, 4.0)) / 5.0) - 4 * pow(w, 6.0)) / (pow(h, 3.0) * pow(w, 5.0));
	Kp(7, 8) = -((4.0 * nu) / 5.0 + 1.0 / 5.0) / w - (2.0 * w) / pow(h, 2.0);
	Kp(7, 9) = ((4.0 * nu) / 5.0 + 1.0 / 5.0) / h + (2.0 * h) / pow(w, 2.0);
	Kp(7, 10) = (pow(h, 2.0) * ((4.0 * nu * pow(w, 2.0)) / 5.0 - (14.0 * pow(w, 2.0)) / 5.0) + 2.0 * pow(w, 4.0)) / (pow(h, 3.0) * pow(w, 3.0)) - (4.0 * h) / pow(w, 3.0);
	Kp(7, 11) = ((4.0 * nu) / 5.0 + 1.0 / 5.0) / w - w / pow(h, 2.0);
	Kp(7, 12) = (2.0 * h) / pow(w, 2.0) - (nu / 5.0 - 1.0 / 5.0) / h;

	Kp(8, 1) = Kp(1, 8);
	Kp(8, 2) = Kp(2, 8);
	Kp(8, 3) = Kp(3, 8);
	Kp(8, 4) = Kp(4, 8);
	Kp(8, 5) = Kp(5, 8);
	Kp(8, 6) = Kp(6, 8);
	Kp(8, 7) = Kp(7, 8);
	Kp(8, 8) = (4.0 * w) / (3.0 * h) - (h * ((4.0 * nu) / 15.0 - 4.0 / 15.0)) / w;
	Kp(8, 9) = -nu;
	Kp(8, 10) = ((4.0 * nu) / 5.0 + 1.0 / 5.0) / w - w / pow(h, 2.0);
	Kp(8, 11) = (2.0 * w) / (3.0 * h) + (h * ((4.0 * nu) / 15.0 - 4.0 / 15.0)) / w;
	Kp(8, 12) = 0.0;

	Kp(9, 1) = Kp(1, 9);
	Kp(9, 2) = Kp(2, 9);
	Kp(9, 3) = Kp(3, 9);
	Kp(9, 4) = Kp(4, 9);
	Kp(9, 5) = Kp(5, 9);
	Kp(9, 6) = Kp(6, 9);
	Kp(9, 7) = Kp(7, 9);
	Kp(9, 8) = Kp(8, 9);
	Kp(9, 9) = (4.0 * h) / (3.0 * w) - (8.0 * w * (nu / 2.0 - 1.0 / 2.0)) / (15.0 * h);
	Kp(9, 10) = (nu / 5.0 - 1.0 / 5.0) / h - (2.0 * h) / pow(w, 2.0);
	Kp(9, 11) = 0.0;
	Kp(9, 12) = (2.0 * h) / (3.0 * w) + (w * (nu / 15.0 - 1.0 / 15.0)) / h;

	Kp(10, 1) = Kp(1, 10);
	Kp(10, 2) = Kp(2, 10);
	Kp(10, 3) = Kp(3, 10);
	Kp(10, 4) = Kp(4, 10);
	Kp(10, 5) = Kp(5, 10);
	Kp(10, 6) = Kp(6, 10);
	Kp(10, 7) = Kp(7, 10);
	Kp(10, 8) = Kp(8, 10);
	Kp(10, 9) = Kp(9, 10);
	Kp(10, 10) = (4.0 * h) / pow(w, 3.0) - (pow(h, 2.0) * ((4.0 * nu * pow(w, 4.0)) / 5.0 - (14.0 * pow(w, 4.0)) / 5.0) - 4.0 * pow(w, 6.0)) / (pow(h, 3.0) * pow(w, 5.0));
	Kp(10, 11) = -((4.0 * nu) / 5.0 + 1.0 / 5.0) / w - (2.0 * w) / pow(h, 2.0);
	Kp(10, 12) = -((4.0 * nu) / 5.0 + 1.0 / 5.0) / h - (2.0 * h) / pow(w, 2.0);

	Kp(11, 1) = Kp(1, 11);
	Kp(11, 2) = Kp(2, 11);
	Kp(11, 3) = Kp(3, 11);
	Kp(11, 4) = Kp(4, 11);
	Kp(11, 5) = Kp(5, 11);
	Kp(11, 6) = Kp(6, 11);
	Kp(11, 7) = Kp(7, 11);
	Kp(11, 8) = Kp(8, 11);
	Kp(11, 9) = Kp(9, 11);
	Kp(11, 10) = Kp(10, 11);
	Kp(11, 11) = (4.0 * w) / (3.0 * h) - (h * ((4.0 * nu) / 15.0 - 4.0 / 15.0)) / w;
	Kp(11, 12) = nu;

	Kp(12, 1) = Kp(1, 12);
	Kp(12, 2) = Kp(2, 12);
	Kp(12, 3) = Kp(3, 12);
	Kp(12, 4) = Kp(4, 12);
	Kp(12, 5) = Kp(5, 12);
	Kp(12, 6) = Kp(6, 12);
	Kp(12, 7) = Kp(7, 12);
	Kp(12, 8) = Kp(8, 12);
	Kp(12, 9) = Kp(9, 12);
	Kp(12, 10) = Kp(10, 12);
	Kp(12, 11) = Kp(11, 12);
	Kp(12, 12) = (4.0 * h) / (3.0 * w) - (8.0 * w * (nu / 2.0 - 1.0 / 2.0)) / (15.0 * h);

	// Multiply with factor and assign to a matrix
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			KPlBend(i, j) = factor * Kp(i + 1, j + 1);
		}
	}

	/*// DEBUG
	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			opserr << KPlBend(i,j) << "\t";
		}
		opserr << "\n";opserr << "\n";
	} // */

}

// Membrane (in-plane) shape functions. Bi-linear interpolation (textbook) - OK by KK
double FourNodeQuadWall3d::shapeFunction(double xi, double eta) {

	double oneMinuseta = 1.0 - eta;
	double onePluseta = 1.0 + eta;
	double oneMinusxi = 1.0 - xi;
	double onePlusxi = 1.0 + xi;

	shp[2][0] = 0.25 * oneMinusxi * oneMinuseta;	// N_1
	shp[2][1] = 0.25 * onePlusxi * oneMinuseta;		// N_2
	shp[2][2] = 0.25 * onePlusxi * onePluseta;		// N_3
	shp[2][3] = 0.25 * oneMinusxi * onePluseta;		// N_4

	double J[2][2];

	// Jacobian (based on node locations in local cs)
	J[0][0] = 0.25 * (-nd1Crds(0) * oneMinuseta + nd2Crds(0) * oneMinuseta +
		nd3Crds(0) * (onePluseta)-nd4Crds(0) * (onePluseta));

	J[0][1] = 0.25 * (-nd1Crds(0) * oneMinusxi - nd2Crds(0) * onePlusxi +
		nd3Crds(0) * onePlusxi + nd4Crds(0) * oneMinusxi);

	J[1][0] = 0.25 * (-nd1Crds(1) * oneMinuseta + nd2Crds(1) * oneMinuseta +
		nd3Crds(1) * onePluseta - nd4Crds(1) * onePluseta);

	J[1][1] = 0.25 * (-nd1Crds(1) * oneMinusxi - nd2Crds(1) * onePlusxi +
		nd3Crds(1) * onePlusxi + nd4Crds(1) * oneMinusxi);

	double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
	double oneOverdetJ = 1.0 / detJ;
	double L[2][2];

	// L = inv(J)
	L[0][0] = J[1][1] * oneOverdetJ;
	L[1][0] = -J[0][1] * oneOverdetJ;
	L[0][1] = -J[1][0] * oneOverdetJ;
	L[1][1] = J[0][0] * oneOverdetJ;

	double L00 = 0.25 * L[0][0];
	double L10 = 0.25 * L[1][0];
	double L01 = 0.25 * L[0][1];
	double L11 = 0.25 * L[1][1];

	double L00oneMinuseta = L00 * oneMinuseta;
	double L00onePluseta = L00 * onePluseta;
	double L01oneMinusxi = L01 * oneMinusxi;
	double L01onePlusxi = L01 * onePlusxi;

	double L10oneMinuseta = L10 * oneMinuseta;
	double L10onePluseta = L10 * onePluseta;
	double L11oneMinusxi = L11 * oneMinusxi;
	double L11onePlusxi = L11 * onePlusxi;

	// See Cook, Malkus, Plesha p. 169 for the derivation of these terms
	shp[0][0] = -L00oneMinuseta - L01oneMinusxi;	// N_1,1
	shp[0][1] = L00oneMinuseta - L01onePlusxi;		// N_2,1
	shp[0][2] = L00onePluseta + L01onePlusxi;		// N_3,1
	shp[0][3] = -L00onePluseta + L01oneMinusxi;	// N_4,1

	shp[1][0] = -L10oneMinuseta - L11oneMinusxi;	// N_1,2
	shp[1][1] = L10oneMinuseta - L11onePlusxi;		// N_2,2
	shp[1][2] = L10onePluseta + L11onePlusxi;		// N_3,2
	shp[1][3] = -L10onePluseta + L11oneMinusxi;	// N_4,2

	return detJ;
}

// Initial stiffness - Ok by KK
const Matrix& FourNodeQuadWall3d::getInitialStiff() {

	Kel.Zero(); // element stiffness matrix in local coordinate system
	Kgl.Zero(); // element stiffness matrix in global coordinate system

	// Obtain membrane (in-lane) matrix ........................................
	Matrix Km(8, 8);
	Km.Zero();

	double dvol;
	double DB[3][2];

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1]);
		dvol *= (thickness * wts[i]);

		// Get the material tangent
		const Matrix& D = theMaterial[i]->getInitialTangent();

		// Perform numerical integration
		//K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
		//K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);

		double D00 = D(0, 0); double D01 = D(0, 1); double D02 = D(0, 2);
		double D10 = D(1, 0); double D11 = D(1, 1); double D12 = D(1, 2);
		double D20 = D(2, 0); double D21 = D(2, 1); double D22 = D(2, 2);

		/*// DEBUG
		opserr << D00 << "\t" << D01 <<  "\t" << D02 << "\n";
		opserr << D10 << "\t" << D11 <<  "\t" << D12 << "\n";
		opserr << D20 << "\t" << D21 << "\t" << D22 << "\n";
		opserr <<  "\n"; // */

		//	 for (int beta = 0, ib = 0, colIb =0, colIbP1 = 8; 
		//   beta < 4; 
		//   beta++, ib += 2, colIb += 16, colIbP1 += 16) {

		for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
			for (int beta = 0, ib = 0; beta < 4; beta++, ib += 2) {

				// Local c.s. = global c.s.
				DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
				DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
				DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
				DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
				DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
				DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);

				Km(ia, ib) += shp[0][alpha] * DB[0][0] + shp[1][alpha] * DB[2][0];
				Km(ia, ib + 1) += shp[0][alpha] * DB[0][1] + shp[1][alpha] * DB[2][1];
				Km(ia + 1, ib) += shp[1][alpha] * DB[1][0] + shp[0][alpha] * DB[2][0];
				Km(ia + 1, ib + 1) += shp[1][alpha] * DB[1][1] + shp[0][alpha] * DB[2][1];

			}
		}
	}

	/*// DEBUG
	opserr << "Kmembrane = \n";
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			opserr << Km(i,j) << "\t";
		}
		opserr << "\n";opserr << "\n";
	} // */

	// End membrane (in-lane) matrix ..........................................

	// Assemble global K
	// Membrane (in-plane) terms
	int tm[8]; tm[0] = 0; tm[1] = 1; tm[2] = 6; tm[3] = 7; tm[4] = 12; tm[5] = 13; tm[6] = 18; tm[7] = 19;

	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			Kel(tm[i], tm[j]) = Km(i, j);
		}
	}



	// Plate (out-of-plane) terms
	int tp[12]; tp[0] = 2; tp[1] = 3; tp[2] = 4; tp[3] = 8; tp[4] = 9; tp[5] = 10; tp[6] = 14; tp[7] = 15; tp[8] = 16; tp[9] = 20; tp[10] = 21; tp[11] = 22;

	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			Kel(tp[i], tp[j]) = KPlBend(i, j);
		}
	}

	// Drill temrs (not considered currently) - ADD AS INPUT PARAMETER !!!
	Kel(5, 5) = 1.0;
	Kel(11, 11) = 1.0;
	Kel(17, 17) = 1.0;
	Kel(23, 23) = 1.0;

	/*// DEBUG
	for (int i = 0; i < 24; i++) {
		for (int j = 0; j < 24; j++) {
			opserr << Kel(i,j) << "\t";
		}
		opserr << "\n";opserr << "\n";
	} // */

	// Perform coordinate transformation - local to global
	// Kgl = Kgl + (T^ Kel * T);
	Kgl.addMatrixTripleProduct(1.0, T, Kel, 1.0);

	return Kgl; // return element stiffness matrix in global cs
}

// Tangent stiffness - OK by KK
const Matrix& FourNodeQuadWall3d::getTangentStiff() {

	Kel.Zero();
	Kgl.Zero();

	// Obtain membrane (in-plane) matrix ........................................
	Matrix Km(8, 8);
	Km.Zero();

	double dvol;
	double DB[3][2];

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1]);
		dvol *= (thickness * wts[i]);

		// Get the material tangent
		const Matrix& D = theMaterial[i]->getTangent();

		// Perform numerical integration
		//K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
		//K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);

		double D00 = D(0, 0); double D01 = D(0, 1); double D02 = D(0, 2);
		double D10 = D(1, 0); double D11 = D(1, 1); double D12 = D(1, 2);
		double D20 = D(2, 0); double D21 = D(2, 1); double D22 = D(2, 2);

		/*// DEBUG
		opserr << "D = \n";
		opserr << D00 << "\t" << D01 <<  "\t" << D02 << "\n";
		opserr << D10 << "\t" << D11 <<  "\t" << D12 << "\n";
		opserr << D20 << "\t" << D21 << "\t" << D22 << "\n";
		opserr <<  "\n"; // */

		for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {
			for (int beta = 0, ib = 0; beta < 4; beta++, ib += 2) {

				// Current assumption: Local c.s. = global c.s.
				DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
				DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
				DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
				DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
				DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
				DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);

				Km(ia, ib) += shp[0][alpha] * DB[0][0] + shp[1][alpha] * DB[2][0];
				Km(ia, ib + 1) += shp[0][alpha] * DB[0][1] + shp[1][alpha] * DB[2][1];
				Km(ia + 1, ib) += shp[1][alpha] * DB[1][0] + shp[0][alpha] * DB[2][0];
				Km(ia + 1, ib + 1) += shp[1][alpha] * DB[1][1] + shp[0][alpha] * DB[2][1];

			}
		}
	}

	/*// DEBUG
	opserr << "Kmembrane = \n";
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			opserr << Km(i,j) << "\t";
		}
		opserr << "\n";opserr << "\n";
	} // */

	// End membrane (in-plane) matrix ..........................................

	// Assemble global K ......................
	// Membrane (in-plane) terms
	int tm[8]; tm[0] = 0; tm[1] = 1; tm[2] = 6; tm[3] = 7; tm[4] = 12; tm[5] = 13; tm[6] = 18; tm[7] = 19;

	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			Kel(tm[i], tm[j]) = Km(i, j);
		}
	}

	// Plate (out-of-plane) terms
	int tp[12]; tp[0] = 2; tp[1] = 3; tp[2] = 4; tp[3] = 8; tp[4] = 9; tp[5] = 10; tp[6] = 14; tp[7] = 15; tp[8] = 16; tp[9] = 20; tp[10] = 21; tp[11] = 22;

	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			Kel(tp[i], tp[j]) = KPlBend(i, j);
		}
	}

	// Drill temrs (not considered currently) - ADD AS INPUT PARAMETER !!!
	Kel(5, 5) = 1.0;
	Kel(11, 11) = 1.0;
	Kel(17, 17) = 1.0;
	Kel(23, 23) = 1.0;

	/*// DEBUG
	for (int i = 0; i < 24; i++) {
		for (int j = 0; j < 24; j++) {
			opserr << Kel(i,j) << "\t";
		}
		opserr << "\n";opserr << "\n";
	} // */

	// Perform coordinate transformation - local to global
	// Kgl = Kgl + (T^ Kel * T);
	Kgl.addMatrixTripleProduct(1.0, T, Kel, 1.0);

	return Kgl; // return element stiffness matrix in global cs
}

// Force vector - OK by KK
const Vector& FourNodeQuadWall3d::getResistingForce() {

	Pel.Zero();
	Pgl.Zero();

	// Compute membrane (in-plane) forces
	Vector Pm(8);
	Pm.Zero();

	double dvol;

	// Loop over the integration points
	for (int i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		dvol = this->shapeFunction(pts[i][0], pts[i][1]);
		dvol *= (thickness * wts[i]);

		// Get material stress response
		const Vector& sigma = theMaterial[i]->getStress();

		// Perform numerical integration on internal force
		//P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
		//P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);

		for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 2) {

			Pm(ia) += dvol * (shp[0][alpha] * sigma(0) + shp[1][alpha] * sigma(2));
			Pm(ia + 1) += dvol * (shp[1][alpha] * sigma(1) + shp[0][alpha] * sigma(2));

			// Subtract equiv. body forces from the nodes
			//P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
			//P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
			if (applyLoad == 0) {
				Pm(ia) -= dvol * (shp[2][alpha] * b[0]);
				Pm(ia + 1) -= dvol * (shp[2][alpha] * b[1]);
			}
			else {
				Pm(ia) -= dvol * (shp[2][alpha] * appliedB[0]);
				Pm(ia + 1) -= dvol * (shp[2][alpha] * appliedB[1]);
			}
		}

	}

	/*// DEBUG
	for (int i = 0; i < 8; i++) {
			opserr << Pm(i) << "\t";
		opserr << "\n";opserr << "\n";
	} // */

	// Obtain plate (out-of-plane) trial displacements
	Matrix Dplate(12, 1);
	Dplate.Zero();

	const Vector& dispGL1 = theNodes[0]->getTrialDisp();
	const Vector& dispGL2 = theNodes[1]->getTrialDisp();
	const Vector& dispGL3 = theNodes[2]->getTrialDisp();
	const Vector& dispGL4 = theNodes[3]->getTrialDisp();

	// Pack all displacements into one vector
	Vector DispGL(24); // Vector of displacements in global cs
	DispGL.Zero();

	for (int i = 0; i < 6; i++) {
		DispGL(i) = dispGL1(i);
		DispGL(i + 6) = dispGL2(i);
		DispGL(i + 12) = dispGL3(i);
		DispGL(i + 18) = dispGL4(i);
	}

	// Coordinate transformation - from global to local
	Vector DispLOC(24); // Vector of displacements in local cs
	DispLOC.Zero();

	// Dloc = T Dgl
	DispLOC.addMatrixVector(1.0, T, DispGL, 1.0);

	Dplate(0, 0) = DispLOC(2); // DOF 3
	Dplate(1, 0) = DispLOC(3); // DOF 4
	Dplate(2, 0) = DispLOC(4); // DOF 5
	Dplate(3, 0) = DispLOC(8); // DOF 9
	Dplate(4, 0) = DispLOC(9); // DOF 10
	Dplate(5, 0) = DispLOC(10); // DOF 11
	Dplate(6, 0) = DispLOC(14); // DOF 15
	Dplate(7, 0) = DispLOC(15); // DOF 16
	Dplate(8, 0) = DispLOC(16); // DOF 17
	Dplate(9, 0) = DispLOC(20); // DOF 21
	Dplate(10, 0) = DispLOC(21); // DOF 22
	Dplate(11, 0) = DispLOC(22); // DOF 23

	/*// DEBUG
	for (int i = 0; i<12; i++) {
		opserr << "Dplate" << i+1 << " = " << Dplate(i,0) << "\n";
	} ; opserr << "\n"; // */

	// Compute plate (out-of-plane) forces Pp = Kplate*Dplate
	Matrix Pp(12, 1);
	Pp.Zero();

	Pp.addMatrixProduct(0.0, KPlBend, Dplate, 1.0);

	/*// DEBUG
	for (int i = 0; i<12; i++) {
		opserr << Pp(i,0) << "\n";
	} ; opserr << "\n"; // */

	// Assemble global force vector
	// Membrane (in-plane) terms
	int tm[8]; tm[0] = 0; tm[1] = 1; tm[2] = 6; tm[3] = 7; tm[4] = 12; tm[5] = 13; tm[6] = 18; tm[7] = 19;

	for (int i = 0; i < 8; i++) {
		Pel(tm[i]) = Pm(i);
	}

	// Plate (out-of-plane) terms
	int tp[12]; tp[0] = 2; tp[1] = 3; tp[2] = 4; tp[3] = 8; tp[4] = 9; tp[5] = 10; tp[6] = 14; tp[7] = 15; tp[8] = 16; tp[9] = 20; tp[10] = 21; tp[11] = 22;

	for (int i = 0; i < 12; i++) {
		Pel(tp[i]) = Pp(i, 0);
	}

	// Drill temrs(not considered currently) - ADD Pdrill = Kdrill x DOF 
	Pel(5) = 0.0;
	Pel(11) = 0.0;
	Pel(17) = 0.0;
	Pel(23) = 0.0;

	// Subtract pressure loading from resisting force
	if (pressure != 0.0) {
		//P = P - pressureLoad;
		Pel.addVector(1.0, pressureLoad, -1.0);
	}

	// Subtract other external nodal loads ... P_res = P_int - P_ext
	//P = P - Q;
	Pel.addVector(1.0, Q, -1.0);

	/*// DEBUG
	for (int i = 0; i<24; i++) {
		opserr << Pel(i) << "\n";
	} ; opserr << "\n"; // */

	// Perform coordinate transformation - local to global
	// Pgl = Pgl + T^ Pel;
	Pgl.addMatrixTransposeVector(1.0, T, Pel, 1.0);

	return Pgl; // return element force vector in global cs
}


// Transformation matrix using unit vectors
void  FourNodeQuadWall3d::setTransformationMatrix(void) {

	T.Zero();			// element transformation matrix
	Tt.Zero();			// node transformation matrix
	Ttstrain.Zero();	// strain/stress output transformation matrix

	const Vector& Nd1Crd = theNodes[0]->getCrds();
	const Vector& Nd2Crd = theNodes[1]->getCrds();
	const Vector& Nd3Crd = theNodes[2]->getCrds();
	const Vector& Nd4Crd = theNodes[3]->getCrds();

	// Define local axis: 
	// x:  Nd1 -> Nd2
	// y:  Nd1 -> Nd4
	// z:  (x) x (y)

	// Vector components, magnitudes and unit vectors
	double Xx, Xy, Xz, X_, Xex, Xey, Xez;
	double Yx, Yy, Yz, Y_, Yex, Yey, Yez;
	double Zex, Zey, Zez;

	// Components of local X axis
	Xx = Nd2Crd(0) - Nd1Crd(0);
	Xy = Nd2Crd(1) - Nd1Crd(1);
	Xz = Nd2Crd(2) - Nd1Crd(2);

	// Magnitude
	X_ = pow(pow(Xx, 2.0) + pow(Xy, 2.0) + pow(Xz, 2.0), 0.5);

	w = X_; // Element width

	// unit x components
	Xex = Xx / X_;
	Xey = Xy / X_;
	Xez = Xz / X_;

	// Components of local Y axis
	Yx = Nd4Crd(0) - Nd1Crd(0);
	Yy = Nd4Crd(1) - Nd1Crd(1);
	Yz = Nd4Crd(2) - Nd1Crd(2);

	// Magnitude
	Y_ = pow(pow(Yx, 2) + pow(Yy, 2) + pow(Yz, 2), 0.5);

	h = Y_;  // Element height

	// unit y components
	Yex = Yx / Y_;
	Yey = Yy / Y_;
	Yez = Yz / Y_;

	// (Ze) = (Xe) x (Ye)
	Zex = Xey * Yez - Xez * Yey;
	Zey = -(Xex * Yez - Xez * Yex);
	Zez = Xex * Yey - Xey * Yex;

	// Node Transformation Matrix
	Tt(0, 0) = Xex;
	Tt(1, 0) = Yex;
	Tt(2, 0) = Zex;
	Tt(0, 1) = Xey;
	Tt(1, 1) = Yey;
	Tt(2, 1) = Zey;
	Tt(0, 2) = Xez;
	Tt(1, 2) = Yez;
	Tt(2, 2) = Zez;

	/* // DEBUG
	opserr << "NODE Transformation Matrix: \n";
	for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	opserr << Tt(i,j) << "\t";
	}
	opserr << "\n";opserr << "\n";

	} // */

	// Element Transformation Matrix
	for (int j = 0; j < 24; j += 3) {

		T(j + 0, j + 0) = Xex;
		T(j + 1, j + 0) = Yex;
		T(j + 2, j + 0) = Zex;
		T(j + 0, j + 1) = Xey;
		T(j + 1, j + 1) = Yey;
		T(j + 2, j + 1) = Zey;
		T(j + 0, j + 2) = Xez;
		T(j + 1, j + 2) = Yez;
		T(j + 2, j + 2) = Zez;

	}

	/*// DEBUG
	opserr << "ELEMENT Transformation Matrix: \n";
	for (int i = 0; i < 24; i++) {
		for (int j = 0; j < 24; j++) {
			opserr << T(i,j) << "\t";
		}
		opserr << "\n";opserr << "\n";

	} // */

	// Obtain node coordinates in local cs
	// Get node coordinates in global cs
	const Vector& nd1CrdsGL = theNodes[0]->getCrds();
	const Vector& nd2CrdsGL = theNodes[1]->getCrds();
	const Vector& nd3CrdsGL = theNodes[2]->getCrds();
	const Vector& nd4CrdsGL = theNodes[3]->getCrds();

	/*// DEBUG - global node coordinates
	opserr << "GLOBAL Node Coordinates: \n";
	for (int i = 0; i < 3; i++)
	opserr << nd1CrdsGL(i) << "\t";

	opserr << "\n";opserr << "\n";

	for (int i = 0; i < 3; i++)
	opserr << nd2CrdsGL(i) << "\t";

	opserr << "\n";opserr << "\n";

	for (int i = 0; i < 3; i++)
	opserr << nd3CrdsGL(i) << "\t";

	opserr << "\n";opserr << "\n";

	for (int i = 0; i < 3; i++)
	opserr << nd4CrdsGL(i) << "\t";

	opserr << "\n";opserr << "\n";
	opserr << "\n";opserr << "\n";
	// */

	// Coordinate transformation - from global to local
	nd1Crds.addMatrixVector(1.0, Tt, nd1CrdsGL, 1.0);
	nd2Crds.addMatrixVector(1.0, Tt, nd2CrdsGL, 1.0);
	nd3Crds.addMatrixVector(1.0, Tt, nd3CrdsGL, 1.0);
	nd4Crds.addMatrixVector(1.0, Tt, nd4CrdsGL, 1.0);

	/*// DEBUG - local node coordinates
	opserr << "LOCAL Node Coordinates: \n";
	for (int i = 0; i < 3; i++)
	opserr << nd1Crds(i) << "\t";

	opserr << "\n";opserr << "\n";

	for (int i = 0; i < 3; i++)
	opserr << nd2Crds(i) << "\t";

	opserr << "\n";opserr << "\n";

	for (int i = 0; i < 3; i++)
	opserr << nd3Crds(i) << "\t";

	opserr << "\n";opserr << "\n";

	for (int i = 0; i < 3; i++)
	opserr << nd4Crds(i) << "\t";

	opserr << "\n";opserr << "\n";
	// */

	// Transformation Matrix for Strain/Stress output only - works only if elements are in planes of the coordinate system !!!
	if (dirns[2] == 0) { // xy plane
		Ttstrain(0, 0) = 1.0;
		Ttstrain(1, 1) = 1.0;
		Ttstrain(3, 2) = 1.0;
	}
	else if (dirns[0] == 0) { // yz plane
		Ttstrain(1, 0) = 1.0;
		Ttstrain(2, 1) = 1.0;
		Ttstrain(4, 2) = 1.0;
	}
	else if (dirns[1] == 0) { // zx plane
		Ttstrain(0, 1) = 1.0;
		Ttstrain(2, 0) = 1.0;
		Ttstrain(5, 2) = 1.0;
	}
	else {

	}
}


/* FOR DYNAMIC ANALYSIS *//////////////////////////////////////////////////////////////////////////
// Element Mass Matrix - OK by KK, comment: includes vertical mass !!!
const Matrix& FourNodeQuadWall3d::getMass() {

	Mel.Zero();
	Mgl.Zero();

	int i;
	static double rhoi[4];
	double sum = 0.0;
	for (i = 0; i < 4; i++) {
		if (rho == 0)
			rhoi[i] = theMaterial[i]->getRho();
		else
			rhoi[i] = rho;
		sum += rhoi[i];
	}

	if (sum == 0.0)
		return Mel;

	double rhodvol, Nrho;

	// Compute a lumped mass matrix
	for (i = 0; i < 4; i++) {

		// Determine Jacobian for this integration point
		rhodvol = this->shapeFunction(pts[i][0], pts[i][1]);

		// Element plus material density ... MAY WANT TO REMOVE ELEMENT DENSITY
		rhodvol *= (rhoi[i] * thickness * wts[i]);

		for (int alpha = 0, ia = 0; alpha < 4; alpha++, ia += 6) {
			Nrho = shp[2][alpha] * rhodvol;

			Mel(ia, ia) += Nrho;
			Mel(ia + 1, ia + 1) += Nrho;
			Mel(ia + 2, ia + 2) += Nrho;
		}

	}

	/*// DEBUG
	opserr << "Element Mass Matrix in LOCAL cs: \n";
	for (int i = 0; i < 24; i++) {
	for (int j = 0; j < 24; j++) {
	opserr << Mel(i,j) << "\t";
	}
	opserr << "\n";opserr << "\n";

	} // */

	// Perform coordinate transformation - local to global
	// Mgl = Mgl + (T^ Mel * T);
	Mgl.addMatrixTripleProduct(1.0, T, Mel, 1.0);

	/*// DEBUG
	opserr << "Element Mass Matrix in GLOBAL cs: \n";
	for (int i = 0; i < 24; i++) {
		for (int j = 0; j < 24; j++) {
			opserr << Mgl(i,j) << "\t";
		}
		opserr << "\n";opserr << "\n";

	} // */

	return Mgl; // return element mass matrix in global cs
}

// !!! WORK ON THIS / CHECK (check dimensions) !!! used only if assigned material density
int FourNodeQuadWall3d::addInertiaLoadToUnbalance(const Vector& accel) {

	int i;
	static double rhoi[4];
	double sum = 0.0;
	for (i = 0; i < 4; i++) {
		rhoi[i] = theMaterial[i]->getRho();
		sum += rhoi[i];
	}

	if (sum == 0.0)
		return 0;

	// Get R * accel from the nodes
	const Vector& Raccel1 = theNodes[0]->getRV(accel);
	const Vector& Raccel2 = theNodes[1]->getRV(accel);
	const Vector& Raccel3 = theNodes[2]->getRV(accel);
	const Vector& Raccel4 = theNodes[3]->getRV(accel);

	static double ra[12];

	// Translational accelerations only (no rotational masses)
	ra[0] = Raccel1(0);
	ra[1] = Raccel1(1);
	ra[2] = Raccel1(2);
	ra[3] = Raccel2(0);
	ra[4] = Raccel2(1);
	ra[5] = Raccel2(2);
	ra[6] = Raccel3(0);
	ra[7] = Raccel3(1);
	ra[8] = Raccel3(2);
	ra[9] = Raccel4(0);
	ra[10] = Raccel4(1);
	ra[11] = Raccel4(2);

	// Compute mass matrix
	this->getMass();

	// FIGURE OUT THESE TIMENSIONS !!! check with other elements !!!
	// Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	for (i = 0; i < 12; i++)
		Q(i) += -Kel(i, i) * ra[i];

	return 0;
}


// !!! WORK ON THIS / CHECK (check dimensions) !!!
const Vector& FourNodeQuadWall3d::getResistingForceIncInertia() {

	int i;
	static double rhoi[4];
	double sum = 0.0;
	for (i = 0; i < 4; i++) {
		rhoi[i] = theMaterial[i]->getRho();
		sum += rhoi[i];
	}

	// if no mass terms .. just add damping terms
	if (sum == 0.0) {
		this->getResistingForce();

		// add the damping forces if rayleigh damping
		if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			Pgl += this->getRayleighDampingForces();

		return Pgl;
	}

	// BELOW IS USED ONLY IF MATERIAL DENSITY ASSIGNED - DIMENSIONS ARE NOT CONSISTENT - CODE AND TEST !!!
	const Vector& accel1 = theNodes[0]->getTrialAccel();
	const Vector& accel2 = theNodes[1]->getTrialAccel();
	const Vector& accel3 = theNodes[2]->getTrialAccel();
	const Vector& accel4 = theNodes[3]->getTrialAccel();

	static double a[12];

	a[0] = accel1(0);
	a[1] = accel1(1);
	a[2] = accel1(2);
	a[3] = accel2(0);
	a[4] = accel2(1);
	a[5] = accel2(2);
	a[6] = accel3(0);
	a[7] = accel3(1);
	a[8] = accel3(2);
	a[9] = accel4(0);
	a[10] = accel4(1);
	a[11] = accel4(2);

	// Compute the current resisting force
	this->getResistingForce();

	// Compute the mass matrix
	this->getMass();

	// Take advantage of lumped mass matrix
	for (i = 0; i < 12; i++)
		Pel(i) += Kel(i, i) * a[i];

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
		Pel += this->getRayleighDampingForces();

	return Pel;
}


/* FOR RECORDERS *//////////////////////////////////////////////////////////////////////////
// Define recorders - OK by KK
Response* FourNodeQuadWall3d::setResponse(const char** argv, int argc, OPS_Stream& output) {

	Response* theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", "FourNodeQuadWall3d");
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);
	output.attr("node3", connectedExternalNodes[2]);
	output.attr("node4", connectedExternalNodes[3]);

	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) { // 1: Nodal forces in local cs - needs to be tested !!!

		char dataOut[26];
		for (int i = 1; i <= 4; i++) {
			sprintf(dataOut, "P1_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "P2_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "P3_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "M1_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "M2_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "M3_%d", i);
			output.tag("ResponseType", dataOut);
		}

		theResponse = new ElementResponse(this, 1, Vector(24));
	}

	else if (strcmp(argv[0], "force_global") == 0 || strcmp(argv[0], "forces_global") == 0) { // 2: Nodal forces in global cs - needs to be tested !!!

		char dataOut[26];
		for (int i = 1; i <= 4; i++) {
			sprintf(dataOut, "P1_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "P2_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "P3_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "M1_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "M2_%d", i);
			output.tag("ResponseType", dataOut);
			sprintf(dataOut, "M3_%d", i);
			output.tag("ResponseType", dataOut);
		}

		theResponse = new ElementResponse(this, 2, Vector(24));
	}

	else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) { // 3: nD material response at a integration point in local cs - needs to be tested !!!

		int pointNum = atoi(argv[1]);
		if (pointNum > 0 && pointNum <= 4) {

			output.tag("GaussPoint");
			output.attr("number", pointNum);
			output.attr("eta", pts[pointNum - 1][0]);
			output.attr("neta", pts[pointNum - 1][1]);

			theResponse = theMaterial[pointNum - 1]->setResponse(&argv[2], argc - 2, output);

			output.endTag();

		}
	}

	else if ((strcmp(argv[0], "strains") == 0) || (strcmp(argv[0], "strain") == 0)) { // 4: strains at all integration points in local cs

		for (int i = 0; i < 4; i++) {
			output.tag("GaussPoint");
			output.attr("number", i + 1);
			output.attr("eta", pts[i][0]);
			output.attr("neta", pts[i][1]);

			output.tag("NdMaterialOutput");
			output.attr("classType", theMaterial[i]->getClassTag());
			output.attr("tag", theMaterial[i]->getTag());

			output.tag("ResponseType", "eps11");
			output.tag("ResponseType", "eps22");
			output.tag("ResponseType", "gamma12");

			output.endTag(); // GaussPoint
			output.endTag(); // NdMaterialOutput
		}

		theResponse = new ElementResponse(this, 4, Vector(12));
	}

	else if ((strcmp(argv[0], "average_strains") == 0) || (strcmp(argv[0], "ave_strains") == 0)) { // 5: average strains from all integration points in local cs

		output.tag("NdMaterialOutput");

		output.tag("ResponseType", "eps11");
		output.tag("ResponseType", "eps22");
		output.tag("ResponseType", "gamma12");

		theResponse = new ElementResponse(this, 5, Vector(3));
	}

	else if ((strcmp(argv[0], "average_strains_global") == 0) || (strcmp(argv[0], "ave_strains_global") == 0)) { // 6: average strains from all integration points in global cs

		output.tag("NdMaterialOutput");

		output.tag("ResponseType", "EPS11");
		output.tag("ResponseType", "EPS22");
		output.tag("ResponseType", "GAMMA12");

		theResponse = new ElementResponse(this, 6, Vector(6));
	}

	else if ((strcmp(argv[0], "stresses") == 0) || (strcmp(argv[0], "stress") == 0)) { // 7: stresses at all integration points in local cs

		for (int i = 0; i < 4; i++) {
			output.tag("GaussPoint");
			output.attr("number", i + 1);
			output.attr("eta", pts[i][0]);
			output.attr("neta", pts[i][1]);

			output.tag("NdMaterialOutput");
			output.attr("classType", theMaterial[i]->getClassTag());
			output.attr("tag", theMaterial[i]->getTag());

			output.tag("ResponseType", "sigma11");
			output.tag("ResponseType", "sigma22");
			output.tag("ResponseType", "sigma12");

			output.endTag(); // GaussPoint
			output.endTag(); // NdMaterialOutput
		}

		theResponse = new ElementResponse(this, 7, Vector(12));
	}

	else if ((strcmp(argv[0], "average_stresses") == 0) || (strcmp(argv[0], "ave_stress") == 0)) { // 8: average stresses from all integration points in local cs

		output.tag("NdMaterialOutput");

		output.tag("ResponseType", "sigma11");
		output.tag("ResponseType", "sigma22");
		output.tag("ResponseType", "tau12");

		theResponse = new ElementResponse(this, 8, Vector(3));
	}

	else if ((strcmp(argv[0], "average_stresses_global") == 0) || (strcmp(argv[0], "ave_stress_global") == 0)) { // 9: average stresses from all integration points in global cs

		output.tag("NdMaterialOutput");

		output.tag("ResponseType", "SIGMA11");
		output.tag("ResponseType", "SIGMA22");
		output.tag("ResponseType", "TAU12");

		theResponse = new ElementResponse(this, 9, Vector(6));
	}

	output.endTag(); // ElementOutput

	return theResponse;
}

// Record values - OK by KK
int  FourNodeQuadWall3d::getResponse(int responseID, Information& eleInfo) {

	if (responseID == 1) { // 1: Nodal forces in local cs - needs to be tested !!!

		return eleInfo.setVector(Pel);

	}
	else if (responseID == 2) { // 2: Nodal forces in global cs - needs to be tested !!!

		return eleInfo.setVector(this->getResistingForce());

	}
	else if (responseID == 4) { // 4: strains at all integration points in local cs

	 // Loop over the integration points
		static Vector strains(12);
		int cnt = 0;

		for (int i = 0; i < 4; i++) {

			// Get material stress response
			const Vector& epsilon = theMaterial[i]->getStrain();
			strains(cnt) = epsilon(0);
			strains(cnt + 1) = epsilon(1);
			strains(cnt + 2) = epsilon(2);
			cnt += 3;
		}
		return eleInfo.setVector(strains);

	}
	else if (responseID == 5) { // 5: average strains from all integration points in local cs

	 // Loop over the integration points
		static Vector avestrains(3);

		double Dx, Dy, Dxy;

		Dx = 0.0;
		Dy = 0.0;
		Dxy = 0.0;

		for (int i = 0; i < 4; i++) {

			// Get material strain response
			const Vector& epsilon = theMaterial[i]->getStrain();

			Dx += epsilon(0);
			Dy += epsilon(1);
			Dxy += epsilon(2);

		}

		avestrains(0) = Dx / 4.0;
		avestrains(1) = Dy / 4.0;
		avestrains(2) = Dxy / 4.0;

		return eleInfo.setVector(avestrains);

	}
	else if (responseID == 6) { // 6: average strains from all integration points in global cs

	 // Loop over the integration points
		static Vector avestrainslocal(3);
		static Vector avestrainsglobal(6);
		avestrainsglobal.Zero();

		double Dx, Dy, Dxy;

		Dx = 0.0;
		Dy = 0.0;
		Dxy = 0.0;

		for (int i = 0; i < 4; i++) {

			// Get material strain response
			const Vector& epsilon = theMaterial[i]->getStrain();

			Dx += epsilon(0);
			Dy += epsilon(1);
			Dxy += epsilon(2);

		}

		avestrainslocal(0) = Dx / 4.0;
		avestrainslocal(1) = Dy / 4.0;
		avestrainslocal(2) = Dxy / 4.0;

		/*// DEBUG
		opserr << "EpsXLOC = " << avestrainslocal(0) << "\n";
		opserr << "EpsYLOC = " << avestrainslocal(1) << "\n";
		opserr << "EpsXYLOC = " << avestrainslocal(2) << "\n\n";
		// */

		avestrainsglobal.addMatrixVector(1.0, Ttstrain, avestrainslocal, 1.0);

		/*// DEBUG
		opserr << "EpsXGL = " << avestrainsglobal(0) << "\n";
		opserr << "EpsYGL = " << avestrainsglobal(1) << "\n";
		opserr << "EpsZGL = " << avestrainsglobal(2) << "\n";
		opserr << "EpsXYGL = " << avestrainsglobal(3) << "\n";
		opserr << "EpsYZGL = " << avestrainsglobal(4) << "\n";
		opserr << "EpsZXGL = " << avestrainsglobal(5) << "\n\n\n";
		// */

		return eleInfo.setVector(avestrainsglobal);

	}
	else if (responseID == 7) { // 7: stresses at all integration points in local cs

	 // Loop over the integration points
		static Vector stresses(12);
		int cnt = 0;
		for (int i = 0; i < 4; i++) {

			// Get material stress response
			const Vector& sigma = theMaterial[i]->getStress();
			stresses(cnt) = sigma(0);
			stresses(cnt + 1) = sigma(1);
			stresses(cnt + 2) = sigma(2);
			cnt += 3;
		}
		return eleInfo.setVector(stresses);

	}
	else if (responseID == 8) { // 8: average stresses from all integration points in local cs

// Loop over the integration points
		static Vector avestress(3);

		double Sx, Sy, Sxy;

		Sx = 0.0;
		Sy = 0.0;
		Sxy = 0.0;

		for (int i = 0; i < 4; i++) {

			// Get material stress response
			const Vector& epsilon = theMaterial[i]->getStress();

			Sx += epsilon(0);
			Sy += epsilon(1);
			Sxy += epsilon(2);

		}

		avestress(0) = Sx / 4.0;
		avestress(1) = Sy / 4.0;
		avestress(2) = Sxy / 4.0;

		return eleInfo.setVector(avestress);

	}
	else if (responseID == 9) { // 9: average stresses from all integration points in global cs

	 // Loop over the integration points
		static Vector avestresslocal(3);
		static Vector avestressglobal(6);
		avestressglobal.Zero();

		double Sx, Sy, Sxy;

		Sx = 0.0;
		Sy = 0.0;
		Sxy = 0.0;

		for (int i = 0; i < 4; i++) {

			// Get material stress response
			const Vector& epsilon = theMaterial[i]->getStress();

			Sx += epsilon(0);
			Sy += epsilon(1);
			Sxy += epsilon(2);

		}

		avestresslocal(0) = Sx / 4.0;
		avestresslocal(1) = Sy / 4.0;
		avestresslocal(2) = Sxy / 4.0;

		avestressglobal.addMatrixVector(1.0, Ttstrain, avestresslocal, 1.0);

		return eleInfo.setVector(avestressglobal);

	}

	else

		return -1;
}


/* FOR PARALLEL COMPUTING *//////////////////////////////////////////////////////////////////////////
// !!! WORK ON THIS / CHECK !!!
int FourNodeQuadWall3d::sendSelf(int commitTag, Channel& theChannel) {

	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// Quad packs its data into a Vector and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
	static Vector data(10);
	data(0) = this->getTag();
	data(1) = thickness;
	data(3) = b[0];
	data(4) = b[1];
	data(5) = pressure;

	data(6) = alphaM;
	data(7) = betaK;
	data(8) = betaK0;
	data(9) = betaKc;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING FourNodeQuadWall3d::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}


	// Now quad sends the ids of its materials
	int matDbTag;

	static ID idData(12);

	int i;
	for (i = 0; i < 4; i++) {
		idData(i) = theMaterial[i]->getClassTag();
		matDbTag = theMaterial[i]->getDbTag();
		// NOTE: we do have to ensure that the material has a database
		// tag if we are sending to a database channel.
		if (matDbTag == 0) {
			matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
				theMaterial[i]->setDbTag(matDbTag);
		}
		idData(i + 4) = matDbTag;
	}

	idData(8) = connectedExternalNodes(0);
	idData(9) = connectedExternalNodes(1);
	idData(10) = connectedExternalNodes(2);
	idData(11) = connectedExternalNodes(3);

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING FourNodeQuadWall3d::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	// Finally, quad asks its material objects to send themselves
	for (i = 0; i < 4; i++) {
		res += theMaterial[i]->sendSelf(commitTag, theChannel);
		if (res < 0) {
			opserr << "WARNING FourNodeQuadWall3d::sendSelf() - " << this->getTag() << " failed to send its Material\n";
			return res;
		}
	}

	return res;
}

// !!! WORK ON THIS / CHECK !!!
int FourNodeQuadWall3d::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {

	int res = 0;

	int dataTag = this->getDbTag();

	// Quad creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(10);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING FourNodeQuadWall3d::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	thickness = data(1);
	b[0] = data(3);
	b[1] = data(4);
	pressure = data(5);

	alphaM = data(6);
	betaK = data(7);
	betaK0 = data(8);
	betaKc = data(9);

	static ID idData(12);
	// Quad now receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING FourNodeQuadWall3d::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	connectedExternalNodes(0) = idData(8);
	connectedExternalNodes(1) = idData(9);
	connectedExternalNodes(2) = idData(10);
	connectedExternalNodes(3) = idData(11);


	if (theMaterial == 0) {
		// Allocate new materials
		theMaterial = new NDMaterial * [4];
		if (theMaterial == 0) {
			opserr << "FourNodeQuadWall3d::recvSelf() - Could not allocate NDMaterial* array\n";
			return -1;
		}
		for (int i = 0; i < 4; i++) {
			int matClassTag = idData(i);
			int matDbTag = idData(i + 4);
			// Allocate new material with the sent class tag
			theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
			if (theMaterial[i] == 0) {
				opserr << "FourNodeQuadWall3d::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
				return -1;
			}
			// Now receive materials into the newly allocated space
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to recv itself\n";
				return res;
			}
		}
	}

	// materials exist , ensure materials of correct type and recvSelf on them
	else {
		for (int i = 0; i < 4; i++) {
			int matClassTag = idData(i);
			int matDbTag = idData(i + 4);
			// Check that material is of the right type; if not,
			// delete it and create a new one of the right type
			if (theMaterial[i]->getClassTag() != matClassTag) {
				delete theMaterial[i];
				theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
				if (theMaterial[i] == 0) {
					opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to create\n";

					return -1;
				}
			}
			// Receive the material
			theMaterial[i]->setDbTag(matDbTag);
			res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "NLBeamColumn3d::recvSelf() - material " << i << "failed to recv itself\n";
				return res;
			}
		}
	}

	return res;
}

/* MISC *//////////////////////////////////////////////////////////////////////////
// Print element state on the screen - OK by KK
void FourNodeQuadWall3d::Print(OPS_Stream& s, int flag) {

	if (flag == 2) {

		s << "QuadWall3d\n";

		int i;
		const int numNodes = 4;
		const int nstress = 3;

		for (i = 0; i < numNodes; i++) {
			const Vector& nodeCrd = theNodes[i]->getCrds();
			const Vector& nodeDisp = theNodes[i]->getDisp();
			s << "Node: " << theNodes[i]->getTag() << ": " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2) << " " << endln;
		}

		// spit out the section location & invoke print on the scetion
		const int numMaterials = 4;

		static Vector avgStress(nstress);
		static Vector avgStrain(nstress);
		avgStress.Zero();
		avgStrain.Zero();
		for (i = 0; i < numMaterials; i++) {
			avgStress += theMaterial[i]->getStress();
			avgStrain += theMaterial[i]->getStrain();
		}
		avgStress /= numMaterials;
		avgStrain /= numMaterials;

		s << "Average_strain ";
		for (i = 0; i < nstress; i++)
			s << avgStress(i) << " ";
		s << endln;

		s << "Average_Stress";
		for (i = 0; i < nstress; i++)
			s << avgStrain(i) << " ";
		s << endln;

	}
	else {
		s << "\nQuadWall3d, Tag:  " << this->getTag() << endln;
		s << "\tConnected external nodes:  " << connectedExternalNodes;
		s << "\tthickness:  " << thickness << endln;
		s << "\tsurface pressure:  " << pressure << endln; // REMOVE IF NOT USED !!!
		s << "\tbody forces:  " << b[0] << " " << b[1] << endln; // REMOVE IF NOT USED !!!
		theMaterial[0]->Print(s, flag);
		s << "\tStrains (xx yy xy)" << endln;
		for (int i = 0; i < 4; i++)
			s << "\t\tGauss point " << i + 1 << ": " << theMaterial[i]->getStrain();
	}
}

// Display elements on the screen - OK by KK - add plotting of cracks !!!
int FourNodeQuadWall3d::displaySelf(Renderer& theViewer, int displayMode, float fact) {

	// first set the quantity to be displayed at the nodes;
	// if displayMode is 1 through 3 we will plot material stresses otherwise 0.0

	static Vector values(4);

	for (int j = 0; j < 4; j++)
		values(j) = 0.0;

	if (displayMode < 4 && displayMode > 0) {
		for (int i = 0; i < 4; i++) {
			const Vector& stress = theMaterial[i]->getStress();
			values(i) = stress(displayMode - 1);
		}
	}

	// now  determine the end points of the quad based on
	// the display factor (a measure of the distorted image)
	// store this information in 4 3d vectors v1 through v4
	const Vector& end1Crd = theNodes[0]->getCrds();
	const Vector& end2Crd = theNodes[1]->getCrds();
	const Vector& end3Crd = theNodes[2]->getCrds();
	const Vector& end4Crd = theNodes[3]->getCrds();

	static Matrix coords(4, 3);

	if (displayMode >= 0) {

		const Vector& end1Disp = theNodes[0]->getDisp();
		const Vector& end2Disp = theNodes[1]->getDisp();
		const Vector& end3Disp = theNodes[2]->getDisp();
		const Vector& end4Disp = theNodes[3]->getDisp();

		for (int i = 0; i < 3; i++) {
			coords(0, i) = end1Crd(i) + end1Disp(i) * fact;
			coords(1, i) = end2Crd(i) + end2Disp(i) * fact;
			coords(2, i) = end3Crd(i) + end3Disp(i) * fact;
			coords(3, i) = end4Crd(i) + end4Disp(i) * fact;
		}
	}
	else {
		int mode = displayMode * -1;
		const Matrix& eigen1 = theNodes[0]->getEigenvectors();
		const Matrix& eigen2 = theNodes[1]->getEigenvectors();
		const Matrix& eigen3 = theNodes[2]->getEigenvectors();
		const Matrix& eigen4 = theNodes[3]->getEigenvectors();
		if (eigen1.noCols() >= mode) {
			for (int i = 0; i < 3; i++) {
				coords(0, i) = end1Crd(i) + eigen1(i, mode - 1) * fact;
				coords(1, i) = end2Crd(i) + eigen2(i, mode - 1) * fact;
				coords(2, i) = end3Crd(i) + eigen3(i, mode - 1) * fact;
				coords(3, i) = end4Crd(i) + eigen4(i, mode - 1) * fact;
			}
		}
		else {
			for (int i = 0; i < 3; i++) {
				coords(0, i) = end1Crd(i);
				coords(1, i) = end2Crd(i);
				coords(2, i) = end3Crd(i);
				coords(3, i) = end4Crd(i);
			}
		}
	}

	int error = 0;

	// finally we draw the element using drawPolygon
	error += theViewer.drawPolygon(coords, values);

	return error;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// IMPLEMENT INTO FINAL VERSION /////////////////////////////////////////
////////////////////////////////////// think about pros/cons /////////////////////////////////////////////

/* FOR BODY FORCES AND PRESSURE LOAD - NOT USED SO FAR - CONSIDER REMOVING *//////////////////////////////
// No body force applied - NOT USED SO FAR !!!
void FourNodeQuadWall3d::zeroLoad(void) {

	Q.Zero();

	applyLoad = 0;

	appliedB[0] = 0.0;
	appliedB[1] = 0.0;

	return;
}

// Apply body forces - NOT USED SO FAR !!!
int  FourNodeQuadWall3d::addLoad(ElementalLoad* theLoad, double loadFactor) {

	// Added option for applying body forces in load pattern: C.McGann, U.Washington
	int type;
	const Vector& data = theLoad->getData(type, loadFactor);

	if (type == LOAD_TAG_SelfWeight) {
		applyLoad = 1;
		appliedB[0] += loadFactor * data(0) * b[0];
		appliedB[1] += loadFactor * data(1) * b[1];
		return 0;
	}
	else {
		opserr << "FourNodeQuadWall3d::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
		return -1;
	}

	return -1;
}

// Apply pressure load - NOT USED SO FAR !!!
void  FourNodeQuadWall3d::setPressureLoadAtNodes(void) {

	pressureLoad.Zero();

	if (pressure == 0.0)
		return;

	const Vector& node1 = theNodes[0]->getCrds();
	const Vector& node2 = theNodes[1]->getCrds();
	const Vector& node3 = theNodes[2]->getCrds();
	const Vector& node4 = theNodes[3]->getCrds();

	double x1 = node1(dirn[0]);
	double y1 = node1(dirn[1]);
	double x2 = node2(dirn[0]);
	double y2 = node2(dirn[1]);
	double x3 = node3(dirn[0]);
	double y3 = node3(dirn[1]);
	double x4 = node4(dirn[0]);
	double y4 = node4(dirn[1]);

	double dx12 = x2 - x1;
	double dy12 = y2 - y1;
	double dx23 = x3 - x2;
	double dy23 = y3 - y2;
	double dx34 = x4 - x3;
	double dy34 = y4 - y3;
	double dx41 = x1 - x4;
	double dy41 = y1 - y4;

	double pressureOver2 = pressure / 2.0;

	// Contribution from side 12
	pressureLoad(dirn[0]) += pressureOver2 * dy12;
	pressureLoad(dirn[0] + 3) += pressureOver2 * dy12;
	pressureLoad(dirn[1]) += pressureOver2 * -dx12;
	pressureLoad(dirn[1] + 3) += pressureOver2 * -dx12;

	// Contribution from side 23
	pressureLoad(dirn[0] + 3) += pressureOver2 * dy23;
	pressureLoad(dirn[0] + 6) += pressureOver2 * dy23;
	pressureLoad(dirn[1] + 3) += pressureOver2 * -dx23;
	pressureLoad(dirn[1] + 6) += pressureOver2 * -dx23;

	// Contribution from side 34
	pressureLoad(dirn[0] + 6) += pressureOver2 * dy34;
	pressureLoad(dirn[0] + 9) += pressureOver2 * dy34;
	pressureLoad(dirn[1] + 6) += pressureOver2 * -dx34;
	pressureLoad(dirn[1] + 9) += pressureOver2 * -dx34;

	// Contribution from side 41
	pressureLoad(dirn[0] + 9) += pressureOver2 * dy41;
	pressureLoad(dirn[0]) += pressureOver2 * dy41;
	pressureLoad(dirn[1] + 9) += pressureOver2 * -dx41;
	pressureLoad(dirn[1]) += pressureOver2 * -dx41;

	//pressureLoad = pressureLoad*thickness;
}

/* FOR MONTE CARLO SIMULATIONS *//////////////////////////////////////////////////////////////////////////
// !!! WORK ON THIS / CHECK !!!
int FourNodeQuadWall3d::setParameter(const char** argv, int argc, Parameter& param) {

	if (argc < 1)
		return -1;

	int res = -1;

	// quad pressure loading
	if (strcmp(argv[0], "pressure") == 0) {
		return param.addObject(2, this);
	}
	// a material parameter
	else if ((strstr(argv[0], "material") != 0) && (strcmp(argv[0], "materialState") != 0)) {

		if (argc < 3)
			return -1;

		int pointNum = atoi(argv[1]);
		if (pointNum > 0 && pointNum <= 4)
			return theMaterial[pointNum - 1]->setParameter(&argv[2], argc - 2, param);
		else
			return -1;
	}

	// otherwise it could be just a forall material parameter
	else {

		int matRes = res;
		for (int i = 0; i < 4; i++) {

			matRes = theMaterial[i]->setParameter(argv, argc, param);

			if (matRes != -1)
				res = matRes;
		}
	}

	return res;
}

// !!! WORK ON THIS / CHECK !!!
int FourNodeQuadWall3d::updateParameter(int parameterID, Information& info) {

	int res = -1;
	int matRes = res;
	switch (parameterID) {
	case -1:
		return -1;

	case 1:

		for (int i = 0; i < 4; i++) {
			matRes = theMaterial[i]->updateParameter(parameterID, info);
		}
		if (matRes != -1) {
			res = matRes;
		}
		return res;

	case 2:
		pressure = info.theDouble;
		this->setPressureLoadAtNodes();	// update consistent nodal loads
		return 0;

	default:
		/*
		if (parameterID >= 100) { // material parameter
		int pointNum = parameterID/100;
		if (pointNum > 0 && pointNum <= 4)
		return theMaterial[pointNum-1]->updateParameter(parameterID-100*pointNum, info);
		else
		return -1;
		} else // unknown
		*/
		return -1;
	}
}