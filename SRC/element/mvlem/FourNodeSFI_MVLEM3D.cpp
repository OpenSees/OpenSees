// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								Kamiar Kalbasi
//								California State University, Fullerton 
//
// Created:
//
// Description: 
//
// References:
//
// Source: /usr/local/cvs/OpenSees/SRC/element/FourNodeMVLEM3D/FourNodeMVLEM3D.h
//
// Rev: 1

#include <elementAPI.h>
#include <G3Globals.h>
#include <UniaxialMaterial.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>
#include "FourNodeSFI_MVLEM3D.h"
#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <NDMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <DummyStream.h>

#include <TclModelBuilder.h> // for creating/adding internal nodes (theNodesX) to the domain 

// Typical constructor
FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D(int tag,
	int Nd1, int Nd2, int Nd3, int Nd4,
	Node * theNd1, Node * theNd2, Node * theNd3, Node * theNd4,
	NDMaterial **materials,
	double *thickness,
	double *width,
	Domain *theTclDomain,
	int mm = 0,
	double cc = 0.0,
	double nn = 0.0,
	double tf = 0.0)

	:Element(tag, ELE_TAG_FourNodeSFI_MVLEM3D),
	externalNodes(4 + mm),
	theNd1(0),
	theNd2(0),
	theNd3(0),
	theNd4(0),
	theNodesX(0),
	theNodesALL(0),
	theMaterial(0), theLoad(0),
	FourNodeSFI_MVLEM3DStrainX(0), FourNodeSFI_MVLEM3DStrainY(0), FourNodeSFI_MVLEM3DStrainXY(0), FourNodeSFI_MVLEM3DStrain(0),
	x(0), b(0), AcX(0), AcY(0), kx(0), ky(0), Kh(0), Fx(0), Fy(0), Fxy(0), Dens(0), Dx(0), Dy(0), Dxy(0),
	FourNodeSFI_MVLEM3DK(24 + m, 24 + m), FourNodeSFI_MVLEM3DR(24 + m), FourNodeSFI_MVLEM3DD(24 + m, 24 + m), FourNodeSFI_MVLEM3DM(24 + m, 24 + m),
	FourNodeSFI_MVLEM3DKlocal(24 + m, 24 + m), FourNodeSFI_MVLEM3DDlocal(24 + m, 24 + m), FourNodeSFI_MVLEM3DRlocal(24 + m), FourNodeSFI_MVLEM3DMlocal(24 + m, 24 + m),
	P_24DOF(24), P_24DOF_local(24),
	m(mm), c(cc), NUelastic(nn), Tfactor(tf),
	T(24 + m, 24 + m), Tt(3, 3), T6(6, 6),
	nd1Crds(3), nd2Crds(3), nd3Crds(3), nd4Crds(3), modifiedT(0), t(0)
{
	
	TotalMass = 0.0;
	NodeMass = 0.0;
	h = 0.0;
	d = 0.0;
	Lw = 0.0;
	K_drilling = 0.0;

	// Out of Plane parameters
	Eave = 0.0;
	Tave = 0.0;

	// Imaginary beam properties
	Eim = 0.0;
	Him = 0.0;
	Iim = 0.0;
	Aim = 0.0;

	// Check number of fibers - max is 999 to avoid overlapping in internal node tags
	if (m > 999) {
		opserr << "WARNING: Number of fibers assigned is " << m << ". Maximum allowed number of fibers is 999!\n";
		exit(-1);
	}

	// Fill in the ID containing external node info with node id's 
	if (externalNodes.Size() != 4 + m)
		opserr << "FATAL FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D() - out of memory, could not create an ID of size 2+m\n";

	externalNodes(0) = Nd1;
	externalNodes(1) = Nd2;
	externalNodes(2) = Nd3;
	externalNodes(3) = Nd4;

	// Create a internal node tag
	for (int i = 0; i < m; i++) { // Large NEGATIVE integer starting with tag of the element
		externalNodes(i + 4) = -(tag * 1000 + i + 1); // Max fibers is 999 to avoid overlap
	}

	// Set external node pointers to NULL - external nodes
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;

	// Allocate memory for the m internal nodes
	theNodesX = new Node*[m];
	theNodesALL = new Node*[m + 4];

	// Set NodeX pointers to NULL - m internal nodes
	for (int i = 0; i < m; i++) {
		theNodesX[i] = 0;
	}

	// Set theNodesALL pointers to NULL - m internal nodes
	for (int i = 0; i < m + 4; i++) {
		theNodesALL[i] = 0;
	}

	// Get coordinates of end nodes
	nd1Crds = theNd1->getCrds();
	nd2Crds = theNd2->getCrds();
	nd3Crds = theNd3->getCrds();
	nd4Crds = theNd4->getCrds();

	// Check thickness and width input
	if (thickness == 0) {
		opserr << "FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D() - "
			<< "Null thickness array passed.\n";
		exit(-1);
	}

	if (width == 0) {
		opserr << "FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D() - "
			<< "Null width array passed.\n";
		exit(-1);
	}

	// Allocate memory for the thickness and width
	// Input parameters
	t = new double[m];
	b = new double[m];

	for (int i = 0; i<m; i++) {
		t[i] = thickness[i];
		b[i] = width[i];
		Lw += b[i];		// Total length of the wall
	}

	// Calculate locations of concrete macro-fibers in the cross-section (centerline - x = 0.0)
	x = new double[m];
	for (int i = 0; i < m; i++)
		x[i] = 0.0;

	for (int i = 0; i < m; i++) {
		double sumb_i = 0.0;
		for (int j = 0; j<i + 1; j++)
			sumb_i += b[j];

		x[i] = (sumb_i - b[i] / 2.0) - Lw / 2.0;
	}

	// Compute coordinate transformation matrix
	setTransformationMatrix();
	
	// Define internal dummy nodes
	// Node coordinates in local coordinate system
	Vector nd1CrdL(3);
	Vector nd2CrdL(3);
	Vector nd3CrdL(3);
	Vector nd4CrdL(3);

	nd1CrdL.addMatrixVector(0.0, Tt, nd1Crds, 1.0);
	nd2CrdL.addMatrixVector(0.0, Tt, nd2Crds, 1.0);
	nd3CrdL.addMatrixVector(0.0, Tt, nd3Crds, 1.0);
	nd4CrdL.addMatrixVector(0.0, Tt, nd4Crds, 1.0);

	// Build m internal nodes (NodesX) and add them to the domain
	for (int i = 0; i < m; i++) {

		int nodeId_temp = externalNodes(i + 4); // Store node tag temporarily 

		// Create coordinates wrt top and bottom element node 
		double xLoc_temp = nd1CrdL(0) + x[i];
		double yLoc_temp = 0.5*(nd1CrdL(1) + nd3CrdL(1)); // Mid-height
		double zLoc_temp = nd1CrdL(2);

		Vector Dummy_LCcs(3);
		Dummy_LCcs(0) = xLoc_temp;
		Dummy_LCcs(1) = yLoc_temp;
		Dummy_LCcs(2) = zLoc_temp;

		Vector Dummy_GLcs(3);
		Dummy_GLcs.addMatrixTransposeVector(0.0, Tt, Dummy_LCcs, 1.0);

		// Create Node and add it to the domain
		Node *theNode = 0;

		theNode = new Node(nodeId_temp, 1, Dummy_GLcs(0), Dummy_GLcs(1), Dummy_GLcs(2)); // create internal node with 1 DOF	

		if (theNode == 0) {
			opserr << "WARNING ran out of memory creating node\n";
			opserr << "node: " << nodeId_temp << " in FourNodeSFI_MVLEM3D." << endln; endln;
			exit(-1);
		}

		if (theTclDomain->addNode(theNode) == false) { // add internal node to the domain
			opserr << "WARNING failed to add node to the domain\n";
			opserr << "node: " << nodeId_temp << " in FourNodeSFI_MVLEM3D." << endln;
			delete theNode; // otherwise memory leak
			exit(-1);
		}
	} // END create/add internal nodes

	// Check material input
	if (materials == 0) {
		opserr << "FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D() - "
			<< "Null material array passed.\n";
		exit(-1);
	}

	// Allocate memory for the ND materials
	theMaterial = new NDMaterial*[m];

	if (theMaterial == 0) {
		opserr << "FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D() - "
			<< "Failed to allocate pointers for uniaxial materials.\n";
		exit(-1);
	}

	// Get copies of the ND materials
	for (int i = 0; i < m; i++) {
		if (materials[i] == 0) {
			opserr << "FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D() - "
				"Null ND material pointer passed.\n";
			exit(-1);
		}

		theMaterial[i] = materials[i]->getCopy();

		if (theMaterial[i] == 0) {
			opserr << "FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D() - "
				<< "Failed to copy ND material.\n";
			exit(-1);
		}
	}

	// Allocate memory for element arrays
	// Area of concrete fibers
	AcX = new double[m];
	AcY = new double[m];

	// Panel stiffness (trial)
	kx = new double[m];
	ky = new double[m];
	//Kh = new double[1];

	// Panel force (trial)
	Fx = new double[m];
	Fy = new double[m];
	Fxy = new double[m];

	// Panel stiffness (trial)
	Dx = new double[m];
	Dy = new double[m];
	Dxy = new double[m];

	// Panel strains
	FourNodeSFI_MVLEM3DStrainX = new double[m];
	FourNodeSFI_MVLEM3DStrainY = new double[m];
	FourNodeSFI_MVLEM3DStrainXY = new double[m];
	FourNodeSFI_MVLEM3DStrain = new double[3 * m];

	// Density
	Dens = new double[m];

	// Assign zero to element arrays
	for (int i = 0; i < m; i++) {

		AcX[i] = 0.0;
		AcY[i] = 0.0;

		kx[i] = 0.0;
		ky[i] = 0.0;

		Fx[i] = 0.0;
		Fy[i] = 0.0;
		Fxy[i] = 0.0;

		Dx[i] = 0.0;
		Dy[i] = 0.0;
		Dxy[i] = 0.0;

		FourNodeSFI_MVLEM3DStrainX[i] = 0.0;
		FourNodeSFI_MVLEM3DStrainY[i] = 0.0;
		FourNodeSFI_MVLEM3DStrainXY[i] = 0.0;

		FourNodeSFI_MVLEM3DStrain[i] = 0.0;
		FourNodeSFI_MVLEM3DStrain[i + m] = 0.0;
		FourNodeSFI_MVLEM3DStrain[i + 2 * m] = 0.0;

		Dens[i] = 0.0;
	}

	Kh = 0.0;

	//opserr << "  1		Typical constructor\n";

	// Revert to start
	this->revertToStart();
}

// Constructor which should be invoked by an FE_ObjectBroker only
FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D()
	:Element(0, ELE_TAG_FourNodeSFI_MVLEM3D),
	externalNodes(4 + m),
	theNd1(0),
	theNd2(0),
	theNd3(0),
	theNd4(0),
	theNodesX(0),
	theNodesALL(0),
	theMaterial(0), theLoad(0),
	FourNodeSFI_MVLEM3DStrainX(0), FourNodeSFI_MVLEM3DStrainY(0), FourNodeSFI_MVLEM3DStrainXY(0), FourNodeSFI_MVLEM3DStrain(0),
	x(0), b(0), AcX(0), AcY(0), kx(0), ky(0), Kh(0), Fx(0), Fy(0), Fxy(0), Dens(0), Dx(0), Dy(0), Dxy(0),
	FourNodeSFI_MVLEM3DK(24 + m, 24 + m), FourNodeSFI_MVLEM3DR(24 + m), FourNodeSFI_MVLEM3DD(24 + m, 24 + m), FourNodeSFI_MVLEM3DM(24 + m, 24 + m),
	FourNodeSFI_MVLEM3DKlocal(24 + m, 24 + m), FourNodeSFI_MVLEM3DDlocal(24 + m, 24 + m), FourNodeSFI_MVLEM3DRlocal(24 + m), FourNodeSFI_MVLEM3DMlocal(24 + m, 24 + m),
	P_24DOF(24), P_24DOF_local(24),
	m(m), c(c), NUelastic(NUelastic), Tfactor(Tfactor),
	T(24 + m, 24 + m), Tt(3, 3), T6(6, 6),
	nd1Crds(3), nd2Crds(3), nd3Crds(3), nd4Crds(3), modifiedT(0), t(0)
{
	if (externalNodes.Size() != 4 + m)
		opserr << "FATAL FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D() - out of memory, could not create an ID of size 2\n";

	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;

	// Allocate memory for all nodes (internal and external)
	theNodesX = new Node*[m];		// internal nodes
	theNodesALL = new Node*[m + 4];	// all nodes

	// Set NodeX pointers to zero
	for (int i = 0; i < m; i++)
	{
		theNodesX[i] = 0;
	}

	// Set theNodesALL pointers to zero
	for (int i = 0; i < m + 4; i++)
	{
		theNodesALL[i] = 0;
	}

	//opserr << "  2		FE_ObjectBroker constructor\n";

}

//  Destructor - provided to clean up any memory 
FourNodeSFI_MVLEM3D::~FourNodeSFI_MVLEM3D()
{
	// clean up the memory associated with the element, this is
	// memory the FourNodeSFI_MVLEM3D objects allocates and memory allocated 
	// by other objects that the FourNodeSFI_MVLEM3D object is responsible for 
	// cleaning up, i.e. the MaterialObject.

	if (theMaterial != 0) {
		for (int i = 0; i < m; i++)
			if (theMaterial[i] != 0)
				delete theMaterial[i];
		delete[] theMaterial;
	}

	if (theLoad != 0)
		delete theLoad;
	if (x != 0)
		delete x;
	if (b != 0)
		delete b;
	if (AcX != 0)
		delete AcX;
	if (AcY != 0)
		delete AcY;
	if (kx != 0)
		delete kx;
	if (ky != 0)
		delete ky;
	if (Fx != 0)
		delete Fx;
	if (Fy != 0)
		delete Fy;
	if (Fxy != 0)
		delete Fxy;
	if (Dens != 0)
		delete Dens;
	if (Dx != 0)
		delete Dx;
	if (Dy != 0)
		delete Dy;
	if (Dxy != 0)
		delete Dxy;
	if (FourNodeSFI_MVLEM3DStrainX != 0)
		delete FourNodeSFI_MVLEM3DStrainX;
	if (FourNodeSFI_MVLEM3DStrainY != 0)
		delete FourNodeSFI_MVLEM3DStrainY;
	if (FourNodeSFI_MVLEM3DStrainXY != 0)
		delete FourNodeSFI_MVLEM3DStrainXY;
	if (FourNodeSFI_MVLEM3DStrain != 0)
		delete FourNodeSFI_MVLEM3DStrain;
	if (theNodesX != 0)
		delete theNodesX;
	if (theNodesALL != 0)
		delete theNodesALL;
	if (modifiedT != 0)
		delete modifiedT;
	if (t != 0)
		delete t;
	// !!! tripple check that all global arrays generated in constructor are cleaned up here @@@ DONE
	//opserr << "  3		Destructor\n";
}

// Get number of nodes (external + internal)
int FourNodeSFI_MVLEM3D::getNumExternalNodes(void) const
{
	return 4 + m;
	//opserr << "  4		getNumExternalNodes\n";
}

// Get node tags
const ID & FourNodeSFI_MVLEM3D::getExternalNodes(void)
{
	return externalNodes;
	//opserr << "  5		getExternalNodes\n";
}

// Get node pointers
Node ** FourNodeSFI_MVLEM3D::getNodePtrs(void)
{

	// Pack external and internal node pointers into one array
	for (int i = 0; i < 4; i++) {
		theNodesALL[i] = theNodes[i];
	}

	for (int i = 4; i < m + 4; i++) {
		theNodesALL[i] = theNodesX[i - 4];
	}
	//opserr << "  6		getNodePtrs\n";
	return theNodesALL;
}

// Get number of DOFs 
int FourNodeSFI_MVLEM3D::getNumDOF(void) {

	int NumDOF = 24 + m; // 6 DOFs per external node x 4 external nodes x 1 DOF per m internal nodes
	//opserr << "  7		getNumDOF\n";
	return NumDOF;
}

// Set Domain
void FourNodeSFI_MVLEM3D::setDomain(Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0)
	{
		return;
	}

	// First ensure external nodes exist in Domain and set the node pointers
	Node *end1Ptr, *end2Ptr, *end3Ptr, *end4Ptr;

	int Nd1 = externalNodes(0);
	int Nd2 = externalNodes(1);
	int Nd3 = externalNodes(2);
	int Nd4 = externalNodes(3);
	end1Ptr = theDomain->getNode(Nd1);
	end2Ptr = theDomain->getNode(Nd2);
	end3Ptr = theDomain->getNode(Nd3);
	end4Ptr = theDomain->getNode(Nd4);

	theNodes[0] = end1Ptr;
	theNodes[1] = end2Ptr;
	theNodes[2] = end3Ptr;
	theNodes[3] = end4Ptr;

	if (theNodes[0] == 0) {
		opserr << "WARNING FourNodeSFI_MVLEM3D::setDomain() - at FourNodeSFI_MVLEM3D " << this->getTag() << " node " <<
			Nd1 << " does not exist in domain\n";
		return;  // Don't go any further - otherwise segemntation fault
	}

	if (theNodes[1] == 0) {
		opserr << "WARNING FourNodeSFI_MVLEM3D::setDomain() - at FourNodeSFI_MVLEM3D " << this->getTag() << " node " <<
			Nd2 << " does not exist in domain\n";
		return;
	}

	if (theNodes[2] == 0) {
		opserr << "WARNING FourNodeSFI_MVLEM3D::setDomain() - at FourNodeSFI_MVLEM3D " << this->getTag() << " node " <<
			Nd3 << " does not exist in domain\n";
		return;  // Don't go any further - otherwise segemntation fault
	}
	if (theNodes[3] == 0) {
		opserr << "WARNING FourNodeSFI_MVLEM3D::setDomain() - at FourNodeSFI_MVLEM3D " << this->getTag() << " node " <<
			Nd4 << " does not exist in domain\n";
		return;
	}

	// Then ensure internal NodesX exist in Domain and set the node pointers
	for (int i = 0; i<m; i++) {

		int NdX_temp1 = externalNodes(i + 4);

		theNodesX[i] = theDomain->getNode(NdX_temp1);

		if (theNodesX[i] == 0) {
			opserr << "WARNING FourNodeSFI_MVLEM3D::setDomain() - at FourNodeSFI_MVLEM3D " << this->getTag() << " node " <<
				NdX_temp1 << " does not exist in domain\n";
			return;  // Don't go any further - otherwise segemntation fault
		}
	}

	// Call the DomainComponent class method 
	this->DomainComponent::setDomain(theDomain);

	// Ensure conected nodes have correct number of dof's
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();
	int dofNd3 = theNodes[2]->getNumberDOF();
	int dofNd4 = theNodes[3]->getNumberDOF();

	if ((dofNd1 != 6) || (dofNd2 != 6) || (dofNd3 != 6) || (dofNd4 != 6)) {
		opserr << "FourNodeSFI_MVLEM3D::setDomain(): 6 dof required at all nodes. " << dofNd1 << " provided at node 1, " << dofNd2 << " provided at node 2, "
			<< dofNd3 << " provided at node 4, " << dofNd4 << " provided at node 3";
	}

	for (int i = 0; i < m; i++) {
		int dofNdXi = theNodesX[i]->getNumberDOF();
		if (dofNdXi != 1)				// 1 DOF at internal nodes 
		{
			opserr << "FourNodeSFI_MVLEM3D::setDomain(): 1 dof required at internal nodes, " << dofNdXi << " provided\n";
		}
	}

	// Compute coordinate transformation matrix
	//setTransformationMatrix();

	// Calculate the element height and perform checks
	double h1 = pow(pow(nd3Crds(0) - nd1Crds(0), 2.0) + pow(nd3Crds(1) - nd1Crds(1), 2.0) + pow(nd3Crds(2) - nd1Crds(2), 2.0), 0.5);
	double h2 = pow(pow(nd4Crds(0) - nd2Crds(0), 2.0) + pow(nd4Crds(1) - nd2Crds(1), 2.0) + pow(nd4Crds(2) - nd2Crds(2), 2.0), 0.5);

	// Check if element height is zero
	if ((h1 == 0.0) || (h2 == 0.0)) {
		opserr << "WARNING: One of the element sides is ZERO. Check geometry";
		exit(-1);
	}

	// Check if element has constant height
	if ((h1 / h2 > 1.01) || (h1 / h2 < 0.99)) {
		opserr << "WARNING: Element does not have constant height. Check geometry.";
		exit(-1);
	}

	// Element height
	h = (h1 + h2) / 2.0;

	// Calculate average wall thickness
	for (int i = 0; i<m; i++) {
		Tave += t[i] * b[i] / Lw;
	}

	// Calculate the element width and perform checks
	double b1 = pow(pow(nd1Crds(0) - nd2Crds(0), 2.0) + pow(nd1Crds(1) - nd2Crds(1), 2.0) + pow(nd1Crds(2) - nd2Crds(2), 2.0), 0.5);
	double b2 = pow(pow(nd4Crds(0) - nd3Crds(0), 2.0) + pow(nd4Crds(1) - nd3Crds(1), 2.0) + pow(nd4Crds(2) - nd3Crds(2), 2.0), 0.5);

	// Check width of element
	if ((Lw / b1 > 1.01) || (Lw / b1 < 0.99)) {
		opserr << "WARNING: Element nodes coordinates are not matched with fibers width. Check geometry.";
		exit(-1);
	}

	if ((Lw / b2 > 1.01) || (Lw / b2 < 0.99)) {
		opserr << "WARNING: Element nodes coordinates are not matched with fibers width. Check geometry.";
		exit(-1);
	}

	// Calculate distance of corner nodes from center axis
	d = Lw / 2.0;

	// Calculate concrete areas in X and Y directions
	double A = 0.0;
	for (int i = 0; i < m; i++) {
		AcX[i] = h * t[i];
		AcY[i] = b[i] * t[i];

		A += AcY[i];
	}

	// Get panel density from 2-D materials
	for (int i = 0; i < m; i++) {
		Dens[i] = theMaterial[i]->getRho();
	}

	// Calculate the nodal mass (external nodes only) for lumped mass approach
	for (int i = 0; i < m; i++) {
		TotalMass += Dens[i] * AcY[i] * h;
	}

	NodeMass = TotalMass / 4.0;

	// Get Concrete Young's Modulus
	theResponses = new Response *[1];
	if (theResponses == 0) {
		opserr << " FourNodeSFI_MVLEM3D::FourNodeSFI_MVLEM3D - failed allocate responses array\n";
		exit(-1);
	}

	OPS_Stream *theDummyStream = new DummyStream();
	const char **argv = new const char *[1];

	argv[0] = "getInputParameters"; // to get input parameters from concrete material
	for (int i = 0; i < m; i++)
	{
		theResponses[0] = theMaterial[i]->setResponse(argv, 1, *theDummyStream);

		if (theResponses[0] == 0) {
			opserr << " FSAM::FSAM - failed to set appropriate materials tag: " << this->getTag() << "\n";
			exit(-1);
		}

		// Get FSAM material input variables
		theResponses[0]->getResponse();
		Information &theInfoInput = theResponses[0]->getInformation();
		const Vector InputNDMat = theInfoInput.getData();

		Vector InputNDMaterial(InputNDMat.Size());

		for (int j = 0; j < InputNDMat.Size(); j++)
			InputNDMaterial[j] = InputNDMat[j];

		// Calculate out-of-plane modulus of elasticity (average modulus)
		Eave += AcY[i] * InputNDMaterial[9] / A; // !!! FIX THIS @@@ DONE, look at line 422

		// Drilling DOF
		K_drilling += ((InputNDMaterial[9] * AcY[i] * (1.0 - InputNDMaterial[4]) + InputNDMaterial[10] * AcY[i] * InputNDMaterial[4]) / h) * (0.5 * Lw + x[i]) * (0.5 + x[i] / Lw) * (0.5 * Lw + x[i]) * (0.5 + x[i] / Lw);
	}
	
	//opserr << "K_drilling: "<< K_drilling << "\n";
	// Imaginary Beams Calculation
	Eim = 1.0e4 * Eave;
	Him = 0.5 * h; // !!! FIX THIS - USE SOMTHING MEANINGFUL - Change to Hwall/2 = Him @@@ DONE
	Aim = Tave * Him;
	Iim = Tave * Him * Him * Him / 12.0;
	//Eim = K_drilling * Lw / (4.0 * Iim);

	//opserr << "Eim: " << Eim << "\n";

	// Create a vector to hop applied loads - NOT used in the current model formulation (no element loads)
	if (theLoad == 0)
		theLoad = new Vector(24 + m);
	if (theLoad == 0) {
		opserr << "FourNodeSFI_MVLEM3D::setDomain() - element: " << this->getTag()
			<< " out of memory creating vector of size: " << 24 + m << endln;
		return;
	}
	//opserr << "  8		setDomain\n";
}

// Commit state of the materials
int FourNodeSFI_MVLEM3D::commitState()
{
	int errCode = 0;

	// Commit material models
	for (int i = 0; i < m; i++) {
		errCode += theMaterial[i]->commitState();
	}

	//opserr << "  9		commitState\n";
	return errCode;
}

// Revert to last commited state (if convergence is not achieved)
int FourNodeSFI_MVLEM3D::revertToLastCommit()
{
	int errCode = 0;

	// Revert material models
	for (int i = 0; i < m; i++) {
		errCode += theMaterial[i]->revertToLastCommit();
	}
	//opserr << "  10		revertToLastCommit\n";
	return errCode;
}

// Revert to start
int FourNodeSFI_MVLEM3D::revertToStart()
{

	int errCode = 0;

	// Revert material models
	for (int i = 0; i < m; i++)
		errCode += theMaterial[i]->revertToStart();

	// Compute initial stiffness
	this->getInitialStiff();
	//opserr << "  11		revertToStart\n";
	return errCode;

}

// Update state
int FourNodeSFI_MVLEM3D::update()
{
	// setTransformationMatrix(); !!! for implementing P-Delta effects - implement eventually 

	// Get the current strain given trial displacements at nodes 
	//FourNodeSFI_MVLEM3DStrain = this->computeCurrentStrain();
	this->computeCurrentStrain();

	// Set the strain in the materials
	int errCode = 0;

	for (int i = 0; i < m; i++) {

		Vector strain(3);

		strain(0) = FourNodeSFI_MVLEM3DStrain[i];
		strain(1) = FourNodeSFI_MVLEM3DStrain[i + m];
		strain(2) = FourNodeSFI_MVLEM3DStrain[i + 2 * m];

		// Set trial response for material models
		errCode += theMaterial[i]->setTrialStrain(strain);

	}
	//opserr << "  12		update\n";
	return errCode;
}

// Get current strains at RC panels (macro-fibers)
double *FourNodeSFI_MVLEM3D::computeCurrentStrain(void)
{

	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();

	// obtain displacements in X direction from the internal nodes
	for (int i = 0; i < m; i++) {
		const Vector &dispXi = theNodesX[i]->getTrialDisp(); // 1 DOF at theNodesX - Vector of size 1
		Dx[i] = dispXi(0); // get displacements in X direction from the nodes
	}

	Vector dispG(24 + m); // Vector of total 24+m displacemets in global coordinates
	//dispG.Zero();
	Vector dispL(24 + m); // Vector of total 24+m displacemets in local coordinates
	//dispL.Zero();
	Vector dispL_inPlan2N(6); // in-plane displacements of equlivalent 2D, 2N SFI model
	//dispL_inPlan2N.Zero();

	// store nodal displacemnts in global vs
	for (int i = 0; i < 6; i++) {
		dispG(i) = disp1(i);
		dispG(i + 6) = disp2(i);
		dispG(i + 12) = disp3(i);
		dispG(i + 18) = disp4(i);
	}

	for (int i = 0; i < m; ++i) {
		dispG(i + 24) = Dx[i];
	}

	// tranform nodal displacements from global to local cs
	dispL.addMatrixVector(0.0, T, dispG, 1.0);

	dispL_inPlan2N(0) = dispL(0) / 2.0 + dispL(6) / 2.0;
	dispL_inPlan2N(1) = dispL(1) / 2.0 + dispL(7) / 2.0;
	dispL_inPlan2N(2) = dispL(5) / (2.0 * (d*d) + 2.0) + dispL(11) / (2.0 * (d*d) + 2.0) - (dispL(1)*d) / (2.0 * (d*d) + 2.0) + (dispL(7)*d) / (2.0 * (d*d) + 2.0);
	dispL_inPlan2N(3) = dispL(12) / 2.0 + dispL(18) / 2.0;
	dispL_inPlan2N(4) = dispL(13) / 2.0 + dispL(19) / 2.0;
	dispL_inPlan2N(5) = dispL(17) / (2.0 * (d*d) + 2.0) + dispL(23) / (2.0 * (d*d) + 2.0) - (dispL(13)*d) / (2.0 * (d*d) + 2.0) + (dispL(19)*d) / (2.0 * (d*d) + 2.0);

	// Deformations at each RC panel (macro-fiber) - MVLEM formulation
	for (int i = 0; i < m; i++) {
		Dy[i] = -dispL_inPlan2N(1) - x[i] * dispL_inPlan2N(2) + dispL_inPlan2N(4) + x[i] * dispL_inPlan2N(5);
		Dxy[i] = dispL_inPlan2N(0) - dispL_inPlan2N(3) - c * h*dispL_inPlan2N(2) - (1.0 - c)*h*dispL_inPlan2N(5);
	}

	Dsh = -Dxy[0]; // Store shear deformations for the recorder

	// Strains at each RC panel (macro-fiber)
	for (int i = 0; i < m; i++) {
		FourNodeSFI_MVLEM3DStrainX[i] = Dx[i] / b[i];
		FourNodeSFI_MVLEM3DStrainY[i] = Dy[i] / h;
		FourNodeSFI_MVLEM3DStrainXY[i] = -Dxy[i] / h;
	}

	// Store strains into a single vector
	for (int i = 0; i < m; i++) {
		FourNodeSFI_MVLEM3DStrain[i] = FourNodeSFI_MVLEM3DStrainX[i];
		FourNodeSFI_MVLEM3DStrain[i + m] = FourNodeSFI_MVLEM3DStrainY[i];
		FourNodeSFI_MVLEM3DStrain[i + 2 * m] = FourNodeSFI_MVLEM3DStrainXY[i];
	}
	//opserr << "  13		computeCurrentStrain\n";
	// Return strain vector
	return FourNodeSFI_MVLEM3DStrain;

}

// Get the element intial element tangent matrix
const Matrix & FourNodeSFI_MVLEM3D::getInitialStiff(void)
{

	//FourNodeSFI_MVLEM3DK.Zero();		// Global stiffness matrix
	//FourNodeSFI_MVLEM3DKlocal.Zero();	// Local stiffness matrix

	Kh = 0.0;

	for (int i = 0; i < m; i++)
	{
		// Get material initial tangent
		const Matrix &D = theMaterial[i]->getInitialTangent();

		double D00 = D(0, 0); double D01 = D(0, 1); double D02 = D(0, 2);
		double D10 = D(1, 0); double D11 = D(1, 1); double D12 = D(1, 2);
		double D20 = D(2, 0); double D21 = D(2, 1); double D22 = D(2, 2);

		kx[i] = D00 * h*t[i] / b[i];
		ky[i] = D11 * b[i] * t[i] / h;
		Kh += D22 * b[i] * t[i] / h;

	}

	//opserr << "K_drilling in initialStiff: " << K_drilling << "\n";
	/*opserr << "x in initialStiff: ";
	opserr << "\n";
	for (int i = 0; i < m; i++) {
		opserr << x[i] << ",";
		opserr << "\t";
	}
	opserr << "\n" << "\n";*/

	// Build the initial stiffness matrix
	double Kv = 0.0; double Km = 0.0; double e = 0.0; // double ex = 0.0;

	for (int i = 0; i<m; ++i)
	{
		Kv += ky[i];
		Km += ky[i] * x[i] * x[i];
		e += ky[i] * x[i];

		FourNodeSFI_MVLEM3DKlocal(24 + i, 24 + i) = kx[i]; // Diagonal terms accounting for horizontal stiffness
	}

	/*opserr << "Initial Stiff: ";
	opserr << "kh: " << Kh << ", " << "km: " << Km << ", " << "kv: " << Kv << ", " << "e: " << e << ";";
	opserr << "\n" << "\n";*/

	// Assemble element stiffness matrix
	FourNodeSFI_MVLEM3DKlocal(0, 0) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(0, 1) = (Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 2) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 3) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 4) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 5) = -(Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 6) = Kh / 4.0 - (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(0, 7) = -(Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 11) = -(Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 12) = -Kh / 4.0;
	FourNodeSFI_MVLEM3DKlocal(0, 13) = -(Kh*d*h*(c - 1)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 17) = (Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 18) = -Kh / 4;
	FourNodeSFI_MVLEM3DKlocal(0, 19) = (Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 23) = (Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(1, 0) = FourNodeSFI_MVLEM3DKlocal(0, 1);
	FourNodeSFI_MVLEM3DKlocal(1, 1) = Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - (d*(e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(1, 2) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 3) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 4) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 5) = (e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(1, 6) = (Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(1, 7) = Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) + (d*(e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) - (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(1, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 11) = (e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(1, 12) = -(Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(1, 13) = (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - Kv / 4.0 + (d*(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(1, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 17) = -(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(1, 18) = -(Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(1, 19) = (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - Kv / 4.0 - (d*(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(1, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 23) = -(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);

	FourNodeSFI_MVLEM3DKlocal(2, 0) = FourNodeSFI_MVLEM3DKlocal(0, 2);
	FourNodeSFI_MVLEM3DKlocal(2, 1) = FourNodeSFI_MVLEM3DKlocal(1, 2);
	FourNodeSFI_MVLEM3DKlocal(2, 2) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 3) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 4) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 5) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 6) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 7) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 9) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 10) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 16) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(3, 0) = FourNodeSFI_MVLEM3DKlocal(0, 3);
	FourNodeSFI_MVLEM3DKlocal(3, 1) = FourNodeSFI_MVLEM3DKlocal(1, 3);
	FourNodeSFI_MVLEM3DKlocal(3, 2) = FourNodeSFI_MVLEM3DKlocal(2, 3);
	FourNodeSFI_MVLEM3DKlocal(3, 3) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 4) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeSFI_MVLEM3DKlocal(3, 5) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 6) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 7) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(4, 0) = FourNodeSFI_MVLEM3DKlocal(0, 4);
	FourNodeSFI_MVLEM3DKlocal(4, 1) = FourNodeSFI_MVLEM3DKlocal(1, 4);
	FourNodeSFI_MVLEM3DKlocal(4, 2) = FourNodeSFI_MVLEM3DKlocal(2, 4);
	FourNodeSFI_MVLEM3DKlocal(4, 3) = FourNodeSFI_MVLEM3DKlocal(3, 4);
	FourNodeSFI_MVLEM3DKlocal(4, 4) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 5) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 6) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 7) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 8) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 10) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 16) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(5, 0) = FourNodeSFI_MVLEM3DKlocal(0, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 1) = FourNodeSFI_MVLEM3DKlocal(1, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 2) = FourNodeSFI_MVLEM3DKlocal(2, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 3) = FourNodeSFI_MVLEM3DKlocal(3, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 4) = FourNodeSFI_MVLEM3DKlocal(4, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 5) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw + 0*K_drilling;
	FourNodeSFI_MVLEM3DKlocal(5, 6) = -(Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 7) = e / (2.0 * (2.0 * (d*d) + 2.0)) + (d*(Km + Kh*(c*c)*(h*h))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(5, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 11) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (2.0 * Eim*Iim) / Lw + 0*K_drilling;
	FourNodeSFI_MVLEM3DKlocal(5, 12) = (Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 13) = (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - e / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 17) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 18) = (Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 19) = -e / (2.0 * (2.0 * (d*d) + 2.0)) - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 23) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(6, 0) = FourNodeSFI_MVLEM3DKlocal(0, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 1) = FourNodeSFI_MVLEM3DKlocal(1, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 2) = FourNodeSFI_MVLEM3DKlocal(2, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 3) = FourNodeSFI_MVLEM3DKlocal(3, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 4) = FourNodeSFI_MVLEM3DKlocal(4, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 5) = FourNodeSFI_MVLEM3DKlocal(5, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 6) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(6, 7) = -(Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 11) = -(Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 12) = -Kh / 4;
	FourNodeSFI_MVLEM3DKlocal(6, 13) = -(Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 17) = (Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 18) = -Kh / 4.0;
	FourNodeSFI_MVLEM3DKlocal(6, 19) = (Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 23) = (Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(7, 0) = FourNodeSFI_MVLEM3DKlocal(0, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 1) = FourNodeSFI_MVLEM3DKlocal(1, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 2) = FourNodeSFI_MVLEM3DKlocal(2, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 3) = FourNodeSFI_MVLEM3DKlocal(3, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 4) = FourNodeSFI_MVLEM3DKlocal(4, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 5) = FourNodeSFI_MVLEM3DKlocal(5, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 6) = FourNodeSFI_MVLEM3DKlocal(6, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 7) = Kv / 4.0 + (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) + (d*(e / 2.0 + (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(7, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 11) = (e / 2.0 + (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) - (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(7, 12) = (Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(7, 13) = (d*(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - Kv / 4.0;
	FourNodeSFI_MVLEM3DKlocal(7, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 17) = -(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(7, 18) = (Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(7, 19) = -Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - (d*(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(7, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 23) = -(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);

	FourNodeSFI_MVLEM3DKlocal(8, 0) = FourNodeSFI_MVLEM3DKlocal(0, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 1) = FourNodeSFI_MVLEM3DKlocal(1, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 2) = FourNodeSFI_MVLEM3DKlocal(2, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 3) = FourNodeSFI_MVLEM3DKlocal(3, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 4) = FourNodeSFI_MVLEM3DKlocal(4, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 5) = FourNodeSFI_MVLEM3DKlocal(5, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 6) = FourNodeSFI_MVLEM3DKlocal(6, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 7) = FourNodeSFI_MVLEM3DKlocal(7, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 8) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 10) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 16) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(9, 0) = FourNodeSFI_MVLEM3DKlocal(0, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 1) = FourNodeSFI_MVLEM3DKlocal(1, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 2) = FourNodeSFI_MVLEM3DKlocal(2, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 3) = FourNodeSFI_MVLEM3DKlocal(3, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 4) = FourNodeSFI_MVLEM3DKlocal(4, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 5) = FourNodeSFI_MVLEM3DKlocal(5, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 6) = FourNodeSFI_MVLEM3DKlocal(6, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 7) = FourNodeSFI_MVLEM3DKlocal(7, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 8) = FourNodeSFI_MVLEM3DKlocal(8, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 10) = -(Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeSFI_MVLEM3DKlocal(9, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(10, 0) = FourNodeSFI_MVLEM3DKlocal(0, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 1) = FourNodeSFI_MVLEM3DKlocal(1, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 2) = FourNodeSFI_MVLEM3DKlocal(2, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 3) = FourNodeSFI_MVLEM3DKlocal(3, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 4) = FourNodeSFI_MVLEM3DKlocal(4, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 5) = FourNodeSFI_MVLEM3DKlocal(5, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 6) = FourNodeSFI_MVLEM3DKlocal(6, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 7) = FourNodeSFI_MVLEM3DKlocal(7, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 8) = FourNodeSFI_MVLEM3DKlocal(8, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 9) = FourNodeSFI_MVLEM3DKlocal(9, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 10) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 22) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(11, 0) = FourNodeSFI_MVLEM3DKlocal(0, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 1) = FourNodeSFI_MVLEM3DKlocal(1, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 2) = FourNodeSFI_MVLEM3DKlocal(2, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 3) = FourNodeSFI_MVLEM3DKlocal(3, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 4) = FourNodeSFI_MVLEM3DKlocal(4, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 5) = FourNodeSFI_MVLEM3DKlocal(5, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 6) = FourNodeSFI_MVLEM3DKlocal(6, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 7) = FourNodeSFI_MVLEM3DKlocal(7, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 8) = FourNodeSFI_MVLEM3DKlocal(8, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 9) = FourNodeSFI_MVLEM3DKlocal(9, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 10) = FourNodeSFI_MVLEM3DKlocal(10, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 11) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw + 0*K_drilling;
	FourNodeSFI_MVLEM3DKlocal(11, 12) = (Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 13) = (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - e / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 17) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 18) = (Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 19) = -e / (2.0 * (2.0 * (d*d) + 2.0)) - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 23) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(12, 0) = FourNodeSFI_MVLEM3DKlocal(0, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 1) = FourNodeSFI_MVLEM3DKlocal(1, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 2) = FourNodeSFI_MVLEM3DKlocal(2, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 3) = FourNodeSFI_MVLEM3DKlocal(3, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 4) = FourNodeSFI_MVLEM3DKlocal(4, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 5) = FourNodeSFI_MVLEM3DKlocal(5, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 6) = FourNodeSFI_MVLEM3DKlocal(6, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 7) = FourNodeSFI_MVLEM3DKlocal(7, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 8) = FourNodeSFI_MVLEM3DKlocal(8, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 9) = FourNodeSFI_MVLEM3DKlocal(9, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 10) = FourNodeSFI_MVLEM3DKlocal(10, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 11) = FourNodeSFI_MVLEM3DKlocal(11, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 12) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(12, 13) = (Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(12, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 17) = -(Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(12, 18) = Kh / 4.0 - (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(12, 19) = -(Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(12, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 23) = -(Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(13, 0) = FourNodeSFI_MVLEM3DKlocal(0, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 1) = FourNodeSFI_MVLEM3DKlocal(1, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 2) = FourNodeSFI_MVLEM3DKlocal(2, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 3) = FourNodeSFI_MVLEM3DKlocal(3, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 4) = FourNodeSFI_MVLEM3DKlocal(4, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 5) = FourNodeSFI_MVLEM3DKlocal(5, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 6) = FourNodeSFI_MVLEM3DKlocal(6, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 7) = FourNodeSFI_MVLEM3DKlocal(7, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 8) = FourNodeSFI_MVLEM3DKlocal(8, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 9) = FourNodeSFI_MVLEM3DKlocal(9, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 10) = FourNodeSFI_MVLEM3DKlocal(10, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 11) = FourNodeSFI_MVLEM3DKlocal(11, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 12) = FourNodeSFI_MVLEM3DKlocal(12, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 13) = Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) + (12.0 * Eim*Iim) / (Lw*Lw*Lw) - (d*(e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(13, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 17) = (e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(13, 18) = (Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(13, 19) = Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - (12.0 * Eim*Iim) / (Lw*Lw*Lw) + (d*(e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(13, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 23) = (e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);

	FourNodeSFI_MVLEM3DKlocal(14, 0) = FourNodeSFI_MVLEM3DKlocal(0, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 1) = FourNodeSFI_MVLEM3DKlocal(1, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 2) = FourNodeSFI_MVLEM3DKlocal(2, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 3) = FourNodeSFI_MVLEM3DKlocal(3, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 4) = FourNodeSFI_MVLEM3DKlocal(4, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 5) = FourNodeSFI_MVLEM3DKlocal(5, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 6) = FourNodeSFI_MVLEM3DKlocal(6, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 7) = FourNodeSFI_MVLEM3DKlocal(7, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 8) = FourNodeSFI_MVLEM3DKlocal(8, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 9) = FourNodeSFI_MVLEM3DKlocal(9, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 10) = FourNodeSFI_MVLEM3DKlocal(10, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 11) = FourNodeSFI_MVLEM3DKlocal(11, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 12) = FourNodeSFI_MVLEM3DKlocal(12, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 13) = FourNodeSFI_MVLEM3DKlocal(13, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 15) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(14, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(14, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(14, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(15, 0) = FourNodeSFI_MVLEM3DKlocal(0, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 1) = FourNodeSFI_MVLEM3DKlocal(1, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 2) = FourNodeSFI_MVLEM3DKlocal(2, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 3) = FourNodeSFI_MVLEM3DKlocal(3, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 4) = FourNodeSFI_MVLEM3DKlocal(4, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 5) = FourNodeSFI_MVLEM3DKlocal(5, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 6) = FourNodeSFI_MVLEM3DKlocal(6, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 7) = FourNodeSFI_MVLEM3DKlocal(7, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 8) = FourNodeSFI_MVLEM3DKlocal(8, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 9) = FourNodeSFI_MVLEM3DKlocal(9, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 10) = FourNodeSFI_MVLEM3DKlocal(10, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 11) = FourNodeSFI_MVLEM3DKlocal(11, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 12) = FourNodeSFI_MVLEM3DKlocal(12, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 13) = FourNodeSFI_MVLEM3DKlocal(13, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 14) = FourNodeSFI_MVLEM3DKlocal(14, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(15, 16) = -(Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeSFI_MVLEM3DKlocal(15, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(15, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(15, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(15, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(15, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(15, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(15, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(16, 0) = FourNodeSFI_MVLEM3DKlocal(0, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 1) = FourNodeSFI_MVLEM3DKlocal(1, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 2) = FourNodeSFI_MVLEM3DKlocal(2, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 3) = FourNodeSFI_MVLEM3DKlocal(3, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 4) = FourNodeSFI_MVLEM3DKlocal(4, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 5) = FourNodeSFI_MVLEM3DKlocal(5, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 6) = FourNodeSFI_MVLEM3DKlocal(6, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 7) = FourNodeSFI_MVLEM3DKlocal(7, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 8) = FourNodeSFI_MVLEM3DKlocal(8, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 9) = FourNodeSFI_MVLEM3DKlocal(9, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 10) = FourNodeSFI_MVLEM3DKlocal(10, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 11) = FourNodeSFI_MVLEM3DKlocal(11, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 12) = FourNodeSFI_MVLEM3DKlocal(12, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 13) = FourNodeSFI_MVLEM3DKlocal(13, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 14) = FourNodeSFI_MVLEM3DKlocal(14, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 15) = FourNodeSFI_MVLEM3DKlocal(15, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(16, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(16, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(16, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(16, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(16, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(16, 22) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(16, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(17, 0) = FourNodeSFI_MVLEM3DKlocal(0, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 1) = FourNodeSFI_MVLEM3DKlocal(1, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 2) = FourNodeSFI_MVLEM3DKlocal(2, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 3) = FourNodeSFI_MVLEM3DKlocal(3, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 4) = FourNodeSFI_MVLEM3DKlocal(4, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 5) = FourNodeSFI_MVLEM3DKlocal(5, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 6) = FourNodeSFI_MVLEM3DKlocal(6, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 7) = FourNodeSFI_MVLEM3DKlocal(7, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 8) = FourNodeSFI_MVLEM3DKlocal(8, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 9) = FourNodeSFI_MVLEM3DKlocal(9, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 10) = FourNodeSFI_MVLEM3DKlocal(10, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 11) = FourNodeSFI_MVLEM3DKlocal(11, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 12) = FourNodeSFI_MVLEM3DKlocal(12, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 13) = FourNodeSFI_MVLEM3DKlocal(13, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 14) = FourNodeSFI_MVLEM3DKlocal(14, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 15) = FourNodeSFI_MVLEM3DKlocal(15, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 16) = FourNodeSFI_MVLEM3DKlocal(16, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 17) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw + 0*K_drilling;
	FourNodeSFI_MVLEM3DKlocal(17, 18) = -(Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(17, 19) = e / (2.0 * (2.0 * (d*d) + 2.0)) - (6.0 * Eim*Iim) / (Lw*Lw) + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(17, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(17, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(17, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(17, 23) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (2.0 * Eim*Iim) / Lw + 0*K_drilling;

	FourNodeSFI_MVLEM3DKlocal(18, 0) = FourNodeSFI_MVLEM3DKlocal(0, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 1) = FourNodeSFI_MVLEM3DKlocal(1, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 2) = FourNodeSFI_MVLEM3DKlocal(2, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 3) = FourNodeSFI_MVLEM3DKlocal(3, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 4) = FourNodeSFI_MVLEM3DKlocal(4, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 5) = FourNodeSFI_MVLEM3DKlocal(5, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 6) = FourNodeSFI_MVLEM3DKlocal(6, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 7) = FourNodeSFI_MVLEM3DKlocal(7, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 8) = FourNodeSFI_MVLEM3DKlocal(8, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 9) = FourNodeSFI_MVLEM3DKlocal(9, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 10) = FourNodeSFI_MVLEM3DKlocal(10, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 11) = FourNodeSFI_MVLEM3DKlocal(11, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 12) = FourNodeSFI_MVLEM3DKlocal(12, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 13) = FourNodeSFI_MVLEM3DKlocal(13, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 14) = FourNodeSFI_MVLEM3DKlocal(14, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 15) = FourNodeSFI_MVLEM3DKlocal(15, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 16) = FourNodeSFI_MVLEM3DKlocal(16, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 17) = FourNodeSFI_MVLEM3DKlocal(17, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 18) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(18, 19) = -(Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(18, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(18, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(18, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(18, 23) = -(Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(19, 0) = FourNodeSFI_MVLEM3DKlocal(0, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 1) = FourNodeSFI_MVLEM3DKlocal(1, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 2) = FourNodeSFI_MVLEM3DKlocal(2, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 3) = FourNodeSFI_MVLEM3DKlocal(3, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 4) = FourNodeSFI_MVLEM3DKlocal(4, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 5) = FourNodeSFI_MVLEM3DKlocal(5, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 6) = FourNodeSFI_MVLEM3DKlocal(6, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 7) = FourNodeSFI_MVLEM3DKlocal(7, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 8) = FourNodeSFI_MVLEM3DKlocal(8, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 9) = FourNodeSFI_MVLEM3DKlocal(9, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 10) = FourNodeSFI_MVLEM3DKlocal(10, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 11) = FourNodeSFI_MVLEM3DKlocal(11, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 12) = FourNodeSFI_MVLEM3DKlocal(12, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 13) = FourNodeSFI_MVLEM3DKlocal(13, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 14) = FourNodeSFI_MVLEM3DKlocal(14, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 15) = FourNodeSFI_MVLEM3DKlocal(15, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 16) = FourNodeSFI_MVLEM3DKlocal(16, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 17) = FourNodeSFI_MVLEM3DKlocal(17, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 18) = FourNodeSFI_MVLEM3DKlocal(18, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 19) = Kv / 4.0 + (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) + (12.0 * Eim*Iim) / (Lw*Lw*Lw) + (d*(e / 2.0 + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(19, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(19, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(19, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(19, 23) = (e / 2.0 + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) - (6.0 * Eim*Iim) / (Lw*Lw);

	FourNodeSFI_MVLEM3DKlocal(20, 0) = FourNodeSFI_MVLEM3DKlocal(0, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 1) = FourNodeSFI_MVLEM3DKlocal(1, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 2) = FourNodeSFI_MVLEM3DKlocal(2, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 3) = FourNodeSFI_MVLEM3DKlocal(3, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 4) = FourNodeSFI_MVLEM3DKlocal(4, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 5) = FourNodeSFI_MVLEM3DKlocal(5, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 6) = FourNodeSFI_MVLEM3DKlocal(6, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 7) = FourNodeSFI_MVLEM3DKlocal(7, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 8) = FourNodeSFI_MVLEM3DKlocal(8, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 9) = FourNodeSFI_MVLEM3DKlocal(9, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 10) = FourNodeSFI_MVLEM3DKlocal(10, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 11) = FourNodeSFI_MVLEM3DKlocal(11, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 12) = FourNodeSFI_MVLEM3DKlocal(12, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 13) = FourNodeSFI_MVLEM3DKlocal(13, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 14) = FourNodeSFI_MVLEM3DKlocal(14, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 15) = FourNodeSFI_MVLEM3DKlocal(15, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 16) = FourNodeSFI_MVLEM3DKlocal(16, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 17) = FourNodeSFI_MVLEM3DKlocal(17, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 18) = FourNodeSFI_MVLEM3DKlocal(18, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 19) = FourNodeSFI_MVLEM3DKlocal(19, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(20, 21) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(20, 22) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(20, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(21, 0) = FourNodeSFI_MVLEM3DKlocal(0, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 1) = FourNodeSFI_MVLEM3DKlocal(1, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 2) = FourNodeSFI_MVLEM3DKlocal(2, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 3) = FourNodeSFI_MVLEM3DKlocal(3, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 4) = FourNodeSFI_MVLEM3DKlocal(4, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 5) = FourNodeSFI_MVLEM3DKlocal(5, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 6) = FourNodeSFI_MVLEM3DKlocal(6, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 7) = FourNodeSFI_MVLEM3DKlocal(7, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 8) = FourNodeSFI_MVLEM3DKlocal(8, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 9) = FourNodeSFI_MVLEM3DKlocal(9, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 10) = FourNodeSFI_MVLEM3DKlocal(10, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 11) = FourNodeSFI_MVLEM3DKlocal(11, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 12) = FourNodeSFI_MVLEM3DKlocal(12, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 13) = FourNodeSFI_MVLEM3DKlocal(13, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 14) = FourNodeSFI_MVLEM3DKlocal(14, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 15) = FourNodeSFI_MVLEM3DKlocal(15, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 16) = FourNodeSFI_MVLEM3DKlocal(16, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 17) = FourNodeSFI_MVLEM3DKlocal(17, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 18) = FourNodeSFI_MVLEM3DKlocal(18, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 19) = FourNodeSFI_MVLEM3DKlocal(19, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 20) = FourNodeSFI_MVLEM3DKlocal(20, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(21, 22) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeSFI_MVLEM3DKlocal(21, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(22, 0) = FourNodeSFI_MVLEM3DKlocal(0, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 1) = FourNodeSFI_MVLEM3DKlocal(1, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 2) = FourNodeSFI_MVLEM3DKlocal(2, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 3) = FourNodeSFI_MVLEM3DKlocal(3, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 4) = FourNodeSFI_MVLEM3DKlocal(4, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 5) = FourNodeSFI_MVLEM3DKlocal(5, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 6) = FourNodeSFI_MVLEM3DKlocal(6, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 7) = FourNodeSFI_MVLEM3DKlocal(7, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 8) = FourNodeSFI_MVLEM3DKlocal(8, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 9) = FourNodeSFI_MVLEM3DKlocal(9, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 10) = FourNodeSFI_MVLEM3DKlocal(10, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 11) = FourNodeSFI_MVLEM3DKlocal(11, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 12) = FourNodeSFI_MVLEM3DKlocal(12, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 13) = FourNodeSFI_MVLEM3DKlocal(13, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 14) = FourNodeSFI_MVLEM3DKlocal(14, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 15) = FourNodeSFI_MVLEM3DKlocal(15, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 16) = FourNodeSFI_MVLEM3DKlocal(16, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 17) = FourNodeSFI_MVLEM3DKlocal(17, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 18) = FourNodeSFI_MVLEM3DKlocal(18, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 19) = FourNodeSFI_MVLEM3DKlocal(19, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 20) = FourNodeSFI_MVLEM3DKlocal(20, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 21) = FourNodeSFI_MVLEM3DKlocal(21, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(22, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(23, 0) = FourNodeSFI_MVLEM3DKlocal(0, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 1) = FourNodeSFI_MVLEM3DKlocal(1, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 2) = FourNodeSFI_MVLEM3DKlocal(2, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 3) = FourNodeSFI_MVLEM3DKlocal(3, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 4) = FourNodeSFI_MVLEM3DKlocal(4, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 5) = FourNodeSFI_MVLEM3DKlocal(5, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 6) = FourNodeSFI_MVLEM3DKlocal(6, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 7) = FourNodeSFI_MVLEM3DKlocal(7, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 8) = FourNodeSFI_MVLEM3DKlocal(8, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 9) = FourNodeSFI_MVLEM3DKlocal(9, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 10) = FourNodeSFI_MVLEM3DKlocal(10, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 11) = FourNodeSFI_MVLEM3DKlocal(11, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 12) = FourNodeSFI_MVLEM3DKlocal(12, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 13) = FourNodeSFI_MVLEM3DKlocal(13, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 14) = FourNodeSFI_MVLEM3DKlocal(14, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 15) = FourNodeSFI_MVLEM3DKlocal(15, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 16) = FourNodeSFI_MVLEM3DKlocal(16, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 17) = FourNodeSFI_MVLEM3DKlocal(17, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 18) = FourNodeSFI_MVLEM3DKlocal(18, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 19) = FourNodeSFI_MVLEM3DKlocal(19, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 20) = FourNodeSFI_MVLEM3DKlocal(20, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 21) = FourNodeSFI_MVLEM3DKlocal(21, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 22) = FourNodeSFI_MVLEM3DKlocal(22, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 23) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw + 0 * K_drilling;

	FourNodeSFI_MVLEM3DK.addMatrixTripleProduct(0.0, T, FourNodeSFI_MVLEM3DKlocal, 1.0);  // Convert matrix from local to global cs

	/*opserr << "FourNodeSFI_MVLEM3DK :";
	opserr << "\n";
	for (int j = 0; j < 24+m; j++) {
		for (int i = 0; i < 24+m; i++) {
			opserr << FourNodeSFI_MVLEM3DK(i, j) << " , ";
		}
		opserr << "\n";
	}
	opserr << "\n" << "\n";*/
	//opserr << "  14		getInitialStiff\n";
	// Return element stiffness matrix
	return FourNodeSFI_MVLEM3DK;

}

// Get current element tangent stiffness matrix from the material for the last updated strain
const Matrix & FourNodeSFI_MVLEM3D::getTangentStiff(void)
{

	//FourNodeSFI_MVLEM3DK.Zero();		// Global stiffness matrix
	//FourNodeSFI_MVLEM3DKlocal.Zero();	// Local stiffness matrix

	Kh = 0.0;

	for (int i = 0; i < m; i++)
	{
		// Get the material tangent
		const Matrix &D = theMaterial[i]->getTangent();

		double D00 = D(0, 0); double D01 = D(0, 1); double D02 = D(0, 2);
		double D10 = D(1, 0); double D11 = D(1, 1); double D12 = D(1, 2);
		double D20 = D(2, 0); double D21 = D(2, 1); double D22 = D(2, 2);

		kx[i] = D00 * h*t[i] / b[i];
		ky[i] = D11 * b[i] * t[i] / h;
		Kh += D22 * b[i] * t[i] / h;

	}

	//opserr << "K_drilling in TangentStiff: " << K_drilling << "\n";
	// Build the tangent stiffness matrix
	double Kv = 0.0; double Km = 0.0; double e = 0.0; // double ex = 0.0;

	for (int i = 0; i<m; ++i)
	{
		Kv += ky[i];
		Km += ky[i] * x[i] * x[i];
		e += ky[i] * x[i];

		FourNodeSFI_MVLEM3DKlocal(24 + i, 24 + i) = kx[i]; // Diagonal terms accounting for horizontal stiffness
	}

	/*opserr << "Tangent Stiff: ";
	opserr << "kh: " << Kh << ", " << "km: " << Km << ", " << "kv: " << Kv << ", " << "e: " << e << ";";
	opserr << "\n" << "\n";*/

	// Assemble element stiffness matrix
	FourNodeSFI_MVLEM3DKlocal(0, 0) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(0, 1) = (Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 2) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 3) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 4) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 5) = -(Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 6) = Kh / 4.0 - (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(0, 7) = -(Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 11) = -(Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 12) = -Kh / 4.0;
	FourNodeSFI_MVLEM3DKlocal(0, 13) = -(Kh*d*h*(c - 1)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 17) = (Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 18) = -Kh / 4;
	FourNodeSFI_MVLEM3DKlocal(0, 19) = (Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(0, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(0, 23) = (Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(1, 0) = FourNodeSFI_MVLEM3DKlocal(0, 1);
	FourNodeSFI_MVLEM3DKlocal(1, 1) = Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - (d*(e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(1, 2) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 3) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 4) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 5) = (e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(1, 6) = (Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(1, 7) = Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) + (d*(e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) - (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(1, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 11) = (e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(1, 12) = -(Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(1, 13) = (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - Kv / 4.0 + (d*(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(1, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 17) = -(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(1, 18) = -(Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(1, 19) = (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - Kv / 4.0 - (d*(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(1, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(1, 23) = -(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);

	FourNodeSFI_MVLEM3DKlocal(2, 0) = FourNodeSFI_MVLEM3DKlocal(0, 2);
	FourNodeSFI_MVLEM3DKlocal(2, 1) = FourNodeSFI_MVLEM3DKlocal(1, 2);
	FourNodeSFI_MVLEM3DKlocal(2, 2) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 3) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 4) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 5) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 6) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 7) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 9) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 10) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 16) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(2, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(2, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(3, 0) = FourNodeSFI_MVLEM3DKlocal(0, 3);
	FourNodeSFI_MVLEM3DKlocal(3, 1) = FourNodeSFI_MVLEM3DKlocal(1, 3);
	FourNodeSFI_MVLEM3DKlocal(3, 2) = FourNodeSFI_MVLEM3DKlocal(2, 3);
	FourNodeSFI_MVLEM3DKlocal(3, 3) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 4) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeSFI_MVLEM3DKlocal(3, 5) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 6) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 7) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(3, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(3, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(4, 0) = FourNodeSFI_MVLEM3DKlocal(0, 4);
	FourNodeSFI_MVLEM3DKlocal(4, 1) = FourNodeSFI_MVLEM3DKlocal(1, 4);
	FourNodeSFI_MVLEM3DKlocal(4, 2) = FourNodeSFI_MVLEM3DKlocal(2, 4);
	FourNodeSFI_MVLEM3DKlocal(4, 3) = FourNodeSFI_MVLEM3DKlocal(3, 4);
	FourNodeSFI_MVLEM3DKlocal(4, 4) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 5) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 6) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 7) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 8) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 10) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 16) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(4, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(4, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(5, 0) = FourNodeSFI_MVLEM3DKlocal(0, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 1) = FourNodeSFI_MVLEM3DKlocal(1, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 2) = FourNodeSFI_MVLEM3DKlocal(2, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 3) = FourNodeSFI_MVLEM3DKlocal(3, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 4) = FourNodeSFI_MVLEM3DKlocal(4, 5);
	FourNodeSFI_MVLEM3DKlocal(5, 5) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw + 0*K_drilling;
	FourNodeSFI_MVLEM3DKlocal(5, 6) = -(Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 7) = e / (2.0 * (2.0 * (d*d) + 2.0)) + (d*(Km + Kh*(c*c)*(h*h))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(5, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 11) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (2.0 * Eim*Iim) / Lw + 0*K_drilling;
	FourNodeSFI_MVLEM3DKlocal(5, 12) = (Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 13) = (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - e / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 17) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 18) = (Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 19) = -e / (2.0 * (2.0 * (d*d) + 2.0)) - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(5, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(5, 23) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(6, 0) = FourNodeSFI_MVLEM3DKlocal(0, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 1) = FourNodeSFI_MVLEM3DKlocal(1, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 2) = FourNodeSFI_MVLEM3DKlocal(2, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 3) = FourNodeSFI_MVLEM3DKlocal(3, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 4) = FourNodeSFI_MVLEM3DKlocal(4, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 5) = FourNodeSFI_MVLEM3DKlocal(5, 6);
	FourNodeSFI_MVLEM3DKlocal(6, 6) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(6, 7) = -(Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 11) = -(Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 12) = -Kh / 4;
	FourNodeSFI_MVLEM3DKlocal(6, 13) = -(Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 17) = (Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 18) = -Kh / 4.0;
	FourNodeSFI_MVLEM3DKlocal(6, 19) = (Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(6, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(6, 23) = (Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(7, 0) = FourNodeSFI_MVLEM3DKlocal(0, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 1) = FourNodeSFI_MVLEM3DKlocal(1, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 2) = FourNodeSFI_MVLEM3DKlocal(2, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 3) = FourNodeSFI_MVLEM3DKlocal(3, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 4) = FourNodeSFI_MVLEM3DKlocal(4, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 5) = FourNodeSFI_MVLEM3DKlocal(5, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 6) = FourNodeSFI_MVLEM3DKlocal(6, 7);
	FourNodeSFI_MVLEM3DKlocal(7, 7) = Kv / 4.0 + (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) + (d*(e / 2.0 + (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(7, 8) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 9) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 10) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 11) = (e / 2.0 + (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) - (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(7, 12) = (Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(7, 13) = (d*(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - Kv / 4.0;
	FourNodeSFI_MVLEM3DKlocal(7, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 17) = -(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(7, 18) = (Kh*c*d*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(7, 19) = -Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - (d*(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(7, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(7, 23) = -(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);

	FourNodeSFI_MVLEM3DKlocal(8, 0) = FourNodeSFI_MVLEM3DKlocal(0, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 1) = FourNodeSFI_MVLEM3DKlocal(1, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 2) = FourNodeSFI_MVLEM3DKlocal(2, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 3) = FourNodeSFI_MVLEM3DKlocal(3, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 4) = FourNodeSFI_MVLEM3DKlocal(4, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 5) = FourNodeSFI_MVLEM3DKlocal(5, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 6) = FourNodeSFI_MVLEM3DKlocal(6, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 7) = FourNodeSFI_MVLEM3DKlocal(7, 8);
	FourNodeSFI_MVLEM3DKlocal(8, 8) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 10) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 16) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(8, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(8, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(9, 0) = FourNodeSFI_MVLEM3DKlocal(0, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 1) = FourNodeSFI_MVLEM3DKlocal(1, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 2) = FourNodeSFI_MVLEM3DKlocal(2, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 3) = FourNodeSFI_MVLEM3DKlocal(3, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 4) = FourNodeSFI_MVLEM3DKlocal(4, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 5) = FourNodeSFI_MVLEM3DKlocal(5, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 6) = FourNodeSFI_MVLEM3DKlocal(6, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 7) = FourNodeSFI_MVLEM3DKlocal(7, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 8) = FourNodeSFI_MVLEM3DKlocal(8, 9);
	FourNodeSFI_MVLEM3DKlocal(9, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 10) = -(Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeSFI_MVLEM3DKlocal(9, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(9, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(9, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(10, 0) = FourNodeSFI_MVLEM3DKlocal(0, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 1) = FourNodeSFI_MVLEM3DKlocal(1, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 2) = FourNodeSFI_MVLEM3DKlocal(2, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 3) = FourNodeSFI_MVLEM3DKlocal(3, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 4) = FourNodeSFI_MVLEM3DKlocal(4, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 5) = FourNodeSFI_MVLEM3DKlocal(5, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 6) = FourNodeSFI_MVLEM3DKlocal(6, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 7) = FourNodeSFI_MVLEM3DKlocal(7, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 8) = FourNodeSFI_MVLEM3DKlocal(8, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 9) = FourNodeSFI_MVLEM3DKlocal(9, 10);
	FourNodeSFI_MVLEM3DKlocal(10, 10) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 11) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 12) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 13) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(10, 22) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(10, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(11, 0) = FourNodeSFI_MVLEM3DKlocal(0, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 1) = FourNodeSFI_MVLEM3DKlocal(1, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 2) = FourNodeSFI_MVLEM3DKlocal(2, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 3) = FourNodeSFI_MVLEM3DKlocal(3, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 4) = FourNodeSFI_MVLEM3DKlocal(4, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 5) = FourNodeSFI_MVLEM3DKlocal(5, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 6) = FourNodeSFI_MVLEM3DKlocal(6, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 7) = FourNodeSFI_MVLEM3DKlocal(7, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 8) = FourNodeSFI_MVLEM3DKlocal(8, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 9) = FourNodeSFI_MVLEM3DKlocal(9, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 10) = FourNodeSFI_MVLEM3DKlocal(10, 11);
	FourNodeSFI_MVLEM3DKlocal(11, 11) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw + 0*K_drilling;
	FourNodeSFI_MVLEM3DKlocal(11, 12) = (Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 13) = (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - e / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 17) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 18) = (Kh*c*h) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 19) = -e / (2.0 * (2.0 * (d*d) + 2.0)) - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(11, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(11, 23) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(12, 0) = FourNodeSFI_MVLEM3DKlocal(0, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 1) = FourNodeSFI_MVLEM3DKlocal(1, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 2) = FourNodeSFI_MVLEM3DKlocal(2, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 3) = FourNodeSFI_MVLEM3DKlocal(3, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 4) = FourNodeSFI_MVLEM3DKlocal(4, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 5) = FourNodeSFI_MVLEM3DKlocal(5, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 6) = FourNodeSFI_MVLEM3DKlocal(6, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 7) = FourNodeSFI_MVLEM3DKlocal(7, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 8) = FourNodeSFI_MVLEM3DKlocal(8, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 9) = FourNodeSFI_MVLEM3DKlocal(9, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 10) = FourNodeSFI_MVLEM3DKlocal(10, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 11) = FourNodeSFI_MVLEM3DKlocal(11, 12);
	FourNodeSFI_MVLEM3DKlocal(12, 12) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(12, 13) = (Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(12, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 17) = -(Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(12, 18) = Kh / 4.0 - (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(12, 19) = -(Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(12, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(12, 23) = -(Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(13, 0) = FourNodeSFI_MVLEM3DKlocal(0, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 1) = FourNodeSFI_MVLEM3DKlocal(1, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 2) = FourNodeSFI_MVLEM3DKlocal(2, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 3) = FourNodeSFI_MVLEM3DKlocal(3, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 4) = FourNodeSFI_MVLEM3DKlocal(4, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 5) = FourNodeSFI_MVLEM3DKlocal(5, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 6) = FourNodeSFI_MVLEM3DKlocal(6, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 7) = FourNodeSFI_MVLEM3DKlocal(7, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 8) = FourNodeSFI_MVLEM3DKlocal(8, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 9) = FourNodeSFI_MVLEM3DKlocal(9, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 10) = FourNodeSFI_MVLEM3DKlocal(10, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 11) = FourNodeSFI_MVLEM3DKlocal(11, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 12) = FourNodeSFI_MVLEM3DKlocal(12, 13);
	FourNodeSFI_MVLEM3DKlocal(13, 13) = Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) + (12.0 * Eim*Iim) / (Lw*Lw*Lw) - (d*(e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(13, 14) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 15) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 16) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 17) = (e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeSFI_MVLEM3DKlocal(13, 18) = (Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(13, 19) = Kv / 4.0 - (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) - (12.0 * Eim*Iim) / (Lw*Lw*Lw) + (d*(e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(13, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(13, 23) = (e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);

	FourNodeSFI_MVLEM3DKlocal(14, 0) = FourNodeSFI_MVLEM3DKlocal(0, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 1) = FourNodeSFI_MVLEM3DKlocal(1, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 2) = FourNodeSFI_MVLEM3DKlocal(2, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 3) = FourNodeSFI_MVLEM3DKlocal(3, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 4) = FourNodeSFI_MVLEM3DKlocal(4, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 5) = FourNodeSFI_MVLEM3DKlocal(5, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 6) = FourNodeSFI_MVLEM3DKlocal(6, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 7) = FourNodeSFI_MVLEM3DKlocal(7, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 8) = FourNodeSFI_MVLEM3DKlocal(8, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 9) = FourNodeSFI_MVLEM3DKlocal(9, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 10) = FourNodeSFI_MVLEM3DKlocal(10, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 11) = FourNodeSFI_MVLEM3DKlocal(11, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 12) = FourNodeSFI_MVLEM3DKlocal(12, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 13) = FourNodeSFI_MVLEM3DKlocal(13, 14);
	FourNodeSFI_MVLEM3DKlocal(14, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 15) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(14, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(14, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(14, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(14, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(15, 0) = FourNodeSFI_MVLEM3DKlocal(0, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 1) = FourNodeSFI_MVLEM3DKlocal(1, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 2) = FourNodeSFI_MVLEM3DKlocal(2, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 3) = FourNodeSFI_MVLEM3DKlocal(3, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 4) = FourNodeSFI_MVLEM3DKlocal(4, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 5) = FourNodeSFI_MVLEM3DKlocal(5, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 6) = FourNodeSFI_MVLEM3DKlocal(6, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 7) = FourNodeSFI_MVLEM3DKlocal(7, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 8) = FourNodeSFI_MVLEM3DKlocal(8, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 9) = FourNodeSFI_MVLEM3DKlocal(9, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 10) = FourNodeSFI_MVLEM3DKlocal(10, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 11) = FourNodeSFI_MVLEM3DKlocal(11, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 12) = FourNodeSFI_MVLEM3DKlocal(12, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 13) = FourNodeSFI_MVLEM3DKlocal(13, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 14) = FourNodeSFI_MVLEM3DKlocal(14, 15);
	FourNodeSFI_MVLEM3DKlocal(15, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(15, 16) = -(Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeSFI_MVLEM3DKlocal(15, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(15, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(15, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(15, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(15, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(15, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(15, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(16, 0) = FourNodeSFI_MVLEM3DKlocal(0, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 1) = FourNodeSFI_MVLEM3DKlocal(1, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 2) = FourNodeSFI_MVLEM3DKlocal(2, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 3) = FourNodeSFI_MVLEM3DKlocal(3, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 4) = FourNodeSFI_MVLEM3DKlocal(4, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 5) = FourNodeSFI_MVLEM3DKlocal(5, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 6) = FourNodeSFI_MVLEM3DKlocal(6, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 7) = FourNodeSFI_MVLEM3DKlocal(7, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 8) = FourNodeSFI_MVLEM3DKlocal(8, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 9) = FourNodeSFI_MVLEM3DKlocal(9, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 10) = FourNodeSFI_MVLEM3DKlocal(10, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 11) = FourNodeSFI_MVLEM3DKlocal(11, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 12) = FourNodeSFI_MVLEM3DKlocal(12, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 13) = FourNodeSFI_MVLEM3DKlocal(13, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 14) = FourNodeSFI_MVLEM3DKlocal(14, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 15) = FourNodeSFI_MVLEM3DKlocal(15, 16);
	FourNodeSFI_MVLEM3DKlocal(16, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(16, 17) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(16, 18) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(16, 19) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(16, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(16, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(16, 22) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(16, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(17, 0) = FourNodeSFI_MVLEM3DKlocal(0, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 1) = FourNodeSFI_MVLEM3DKlocal(1, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 2) = FourNodeSFI_MVLEM3DKlocal(2, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 3) = FourNodeSFI_MVLEM3DKlocal(3, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 4) = FourNodeSFI_MVLEM3DKlocal(4, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 5) = FourNodeSFI_MVLEM3DKlocal(5, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 6) = FourNodeSFI_MVLEM3DKlocal(6, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 7) = FourNodeSFI_MVLEM3DKlocal(7, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 8) = FourNodeSFI_MVLEM3DKlocal(8, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 9) = FourNodeSFI_MVLEM3DKlocal(9, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 10) = FourNodeSFI_MVLEM3DKlocal(10, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 11) = FourNodeSFI_MVLEM3DKlocal(11, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 12) = FourNodeSFI_MVLEM3DKlocal(12, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 13) = FourNodeSFI_MVLEM3DKlocal(13, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 14) = FourNodeSFI_MVLEM3DKlocal(14, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 15) = FourNodeSFI_MVLEM3DKlocal(15, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 16) = FourNodeSFI_MVLEM3DKlocal(16, 17);
	FourNodeSFI_MVLEM3DKlocal(17, 17) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw + 0*K_drilling;
	FourNodeSFI_MVLEM3DKlocal(17, 18) = -(Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(17, 19) = e / (2.0 * (2.0 * (d*d) + 2.0)) - (6.0 * Eim*Iim) / (Lw*Lw) + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(17, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(17, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(17, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(17, 23) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (2.0 * Eim*Iim) / Lw + 0*K_drilling;

	FourNodeSFI_MVLEM3DKlocal(18, 0) = FourNodeSFI_MVLEM3DKlocal(0, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 1) = FourNodeSFI_MVLEM3DKlocal(1, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 2) = FourNodeSFI_MVLEM3DKlocal(2, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 3) = FourNodeSFI_MVLEM3DKlocal(3, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 4) = FourNodeSFI_MVLEM3DKlocal(4, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 5) = FourNodeSFI_MVLEM3DKlocal(5, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 6) = FourNodeSFI_MVLEM3DKlocal(6, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 7) = FourNodeSFI_MVLEM3DKlocal(7, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 8) = FourNodeSFI_MVLEM3DKlocal(8, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 9) = FourNodeSFI_MVLEM3DKlocal(9, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 10) = FourNodeSFI_MVLEM3DKlocal(10, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 11) = FourNodeSFI_MVLEM3DKlocal(11, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 12) = FourNodeSFI_MVLEM3DKlocal(12, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 13) = FourNodeSFI_MVLEM3DKlocal(13, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 14) = FourNodeSFI_MVLEM3DKlocal(14, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 15) = FourNodeSFI_MVLEM3DKlocal(15, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 16) = FourNodeSFI_MVLEM3DKlocal(16, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 17) = FourNodeSFI_MVLEM3DKlocal(17, 18);
	FourNodeSFI_MVLEM3DKlocal(18, 18) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeSFI_MVLEM3DKlocal(18, 19) = -(Kh*d*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));
	FourNodeSFI_MVLEM3DKlocal(18, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(18, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(18, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(18, 23) = -(Kh*h*(c - 1.0)) / (2.0 * (2.0 * (d*d) + 2.0));

	FourNodeSFI_MVLEM3DKlocal(19, 0) = FourNodeSFI_MVLEM3DKlocal(0, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 1) = FourNodeSFI_MVLEM3DKlocal(1, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 2) = FourNodeSFI_MVLEM3DKlocal(2, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 3) = FourNodeSFI_MVLEM3DKlocal(3, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 4) = FourNodeSFI_MVLEM3DKlocal(4, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 5) = FourNodeSFI_MVLEM3DKlocal(5, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 6) = FourNodeSFI_MVLEM3DKlocal(6, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 7) = FourNodeSFI_MVLEM3DKlocal(7, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 8) = FourNodeSFI_MVLEM3DKlocal(8, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 9) = FourNodeSFI_MVLEM3DKlocal(9, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 10) = FourNodeSFI_MVLEM3DKlocal(10, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 11) = FourNodeSFI_MVLEM3DKlocal(11, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 12) = FourNodeSFI_MVLEM3DKlocal(12, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 13) = FourNodeSFI_MVLEM3DKlocal(13, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 14) = FourNodeSFI_MVLEM3DKlocal(14, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 15) = FourNodeSFI_MVLEM3DKlocal(15, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 16) = FourNodeSFI_MVLEM3DKlocal(16, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 17) = FourNodeSFI_MVLEM3DKlocal(17, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 18) = FourNodeSFI_MVLEM3DKlocal(18, 19);
	FourNodeSFI_MVLEM3DKlocal(19, 19) = Kv / 4.0 + (d*e) / (2.0 * (2.0 * (d*d) + 2.0)) + (12.0 * Eim*Iim) / (Lw*Lw*Lw) + (d*(e / 2.0 + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeSFI_MVLEM3DKlocal(19, 20) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(19, 21) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(19, 22) = 0.0;
	FourNodeSFI_MVLEM3DKlocal(19, 23) = (e / 2.0 + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) - (6.0 * Eim*Iim) / (Lw*Lw);

	FourNodeSFI_MVLEM3DKlocal(20, 0) = FourNodeSFI_MVLEM3DKlocal(0, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 1) = FourNodeSFI_MVLEM3DKlocal(1, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 2) = FourNodeSFI_MVLEM3DKlocal(2, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 3) = FourNodeSFI_MVLEM3DKlocal(3, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 4) = FourNodeSFI_MVLEM3DKlocal(4, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 5) = FourNodeSFI_MVLEM3DKlocal(5, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 6) = FourNodeSFI_MVLEM3DKlocal(6, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 7) = FourNodeSFI_MVLEM3DKlocal(7, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 8) = FourNodeSFI_MVLEM3DKlocal(8, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 9) = FourNodeSFI_MVLEM3DKlocal(9, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 10) = FourNodeSFI_MVLEM3DKlocal(10, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 11) = FourNodeSFI_MVLEM3DKlocal(11, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 12) = FourNodeSFI_MVLEM3DKlocal(12, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 13) = FourNodeSFI_MVLEM3DKlocal(13, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 14) = FourNodeSFI_MVLEM3DKlocal(14, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 15) = FourNodeSFI_MVLEM3DKlocal(15, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 16) = FourNodeSFI_MVLEM3DKlocal(16, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 17) = FourNodeSFI_MVLEM3DKlocal(17, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 18) = FourNodeSFI_MVLEM3DKlocal(18, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 19) = FourNodeSFI_MVLEM3DKlocal(19, 20);
	FourNodeSFI_MVLEM3DKlocal(20, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(20, 21) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(20, 22) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(20, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(21, 0) = FourNodeSFI_MVLEM3DKlocal(0, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 1) = FourNodeSFI_MVLEM3DKlocal(1, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 2) = FourNodeSFI_MVLEM3DKlocal(2, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 3) = FourNodeSFI_MVLEM3DKlocal(3, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 4) = FourNodeSFI_MVLEM3DKlocal(4, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 5) = FourNodeSFI_MVLEM3DKlocal(5, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 6) = FourNodeSFI_MVLEM3DKlocal(6, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 7) = FourNodeSFI_MVLEM3DKlocal(7, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 8) = FourNodeSFI_MVLEM3DKlocal(8, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 9) = FourNodeSFI_MVLEM3DKlocal(9, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 10) = FourNodeSFI_MVLEM3DKlocal(10, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 11) = FourNodeSFI_MVLEM3DKlocal(11, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 12) = FourNodeSFI_MVLEM3DKlocal(12, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 13) = FourNodeSFI_MVLEM3DKlocal(13, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 14) = FourNodeSFI_MVLEM3DKlocal(14, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 15) = FourNodeSFI_MVLEM3DKlocal(15, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 16) = FourNodeSFI_MVLEM3DKlocal(16, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 17) = FourNodeSFI_MVLEM3DKlocal(17, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 18) = FourNodeSFI_MVLEM3DKlocal(18, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 19) = FourNodeSFI_MVLEM3DKlocal(19, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 20) = FourNodeSFI_MVLEM3DKlocal(20, 21);
	FourNodeSFI_MVLEM3DKlocal(21, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(21, 22) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeSFI_MVLEM3DKlocal(21, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(22, 0) = FourNodeSFI_MVLEM3DKlocal(0, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 1) = FourNodeSFI_MVLEM3DKlocal(1, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 2) = FourNodeSFI_MVLEM3DKlocal(2, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 3) = FourNodeSFI_MVLEM3DKlocal(3, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 4) = FourNodeSFI_MVLEM3DKlocal(4, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 5) = FourNodeSFI_MVLEM3DKlocal(5, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 6) = FourNodeSFI_MVLEM3DKlocal(6, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 7) = FourNodeSFI_MVLEM3DKlocal(7, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 8) = FourNodeSFI_MVLEM3DKlocal(8, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 9) = FourNodeSFI_MVLEM3DKlocal(9, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 10) = FourNodeSFI_MVLEM3DKlocal(10, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 11) = FourNodeSFI_MVLEM3DKlocal(11, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 12) = FourNodeSFI_MVLEM3DKlocal(12, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 13) = FourNodeSFI_MVLEM3DKlocal(13, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 14) = FourNodeSFI_MVLEM3DKlocal(14, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 15) = FourNodeSFI_MVLEM3DKlocal(15, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 16) = FourNodeSFI_MVLEM3DKlocal(16, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 17) = FourNodeSFI_MVLEM3DKlocal(17, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 18) = FourNodeSFI_MVLEM3DKlocal(18, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 19) = FourNodeSFI_MVLEM3DKlocal(19, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 20) = FourNodeSFI_MVLEM3DKlocal(20, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 21) = FourNodeSFI_MVLEM3DKlocal(21, 22);
	FourNodeSFI_MVLEM3DKlocal(22, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DKlocal(22, 23) = 0.0;

	FourNodeSFI_MVLEM3DKlocal(23, 0) = FourNodeSFI_MVLEM3DKlocal(0, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 1) = FourNodeSFI_MVLEM3DKlocal(1, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 2) = FourNodeSFI_MVLEM3DKlocal(2, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 3) = FourNodeSFI_MVLEM3DKlocal(3, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 4) = FourNodeSFI_MVLEM3DKlocal(4, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 5) = FourNodeSFI_MVLEM3DKlocal(5, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 6) = FourNodeSFI_MVLEM3DKlocal(6, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 7) = FourNodeSFI_MVLEM3DKlocal(7, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 8) = FourNodeSFI_MVLEM3DKlocal(8, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 9) = FourNodeSFI_MVLEM3DKlocal(9, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 10) = FourNodeSFI_MVLEM3DKlocal(10, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 11) = FourNodeSFI_MVLEM3DKlocal(11, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 12) = FourNodeSFI_MVLEM3DKlocal(12, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 13) = FourNodeSFI_MVLEM3DKlocal(13, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 14) = FourNodeSFI_MVLEM3DKlocal(14, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 15) = FourNodeSFI_MVLEM3DKlocal(15, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 16) = FourNodeSFI_MVLEM3DKlocal(16, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 17) = FourNodeSFI_MVLEM3DKlocal(17, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 18) = FourNodeSFI_MVLEM3DKlocal(18, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 19) = FourNodeSFI_MVLEM3DKlocal(19, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 20) = FourNodeSFI_MVLEM3DKlocal(20, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 21) = FourNodeSFI_MVLEM3DKlocal(21, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 22) = FourNodeSFI_MVLEM3DKlocal(22, 23);
	FourNodeSFI_MVLEM3DKlocal(23, 23) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw + 0 * K_drilling;

	FourNodeSFI_MVLEM3DK.addMatrixTripleProduct(0.0, T, FourNodeSFI_MVLEM3DKlocal, 1.0); // Convert matrix from local to global cs
	//opserr << "  15		getTangentStiff\n";
	// Return element Global stiffness matrix
	return FourNodeSFI_MVLEM3DK;

}

// Get element mass matrix assuming lumped mass
const Matrix & FourNodeSFI_MVLEM3D::getMass(void)
{

	//FourNodeSFI_MVLEM3DM.Zero();
	//FourNodeSFI_MVLEM3DMlocal.Zero();

	// No rotational mass, no mass at internal (dummy) nodes
	FourNodeSFI_MVLEM3DMlocal(0, 0) = NodeMass;
	FourNodeSFI_MVLEM3DMlocal(1, 1) = NodeMass;
	FourNodeSFI_MVLEM3DMlocal(2, 2) = NodeMass;

	FourNodeSFI_MVLEM3DMlocal(6, 6) = NodeMass;
	FourNodeSFI_MVLEM3DMlocal(7, 7) = NodeMass;
	FourNodeSFI_MVLEM3DMlocal(8, 8) = NodeMass;

	FourNodeSFI_MVLEM3DMlocal(12, 12) = NodeMass;
	FourNodeSFI_MVLEM3DMlocal(13, 13) = NodeMass;
	FourNodeSFI_MVLEM3DMlocal(14, 14) = NodeMass;

	FourNodeSFI_MVLEM3DMlocal(18, 18) = NodeMass;
	FourNodeSFI_MVLEM3DMlocal(19, 19) = NodeMass;
	FourNodeSFI_MVLEM3DMlocal(20, 20) = NodeMass;

	FourNodeSFI_MVLEM3DM.addMatrixTripleProduct(0.0, T, FourNodeSFI_MVLEM3DMlocal, 1.0); // Convert matrix from local to global cs
	//opserr << "  16		getMass\n";
	// Return element mass matrix
	return FourNodeSFI_MVLEM3DM;
}

// Get element damping matrix
const Matrix & FourNodeSFI_MVLEM3D::getDamp(void)
{
	//FourNodeSFI_MVLEM3DD.Zero();

	FourNodeSFI_MVLEM3DD = this->Element::getDamp();
	//opserr << "  17		getDamp\n";
	// Return element damping matrix
	return FourNodeSFI_MVLEM3DD; // !!! Is this already in global cs? I dont see you are doing a transformation. CHECK WITH OTHER ELEMENTS
}

// N/A to this model - no element loads
void FourNodeSFI_MVLEM3D::zeroLoad(void)
{
	//opserr << "  18		zeroLoad\n";
	// does nothing 
}

// N/A to this model - no element loads
int FourNodeSFI_MVLEM3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	//opserr << "  19		addLoad\n";
	return 0;
}

// N/A to this model - no element loads
int FourNodeSFI_MVLEM3D::addInertiaLoadToUnbalance(const Vector &accel)
{
	//opserr << "  20		addInertiaLoadToUnbalance\n";
	return 0;
}

// Get element force vector
const Vector & FourNodeSFI_MVLEM3D::getResistingForce()
{

	//FourNodeSFI_MVLEM3DR.Zero();
	//FourNodeSFI_MVLEM3DRlocal.Zero();

	// Get Trial Displacements
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();

	Vector dispG(24 + m); // Vector of total 24+m displacemets in global coordinates
	dispG.Zero();
	Vector dispL(24 + m); // Vector of total 24+m displacemets in local coordinates
	//dispL.Zero();

	// Assigning all displacements in global CS into one vector
	for (int i = 0; i < 6; i++) {
		dispG(i) = disp1(i);
		dispG(i + 6) = disp2(i);
		dispG(i + 12) = disp3(i);
		dispG(i + 18) = disp4(i);
	}

	// Convert nodal displacements from global to local cs
	dispL.addMatrixVector(0.0, T, dispG, 1.0);

	// In-plane forces from 2-node 6DOF MVLEM formulation
	double R1 = 0.0;
	double R2 = 0.0;
	double R3 = 0.0;
	double R4 = 0.0;
	double R5 = 0.0;
	double R6 = 0.0;

	// Get the current force matrix from panel stresses
	for (int i = 0; i < m; i++)
	{
		// Get the material stress
		const Vector &Stress = theMaterial[i]->getStress();

		double fx = Stress(0);
		double fy = Stress(1);
		double tauxy = Stress(2);

		Fx[i] = fx * AcX[i];
		Fy[i] = fy * AcY[i];
		Fxy[i] = tauxy * AcY[i];

	}

	// Build force vector 
	double Fh = 0.0; // Force in horizontal spring (at location c*h)
	double Fysum = 0.0; // Sum of vertical forces

	for (int i = 0; i < m; i++)
	{
		Fh += -1.0*Fxy[i];
		Fysum += Fy[i];
		FourNodeSFI_MVLEM3DRlocal[24 + i] = Fx[i]; // Force on internal (dummy) DOFs
	}

	R1 = Fh;
	R2 = -Fysum;
	R3 = -Fh*c*h;
	R4 = -Fh;
	R5 = Fysum;
	R6 = -Fh*(1.0 - c)*h;

	for (int i = 0; i<m; i++) {
		R3 -= Fy[i] * x[i];
		R6 += +Fy[i] * x[i];
	}

	//opserr << "K_drilling in R: " << K_drilling << "\n";
	//opserr << "dispL 5: " << dispL(5) << "\n";
	//opserr << "dispL 11: " << dispL(11) << "\n";
	// Calculate force vector in local cs
	FourNodeSFI_MVLEM3DRlocal(0) = R1 / 2.0 + (Aim*Eim*dispL(0)) / Lw - (Aim*Eim*dispL(6)) / Lw;
	FourNodeSFI_MVLEM3DRlocal(1) = R2 / 2.0 - (R3*d) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim*dispL(1)) / (Lw*Lw*Lw) + (6.0 * Eim*Iim*dispL(5)) / (Lw*Lw) - (12.0 * Eim*Iim*dispL(7)) / (Lw*Lw*Lw) + (6.0 * Eim*Iim*dispL(11)) / (Lw*Lw);
	FourNodeSFI_MVLEM3DRlocal(2) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(3) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h) - (h*h)*NUelastic + 5 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(4) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (dispL(4))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(10))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(22))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(16))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(5) = 0*K_drilling * (dispL(5) + 0.5*dispL(11)) + R3 / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim*dispL(1)) / (Lw*Lw) + (4.0 * Eim*Iim*dispL(5)) / Lw - (6.0 * Eim*Iim*dispL(7)) / (Lw*Lw) + (2.0 * Eim*Iim*dispL(11)) / Lw;
	FourNodeSFI_MVLEM3DRlocal(6) = R1 / 2.0 - (Aim*Eim*dispL(0)) / Lw + (Aim*Eim*dispL(6)) / Lw;
	FourNodeSFI_MVLEM3DRlocal(7) = R2 / 2.0 + (R3*d) / (2.0 * (d*d) + 2.0) - (12.0 * Eim*Iim*dispL(1)) / (Lw*Lw*Lw) - (6.0 * Eim*Iim*dispL(5)) / (Lw*Lw) + (12.0 * Eim*Iim*dispL(7)) / (Lw*Lw*Lw) - (6.0 * Eim*Iim*dispL(11)) / (Lw*Lw);
	FourNodeSFI_MVLEM3DRlocal(8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*(4.0 * (h*h)*NUelastic + (h*h) + 10 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(9) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(10) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (dispL(4))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(16))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (dispL(22))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(5.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(11) = 0*K_drilling * (dispL(11) + 0.5*dispL(5)) + R3 / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim*dispL(1)) / (Lw*Lw) + (2.0 * Eim*Iim*dispL(5)) / Lw - (6.0 * Eim*Iim*dispL(7)) / (Lw*Lw) + (4.0 * Eim*Iim*dispL(11)) / Lw;
	FourNodeSFI_MVLEM3DRlocal(12) = R4 / 2.0 + (Aim*Eim*dispL(12)) / Lw - (Aim*Eim*dispL(18)) / Lw;
	FourNodeSFI_MVLEM3DRlocal(13) = R5 / 2.0 - (R6*d) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim*dispL(13)) / (Lw*Lw*Lw) + (6.0 * Eim*Iim*dispL(17)) / (Lw*Lw) - (12.0 * Eim*Iim*dispL(19)) / (Lw*Lw*Lw) + (6.0 * Eim*Iim*dispL(23)) / (Lw*Lw);
	FourNodeSFI_MVLEM3DRlocal(14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(15) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (dispL(16))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(22))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(10))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (dispL(4))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(17) = 0*K_drilling * (dispL(17) + 0.5*dispL(23)) + R6 / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim*dispL(13)) / (Lw*Lw) + (4.0 * Eim*Iim*dispL(17)) / Lw - (6.0 * Eim*Iim*dispL(19)) / (Lw*Lw) + (2.0 * Eim*Iim*dispL(23)) / Lw;
	FourNodeSFI_MVLEM3DRlocal(18) = R4 / 2.0 - (Aim*Eim*dispL(12)) / Lw + (Aim*Eim*dispL(18)) / Lw;
	FourNodeSFI_MVLEM3DRlocal(19) = R5 / 2.0 + (R6*d) / (2.0 * (d*d) + 2.0) - (12.0 * Eim*Iim*dispL(13)) / (Lw*Lw*Lw) - (6.0 * Eim*Iim*dispL(17)) / (Lw*Lw) + (12.0 * Eim*Iim*dispL(19)) / (Lw*Lw*Lw) - (6.0 * Eim*Iim*dispL(23)) / (Lw*Lw);
	FourNodeSFI_MVLEM3DRlocal(20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(21) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(22) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (dispL(22))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(16))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(4))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(10))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeSFI_MVLEM3DRlocal(23) = 0*K_drilling * (dispL(23) + 0.5*dispL(17)) + R6 / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim*dispL(13)) / (Lw*Lw) + (2.0 * Eim*Iim*dispL(17)) / Lw - (6.0 * Eim*Iim*dispL(19)) / (Lw*Lw) + (4.0 * Eim*Iim*dispL(23)) / Lw;
	
	// Convert force vector from local to global cs
	FourNodeSFI_MVLEM3DR.addMatrixTransposeVector(0.0, T, FourNodeSFI_MVLEM3DRlocal, 1.0);

	/*opserr << "FourNodeSFI_MVLEM3DR :";
	opserr << "\n";
	for (int i = 0; i < 24+m; i++) {
	opserr << FourNodeSFI_MVLEM3DR(i) << " , ";
	opserr << "\n";
	}
	opserr << "\n" << "\n";*/
	//opserr << "  21		getResistingForce\n";
	// Return element force vector
	return FourNodeSFI_MVLEM3DR;
}

// Get resisting force incremenet from inertial forces
const Vector & FourNodeSFI_MVLEM3D::getResistingForceIncInertia()
{

	// compute the current resisting force
	this->getResistingForce();

	if (TotalMass != 0.0) {

		// Get nodal accelerations
		const Vector &accel1 = theNodes[0]->getTrialAccel();
		const Vector &accel2 = theNodes[1]->getTrialAccel();
		const Vector &accel3 = theNodes[2]->getTrialAccel();
		const Vector &accel4 = theNodes[3]->getTrialAccel();

		Vector accelG(24);
		accelG.Zero();
		Vector accelL(24);
		accelL.Zero();

		// pack all nodal accelerations in one vector
		for (int i = 0; i < 6; i++) {
			accelG(i) = accel1(i);
			accelG(i + 6) = accel2(i);
			accelG(i + 12) = accel3(i);
			accelG(i + 18) = accel4(i);
		}

		// transform nodal accelerations from global to local cs
		accelL.addMatrixVector(0.0, T, accelG, 1.0);

		// !!! SAME COMMENT AS ABOVE REGARDING NodeMass @@@ See derivation MATLAB code
		FourNodeSFI_MVLEM3DRlocal(0) += NodeMass / 4.0*accelL(0) + NodeMass / 4.0*accelL(6);
		FourNodeSFI_MVLEM3DRlocal(1) += NodeMass / 4.0*accelL(1) + NodeMass / 4.0*accelL(7);
		FourNodeSFI_MVLEM3DRlocal(2) += NodeMass / 4.0*accelL(2) + NodeMass / 4.0*accelL(8);
		FourNodeSFI_MVLEM3DRlocal(6) += NodeMass / 4.0*accelL(0) + NodeMass / 4.0*accelL(6);
		FourNodeSFI_MVLEM3DRlocal(7) += NodeMass / 4.0*accelL(1) + NodeMass / 4.0*accelL(7);
		FourNodeSFI_MVLEM3DRlocal(8) += NodeMass / 4.0*accelL(2) + NodeMass / 4.0*accelL(8);
		FourNodeSFI_MVLEM3DRlocal(12) += NodeMass / 4.0*accelL(12) + NodeMass / 4.0*accelL(18);
		FourNodeSFI_MVLEM3DRlocal(13) += NodeMass / 4.0*accelL(13) + NodeMass / 4.0*accelL(19);
		FourNodeSFI_MVLEM3DRlocal(14) += NodeMass / 4.0*accelL(14) + NodeMass / 4.0*accelL(20);
		FourNodeSFI_MVLEM3DRlocal(18) += NodeMass / 4.0*accelL(12) + NodeMass / 4.0*accelL(18);
		FourNodeSFI_MVLEM3DRlocal(19) += NodeMass / 4.0*accelL(13) + NodeMass / 4.0*accelL(19);
		FourNodeSFI_MVLEM3DRlocal(20) += NodeMass / 4.0*accelL(14) + NodeMass / 4.0*accelL(20);

		// convert forces from local to global cs
		FourNodeSFI_MVLEM3DR.addMatrixTransposeVector(1.0, T, FourNodeSFI_MVLEM3DRlocal, 1.0);

		// Add the damping forces if rayleigh damping
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			FourNodeSFI_MVLEM3DR += this->getRayleighDampingForces();

	}
	else {

		// Add the damping forces if rayleigh damping
		if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			FourNodeSFI_MVLEM3DR += this->getRayleighDampingForces();
	}
	//opserr << "  22		getResistingForceIncInertia\n";
	return FourNodeSFI_MVLEM3DR;
}

// !!! this function is not up to date - correct it @@@ DONE
// Send Self
int FourNodeSFI_MVLEM3D::sendSelf(int commitTag, Channel &theChannel)
{
	int res;
	int dataTag = this->getDbTag();

	static Vector data(5);  // One bigger than needed so no clash later

	data(0) = this->getTag();
	data(1) = m;
	data(2) = c;
	data(3) = NUelastic;
	data(4) = Tfactor;

	// FourNodeSFI_MVLEM3D then sends the tags of it's nodes
	res = theChannel.sendID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING FourNodeSFI_MVLEM3D::sendSelf() - failed to send ID\n";
		return -2;
	}

	// Send the material class tags
	ID matClassTags(m);
	for (int i = 0; i < m; i++)
		matClassTags(i) = theMaterial[i]->getClassTag();
	res = theChannel.sendID(0, commitTag, matClassTags);

	// Send the material models
	for (int i = 0; i < m; i++)
		theMaterial[i]->sendSelf(commitTag, theChannel);
	//opserr << "  23		sendSelf\n";
	return 0;
}

// !!! this function should be updated. @@@ DONE
// Receive Self
int FourNodeSFI_MVLEM3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res;
	int dataTag = this->getDbTag();

	// FourNodeSFI_MVLEM3D creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	// delete dynamic memory
	if (theMaterial != 0) {
		for (int i = 0; i < m; i++)
			if (theMaterial[i] != 0)
				delete theMaterial[i];
		delete[] theMaterial;
	}

	Vector data(5); // One bigger than needed so no clash later
	res = theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING FourNodeSFI_MVLEM3D::recvSelf() - failed to receive Vector\n";
		return -1;
	}

	this->setTag((int)data(0));
	m = data(1);
	c = data(2);
	NUelastic = data(3);
	Tfactor = data(4);

	// FourNodeSFI_MVLEM3D now receives the tags of it's four external nodes
	res = theChannel.recvID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING FourNodeSFI_MVLEM3D::recvSelf() - failed to receive ID\n";
		return -2;
	}

	// Receive the material class tags
	ID matClassTags(m);
	res = theChannel.recvID(0, commitTag, matClassTags);

	// Allocate memory for the uniaxial materials
	theMaterial = new NDMaterial*[m];
	if (theMaterial == 0) {
		opserr << "FourNodeSFI_MVLEM3D::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Receive the material models
	for (int i = 0; i < m; i++) {
		theMaterial[i] = theBroker.getNewNDMaterial(matClassTags(i));
		if (theMaterial[i] == 0) {
			opserr << "FourNodeSFI_MVLEM3D::recvSelf() - "
				<< "failed to get blank uniaxial material.\n";
			return -3;
		}
		theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
	}
	//opserr << "  24		recvSelf\n";
	return 0;
}

// !!! tripple check and demonstrate the it works @@@ DONE
// Display model
int FourNodeSFI_MVLEM3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	// First get the end points of the beam based on
	// the display factor (a measure of the distorted image)
	// get location of nodes
	/*const Vector &nd1Crds = theNodes[0]->getCrds();
	const Vector &nd2Crds = theNodes[1]->getCrds();
	const Vector &nd3Crds = theNodes[2]->getCrds();
	const Vector &nd4Crds = theNodes[3]->getCrds();*/

	static Vector Gv1(3);
	static Vector Gv2(3);
	static Vector Gv3(3);
	static Vector Gv4(3);
	Gv1.Zero();
	Gv2.Zero();
	Gv3.Zero();
	Gv4.Zero();

	// scalse the up based on thedisplay factor
	if (displayMode >= 0) {

		const Vector &end1Disp = theNodes[0]->getDisp();
		const Vector &end2Disp = theNodes[1]->getDisp();
		const Vector &end3Disp = theNodes[2]->getDisp();
		const Vector &end4Disp = theNodes[3]->getDisp();

		for (int i = 0; i < 3; i++) { // loop over coordinates (3 for 3D elements)

	// add displacement (multiplied with the displacement facotr) to the original node location to obtain current node location
			Gv1(i) = nd1Crds(i) + end1Disp(i)*fact;
			Gv2(i) = nd2Crds(i) + end2Disp(i)*fact;
			Gv3(i) = nd3Crds(i) + end3Disp(i)*fact;
			Gv4(i) = nd4Crds(i) + end4Disp(i)*fact;

		}

	}
	else {

		int mode = displayMode * -1;

		const Matrix &eigen1 = theNodes[0]->getEigenvectors();
		const Matrix &eigen2 = theNodes[1]->getEigenvectors();
		const Matrix &eigen3 = theNodes[2]->getEigenvectors();
		const Matrix &eigen4 = theNodes[3]->getEigenvectors();

		if (eigen1.noCols() >= mode) {

			for (int i = 0; i < 3; i++) {

				Gv1(i) = nd1Crds(i) + eigen1(i, mode - 1)*fact;
				Gv2(i) = nd2Crds(i) + eigen2(i, mode - 1)*fact;
				Gv3(i) = nd3Crds(i) + eigen3(i, mode - 1)*fact;
				Gv4(i) = nd4Crds(i) + eigen4(i, mode - 1)*fact;

			}

		}
		else {

			for (int i = 0; i < 3; i++) {

				Gv1(i) = nd1Crds(i);
				Gv2(i) = nd2Crds(i);
				Gv3(i) = nd3Crds(i);
				Gv4(i) = nd4Crds(i);

			}
		}
	}

	int error = 0;

	Vector RGB(3); // setiing up the colors
	RGB(0) = 0.0;
	RGB(1) = 1.0;
	RGB(2) = 1.0;

	// Add 2 vectors for top and bottom middle nodes
	Vector Gv1_(3); // Centrline node at the bottom
	Vector Gv2_(3); // Cetnrline node at the top
	Gv1_.Zero();
	Gv2_.Zero();

	// Calculate x, y and z coordinates of V1_1 and v1_2 based on x,y,z coordinates
	Gv1_(0) = 0.5 * (Gv1(0) + Gv2(0));
	Gv1_(1) = 0.5 * (Gv1(1) + Gv2(1));
	Gv1_(2) = 0.5 * (Gv1(2) + Gv2(2));

	Gv2_(0) = 0.5 * (Gv3(0) + Gv4(0));
	Gv2_(1) = 0.5 * (Gv3(1) + Gv4(1));
	Gv2_(2) = 0.5 * (Gv3(2) + Gv4(2));

	Vector Lv1_(3); // Centrline node at the bottom
	Vector Lv2_(3); // Cetnrline node at the top
	Lv1_.Zero();
	Lv2_.Zero();
	Lv1_.addMatrixVector(0.0, Tt, Gv1_, 1.0);
	Lv2_.addMatrixVector(0.0, Tt, Gv2_, 1.0);

	// Displaying Fibers
	for (int panel = 0; panel < m; panel++) // loop over m panels
	{

		Matrix NodePLotCrds(m, 13); // (panel id, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)

		// First set the quantity to be displayed at the nodes;
		// if displayMode is 1 through 3 we will plot material stresses otherwise 0.0
		static Vector values(1); // values of epsX to be plotted 

		values(0) = 0.0;

		if (displayMode < 4 && displayMode > 0) {
			const Vector &stress = theMaterial[panel]->getStrain();
			values(0) = stress(displayMode - 1);
		}

		// Determine the deformation - rotation - other is taken from v1, v2
		const Vector &end1Disp4G = theNodes[0]->getDisp();
		const Vector &end2Disp4G = theNodes[1]->getDisp();
		const Vector &end3Disp4G = theNodes[2]->getDisp();
		const Vector &end4Disp4G = theNodes[3]->getDisp();

		static Vector end1Disp4(6); end1Disp4.Zero();
		static Vector end2Disp4(6); end2Disp4.Zero();
		static Vector end3Disp4(6); end3Disp4.Zero();
		static Vector end4Disp4(6); end4Disp4.Zero();

		end1Disp4.addMatrixVector(0.0, T6, end1Disp4G, 1.0);
		end2Disp4.addMatrixVector(0.0, T6, end2Disp4G, 1.0);
		end3Disp4.addMatrixVector(0.0, T6, end3Disp4G, 1.0);
		end4Disp4.addMatrixVector(0.0, T6, end4Disp4G, 1.0);

		static Vector end1Disp(6);
		static Vector end2Disp(6);
		end1Disp.Zero();
		end2Disp.Zero();

		for (int i = 0; i < 4; i++) {
			end1Disp(i) = 0.5 * (end1Disp4(i) + end2Disp4(i));
			end2Disp(i) = 0.5 * (end3Disp4(i) + end4Disp4(i));
		}
		end1Disp(4) = end1Disp4(4) / (2.0 * (d*d) + 2.0) + end2Disp4(4) / (2.0 * (d*d) + 2.0) + (end1Disp4(2)*d) / (2.0 * (d*d) + 2.0) - (end2Disp4(2)*d) / (2.0 * (d*d) + 2.0);
		end1Disp(5) = end1Disp4(5) / (2.0 * (d*d) + 2.0) + end2Disp4(5) / (2.0 * (d*d) + 2.0) - (end1Disp4(1)*d) / (2.0 * (d*d) + 2.0) + (end2Disp4(1)*d) / (2.0 * (d*d) + 2.0);

		end2Disp(4) = end3Disp4(4) / (2.0 * (d*d) + 2.0) + end4Disp4(4) / (2.0 * (d*d) + 2.0) + (end3Disp4(2)*d) / (2.0 * (d*d) + 2.0) - (end4Disp4(2)*d) / (2.0 * (d*d) + 2.0);
		end2Disp(5) = end3Disp4(5) / (2.0 * (d*d) + 2.0) + end4Disp4(5) / (2.0 * (d*d) + 2.0) - (end3Disp4(1)*d) / (2.0 * (d*d) + 2.0) + (end4Disp4(1)*d) / (2.0 * (d*d) + 2.0);

		// Fiber nodes
		NodePLotCrds(panel, 0) = panel + 1; // panel id

		Vector LocCoord(3); LocCoord.Zero();
		Vector GlCoord(3); GlCoord.Zero();
		// Local node 1 - bottom left
		LocCoord(0) = Lv1_(0) + x[panel] - b[panel] / 2.0; // x 
		LocCoord(1) = Lv1_(1) + (x[panel] - b[panel] / 2.0)*end1Disp(5)*fact; // y
		LocCoord(2) = Lv1_(2) - (x[panel] - b[panel] / 2.0)*end1Disp(4)*fact; // z
		GlCoord.addMatrixTransposeVector(0.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 1) = GlCoord(0);
		NodePLotCrds(panel, 2) = GlCoord(1);
		NodePLotCrds(panel, 3) = GlCoord(2);
		LocCoord.Zero();
		GlCoord.Zero();
		// Local node 2 - bottom right
		LocCoord(0) = Lv1_(0) + x[panel] + b[panel] / 2.0; // x
		LocCoord(1) = Lv1_(1) + (x[panel] + b[panel] / 2.0)*end1Disp(5)*fact; // y
		LocCoord(2) = Lv1_(2) - (x[panel] + b[panel] / 2.0)*end1Disp(4)*fact; // z
		GlCoord.addMatrixTransposeVector(0.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 4) = GlCoord(0);
		NodePLotCrds(panel, 5) = GlCoord(1);
		NodePLotCrds(panel, 6) = GlCoord(2);
		LocCoord.Zero();
		GlCoord.Zero();
		// Local node 3 - top left
		LocCoord(0) = Lv2_(0) + x[panel] + b[panel] / 2.0; // x
		LocCoord(1) = Lv2_(1) + (x[panel] + b[panel] / 2.0)*end2Disp(5)*fact; // y
		LocCoord(2) = Lv2_(2) - (x[panel] + b[panel] / 2.0)*end2Disp(4)*fact; // z
		GlCoord.addMatrixTransposeVector(0.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 7) = GlCoord(0);
		NodePLotCrds(panel, 8) = GlCoord(1);
		NodePLotCrds(panel, 9) = GlCoord(2);
		LocCoord.Zero();
		GlCoord.Zero();
		// Local node 4 - top rigth
		LocCoord(0) = Lv2_(0) + x[panel] - b[panel] / 2.0; // x
		LocCoord(1) = Lv2_(1) + (x[panel] - b[panel] / 2.0)*end2Disp(5)*fact; // y
		LocCoord(2) = Lv2_(2) - (x[panel] - b[panel] / 2.0)*end2Disp(4)*fact; // z
		GlCoord.addMatrixTransposeVector(0.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 10) = GlCoord(0);
		NodePLotCrds(panel, 11) = GlCoord(1);
		NodePLotCrds(panel, 12) = GlCoord(2);

		Matrix coords(4, 3); // Temporary coordinates for plotting

		coords(0, 0) = NodePLotCrds(panel, 1); // node 1 x
		coords(1, 0) = NodePLotCrds(panel, 4); // node 2 x
		coords(2, 0) = NodePLotCrds(panel, 7); // node 3 x
		coords(3, 0) = NodePLotCrds(panel, 10);// node 4 x

		coords(0, 1) = NodePLotCrds(panel, 2); // node 1 y
		coords(1, 1) = NodePLotCrds(panel, 5); // node 2 y
		coords(2, 1) = NodePLotCrds(panel, 8); // node 3 y
		coords(3, 1) = NodePLotCrds(panel, 11); // node 4 y

		coords(0, 2) = NodePLotCrds(panel, 3); // node 1 z
		coords(1, 2) = NodePLotCrds(panel, 6); // node 2 z
		coords(2, 2) = NodePLotCrds(panel, 9); // node 3 z
		coords(3, 2) = NodePLotCrds(panel, 12); // node 4 z

		error += theViewer.drawPolygon(coords, values);

	}
	//opserr << "  25		displaySelf\n";
	return error;

}

// Print Element Information
void FourNodeSFI_MVLEM3D::Print(OPS_Stream &s, int flag)
{
	if (flag == 0) {
		s << "FourNodeSFI_MVLEM3D Element tag: " << this->getTag() << endln;
		s << "iNode: " << externalNodes(0) << ", jNode: " << externalNodes(1) << "lNode: " << externalNodes(2) << ", kNode: " << externalNodes(3) << endln;
		s << "Element height: " << h << endln;
		s << "Number of RC panel elements: " << m << endln;

		// get resisting forces in global system
		s << "Global resisting forces: " << this->getResistingForce_24DOF();

		for (int i = 0; i < m; i++) {
			s << "\nPanel #: " << i + 1 << endln;
			theMaterial[i]->Print(s, flag);
		}

	}
	else if (flag == 1) {
		// does nothing
	}
	//opserr << "  26		Print\n";
}

// Set element responses
Response *FourNodeSFI_MVLEM3D::setResponse(const char **argv, int argc, OPS_Stream &s)
{
	// !!! I thought we added in-plane forces only, or local forces or something like that. It needs to be added here and in MVLEM for easier processing @@@ DONE
	Response *theResponse = 0;

	s.tag("ElementOutput");
	s.attr("eleType", "FourNodeSFI_MVLEM3D");
	s.attr("eleTag", this->getTag());
	s.attr("node1", externalNodes[0]);
	s.attr("node2", externalNodes[1]);
	s.attr("node3", externalNodes[3]);
	s.attr("node4", externalNodes[2]);

	// Nodal forces in global cs
	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

		s.tag("ResponseType", "Fx_i");
		s.tag("ResponseType", "Fy_i");
		s.tag("ResponseType", "Fz_i");
		s.tag("ResponseType", "Mx_i");
		s.tag("ResponseType", "My_i");
		s.tag("ResponseType", "Mz_i");
		s.tag("ResponseType", "Fx_j");
		s.tag("ResponseType", "Fy_j");
		s.tag("ResponseType", "Fz_j");
		s.tag("ResponseType", "Mx_j");
		s.tag("ResponseType", "My_j");
		s.tag("ResponseType", "Mz_j");
		s.tag("ResponseType", "Fx_k");
		s.tag("ResponseType", "Fy_k");
		s.tag("ResponseType", "Fz_k");
		s.tag("ResponseType", "Mx_k");
		s.tag("ResponseType", "My_k");
		s.tag("ResponseType", "Mz_k");
		s.tag("ResponseType", "Fx_l");
		s.tag("ResponseType", "Fy_l");
		s.tag("ResponseType", "Fz_l");
		s.tag("ResponseType", "Mx_l");
		s.tag("ResponseType", "My_l");
		s.tag("ResponseType", "Mz_l");

		return theResponse = new ElementResponse(this, 1, Vector(24));

	}
	// Nodal forces in local cs
	else if (strcmp(argv[0], "forceL") == 0 || strcmp(argv[0], "forcesL") == 0 ||
		strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

		s.tag("ResponseType", "Fx_i");
		s.tag("ResponseType", "Fy_i");
		s.tag("ResponseType", "Fz_i");
		s.tag("ResponseType", "Mx_i");
		s.tag("ResponseType", "My_i");
		s.tag("ResponseType", "Mz_i");
		s.tag("ResponseType", "Fx_j");
		s.tag("ResponseType", "Fy_j");
		s.tag("ResponseType", "Fz_j");
		s.tag("ResponseType", "Mx_j");
		s.tag("ResponseType", "My_j");
		s.tag("ResponseType", "Mz_j");
		s.tag("ResponseType", "Fx_k");
		s.tag("ResponseType", "Fy_k");
		s.tag("ResponseType", "Fz_k");
		s.tag("ResponseType", "Mx_k");
		s.tag("ResponseType", "My_k");
		s.tag("ResponseType", "Mz_k");
		s.tag("ResponseType", "Fx_l");
		s.tag("ResponseType", "Fy_l");
		s.tag("ResponseType", "Fz_l");
		s.tag("ResponseType", "Mx_l");
		s.tag("ResponseType", "My_l");
		s.tag("ResponseType", "Mz_l");

		return theResponse = new ElementResponse(this, 2, Vector(24));

	}

	// Shear deformation
	else if (strcmp(argv[0], "ShearDef") == 0 || strcmp(argv[0], "sheardef") == 0) {

		s.tag("ResponseType", "Dsh");

		return theResponse = new ElementResponse(this, 3, 0.0);

	}

	// Element curvature
	else if (strcmp(argv[0], "Curvature") == 0 || strcmp(argv[0], "curvature") == 0) {

		s.tag("ResponseType", "fi");

		return theResponse = new ElementResponse(this, 4, 0.0);
	}

	// Material output
	else if (strcmp(argv[0], "RCpanel") == 0 || strcmp(argv[0], "RCPanel")
		|| strcmp(argv[0], "RC_panel") || strcmp(argv[0], "RC_Panel") == 0)
	{

		// Check if correct # of arguments passed
		if (argc != 3) {
			opserr << "WARNING: Number of recorder input for RC Panel is: " << argc - 1 << "; should be 2: panTag (one panel only: 1 to m) and $Response_Type.\n";
			return 0;
		}

		int matNum = atoi(argv[1]);

		s.tag("Material");
		s.attr("number", matNum);

		return theResponse = theMaterial[matNum - 1]->setResponse(&argv[argc - 1], argc - 2, s);

	}

	s.endTag();
	//opserr << "  27		setResponse\n";
	return 0;
}

// Get shear deformation
double FourNodeSFI_MVLEM3D::getShearDef(void)
{
	//opserr << "  28		getShearDef\n";
	return Dsh;
}

// Get curvature (from vertical strains)
double FourNodeSFI_MVLEM3D::getCurvature(void)
{
	double Curv;

	Curv = (FourNodeSFI_MVLEM3DStrainY[0] - FourNodeSFI_MVLEM3DStrainY[m - 1]) / (x[0] - x[m - 1]);
	//opserr << "  29		getCurvature\n";
	return Curv;
}

// Get global forces at 24 DOFs (top and bottom node)
Vector FourNodeSFI_MVLEM3D::getResistingForce_24DOF(void)
{

	for (int i = 0; i < 24; i++) {
		P_24DOF(i) = FourNodeSFI_MVLEM3DR(i);
	}
	//opserr << "  30		getResistingForce_24DOF\n";
	return P_24DOF;
}

// Obtain element responses
int FourNodeSFI_MVLEM3D::getResponse(int responseID, Information &eleInfo)
{
	// !!! modify accordingly @@@ DONE
	switch (responseID)
	{
	case 1:  // Global forces
		return eleInfo.setVector(this->getResistingForce_24DOF());

	case 2:  // Local forces
		return eleInfo.setVector(this->getResistingForce_24DOF_local());

	case 3:  // Shear deformation
		return eleInfo.setDouble(this->getShearDef());

	case 4:  // Curvature
		return eleInfo.setDouble(this->getCurvature());

	default:

		return 0;

	}
	//opserr << "  31		getResponse\n";
}

// Return element local forces 
Vector FourNodeSFI_MVLEM3D::getResistingForce_24DOF_local(void)
{
	for (int i = 0; i < 24; i++) {
		P_24DOF_local(i) = FourNodeSFI_MVLEM3DRlocal(i);
	}
	//opserr << "  32		getResistingForce_24DOF_local\n";
	return P_24DOF_local;
}

// Compute element transformation matrix
void  FourNodeSFI_MVLEM3D::setTransformationMatrix(void) {

	//T(24 + m, 24 + m); // !!! why is the size not defined for all matrices? same for MVLEM! @@@ because of internal nodes
	//T.Zero(); // element transformation matrix
	//Tt.Zero(); // 3 x 3 - basic transformation matrix
	//T6.Zero(); // 6 x 6 - nodal transformation matrix !!! is this correct? @@@ YES

	// !!! Check if you can store nodal coordinates as global variables and just use them - this way it is slowing down a lot @@@ DONE
	/*const Vector &nd1Crds = theNodes[0]->getCrds();
	const Vector &nd2Crds = theNodes[1]->getCrds();
	const Vector &nd3Crds = theNodes[2]->getCrds();
	const Vector &nd4Crds = theNodes[3]->getCrds();*/

	// Define local axis: 
	// x: user input [x1, x2, x3]
	// y: Nd1 -> Nd2
	// z: (x) x (y)

	// Vector components, magnitudes and iunit vectors
	double Xx, Xy, Xz, X_, Xex, Xey, Xez;
	double Yx, Yy, Yz, Y_, Yex, Yey, Yez;
	double Zex, Zey, Zez;

	Xx = nd2Crds(0) - nd1Crds(0);
	Xy = nd2Crds(1) - nd1Crds(1);
	Xz = nd2Crds(2) - nd1Crds(2);

	//opserr << " Xx: " << Xx;
	//opserr << " Xy: " << Xy;
	//opserr << " Xz: " << Xz;

	// Magnitude
	X_ = pow(pow(Xx, 2) + pow(Xy, 2) + pow(Xz, 2), 0.5);

	// unit x components
	Xex = Xx / X_;
	Xey = Xy / X_;
	Xez = Xz / X_;

	// !!! are you sure this is correct considering your node ordering? It should be perpendicular to x @@@ YES
	// k------l
	// |      |
	// |      |
	// i------j
	// Components of local Y axis
	Yx = nd3Crds(0) - nd1Crds(0);
	Yy = nd3Crds(1) - nd1Crds(1);
	Yz = nd3Crds(2) - nd1Crds(2);

	// Magnitude
	Y_ = pow(pow(Yx, 2) + pow(Yy, 2) + pow(Yz, 2), 0.5);

	// unit y components
	Yex = Yx / Y_;
	Yey = Yy / Y_;
	Yez = Yz / Y_;

	// (Ze) = (Xe) x (Ye)
	Zex = Xey * Yez - Xez * Yey;
	Zey = -(Xex*Yez - Xez * Yex);
	Zez = Xex * Yey - Xey * Yex;

	// Fill in transformation matrices 
	// 3 x 3 - basic matrix
	Tt(0, 0) = Xex;
	Tt(1, 0) = Yex;
	Tt(2, 0) = Zex;
	Tt(0, 1) = Xey;
	Tt(1, 1) = Yey;
	Tt(2, 1) = Zey;
	Tt(0, 2) = Xez;
	Tt(1, 2) = Yez;
	Tt(2, 2) = Zez;

	// 6 x 6
	for (int j = 0; j < 6; j += 3) {

		T6(j + 0, j + 0) = Xex;
		T6(j + 1, j + 0) = Yex;
		T6(j + 2, j + 0) = Zex;
		T6(j + 0, j + 1) = Xey;
		T6(j + 1, j + 1) = Yey;
		T6(j + 2, j + 1) = Zey;
		T6(j + 0, j + 2) = Xez;
		T6(j + 1, j + 2) = Yez;
		T6(j + 2, j + 2) = Zez;
	}

	// 24 x 24
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
	for (int i = 0; i < m; ++i)
		T(24 + i, 24 + i) = 1.0; // Diagonal terms accounting for horizontal stiffness
	//opserr << "  33		setTransformationMatrix\n";
};