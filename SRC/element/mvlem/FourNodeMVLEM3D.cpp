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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>
#include "FourNodeMVLEM3D.h"
#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <Message.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

// initialize the class wide variables
// global coordinate system
Matrix FourNodeMVLEM3D::FourNodeMVLEM3DK(24, 24);
Matrix FourNodeMVLEM3D::FourNodeMVLEM3DM(24, 24);
Matrix FourNodeMVLEM3D::FourNodeMVLEM3DD(24, 24);
Vector FourNodeMVLEM3D::FourNodeMVLEM3DR(24);

// local coordinate system
Matrix FourNodeMVLEM3D::FourNodeMVLEM3DKlocal(24, 24);
Matrix FourNodeMVLEM3D::FourNodeMVLEM3DDlocal(24, 24);
// Vector FourNodeMVLEM3D::FourNodeMVLEM3DRlocal(24); // !!! why is this one not defined and initiated in the same way as others? Why not static?
Matrix FourNodeMVLEM3D::FourNodeMVLEM3DMlocal(24, 24);

// typical constructor - !!! add optional parameters for stiffness of imaginary beam
FourNodeMVLEM3D::FourNodeMVLEM3D(int tag,
	double Dens,
	int Nd1, int Nd2, int Nd3, int Nd4,
	Node *theNd1, Node *theNd2, Node *theNd3, Node *theNd4,
	UniaxialMaterial **materialsConcrete,
	UniaxialMaterial **materialsSteel,
	UniaxialMaterial **materialsShear,
	double *Rho,
	double *thickness,
	double *width,
	int mm = 0,
	double cc = 0.0,
	double nn = 0.0,
	double tf = 0.0)

	:Element(tag, ELE_TAG_FourNodeMVLEM3D),
	density(Dens),
	externalNodes(4),
	theNd1(0),
	theNd2(0),
	theNd3(0),
	theNd4(0),
	theMaterialsConcrete(0), theMaterialsSteel(0), theMaterialsShear(0),
	theLoad(0), FourNodeMVLEM3DStrain(0),
	c(cc), m(mm), NUelastic(nn), Tfactor(tf), Eave(0.0),
	T(24, 24), T6(6, 6), Tt(3, 3), FourNodeMVLEM3DRlocal(24)
{
	// Fill with ZEROs all element matrices
	FourNodeMVLEM3DK.Zero();
	FourNodeMVLEM3DR.Zero();
	FourNodeMVLEM3DD.Zero();
	FourNodeMVLEM3DM.Zero();

	FourNodeMVLEM3DKlocal.Zero();
	FourNodeMVLEM3DDlocal.Zero();
	FourNodeMVLEM3DMlocal.Zero();

	NodeMass = 0.0;
	h = 0.0;
	d = 0.0;
	Lw = 0.0;
	Tave = 0.0;
	Eave = 0.0;

	// Fill in the ID containing external node info with node id's    
	if (externalNodes.Size() != 4)
		opserr << "FATAL FourNodeMVLEM3D::FourNodeMVLEM3D() - out of memory, could not create an ID of size 4\n";

	// Assign node tags to external nodes
	externalNodes(0) = Nd1;
	externalNodes(1) = Nd2;
	externalNodes(2) = Nd3;
	externalNodes(3) = Nd4;

	//Set node pointers to NULL
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;

	// Get coordinates of end nodes 
	nd1Crds = theNd1->getCrds();
	nd2Crds = theNd2->getCrds();
	nd3Crds = theNd3->getCrds();
	nd4Crds = theNd4->getCrds();

	// Check thickness and width input
	if (thickness == 0) {
		opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
			<< "Null thickness array passed.\n";
		exit(-1);
	}

	if (width == 0) {
		opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
			<< "Null width array passed.\n";
		exit(-1);
	}

	// Allocate memory for element arrays
	// Input parameters
	t = new double[m];
	b = new double[m];
	rho = new double[m];
	
	// Assign values from input
	for (int i = 0; i<m; i++) {
		t[i] = thickness[i];
		b[i] = width[i];
		rho[i] = Rho[i];
		Lw += b[i];
	}

	// Area of concrete and steel fibers
	Ac = new double[m];
	As = new double[m];

	// Stiffness of concrete and steel fibers
	Ec = new double[m];
	Es = new double[m];

	// Fiber stiffness (trial)
	ky = new double[m];
	kh = new double[1];

	stressC = new double[m];
	stressS = new double[m];

	// Fiber strains
	FourNodeMVLEM3DStrain = new double[m + 1];

	// Assign zero to element arrays
	for (int i = 0; i < m; i++) {

		Ac[i] = 0.0;
		As[i] = 0.0;

		ky[i] = 0.0;

		stressC[i] = 0.0;
		stressS[i] = 0.0;

		Ec[i] = 0.0;
		Es[i] = 0.0;

		FourNodeMVLEM3DStrain[i] = 0.0;
	}

	FourNodeMVLEM3DStrain[m] = 0.0;

	kh[0] = 0.0;

	// Check Concrete material input
	if (materialsConcrete == 0) {
		opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
			<< "null Concrete material array passed.\n";
		exit(-1);
	}

	// Check Steel material input
	if (materialsSteel == 0) {
		opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
			<< "null Steel material array passed.\n";
		exit(-1);
	}

	// Check Shear material input
	if (materialsShear == 0) {
		opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
			<< "null Shear material passed.\n";
		exit(-1);
	}

	// Allocate memory for the Concrete uniaxial materials
	theMaterialsConcrete = new UniaxialMaterial*[m];
	if (theMaterialsConcrete == 0) {
		opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
			<< "failed to allocate pointers for Concrete uniaxial materials.\n";
		exit(-1);
	}

	// Get copies of the Concrete uniaxial materials
	for (int i = 0; i < m; i++) {
		if (materialsConcrete[i] == 0) {
			opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
				"null uniaxial Concrete material pointer passed.\n";
			exit(-1);
		}

		theMaterialsConcrete[i] = materialsConcrete[i]->getCopy();

		if (theMaterialsConcrete[i] == 0) {
			opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
				<< "failed to copy Concrete uniaxial material.\n";
			exit(-1);
		}
	}

	// Allocate memory for the Steel uniaxial materials
	theMaterialsSteel = new UniaxialMaterial*[m];
	if (theMaterialsSteel == 0) {
		opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
			<< "failed to allocate pointers for Steel uniaxial materials.\n";
		exit(-1);
	}

	// Get copies of the uniaxial materials
	for (int i = 0; i < m; i++) {
		if (materialsSteel[i] == 0) {
			opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
				"null uniaxial Steel material pointer passed.\n";
			exit(-1);
		}

		theMaterialsSteel[i] = materialsSteel[i]->getCopy();

		if (theMaterialsSteel[i] == 0) {
			opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
				<< "failed to copy Steel uniaxial material.\n";
			exit(-1);
		}
	}

	// Allocate memory for the Shear uniaxial materials
	theMaterialsShear = new UniaxialMaterial*[1];
	if (theMaterialsShear == 0) {
		opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
			<< "failed to allocate pointers for Shear uniaxial materials.\n";
		exit(-1);
	}

	// Get copies of the uniaxial materials
	for (int i = 0; i < 1; i++) {
		if (materialsShear[i] == 0) {
			opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
				"null uniaxial Shear material pointer passed.\n";
			exit(-1);
		}

		theMaterialsShear[i] = materialsShear[i]->getCopy();

		if (theMaterialsShear[i] == 0) {
			opserr << "FourNodeMVLEM3D::FourNodeMVLEM3D() - "
				<< "failed to copy Shear uniaxial material.\n";
			exit(-1);
		}
	}

	// Revert to start
	this->revertToStart();
}

// Constructor which should be invoked by an FE_ObjectBroker only
FourNodeMVLEM3D::FourNodeMVLEM3D()
	:Element(0, ELE_TAG_FourNodeMVLEM3D),
	density(0.0),
	externalNodes(4),
	theMaterialsConcrete(0), theMaterialsSteel(0), theMaterialsShear(0), theLoad(0), FourNodeMVLEM3DStrain(0),
	h(0.0), c(0.0), m(0), NUelastic(0.0), Tfactor(0.0) 

{
	if (externalNodes.Size() != 4)
		opserr << "FATAL FourNodeMVLEM3D::FourNodeMVLEM3D() - out of memory, could not create an ID of size 2\n";
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;
}

//  Destructor - provided to clean up any memory
FourNodeMVLEM3D::~FourNodeMVLEM3D()
{
	// clean up the memory associated with the element, this is
	// memory the FourNodeMVLEM3D objects allocates and memory allocated 
	// by other objects that the FourNodeMVLEM3D object is responsible for 
	// cleaning up, i.e. the MaterialObject.
	if (theMaterialsConcrete != 0) {
		for (int i = 0; i < m; i++)
			if (theMaterialsConcrete[i] != 0)
				delete theMaterialsConcrete[i];
		delete[] theMaterialsConcrete;
	}

	if (theMaterialsSteel != 0) {
		for (int i = 0; i < m; i++)
			if (theMaterialsSteel[i] != 0)
				delete theMaterialsSteel[i];
		delete[] theMaterialsSteel;
	}

	if (theMaterialsShear != 0) {
		for (int i = 0; i < 1; i++)
			if (theMaterialsShear[i] != 0)
				delete theMaterialsShear[i];
		delete[] theMaterialsShear;
	}

	if (theLoad != 0)
		delete theLoad;

	if (x != 0)
		delete x;
	if (t != 0)
		delete t;
	if (b != 0)
		delete b;
	if (rho != 0)
		delete rho;
	if (Ac != 0)
		delete Ac;
	if (As != 0)
		delete As;
	if (ky != 0)
		delete ky;
	if (kh != 0)
		delete kh;
	if (Ec != 0)
		delete Ec;
	if (Es != 0)
		delete Es;
	if (stressC != 0)
		delete stressC;
	if (stressS != 0)
		delete stressS;
	if (FourNodeMVLEM3DStrain != 0)
		delete FourNodeMVLEM3DStrain;
	// !!! CHECK IF ALL DELETED
}

int FourNodeMVLEM3D::getNumExternalNodes(void) const
{
	return 4;
}

const ID &FourNodeMVLEM3D::getExternalNodes(void)
{
	return externalNodes;
}

Node **FourNodeMVLEM3D::getNodePtrs(void)
{
	return theNodes;
}

int FourNodeMVLEM3D::getNumDOF(void) {
	return 24;
}

void FourNodeMVLEM3D::setDomain(Domain *theDomain)
{
	// check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		return;
	}

	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;

	// first ensure nodes exist in Domain and set the node pointers
	int Nd1 = externalNodes(0);
	int Nd2 = externalNodes(1);
	int Nd3 = externalNodes(2);
	int Nd4 = externalNodes(3);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);
	theNodes[2] = theDomain->getNode(Nd3);
	theNodes[3] = theDomain->getNode(Nd4);

	if (theNodes[0] == 0) {
		opserr << "WARNING FourNodeMVLEM3D::setDomain() - at FourNodeMVLEM3D " << this->getTag() << " node " <<
			Nd1 << " does not exist in domain\n";
		return;  // don't go any further - otherwise segmentation fault
	}

	if (theNodes[1] == 0) {
		opserr << "WARNING FourNodeMVLEM3D::setDomain() - at FourNodeMVLEM3D " << this->getTag() << " node " <<
			Nd2 << " does not exist in domain\n";
		return;
	}

	if (theNodes[2] == 0) {
		opserr << "WARNING FourNodeMVLEM3D::setDomain() - at FourNodeMVLEM3D " << this->getTag() << " node " <<
			Nd3 << " does not exist in domain\n";
		return;
	}

	if (theNodes[3] == 0) {
		opserr << "WARNING FourNodeMVLEM3D::setDomain() - at FourNodeMVLEM3D " << this->getTag() << " node " <<
			Nd4 << " does not exist in domain\n";
		return;
	}

	// Call the DomainComponent class method THIS IS VERY IMPORTANT
	this->DomainComponent::setDomain(theDomain);

	// Ensure connected nodes have correct number of dof's
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();
	int dofNd3 = theNodes[2]->getNumberDOF();
	int dofNd4 = theNodes[3]->getNumberDOF();

	if ((dofNd1 != 6) || (dofNd2 != 6) || (dofNd3 != 6) || (dofNd4 != 6)) {
		opserr << "FourNodeMVLEM3D::setDomain(): 6 dof required at all nodes. " << dofNd1 << " provided at node 1, " << dofNd2 << " provided at node 2, "
			<< dofNd3 << " provided at node 3, " << dofNd4 << " provided at node 4";

	}

	// Calculate the element height based on distance between top and bottom nodes and perform checks
	double h1 = pow(pow(nd3Crds(0) - nd1Crds(0), 2.0) + pow(nd3Crds(1) - nd1Crds(1), 2.0) + pow(nd3Crds(2) - nd1Crds(2), 2.0), 0.5);
	double h2 = pow(pow(nd4Crds(0) - nd2Crds(0), 2.0) + pow(nd4Crds(1) - nd2Crds(1), 2.0) + pow(nd4Crds(2) - nd2Crds(2), 2.0), 0.5);

	// Check if element height is zero
	if ((h1 == 0.0) || (h2 == 0.0)) {
		opserr << "WARNING: FourNodeMVLEM3D element with tag " << this->getTag() <<
			" has ZERO height. Check geometry.";
		exit(-1);
	}

	// Check if element has constant height
	if ((h1 / h2 > 1.01) || (h1 / h2 < 0.99)) {
		opserr << "WARNING: FourNodeMVLEM3D element with tag " << this->getTag() << 
			" does not have constant height. Heights of the element are " << h1 << " and " << h2 << ". Check geometry.";
		exit(-1);
	}

	// Element height
	h = (h1 + h2) / 2.0;

	// Calculate average wall thickness
	for (int i = 0; i<m; i++) {
		Tave += t[i] * b[i] / (m * Lw); // !!! what is this values in an example. Check formula
	}

	// Calculate the element length based on distance between left and right nodes and perform checks
	double L1 = pow(pow(nd1Crds(0) - nd2Crds(0), 2.0) + pow(nd1Crds(1) - nd2Crds(1), 2.0) + pow(nd1Crds(2) - nd2Crds(2), 2.0), 0.5);
	double L2 = pow(pow(nd4Crds(0) - nd3Crds(0), 2.0) + pow(nd4Crds(1) - nd3Crds(1), 2.0) + pow(nd4Crds(2) - nd3Crds(2), 2.0), 0.5);

	// Check width of element
	if ((L1 / L2 > 1.01) || (L1 / L2 < 0.99)) {
		opserr << "WARNING: FourNodeMVLEM3D element with tag " << this->getTag() <<
			" does not have constant length. Top and bottom lengths of the element are " << L1 << " and " << L2 << ". Check geometry.";
		exit(-1);
	}
	
	if ((Lw / L1 > 1.01) || (Lw / L1 < 0.99)) {
		opserr << "WARNING: Node coordinates do not match sum of fiber widths for FourNodeMVLEM3D element with tag " << this->getTag() <<
			". Element width based on model geometry is " << L1 << " and sum of fiber widths is " << Lw << ". Check input and geometry.";
		exit(-1);
	}

	if ((Lw / L2 > 1.01) || (Lw / L2 < 0.99)) {
		opserr << "WARNING: Node coordinates do not match sum of fiber widths for FourNodeMVLEM3D element with tag " << this->getTag() <<
			". Element width based on model geometry is " << L2 << " and sum of fiber widths is " << Lw << ". Check input and geometry.";
		exit(-1);
	}

	// Calculate concrete and steel areas in Y directions
	for (int i = 0; i < m; i++) {
		As[i] = (b[i] * t[i])*rho[i];
		Ac[i] = (b[i] * t[i]) - As[i];
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

	// Calculate distance of corner nodes from center axis
	d = Lw / 2.0;

	// Determine the nodal mass for lumped mass approach
	A = 0;
	for (int i = 0; i < m; i++) {
		A += b[i] * t[i];
	}

	NodeMass = density * A * h / 4.0;

	// Calculate out-of-plane modulus of elasticity (average modulus)
	for (int i = 0; i < m; ++i) {
		Ec[i] = theMaterialsConcrete[i]->getInitialTangent();
		Eave += Ec[i] * b[i] * t[i] / A;
	}

	// Imaginary beam parameters
	Eim = 1.0e4 * Eave; // !!! implement optional input parameter instead of 1.0e4
	Him = 0.5 * h;		// !!! implement optional input parameter instead of 0.5
	Aim = Tave * Him;
	Iim = Tave * Him*Him*Him / 12.0;

	// Determine the transformation matrix
	setTransformationMatrix();

	// Create a vector to hop applied loads
	if (theLoad == 0)
		theLoad = new Vector(24);
	if (theLoad == 0) {
		opserr << "FourNodeMVLEM3D::setDomain() - element: " << this->getTag()
			<< " out of memory creating vector of size: " << 24 << endln;
		return;
	}
}

// Commit state of the materials
int FourNodeMVLEM3D::commitState() {

	int errCode = 0;

	// Commit Concrete material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsConcrete[i]->commitState();

	// Commit Steel material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsSteel[i]->commitState();

	// Commit Shear material models
	for (int i = 0; i < 1; i++)
		errCode += theMaterialsShear[i]->commitState();

	return errCode;
}

// Revert materials to last commit
int FourNodeMVLEM3D::revertToLastCommit() {

	int errCode = 0;

	// Revert Concrete material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsConcrete[i]->revertToLastCommit();

	// Revert Steel material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsSteel[i]->revertToLastCommit();

	// Revert Shear material model
	for (int i = 0; i < 1; i++)
		errCode += theMaterialsShear[i]->revertToLastCommit();

	return errCode;
}

// Revert materials to start
int FourNodeMVLEM3D::revertToStart() {

	int errCode = 0;

	// Revert Concrete material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsConcrete[i]->revertToStart();

	// Revert Steel material models
	for (int i = 0; i < m; i++)
		errCode += theMaterialsSteel[i]->revertToStart();

	// Revert Shear material model
	for (int i = 0; i < 1; i++)
		errCode += theMaterialsShear[i]->revertToStart();

	return errCode;
}

// Update the element
int FourNodeMVLEM3D::update() {

	// setTransformationMatrix(); !!! add switch for P-Delta effects - research a bit what is proper way to do this

	// Determine the current strain given trial displacements at nodes
	FourNodeMVLEM3DStrain = this->computeCurrentStrain();

	// Set the strain in the materials
	int errCode1 = 0;

	// Set trial response for Concrete material models
	for (int i = 0; i < m; i++)
		errCode1 += theMaterialsConcrete[i]->setTrialStrain(FourNodeMVLEM3DStrain[i]);

	// Set trial response for Steel material models
	for (int i = 0; i < m; i++)
		errCode1 += theMaterialsSteel[i]->setTrialStrain(FourNodeMVLEM3DStrain[i]);

	// Set trial response for Shear material model
	errCode1 += theMaterialsShear[0]->setTrialStrain(FourNodeMVLEM3DStrain[m]);

	return errCode1;
}

double * FourNodeMVLEM3D::computeCurrentStrain(void)
{

	// get nodal dispalcements in global cs
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();

	Vector dispG(24); // global cs
	Vector dispL(24); // local cs
	Vector dispL_inPlan2N(6); // local displacement vector of original 2-node 6DOF MVLEM formulation

	// store nodal displacemnts in global cs in one vector
	for (int i = 0; i < 6; i++) {
		dispG(i) = disp1(i);
		dispG(i + 6) = disp2(i);
		dispG(i + 12) = disp3(i);
		dispG(i + 18) = disp4(i);
	}

	// tranform nodal displacements from global to local cs
	dispL.addMatrixVector(0.0, T, dispG, 1.0);

	// Calculate 2-node 6DOF MVLEM local displacement vector
	dispL_inPlan2N(0) = dispL(0) / 2.0 + dispL(6) / 2.0;
	dispL_inPlan2N(1) = dispL(1) / 2.0 + dispL(7) / 2.0;
	dispL_inPlan2N(2) = dispL(5) / (2.0 * (d*d) + 2.0) + dispL(11) / (2.0 * (d*d) + 2.0) - (dispL(1)*d) / (2.0 * (d*d) + 2.0) + (dispL(7)*d) / (2.0 * (d*d) + 2.0);
	dispL_inPlan2N(3) = dispL(12) / 2.0 + dispL(18) / 2.0;
	dispL_inPlan2N(4) = dispL(13) / 2.0 + dispL(19) / 2.0;
	dispL_inPlan2N(5) = dispL(17) / (2.0 * (d*d) + 2.0) + dispL(23) / (2.0 * (d*d) + 2.0) - (dispL(13)*d) / (2.0 * (d*d) + 2.0) + (dispL(19)*d) / (2.0 * (d*d) + 2.0);

	// Fiber (Flexural) Strains
	for (int i = 0; i < m; i++) {
		FourNodeMVLEM3DStrain[i] = (-dispL_inPlan2N(1) - x[i] * dispL_inPlan2N(2) + dispL_inPlan2N(4) + x[i] * dispL_inPlan2N(5)) / h;	
	}

	// Shear deformation !!! deformation or strain? units? figure out how to make it more practical
	FourNodeMVLEM3DStrain[m] = dispL_inPlan2N(0) - dispL_inPlan2N(3) - c * h*dispL_inPlan2N(2) - (1.0 - c)*h*dispL_inPlan2N(5);

	return FourNodeMVLEM3DStrain;

}

const Matrix & FourNodeMVLEM3D::getInitialStiff(void)
{

	// Get vertical fiber materials initial tangent
	for (int i = 0; i < m; ++i)
	{
		Ec[i] = theMaterialsConcrete[i]->getInitialTangent();
		Es[i] = theMaterialsSteel[i]->getInitialTangent();
		ky[i] = Ec[i] * Ac[i] / h + Es[i] * As[i] / h;
	}

	// Build the initial stiffness matrix
	double Kv = 0.0; double Kh = 0.0; double Km = 0.0; double e = 0.0; double ex = 0.0;

	for (int i = 0; i < m; ++i)
	{
		Kv += ky[i];
		Km += ky[i] * x[i] * x[i];
		e += ky[i] * x[i];
	}

	// Get shear stiffness from shear material
	Kh = theMaterialsShear[0]->getInitialTangent();

	// Assemble element stiffness matrix !!! check these formulas in details against Matlab derivation
	FourNodeMVLEM3DKlocal(0, 0) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(0, 1) = (Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 2) = 0.0;
	FourNodeMVLEM3DKlocal(0, 3) = 0.0;
	FourNodeMVLEM3DKlocal(0, 4) = 0.0;
	FourNodeMVLEM3DKlocal(0, 5) = -(Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 6) = Kh / 4.0 - (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(0, 7) = -(Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 8) = 0.0;
	FourNodeMVLEM3DKlocal(0, 9) = 0.0;
	FourNodeMVLEM3DKlocal(0, 10) = 0.0;
	FourNodeMVLEM3DKlocal(0, 11) = -(Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 12) = -Kh / 4.0;
	FourNodeMVLEM3DKlocal(0, 13) = -(Kh*d*h*(c - 1)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 14) = 0.0;
	FourNodeMVLEM3DKlocal(0, 15) = 0.0;
	FourNodeMVLEM3DKlocal(0, 16) = 0.0;
	FourNodeMVLEM3DKlocal(0, 17) = (Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 18) = -Kh / 4.0;
	FourNodeMVLEM3DKlocal(0, 19) = (Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 20) = 0.0;
	FourNodeMVLEM3DKlocal(0, 21) = 0.0;
	FourNodeMVLEM3DKlocal(0, 22) = 0.0;
	FourNodeMVLEM3DKlocal(0, 23) = (Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);

	FourNodeMVLEM3DKlocal(1, 0) = FourNodeMVLEM3DKlocal(0, 1);
	FourNodeMVLEM3DKlocal(1, 1) = Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) - (d*(e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeMVLEM3DKlocal(1, 2) = 0.0;
	FourNodeMVLEM3DKlocal(1, 3) = 0.0;
	FourNodeMVLEM3DKlocal(1, 4) = 0.0;
	FourNodeMVLEM3DKlocal(1, 5) = (e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(1, 6) = (Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(1, 7) = Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) + (d*(e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) - (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeMVLEM3DKlocal(1, 8) = 0.0;
	FourNodeMVLEM3DKlocal(1, 9) = 0.0;
	FourNodeMVLEM3DKlocal(1, 10) = 0.0;
	FourNodeMVLEM3DKlocal(1, 11) = (e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(1, 12) = -(Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(1, 13) = (d*e) / (4.0 * (d*d) + 4.0) - Kv / 4.0 + (d*(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(1, 14) = 0.0;
	FourNodeMVLEM3DKlocal(1, 15) = 0.0;
	FourNodeMVLEM3DKlocal(1, 16) = 0.0;
	FourNodeMVLEM3DKlocal(1, 17) = -(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(1, 18) = -(Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(1, 19) = (d*e) / (4.0 * (d*d) + 4.0) - Kv / 4.0 - (d*(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(1, 20) = 0.0;
	FourNodeMVLEM3DKlocal(1, 21) = 0.0;
	FourNodeMVLEM3DKlocal(1, 22) = 0.0;
	FourNodeMVLEM3DKlocal(1, 23) = -(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);

	FourNodeMVLEM3DKlocal(2, 0) = FourNodeMVLEM3DKlocal(0, 2);
	FourNodeMVLEM3DKlocal(2, 1) = FourNodeMVLEM3DKlocal(1, 2);
	FourNodeMVLEM3DKlocal(2, 2) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 3) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 4) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 5) = 0.0;
	FourNodeMVLEM3DKlocal(2, 6) = 0.0;
	FourNodeMVLEM3DKlocal(2, 7) = 0.0;
	FourNodeMVLEM3DKlocal(2, 8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 9) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 10) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 11) = 0.0;
	FourNodeMVLEM3DKlocal(2, 12) = 0.0;
	FourNodeMVLEM3DKlocal(2, 13) = 0.0;
	FourNodeMVLEM3DKlocal(2, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 16) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 17) = 0.0;
	FourNodeMVLEM3DKlocal(2, 18) = 0.0;
	FourNodeMVLEM3DKlocal(2, 19) = 0.0;
	FourNodeMVLEM3DKlocal(2, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 23) = 0.0;

	FourNodeMVLEM3DKlocal(3, 0) = FourNodeMVLEM3DKlocal(0, 3);
	FourNodeMVLEM3DKlocal(3, 1) = FourNodeMVLEM3DKlocal(1, 3);
	FourNodeMVLEM3DKlocal(3, 2) = FourNodeMVLEM3DKlocal(2, 3);
	FourNodeMVLEM3DKlocal(3, 3) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 4) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeMVLEM3DKlocal(3, 5) = 0.0;
	FourNodeMVLEM3DKlocal(3, 6) = 0.0;
	FourNodeMVLEM3DKlocal(3, 7) = 0.0;
	FourNodeMVLEM3DKlocal(3, 8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 10) = 0.0;
	FourNodeMVLEM3DKlocal(3, 11) = 0.0;
	FourNodeMVLEM3DKlocal(3, 12) = 0.0;
	FourNodeMVLEM3DKlocal(3, 13) = 0.0;
	FourNodeMVLEM3DKlocal(3, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 16) = 0.0;
	FourNodeMVLEM3DKlocal(3, 17) = 0.0;
	FourNodeMVLEM3DKlocal(3, 18) = 0.0;
	FourNodeMVLEM3DKlocal(3, 19) = 0.0;
	FourNodeMVLEM3DKlocal(3, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 22) = 0.0;
	FourNodeMVLEM3DKlocal(3, 23) = 0.0;

	FourNodeMVLEM3DKlocal(4, 0) = FourNodeMVLEM3DKlocal(0, 4);
	FourNodeMVLEM3DKlocal(4, 1) = FourNodeMVLEM3DKlocal(1, 4);
	FourNodeMVLEM3DKlocal(4, 2) = FourNodeMVLEM3DKlocal(2, 4);
	FourNodeMVLEM3DKlocal(4, 3) = FourNodeMVLEM3DKlocal(3, 4);
	FourNodeMVLEM3DKlocal(4, 4) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 5) = 0.0;
	FourNodeMVLEM3DKlocal(4, 6) = 0.0;
	FourNodeMVLEM3DKlocal(4, 7) = 0.0;
	FourNodeMVLEM3DKlocal(4, 8) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 9) = 0.0;
	FourNodeMVLEM3DKlocal(4, 10) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 11) = 0.0;
	FourNodeMVLEM3DKlocal(4, 12) = 0.0;
	FourNodeMVLEM3DKlocal(4, 13) = 0.0;
	FourNodeMVLEM3DKlocal(4, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 15) = 0.0;
	FourNodeMVLEM3DKlocal(4, 16) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 17) = 0.0;
	FourNodeMVLEM3DKlocal(4, 18) = 0.0;
	FourNodeMVLEM3DKlocal(4, 19) = 0.0;
	FourNodeMVLEM3DKlocal(4, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 21) = 0.0;
	FourNodeMVLEM3DKlocal(4, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 23) = 0.0;

	FourNodeMVLEM3DKlocal(5, 0) = FourNodeMVLEM3DKlocal(0, 5);
	FourNodeMVLEM3DKlocal(5, 1) = FourNodeMVLEM3DKlocal(1, 5);
	FourNodeMVLEM3DKlocal(5, 2) = FourNodeMVLEM3DKlocal(2, 5);
	FourNodeMVLEM3DKlocal(5, 3) = FourNodeMVLEM3DKlocal(3, 5);
	FourNodeMVLEM3DKlocal(5, 4) = FourNodeMVLEM3DKlocal(4, 5);
	FourNodeMVLEM3DKlocal(5, 5) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw;
	FourNodeMVLEM3DKlocal(5, 6) = -(Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(5, 7) = e / (4.0 * (d*d) + 4.0) + (d*(Km + Kh*(c*c)*(h*h))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(5, 8) = 0.0;
	FourNodeMVLEM3DKlocal(5, 9) = 0.0;
	FourNodeMVLEM3DKlocal(5, 10) = 0.0;
	FourNodeMVLEM3DKlocal(5, 11) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (2.0 * Eim*Iim) / Lw;
	FourNodeMVLEM3DKlocal(5, 12) = (Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(5, 13) = (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - e / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(5, 14) = 0.0;
	FourNodeMVLEM3DKlocal(5, 15) = 0.0;
	FourNodeMVLEM3DKlocal(5, 16) = 0.0;
	FourNodeMVLEM3DKlocal(5, 17) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(5, 18) = (Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(5, 19) = -e / (4.0 * (d*d) + 4.0) - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(5, 20) = 0.0;
	FourNodeMVLEM3DKlocal(5, 21) = 0.0;
	FourNodeMVLEM3DKlocal(5, 22) = 0.0;
	FourNodeMVLEM3DKlocal(5, 23) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));

	FourNodeMVLEM3DKlocal(6, 0) = FourNodeMVLEM3DKlocal(0, 6);
	FourNodeMVLEM3DKlocal(6, 1) = FourNodeMVLEM3DKlocal(1, 6);
	FourNodeMVLEM3DKlocal(6, 2) = FourNodeMVLEM3DKlocal(2, 6);
	FourNodeMVLEM3DKlocal(6, 3) = FourNodeMVLEM3DKlocal(3, 6);
	FourNodeMVLEM3DKlocal(6, 4) = FourNodeMVLEM3DKlocal(4, 6);
	FourNodeMVLEM3DKlocal(6, 5) = FourNodeMVLEM3DKlocal(5, 6);
	FourNodeMVLEM3DKlocal(6, 6) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(6, 7) = -(Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 8) = 0.0;
	FourNodeMVLEM3DKlocal(6, 9) = 0.0;
	FourNodeMVLEM3DKlocal(6, 10) = 0.0;
	FourNodeMVLEM3DKlocal(6, 11) = -(Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 12) = -Kh / 4.0;
	FourNodeMVLEM3DKlocal(6, 13) = -(Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 14) = 0.0;
	FourNodeMVLEM3DKlocal(6, 15) = 0.0;
	FourNodeMVLEM3DKlocal(6, 16) = 0.0;
	FourNodeMVLEM3DKlocal(6, 17) = (Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 18) = -Kh / 4;
	FourNodeMVLEM3DKlocal(6, 19) = (Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 20) = 0.0;
	FourNodeMVLEM3DKlocal(6, 21) = 0.0;
	FourNodeMVLEM3DKlocal(6, 22) = 0.0;
	FourNodeMVLEM3DKlocal(6, 23) = (Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);

	FourNodeMVLEM3DKlocal(7, 0) = FourNodeMVLEM3DKlocal(0, 7);
	FourNodeMVLEM3DKlocal(7, 1) = FourNodeMVLEM3DKlocal(1, 7);
	FourNodeMVLEM3DKlocal(7, 2) = FourNodeMVLEM3DKlocal(2, 7);
	FourNodeMVLEM3DKlocal(7, 3) = FourNodeMVLEM3DKlocal(3, 7);
	FourNodeMVLEM3DKlocal(7, 4) = FourNodeMVLEM3DKlocal(4, 7);
	FourNodeMVLEM3DKlocal(7, 5) = FourNodeMVLEM3DKlocal(5, 7);
	FourNodeMVLEM3DKlocal(7, 6) = FourNodeMVLEM3DKlocal(6, 7);
	FourNodeMVLEM3DKlocal(7, 7) = Kv / 4.0 + (d*e) / (4.0 * (d*d) + 4.0) + (d*(e / 2.0 + (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeMVLEM3DKlocal(7, 8) = 0.0;
	FourNodeMVLEM3DKlocal(7, 9) = 0.0;
	FourNodeMVLEM3DKlocal(7, 10) = 0.0;
	FourNodeMVLEM3DKlocal(7, 11) = (e / 2.0 + (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) - (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(7, 12) = (Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(7, 13) = (d*(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) - (d*e) / (4.0 * (d*d) + 4.0) - Kv / 4.0;
	FourNodeMVLEM3DKlocal(7, 14) = 0.0;
	FourNodeMVLEM3DKlocal(7, 15) = 0.0;
	FourNodeMVLEM3DKlocal(7, 16) = 0.0;
	FourNodeMVLEM3DKlocal(7, 17) = -(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(7, 18) = (Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(7, 19) = -Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) - (d*(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(7, 20) = 0.0;
	FourNodeMVLEM3DKlocal(7, 21) = 0.0;
	FourNodeMVLEM3DKlocal(7, 22) = 0.0;
	FourNodeMVLEM3DKlocal(7, 23) = -(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);

	FourNodeMVLEM3DKlocal(8, 0) = FourNodeMVLEM3DKlocal(0, 8);
	FourNodeMVLEM3DKlocal(8, 1) = FourNodeMVLEM3DKlocal(1, 8);
	FourNodeMVLEM3DKlocal(8, 2) = FourNodeMVLEM3DKlocal(2, 8);
	FourNodeMVLEM3DKlocal(8, 3) = FourNodeMVLEM3DKlocal(3, 8);
	FourNodeMVLEM3DKlocal(8, 4) = FourNodeMVLEM3DKlocal(4, 8);
	FourNodeMVLEM3DKlocal(8, 5) = FourNodeMVLEM3DKlocal(5, 8);
	FourNodeMVLEM3DKlocal(8, 6) = FourNodeMVLEM3DKlocal(6, 8);
	FourNodeMVLEM3DKlocal(8, 7) = FourNodeMVLEM3DKlocal(7, 8);
	FourNodeMVLEM3DKlocal(8, 8) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 10) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 11) = 0.0;
	FourNodeMVLEM3DKlocal(8, 12) = 0.0;
	FourNodeMVLEM3DKlocal(8, 13) = 0.0;
	FourNodeMVLEM3DKlocal(8, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 16) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 17) = 0.0;
	FourNodeMVLEM3DKlocal(8, 18) = 0.0;
	FourNodeMVLEM3DKlocal(8, 19) = 0.0;
	FourNodeMVLEM3DKlocal(8, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 23) = 0.0;

	FourNodeMVLEM3DKlocal(9, 0) = FourNodeMVLEM3DKlocal(0, 9);
	FourNodeMVLEM3DKlocal(9, 1) = FourNodeMVLEM3DKlocal(1, 9);
	FourNodeMVLEM3DKlocal(9, 2) = FourNodeMVLEM3DKlocal(2, 9);
	FourNodeMVLEM3DKlocal(9, 3) = FourNodeMVLEM3DKlocal(3, 9);
	FourNodeMVLEM3DKlocal(9, 4) = FourNodeMVLEM3DKlocal(4, 9);
	FourNodeMVLEM3DKlocal(9, 5) = FourNodeMVLEM3DKlocal(5, 9);
	FourNodeMVLEM3DKlocal(9, 6) = FourNodeMVLEM3DKlocal(6, 9);
	FourNodeMVLEM3DKlocal(9, 7) = FourNodeMVLEM3DKlocal(7, 9);
	FourNodeMVLEM3DKlocal(9, 8) = FourNodeMVLEM3DKlocal(8, 9);
	FourNodeMVLEM3DKlocal(9, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 10) = -(Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeMVLEM3DKlocal(9, 11) = 0.0;
	FourNodeMVLEM3DKlocal(9, 12) = 0.0;
	FourNodeMVLEM3DKlocal(9, 13) = 0.0;
	FourNodeMVLEM3DKlocal(9, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 16) = 0.0;
	FourNodeMVLEM3DKlocal(9, 17) = 0.0;
	FourNodeMVLEM3DKlocal(9, 18) = 0.0;
	FourNodeMVLEM3DKlocal(9, 19) = 0.0;
	FourNodeMVLEM3DKlocal(9, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 22) = 0.0;
	FourNodeMVLEM3DKlocal(9, 23) = 0.0;

	FourNodeMVLEM3DKlocal(10, 0) = FourNodeMVLEM3DKlocal(0, 10);
	FourNodeMVLEM3DKlocal(10, 1) = FourNodeMVLEM3DKlocal(1, 10);
	FourNodeMVLEM3DKlocal(10, 2) = FourNodeMVLEM3DKlocal(2, 10);
	FourNodeMVLEM3DKlocal(10, 3) = FourNodeMVLEM3DKlocal(3, 10);
	FourNodeMVLEM3DKlocal(10, 4) = FourNodeMVLEM3DKlocal(4, 10);
	FourNodeMVLEM3DKlocal(10, 5) = FourNodeMVLEM3DKlocal(5, 10);
	FourNodeMVLEM3DKlocal(10, 6) = FourNodeMVLEM3DKlocal(6, 10);
	FourNodeMVLEM3DKlocal(10, 7) = FourNodeMVLEM3DKlocal(7, 10);
	FourNodeMVLEM3DKlocal(10, 8) = FourNodeMVLEM3DKlocal(8, 10);
	FourNodeMVLEM3DKlocal(10, 9) = FourNodeMVLEM3DKlocal(9, 10);
	FourNodeMVLEM3DKlocal(10, 10) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 11) = 0.0;
	FourNodeMVLEM3DKlocal(10, 12) = 0.0;
	FourNodeMVLEM3DKlocal(10, 13) = 0.0;
	FourNodeMVLEM3DKlocal(10, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 15) = 0.0;
	FourNodeMVLEM3DKlocal(10, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 17) = 0.0;
	FourNodeMVLEM3DKlocal(10, 18) = 0.0;
	FourNodeMVLEM3DKlocal(10, 19) = 0.0;
	FourNodeMVLEM3DKlocal(10, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 21) = 0.0;
	FourNodeMVLEM3DKlocal(10, 22) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 23) = 0.0;

	FourNodeMVLEM3DKlocal(11, 0) = FourNodeMVLEM3DKlocal(0, 11);
	FourNodeMVLEM3DKlocal(11, 1) = FourNodeMVLEM3DKlocal(1, 11);
	FourNodeMVLEM3DKlocal(11, 2) = FourNodeMVLEM3DKlocal(2, 11);
	FourNodeMVLEM3DKlocal(11, 3) = FourNodeMVLEM3DKlocal(3, 11);
	FourNodeMVLEM3DKlocal(11, 4) = FourNodeMVLEM3DKlocal(4, 11);
	FourNodeMVLEM3DKlocal(11, 5) = FourNodeMVLEM3DKlocal(5, 11);
	FourNodeMVLEM3DKlocal(11, 6) = FourNodeMVLEM3DKlocal(6, 11);
	FourNodeMVLEM3DKlocal(11, 7) = FourNodeMVLEM3DKlocal(7, 11);
	FourNodeMVLEM3DKlocal(11, 8) = FourNodeMVLEM3DKlocal(8, 11);
	FourNodeMVLEM3DKlocal(11, 9) = FourNodeMVLEM3DKlocal(9, 11);
	FourNodeMVLEM3DKlocal(11, 10) = FourNodeMVLEM3DKlocal(10, 11);
	FourNodeMVLEM3DKlocal(11, 11) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw;
	FourNodeMVLEM3DKlocal(11, 12) = (Kh*c*h) / (4 * (d*d) + 4);
	FourNodeMVLEM3DKlocal(11, 13) = (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - e / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(11, 14) = 0.0;
	FourNodeMVLEM3DKlocal(11, 15) = 0.0;
	FourNodeMVLEM3DKlocal(11, 16) = 0.0;
	FourNodeMVLEM3DKlocal(11, 17) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(11, 18) = (Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(11, 19) = -e / (4.0 * (d*d) + 4.0) - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(11, 20) = 0.0;
	FourNodeMVLEM3DKlocal(11, 21) = 0.0;
	FourNodeMVLEM3DKlocal(11, 22) = 0.0;
	FourNodeMVLEM3DKlocal(11, 23) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));

	FourNodeMVLEM3DKlocal(12, 0) = FourNodeMVLEM3DKlocal(0, 12);
	FourNodeMVLEM3DKlocal(12, 1) = FourNodeMVLEM3DKlocal(1, 12);
	FourNodeMVLEM3DKlocal(12, 2) = FourNodeMVLEM3DKlocal(2, 12);
	FourNodeMVLEM3DKlocal(12, 3) = FourNodeMVLEM3DKlocal(3, 12);
	FourNodeMVLEM3DKlocal(12, 4) = FourNodeMVLEM3DKlocal(4, 12);
	FourNodeMVLEM3DKlocal(12, 5) = FourNodeMVLEM3DKlocal(5, 12);
	FourNodeMVLEM3DKlocal(12, 6) = FourNodeMVLEM3DKlocal(6, 12);
	FourNodeMVLEM3DKlocal(12, 7) = FourNodeMVLEM3DKlocal(7, 12);
	FourNodeMVLEM3DKlocal(12, 8) = FourNodeMVLEM3DKlocal(8, 12);
	FourNodeMVLEM3DKlocal(12, 9) = FourNodeMVLEM3DKlocal(9, 12);
	FourNodeMVLEM3DKlocal(12, 10) = FourNodeMVLEM3DKlocal(10, 12);
	FourNodeMVLEM3DKlocal(12, 11) = FourNodeMVLEM3DKlocal(11, 12);
	FourNodeMVLEM3DKlocal(12, 12) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(12, 13) = (Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(12, 14) = 0.0;
	FourNodeMVLEM3DKlocal(12, 15) = 0.0;
	FourNodeMVLEM3DKlocal(12, 16) = 0.0;
	FourNodeMVLEM3DKlocal(12, 17) = -(Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(12, 18) = Kh / 4.0 - (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(12, 19) = -(Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(12, 20) = 0.0;
	FourNodeMVLEM3DKlocal(12, 21) = 0.0;
	FourNodeMVLEM3DKlocal(12, 22) = 0.0;
	FourNodeMVLEM3DKlocal(12, 23) = -(Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);

	FourNodeMVLEM3DKlocal(13, 0) = FourNodeMVLEM3DKlocal(0, 13);
	FourNodeMVLEM3DKlocal(13, 1) = FourNodeMVLEM3DKlocal(1, 13);
	FourNodeMVLEM3DKlocal(13, 2) = FourNodeMVLEM3DKlocal(2, 13);
	FourNodeMVLEM3DKlocal(13, 3) = FourNodeMVLEM3DKlocal(3, 13);
	FourNodeMVLEM3DKlocal(13, 4) = FourNodeMVLEM3DKlocal(4, 13);
	FourNodeMVLEM3DKlocal(13, 5) = FourNodeMVLEM3DKlocal(5, 13);
	FourNodeMVLEM3DKlocal(13, 6) = FourNodeMVLEM3DKlocal(6, 13);
	FourNodeMVLEM3DKlocal(13, 7) = FourNodeMVLEM3DKlocal(7, 13);
	FourNodeMVLEM3DKlocal(13, 8) = FourNodeMVLEM3DKlocal(8, 13);
	FourNodeMVLEM3DKlocal(13, 9) = FourNodeMVLEM3DKlocal(9, 13);
	FourNodeMVLEM3DKlocal(13, 10) = FourNodeMVLEM3DKlocal(10, 13);
	FourNodeMVLEM3DKlocal(13, 11) = FourNodeMVLEM3DKlocal(11, 13);
	FourNodeMVLEM3DKlocal(13, 12) = FourNodeMVLEM3DKlocal(12, 13);
	FourNodeMVLEM3DKlocal(13, 13) = Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw) - (d*(e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(13, 14) = 0.0;
	FourNodeMVLEM3DKlocal(13, 15) = 0.0;
	FourNodeMVLEM3DKlocal(13, 16) = 0.0;
	FourNodeMVLEM3DKlocal(13, 17) = (e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(13, 18) = (Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(13, 19) = Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) - (12.0 * Eim*Iim) / (Lw*Lw*Lw) + (d*(e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(13, 20) = 0.0;
	FourNodeMVLEM3DKlocal(13, 21) = 0.0;
	FourNodeMVLEM3DKlocal(13, 22) = 0.0;
	FourNodeMVLEM3DKlocal(13, 23) = (e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);

	FourNodeMVLEM3DKlocal(14, 0) = FourNodeMVLEM3DKlocal(0, 14);
	FourNodeMVLEM3DKlocal(14, 1) = FourNodeMVLEM3DKlocal(1, 14);
	FourNodeMVLEM3DKlocal(14, 2) = FourNodeMVLEM3DKlocal(2, 14);
	FourNodeMVLEM3DKlocal(14, 3) = FourNodeMVLEM3DKlocal(3, 14);
	FourNodeMVLEM3DKlocal(14, 4) = FourNodeMVLEM3DKlocal(4, 14);
	FourNodeMVLEM3DKlocal(14, 5) = FourNodeMVLEM3DKlocal(5, 14);
	FourNodeMVLEM3DKlocal(14, 6) = FourNodeMVLEM3DKlocal(6, 14);
	FourNodeMVLEM3DKlocal(14, 7) = FourNodeMVLEM3DKlocal(7, 14);
	FourNodeMVLEM3DKlocal(14, 8) = FourNodeMVLEM3DKlocal(8, 14);
	FourNodeMVLEM3DKlocal(14, 9) = FourNodeMVLEM3DKlocal(9, 14);
	FourNodeMVLEM3DKlocal(14, 10) = FourNodeMVLEM3DKlocal(10, 14);
	FourNodeMVLEM3DKlocal(14, 11) = FourNodeMVLEM3DKlocal(11, 14);
	FourNodeMVLEM3DKlocal(14, 12) = FourNodeMVLEM3DKlocal(12, 14);
	FourNodeMVLEM3DKlocal(14, 13) = FourNodeMVLEM3DKlocal(13, 14);
	FourNodeMVLEM3DKlocal(14, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 15) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 17) = 0.0;
	FourNodeMVLEM3DKlocal(14, 18) = 0.0;
	FourNodeMVLEM3DKlocal(14, 19) = 0.0;
	FourNodeMVLEM3DKlocal(14, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 23) = 0.0;

	FourNodeMVLEM3DKlocal(15, 0) = FourNodeMVLEM3DKlocal(0, 15);
	FourNodeMVLEM3DKlocal(15, 1) = FourNodeMVLEM3DKlocal(1, 15);
	FourNodeMVLEM3DKlocal(15, 2) = FourNodeMVLEM3DKlocal(2, 15);
	FourNodeMVLEM3DKlocal(15, 3) = FourNodeMVLEM3DKlocal(3, 15);
	FourNodeMVLEM3DKlocal(15, 4) = FourNodeMVLEM3DKlocal(4, 15);
	FourNodeMVLEM3DKlocal(15, 5) = FourNodeMVLEM3DKlocal(5, 15);
	FourNodeMVLEM3DKlocal(15, 6) = FourNodeMVLEM3DKlocal(6, 15);
	FourNodeMVLEM3DKlocal(15, 7) = FourNodeMVLEM3DKlocal(7, 15);
	FourNodeMVLEM3DKlocal(15, 8) = FourNodeMVLEM3DKlocal(8, 15);
	FourNodeMVLEM3DKlocal(15, 9) = FourNodeMVLEM3DKlocal(9, 15);
	FourNodeMVLEM3DKlocal(15, 10) = FourNodeMVLEM3DKlocal(10, 15);
	FourNodeMVLEM3DKlocal(15, 11) = FourNodeMVLEM3DKlocal(11, 15);
	FourNodeMVLEM3DKlocal(15, 12) = FourNodeMVLEM3DKlocal(12, 15);
	FourNodeMVLEM3DKlocal(15, 13) = FourNodeMVLEM3DKlocal(13, 15);
	FourNodeMVLEM3DKlocal(15, 14) = FourNodeMVLEM3DKlocal(14, 15);
	FourNodeMVLEM3DKlocal(15, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(15, 16) = -(Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeMVLEM3DKlocal(15, 17) = 0.0;
	FourNodeMVLEM3DKlocal(15, 18) = 0.0;
	FourNodeMVLEM3DKlocal(15, 19) = 0.0;
	FourNodeMVLEM3DKlocal(15, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(15, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(15, 22) = 0.0;
	FourNodeMVLEM3DKlocal(15, 23) = 0.0;

	FourNodeMVLEM3DKlocal(16, 0) = FourNodeMVLEM3DKlocal(0, 16);
	FourNodeMVLEM3DKlocal(16, 1) = FourNodeMVLEM3DKlocal(1, 16);
	FourNodeMVLEM3DKlocal(16, 2) = FourNodeMVLEM3DKlocal(2, 16);
	FourNodeMVLEM3DKlocal(16, 3) = FourNodeMVLEM3DKlocal(3, 16);
	FourNodeMVLEM3DKlocal(16, 4) = FourNodeMVLEM3DKlocal(4, 16);
	FourNodeMVLEM3DKlocal(16, 5) = FourNodeMVLEM3DKlocal(5, 16);
	FourNodeMVLEM3DKlocal(16, 6) = FourNodeMVLEM3DKlocal(6, 16);
	FourNodeMVLEM3DKlocal(16, 7) = FourNodeMVLEM3DKlocal(7, 16);
	FourNodeMVLEM3DKlocal(16, 8) = FourNodeMVLEM3DKlocal(8, 16);
	FourNodeMVLEM3DKlocal(16, 9) = FourNodeMVLEM3DKlocal(9, 16);
	FourNodeMVLEM3DKlocal(16, 10) = FourNodeMVLEM3DKlocal(10, 16);
	FourNodeMVLEM3DKlocal(16, 11) = FourNodeMVLEM3DKlocal(11, 16);
	FourNodeMVLEM3DKlocal(16, 12) = FourNodeMVLEM3DKlocal(12, 16);
	FourNodeMVLEM3DKlocal(16, 13) = FourNodeMVLEM3DKlocal(13, 16);
	FourNodeMVLEM3DKlocal(16, 14) = FourNodeMVLEM3DKlocal(14, 16);
	FourNodeMVLEM3DKlocal(16, 15) = FourNodeMVLEM3DKlocal(15, 16);
	FourNodeMVLEM3DKlocal(16, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(16, 17) = 0.0;
	FourNodeMVLEM3DKlocal(16, 18) = 0.0;
	FourNodeMVLEM3DKlocal(16, 19) = 0.0;
	FourNodeMVLEM3DKlocal(16, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(16, 21) = 0.0;
	FourNodeMVLEM3DKlocal(16, 22) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(16, 23) = 0.0;

	FourNodeMVLEM3DKlocal(17, 0) = FourNodeMVLEM3DKlocal(0, 17);
	FourNodeMVLEM3DKlocal(17, 1) = FourNodeMVLEM3DKlocal(1, 17);
	FourNodeMVLEM3DKlocal(17, 2) = FourNodeMVLEM3DKlocal(2, 17);
	FourNodeMVLEM3DKlocal(17, 3) = FourNodeMVLEM3DKlocal(3, 17);
	FourNodeMVLEM3DKlocal(17, 4) = FourNodeMVLEM3DKlocal(4, 17);
	FourNodeMVLEM3DKlocal(17, 5) = FourNodeMVLEM3DKlocal(5, 17);
	FourNodeMVLEM3DKlocal(17, 6) = FourNodeMVLEM3DKlocal(6, 17);
	FourNodeMVLEM3DKlocal(17, 7) = FourNodeMVLEM3DKlocal(7, 17);
	FourNodeMVLEM3DKlocal(17, 8) = FourNodeMVLEM3DKlocal(8, 17);
	FourNodeMVLEM3DKlocal(17, 9) = FourNodeMVLEM3DKlocal(9, 17);
	FourNodeMVLEM3DKlocal(17, 10) = FourNodeMVLEM3DKlocal(10, 17);
	FourNodeMVLEM3DKlocal(17, 11) = FourNodeMVLEM3DKlocal(11, 17);
	FourNodeMVLEM3DKlocal(17, 12) = FourNodeMVLEM3DKlocal(12, 17);
	FourNodeMVLEM3DKlocal(17, 13) = FourNodeMVLEM3DKlocal(13, 17);
	FourNodeMVLEM3DKlocal(17, 14) = FourNodeMVLEM3DKlocal(14, 17);
	FourNodeMVLEM3DKlocal(17, 15) = FourNodeMVLEM3DKlocal(15, 17);
	FourNodeMVLEM3DKlocal(17, 16) = FourNodeMVLEM3DKlocal(16, 17);
	FourNodeMVLEM3DKlocal(17, 17) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw;
	FourNodeMVLEM3DKlocal(17, 18) = -(Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(17, 19) = e / (4.0 * (d*d) + 4.0) - (6.0 * Eim*Iim) / (Lw*Lw) + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(17, 20) = 0.0;
	FourNodeMVLEM3DKlocal(17, 21) = 0.0;
	FourNodeMVLEM3DKlocal(17, 22) = 0.0;
	FourNodeMVLEM3DKlocal(17, 23) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (2.0 * Eim*Iim) / Lw;

	FourNodeMVLEM3DKlocal(18, 0) = FourNodeMVLEM3DKlocal(0, 18);
	FourNodeMVLEM3DKlocal(18, 1) = FourNodeMVLEM3DKlocal(1, 18);
	FourNodeMVLEM3DKlocal(18, 2) = FourNodeMVLEM3DKlocal(2, 18);
	FourNodeMVLEM3DKlocal(18, 3) = FourNodeMVLEM3DKlocal(3, 18);
	FourNodeMVLEM3DKlocal(18, 4) = FourNodeMVLEM3DKlocal(4, 18);
	FourNodeMVLEM3DKlocal(18, 5) = FourNodeMVLEM3DKlocal(5, 18);
	FourNodeMVLEM3DKlocal(18, 6) = FourNodeMVLEM3DKlocal(6, 18);
	FourNodeMVLEM3DKlocal(18, 7) = FourNodeMVLEM3DKlocal(7, 18);
	FourNodeMVLEM3DKlocal(18, 8) = FourNodeMVLEM3DKlocal(8, 18);
	FourNodeMVLEM3DKlocal(18, 9) = FourNodeMVLEM3DKlocal(9, 18);
	FourNodeMVLEM3DKlocal(18, 10) = FourNodeMVLEM3DKlocal(10, 18);
	FourNodeMVLEM3DKlocal(18, 11) = FourNodeMVLEM3DKlocal(11, 18);
	FourNodeMVLEM3DKlocal(18, 12) = FourNodeMVLEM3DKlocal(12, 18);
	FourNodeMVLEM3DKlocal(18, 13) = FourNodeMVLEM3DKlocal(13, 18);
	FourNodeMVLEM3DKlocal(18, 14) = FourNodeMVLEM3DKlocal(14, 18);
	FourNodeMVLEM3DKlocal(18, 15) = FourNodeMVLEM3DKlocal(15, 18);
	FourNodeMVLEM3DKlocal(18, 16) = FourNodeMVLEM3DKlocal(16, 18);
	FourNodeMVLEM3DKlocal(18, 17) = FourNodeMVLEM3DKlocal(17, 18);
	FourNodeMVLEM3DKlocal(18, 18) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(18, 19) = -(Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(18, 20) = 0.0;
	FourNodeMVLEM3DKlocal(18, 21) = 0.0;
	FourNodeMVLEM3DKlocal(18, 22) = 0.0;
	FourNodeMVLEM3DKlocal(18, 23) = -(Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);

	FourNodeMVLEM3DKlocal(19, 0) = FourNodeMVLEM3DKlocal(0, 19);
	FourNodeMVLEM3DKlocal(19, 1) = FourNodeMVLEM3DKlocal(1, 19);
	FourNodeMVLEM3DKlocal(19, 2) = FourNodeMVLEM3DKlocal(2, 19);
	FourNodeMVLEM3DKlocal(19, 3) = FourNodeMVLEM3DKlocal(3, 19);
	FourNodeMVLEM3DKlocal(19, 4) = FourNodeMVLEM3DKlocal(4, 19);
	FourNodeMVLEM3DKlocal(19, 5) = FourNodeMVLEM3DKlocal(5, 19);
	FourNodeMVLEM3DKlocal(19, 6) = FourNodeMVLEM3DKlocal(6, 19);
	FourNodeMVLEM3DKlocal(19, 7) = FourNodeMVLEM3DKlocal(7, 19);
	FourNodeMVLEM3DKlocal(19, 8) = FourNodeMVLEM3DKlocal(8, 19);
	FourNodeMVLEM3DKlocal(19, 9) = FourNodeMVLEM3DKlocal(9, 19);
	FourNodeMVLEM3DKlocal(19, 10) = FourNodeMVLEM3DKlocal(10, 19);
	FourNodeMVLEM3DKlocal(19, 11) = FourNodeMVLEM3DKlocal(11, 19);
	FourNodeMVLEM3DKlocal(19, 12) = FourNodeMVLEM3DKlocal(12, 19);
	FourNodeMVLEM3DKlocal(19, 13) = FourNodeMVLEM3DKlocal(13, 19);
	FourNodeMVLEM3DKlocal(19, 14) = FourNodeMVLEM3DKlocal(14, 19);
	FourNodeMVLEM3DKlocal(19, 15) = FourNodeMVLEM3DKlocal(15, 19);
	FourNodeMVLEM3DKlocal(19, 16) = FourNodeMVLEM3DKlocal(16, 19);
	FourNodeMVLEM3DKlocal(19, 17) = FourNodeMVLEM3DKlocal(17, 19);
	FourNodeMVLEM3DKlocal(19, 18) = FourNodeMVLEM3DKlocal(18, 19);
	FourNodeMVLEM3DKlocal(19, 19) = Kv / 4.0 + (d*e) / (4.0 * (d*d) + 4.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw) + (d*(e / 2.0 + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(19, 20) = 0.0;
	FourNodeMVLEM3DKlocal(19, 21) = 0.0;
	FourNodeMVLEM3DKlocal(19, 22) = 0.0;
	FourNodeMVLEM3DKlocal(19, 23) = (e / 2.0 + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) - (6.0 * Eim*Iim) / (Lw*Lw);

	FourNodeMVLEM3DKlocal(20, 0) = FourNodeMVLEM3DKlocal(0, 20);
	FourNodeMVLEM3DKlocal(20, 1) = FourNodeMVLEM3DKlocal(1, 20);
	FourNodeMVLEM3DKlocal(20, 2) = FourNodeMVLEM3DKlocal(2, 20);
	FourNodeMVLEM3DKlocal(20, 3) = FourNodeMVLEM3DKlocal(3, 20);
	FourNodeMVLEM3DKlocal(20, 4) = FourNodeMVLEM3DKlocal(4, 20);
	FourNodeMVLEM3DKlocal(20, 5) = FourNodeMVLEM3DKlocal(5, 20);
	FourNodeMVLEM3DKlocal(20, 6) = FourNodeMVLEM3DKlocal(6, 20);
	FourNodeMVLEM3DKlocal(20, 7) = FourNodeMVLEM3DKlocal(7, 20);
	FourNodeMVLEM3DKlocal(20, 8) = FourNodeMVLEM3DKlocal(8, 20);
	FourNodeMVLEM3DKlocal(20, 9) = FourNodeMVLEM3DKlocal(9, 20);
	FourNodeMVLEM3DKlocal(20, 10) = FourNodeMVLEM3DKlocal(10, 20);
	FourNodeMVLEM3DKlocal(20, 11) = FourNodeMVLEM3DKlocal(11, 20);
	FourNodeMVLEM3DKlocal(20, 12) = FourNodeMVLEM3DKlocal(12, 20);
	FourNodeMVLEM3DKlocal(20, 13) = FourNodeMVLEM3DKlocal(13, 20);
	FourNodeMVLEM3DKlocal(20, 14) = FourNodeMVLEM3DKlocal(14, 20);
	FourNodeMVLEM3DKlocal(20, 15) = FourNodeMVLEM3DKlocal(15, 20);
	FourNodeMVLEM3DKlocal(20, 16) = FourNodeMVLEM3DKlocal(16, 20);
	FourNodeMVLEM3DKlocal(20, 17) = FourNodeMVLEM3DKlocal(17, 20);
	FourNodeMVLEM3DKlocal(20, 18) = FourNodeMVLEM3DKlocal(18, 20);
	FourNodeMVLEM3DKlocal(20, 19) = FourNodeMVLEM3DKlocal(19, 20);
	FourNodeMVLEM3DKlocal(20, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(20, 21) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(20, 22) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(20, 23) = 0.0;

	FourNodeMVLEM3DKlocal(21, 0) = FourNodeMVLEM3DKlocal(0, 21);
	FourNodeMVLEM3DKlocal(21, 1) = FourNodeMVLEM3DKlocal(1, 21);
	FourNodeMVLEM3DKlocal(21, 2) = FourNodeMVLEM3DKlocal(2, 21);
	FourNodeMVLEM3DKlocal(21, 3) = FourNodeMVLEM3DKlocal(3, 21);
	FourNodeMVLEM3DKlocal(21, 4) = FourNodeMVLEM3DKlocal(4, 21);
	FourNodeMVLEM3DKlocal(21, 5) = FourNodeMVLEM3DKlocal(5, 21);
	FourNodeMVLEM3DKlocal(21, 6) = FourNodeMVLEM3DKlocal(6, 21);
	FourNodeMVLEM3DKlocal(21, 7) = FourNodeMVLEM3DKlocal(7, 21);
	FourNodeMVLEM3DKlocal(21, 8) = FourNodeMVLEM3DKlocal(8, 21);
	FourNodeMVLEM3DKlocal(21, 9) = FourNodeMVLEM3DKlocal(9, 21);
	FourNodeMVLEM3DKlocal(21, 10) = FourNodeMVLEM3DKlocal(10, 21);
	FourNodeMVLEM3DKlocal(21, 11) = FourNodeMVLEM3DKlocal(11, 21);
	FourNodeMVLEM3DKlocal(21, 12) = FourNodeMVLEM3DKlocal(12, 21);
	FourNodeMVLEM3DKlocal(21, 13) = FourNodeMVLEM3DKlocal(13, 21);
	FourNodeMVLEM3DKlocal(21, 14) = FourNodeMVLEM3DKlocal(14, 21);
	FourNodeMVLEM3DKlocal(21, 15) = FourNodeMVLEM3DKlocal(15, 21);
	FourNodeMVLEM3DKlocal(21, 16) = FourNodeMVLEM3DKlocal(16, 21);
	FourNodeMVLEM3DKlocal(21, 17) = FourNodeMVLEM3DKlocal(17, 21);
	FourNodeMVLEM3DKlocal(21, 18) = FourNodeMVLEM3DKlocal(18, 21);
	FourNodeMVLEM3DKlocal(21, 19) = FourNodeMVLEM3DKlocal(19, 21);
	FourNodeMVLEM3DKlocal(21, 20) = FourNodeMVLEM3DKlocal(20, 21);
	FourNodeMVLEM3DKlocal(21, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(21, 22) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeMVLEM3DKlocal(21, 23) = 0.0;

	FourNodeMVLEM3DKlocal(22, 0) = FourNodeMVLEM3DKlocal(0, 22);
	FourNodeMVLEM3DKlocal(22, 1) = FourNodeMVLEM3DKlocal(1, 22);
	FourNodeMVLEM3DKlocal(22, 2) = FourNodeMVLEM3DKlocal(2, 22);
	FourNodeMVLEM3DKlocal(22, 3) = FourNodeMVLEM3DKlocal(3, 22);
	FourNodeMVLEM3DKlocal(22, 4) = FourNodeMVLEM3DKlocal(4, 22);
	FourNodeMVLEM3DKlocal(22, 5) = FourNodeMVLEM3DKlocal(5, 22);
	FourNodeMVLEM3DKlocal(22, 6) = FourNodeMVLEM3DKlocal(6, 22);
	FourNodeMVLEM3DKlocal(22, 7) = FourNodeMVLEM3DKlocal(7, 22);
	FourNodeMVLEM3DKlocal(22, 8) = FourNodeMVLEM3DKlocal(8, 22);
	FourNodeMVLEM3DKlocal(22, 9) = FourNodeMVLEM3DKlocal(9, 22);
	FourNodeMVLEM3DKlocal(22, 10) = FourNodeMVLEM3DKlocal(10, 22);
	FourNodeMVLEM3DKlocal(22, 11) = FourNodeMVLEM3DKlocal(11, 22);
	FourNodeMVLEM3DKlocal(22, 12) = FourNodeMVLEM3DKlocal(12, 22);
	FourNodeMVLEM3DKlocal(22, 13) = FourNodeMVLEM3DKlocal(13, 22);
	FourNodeMVLEM3DKlocal(22, 14) = FourNodeMVLEM3DKlocal(14, 22);
	FourNodeMVLEM3DKlocal(22, 15) = FourNodeMVLEM3DKlocal(15, 22);
	FourNodeMVLEM3DKlocal(22, 16) = FourNodeMVLEM3DKlocal(16, 22);
	FourNodeMVLEM3DKlocal(22, 17) = FourNodeMVLEM3DKlocal(17, 22);
	FourNodeMVLEM3DKlocal(22, 18) = FourNodeMVLEM3DKlocal(18, 22);
	FourNodeMVLEM3DKlocal(22, 19) = FourNodeMVLEM3DKlocal(19, 22);
	FourNodeMVLEM3DKlocal(22, 20) = FourNodeMVLEM3DKlocal(20, 22);
	FourNodeMVLEM3DKlocal(22, 21) = FourNodeMVLEM3DKlocal(21, 22);
	FourNodeMVLEM3DKlocal(22, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(22, 23) = 0.0;

	FourNodeMVLEM3DKlocal(23, 0) = FourNodeMVLEM3DKlocal(0, 23);
	FourNodeMVLEM3DKlocal(23, 1) = FourNodeMVLEM3DKlocal(1, 23);
	FourNodeMVLEM3DKlocal(23, 2) = FourNodeMVLEM3DKlocal(2, 23);
	FourNodeMVLEM3DKlocal(23, 3) = FourNodeMVLEM3DKlocal(3, 23);
	FourNodeMVLEM3DKlocal(23, 4) = FourNodeMVLEM3DKlocal(4, 23);
	FourNodeMVLEM3DKlocal(23, 5) = FourNodeMVLEM3DKlocal(5, 23);
	FourNodeMVLEM3DKlocal(23, 6) = FourNodeMVLEM3DKlocal(6, 23);
	FourNodeMVLEM3DKlocal(23, 7) = FourNodeMVLEM3DKlocal(7, 23);
	FourNodeMVLEM3DKlocal(23, 8) = FourNodeMVLEM3DKlocal(8, 23);
	FourNodeMVLEM3DKlocal(23, 9) = FourNodeMVLEM3DKlocal(9, 23);
	FourNodeMVLEM3DKlocal(23, 10) = FourNodeMVLEM3DKlocal(10, 23);
	FourNodeMVLEM3DKlocal(23, 11) = FourNodeMVLEM3DKlocal(11, 23);
	FourNodeMVLEM3DKlocal(23, 12) = FourNodeMVLEM3DKlocal(12, 23);
	FourNodeMVLEM3DKlocal(23, 13) = FourNodeMVLEM3DKlocal(13, 23);
	FourNodeMVLEM3DKlocal(23, 14) = FourNodeMVLEM3DKlocal(14, 23);
	FourNodeMVLEM3DKlocal(23, 15) = FourNodeMVLEM3DKlocal(15, 23);
	FourNodeMVLEM3DKlocal(23, 16) = FourNodeMVLEM3DKlocal(16, 23);
	FourNodeMVLEM3DKlocal(23, 17) = FourNodeMVLEM3DKlocal(17, 23);
	FourNodeMVLEM3DKlocal(23, 18) = FourNodeMVLEM3DKlocal(18, 23);
	FourNodeMVLEM3DKlocal(23, 19) = FourNodeMVLEM3DKlocal(19, 23);
	FourNodeMVLEM3DKlocal(23, 20) = FourNodeMVLEM3DKlocal(20, 23);
	FourNodeMVLEM3DKlocal(23, 21) = FourNodeMVLEM3DKlocal(21, 23);
	FourNodeMVLEM3DKlocal(23, 22) = FourNodeMVLEM3DKlocal(22, 23);
	FourNodeMVLEM3DKlocal(23, 23) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw;

	// Convert matrix from local to global cs !!! check 0.0 versus 1.0 - make consistent everywhere
	FourNodeMVLEM3DK.addMatrixTripleProduct(0.0, T, FourNodeMVLEM3DKlocal, 1.0); 

	// Return element stiffness matrix
	return FourNodeMVLEM3DK;

}


const Matrix & FourNodeMVLEM3D::getTangentStiff(void)
{

	for (int i = 0; i < m; ++i)
	{
		Ec[i] = theMaterialsConcrete[i]->getTangent();
		Es[i] = theMaterialsSteel[i]->getTangent();
		ky[i] = Ec[i] * Ac[i] / h + Es[i] * As[i] / h;
	}

	// Build the initial stiffness matrix
	double Kv = 0.0; double Kh = 0.0; double Km = 0.0; double e = 0.0; double ex = 0.0;

	for (int i = 0; i < m; ++i)
	{
		Kv += ky[i];
		Km += ky[i] * x[i] * x[i];
		e += ky[i] * x[i];

	}

	// Get shear stiffness from shear material
	Kh = theMaterialsShear[0]->getTangent();

	// Assemble element stiffness matrix
	FourNodeMVLEM3DKlocal(0, 0) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(0, 1) = (Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 2) = 0.0;
	FourNodeMVLEM3DKlocal(0, 3) = 0.0;
	FourNodeMVLEM3DKlocal(0, 4) = 0.0;
	FourNodeMVLEM3DKlocal(0, 5) = -(Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 6) = Kh / 4.0 - (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(0, 7) = -(Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 8) = 0.0;
	FourNodeMVLEM3DKlocal(0, 9) = 0.0;
	FourNodeMVLEM3DKlocal(0, 10) = 0.0;
	FourNodeMVLEM3DKlocal(0, 11) = -(Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 12) = -Kh / 4.0;
	FourNodeMVLEM3DKlocal(0, 13) = -(Kh*d*h*(c - 1)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 14) = 0.0;
	FourNodeMVLEM3DKlocal(0, 15) = 0.0;
	FourNodeMVLEM3DKlocal(0, 16) = 0.0;
	FourNodeMVLEM3DKlocal(0, 17) = (Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 18) = -Kh / 4.0;
	FourNodeMVLEM3DKlocal(0, 19) = (Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(0, 20) = 0.0;
	FourNodeMVLEM3DKlocal(0, 21) = 0.0;
	FourNodeMVLEM3DKlocal(0, 22) = 0.0;
	FourNodeMVLEM3DKlocal(0, 23) = (Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);

	FourNodeMVLEM3DKlocal(1, 0) = FourNodeMVLEM3DKlocal(0, 1);
	FourNodeMVLEM3DKlocal(1, 1) = Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) - (d*(e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeMVLEM3DKlocal(1, 2) = 0.0;
	FourNodeMVLEM3DKlocal(1, 3) = 0.0;
	FourNodeMVLEM3DKlocal(1, 4) = 0.0;
	FourNodeMVLEM3DKlocal(1, 5) = (e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(1, 6) = (Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(1, 7) = Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) + (d*(e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) - (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeMVLEM3DKlocal(1, 8) = 0.0;
	FourNodeMVLEM3DKlocal(1, 9) = 0.0;
	FourNodeMVLEM3DKlocal(1, 10) = 0.0;
	FourNodeMVLEM3DKlocal(1, 11) = (e / 2.0 - (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(1, 12) = -(Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(1, 13) = (d*e) / (4.0 * (d*d) + 4.0) - Kv / 4.0 + (d*(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(1, 14) = 0.0;
	FourNodeMVLEM3DKlocal(1, 15) = 0.0;
	FourNodeMVLEM3DKlocal(1, 16) = 0.0;
	FourNodeMVLEM3DKlocal(1, 17) = -(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(1, 18) = -(Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(1, 19) = (d*e) / (4.0 * (d*d) + 4.0) - Kv / 4.0 - (d*(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(1, 20) = 0.0;
	FourNodeMVLEM3DKlocal(1, 21) = 0.0;
	FourNodeMVLEM3DKlocal(1, 22) = 0.0;
	FourNodeMVLEM3DKlocal(1, 23) = -(e / 2.0 - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);

	FourNodeMVLEM3DKlocal(2, 0) = FourNodeMVLEM3DKlocal(0, 2);
	FourNodeMVLEM3DKlocal(2, 1) = FourNodeMVLEM3DKlocal(1, 2);
	FourNodeMVLEM3DKlocal(2, 2) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 3) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 4) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 5) = 0.0;
	FourNodeMVLEM3DKlocal(2, 6) = 0.0;
	FourNodeMVLEM3DKlocal(2, 7) = 0.0;
	FourNodeMVLEM3DKlocal(2, 8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 9) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 10) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 11) = 0.0;
	FourNodeMVLEM3DKlocal(2, 12) = 0.0;
	FourNodeMVLEM3DKlocal(2, 13) = 0.0;
	FourNodeMVLEM3DKlocal(2, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 16) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 17) = 0.0;
	FourNodeMVLEM3DKlocal(2, 18) = 0.0;
	FourNodeMVLEM3DKlocal(2, 19) = 0.0;
	FourNodeMVLEM3DKlocal(2, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(2, 23) = 0.0;

	FourNodeMVLEM3DKlocal(3, 0) = FourNodeMVLEM3DKlocal(0, 3);
	FourNodeMVLEM3DKlocal(3, 1) = FourNodeMVLEM3DKlocal(1, 3);
	FourNodeMVLEM3DKlocal(3, 2) = FourNodeMVLEM3DKlocal(2, 3);
	FourNodeMVLEM3DKlocal(3, 3) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 4) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeMVLEM3DKlocal(3, 5) = 0.0;
	FourNodeMVLEM3DKlocal(3, 6) = 0.0;
	FourNodeMVLEM3DKlocal(3, 7) = 0.0;
	FourNodeMVLEM3DKlocal(3, 8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 10) = 0.0;
	FourNodeMVLEM3DKlocal(3, 11) = 0.0;
	FourNodeMVLEM3DKlocal(3, 12) = 0.0;
	FourNodeMVLEM3DKlocal(3, 13) = 0.0;
	FourNodeMVLEM3DKlocal(3, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 16) = 0.0;
	FourNodeMVLEM3DKlocal(3, 17) = 0.0;
	FourNodeMVLEM3DKlocal(3, 18) = 0.0;
	FourNodeMVLEM3DKlocal(3, 19) = 0.0;
	FourNodeMVLEM3DKlocal(3, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(3, 22) = 0.0;
	FourNodeMVLEM3DKlocal(3, 23) = 0.0;

	FourNodeMVLEM3DKlocal(4, 0) = FourNodeMVLEM3DKlocal(0, 4);
	FourNodeMVLEM3DKlocal(4, 1) = FourNodeMVLEM3DKlocal(1, 4);
	FourNodeMVLEM3DKlocal(4, 2) = FourNodeMVLEM3DKlocal(2, 4);
	FourNodeMVLEM3DKlocal(4, 3) = FourNodeMVLEM3DKlocal(3, 4);
	FourNodeMVLEM3DKlocal(4, 4) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 5) = 0.0;
	FourNodeMVLEM3DKlocal(4, 6) = 0.0;
	FourNodeMVLEM3DKlocal(4, 7) = 0.0;
	FourNodeMVLEM3DKlocal(4, 8) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 9) = 0.0;
	FourNodeMVLEM3DKlocal(4, 10) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 11) = 0.0;
	FourNodeMVLEM3DKlocal(4, 12) = 0.0;
	FourNodeMVLEM3DKlocal(4, 13) = 0.0;
	FourNodeMVLEM3DKlocal(4, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 15) = 0.0;
	FourNodeMVLEM3DKlocal(4, 16) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 17) = 0.0;
	FourNodeMVLEM3DKlocal(4, 18) = 0.0;
	FourNodeMVLEM3DKlocal(4, 19) = 0.0;
	FourNodeMVLEM3DKlocal(4, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 21) = 0.0;
	FourNodeMVLEM3DKlocal(4, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(4, 23) = 0.0;

	FourNodeMVLEM3DKlocal(5, 0) = FourNodeMVLEM3DKlocal(0, 5);
	FourNodeMVLEM3DKlocal(5, 1) = FourNodeMVLEM3DKlocal(1, 5);
	FourNodeMVLEM3DKlocal(5, 2) = FourNodeMVLEM3DKlocal(2, 5);
	FourNodeMVLEM3DKlocal(5, 3) = FourNodeMVLEM3DKlocal(3, 5);
	FourNodeMVLEM3DKlocal(5, 4) = FourNodeMVLEM3DKlocal(4, 5);
	FourNodeMVLEM3DKlocal(5, 5) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw;
	FourNodeMVLEM3DKlocal(5, 6) = -(Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(5, 7) = e / (4.0 * (d*d) + 4.0) + (d*(Km + Kh*(c*c)*(h*h))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(5, 8) = 0.0;
	FourNodeMVLEM3DKlocal(5, 9) = 0.0;
	FourNodeMVLEM3DKlocal(5, 10) = 0.0;
	FourNodeMVLEM3DKlocal(5, 11) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (2.0 * Eim*Iim) / Lw;
	FourNodeMVLEM3DKlocal(5, 12) = (Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(5, 13) = (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - e / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(5, 14) = 0.0;
	FourNodeMVLEM3DKlocal(5, 15) = 0.0;
	FourNodeMVLEM3DKlocal(5, 16) = 0.0;
	FourNodeMVLEM3DKlocal(5, 17) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(5, 18) = (Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(5, 19) = -e / (4.0 * (d*d) + 4.0) - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(5, 20) = 0.0;
	FourNodeMVLEM3DKlocal(5, 21) = 0.0;
	FourNodeMVLEM3DKlocal(5, 22) = 0.0;
	FourNodeMVLEM3DKlocal(5, 23) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));

	FourNodeMVLEM3DKlocal(6, 0) = FourNodeMVLEM3DKlocal(0, 6);
	FourNodeMVLEM3DKlocal(6, 1) = FourNodeMVLEM3DKlocal(1, 6);
	FourNodeMVLEM3DKlocal(6, 2) = FourNodeMVLEM3DKlocal(2, 6);
	FourNodeMVLEM3DKlocal(6, 3) = FourNodeMVLEM3DKlocal(3, 6);
	FourNodeMVLEM3DKlocal(6, 4) = FourNodeMVLEM3DKlocal(4, 6);
	FourNodeMVLEM3DKlocal(6, 5) = FourNodeMVLEM3DKlocal(5, 6);
	FourNodeMVLEM3DKlocal(6, 6) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(6, 7) = -(Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 8) = 0.0;
	FourNodeMVLEM3DKlocal(6, 9) = 0.0;
	FourNodeMVLEM3DKlocal(6, 10) = 0.0;
	FourNodeMVLEM3DKlocal(6, 11) = -(Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 12) = -Kh / 4.0;
	FourNodeMVLEM3DKlocal(6, 13) = -(Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 14) = 0.0;
	FourNodeMVLEM3DKlocal(6, 15) = 0.0;
	FourNodeMVLEM3DKlocal(6, 16) = 0.0;
	FourNodeMVLEM3DKlocal(6, 17) = (Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 18) = -Kh / 4;
	FourNodeMVLEM3DKlocal(6, 19) = (Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(6, 20) = 0.0;
	FourNodeMVLEM3DKlocal(6, 21) = 0.0;
	FourNodeMVLEM3DKlocal(6, 22) = 0.0;
	FourNodeMVLEM3DKlocal(6, 23) = (Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);

	FourNodeMVLEM3DKlocal(7, 0) = FourNodeMVLEM3DKlocal(0, 7);
	FourNodeMVLEM3DKlocal(7, 1) = FourNodeMVLEM3DKlocal(1, 7);
	FourNodeMVLEM3DKlocal(7, 2) = FourNodeMVLEM3DKlocal(2, 7);
	FourNodeMVLEM3DKlocal(7, 3) = FourNodeMVLEM3DKlocal(3, 7);
	FourNodeMVLEM3DKlocal(7, 4) = FourNodeMVLEM3DKlocal(4, 7);
	FourNodeMVLEM3DKlocal(7, 5) = FourNodeMVLEM3DKlocal(5, 7);
	FourNodeMVLEM3DKlocal(7, 6) = FourNodeMVLEM3DKlocal(6, 7);
	FourNodeMVLEM3DKlocal(7, 7) = Kv / 4.0 + (d*e) / (4.0 * (d*d) + 4.0) + (d*(e / 2.0 + (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw);
	FourNodeMVLEM3DKlocal(7, 8) = 0.0;
	FourNodeMVLEM3DKlocal(7, 9) = 0.0;
	FourNodeMVLEM3DKlocal(7, 10) = 0.0;
	FourNodeMVLEM3DKlocal(7, 11) = (e / 2.0 + (d*(Km + Kh*(c*c)*(h*h))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) - (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(7, 12) = (Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(7, 13) = (d*(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0) - (d*e) / (4.0 * (d*d) + 4.0) - Kv / 4.0;
	FourNodeMVLEM3DKlocal(7, 14) = 0.0;
	FourNodeMVLEM3DKlocal(7, 15) = 0.0;
	FourNodeMVLEM3DKlocal(7, 16) = 0.0;
	FourNodeMVLEM3DKlocal(7, 17) = -(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(7, 18) = (Kh*c*d*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(7, 19) = -Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) - (d*(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(7, 20) = 0.0;
	FourNodeMVLEM3DKlocal(7, 21) = 0.0;
	FourNodeMVLEM3DKlocal(7, 22) = 0.0;
	FourNodeMVLEM3DKlocal(7, 23) = -(e / 2.0 + (d*(Km + Kh*c*(h*h)*(c - 1.0))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0);

	FourNodeMVLEM3DKlocal(8, 0) = FourNodeMVLEM3DKlocal(0, 8);
	FourNodeMVLEM3DKlocal(8, 1) = FourNodeMVLEM3DKlocal(1, 8);
	FourNodeMVLEM3DKlocal(8, 2) = FourNodeMVLEM3DKlocal(2, 8);
	FourNodeMVLEM3DKlocal(8, 3) = FourNodeMVLEM3DKlocal(3, 8);
	FourNodeMVLEM3DKlocal(8, 4) = FourNodeMVLEM3DKlocal(4, 8);
	FourNodeMVLEM3DKlocal(8, 5) = FourNodeMVLEM3DKlocal(5, 8);
	FourNodeMVLEM3DKlocal(8, 6) = FourNodeMVLEM3DKlocal(6, 8);
	FourNodeMVLEM3DKlocal(8, 7) = FourNodeMVLEM3DKlocal(7, 8);
	FourNodeMVLEM3DKlocal(8, 8) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 10) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 11) = 0.0;
	FourNodeMVLEM3DKlocal(8, 12) = 0.0;
	FourNodeMVLEM3DKlocal(8, 13) = 0.0;
	FourNodeMVLEM3DKlocal(8, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 16) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 17) = 0.0;
	FourNodeMVLEM3DKlocal(8, 18) = 0.0;
	FourNodeMVLEM3DKlocal(8, 19) = 0.0;
	FourNodeMVLEM3DKlocal(8, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(8, 23) = 0.0;

	FourNodeMVLEM3DKlocal(9, 0) = FourNodeMVLEM3DKlocal(0, 9);
	FourNodeMVLEM3DKlocal(9, 1) = FourNodeMVLEM3DKlocal(1, 9);
	FourNodeMVLEM3DKlocal(9, 2) = FourNodeMVLEM3DKlocal(2, 9);
	FourNodeMVLEM3DKlocal(9, 3) = FourNodeMVLEM3DKlocal(3, 9);
	FourNodeMVLEM3DKlocal(9, 4) = FourNodeMVLEM3DKlocal(4, 9);
	FourNodeMVLEM3DKlocal(9, 5) = FourNodeMVLEM3DKlocal(5, 9);
	FourNodeMVLEM3DKlocal(9, 6) = FourNodeMVLEM3DKlocal(6, 9);
	FourNodeMVLEM3DKlocal(9, 7) = FourNodeMVLEM3DKlocal(7, 9);
	FourNodeMVLEM3DKlocal(9, 8) = FourNodeMVLEM3DKlocal(8, 9);
	FourNodeMVLEM3DKlocal(9, 9) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 10) = -(Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeMVLEM3DKlocal(9, 11) = 0.0;
	FourNodeMVLEM3DKlocal(9, 12) = 0.0;
	FourNodeMVLEM3DKlocal(9, 13) = 0.0;
	FourNodeMVLEM3DKlocal(9, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 16) = 0.0;
	FourNodeMVLEM3DKlocal(9, 17) = 0.0;
	FourNodeMVLEM3DKlocal(9, 18) = 0.0;
	FourNodeMVLEM3DKlocal(9, 19) = 0.0;
	FourNodeMVLEM3DKlocal(9, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(9, 22) = 0.0;
	FourNodeMVLEM3DKlocal(9, 23) = 0.0;

	FourNodeMVLEM3DKlocal(10, 0) = FourNodeMVLEM3DKlocal(0, 10);
	FourNodeMVLEM3DKlocal(10, 1) = FourNodeMVLEM3DKlocal(1, 10);
	FourNodeMVLEM3DKlocal(10, 2) = FourNodeMVLEM3DKlocal(2, 10);
	FourNodeMVLEM3DKlocal(10, 3) = FourNodeMVLEM3DKlocal(3, 10);
	FourNodeMVLEM3DKlocal(10, 4) = FourNodeMVLEM3DKlocal(4, 10);
	FourNodeMVLEM3DKlocal(10, 5) = FourNodeMVLEM3DKlocal(5, 10);
	FourNodeMVLEM3DKlocal(10, 6) = FourNodeMVLEM3DKlocal(6, 10);
	FourNodeMVLEM3DKlocal(10, 7) = FourNodeMVLEM3DKlocal(7, 10);
	FourNodeMVLEM3DKlocal(10, 8) = FourNodeMVLEM3DKlocal(8, 10);
	FourNodeMVLEM3DKlocal(10, 9) = FourNodeMVLEM3DKlocal(9, 10);
	FourNodeMVLEM3DKlocal(10, 10) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(5.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 11) = 0.0;
	FourNodeMVLEM3DKlocal(10, 12) = 0.0;
	FourNodeMVLEM3DKlocal(10, 13) = 0.0;
	FourNodeMVLEM3DKlocal(10, 14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 15) = 0.0;
	FourNodeMVLEM3DKlocal(10, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 17) = 0.0;
	FourNodeMVLEM3DKlocal(10, 18) = 0.0;
	FourNodeMVLEM3DKlocal(10, 19) = 0.0;
	FourNodeMVLEM3DKlocal(10, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 21) = 0.0;
	FourNodeMVLEM3DKlocal(10, 22) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(10, 23) = 0.0;

	FourNodeMVLEM3DKlocal(11, 0) = FourNodeMVLEM3DKlocal(0, 11);
	FourNodeMVLEM3DKlocal(11, 1) = FourNodeMVLEM3DKlocal(1, 11);
	FourNodeMVLEM3DKlocal(11, 2) = FourNodeMVLEM3DKlocal(2, 11);
	FourNodeMVLEM3DKlocal(11, 3) = FourNodeMVLEM3DKlocal(3, 11);
	FourNodeMVLEM3DKlocal(11, 4) = FourNodeMVLEM3DKlocal(4, 11);
	FourNodeMVLEM3DKlocal(11, 5) = FourNodeMVLEM3DKlocal(5, 11);
	FourNodeMVLEM3DKlocal(11, 6) = FourNodeMVLEM3DKlocal(6, 11);
	FourNodeMVLEM3DKlocal(11, 7) = FourNodeMVLEM3DKlocal(7, 11);
	FourNodeMVLEM3DKlocal(11, 8) = FourNodeMVLEM3DKlocal(8, 11);
	FourNodeMVLEM3DKlocal(11, 9) = FourNodeMVLEM3DKlocal(9, 11);
	FourNodeMVLEM3DKlocal(11, 10) = FourNodeMVLEM3DKlocal(10, 11);
	FourNodeMVLEM3DKlocal(11, 11) = (Km + Kh*(c*c)*(h*h)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw;
	FourNodeMVLEM3DKlocal(11, 12) = (Kh*c*h) / (4 * (d*d) + 4);
	FourNodeMVLEM3DKlocal(11, 13) = (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) - e / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(11, 14) = 0.0;
	FourNodeMVLEM3DKlocal(11, 15) = 0.0;
	FourNodeMVLEM3DKlocal(11, 16) = 0.0;
	FourNodeMVLEM3DKlocal(11, 17) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(11, 18) = (Kh*c*h) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(11, 19) = -e / (4.0 * (d*d) + 4.0) - (d*(Km + Kh*c*(h*h)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(11, 20) = 0.0;
	FourNodeMVLEM3DKlocal(11, 21) = 0.0;
	FourNodeMVLEM3DKlocal(11, 22) = 0.0;
	FourNodeMVLEM3DKlocal(11, 23) = -(Km + Kh*c*(h*h)*(c - 1.0)) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));

	FourNodeMVLEM3DKlocal(12, 0) = FourNodeMVLEM3DKlocal(0, 12);
	FourNodeMVLEM3DKlocal(12, 1) = FourNodeMVLEM3DKlocal(1, 12);
	FourNodeMVLEM3DKlocal(12, 2) = FourNodeMVLEM3DKlocal(2, 12);
	FourNodeMVLEM3DKlocal(12, 3) = FourNodeMVLEM3DKlocal(3, 12);
	FourNodeMVLEM3DKlocal(12, 4) = FourNodeMVLEM3DKlocal(4, 12);
	FourNodeMVLEM3DKlocal(12, 5) = FourNodeMVLEM3DKlocal(5, 12);
	FourNodeMVLEM3DKlocal(12, 6) = FourNodeMVLEM3DKlocal(6, 12);
	FourNodeMVLEM3DKlocal(12, 7) = FourNodeMVLEM3DKlocal(7, 12);
	FourNodeMVLEM3DKlocal(12, 8) = FourNodeMVLEM3DKlocal(8, 12);
	FourNodeMVLEM3DKlocal(12, 9) = FourNodeMVLEM3DKlocal(9, 12);
	FourNodeMVLEM3DKlocal(12, 10) = FourNodeMVLEM3DKlocal(10, 12);
	FourNodeMVLEM3DKlocal(12, 11) = FourNodeMVLEM3DKlocal(11, 12);
	FourNodeMVLEM3DKlocal(12, 12) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(12, 13) = (Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(12, 14) = 0.0;
	FourNodeMVLEM3DKlocal(12, 15) = 0.0;
	FourNodeMVLEM3DKlocal(12, 16) = 0.0;
	FourNodeMVLEM3DKlocal(12, 17) = -(Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(12, 18) = Kh / 4.0 - (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(12, 19) = -(Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(12, 20) = 0.0;
	FourNodeMVLEM3DKlocal(12, 21) = 0.0;
	FourNodeMVLEM3DKlocal(12, 22) = 0.0;
	FourNodeMVLEM3DKlocal(12, 23) = -(Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);

	FourNodeMVLEM3DKlocal(13, 0) = FourNodeMVLEM3DKlocal(0, 13);
	FourNodeMVLEM3DKlocal(13, 1) = FourNodeMVLEM3DKlocal(1, 13);
	FourNodeMVLEM3DKlocal(13, 2) = FourNodeMVLEM3DKlocal(2, 13);
	FourNodeMVLEM3DKlocal(13, 3) = FourNodeMVLEM3DKlocal(3, 13);
	FourNodeMVLEM3DKlocal(13, 4) = FourNodeMVLEM3DKlocal(4, 13);
	FourNodeMVLEM3DKlocal(13, 5) = FourNodeMVLEM3DKlocal(5, 13);
	FourNodeMVLEM3DKlocal(13, 6) = FourNodeMVLEM3DKlocal(6, 13);
	FourNodeMVLEM3DKlocal(13, 7) = FourNodeMVLEM3DKlocal(7, 13);
	FourNodeMVLEM3DKlocal(13, 8) = FourNodeMVLEM3DKlocal(8, 13);
	FourNodeMVLEM3DKlocal(13, 9) = FourNodeMVLEM3DKlocal(9, 13);
	FourNodeMVLEM3DKlocal(13, 10) = FourNodeMVLEM3DKlocal(10, 13);
	FourNodeMVLEM3DKlocal(13, 11) = FourNodeMVLEM3DKlocal(11, 13);
	FourNodeMVLEM3DKlocal(13, 12) = FourNodeMVLEM3DKlocal(12, 13);
	FourNodeMVLEM3DKlocal(13, 13) = Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw) - (d*(e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(13, 14) = 0.0;
	FourNodeMVLEM3DKlocal(13, 15) = 0.0;
	FourNodeMVLEM3DKlocal(13, 16) = 0.0;
	FourNodeMVLEM3DKlocal(13, 17) = (e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);
	FourNodeMVLEM3DKlocal(13, 18) = (Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(13, 19) = Kv / 4.0 - (d*e) / (4.0 * (d*d) + 4.0) - (12.0 * Eim*Iim) / (Lw*Lw*Lw) + (d*(e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(13, 20) = 0.0;
	FourNodeMVLEM3DKlocal(13, 21) = 0.0;
	FourNodeMVLEM3DKlocal(13, 22) = 0.0;
	FourNodeMVLEM3DKlocal(13, 23) = (e / 2.0 - (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim) / (Lw*Lw);

	FourNodeMVLEM3DKlocal(14, 0) = FourNodeMVLEM3DKlocal(0, 14);
	FourNodeMVLEM3DKlocal(14, 1) = FourNodeMVLEM3DKlocal(1, 14);
	FourNodeMVLEM3DKlocal(14, 2) = FourNodeMVLEM3DKlocal(2, 14);
	FourNodeMVLEM3DKlocal(14, 3) = FourNodeMVLEM3DKlocal(3, 14);
	FourNodeMVLEM3DKlocal(14, 4) = FourNodeMVLEM3DKlocal(4, 14);
	FourNodeMVLEM3DKlocal(14, 5) = FourNodeMVLEM3DKlocal(5, 14);
	FourNodeMVLEM3DKlocal(14, 6) = FourNodeMVLEM3DKlocal(6, 14);
	FourNodeMVLEM3DKlocal(14, 7) = FourNodeMVLEM3DKlocal(7, 14);
	FourNodeMVLEM3DKlocal(14, 8) = FourNodeMVLEM3DKlocal(8, 14);
	FourNodeMVLEM3DKlocal(14, 9) = FourNodeMVLEM3DKlocal(9, 14);
	FourNodeMVLEM3DKlocal(14, 10) = FourNodeMVLEM3DKlocal(10, 14);
	FourNodeMVLEM3DKlocal(14, 11) = FourNodeMVLEM3DKlocal(11, 14);
	FourNodeMVLEM3DKlocal(14, 12) = FourNodeMVLEM3DKlocal(12, 14);
	FourNodeMVLEM3DKlocal(14, 13) = FourNodeMVLEM3DKlocal(13, 14);
	FourNodeMVLEM3DKlocal(14, 14) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 15) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 17) = 0.0;
	FourNodeMVLEM3DKlocal(14, 18) = 0.0;
	FourNodeMVLEM3DKlocal(14, 19) = 0.0;
	FourNodeMVLEM3DKlocal(14, 20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(14, 23) = 0.0;

	FourNodeMVLEM3DKlocal(15, 0) = FourNodeMVLEM3DKlocal(0, 15);
	FourNodeMVLEM3DKlocal(15, 1) = FourNodeMVLEM3DKlocal(1, 15);
	FourNodeMVLEM3DKlocal(15, 2) = FourNodeMVLEM3DKlocal(2, 15);
	FourNodeMVLEM3DKlocal(15, 3) = FourNodeMVLEM3DKlocal(3, 15);
	FourNodeMVLEM3DKlocal(15, 4) = FourNodeMVLEM3DKlocal(4, 15);
	FourNodeMVLEM3DKlocal(15, 5) = FourNodeMVLEM3DKlocal(5, 15);
	FourNodeMVLEM3DKlocal(15, 6) = FourNodeMVLEM3DKlocal(6, 15);
	FourNodeMVLEM3DKlocal(15, 7) = FourNodeMVLEM3DKlocal(7, 15);
	FourNodeMVLEM3DKlocal(15, 8) = FourNodeMVLEM3DKlocal(8, 15);
	FourNodeMVLEM3DKlocal(15, 9) = FourNodeMVLEM3DKlocal(9, 15);
	FourNodeMVLEM3DKlocal(15, 10) = FourNodeMVLEM3DKlocal(10, 15);
	FourNodeMVLEM3DKlocal(15, 11) = FourNodeMVLEM3DKlocal(11, 15);
	FourNodeMVLEM3DKlocal(15, 12) = FourNodeMVLEM3DKlocal(12, 15);
	FourNodeMVLEM3DKlocal(15, 13) = FourNodeMVLEM3DKlocal(13, 15);
	FourNodeMVLEM3DKlocal(15, 14) = FourNodeMVLEM3DKlocal(14, 15);
	FourNodeMVLEM3DKlocal(15, 15) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(15, 16) = -(Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeMVLEM3DKlocal(15, 17) = 0.0;
	FourNodeMVLEM3DKlocal(15, 18) = 0.0;
	FourNodeMVLEM3DKlocal(15, 19) = 0.0;
	FourNodeMVLEM3DKlocal(15, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(15, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(15, 22) = 0.0;
	FourNodeMVLEM3DKlocal(15, 23) = 0.0;

	FourNodeMVLEM3DKlocal(16, 0) = FourNodeMVLEM3DKlocal(0, 16);
	FourNodeMVLEM3DKlocal(16, 1) = FourNodeMVLEM3DKlocal(1, 16);
	FourNodeMVLEM3DKlocal(16, 2) = FourNodeMVLEM3DKlocal(2, 16);
	FourNodeMVLEM3DKlocal(16, 3) = FourNodeMVLEM3DKlocal(3, 16);
	FourNodeMVLEM3DKlocal(16, 4) = FourNodeMVLEM3DKlocal(4, 16);
	FourNodeMVLEM3DKlocal(16, 5) = FourNodeMVLEM3DKlocal(5, 16);
	FourNodeMVLEM3DKlocal(16, 6) = FourNodeMVLEM3DKlocal(6, 16);
	FourNodeMVLEM3DKlocal(16, 7) = FourNodeMVLEM3DKlocal(7, 16);
	FourNodeMVLEM3DKlocal(16, 8) = FourNodeMVLEM3DKlocal(8, 16);
	FourNodeMVLEM3DKlocal(16, 9) = FourNodeMVLEM3DKlocal(9, 16);
	FourNodeMVLEM3DKlocal(16, 10) = FourNodeMVLEM3DKlocal(10, 16);
	FourNodeMVLEM3DKlocal(16, 11) = FourNodeMVLEM3DKlocal(11, 16);
	FourNodeMVLEM3DKlocal(16, 12) = FourNodeMVLEM3DKlocal(12, 16);
	FourNodeMVLEM3DKlocal(16, 13) = FourNodeMVLEM3DKlocal(13, 16);
	FourNodeMVLEM3DKlocal(16, 14) = FourNodeMVLEM3DKlocal(14, 16);
	FourNodeMVLEM3DKlocal(16, 15) = FourNodeMVLEM3DKlocal(15, 16);
	FourNodeMVLEM3DKlocal(16, 16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(16, 17) = 0.0;
	FourNodeMVLEM3DKlocal(16, 18) = 0.0;
	FourNodeMVLEM3DKlocal(16, 19) = 0.0;
	FourNodeMVLEM3DKlocal(16, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(16, 21) = 0.0;
	FourNodeMVLEM3DKlocal(16, 22) = -(Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(16, 23) = 0.0;

	FourNodeMVLEM3DKlocal(17, 0) = FourNodeMVLEM3DKlocal(0, 17);
	FourNodeMVLEM3DKlocal(17, 1) = FourNodeMVLEM3DKlocal(1, 17);
	FourNodeMVLEM3DKlocal(17, 2) = FourNodeMVLEM3DKlocal(2, 17);
	FourNodeMVLEM3DKlocal(17, 3) = FourNodeMVLEM3DKlocal(3, 17);
	FourNodeMVLEM3DKlocal(17, 4) = FourNodeMVLEM3DKlocal(4, 17);
	FourNodeMVLEM3DKlocal(17, 5) = FourNodeMVLEM3DKlocal(5, 17);
	FourNodeMVLEM3DKlocal(17, 6) = FourNodeMVLEM3DKlocal(6, 17);
	FourNodeMVLEM3DKlocal(17, 7) = FourNodeMVLEM3DKlocal(7, 17);
	FourNodeMVLEM3DKlocal(17, 8) = FourNodeMVLEM3DKlocal(8, 17);
	FourNodeMVLEM3DKlocal(17, 9) = FourNodeMVLEM3DKlocal(9, 17);
	FourNodeMVLEM3DKlocal(17, 10) = FourNodeMVLEM3DKlocal(10, 17);
	FourNodeMVLEM3DKlocal(17, 11) = FourNodeMVLEM3DKlocal(11, 17);
	FourNodeMVLEM3DKlocal(17, 12) = FourNodeMVLEM3DKlocal(12, 17);
	FourNodeMVLEM3DKlocal(17, 13) = FourNodeMVLEM3DKlocal(13, 17);
	FourNodeMVLEM3DKlocal(17, 14) = FourNodeMVLEM3DKlocal(14, 17);
	FourNodeMVLEM3DKlocal(17, 15) = FourNodeMVLEM3DKlocal(15, 17);
	FourNodeMVLEM3DKlocal(17, 16) = FourNodeMVLEM3DKlocal(16, 17);
	FourNodeMVLEM3DKlocal(17, 17) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw;
	FourNodeMVLEM3DKlocal(17, 18) = -(Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(17, 19) = e / (4.0 * (d*d) + 4.0) - (6.0 * Eim*Iim) / (Lw*Lw) + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0));
	FourNodeMVLEM3DKlocal(17, 20) = 0.0;
	FourNodeMVLEM3DKlocal(17, 21) = 0.0;
	FourNodeMVLEM3DKlocal(17, 22) = 0.0;
	FourNodeMVLEM3DKlocal(17, 23) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (2.0 * Eim*Iim) / Lw;

	FourNodeMVLEM3DKlocal(18, 0) = FourNodeMVLEM3DKlocal(0, 18);
	FourNodeMVLEM3DKlocal(18, 1) = FourNodeMVLEM3DKlocal(1, 18);
	FourNodeMVLEM3DKlocal(18, 2) = FourNodeMVLEM3DKlocal(2, 18);
	FourNodeMVLEM3DKlocal(18, 3) = FourNodeMVLEM3DKlocal(3, 18);
	FourNodeMVLEM3DKlocal(18, 4) = FourNodeMVLEM3DKlocal(4, 18);
	FourNodeMVLEM3DKlocal(18, 5) = FourNodeMVLEM3DKlocal(5, 18);
	FourNodeMVLEM3DKlocal(18, 6) = FourNodeMVLEM3DKlocal(6, 18);
	FourNodeMVLEM3DKlocal(18, 7) = FourNodeMVLEM3DKlocal(7, 18);
	FourNodeMVLEM3DKlocal(18, 8) = FourNodeMVLEM3DKlocal(8, 18);
	FourNodeMVLEM3DKlocal(18, 9) = FourNodeMVLEM3DKlocal(9, 18);
	FourNodeMVLEM3DKlocal(18, 10) = FourNodeMVLEM3DKlocal(10, 18);
	FourNodeMVLEM3DKlocal(18, 11) = FourNodeMVLEM3DKlocal(11, 18);
	FourNodeMVLEM3DKlocal(18, 12) = FourNodeMVLEM3DKlocal(12, 18);
	FourNodeMVLEM3DKlocal(18, 13) = FourNodeMVLEM3DKlocal(13, 18);
	FourNodeMVLEM3DKlocal(18, 14) = FourNodeMVLEM3DKlocal(14, 18);
	FourNodeMVLEM3DKlocal(18, 15) = FourNodeMVLEM3DKlocal(15, 18);
	FourNodeMVLEM3DKlocal(18, 16) = FourNodeMVLEM3DKlocal(16, 18);
	FourNodeMVLEM3DKlocal(18, 17) = FourNodeMVLEM3DKlocal(17, 18);
	FourNodeMVLEM3DKlocal(18, 18) = Kh / 4.0 + (Aim*Eim) / Lw;
	FourNodeMVLEM3DKlocal(18, 19) = -(Kh*d*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);
	FourNodeMVLEM3DKlocal(18, 20) = 0.0;
	FourNodeMVLEM3DKlocal(18, 21) = 0.0;
	FourNodeMVLEM3DKlocal(18, 22) = 0.0;
	FourNodeMVLEM3DKlocal(18, 23) = -(Kh*h*(c - 1.0)) / (4.0 * (d*d) + 4.0);

	FourNodeMVLEM3DKlocal(19, 0) = FourNodeMVLEM3DKlocal(0, 19);
	FourNodeMVLEM3DKlocal(19, 1) = FourNodeMVLEM3DKlocal(1, 19);
	FourNodeMVLEM3DKlocal(19, 2) = FourNodeMVLEM3DKlocal(2, 19);
	FourNodeMVLEM3DKlocal(19, 3) = FourNodeMVLEM3DKlocal(3, 19);
	FourNodeMVLEM3DKlocal(19, 4) = FourNodeMVLEM3DKlocal(4, 19);
	FourNodeMVLEM3DKlocal(19, 5) = FourNodeMVLEM3DKlocal(5, 19);
	FourNodeMVLEM3DKlocal(19, 6) = FourNodeMVLEM3DKlocal(6, 19);
	FourNodeMVLEM3DKlocal(19, 7) = FourNodeMVLEM3DKlocal(7, 19);
	FourNodeMVLEM3DKlocal(19, 8) = FourNodeMVLEM3DKlocal(8, 19);
	FourNodeMVLEM3DKlocal(19, 9) = FourNodeMVLEM3DKlocal(9, 19);
	FourNodeMVLEM3DKlocal(19, 10) = FourNodeMVLEM3DKlocal(10, 19);
	FourNodeMVLEM3DKlocal(19, 11) = FourNodeMVLEM3DKlocal(11, 19);
	FourNodeMVLEM3DKlocal(19, 12) = FourNodeMVLEM3DKlocal(12, 19);
	FourNodeMVLEM3DKlocal(19, 13) = FourNodeMVLEM3DKlocal(13, 19);
	FourNodeMVLEM3DKlocal(19, 14) = FourNodeMVLEM3DKlocal(14, 19);
	FourNodeMVLEM3DKlocal(19, 15) = FourNodeMVLEM3DKlocal(15, 19);
	FourNodeMVLEM3DKlocal(19, 16) = FourNodeMVLEM3DKlocal(16, 19);
	FourNodeMVLEM3DKlocal(19, 17) = FourNodeMVLEM3DKlocal(17, 19);
	FourNodeMVLEM3DKlocal(19, 18) = FourNodeMVLEM3DKlocal(18, 19);
	FourNodeMVLEM3DKlocal(19, 19) = Kv / 4.0 + (d*e) / (4.0 * (d*d) + 4.0) + (12.0 * Eim*Iim) / (Lw*Lw*Lw) + (d*(e / 2.0 + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0))) / (2.0 * (d*d) + 2.0);
	FourNodeMVLEM3DKlocal(19, 20) = 0.0;
	FourNodeMVLEM3DKlocal(19, 21) = 0.0;
	FourNodeMVLEM3DKlocal(19, 22) = 0.0;
	FourNodeMVLEM3DKlocal(19, 23) = (e / 2.0 + (d*(Km + Kh*(h*h)*((c - 1.0)*(c - 1.0)))) / (2.0 * (d*d) + 2.0)) / (2.0 * (d*d) + 2.0) - (6.0 * Eim*Iim) / (Lw*Lw);

	FourNodeMVLEM3DKlocal(20, 0) = FourNodeMVLEM3DKlocal(0, 20);
	FourNodeMVLEM3DKlocal(20, 1) = FourNodeMVLEM3DKlocal(1, 20);
	FourNodeMVLEM3DKlocal(20, 2) = FourNodeMVLEM3DKlocal(2, 20);
	FourNodeMVLEM3DKlocal(20, 3) = FourNodeMVLEM3DKlocal(3, 20);
	FourNodeMVLEM3DKlocal(20, 4) = FourNodeMVLEM3DKlocal(4, 20);
	FourNodeMVLEM3DKlocal(20, 5) = FourNodeMVLEM3DKlocal(5, 20);
	FourNodeMVLEM3DKlocal(20, 6) = FourNodeMVLEM3DKlocal(6, 20);
	FourNodeMVLEM3DKlocal(20, 7) = FourNodeMVLEM3DKlocal(7, 20);
	FourNodeMVLEM3DKlocal(20, 8) = FourNodeMVLEM3DKlocal(8, 20);
	FourNodeMVLEM3DKlocal(20, 9) = FourNodeMVLEM3DKlocal(9, 20);
	FourNodeMVLEM3DKlocal(20, 10) = FourNodeMVLEM3DKlocal(10, 20);
	FourNodeMVLEM3DKlocal(20, 11) = FourNodeMVLEM3DKlocal(11, 20);
	FourNodeMVLEM3DKlocal(20, 12) = FourNodeMVLEM3DKlocal(12, 20);
	FourNodeMVLEM3DKlocal(20, 13) = FourNodeMVLEM3DKlocal(13, 20);
	FourNodeMVLEM3DKlocal(20, 14) = FourNodeMVLEM3DKlocal(14, 20);
	FourNodeMVLEM3DKlocal(20, 15) = FourNodeMVLEM3DKlocal(15, 20);
	FourNodeMVLEM3DKlocal(20, 16) = FourNodeMVLEM3DKlocal(16, 20);
	FourNodeMVLEM3DKlocal(20, 17) = FourNodeMVLEM3DKlocal(17, 20);
	FourNodeMVLEM3DKlocal(20, 18) = FourNodeMVLEM3DKlocal(18, 20);
	FourNodeMVLEM3DKlocal(20, 19) = FourNodeMVLEM3DKlocal(19, 20);
	FourNodeMVLEM3DKlocal(20, 20) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(20, 21) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(20, 22) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(20, 23) = 0.0;

	FourNodeMVLEM3DKlocal(21, 0) = FourNodeMVLEM3DKlocal(0, 21);
	FourNodeMVLEM3DKlocal(21, 1) = FourNodeMVLEM3DKlocal(1, 21);
	FourNodeMVLEM3DKlocal(21, 2) = FourNodeMVLEM3DKlocal(2, 21);
	FourNodeMVLEM3DKlocal(21, 3) = FourNodeMVLEM3DKlocal(3, 21);
	FourNodeMVLEM3DKlocal(21, 4) = FourNodeMVLEM3DKlocal(4, 21);
	FourNodeMVLEM3DKlocal(21, 5) = FourNodeMVLEM3DKlocal(5, 21);
	FourNodeMVLEM3DKlocal(21, 6) = FourNodeMVLEM3DKlocal(6, 21);
	FourNodeMVLEM3DKlocal(21, 7) = FourNodeMVLEM3DKlocal(7, 21);
	FourNodeMVLEM3DKlocal(21, 8) = FourNodeMVLEM3DKlocal(8, 21);
	FourNodeMVLEM3DKlocal(21, 9) = FourNodeMVLEM3DKlocal(9, 21);
	FourNodeMVLEM3DKlocal(21, 10) = FourNodeMVLEM3DKlocal(10, 21);
	FourNodeMVLEM3DKlocal(21, 11) = FourNodeMVLEM3DKlocal(11, 21);
	FourNodeMVLEM3DKlocal(21, 12) = FourNodeMVLEM3DKlocal(12, 21);
	FourNodeMVLEM3DKlocal(21, 13) = FourNodeMVLEM3DKlocal(13, 21);
	FourNodeMVLEM3DKlocal(21, 14) = FourNodeMVLEM3DKlocal(14, 21);
	FourNodeMVLEM3DKlocal(21, 15) = FourNodeMVLEM3DKlocal(15, 21);
	FourNodeMVLEM3DKlocal(21, 16) = FourNodeMVLEM3DKlocal(16, 21);
	FourNodeMVLEM3DKlocal(21, 17) = FourNodeMVLEM3DKlocal(17, 21);
	FourNodeMVLEM3DKlocal(21, 18) = FourNodeMVLEM3DKlocal(18, 21);
	FourNodeMVLEM3DKlocal(21, 19) = FourNodeMVLEM3DKlocal(19, 21);
	FourNodeMVLEM3DKlocal(21, 20) = FourNodeMVLEM3DKlocal(20, 21);
	FourNodeMVLEM3DKlocal(21, 21) = -(Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(21, 22) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (12.0 * (NUelastic*NUelastic) - 12.0);
	FourNodeMVLEM3DKlocal(21, 23) = 0.0;

	FourNodeMVLEM3DKlocal(22, 0) = FourNodeMVLEM3DKlocal(0, 22);
	FourNodeMVLEM3DKlocal(22, 1) = FourNodeMVLEM3DKlocal(1, 22);
	FourNodeMVLEM3DKlocal(22, 2) = FourNodeMVLEM3DKlocal(2, 22);
	FourNodeMVLEM3DKlocal(22, 3) = FourNodeMVLEM3DKlocal(3, 22);
	FourNodeMVLEM3DKlocal(22, 4) = FourNodeMVLEM3DKlocal(4, 22);
	FourNodeMVLEM3DKlocal(22, 5) = FourNodeMVLEM3DKlocal(5, 22);
	FourNodeMVLEM3DKlocal(22, 6) = FourNodeMVLEM3DKlocal(6, 22);
	FourNodeMVLEM3DKlocal(22, 7) = FourNodeMVLEM3DKlocal(7, 22);
	FourNodeMVLEM3DKlocal(22, 8) = FourNodeMVLEM3DKlocal(8, 22);
	FourNodeMVLEM3DKlocal(22, 9) = FourNodeMVLEM3DKlocal(9, 22);
	FourNodeMVLEM3DKlocal(22, 10) = FourNodeMVLEM3DKlocal(10, 22);
	FourNodeMVLEM3DKlocal(22, 11) = FourNodeMVLEM3DKlocal(11, 22);
	FourNodeMVLEM3DKlocal(22, 12) = FourNodeMVLEM3DKlocal(12, 22);
	FourNodeMVLEM3DKlocal(22, 13) = FourNodeMVLEM3DKlocal(13, 22);
	FourNodeMVLEM3DKlocal(22, 14) = FourNodeMVLEM3DKlocal(14, 22);
	FourNodeMVLEM3DKlocal(22, 15) = FourNodeMVLEM3DKlocal(15, 22);
	FourNodeMVLEM3DKlocal(22, 16) = FourNodeMVLEM3DKlocal(16, 22);
	FourNodeMVLEM3DKlocal(22, 17) = FourNodeMVLEM3DKlocal(17, 22);
	FourNodeMVLEM3DKlocal(22, 18) = FourNodeMVLEM3DKlocal(18, 22);
	FourNodeMVLEM3DKlocal(22, 19) = FourNodeMVLEM3DKlocal(19, 22);
	FourNodeMVLEM3DKlocal(22, 20) = FourNodeMVLEM3DKlocal(20, 22);
	FourNodeMVLEM3DKlocal(22, 21) = FourNodeMVLEM3DKlocal(21, 22);
	FourNodeMVLEM3DKlocal(22, 22) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0)) - (Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DKlocal(22, 23) = 0.0;

	FourNodeMVLEM3DKlocal(23, 0) = FourNodeMVLEM3DKlocal(0, 23);
	FourNodeMVLEM3DKlocal(23, 1) = FourNodeMVLEM3DKlocal(1, 23);
	FourNodeMVLEM3DKlocal(23, 2) = FourNodeMVLEM3DKlocal(2, 23);
	FourNodeMVLEM3DKlocal(23, 3) = FourNodeMVLEM3DKlocal(3, 23);
	FourNodeMVLEM3DKlocal(23, 4) = FourNodeMVLEM3DKlocal(4, 23);
	FourNodeMVLEM3DKlocal(23, 5) = FourNodeMVLEM3DKlocal(5, 23);
	FourNodeMVLEM3DKlocal(23, 6) = FourNodeMVLEM3DKlocal(6, 23);
	FourNodeMVLEM3DKlocal(23, 7) = FourNodeMVLEM3DKlocal(7, 23);
	FourNodeMVLEM3DKlocal(23, 8) = FourNodeMVLEM3DKlocal(8, 23);
	FourNodeMVLEM3DKlocal(23, 9) = FourNodeMVLEM3DKlocal(9, 23);
	FourNodeMVLEM3DKlocal(23, 10) = FourNodeMVLEM3DKlocal(10, 23);
	FourNodeMVLEM3DKlocal(23, 11) = FourNodeMVLEM3DKlocal(11, 23);
	FourNodeMVLEM3DKlocal(23, 12) = FourNodeMVLEM3DKlocal(12, 23);
	FourNodeMVLEM3DKlocal(23, 13) = FourNodeMVLEM3DKlocal(13, 23);
	FourNodeMVLEM3DKlocal(23, 14) = FourNodeMVLEM3DKlocal(14, 23);
	FourNodeMVLEM3DKlocal(23, 15) = FourNodeMVLEM3DKlocal(15, 23);
	FourNodeMVLEM3DKlocal(23, 16) = FourNodeMVLEM3DKlocal(16, 23);
	FourNodeMVLEM3DKlocal(23, 17) = FourNodeMVLEM3DKlocal(17, 23);
	FourNodeMVLEM3DKlocal(23, 18) = FourNodeMVLEM3DKlocal(18, 23);
	FourNodeMVLEM3DKlocal(23, 19) = FourNodeMVLEM3DKlocal(19, 23);
	FourNodeMVLEM3DKlocal(23, 20) = FourNodeMVLEM3DKlocal(20, 23);
	FourNodeMVLEM3DKlocal(23, 21) = FourNodeMVLEM3DKlocal(21, 23);
	FourNodeMVLEM3DKlocal(23, 22) = FourNodeMVLEM3DKlocal(22, 23);
	FourNodeMVLEM3DKlocal(23, 23) = (Km + Kh*(h*h)*((c - 1.0)*(c - 1.0))) / ((2.0 * (d*d) + 2.0)*(2.0 * (d*d) + 2.0)) + (4.0 * Eim*Iim) / Lw;

	// Convert matrix from local to global cs
	FourNodeMVLEM3DK.addMatrixTripleProduct(0.0, T, FourNodeMVLEM3DKlocal, 1.0); 

	// Return element stiffness matrix
	return FourNodeMVLEM3DK;
}

// Get element mass matrix assuming lumped mass
const Matrix & FourNodeMVLEM3D::getMass(void)
{

	// !!! Zeroing matrices at beginning versus not zeroing at the beginning?!

	// No rotational mass
	FourNodeMVLEM3DMlocal(0, 0) = NodeMass;
	FourNodeMVLEM3DMlocal(1, 1) = NodeMass;
	FourNodeMVLEM3DMlocal(2, 2) = NodeMass;

	FourNodeMVLEM3DMlocal(6, 6) = NodeMass;
	FourNodeMVLEM3DMlocal(7, 7) = NodeMass;
	FourNodeMVLEM3DMlocal(8, 8) = NodeMass;

	FourNodeMVLEM3DMlocal(12, 12) = NodeMass;
	FourNodeMVLEM3DMlocal(13, 13) = NodeMass;
	FourNodeMVLEM3DMlocal(14, 14) = NodeMass;

	FourNodeMVLEM3DMlocal(18, 18) = NodeMass;
	FourNodeMVLEM3DMlocal(19, 19) = NodeMass;
	FourNodeMVLEM3DMlocal(20, 20) = NodeMass;

	// Convert matrix from local to global cs !!! why is mass defined in local coordinate system? I think it should be global.
	FourNodeMVLEM3DM.addMatrixTripleProduct(0.0, T, FourNodeMVLEM3DMlocal, 1.0); 

	// Return element mass matrix
	return FourNodeMVLEM3DM;
}

// Get element damping matrix
const Matrix & FourNodeMVLEM3D::getDamp(void)
{
	FourNodeMVLEM3DD.Zero();

	FourNodeMVLEM3DD = this->Element::getDamp();

	// Return element damping matrix
	return FourNodeMVLEM3DD;
}

// zeroLoad
void FourNodeMVLEM3D::zeroLoad(void)
{
	// does nothing - no elemental loads
}

// addLoad
int FourNodeMVLEM3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	return 0;
}

int FourNodeMVLEM3D::addInertiaLoadToUnbalance(const Vector &accel)
{
	return 0;
}

// Get element force vector
const Vector & FourNodeMVLEM3D::getResistingForce()
{

	// Get Trial Displacements
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();

	// Create force vectors and assign zero
	Vector dispG(24); // global cs
	Vector dispL(24); // local cs

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

	// Get force from shear force-deformation relationship
	R1 = theMaterialsShear[0]->getStress(); 

	// Get stresses from uniaxial material models for each fiber
	for (int i = 0; i < m; ++i) {
		stressC[i] = theMaterialsConcrete[i]->getStress();
		stressS[i] = theMaterialsSteel[i]->getStress();
	}

	for (int i = 0; i<m; i++) {
		R2 += -stressC[i] * Ac[i] - stressS[i] * As[i];
		R3 += -stressC[i] * Ac[i] * x[i] - stressS[i] * As[i] * x[i];
		R6 += stressC[i] * Ac[i] * x[i] + stressS[i] * As[i] * x[i];
	}

	R3 += -R1*c*h;
	R4 = -R1;
	R5 = -R2;
	R6 += -R1*(1.0 - c)*h;

	// Calculate force vector in local cs - check Matlab derivations !!!
	FourNodeMVLEM3DRlocal(0) = R1 / 2.0 + (Aim*Eim*dispL(0)) / Lw - (Aim*Eim*dispL(6)) / Lw;
	FourNodeMVLEM3DRlocal(1) = R2 / 2.0 - (R3*d) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim*dispL(1)) / (Lw*Lw*Lw) + (6.0 * Eim*Iim*dispL(5)) / (Lw*Lw) - (12.0 * Eim*Iim*dispL(7)) / (Lw*Lw*Lw) + (6.0 * Eim*Iim*dispL(11)) / (Lw*Lw);
	FourNodeMVLEM3DRlocal(2) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(3) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h) - (h*h)*NUelastic + 5 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(4) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (dispL(4))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(10))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(22))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(16))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(5) = R3 / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim*dispL(1)) / (Lw*Lw) + (4.0 * Eim*Iim*dispL(5)) / Lw - (6.0 * Eim*Iim*dispL(7)) / (Lw*Lw) + (2.0 * Eim*Iim*dispL(11)) / Lw;
	FourNodeMVLEM3DRlocal(6) = R1 / 2.0 - (Aim*Eim*dispL(0)) / Lw + (Aim*Eim*dispL(6)) / Lw;
	FourNodeMVLEM3DRlocal(7) = R2 / 2.0 + (R3*d) / (2.0 * (d*d) + 2.0) - (12.0 * Eim*Iim*dispL(1)) / (Lw*Lw*Lw) - (6.0 * Eim*Iim*dispL(5)) / (Lw*Lw) + (12.0 * Eim*Iim*dispL(7)) / (Lw*Lw*Lw) - (6.0 * Eim*Iim*dispL(11)) / (Lw*Lw);
	FourNodeMVLEM3DRlocal(8) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*(4.0 * (h*h)*NUelastic + (h*h) + 10 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(9) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(10) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (dispL(4))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(16))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (dispL(22))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(5.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(11) = R3 / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim*dispL(1)) / (Lw*Lw) + (2.0 * Eim*Iim*dispL(5)) / Lw - (6.0 * Eim*Iim*dispL(7)) / (Lw*Lw) + (4.0 * Eim*Iim*dispL(11)) / Lw;
	FourNodeMVLEM3DRlocal(12) = R4 / 2.0 + (Aim*Eim*dispL(12)) / Lw - (Aim*Eim*dispL(18)) / Lw;
	FourNodeMVLEM3DRlocal(13) = R5 / 2.0 - (R6*d) / (2.0 * (d*d) + 2.0) + (12.0 * Eim*Iim*dispL(13)) / (Lw*Lw*Lw) + (6.0 * Eim*Iim*dispL(17)) / (Lw*Lw) - (12.0 * Eim*Iim*dispL(19)) / (Lw*Lw*Lw) + (6.0 * Eim*Iim*dispL(23)) / (Lw*Lw);
	FourNodeMVLEM3DRlocal(14) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(15) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(16) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (dispL(16))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(22))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(10))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (dispL(4))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(17) = R6 / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim*dispL(13)) / (Lw*Lw) + (4.0 * Eim*Iim*dispL(17)) / Lw - (6.0 * Eim*Iim*dispL(19)) / (Lw*Lw) + (2.0 * Eim*Iim*dispL(23)) / Lw;
	FourNodeMVLEM3DRlocal(18) = R4 / 2.0 - (Aim*Eim*dispL(12)) / Lw + (Aim*Eim*dispL(18)) / Lw;
	FourNodeMVLEM3DRlocal(19) = R5 / 2.0 + (R6*d) / (2.0 * (d*d) + 2.0) - (12.0 * Eim*Iim*dispL(13)) / (Lw*Lw*Lw) - (6.0 * Eim*Iim*dispL(17)) / (Lw*Lw) + (12.0 * Eim*Iim*dispL(19)) / (Lw*Lw*Lw) - (6.0 * Eim*Iim*dispL(23)) / (Lw*Lw);
	FourNodeMVLEM3DRlocal(20) = (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(10))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(16))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(4))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(21) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(22))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(3))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(9))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(15))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(22) = (Eave*NUelastic*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(21))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (dispL(22))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(16))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(4))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(10))*((Eave*h*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(8))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(14))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(20))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave*Tfactor*Tfactor*Tfactor)*(dispL(2))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	FourNodeMVLEM3DRlocal(23) = R6 / (2.0 * (d*d) + 2.0) + (6.0 * Eim*Iim*dispL(13)) / (Lw*Lw) + (2.0 * Eim*Iim*dispL(17)) / Lw - (6.0 * Eim*Iim*dispL(19)) / (Lw*Lw) + (4.0 * Eim*Iim*dispL(23)) / Lw;

	// Convert force vector from local to global cs
	FourNodeMVLEM3DR.addMatrixTransposeVector(0.0, T, FourNodeMVLEM3DRlocal, 1.0);

	// Return element force vector
	return FourNodeMVLEM3DR;
}

// getResistingForceIncInertia
const Vector & FourNodeMVLEM3D::getResistingForceIncInertia()
{
	this->getResistingForce();

	if (NodeMass != 0.0) {
		const Vector &accel1 = theNodes[0]->getTrialAccel();
		const Vector &accel2 = theNodes[1]->getTrialAccel();
		const Vector &accel3 = theNodes[2]->getTrialAccel();
		const Vector &accel4 = theNodes[3]->getTrialAccel();

		Vector accelG(24);
		accelG.Zero();
		Vector accelL(24);
		accelL.Zero();

		for (int i = 0; i < 6; i++) {
			accelG(i) = accel1(i);
			accelG(i + 6) = accel2(i);
			accelG(i + 12) = accel3(i);
			accelG(i + 18) = accel4(i);
		}

		accelL.addMatrixVector(0.0, T, accelG, 1.0);

		// Compute the current resisting force
		this->getResistingForce();

		// !!! NEEDS TO BE CHECKED - Why local CS? 
		FourNodeMVLEM3DRlocal(0) += NodeMass / 4.0*accelL(0) + NodeMass / 4.0*accelL(6);
		FourNodeMVLEM3DRlocal(1) += NodeMass / 4.0*accelL(1) + NodeMass / 4.0*accelL(7);
		//FourNodeMVLEM3DRlocal(2) += NodeMass / 4.0*accelL(2) + NodeMass / 4.0*accelL(8);
		FourNodeMVLEM3DRlocal(6) += NodeMass / 4.0*accelL(0) + NodeMass / 4.0*accelL(6);
		FourNodeMVLEM3DRlocal(7) += NodeMass / 4.0*accelL(1) + NodeMass / 4.0*accelL(7);
		//FourNodeMVLEM3DRlocal(8) += NodeMass / 4.0*accelL(2) + NodeMass / 4.0*accelL(8);
		FourNodeMVLEM3DRlocal(12) += NodeMass / 4.0*accelL(12) + NodeMass / 4.0*accelL(18);
		FourNodeMVLEM3DRlocal(13) += NodeMass / 4.0*accelL(13) + NodeMass / 4.0*accelL(19);
		//FourNodeMVLEM3DRlocal(14) += NodeMass / 4.0*accelL(14) + NodeMass / 4.0*accelL(20);
		FourNodeMVLEM3DRlocal(18) += NodeMass / 4.0*accelL(12) + NodeMass / 4.0*accelL(18);
		FourNodeMVLEM3DRlocal(19) += NodeMass / 4.0*accelL(13) + NodeMass / 4.0*accelL(19);
		//FourNodeMVLEM3DRlocal(20) += NodeMass / 4.0*accelL(14) + NodeMass / 4.0*accelL(20);

		FourNodeMVLEM3DR.addMatrixTransposeVector(1.0, T, FourNodeMVLEM3DRlocal, 1.0);

		// Add the damping forces if rayleigh damping
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			FourNodeMVLEM3DR += this->getRayleighDampingForces();

	}
	else {

		// Add the damping forces if rayleigh damping
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			FourNodeMVLEM3DR += this->getRayleighDampingForces();
	}

	return FourNodeMVLEM3DR;
}

// sendSelf 
int FourNodeMVLEM3D::sendSelf(int commitTag, Channel &theChannel)
{
	int res;
	int dataTag = this->getDbTag();

	Vector data(6);

	data(0) = this->getTag();
	data(1) = density;
	data(2) = m;
	data(3) = c;
	data(4) = NUelastic;
	data(5) = Tfactor;

	// FourNodeMVLEM3D then sends the tags of it's four end nodes
	res = theChannel.sendID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING FourNodeMVLEM3D::sendSelf() - failed to send ID\n";
		return -2;
	}

	// Send the material class tags
	ID matClassTags(2 * m + 1);
	for (int i = 0; i < m; i++) {
		matClassTags(i) = theMaterialsConcrete[i]->getClassTag();
		matClassTags(i + m) = theMaterialsSteel[i]->getClassTag();
	}

	matClassTags(2 * m) = theMaterialsShear[0]->getClassTag();
	res = theChannel.sendID(0, commitTag, matClassTags);

	// Send the material models
	for (int i = 0; i < m; i++) {
		theMaterialsConcrete[i]->sendSelf(commitTag, theChannel);
		theMaterialsSteel[i]->sendSelf(commitTag, theChannel);
	}
	theMaterialsShear[0]->sendSelf(commitTag, theChannel);

	return 0;

}

// recvSelf
int FourNodeMVLEM3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res;
	int dataTag = this->getDbTag();

	// FourNodeMVLEM3D creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	// delete dynamic memory
	if (theMaterialsConcrete != 0) {
		for (int i = 0; i < m; i++)
			if (theMaterialsConcrete[i] != 0)
				delete theMaterialsConcrete[i];
		delete[] theMaterialsConcrete;
	}

	if (theMaterialsSteel != 0) {
		for (int i = 0; i < m; i++)
			if (theMaterialsSteel[i] != 0)
				delete theMaterialsSteel[i];
		delete[] theMaterialsSteel;
	}

	if (theMaterialsShear != 0) {
		for (int i = 0; i < 1; i++)
			if (theMaterialsShear[i] != 0)
				delete theMaterialsShear[i];
		delete[] theMaterialsShear;
	}

	Vector data(6);
	res = theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING FourNodeMVLEM3D::recvSelf() - failed to receive Vector\n";
		return -1;
	}

	this->setTag((int)data(0));
	density = data(1);
	m = data(2);
	c = data(3);
	NUelastic = data(4);
	Tfactor = data(5);

	// FourNodeMVLEM3D now receives the tags of it's four external nodes
	res = theChannel.recvID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING FourNodeMVLEM3D::recvSelf() - failed to receive ID\n";
		return -2;
	}

	// Receive the material class tags
	ID matClassTags(2 * m + 1);
	res = theChannel.recvID(0, commitTag, matClassTags);

	// Allocate memory for the Concrete uniaxial materials
	theMaterialsConcrete = new UniaxialMaterial*[m];
	if (theMaterialsConcrete == 0) {
		opserr << "FourNodeMVLEM3D::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Allocate memory for the Steel uniaxial materials
	theMaterialsSteel = new UniaxialMaterial*[m];
	if (theMaterialsSteel == 0) {
		opserr << "FourNodeMVLEM3D::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Allocate memory for the Shear uniaxial material
	theMaterialsShear = new UniaxialMaterial*[1];
	if (theMaterialsShear == 0) {
		opserr << "FourNodeMVLEM3D::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Receive the Concrete material models
	for (int i = 0; i < m; i++) {
		theMaterialsConcrete[i] = theBroker.getNewUniaxialMaterial(matClassTags(i));
		if (theMaterialsConcrete[i] == 0) {
			opserr << "FourNodeMVLEM3D::recvSelf() - "
				<< "failed to get blank uniaxial material.\n";
			return -3;
		}
		theMaterialsConcrete[i]->recvSelf(commitTag, theChannel, theBroker);
	}

	// Receive the Steel material models
	for (int i = 0; i < m; i++) {
		theMaterialsSteel[i] = theBroker.getNewUniaxialMaterial(matClassTags(i + m));
		if (theMaterialsSteel[i] == 0) {
			opserr << "FourNodeMVLEM3D::recvSelf() - "
				<< "failed to get blank uniaxial material.\n";
			return -3;
		}
		theMaterialsSteel[i]->recvSelf(commitTag, theChannel, theBroker);
	}

	// Receive the Shear material model
	theMaterialsShear[0] = theBroker.getNewUniaxialMaterial(matClassTags(2 * m));
	if (theMaterialsShear[0] == 0) {
		opserr << "FourNodeMVLEM3D::recvSelf() - "
			<< "failed to get blank uniaxial material.\n";
		return -3;
	}
	theMaterialsShear[0]->recvSelf(commitTag, theChannel, theBroker);

	return 0;
}

// Display model
int FourNodeMVLEM3D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	// First get the end points of the beam based on
	// the display factor (a measure of the distorted image)
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

		int mode = displayMode  *  -1;

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

	// Calculate x, y and z coordinates of V1_1 and v1_2 based on x,y,z coordinates of v1, v2, v3, v4 (take the average)
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
	Lv1_.addMatrixVector(1.0, Tt, Gv1_, 1.0);
	Lv2_.addMatrixVector(1.0, Tt, Gv2_, 1.0);

	// Displaying Fibers
	for (int panel = 0; panel < m; panel++) // loop over m panels
	{

		Matrix NodePLotCrds(m, 13); // (panel id, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)

		// First set the quantity to be displayed at the nodes;
		// if displayMode is 1 through 3 we will plot material stresses otherwise 0.0

		static Vector values(1); // values of epsX to be plotted 

		values(0) = 0.0;

		if (displayMode < 4 && displayMode > 0) {

			values(0) = theMaterialsConcrete[panel]->getStrain();
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

		end1Disp4.addMatrixVector(1.0, T6, end1Disp4G, 1.0);
		end2Disp4.addMatrixVector(1.0, T6, end2Disp4G, 1.0);
		end3Disp4.addMatrixVector(1.0, T6, end3Disp4G, 1.0);
		end4Disp4.addMatrixVector(1.0, T6, end4Disp4G, 1.0);

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
		GlCoord.addMatrixTransposeVector(1.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 1) = GlCoord(0);
		NodePLotCrds(panel, 2) = GlCoord(1);
		NodePLotCrds(panel, 3) = GlCoord(2);
		LocCoord.Zero();
		GlCoord.Zero();
		// Local node 2 - bottom right
		LocCoord(0) = Lv1_(0) + x[panel] + b[panel] / 2.0; // x
		LocCoord(1) = Lv1_(1) + (x[panel] + b[panel] / 2.0)*end1Disp(5)*fact; // y
		LocCoord(2) = Lv1_(2) - (x[panel] + b[panel] / 2.0)*end1Disp(4)*fact; // z
		GlCoord.addMatrixTransposeVector(1.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 4) = GlCoord(0);
		NodePLotCrds(panel, 5) = GlCoord(1);
		NodePLotCrds(panel, 6) = GlCoord(2);
		LocCoord.Zero();
		GlCoord.Zero();
		// Local node 3 - top left
		LocCoord(0) = Lv2_(0) + x[panel] + b[panel] / 2.0; // x
		LocCoord(1) = Lv2_(1) + (x[panel] + b[panel] / 2.0)*end2Disp(5)*fact; // y
		LocCoord(2) = Lv2_(2) - (x[panel] + b[panel] / 2.0)*end2Disp(4)*fact; // z
		GlCoord.addMatrixTransposeVector(1.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 7) = GlCoord(0);
		NodePLotCrds(panel, 8) = GlCoord(1);
		NodePLotCrds(panel, 9) = GlCoord(2);
		LocCoord.Zero();
		GlCoord.Zero();
		// Local node 4 - top rigth
		LocCoord(0) = Lv2_(0) + x[panel] - b[panel] / 2.0; // x
		LocCoord(1) = Lv2_(1) + (x[panel] - b[panel] / 2.0)*end2Disp(5)*fact; // y
		LocCoord(2) = Lv2_(2) - (x[panel] - b[panel] / 2.0)*end2Disp(4)*fact; // z
		GlCoord.addMatrixTransposeVector(1.0, Tt, LocCoord, 1.0);
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

	return error;

}


void FourNodeMVLEM3D::Print(OPS_Stream &s, int flag)
{
	if (flag == 0)
	{
		// Print out element properties
		s << "Element: " << this->getTag() << endln;
		s << "  type: FourNodeMVLEM3D" << endln;
		s << "  iNode: " << externalNodes(0) << ", jNode: " << externalNodes(1) << "  kNode: " << externalNodes(3) << ", lNode: " << externalNodes(2) << endln;
		s << "Element height: " << h << endln;
		s << "Number of uniaxial fibers elements: " << m << endln << endln;

		// determine resisting forces in global system
		s << "  Global resisting force: " << this->getResistingForce() << endln << endln;

		s << "Fiber responses: " << endln;

		for (int i = 0; i < m; i++)
		{
			s << "Fiber #: " << i + 1 << endln;
			s << "Concrete material with tag: " << theMaterialsConcrete[i]->getTag() << endln;
			theMaterialsConcrete[i]->Print(s, flag);

			s << "Steel material with tag: " << theMaterialsSteel[i]->getTag() << endln;
			theMaterialsSteel[i]->Print(s, flag);

		}

		s << "Shear material with tag: " << theMaterialsShear[0]->getTag() << endln;
		theMaterialsShear[0]->Print(s, flag);

	}
	else if (flag == 1) {
		// does nothing
	}
}


// Set recorders
Response *FourNodeMVLEM3D::setResponse(const char **argv, int argc, OPS_Stream &s)
{
	Response *theResponse = 0;

	s.tag("ElementOutput");
	s.attr("eleType", "FourNodeMVLEM3D");
	s.attr("eleTag", this->getTag());
	s.attr("node1", externalNodes[0]);
	s.attr("node2", externalNodes[1]);
	s.attr("node3", externalNodes[3]);
	s.attr("node4", externalNodes[2]);

	// Nodal forces in global cs
	if (strcmp(argv[0], "forceG") == 0 || strcmp(argv[0], "forcesG") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

		s.tag("ResponseType", "globalFx_i");
		s.tag("ResponseType", "globalFy_i");
		s.tag("ResponseType", "globalFz_i");
		s.tag("ResponseType", "globalMx_i");
		s.tag("ResponseType", "globalMy_i");
		s.tag("ResponseType", "globalMz_i");
		s.tag("ResponseType", "globalFx_j");
		s.tag("ResponseType", "globalFy_j");
		s.tag("ResponseType", "globalFz_j");
		s.tag("ResponseType", "globalMx_j");
		s.tag("ResponseType", "globalMy_j");
		s.tag("ResponseType", "globalMz_j");
		s.tag("ResponseType", "globalFx_k");
		s.tag("ResponseType", "globalFy_k");
		s.tag("ResponseType", "globalFz_k");
		s.tag("ResponseType", "globalMx_k");
		s.tag("ResponseType", "globalMy_k");
		s.tag("ResponseType", "globalMz_k");
		s.tag("ResponseType", "globalFx_l");
		s.tag("ResponseType", "globalFy_l");
		s.tag("ResponseType", "globalFz_l");
		s.tag("ResponseType", "globalMx_l");
		s.tag("ResponseType", "globalMy_l");
		s.tag("ResponseType", "globalMz_l");

		return theResponse = new ElementResponse(this, 1, Vector(24));

	}

	// Nodal forces in local cs
	else if (strcmp(argv[0], "forceL") == 0 || strcmp(argv[0], "forcesL") == 0 ||
		strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

		s.tag("ResponseType", "localFx_i");
		s.tag("ResponseType", "localFy_i");
		s.tag("ResponseType", "localFz_i");
		s.tag("ResponseType", "localMx_i");
		s.tag("ResponseType", "localMy_i");
		s.tag("ResponseType", "localMz_i");
		s.tag("ResponseType", "localFx_j");
		s.tag("ResponseType", "localFy_j");
		s.tag("ResponseType", "localFz_j");
		s.tag("ResponseType", "localMx_j");
		s.tag("ResponseType", "localMy_j");
		s.tag("ResponseType", "localMz_j");
		s.tag("ResponseType", "localFx_k");
		s.tag("ResponseType", "localFy_k");
		s.tag("ResponseType", "localFz_k");
		s.tag("ResponseType", "localMx_k");
		s.tag("ResponseType", "localMy_k");
		s.tag("ResponseType", "localMz_k");
		s.tag("ResponseType", "localFx_l");
		s.tag("ResponseType", "localFy_l");
		s.tag("ResponseType", "localFz_l");
		s.tag("ResponseType", "localMx_l");
		s.tag("ResponseType", "localMy_l");
		s.tag("ResponseType", "localMz_l");

		return theResponse = new ElementResponse(this, 2, Vector(24));

	}

	// Element curvature
	else if (strcmp(argv[0], "Curvature") == 0 || strcmp(argv[0], "curvature") == 0) {

		s.tag("ResponseType", "fi");

		return theResponse = new ElementResponse(this, 3, 0.0);
	}

	// Fiber strain
	else if (strcmp(argv[0], "Fiber_Strain") == 0 || strcmp(argv[0], "fiber_strain") == 0) {

		s.tag("ResponseType", "epsy");

		return theResponse = new ElementResponse(this, 4, Vector(m));
	}

	// Fiber concrete stresses
	else if (strcmp(argv[0], "Fiber_Stress_Concrete") == 0 || strcmp(argv[0], "fiber_stress_concrete") == 0) {

		s.tag("ResponseType", "sigmayc");

		return theResponse = new ElementResponse(this, 5, Vector(m));
	}

	// Fiber steel stresses
	else if (strcmp(argv[0], "Fiber_Stress_Steel") == 0 || strcmp(argv[0], "fiber_stress_steel") == 0) {

		s.tag("ResponseType", "sigmays");

		return theResponse = new ElementResponse(this, 6, Vector(m));
	}

	// Shear force deformation
	else if (strcmp(argv[0], "Shear_Force_Deformation") == 0 || strcmp(argv[0], "shear_force_deformation") == 0) {

		s.tag("ResponseType", "shearFD");

		return theResponse = new ElementResponse(this, 7, Vector(2));
	}

	// Shear Deformation
	else if (strcmp(argv[0], "ShearDef") == 0 || strcmp(argv[0], "sheardef") == 0) {

		s.tag("ResponseType", "shearDef");

		return theResponse = new ElementResponse(this, 8, 0.0);
	}

	s.endTag();

	return 0;
}

// get recorders
int FourNodeMVLEM3D::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID)
	{
	case 1:  // Nodal forces in global cs
		return eleInfo.setVector(this->getResistingForce());

	case 2:  // Nodal forces in local cs
		return eleInfo.setVector(this->getResistingForceLocal());

	case 3:  // Curvature
		return eleInfo.setDouble(this->getCurvature());

	case 4:  // Fiber Strains
		return eleInfo.setVector(this->getStrain());

	case 5:  // Fiber Concrete Stress
		return eleInfo.setVector(this->getStressConcrete());

	case 6:  // Fiber Steel Stress
		return eleInfo.setVector(this->getStressSteel());

	case 7:  // Shear Force-Deformtion
		return eleInfo.setVector(this->getShearFD());

	case 8:  // Shear Deformtion
		return eleInfo.setVector(this->getShearDef());

	default:

		return 0;

	}
}

// !!! test out all recorders (particularly force)
// !!! add another recorder or a switch where user can selct which forces at hodes to record
// Return element local forces !!! check if this is OK
Vector FourNodeMVLEM3D::getResistingForceLocal(void)
{
	return FourNodeMVLEM3DRlocal;
}

// Get curvature (from vertical strains)
double FourNodeMVLEM3D::getCurvature(void)
{
	double Curv;

	Curv = (FourNodeMVLEM3DStrain[0] - FourNodeMVLEM3DStrain[m - 1]) / (x[0] - x[m - 1]);

	return Curv;
}

// Get fiber strains
Vector FourNodeMVLEM3D::getStrain(void)
{
	Vector fiberStrain(m);

	for (int i = 0; i<m; i++) {
		fiberStrain(i) = FourNodeMVLEM3DStrain[i];
	}

	return fiberStrain;
}

// Get Concrete Stress 
Vector FourNodeMVLEM3D::getStressConcrete(void)
{
	Vector concreteStress(m);

	for (int i = 0; i<m; i++) {
		concreteStress(i) = theMaterialsConcrete[i]->getStress();
	}

	return concreteStress;
}

// Get Steel Stress 
Vector FourNodeMVLEM3D::getStressSteel(void)
{
	Vector steelStress(m);

	for (int i = 0; i<m; i++) {
		steelStress(i) = theMaterialsSteel[i]->getStress();
	}

	return steelStress;
}

// Get Shear Stress-Strain 
Vector FourNodeMVLEM3D::getShearFD(void)
{
	Vector shearStrainStress(2);

	shearStrainStress(0) = theMaterialsShear[0]->getStrain();
	shearStrainStress(1) = theMaterialsShear[0]->getStress();

	return shearStrainStress;
}

// Get Shear Stress-Strain 
double FourNodeMVLEM3D::getShearDef(void)
{
	double shearDef;

	shearDef = theMaterialsShear[0]->getStrain();

	return shearDef;
}

// Compute element transformation matrix
void  FourNodeMVLEM3D::setTransformationMatrix(void) {

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

	// Magnitude
	X_ = pow(pow(Xx, 2.0) + pow(Xy, 2.0) + pow(Xz, 2.0), 0.5);

	// unit x components
	Xex = Xx / X_;
	Xey = Xy / X_;
	Xez = Xz / X_;

	// k------l
	// |      |
	// |      |
	// i------j
	// Components of local Y axis
	Yx = nd3Crds(0) - nd1Crds(0);
	Yy = nd3Crds(1) - nd1Crds(1);
	Yz = nd3Crds(2) - nd1Crds(2);

	// Magnitude
	Y_ = pow(pow(Yx, 2.0) + pow(Yy, 2.0) + pow(Yz, 2.0), 0.5);

	// unit y components
	Yex = Yx / Y_;
	Yey = Yy / Y_;
	Yez = Yz / Y_;

	// (Ze) = (Xe) x (Ye)
	Zex = Xey*Yez - Xez*Yey;
	Zey = -(Xex*Yez - Xez*Yex);
	Zez = Xex*Yey - Xey*Yex;

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

	// DEBUG
	/*opserr << "T :";
	opserr << "\n";
	for (int j = 0; j < 12; j++) {
	for (int i = 0; i < 12; i++) {
	opserr << T(i, j) << " , ";
	}
	opserr << "\n";
	}
	opserr << "\n" << "\n";*/

	/*opserr << "Tt :";
	opserr << "\n";
	for (int j = 0; j < 3; j++) {
	for (int i = 0; i < 3; i++) {
	opserr << Tt(i, j) << " , ";
	}
	opserr << "\n";
	}
	opserr << "\n" << "\n";*/

}