// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								Carlos N. López
//								Leonardo M. Massone
//								 
// User documentation available at: https://kkolozvari.github.io/E-SFI-MVLEM-3D/
//
// Created: 08/2023
//
// Description: The E-SFI-MVLEM-3D model is a three-dimensional four-node element with 24 DOFs that incorporates axial-flexural-shear interaction and 
// can be used for nonlinear analysis of non-rectangular reinforced concrete walls subjected to multidirectional loading. The E-SFI-MVLEM-3D model is derived by combining two 
// previously available models, a two-dimensional (E-SFI) model, and a three-dimensional (SFI-MVLEM-3D) model. The major enhancement in the model formulation compared to its 
// parent SFI-MVLEM-3D comes from implementing a closed-form solution for calculating horizontal axial strains at fibers of the wall element. This significantly reduced the
// number of element degrees of freedom, which resulted in analysis run-time that is reduced to approximately 25% and a convergence rate that is increased roughly two times.
//
// Notes:
// Nodes should be assigned in counterclockwise direction.
//    4........3 
//    .        .
//    .        .
//    .        . ^ y
//    1........2 |-> x
//
// Reference:
// Kristijan Kolozvari, Carlos N. López, Leonardo M. Massone (2023), “Efficient Three-dimensional Shear-flexure Interaction Model for Reinforced Concrete Walls”, Engineering Structures, Vol. 294, 116700, https://doi.org/10.1016/j.engstruct.2023.116700.
//
// User documentation available at: https://kkolozvari.github.io/E-SFI-MVLEM-3D/
//
// Source: /usr/local/cvs/OpenSees/SRC/element/mvlem/E_SFI_MVLEM_3D.cpp
//
// Rev: 1.0
#include <G3Globals.h>
#include <UniaxialMaterial.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>
#include "E_SFI_MVLEM_3D.h"
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

#include <elementAPI.h>

// Read input parameters and build the element
void* OPS_E_SFI_MVLEM_3D(void)
{
	// Pointer to a uniaxial material that will be returned                       
	Element* theElement = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numArgs < 14) {
		opserr << "Want: element E_SFI_MVLEM_3D eleTag iNode jNode kNode lNode m -thick {Thicknesses} -width {Widths} -mat {Material_tags} <-CoR c> <-ThickMod tMod> <-Poisson Nu> <-Density Dens>\n";
		return 0;
	}

	int iData[6];
	double dData[4];

	// set defaults
	dData[0] = 0.4;		// c
	dData[1] = 0.63;	// tMod (equivalent to cracked out-of-plan stiffness of 0.25Ig)
	dData[2] = 0.25;	// Poisson (concrete)
	dData[3] = 0.0;		// Density

	int numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid tag for element E_SFI_MVLEM_3D" << endln;
		return 0;
	}

	numData = 5;
	if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
		opserr << "WARNING iNode jNode kNode lNode or m for element E_SFI_MVLEM_3D" << iData[0] << endln;
		return 0;
	}

	int m = iData[5];
	const char* str = 0;

	double* theThickness = new double[m];
	double* theWidth = new double[m];
	int* matTags = new int[m];

	NDMaterial** theMaterials = new NDMaterial * [m];
	for (int i = 0; i < m; i++) {
		theThickness[i] = 0.0;
		theWidth[i] = 0.0;
		matTags[i] = 0;
		theMaterials[i] = 0;
	}

	numArgs = OPS_GetNumRemainingInputArgs();

	while (numArgs > 0) {
		//OPS_GetStringCopy(&str);
		str = OPS_GetString();
		if (strcmp(str, "-thick") == 0) {
			numData = m;
			if (OPS_GetDoubleInput(&numData, theThickness) != 0) {
				opserr << "Invalid thick parameter for E_SFI_MVLEM_3D   " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-width") == 0) {
			numData = m;
			if (OPS_GetDoubleInput(&numData, theWidth) != 0) {
				opserr << "Invalid width value for E_SFI_MVLEM_3D  " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-mat") == 0) {
			numData = m;
			if (OPS_GetIntInput(&numData, matTags) != 0) {
				opserr << "Invalid mat tags for E_SFI_MVLEM_3D  " << iData[0] << endln;
				return 0;
			}
			for (int i = 0; i < m; i++) {
				theMaterials[i] = 0;
				theMaterials[i] = OPS_getNDMaterial(matTags[i]);
				if (theMaterials[i] == 0) {
					opserr << "Invalid material tag " << matTags[i] << "  for E_SFI_MVLEM_3D  " << iData[0] << endln;
					return 0;
				}
			}
		}

		// optional parameters
		else if (strcmp(str, "-CoR") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
				opserr << "Invalid CoR parameter for E_SFI_MVLEM_3D   " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-ThickMod") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
				opserr << "Invalid ThickMod parameter for E_SFI_MVLEM_3D   " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-Poisson") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[2]) != 0) {
				opserr << "Invalid Poisson parameter for E_SFI_MVLEM_3D   " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-Density") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[3]) != 0) {
				opserr << "Invalid Dens parameter for E_SFI_MVLEM_3D   " << iData[0] << endln;
				return 0;
			}
		}
		numArgs = OPS_GetNumRemainingInputArgs();

	}

	theElement = new E_SFI_MVLEM_3D(iData[0], dData[3],
		iData[1], iData[2], iData[3], iData[4],
		theMaterials,
		theThickness,
		theWidth,
		iData[5], dData[0], dData[2], dData[1]);

	// Cleanup dynamic memory
	if (theThickness != 0)
		delete[] theThickness;
	if (theWidth != 0)
		delete[] theWidth;
	if (matTags != 0)
		delete[] matTags;

	if (theMaterials != 0)
		delete[] theMaterials;

	return theElement;
}

// Typical constructor
E_SFI_MVLEM_3D::E_SFI_MVLEM_3D(int tag,
	double Dens,
	int Nd1, int Nd2, int Nd3, int Nd4,
	NDMaterial** materials,
	double* thickness,
	double* width,
	int mm = 0,
	double cc = 0.0,
	double nn = 0.0,
	double tf = 0.0)

	:Element(tag, ELE_TAG_E_SFI_MVLEM_3D),
	density(Dens),
	externalNodes(4),
	theMaterial(0), theLoad(0),
	E_SFI_MVLEM_3DStrainX(0), E_SFI_MVLEM_3DStrainY(0), E_SFI_MVLEM_3DStrainXY(0), E_SFI_MVLEM_3DStrain(0),
	x(0), b(0), AcX(0), AcY(0), kx(0), ky(0), Kh(0), Fx(0), Fy(0), Fxy(0), Dx(0), Dy(0), Dxy(0),
	rhox(0), // for epsX expression
	E_SFI_MVLEM_3DK(24, 24), E_SFI_MVLEM_3DR(24), E_SFI_MVLEM_3DD(24, 24), E_SFI_MVLEM_3DM(24, 24),
	E_SFI_MVLEM_3DKlocal(24, 24), E_SFI_MVLEM_3DDlocal(24, 24), E_SFI_MVLEM_3DRlocal(24), E_SFI_MVLEM_3DMlocal(24, 24),
	P_24DOF(24), P_24DOF_local(24),
	m(mm), c(cc), NUelastic(nn), Tfactor(tf),
	T(24, 24), Tt(3, 3), T6(6, 6),
	nd1Crds(3), nd2Crds(3), nd3Crds(3), nd4Crds(3), modifiedT(0), t(0)
{

	TotalMass = 0.0;
	NodeMass = 0.0;
	h = 0.0;
	d = 0.0;
	Lw = 0.0;

	// Out of Plane parameters
	Eave = 0.0;
	Tave = 0.0;

	// Imaginary beam properties
	Eib = 0.0;
	Hib = 0.0;
	Iib = 0.0;
	Aib = 0.0;

	// Check number of fibers - max is 999
	if (m > 999) {
		opserr << "WARNING: Number of fibers assigned is " << m << ". Maximum allowed number of fibers is 999!\n";
		exit(-1);
	}

	// Fill in the ID containing external node info with node id's 
	if (externalNodes.Size() != 4)
		opserr << "FATAL E_SFI_MVLEM_3D::E_SFI_MVLEM_3D() - out of memory, could not create an ID of size 4\n";

	// Assign node tags to external nodes - Node ordering switched on purpose to match the theoretical derivation
	externalNodes(0) = Nd1;
	externalNodes(1) = Nd2;
	externalNodes(3) = Nd3;
	externalNodes(2) = Nd4;

	// Set external node pointers to NULL - external nodes
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;

	// Check thickness and width input
	if (thickness == 0) {
		opserr << "E_SFI_MVLEM_3D::E_SFI_MVLEM_3D() - "
			<< "Null thickness array passed.\n";
		exit(-1);
	}

	if (width == 0) {
		opserr << "E_SFI_MVLEM_3D::E_SFI_MVLEM_3D() - "
			<< "Null width array passed.\n";
		exit(-1);
	}

	// Allocate memory for the thickness and width
	// Input parameters
	t = new double[m];
	b = new double[m];

	for (int i = 0; i < m; i++) {
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
		for (int j = 0; j < i + 1; j++)
			sumb_i += b[j];

		x[i] = (sumb_i - b[i] / 2.0) - Lw / 2.0;
	}

	// Check material input
	if (materials == 0) {
		opserr << "E_SFI_MVLEM_3D::E_SFI_MVLEM_3D() - "
			<< "Null material array passed.\n";
		exit(-1);
	}

	// Allocate memory for the ND materials
	theMaterial = new NDMaterial * [m];

	if (theMaterial == 0) {
		opserr << "E_SFI_MVLEM_3D::E_SFI_MVLEM_3D() - "
			<< "Failed to allocate pointers for uniaxial materials.\n";
		exit(-1);
	}

	// Get copies of the ND materials
	for (int i = 0; i < m; i++) {
		if (materials[i] == 0) {
			opserr << "E_SFI_MVLEM_3D::E_SFI_MVLEM_3D() - "
				"Null ND material pointer passed.\n";
			exit(-1);
		}

		theMaterial[i] = materials[i]->getCopy("PlaneStress2D");

		if (theMaterial[i] == 0) {
			opserr << "E_SFI_MVLEM_3D::E_SFI_MVLEM_3D() - "
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

	// Panel force (trial)
	Fx = new double[m];
	Fy = new double[m];
	Fxy = new double[m];

	// Panel stiffness (trial)
	Dx = new double[m];
	Dy = new double[m];
	Dxy = new double[m];

	// Panel strains
	E_SFI_MVLEM_3DStrainX = new double[m];
	E_SFI_MVLEM_3DStrainY = new double[m];
	E_SFI_MVLEM_3DStrainXY = new double[m];
	E_SFI_MVLEM_3DStrain = new double[3 * m];

	// for epsX expression
	rhox = new double[m];

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

		E_SFI_MVLEM_3DStrainX[i] = 0.0;
		E_SFI_MVLEM_3DStrainY[i] = 0.0;
		E_SFI_MVLEM_3DStrainXY[i] = 0.0;

		E_SFI_MVLEM_3DStrain[i] = 0.0;
		E_SFI_MVLEM_3DStrain[i + m] = 0.0;
		E_SFI_MVLEM_3DStrain[i + 2 * m] = 0.0;

		rhox[i] = 0.0;

	}

	Kh = 0.0;

	// Revert to start
	this->revertToStart();
}

// Constructor which should be invoked by an FE_ObjectBroker only
E_SFI_MVLEM_3D::E_SFI_MVLEM_3D()
	:Element(0, ELE_TAG_E_SFI_MVLEM_3D),
	externalNodes(4),
	theMaterial(0), theLoad(0),
	E_SFI_MVLEM_3DStrainX(0), E_SFI_MVLEM_3DStrainY(0), E_SFI_MVLEM_3DStrainXY(0), E_SFI_MVLEM_3DStrain(0),
	x(0), b(0), AcX(0), AcY(0), kx(0), ky(0), Kh(0), Fx(0), Fy(0), Fxy(0), Dx(0), Dy(0), Dxy(0),
	rhox(0), // for epsX expression
	E_SFI_MVLEM_3DK(24, 24), E_SFI_MVLEM_3DR(24), E_SFI_MVLEM_3DD(24, 24), E_SFI_MVLEM_3DM(24, 24),
	E_SFI_MVLEM_3DKlocal(24, 24), E_SFI_MVLEM_3DDlocal(24, 24), E_SFI_MVLEM_3DRlocal(24), E_SFI_MVLEM_3DMlocal(24, 24),
	P_24DOF(24), P_24DOF_local(24),
	m(m), c(c), NUelastic(NUelastic), Tfactor(Tfactor),
	T(24, 24), Tt(3, 3), T6(6, 6),
	nd1Crds(3), nd2Crds(3), nd3Crds(3), nd4Crds(3), modifiedT(0), t(0)
{
	if (externalNodes.Size() != 4)
		opserr << "FATAL E_SFI_MVLEM_3D::E_SFI_MVLEM_3D() - out of memory, could not create an ID of size 4\n";

	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;

}

//  Destructor - provided to clean up any memory 
E_SFI_MVLEM_3D::~E_SFI_MVLEM_3D()
{
	// clean up the memory associated with the element, this is
	// memory the E_SFI_MVLEM_3D objects allocates and memory allocated 
	// by other objects that the E_SFI_MVLEM_3D object is responsible for 
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
		delete[]x;
	if (b != 0)
		delete[]b;
	if (AcX != 0)
		delete[]AcX;
	if (AcY != 0)
		delete[]AcY;
	if (kx != 0)
		delete[]kx;
	if (ky != 0)
		delete[]ky;
	if (Fx != 0)
		delete[]Fx;
	if (Fy != 0)
		delete[]Fy;
	if (Fxy != 0)
		delete[]Fxy;
	if (Dx != 0)
		delete[]Dx;
	if (Dy != 0)
		delete[]Dy;
	if (Dxy != 0)
		delete[]Dxy;
	if (rhox != 0)
		delete rhox;
	if (E_SFI_MVLEM_3DStrainX != 0)
		delete[]E_SFI_MVLEM_3DStrainX;
	if (E_SFI_MVLEM_3DStrainY != 0)
		delete[]E_SFI_MVLEM_3DStrainY;
	if (E_SFI_MVLEM_3DStrainXY != 0)
		delete[]E_SFI_MVLEM_3DStrainXY;
	if (E_SFI_MVLEM_3DStrain != 0)
		delete[]E_SFI_MVLEM_3DStrain;
	if (modifiedT != 0)
		delete[]modifiedT;
	if (t != 0)
		delete[]t;

}

// Get number of nodes (external + internal)
int E_SFI_MVLEM_3D::getNumExternalNodes(void) const
{
	return 4;
}

// Get node tags
const ID& E_SFI_MVLEM_3D::getExternalNodes(void)
{
	return externalNodes;
}

// Get node pointers
Node** E_SFI_MVLEM_3D::getNodePtrs(void)
{
	return theNodes;
}

// Get number of DOFs 
int E_SFI_MVLEM_3D::getNumDOF(void) {

	int NumDOF = 24; // 6 DOFs per external node x 4 external nodes

	return NumDOF;
}

// Set Domain
void E_SFI_MVLEM_3D::setDomain(Domain* theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0)
	{
		return;
	}

	// Set node pointers to NULL
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;

	// First ensure external nodes exist in Domain and set the node pointers
	int Nd1 = externalNodes(0);
	int Nd2 = externalNodes(1);
	int Nd3 = externalNodes(2);
	int Nd4 = externalNodes(3);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);
	theNodes[2] = theDomain->getNode(Nd3);
	theNodes[3] = theDomain->getNode(Nd4);

	// Get coordinates of end nodes
	nd1Crds = theNodes[0]->getCrds();
	nd2Crds = theNodes[1]->getCrds();
	nd3Crds = theNodes[2]->getCrds();
	nd4Crds = theNodes[3]->getCrds();

	// Compute coordinate transformation matrix
	setTransformationMatrix();

	// Node coordinates in local coordinate system
	Vector nd1CrdL(3); nd1CrdL.Zero();
	Vector nd2CrdL(3); nd2CrdL.Zero();
	Vector nd3CrdL(3); nd3CrdL.Zero();
	Vector nd4CrdL(3); nd4CrdL.Zero();

	nd1CrdL.addMatrixVector(0.0, Tt, nd1Crds, 1.0);
	nd2CrdL.addMatrixVector(0.0, Tt, nd2Crds, 1.0);
	nd3CrdL.addMatrixVector(0.0, Tt, nd3Crds, 1.0);
	nd4CrdL.addMatrixVector(0.0, Tt, nd4Crds, 1.0);

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
	for (int i = 0; i < m; i++) {
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

	if (theNodes[0] == 0) {
		opserr << "WARNING E_SFI_MVLEM_3D::setDomain() - at E_SFI_MVLEM_3D " << this->getTag() << " node " <<
			Nd1 << " does not exist in domain\n";
		return;  // Don't go any further - otherwise segemntation fault
	}

	if (theNodes[1] == 0) {
		opserr << "WARNING E_SFI_MVLEM_3D::setDomain() - at E_SFI_MVLEM_3D " << this->getTag() << " node " <<
			Nd2 << " does not exist in domain\n";
		return;
	}

	if (theNodes[2] == 0) {
		opserr << "WARNING E_SFI_MVLEM_3D::setDomain() - at E_SFI_MVLEM_3D " << this->getTag() << " node " <<
			Nd3 << " does not exist in domain\n";
		return;  // Don't go any further - otherwise segemntation fault
	}
	if (theNodes[3] == 0) {
		opserr << "WARNING E_SFI_MVLEM_3D::setDomain() - at E_SFI_MVLEM_3D " << this->getTag() << " node " <<
			Nd4 << " does not exist in domain\n";
		return;
	}

	// Call the DomainComponent class method 
	this->DomainComponent::setDomain(theDomain);

	// Ensure connected nodes have correct number of dof's
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();
	int dofNd3 = theNodes[2]->getNumberDOF();
	int dofNd4 = theNodes[3]->getNumberDOF();

	if ((dofNd1 != 6) || (dofNd2 != 6) || (dofNd3 != 6) || (dofNd4 != 6)) {
		opserr << "E_SFI_MVLEM_3D::setDomain(): 6 dof required at all nodes. " << dofNd1 << " provided at node 1, " << dofNd2 << " provided at node 2, "
			<< dofNd3 << " provided at node 4, " << dofNd4 << " provided at node 3";
	}

	// Calculate concrete areas in X and Y directions
	double A = 0.0;
	for (int i = 0; i < m; i++) {
		AcX[i] = h * t[i];
		AcY[i] = b[i] * t[i];

		A += AcY[i];
	}

	// Determine the nodal mass for lumped mass approach
	A = 0;
	for (int i = 0; i < m; i++) {
		A += b[i] * t[i];
	}

	NodeMass = density * A * h / 4.0;

	// Get Concrete Young's Modulus
	//theResponses = new Response *[1];
	//if (theResponses == 0) {
	//	opserr << " E_SFI_MVLEM_3D::E_SFI_MVLEM_3D - failed allocate responses array\n";
	//	exit(-1);
	//}

	//OPS_Stream *theDummyStream = new DummyStream();
	DummyStream theDummyStream;
	//const char **argv = new const char *[1];
	//argv[0] = "getInputParameters"; // to get input parameters from concrete material
	char aa[80] = "getInputParameters";
	const char* argv[1];
	argv[0] = aa;

	for (int i = 0; i < m; i++)
	{
		//theResponses[0] = theMaterial[i]->setResponse(argv, 1, *theDummyStream);
		Response* theResponse = theMaterial[i]->setResponse(argv, 1, theDummyStream);

		//if (theResponses[0] == 0) {
		if (theResponse == 0) {
			opserr << " E_SFI_MVLEM_3D::E_SFI_MVLEM_3D - failed to get input parameters for FSAM material with tag: " << this->getTag() << "\n";
			exit(-1);
		}

		// Get FSAM material input variables
		//theResponses[0]->getResponse();
		//Information &theInfoInput = theResponses[0]->getInformation();		
		theResponse->getResponse();
		Information& theInfoInput = theResponse->getInformation();
		const Vector& InputNDMat = theInfoInput.getData();

		// Calculate out-of-plane modulus of elasticity (average modulus)
		Eave += AcY[i] * InputNDMat[9] / A;

		// Assign parameters for EpsX formula calculsions
		rhox[i] = InputNDMat[3];

		delete theResponse;

	}

	//delete theDummyStream;

	// Internal beam parameters
	Eib = Eave;
	Hib = h;
	Aib = Tave * Hib;
	Iib = 0.5 * (Tave * Hib * Hib * Hib / 12.0);

	Tave *= Tfactor; // multiply Tave with the modification factor

	// Create a vector to hop applied loads - NOT used in the current model formulation (no element loads)
	if (theLoad == 0)
		theLoad = new Vector(24);
	if (theLoad == 0) {
		opserr << "E_SFI_MVLEM_3D::setDomain() - element: " << this->getTag()
			<< " out of memory creating vector of size: " << 24 << endln;
		return;
	}

	// Calculate constant terms of stiffness matrix
	K1 = -(Eave * (Tave * Tave * Tave) * (10.0 * (h * h * h * h) + 10.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	K2 = (Eave * (Tave * Tave * Tave) * (4.0 * (h * h) * NUelastic + (h * h) + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0));
	K3 = (Eave * (Tave * Tave * Tave) * (4.0 * NUelastic * (Lw * Lw) + 10.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	K4 = (Eave * (Tave * Tave * Tave) * (10.0 * (h * h * h * h) - 5.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	K5 = (Eave * (Tave * Tave * Tave) * (4.0 * (h * h) * NUelastic + (h * h) - 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0));
	K6 = (Eave * (Tave * Tave * Tave) * (10.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	K7 = -(Eave * (Tave * Tave * Tave) * (5.0 * (h * h * h * h) - 10.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	K8 = (Eave * (Tave * Tave * Tave) * ((h * h) - (h * h) * NUelastic + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0));
	K9 = (Eave * (Tave * Tave * Tave) * (4.0 * NUelastic * (Lw * Lw) - 5.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	K10 = (Eave * (Tave * Tave * Tave) * (5.0 * (h * h * h * h) + 5.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	K11 = (Eave * (Tave * Tave * Tave) * ((h * h) * NUelastic - (h * h) + 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0));
	K12 = (Eave * (Tave * Tave * Tave) * (NUelastic * (Lw * Lw) + 5.0 * (h * h) - (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	K13 = -(Eave * (Tave * Tave * Tave) * ((h * h) - (h * h) * NUelastic + 5.0 * (Lw * Lw))) / (45.0 * h * Lw * ((NUelastic * NUelastic) - 1.0));
	K14 = (Eave * NUelastic * (Tave * Tave * Tave)) / (12.0 * (NUelastic * NUelastic) - 12.0);
	K15 = -(Eave * (Tave * Tave * Tave) * (2.0 * (h * h) * NUelastic - 2.0 * (h * h) + 5.0 * (Lw * Lw))) / (90.0 * h * Lw * ((NUelastic * NUelastic) - 1.0));
	K16 = -(Eave * (Tave * Tave * Tave) * ((h * h) * NUelastic - (h * h) + 10.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0));
	K17 = -(Eave * (Tave * Tave * Tave) * ((h * h) - (h * h) * NUelastic + 5.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0));
	K18 = (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (45.0 * h * ((NUelastic * NUelastic) - 1.0)) - (Eave * h * (Tave * Tave * Tave)) / (9.0 * Lw * ((NUelastic * NUelastic) - 1.0));
	K19 = -(Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0));
	K20 = -(Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * Lw * (2.0 * NUelastic - 2.0)) / (90.0 * h * ((NUelastic * NUelastic) - 1.0));
	K21 = (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0)) - (Eave * h * (Tave * Tave * Tave)) / (36.0 * Lw * ((NUelastic * NUelastic) - 1.0));
	K22 = -(Eave * (Tave * Tave * Tave) * (5.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (45.0 * h * Lw * ((NUelastic * NUelastic) - 1.0));

}

// Commit state of the materials
int E_SFI_MVLEM_3D::commitState()
{
	int errCode = 0;

	// Commit material models
	for (int i = 0; i < m; i++) {
		errCode += theMaterial[i]->commitState();
	}

	return errCode;
}

// Revert to last committed state (if convergence is not achieved)
int E_SFI_MVLEM_3D::revertToLastCommit()
{
	int errCode = 0;

	// Revert material models
	for (int i = 0; i < m; i++) {
		errCode += theMaterial[i]->revertToLastCommit();
	}

	return errCode;
}

// Revert to start
int E_SFI_MVLEM_3D::revertToStart()
{
	int errCode = 0;

	// Revert material models
	for (int i = 0; i < m; i++)
		errCode += theMaterial[i]->revertToStart();

	return errCode;

}

// Update state
int E_SFI_MVLEM_3D::update()
{
	// Get the current strain given trial displacements at nodes 
	this->computeCurrentStrain();

	// Set the strain in the materials
	int errCode = 0;

	for (int i = 0; i < m; i++) {

		Vector strain(3);

		strain(0) = E_SFI_MVLEM_3DStrain[i];
		strain(1) = E_SFI_MVLEM_3DStrain[i + m];
		strain(2) = E_SFI_MVLEM_3DStrain[i + 2 * m];

		// Set trial response for material models
		errCode += theMaterial[i]->setTrialStrain(strain);

	}

	return errCode;
}

// Get current strains at RC panels (macro-fibers)
double* E_SFI_MVLEM_3D::computeCurrentStrain(void)
{

	const Vector& disp1 = theNodes[0]->getTrialDisp();
	const Vector& disp2 = theNodes[1]->getTrialDisp();
	const Vector& disp3 = theNodes[2]->getTrialDisp();
	const Vector& disp4 = theNodes[3]->getTrialDisp();

	Vector dispG(24); // Vector of total 24 displacemets in global coordinates
	dispG.Zero();
	Vector dispL(24); // Vector of total 24 displacemets in local coordinates
	dispL.Zero();
	Vector dispL_inPlan2N(6); // in-plane displacements of equlivalent 2D, 2N SFI model
	dispL_inPlan2N.Zero();

	// store nodal displacemnts in global vs
	for (int i = 0; i < 6; i++) {
		dispG(i) = disp1(i);
		dispG(i + 6) = disp2(i);
		dispG(i + 12) = disp3(i);
		dispG(i + 18) = disp4(i);
	}

	// transform nodal displacements from global to local cs
	dispL.addMatrixVector(0.0, T, dispG, 1.0);

	dispL_inPlan2N(0) = dispL(0) / 2.0 + dispL(6) / 2.0;
	dispL_inPlan2N(1) = dispL(1) / 2.0 + dispL(7) / 2.0;
	dispL_inPlan2N(2) = dispL(5) / (2.0 * (d * d) + 2.0) + dispL(11) / (2.0 * (d * d) + 2.0) - (dispL(1) * d) / (2.0 * (d * d) + 2.0) + (dispL(7) * d) / (2.0 * (d * d) + 2.0);
	dispL_inPlan2N(3) = dispL(12) / 2.0 + dispL(18) / 2.0;
	dispL_inPlan2N(4) = dispL(13) / 2.0 + dispL(19) / 2.0;
	dispL_inPlan2N(5) = dispL(17) / (2.0 * (d * d) + 2.0) + dispL(23) / (2.0 * (d * d) + 2.0) - (dispL(13) * d) / (2.0 * (d * d) + 2.0) + (dispL(19) * d) / (2.0 * (d * d) + 2.0);

	// Deformations at each RC panel (macro-fiber) - MVLEM formulation
	for (int i = 0; i < m; i++) {
		Dy[i] = -dispL_inPlan2N(1) - x[i] * dispL_inPlan2N(2) + dispL_inPlan2N(4) + x[i] * dispL_inPlan2N(5);
		Dxy[i] = dispL_inPlan2N(0) - dispL_inPlan2N(3) - c * h * dispL_inPlan2N(2) - (1.0 - c) * h * dispL_inPlan2N(5);
	}

	Dsh = -Dxy[0]; // Store shear deformations for the recorder

	// Strains at each RC panel (macro-fiber)
	for (int i = 0; i < m; i++) {
		E_SFI_MVLEM_3DStrainY[i] = Dy[i] / h;
		E_SFI_MVLEM_3DStrainXY[i] = -Dxy[i] / h;
		E_SFI_MVLEM_3DStrainX[i] = 0.55 * pow((1.0 + rhox[i]), -60.0) * (1.0 - pow(3.0, -800.0 * abs(E_SFI_MVLEM_3DStrainXY[i]))) * abs(E_SFI_MVLEM_3DStrainXY[i]); // Massone et al.
	}

	// Store strains into a single vector
	for (int i = 0; i < m; i++) {
		E_SFI_MVLEM_3DStrain[i] = E_SFI_MVLEM_3DStrainX[i];
		E_SFI_MVLEM_3DStrain[i + m] = E_SFI_MVLEM_3DStrainY[i];
		E_SFI_MVLEM_3DStrain[i + 2 * m] = E_SFI_MVLEM_3DStrainXY[i];
	}

	// Return strain vector
	return E_SFI_MVLEM_3DStrain;

}

// Get the element initial element tangent matrix
const Matrix& E_SFI_MVLEM_3D::getInitialStiff(void)
{

	E_SFI_MVLEM_3DK.Zero();		// Global stiffness matrix
	E_SFI_MVLEM_3DKlocal.Zero();	// Local stiffness matrix

	Kh = 0.0;

	for (int i = 0; i < m; i++)
	{
		// Get material initial tangent
		const Matrix& D = theMaterial[i]->getInitialTangent();

		double D00 = D(0, 0); double D01 = D(0, 1); double D02 = D(0, 2);
		double D10 = D(1, 0); double D11 = D(1, 1); double D12 = D(1, 2);
		double D20 = D(2, 0); double D21 = D(2, 1); double D22 = D(2, 2);

		kx[i] = D00 * h * t[i] / b[i];
		ky[i] = D11 * b[i] * t[i] / h;
		Kh += D22 * b[i] * t[i] / h;

	}

	// Build the initial stiffness matrix
	double Kv = 0.0; double Km = 0.0; double e = 0.0; // double ex = 0.0;

	for (int i = 0; i < m; ++i)
	{
		Kv += ky[i];
		Km += ky[i] * x[i] * x[i];
		e += ky[i] * x[i];
	}

	// Assemble element stiffness matrix
	E_SFI_MVLEM_3DKlocal(0, 0) = Kh / 4.0 + (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(0, 1) = (Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 2) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 3) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 4) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 5) = -(Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 6) = Kh / 4.0 - (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(0, 7) = -(Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 11) = -(Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 12) = -Kh / 4.0;
	E_SFI_MVLEM_3DKlocal(0, 13) = -(Kh * d * h * (c - 1)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 17) = (Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 18) = -Kh / 4;
	E_SFI_MVLEM_3DKlocal(0, 19) = (Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 23) = (Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(1, 0) = E_SFI_MVLEM_3DKlocal(0, 1);
	E_SFI_MVLEM_3DKlocal(1, 1) = Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - (d * (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	E_SFI_MVLEM_3DKlocal(1, 2) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 3) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 4) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 5) = (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(1, 6) = (Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(1, 7) = Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) + (d * (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) - (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	E_SFI_MVLEM_3DKlocal(1, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 11) = (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(1, 12) = -(Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(1, 13) = (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - Kv / 4.0 + (d * (e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(1, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 17) = -(e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(1, 18) = -(Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(1, 19) = (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - Kv / 4.0 - (d * (e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(1, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 23) = -(e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);

	E_SFI_MVLEM_3DKlocal(2, 0) = E_SFI_MVLEM_3DKlocal(0, 2);
	E_SFI_MVLEM_3DKlocal(2, 1) = E_SFI_MVLEM_3DKlocal(1, 2);
	E_SFI_MVLEM_3DKlocal(2, 2) = K1;
	E_SFI_MVLEM_3DKlocal(2, 3) = -K2;
	E_SFI_MVLEM_3DKlocal(2, 4) = K3;
	E_SFI_MVLEM_3DKlocal(2, 5) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 6) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 7) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 8) = K4;
	E_SFI_MVLEM_3DKlocal(2, 9) = K5;
	E_SFI_MVLEM_3DKlocal(2, 10) = K6;
	E_SFI_MVLEM_3DKlocal(2, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 14) = K7;
	E_SFI_MVLEM_3DKlocal(2, 15) = -K8;
	E_SFI_MVLEM_3DKlocal(2, 16) = -K9;
	E_SFI_MVLEM_3DKlocal(2, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 20) = K10;
	E_SFI_MVLEM_3DKlocal(2, 21) = -K11;
	E_SFI_MVLEM_3DKlocal(2, 22) = K12;
	E_SFI_MVLEM_3DKlocal(2, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(3, 0) = E_SFI_MVLEM_3DKlocal(0, 3);
	E_SFI_MVLEM_3DKlocal(3, 1) = E_SFI_MVLEM_3DKlocal(1, 3);
	E_SFI_MVLEM_3DKlocal(3, 2) = E_SFI_MVLEM_3DKlocal(2, 3);
	E_SFI_MVLEM_3DKlocal(3, 3) = K13;
	E_SFI_MVLEM_3DKlocal(3, 4) = K14;
	E_SFI_MVLEM_3DKlocal(3, 5) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 6) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 7) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 8) = K5;
	E_SFI_MVLEM_3DKlocal(3, 9) = K15;
	E_SFI_MVLEM_3DKlocal(3, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 14) = K8;
	E_SFI_MVLEM_3DKlocal(3, 15) = K16;
	E_SFI_MVLEM_3DKlocal(3, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 20) = K11;
	E_SFI_MVLEM_3DKlocal(3, 21) = K17;
	E_SFI_MVLEM_3DKlocal(3, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(4, 0) = E_SFI_MVLEM_3DKlocal(0, 4);
	E_SFI_MVLEM_3DKlocal(4, 1) = E_SFI_MVLEM_3DKlocal(1, 4);
	E_SFI_MVLEM_3DKlocal(4, 2) = E_SFI_MVLEM_3DKlocal(2, 4);
	E_SFI_MVLEM_3DKlocal(4, 3) = E_SFI_MVLEM_3DKlocal(3, 4);
	E_SFI_MVLEM_3DKlocal(4, 4) = K18;
	E_SFI_MVLEM_3DKlocal(4, 5) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 6) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 7) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 8) = -K6;
	E_SFI_MVLEM_3DKlocal(4, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 10) = K19;
	E_SFI_MVLEM_3DKlocal(4, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 14) = -K9;
	E_SFI_MVLEM_3DKlocal(4, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 16) = K20;
	E_SFI_MVLEM_3DKlocal(4, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 20) = -K12;
	E_SFI_MVLEM_3DKlocal(4, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 22) = K21;
	E_SFI_MVLEM_3DKlocal(4, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(5, 0) = E_SFI_MVLEM_3DKlocal(0, 5);
	E_SFI_MVLEM_3DKlocal(5, 1) = E_SFI_MVLEM_3DKlocal(1, 5);
	E_SFI_MVLEM_3DKlocal(5, 2) = E_SFI_MVLEM_3DKlocal(2, 5);
	E_SFI_MVLEM_3DKlocal(5, 3) = E_SFI_MVLEM_3DKlocal(3, 5);
	E_SFI_MVLEM_3DKlocal(5, 4) = E_SFI_MVLEM_3DKlocal(4, 5);
	E_SFI_MVLEM_3DKlocal(5, 5) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	E_SFI_MVLEM_3DKlocal(5, 6) = -(Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 7) = e / (2.0 * (2.0 * (d * d) + 2.0)) + (d * (Km + Kh * (c * c) * (h * h))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(5, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 11) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (2.0 * Eib * Iib) / Lw;
	E_SFI_MVLEM_3DKlocal(5, 12) = (Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 13) = (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - e / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 17) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 18) = (Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 19) = -e / (2.0 * (2.0 * (d * d) + 2.0)) - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 23) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(6, 0) = E_SFI_MVLEM_3DKlocal(0, 6);
	E_SFI_MVLEM_3DKlocal(6, 1) = E_SFI_MVLEM_3DKlocal(1, 6);
	E_SFI_MVLEM_3DKlocal(6, 2) = E_SFI_MVLEM_3DKlocal(2, 6);
	E_SFI_MVLEM_3DKlocal(6, 3) = E_SFI_MVLEM_3DKlocal(3, 6);
	E_SFI_MVLEM_3DKlocal(6, 4) = E_SFI_MVLEM_3DKlocal(4, 6);
	E_SFI_MVLEM_3DKlocal(6, 5) = E_SFI_MVLEM_3DKlocal(5, 6);
	E_SFI_MVLEM_3DKlocal(6, 6) = Kh / 4.0 + (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(6, 7) = -(Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 11) = -(Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 12) = -Kh / 4;
	E_SFI_MVLEM_3DKlocal(6, 13) = -(Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 17) = (Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 18) = -Kh / 4.0;
	E_SFI_MVLEM_3DKlocal(6, 19) = (Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 23) = (Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(7, 0) = E_SFI_MVLEM_3DKlocal(0, 7);
	E_SFI_MVLEM_3DKlocal(7, 1) = E_SFI_MVLEM_3DKlocal(1, 7);
	E_SFI_MVLEM_3DKlocal(7, 2) = E_SFI_MVLEM_3DKlocal(2, 7);
	E_SFI_MVLEM_3DKlocal(7, 3) = E_SFI_MVLEM_3DKlocal(3, 7);
	E_SFI_MVLEM_3DKlocal(7, 4) = E_SFI_MVLEM_3DKlocal(4, 7);
	E_SFI_MVLEM_3DKlocal(7, 5) = E_SFI_MVLEM_3DKlocal(5, 7);
	E_SFI_MVLEM_3DKlocal(7, 6) = E_SFI_MVLEM_3DKlocal(6, 7);
	E_SFI_MVLEM_3DKlocal(7, 7) = Kv / 4.0 + (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) + (d * (e / 2.0 + (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	E_SFI_MVLEM_3DKlocal(7, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 11) = (e / 2.0 + (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) - (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(7, 12) = (Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(7, 13) = (d * (e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - Kv / 4.0;
	E_SFI_MVLEM_3DKlocal(7, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 17) = -(e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(7, 18) = (Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(7, 19) = -Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - (d * (e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(7, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 23) = -(e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);

	E_SFI_MVLEM_3DKlocal(8, 0) = E_SFI_MVLEM_3DKlocal(0, 8);
	E_SFI_MVLEM_3DKlocal(8, 1) = E_SFI_MVLEM_3DKlocal(1, 8);
	E_SFI_MVLEM_3DKlocal(8, 2) = E_SFI_MVLEM_3DKlocal(2, 8);
	E_SFI_MVLEM_3DKlocal(8, 3) = E_SFI_MVLEM_3DKlocal(3, 8);
	E_SFI_MVLEM_3DKlocal(8, 4) = E_SFI_MVLEM_3DKlocal(4, 8);
	E_SFI_MVLEM_3DKlocal(8, 5) = E_SFI_MVLEM_3DKlocal(5, 8);
	E_SFI_MVLEM_3DKlocal(8, 6) = E_SFI_MVLEM_3DKlocal(6, 8);
	E_SFI_MVLEM_3DKlocal(8, 7) = E_SFI_MVLEM_3DKlocal(7, 8);
	E_SFI_MVLEM_3DKlocal(8, 8) = K1;
	E_SFI_MVLEM_3DKlocal(8, 9) = -K2;
	E_SFI_MVLEM_3DKlocal(8, 10) = -K3;
	E_SFI_MVLEM_3DKlocal(8, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 14) = K10;
	E_SFI_MVLEM_3DKlocal(8, 15) = -K11;
	E_SFI_MVLEM_3DKlocal(8, 16) = -K12;
	E_SFI_MVLEM_3DKlocal(8, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 20) = K7;
	E_SFI_MVLEM_3DKlocal(8, 21) = -K8;
	E_SFI_MVLEM_3DKlocal(8, 22) = K9;
	E_SFI_MVLEM_3DKlocal(8, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(9, 1) = E_SFI_MVLEM_3DKlocal(1, 9);
	E_SFI_MVLEM_3DKlocal(9, 2) = E_SFI_MVLEM_3DKlocal(2, 9);
	E_SFI_MVLEM_3DKlocal(9, 3) = E_SFI_MVLEM_3DKlocal(3, 9);
	E_SFI_MVLEM_3DKlocal(9, 4) = E_SFI_MVLEM_3DKlocal(4, 9);
	E_SFI_MVLEM_3DKlocal(9, 5) = E_SFI_MVLEM_3DKlocal(5, 9);
	E_SFI_MVLEM_3DKlocal(9, 6) = E_SFI_MVLEM_3DKlocal(6, 9);
	E_SFI_MVLEM_3DKlocal(9, 7) = E_SFI_MVLEM_3DKlocal(7, 9);
	E_SFI_MVLEM_3DKlocal(9, 8) = E_SFI_MVLEM_3DKlocal(8, 9);
	E_SFI_MVLEM_3DKlocal(9, 9) = K13;
	E_SFI_MVLEM_3DKlocal(9, 10) = -K14;
	E_SFI_MVLEM_3DKlocal(9, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 14) = K11;
	E_SFI_MVLEM_3DKlocal(9, 15) = K17;
	E_SFI_MVLEM_3DKlocal(9, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 20) = K8;
	E_SFI_MVLEM_3DKlocal(9, 21) = K16;
	E_SFI_MVLEM_3DKlocal(9, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(10, 0) = E_SFI_MVLEM_3DKlocal(0, 10);
	E_SFI_MVLEM_3DKlocal(10, 1) = E_SFI_MVLEM_3DKlocal(1, 10);
	E_SFI_MVLEM_3DKlocal(10, 2) = E_SFI_MVLEM_3DKlocal(2, 10);
	E_SFI_MVLEM_3DKlocal(10, 3) = E_SFI_MVLEM_3DKlocal(3, 10);
	E_SFI_MVLEM_3DKlocal(10, 4) = E_SFI_MVLEM_3DKlocal(4, 10);
	E_SFI_MVLEM_3DKlocal(10, 5) = E_SFI_MVLEM_3DKlocal(5, 10);
	E_SFI_MVLEM_3DKlocal(10, 6) = E_SFI_MVLEM_3DKlocal(6, 10);
	E_SFI_MVLEM_3DKlocal(10, 7) = E_SFI_MVLEM_3DKlocal(7, 10);
	E_SFI_MVLEM_3DKlocal(10, 8) = E_SFI_MVLEM_3DKlocal(8, 10);
	E_SFI_MVLEM_3DKlocal(10, 9) = E_SFI_MVLEM_3DKlocal(9, 10);
	E_SFI_MVLEM_3DKlocal(10, 10) = K22;
	E_SFI_MVLEM_3DKlocal(10, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 14) = K12;
	E_SFI_MVLEM_3DKlocal(10, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 16) = K21;
	E_SFI_MVLEM_3DKlocal(10, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 20) = K9;
	E_SFI_MVLEM_3DKlocal(10, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 22) = K20;
	E_SFI_MVLEM_3DKlocal(10, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(11, 0) = E_SFI_MVLEM_3DKlocal(0, 11);
	E_SFI_MVLEM_3DKlocal(11, 1) = E_SFI_MVLEM_3DKlocal(1, 11);
	E_SFI_MVLEM_3DKlocal(11, 2) = E_SFI_MVLEM_3DKlocal(2, 11);
	E_SFI_MVLEM_3DKlocal(11, 3) = E_SFI_MVLEM_3DKlocal(3, 11);
	E_SFI_MVLEM_3DKlocal(11, 4) = E_SFI_MVLEM_3DKlocal(4, 11);
	E_SFI_MVLEM_3DKlocal(11, 5) = E_SFI_MVLEM_3DKlocal(5, 11);
	E_SFI_MVLEM_3DKlocal(11, 6) = E_SFI_MVLEM_3DKlocal(6, 11);
	E_SFI_MVLEM_3DKlocal(11, 7) = E_SFI_MVLEM_3DKlocal(7, 11);
	E_SFI_MVLEM_3DKlocal(11, 8) = E_SFI_MVLEM_3DKlocal(8, 11);
	E_SFI_MVLEM_3DKlocal(11, 9) = E_SFI_MVLEM_3DKlocal(9, 11);
	E_SFI_MVLEM_3DKlocal(11, 10) = E_SFI_MVLEM_3DKlocal(10, 11);
	E_SFI_MVLEM_3DKlocal(11, 11) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	E_SFI_MVLEM_3DKlocal(11, 12) = (Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 13) = (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - e / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 17) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 18) = (Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 19) = -e / (2.0 * (2.0 * (d * d) + 2.0)) - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 23) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(12, 0) = E_SFI_MVLEM_3DKlocal(0, 12);
	E_SFI_MVLEM_3DKlocal(12, 1) = E_SFI_MVLEM_3DKlocal(1, 12);
	E_SFI_MVLEM_3DKlocal(12, 2) = E_SFI_MVLEM_3DKlocal(2, 12);
	E_SFI_MVLEM_3DKlocal(12, 3) = E_SFI_MVLEM_3DKlocal(3, 12);
	E_SFI_MVLEM_3DKlocal(12, 4) = E_SFI_MVLEM_3DKlocal(4, 12);
	E_SFI_MVLEM_3DKlocal(12, 5) = E_SFI_MVLEM_3DKlocal(5, 12);
	E_SFI_MVLEM_3DKlocal(12, 6) = E_SFI_MVLEM_3DKlocal(6, 12);
	E_SFI_MVLEM_3DKlocal(12, 7) = E_SFI_MVLEM_3DKlocal(7, 12);
	E_SFI_MVLEM_3DKlocal(12, 8) = E_SFI_MVLEM_3DKlocal(8, 12);
	E_SFI_MVLEM_3DKlocal(12, 9) = E_SFI_MVLEM_3DKlocal(9, 12);
	E_SFI_MVLEM_3DKlocal(12, 10) = E_SFI_MVLEM_3DKlocal(10, 12);
	E_SFI_MVLEM_3DKlocal(12, 11) = E_SFI_MVLEM_3DKlocal(11, 12);
	E_SFI_MVLEM_3DKlocal(12, 12) = Kh / 4.0 + (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(12, 13) = (Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(12, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 17) = -(Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(12, 18) = Kh / 4.0 - (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(12, 19) = -(Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(12, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 23) = -(Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(13, 0) = E_SFI_MVLEM_3DKlocal(0, 13);
	E_SFI_MVLEM_3DKlocal(13, 1) = E_SFI_MVLEM_3DKlocal(1, 13);
	E_SFI_MVLEM_3DKlocal(13, 2) = E_SFI_MVLEM_3DKlocal(2, 13);
	E_SFI_MVLEM_3DKlocal(13, 3) = E_SFI_MVLEM_3DKlocal(3, 13);
	E_SFI_MVLEM_3DKlocal(13, 4) = E_SFI_MVLEM_3DKlocal(4, 13);
	E_SFI_MVLEM_3DKlocal(13, 5) = E_SFI_MVLEM_3DKlocal(5, 13);
	E_SFI_MVLEM_3DKlocal(13, 6) = E_SFI_MVLEM_3DKlocal(6, 13);
	E_SFI_MVLEM_3DKlocal(13, 7) = E_SFI_MVLEM_3DKlocal(7, 13);
	E_SFI_MVLEM_3DKlocal(13, 8) = E_SFI_MVLEM_3DKlocal(8, 13);
	E_SFI_MVLEM_3DKlocal(13, 9) = E_SFI_MVLEM_3DKlocal(9, 13);
	E_SFI_MVLEM_3DKlocal(13, 10) = E_SFI_MVLEM_3DKlocal(10, 13);
	E_SFI_MVLEM_3DKlocal(13, 11) = E_SFI_MVLEM_3DKlocal(11, 13);
	E_SFI_MVLEM_3DKlocal(13, 12) = E_SFI_MVLEM_3DKlocal(12, 13);
	E_SFI_MVLEM_3DKlocal(13, 13) = Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) + (12.0 * Eib * Iib) / (Lw * Lw * Lw) - (d * (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(13, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 17) = (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(13, 18) = (Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(13, 19) = Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - (12.0 * Eib * Iib) / (Lw * Lw * Lw) + (d * (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(13, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 23) = (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);

	E_SFI_MVLEM_3DKlocal(14, 0) = E_SFI_MVLEM_3DKlocal(0, 14);
	E_SFI_MVLEM_3DKlocal(14, 1) = E_SFI_MVLEM_3DKlocal(1, 14);
	E_SFI_MVLEM_3DKlocal(14, 2) = E_SFI_MVLEM_3DKlocal(2, 14);
	E_SFI_MVLEM_3DKlocal(14, 3) = E_SFI_MVLEM_3DKlocal(3, 14);
	E_SFI_MVLEM_3DKlocal(14, 4) = E_SFI_MVLEM_3DKlocal(4, 14);
	E_SFI_MVLEM_3DKlocal(14, 5) = E_SFI_MVLEM_3DKlocal(5, 14);
	E_SFI_MVLEM_3DKlocal(14, 6) = E_SFI_MVLEM_3DKlocal(6, 14);
	E_SFI_MVLEM_3DKlocal(14, 7) = E_SFI_MVLEM_3DKlocal(7, 14);
	E_SFI_MVLEM_3DKlocal(14, 8) = E_SFI_MVLEM_3DKlocal(8, 14);
	E_SFI_MVLEM_3DKlocal(14, 9) = E_SFI_MVLEM_3DKlocal(9, 14);
	E_SFI_MVLEM_3DKlocal(14, 10) = E_SFI_MVLEM_3DKlocal(10, 14);
	E_SFI_MVLEM_3DKlocal(14, 11) = E_SFI_MVLEM_3DKlocal(11, 14);
	E_SFI_MVLEM_3DKlocal(14, 12) = E_SFI_MVLEM_3DKlocal(12, 14);
	E_SFI_MVLEM_3DKlocal(14, 13) = E_SFI_MVLEM_3DKlocal(13, 14);
	E_SFI_MVLEM_3DKlocal(14, 14) = K1;
	E_SFI_MVLEM_3DKlocal(14, 15) = K2;
	E_SFI_MVLEM_3DKlocal(14, 16) = K3;
	E_SFI_MVLEM_3DKlocal(14, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(14, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(14, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(14, 20) = K4;
	E_SFI_MVLEM_3DKlocal(14, 21) = -K5;
	E_SFI_MVLEM_3DKlocal(14, 22) = K6;
	E_SFI_MVLEM_3DKlocal(14, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(15, 0) = E_SFI_MVLEM_3DKlocal(0, 15);
	E_SFI_MVLEM_3DKlocal(15, 1) = E_SFI_MVLEM_3DKlocal(1, 15);
	E_SFI_MVLEM_3DKlocal(15, 2) = E_SFI_MVLEM_3DKlocal(2, 15);
	E_SFI_MVLEM_3DKlocal(15, 3) = E_SFI_MVLEM_3DKlocal(3, 15);
	E_SFI_MVLEM_3DKlocal(15, 4) = E_SFI_MVLEM_3DKlocal(4, 15);
	E_SFI_MVLEM_3DKlocal(15, 5) = E_SFI_MVLEM_3DKlocal(5, 15);
	E_SFI_MVLEM_3DKlocal(15, 6) = E_SFI_MVLEM_3DKlocal(6, 15);
	E_SFI_MVLEM_3DKlocal(15, 7) = E_SFI_MVLEM_3DKlocal(7, 15);
	E_SFI_MVLEM_3DKlocal(15, 8) = E_SFI_MVLEM_3DKlocal(8, 15);
	E_SFI_MVLEM_3DKlocal(15, 9) = E_SFI_MVLEM_3DKlocal(9, 15);
	E_SFI_MVLEM_3DKlocal(15, 10) = E_SFI_MVLEM_3DKlocal(10, 15);
	E_SFI_MVLEM_3DKlocal(15, 11) = E_SFI_MVLEM_3DKlocal(11, 15);
	E_SFI_MVLEM_3DKlocal(15, 12) = E_SFI_MVLEM_3DKlocal(12, 15);
	E_SFI_MVLEM_3DKlocal(15, 13) = E_SFI_MVLEM_3DKlocal(13, 15);
	E_SFI_MVLEM_3DKlocal(15, 14) = E_SFI_MVLEM_3DKlocal(14, 15);
	E_SFI_MVLEM_3DKlocal(15, 15) = K13;
	E_SFI_MVLEM_3DKlocal(15, 16) = -K14;
	E_SFI_MVLEM_3DKlocal(15, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(15, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(15, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(15, 20) = -K5;
	E_SFI_MVLEM_3DKlocal(15, 21) = K15;
	E_SFI_MVLEM_3DKlocal(15, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(15, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(16, 0) = E_SFI_MVLEM_3DKlocal(0, 16);
	E_SFI_MVLEM_3DKlocal(16, 1) = E_SFI_MVLEM_3DKlocal(1, 16);
	E_SFI_MVLEM_3DKlocal(16, 2) = E_SFI_MVLEM_3DKlocal(2, 16);
	E_SFI_MVLEM_3DKlocal(16, 3) = E_SFI_MVLEM_3DKlocal(3, 16);
	E_SFI_MVLEM_3DKlocal(16, 4) = E_SFI_MVLEM_3DKlocal(4, 16);
	E_SFI_MVLEM_3DKlocal(16, 5) = E_SFI_MVLEM_3DKlocal(5, 16);
	E_SFI_MVLEM_3DKlocal(16, 6) = E_SFI_MVLEM_3DKlocal(6, 16);
	E_SFI_MVLEM_3DKlocal(16, 7) = E_SFI_MVLEM_3DKlocal(7, 16);
	E_SFI_MVLEM_3DKlocal(16, 8) = E_SFI_MVLEM_3DKlocal(8, 16);
	E_SFI_MVLEM_3DKlocal(16, 9) = E_SFI_MVLEM_3DKlocal(9, 16);
	E_SFI_MVLEM_3DKlocal(16, 10) = E_SFI_MVLEM_3DKlocal(10, 16);
	E_SFI_MVLEM_3DKlocal(16, 11) = E_SFI_MVLEM_3DKlocal(11, 16);
	E_SFI_MVLEM_3DKlocal(16, 12) = E_SFI_MVLEM_3DKlocal(12, 16);
	E_SFI_MVLEM_3DKlocal(16, 13) = E_SFI_MVLEM_3DKlocal(13, 16);
	E_SFI_MVLEM_3DKlocal(16, 14) = E_SFI_MVLEM_3DKlocal(14, 16);
	E_SFI_MVLEM_3DKlocal(16, 15) = E_SFI_MVLEM_3DKlocal(15, 16);
	E_SFI_MVLEM_3DKlocal(16, 16) = K18;
	E_SFI_MVLEM_3DKlocal(16, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(16, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(16, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(16, 20) = -K6;
	E_SFI_MVLEM_3DKlocal(16, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(16, 22) = K19;
	E_SFI_MVLEM_3DKlocal(16, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(17, 0) = E_SFI_MVLEM_3DKlocal(0, 17);
	E_SFI_MVLEM_3DKlocal(17, 1) = E_SFI_MVLEM_3DKlocal(1, 17);
	E_SFI_MVLEM_3DKlocal(17, 2) = E_SFI_MVLEM_3DKlocal(2, 17);
	E_SFI_MVLEM_3DKlocal(17, 3) = E_SFI_MVLEM_3DKlocal(3, 17);
	E_SFI_MVLEM_3DKlocal(17, 4) = E_SFI_MVLEM_3DKlocal(4, 17);
	E_SFI_MVLEM_3DKlocal(17, 5) = E_SFI_MVLEM_3DKlocal(5, 17);
	E_SFI_MVLEM_3DKlocal(17, 6) = E_SFI_MVLEM_3DKlocal(6, 17);
	E_SFI_MVLEM_3DKlocal(17, 7) = E_SFI_MVLEM_3DKlocal(7, 17);
	E_SFI_MVLEM_3DKlocal(17, 8) = E_SFI_MVLEM_3DKlocal(8, 17);
	E_SFI_MVLEM_3DKlocal(17, 9) = E_SFI_MVLEM_3DKlocal(9, 17);
	E_SFI_MVLEM_3DKlocal(17, 10) = E_SFI_MVLEM_3DKlocal(10, 17);
	E_SFI_MVLEM_3DKlocal(17, 11) = E_SFI_MVLEM_3DKlocal(11, 17);
	E_SFI_MVLEM_3DKlocal(17, 12) = E_SFI_MVLEM_3DKlocal(12, 17);
	E_SFI_MVLEM_3DKlocal(17, 13) = E_SFI_MVLEM_3DKlocal(13, 17);
	E_SFI_MVLEM_3DKlocal(17, 14) = E_SFI_MVLEM_3DKlocal(14, 17);
	E_SFI_MVLEM_3DKlocal(17, 15) = E_SFI_MVLEM_3DKlocal(15, 17);
	E_SFI_MVLEM_3DKlocal(17, 16) = E_SFI_MVLEM_3DKlocal(16, 17);
	E_SFI_MVLEM_3DKlocal(17, 17) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	E_SFI_MVLEM_3DKlocal(17, 18) = -(Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(17, 19) = e / (2.0 * (2.0 * (d * d) + 2.0)) - (6.0 * Eib * Iib) / (Lw * Lw) + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(17, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(17, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(17, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(17, 23) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (2.0 * Eib * Iib) / Lw;

	E_SFI_MVLEM_3DKlocal(18, 0) = E_SFI_MVLEM_3DKlocal(0, 18);
	E_SFI_MVLEM_3DKlocal(18, 1) = E_SFI_MVLEM_3DKlocal(1, 18);
	E_SFI_MVLEM_3DKlocal(18, 2) = E_SFI_MVLEM_3DKlocal(2, 18);
	E_SFI_MVLEM_3DKlocal(18, 3) = E_SFI_MVLEM_3DKlocal(3, 18);
	E_SFI_MVLEM_3DKlocal(18, 4) = E_SFI_MVLEM_3DKlocal(4, 18);
	E_SFI_MVLEM_3DKlocal(18, 5) = E_SFI_MVLEM_3DKlocal(5, 18);
	E_SFI_MVLEM_3DKlocal(18, 6) = E_SFI_MVLEM_3DKlocal(6, 18);
	E_SFI_MVLEM_3DKlocal(18, 7) = E_SFI_MVLEM_3DKlocal(7, 18);
	E_SFI_MVLEM_3DKlocal(18, 8) = E_SFI_MVLEM_3DKlocal(8, 18);
	E_SFI_MVLEM_3DKlocal(18, 9) = E_SFI_MVLEM_3DKlocal(9, 18);
	E_SFI_MVLEM_3DKlocal(18, 10) = E_SFI_MVLEM_3DKlocal(10, 18);
	E_SFI_MVLEM_3DKlocal(18, 11) = E_SFI_MVLEM_3DKlocal(11, 18);
	E_SFI_MVLEM_3DKlocal(18, 12) = E_SFI_MVLEM_3DKlocal(12, 18);
	E_SFI_MVLEM_3DKlocal(18, 13) = E_SFI_MVLEM_3DKlocal(13, 18);
	E_SFI_MVLEM_3DKlocal(18, 14) = E_SFI_MVLEM_3DKlocal(14, 18);
	E_SFI_MVLEM_3DKlocal(18, 15) = E_SFI_MVLEM_3DKlocal(15, 18);
	E_SFI_MVLEM_3DKlocal(18, 16) = E_SFI_MVLEM_3DKlocal(16, 18);
	E_SFI_MVLEM_3DKlocal(18, 17) = E_SFI_MVLEM_3DKlocal(17, 18);
	E_SFI_MVLEM_3DKlocal(18, 18) = Kh / 4.0 + (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(18, 19) = -(Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(18, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(18, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(18, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(18, 23) = -(Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(19, 0) = E_SFI_MVLEM_3DKlocal(0, 19);
	E_SFI_MVLEM_3DKlocal(19, 1) = E_SFI_MVLEM_3DKlocal(1, 19);
	E_SFI_MVLEM_3DKlocal(19, 2) = E_SFI_MVLEM_3DKlocal(2, 19);
	E_SFI_MVLEM_3DKlocal(19, 3) = E_SFI_MVLEM_3DKlocal(3, 19);
	E_SFI_MVLEM_3DKlocal(19, 4) = E_SFI_MVLEM_3DKlocal(4, 19);
	E_SFI_MVLEM_3DKlocal(19, 5) = E_SFI_MVLEM_3DKlocal(5, 19);
	E_SFI_MVLEM_3DKlocal(19, 6) = E_SFI_MVLEM_3DKlocal(6, 19);
	E_SFI_MVLEM_3DKlocal(19, 7) = E_SFI_MVLEM_3DKlocal(7, 19);
	E_SFI_MVLEM_3DKlocal(19, 8) = E_SFI_MVLEM_3DKlocal(8, 19);
	E_SFI_MVLEM_3DKlocal(19, 9) = E_SFI_MVLEM_3DKlocal(9, 19);
	E_SFI_MVLEM_3DKlocal(19, 10) = E_SFI_MVLEM_3DKlocal(10, 19);
	E_SFI_MVLEM_3DKlocal(19, 11) = E_SFI_MVLEM_3DKlocal(11, 19);
	E_SFI_MVLEM_3DKlocal(19, 12) = E_SFI_MVLEM_3DKlocal(12, 19);
	E_SFI_MVLEM_3DKlocal(19, 13) = E_SFI_MVLEM_3DKlocal(13, 19);
	E_SFI_MVLEM_3DKlocal(19, 14) = E_SFI_MVLEM_3DKlocal(14, 19);
	E_SFI_MVLEM_3DKlocal(19, 15) = E_SFI_MVLEM_3DKlocal(15, 19);
	E_SFI_MVLEM_3DKlocal(19, 16) = E_SFI_MVLEM_3DKlocal(16, 19);
	E_SFI_MVLEM_3DKlocal(19, 17) = E_SFI_MVLEM_3DKlocal(17, 19);
	E_SFI_MVLEM_3DKlocal(19, 18) = E_SFI_MVLEM_3DKlocal(18, 19);
	E_SFI_MVLEM_3DKlocal(19, 19) = Kv / 4.0 + (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) + (12.0 * Eib * Iib) / (Lw * Lw * Lw) + (d * (e / 2.0 + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(19, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(19, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(19, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(19, 23) = (e / 2.0 + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) - (6.0 * Eib * Iib) / (Lw * Lw);

	E_SFI_MVLEM_3DKlocal(20, 0) = E_SFI_MVLEM_3DKlocal(0, 20);
	E_SFI_MVLEM_3DKlocal(20, 1) = E_SFI_MVLEM_3DKlocal(1, 20);
	E_SFI_MVLEM_3DKlocal(20, 2) = E_SFI_MVLEM_3DKlocal(2, 20);
	E_SFI_MVLEM_3DKlocal(20, 3) = E_SFI_MVLEM_3DKlocal(3, 20);
	E_SFI_MVLEM_3DKlocal(20, 4) = E_SFI_MVLEM_3DKlocal(4, 20);
	E_SFI_MVLEM_3DKlocal(20, 5) = E_SFI_MVLEM_3DKlocal(5, 20);
	E_SFI_MVLEM_3DKlocal(20, 6) = E_SFI_MVLEM_3DKlocal(6, 20);
	E_SFI_MVLEM_3DKlocal(20, 7) = E_SFI_MVLEM_3DKlocal(7, 20);
	E_SFI_MVLEM_3DKlocal(20, 8) = E_SFI_MVLEM_3DKlocal(8, 20);
	E_SFI_MVLEM_3DKlocal(20, 9) = E_SFI_MVLEM_3DKlocal(9, 20);
	E_SFI_MVLEM_3DKlocal(20, 10) = E_SFI_MVLEM_3DKlocal(10, 20);
	E_SFI_MVLEM_3DKlocal(20, 11) = E_SFI_MVLEM_3DKlocal(11, 20);
	E_SFI_MVLEM_3DKlocal(20, 12) = E_SFI_MVLEM_3DKlocal(12, 20);
	E_SFI_MVLEM_3DKlocal(20, 13) = E_SFI_MVLEM_3DKlocal(13, 20);
	E_SFI_MVLEM_3DKlocal(20, 14) = E_SFI_MVLEM_3DKlocal(14, 20);
	E_SFI_MVLEM_3DKlocal(20, 15) = E_SFI_MVLEM_3DKlocal(15, 20);
	E_SFI_MVLEM_3DKlocal(20, 16) = E_SFI_MVLEM_3DKlocal(16, 20);
	E_SFI_MVLEM_3DKlocal(20, 17) = E_SFI_MVLEM_3DKlocal(17, 20);
	E_SFI_MVLEM_3DKlocal(20, 18) = E_SFI_MVLEM_3DKlocal(18, 20);
	E_SFI_MVLEM_3DKlocal(20, 19) = E_SFI_MVLEM_3DKlocal(19, 20);
	E_SFI_MVLEM_3DKlocal(20, 20) = K1;
	E_SFI_MVLEM_3DKlocal(20, 21) = K2;
	E_SFI_MVLEM_3DKlocal(20, 22) = -K3;
	E_SFI_MVLEM_3DKlocal(20, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(21, 0) = E_SFI_MVLEM_3DKlocal(0, 21);
	E_SFI_MVLEM_3DKlocal(21, 1) = E_SFI_MVLEM_3DKlocal(1, 21);
	E_SFI_MVLEM_3DKlocal(21, 2) = E_SFI_MVLEM_3DKlocal(2, 21);
	E_SFI_MVLEM_3DKlocal(21, 3) = E_SFI_MVLEM_3DKlocal(3, 21);
	E_SFI_MVLEM_3DKlocal(21, 4) = E_SFI_MVLEM_3DKlocal(4, 21);
	E_SFI_MVLEM_3DKlocal(21, 5) = E_SFI_MVLEM_3DKlocal(5, 21);
	E_SFI_MVLEM_3DKlocal(21, 6) = E_SFI_MVLEM_3DKlocal(6, 21);
	E_SFI_MVLEM_3DKlocal(21, 7) = E_SFI_MVLEM_3DKlocal(7, 21);
	E_SFI_MVLEM_3DKlocal(21, 8) = E_SFI_MVLEM_3DKlocal(8, 21);
	E_SFI_MVLEM_3DKlocal(21, 9) = E_SFI_MVLEM_3DKlocal(9, 21);
	E_SFI_MVLEM_3DKlocal(21, 10) = E_SFI_MVLEM_3DKlocal(10, 21);
	E_SFI_MVLEM_3DKlocal(21, 11) = E_SFI_MVLEM_3DKlocal(11, 21);
	E_SFI_MVLEM_3DKlocal(21, 12) = E_SFI_MVLEM_3DKlocal(12, 21);
	E_SFI_MVLEM_3DKlocal(21, 13) = E_SFI_MVLEM_3DKlocal(13, 21);
	E_SFI_MVLEM_3DKlocal(21, 14) = E_SFI_MVLEM_3DKlocal(14, 21);
	E_SFI_MVLEM_3DKlocal(21, 15) = E_SFI_MVLEM_3DKlocal(15, 21);
	E_SFI_MVLEM_3DKlocal(21, 16) = E_SFI_MVLEM_3DKlocal(16, 21);
	E_SFI_MVLEM_3DKlocal(21, 17) = E_SFI_MVLEM_3DKlocal(17, 21);
	E_SFI_MVLEM_3DKlocal(21, 18) = E_SFI_MVLEM_3DKlocal(18, 21);
	E_SFI_MVLEM_3DKlocal(21, 19) = E_SFI_MVLEM_3DKlocal(19, 21);
	E_SFI_MVLEM_3DKlocal(21, 20) = E_SFI_MVLEM_3DKlocal(20, 21);
	E_SFI_MVLEM_3DKlocal(21, 21) = K13;
	E_SFI_MVLEM_3DKlocal(21, 22) = K14;
	E_SFI_MVLEM_3DKlocal(21, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(22, 0) = E_SFI_MVLEM_3DKlocal(0, 22);
	E_SFI_MVLEM_3DKlocal(22, 1) = E_SFI_MVLEM_3DKlocal(1, 22);
	E_SFI_MVLEM_3DKlocal(22, 2) = E_SFI_MVLEM_3DKlocal(2, 22);
	E_SFI_MVLEM_3DKlocal(22, 3) = E_SFI_MVLEM_3DKlocal(3, 22);
	E_SFI_MVLEM_3DKlocal(22, 4) = E_SFI_MVLEM_3DKlocal(4, 22);
	E_SFI_MVLEM_3DKlocal(22, 5) = E_SFI_MVLEM_3DKlocal(5, 22);
	E_SFI_MVLEM_3DKlocal(22, 6) = E_SFI_MVLEM_3DKlocal(6, 22);
	E_SFI_MVLEM_3DKlocal(22, 7) = E_SFI_MVLEM_3DKlocal(7, 22);
	E_SFI_MVLEM_3DKlocal(22, 8) = E_SFI_MVLEM_3DKlocal(8, 22);
	E_SFI_MVLEM_3DKlocal(22, 9) = E_SFI_MVLEM_3DKlocal(9, 22);
	E_SFI_MVLEM_3DKlocal(22, 10) = E_SFI_MVLEM_3DKlocal(10, 22);
	E_SFI_MVLEM_3DKlocal(22, 11) = E_SFI_MVLEM_3DKlocal(11, 22);
	E_SFI_MVLEM_3DKlocal(22, 12) = E_SFI_MVLEM_3DKlocal(12, 22);
	E_SFI_MVLEM_3DKlocal(22, 13) = E_SFI_MVLEM_3DKlocal(13, 22);
	E_SFI_MVLEM_3DKlocal(22, 14) = E_SFI_MVLEM_3DKlocal(14, 22);
	E_SFI_MVLEM_3DKlocal(22, 15) = E_SFI_MVLEM_3DKlocal(15, 22);
	E_SFI_MVLEM_3DKlocal(22, 16) = E_SFI_MVLEM_3DKlocal(16, 22);
	E_SFI_MVLEM_3DKlocal(22, 17) = E_SFI_MVLEM_3DKlocal(17, 22);
	E_SFI_MVLEM_3DKlocal(22, 18) = E_SFI_MVLEM_3DKlocal(18, 22);
	E_SFI_MVLEM_3DKlocal(22, 19) = E_SFI_MVLEM_3DKlocal(19, 22);
	E_SFI_MVLEM_3DKlocal(22, 20) = E_SFI_MVLEM_3DKlocal(20, 22);
	E_SFI_MVLEM_3DKlocal(22, 21) = E_SFI_MVLEM_3DKlocal(21, 22);
	E_SFI_MVLEM_3DKlocal(22, 22) = K18;
	E_SFI_MVLEM_3DKlocal(22, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(23, 0) = E_SFI_MVLEM_3DKlocal(0, 23);
	E_SFI_MVLEM_3DKlocal(23, 1) = E_SFI_MVLEM_3DKlocal(1, 23);
	E_SFI_MVLEM_3DKlocal(23, 2) = E_SFI_MVLEM_3DKlocal(2, 23);
	E_SFI_MVLEM_3DKlocal(23, 3) = E_SFI_MVLEM_3DKlocal(3, 23);
	E_SFI_MVLEM_3DKlocal(23, 4) = E_SFI_MVLEM_3DKlocal(4, 23);
	E_SFI_MVLEM_3DKlocal(23, 5) = E_SFI_MVLEM_3DKlocal(5, 23);
	E_SFI_MVLEM_3DKlocal(23, 6) = E_SFI_MVLEM_3DKlocal(6, 23);
	E_SFI_MVLEM_3DKlocal(23, 7) = E_SFI_MVLEM_3DKlocal(7, 23);
	E_SFI_MVLEM_3DKlocal(23, 8) = E_SFI_MVLEM_3DKlocal(8, 23);
	E_SFI_MVLEM_3DKlocal(23, 9) = E_SFI_MVLEM_3DKlocal(9, 23);
	E_SFI_MVLEM_3DKlocal(23, 10) = E_SFI_MVLEM_3DKlocal(10, 23);
	E_SFI_MVLEM_3DKlocal(23, 11) = E_SFI_MVLEM_3DKlocal(11, 23);
	E_SFI_MVLEM_3DKlocal(23, 12) = E_SFI_MVLEM_3DKlocal(12, 23);
	E_SFI_MVLEM_3DKlocal(23, 13) = E_SFI_MVLEM_3DKlocal(13, 23);
	E_SFI_MVLEM_3DKlocal(23, 14) = E_SFI_MVLEM_3DKlocal(14, 23);
	E_SFI_MVLEM_3DKlocal(23, 15) = E_SFI_MVLEM_3DKlocal(15, 23);
	E_SFI_MVLEM_3DKlocal(23, 16) = E_SFI_MVLEM_3DKlocal(16, 23);
	E_SFI_MVLEM_3DKlocal(23, 17) = E_SFI_MVLEM_3DKlocal(17, 23);
	E_SFI_MVLEM_3DKlocal(23, 18) = E_SFI_MVLEM_3DKlocal(18, 23);
	E_SFI_MVLEM_3DKlocal(23, 19) = E_SFI_MVLEM_3DKlocal(19, 23);
	E_SFI_MVLEM_3DKlocal(23, 20) = E_SFI_MVLEM_3DKlocal(20, 23);
	E_SFI_MVLEM_3DKlocal(23, 21) = E_SFI_MVLEM_3DKlocal(21, 23);
	E_SFI_MVLEM_3DKlocal(23, 22) = E_SFI_MVLEM_3DKlocal(22, 23);
	E_SFI_MVLEM_3DKlocal(23, 23) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;

	E_SFI_MVLEM_3DK.addMatrixTripleProduct(0.0, T, E_SFI_MVLEM_3DKlocal, 1.0);  // Convert matrix from local to global cs

	// Return element stiffness matrix
	return E_SFI_MVLEM_3DK;

}

// Get current element tangent stiffness matrix from the material for the last updated strain
const Matrix& E_SFI_MVLEM_3D::getTangentStiff(void)
{

	E_SFI_MVLEM_3DK.Zero();		// Global stiffness matrix
	E_SFI_MVLEM_3DKlocal.Zero();	// Local stiffness matrix

	Kh = 0.0;

	for (int i = 0; i < m; i++)
	{
		// Get the material tangent
		const Matrix& D = theMaterial[i]->getTangent();

		double D00 = D(0, 0); double D01 = D(0, 1); double D02 = D(0, 2);
		double D10 = D(1, 0); double D11 = D(1, 1); double D12 = D(1, 2);
		double D20 = D(2, 0); double D21 = D(2, 1); double D22 = D(2, 2);

		kx[i] = D00 * h * t[i] / b[i];
		ky[i] = D11 * b[i] * t[i] / h;
		Kh += D22 * b[i] * t[i] / h;

	}

	// Build the tangent stiffness matrix
	double Kv = 0.0; double Km = 0.0; double e = 0.0; // double ex = 0.0;

	for (int i = 0; i < m; ++i)
	{
		Kv += ky[i];
		Km += ky[i] * x[i] * x[i];
		e += ky[i] * x[i];

	}

	// Assemble element stiffness matrix
	E_SFI_MVLEM_3DKlocal(0, 0) = Kh / 4.0 + (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(0, 1) = (Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 2) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 3) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 4) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 5) = -(Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 6) = Kh / 4.0 - (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(0, 7) = -(Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 11) = -(Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 12) = -Kh / 4.0;
	E_SFI_MVLEM_3DKlocal(0, 13) = -(Kh * d * h * (c - 1)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 17) = (Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 18) = -Kh / 4;
	E_SFI_MVLEM_3DKlocal(0, 19) = (Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(0, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(0, 23) = (Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(1, 0) = E_SFI_MVLEM_3DKlocal(0, 1);
	E_SFI_MVLEM_3DKlocal(1, 1) = Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - (d * (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	E_SFI_MVLEM_3DKlocal(1, 2) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 3) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 4) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 5) = (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(1, 6) = (Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(1, 7) = Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) + (d * (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) - (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	E_SFI_MVLEM_3DKlocal(1, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 11) = (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(1, 12) = -(Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(1, 13) = (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - Kv / 4.0 + (d * (e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(1, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 17) = -(e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(1, 18) = -(Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(1, 19) = (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - Kv / 4.0 - (d * (e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(1, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(1, 23) = -(e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);

	E_SFI_MVLEM_3DKlocal(2, 0) = E_SFI_MVLEM_3DKlocal(0, 2);
	E_SFI_MVLEM_3DKlocal(2, 1) = E_SFI_MVLEM_3DKlocal(1, 2);
	E_SFI_MVLEM_3DKlocal(2, 2) = K1;
	E_SFI_MVLEM_3DKlocal(2, 3) = -K2;
	E_SFI_MVLEM_3DKlocal(2, 4) = K3;
	E_SFI_MVLEM_3DKlocal(2, 5) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 6) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 7) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 8) = K4;
	E_SFI_MVLEM_3DKlocal(2, 9) = K5;
	E_SFI_MVLEM_3DKlocal(2, 10) = K6;
	E_SFI_MVLEM_3DKlocal(2, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 14) = K7;
	E_SFI_MVLEM_3DKlocal(2, 15) = -K8;
	E_SFI_MVLEM_3DKlocal(2, 16) = -K9;
	E_SFI_MVLEM_3DKlocal(2, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(2, 20) = K10;
	E_SFI_MVLEM_3DKlocal(2, 21) = -K11;
	E_SFI_MVLEM_3DKlocal(2, 22) = K12;
	E_SFI_MVLEM_3DKlocal(2, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(3, 0) = E_SFI_MVLEM_3DKlocal(0, 3);
	E_SFI_MVLEM_3DKlocal(3, 1) = E_SFI_MVLEM_3DKlocal(1, 3);
	E_SFI_MVLEM_3DKlocal(3, 2) = E_SFI_MVLEM_3DKlocal(2, 3);
	E_SFI_MVLEM_3DKlocal(3, 3) = K13;
	E_SFI_MVLEM_3DKlocal(3, 4) = K14;
	E_SFI_MVLEM_3DKlocal(3, 5) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 6) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 7) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 8) = K5;
	E_SFI_MVLEM_3DKlocal(3, 9) = K15;
	E_SFI_MVLEM_3DKlocal(3, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 14) = K8;
	E_SFI_MVLEM_3DKlocal(3, 15) = K16;
	E_SFI_MVLEM_3DKlocal(3, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 20) = K11;
	E_SFI_MVLEM_3DKlocal(3, 21) = K17;
	E_SFI_MVLEM_3DKlocal(3, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(3, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(4, 0) = E_SFI_MVLEM_3DKlocal(0, 4);
	E_SFI_MVLEM_3DKlocal(4, 1) = E_SFI_MVLEM_3DKlocal(1, 4);
	E_SFI_MVLEM_3DKlocal(4, 2) = E_SFI_MVLEM_3DKlocal(2, 4);
	E_SFI_MVLEM_3DKlocal(4, 3) = E_SFI_MVLEM_3DKlocal(3, 4);
	E_SFI_MVLEM_3DKlocal(4, 4) = K18;
	E_SFI_MVLEM_3DKlocal(4, 5) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 6) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 7) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 8) = -K6;
	E_SFI_MVLEM_3DKlocal(4, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 10) = K19;
	E_SFI_MVLEM_3DKlocal(4, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 14) = -K9;
	E_SFI_MVLEM_3DKlocal(4, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 16) = K20;
	E_SFI_MVLEM_3DKlocal(4, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 20) = -K12;
	E_SFI_MVLEM_3DKlocal(4, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(4, 22) = K21;
	E_SFI_MVLEM_3DKlocal(4, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(5, 0) = E_SFI_MVLEM_3DKlocal(0, 5);
	E_SFI_MVLEM_3DKlocal(5, 1) = E_SFI_MVLEM_3DKlocal(1, 5);
	E_SFI_MVLEM_3DKlocal(5, 2) = E_SFI_MVLEM_3DKlocal(2, 5);
	E_SFI_MVLEM_3DKlocal(5, 3) = E_SFI_MVLEM_3DKlocal(3, 5);
	E_SFI_MVLEM_3DKlocal(5, 4) = E_SFI_MVLEM_3DKlocal(4, 5);
	E_SFI_MVLEM_3DKlocal(5, 5) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	E_SFI_MVLEM_3DKlocal(5, 6) = -(Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 7) = e / (2.0 * (2.0 * (d * d) + 2.0)) + (d * (Km + Kh * (c * c) * (h * h))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(5, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 11) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (2.0 * Eib * Iib) / Lw;
	E_SFI_MVLEM_3DKlocal(5, 12) = (Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 13) = (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - e / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 17) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 18) = (Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 19) = -e / (2.0 * (2.0 * (d * d) + 2.0)) - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(5, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(5, 23) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(6, 0) = E_SFI_MVLEM_3DKlocal(0, 6);
	E_SFI_MVLEM_3DKlocal(6, 1) = E_SFI_MVLEM_3DKlocal(1, 6);
	E_SFI_MVLEM_3DKlocal(6, 2) = E_SFI_MVLEM_3DKlocal(2, 6);
	E_SFI_MVLEM_3DKlocal(6, 3) = E_SFI_MVLEM_3DKlocal(3, 6);
	E_SFI_MVLEM_3DKlocal(6, 4) = E_SFI_MVLEM_3DKlocal(4, 6);
	E_SFI_MVLEM_3DKlocal(6, 5) = E_SFI_MVLEM_3DKlocal(5, 6);
	E_SFI_MVLEM_3DKlocal(6, 6) = Kh / 4.0 + (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(6, 7) = -(Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 11) = -(Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 12) = -Kh / 4;
	E_SFI_MVLEM_3DKlocal(6, 13) = -(Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 17) = (Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 18) = -Kh / 4.0;
	E_SFI_MVLEM_3DKlocal(6, 19) = (Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(6, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(6, 23) = (Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(7, 0) = E_SFI_MVLEM_3DKlocal(0, 7);
	E_SFI_MVLEM_3DKlocal(7, 1) = E_SFI_MVLEM_3DKlocal(1, 7);
	E_SFI_MVLEM_3DKlocal(7, 2) = E_SFI_MVLEM_3DKlocal(2, 7);
	E_SFI_MVLEM_3DKlocal(7, 3) = E_SFI_MVLEM_3DKlocal(3, 7);
	E_SFI_MVLEM_3DKlocal(7, 4) = E_SFI_MVLEM_3DKlocal(4, 7);
	E_SFI_MVLEM_3DKlocal(7, 5) = E_SFI_MVLEM_3DKlocal(5, 7);
	E_SFI_MVLEM_3DKlocal(7, 6) = E_SFI_MVLEM_3DKlocal(6, 7);
	E_SFI_MVLEM_3DKlocal(7, 7) = Kv / 4.0 + (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) + (d * (e / 2.0 + (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	E_SFI_MVLEM_3DKlocal(7, 8) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 9) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 10) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 11) = (e / 2.0 + (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) - (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(7, 12) = (Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(7, 13) = (d * (e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - Kv / 4.0;
	E_SFI_MVLEM_3DKlocal(7, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 17) = -(e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(7, 18) = (Kh * c * d * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(7, 19) = -Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - (d * (e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(7, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(7, 23) = -(e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);

	E_SFI_MVLEM_3DKlocal(8, 0) = E_SFI_MVLEM_3DKlocal(0, 8);
	E_SFI_MVLEM_3DKlocal(8, 1) = E_SFI_MVLEM_3DKlocal(1, 8);
	E_SFI_MVLEM_3DKlocal(8, 2) = E_SFI_MVLEM_3DKlocal(2, 8);
	E_SFI_MVLEM_3DKlocal(8, 3) = E_SFI_MVLEM_3DKlocal(3, 8);
	E_SFI_MVLEM_3DKlocal(8, 4) = E_SFI_MVLEM_3DKlocal(4, 8);
	E_SFI_MVLEM_3DKlocal(8, 5) = E_SFI_MVLEM_3DKlocal(5, 8);
	E_SFI_MVLEM_3DKlocal(8, 6) = E_SFI_MVLEM_3DKlocal(6, 8);
	E_SFI_MVLEM_3DKlocal(8, 7) = E_SFI_MVLEM_3DKlocal(7, 8);
	E_SFI_MVLEM_3DKlocal(8, 8) = K1;
	E_SFI_MVLEM_3DKlocal(8, 9) = -K2;
	E_SFI_MVLEM_3DKlocal(8, 10) = -K3;
	E_SFI_MVLEM_3DKlocal(8, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 14) = K10;
	E_SFI_MVLEM_3DKlocal(8, 15) = -K11;
	E_SFI_MVLEM_3DKlocal(8, 16) = -K12;
	E_SFI_MVLEM_3DKlocal(8, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(8, 20) = K7;
	E_SFI_MVLEM_3DKlocal(8, 21) = -K8;
	E_SFI_MVLEM_3DKlocal(8, 22) = K9;
	E_SFI_MVLEM_3DKlocal(8, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(9, 1) = E_SFI_MVLEM_3DKlocal(1, 9);
	E_SFI_MVLEM_3DKlocal(9, 2) = E_SFI_MVLEM_3DKlocal(2, 9);
	E_SFI_MVLEM_3DKlocal(9, 3) = E_SFI_MVLEM_3DKlocal(3, 9);
	E_SFI_MVLEM_3DKlocal(9, 4) = E_SFI_MVLEM_3DKlocal(4, 9);
	E_SFI_MVLEM_3DKlocal(9, 5) = E_SFI_MVLEM_3DKlocal(5, 9);
	E_SFI_MVLEM_3DKlocal(9, 6) = E_SFI_MVLEM_3DKlocal(6, 9);
	E_SFI_MVLEM_3DKlocal(9, 7) = E_SFI_MVLEM_3DKlocal(7, 9);
	E_SFI_MVLEM_3DKlocal(9, 8) = E_SFI_MVLEM_3DKlocal(8, 9);
	E_SFI_MVLEM_3DKlocal(9, 9) = K13;
	E_SFI_MVLEM_3DKlocal(9, 10) = -K14;
	E_SFI_MVLEM_3DKlocal(9, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 14) = K11;
	E_SFI_MVLEM_3DKlocal(9, 15) = K17;
	E_SFI_MVLEM_3DKlocal(9, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 20) = K8;
	E_SFI_MVLEM_3DKlocal(9, 21) = K16;
	E_SFI_MVLEM_3DKlocal(9, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(9, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(10, 0) = E_SFI_MVLEM_3DKlocal(0, 10);
	E_SFI_MVLEM_3DKlocal(10, 1) = E_SFI_MVLEM_3DKlocal(1, 10);
	E_SFI_MVLEM_3DKlocal(10, 2) = E_SFI_MVLEM_3DKlocal(2, 10);
	E_SFI_MVLEM_3DKlocal(10, 3) = E_SFI_MVLEM_3DKlocal(3, 10);
	E_SFI_MVLEM_3DKlocal(10, 4) = E_SFI_MVLEM_3DKlocal(4, 10);
	E_SFI_MVLEM_3DKlocal(10, 5) = E_SFI_MVLEM_3DKlocal(5, 10);
	E_SFI_MVLEM_3DKlocal(10, 6) = E_SFI_MVLEM_3DKlocal(6, 10);
	E_SFI_MVLEM_3DKlocal(10, 7) = E_SFI_MVLEM_3DKlocal(7, 10);
	E_SFI_MVLEM_3DKlocal(10, 8) = E_SFI_MVLEM_3DKlocal(8, 10);
	E_SFI_MVLEM_3DKlocal(10, 9) = E_SFI_MVLEM_3DKlocal(9, 10);
	E_SFI_MVLEM_3DKlocal(10, 10) = K22;
	E_SFI_MVLEM_3DKlocal(10, 11) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 12) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 13) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 14) = K12;
	E_SFI_MVLEM_3DKlocal(10, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 16) = K21;
	E_SFI_MVLEM_3DKlocal(10, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 20) = K9;
	E_SFI_MVLEM_3DKlocal(10, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(10, 22) = K20;
	E_SFI_MVLEM_3DKlocal(10, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(11, 0) = E_SFI_MVLEM_3DKlocal(0, 11);
	E_SFI_MVLEM_3DKlocal(11, 1) = E_SFI_MVLEM_3DKlocal(1, 11);
	E_SFI_MVLEM_3DKlocal(11, 2) = E_SFI_MVLEM_3DKlocal(2, 11);
	E_SFI_MVLEM_3DKlocal(11, 3) = E_SFI_MVLEM_3DKlocal(3, 11);
	E_SFI_MVLEM_3DKlocal(11, 4) = E_SFI_MVLEM_3DKlocal(4, 11);
	E_SFI_MVLEM_3DKlocal(11, 5) = E_SFI_MVLEM_3DKlocal(5, 11);
	E_SFI_MVLEM_3DKlocal(11, 6) = E_SFI_MVLEM_3DKlocal(6, 11);
	E_SFI_MVLEM_3DKlocal(11, 7) = E_SFI_MVLEM_3DKlocal(7, 11);
	E_SFI_MVLEM_3DKlocal(11, 8) = E_SFI_MVLEM_3DKlocal(8, 11);
	E_SFI_MVLEM_3DKlocal(11, 9) = E_SFI_MVLEM_3DKlocal(9, 11);
	E_SFI_MVLEM_3DKlocal(11, 10) = E_SFI_MVLEM_3DKlocal(10, 11);
	E_SFI_MVLEM_3DKlocal(11, 11) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	E_SFI_MVLEM_3DKlocal(11, 12) = (Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 13) = (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - e / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 17) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 18) = (Kh * c * h) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 19) = -e / (2.0 * (2.0 * (d * d) + 2.0)) - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(11, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(11, 23) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(12, 0) = E_SFI_MVLEM_3DKlocal(0, 12);
	E_SFI_MVLEM_3DKlocal(12, 1) = E_SFI_MVLEM_3DKlocal(1, 12);
	E_SFI_MVLEM_3DKlocal(12, 2) = E_SFI_MVLEM_3DKlocal(2, 12);
	E_SFI_MVLEM_3DKlocal(12, 3) = E_SFI_MVLEM_3DKlocal(3, 12);
	E_SFI_MVLEM_3DKlocal(12, 4) = E_SFI_MVLEM_3DKlocal(4, 12);
	E_SFI_MVLEM_3DKlocal(12, 5) = E_SFI_MVLEM_3DKlocal(5, 12);
	E_SFI_MVLEM_3DKlocal(12, 6) = E_SFI_MVLEM_3DKlocal(6, 12);
	E_SFI_MVLEM_3DKlocal(12, 7) = E_SFI_MVLEM_3DKlocal(7, 12);
	E_SFI_MVLEM_3DKlocal(12, 8) = E_SFI_MVLEM_3DKlocal(8, 12);
	E_SFI_MVLEM_3DKlocal(12, 9) = E_SFI_MVLEM_3DKlocal(9, 12);
	E_SFI_MVLEM_3DKlocal(12, 10) = E_SFI_MVLEM_3DKlocal(10, 12);
	E_SFI_MVLEM_3DKlocal(12, 11) = E_SFI_MVLEM_3DKlocal(11, 12);
	E_SFI_MVLEM_3DKlocal(12, 12) = Kh / 4.0 + (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(12, 13) = (Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(12, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 17) = -(Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(12, 18) = Kh / 4.0 - (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(12, 19) = -(Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(12, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(12, 23) = -(Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(13, 0) = E_SFI_MVLEM_3DKlocal(0, 13);
	E_SFI_MVLEM_3DKlocal(13, 1) = E_SFI_MVLEM_3DKlocal(1, 13);
	E_SFI_MVLEM_3DKlocal(13, 2) = E_SFI_MVLEM_3DKlocal(2, 13);
	E_SFI_MVLEM_3DKlocal(13, 3) = E_SFI_MVLEM_3DKlocal(3, 13);
	E_SFI_MVLEM_3DKlocal(13, 4) = E_SFI_MVLEM_3DKlocal(4, 13);
	E_SFI_MVLEM_3DKlocal(13, 5) = E_SFI_MVLEM_3DKlocal(5, 13);
	E_SFI_MVLEM_3DKlocal(13, 6) = E_SFI_MVLEM_3DKlocal(6, 13);
	E_SFI_MVLEM_3DKlocal(13, 7) = E_SFI_MVLEM_3DKlocal(7, 13);
	E_SFI_MVLEM_3DKlocal(13, 8) = E_SFI_MVLEM_3DKlocal(8, 13);
	E_SFI_MVLEM_3DKlocal(13, 9) = E_SFI_MVLEM_3DKlocal(9, 13);
	E_SFI_MVLEM_3DKlocal(13, 10) = E_SFI_MVLEM_3DKlocal(10, 13);
	E_SFI_MVLEM_3DKlocal(13, 11) = E_SFI_MVLEM_3DKlocal(11, 13);
	E_SFI_MVLEM_3DKlocal(13, 12) = E_SFI_MVLEM_3DKlocal(12, 13);
	E_SFI_MVLEM_3DKlocal(13, 13) = Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) + (12.0 * Eib * Iib) / (Lw * Lw * Lw) - (d * (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(13, 14) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 15) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 16) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 17) = (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	E_SFI_MVLEM_3DKlocal(13, 18) = (Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(13, 19) = Kv / 4.0 - (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) - (12.0 * Eib * Iib) / (Lw * Lw * Lw) + (d * (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(13, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(13, 23) = (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);

	E_SFI_MVLEM_3DKlocal(14, 0) = E_SFI_MVLEM_3DKlocal(0, 14);
	E_SFI_MVLEM_3DKlocal(14, 1) = E_SFI_MVLEM_3DKlocal(1, 14);
	E_SFI_MVLEM_3DKlocal(14, 2) = E_SFI_MVLEM_3DKlocal(2, 14);
	E_SFI_MVLEM_3DKlocal(14, 3) = E_SFI_MVLEM_3DKlocal(3, 14);
	E_SFI_MVLEM_3DKlocal(14, 4) = E_SFI_MVLEM_3DKlocal(4, 14);
	E_SFI_MVLEM_3DKlocal(14, 5) = E_SFI_MVLEM_3DKlocal(5, 14);
	E_SFI_MVLEM_3DKlocal(14, 6) = E_SFI_MVLEM_3DKlocal(6, 14);
	E_SFI_MVLEM_3DKlocal(14, 7) = E_SFI_MVLEM_3DKlocal(7, 14);
	E_SFI_MVLEM_3DKlocal(14, 8) = E_SFI_MVLEM_3DKlocal(8, 14);
	E_SFI_MVLEM_3DKlocal(14, 9) = E_SFI_MVLEM_3DKlocal(9, 14);
	E_SFI_MVLEM_3DKlocal(14, 10) = E_SFI_MVLEM_3DKlocal(10, 14);
	E_SFI_MVLEM_3DKlocal(14, 11) = E_SFI_MVLEM_3DKlocal(11, 14);
	E_SFI_MVLEM_3DKlocal(14, 12) = E_SFI_MVLEM_3DKlocal(12, 14);
	E_SFI_MVLEM_3DKlocal(14, 13) = E_SFI_MVLEM_3DKlocal(13, 14);
	E_SFI_MVLEM_3DKlocal(14, 14) = K1;
	E_SFI_MVLEM_3DKlocal(14, 15) = K2;
	E_SFI_MVLEM_3DKlocal(14, 16) = K3;
	E_SFI_MVLEM_3DKlocal(14, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(14, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(14, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(14, 20) = K4;
	E_SFI_MVLEM_3DKlocal(14, 21) = -K5;
	E_SFI_MVLEM_3DKlocal(14, 22) = K6;
	E_SFI_MVLEM_3DKlocal(14, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(15, 0) = E_SFI_MVLEM_3DKlocal(0, 15);
	E_SFI_MVLEM_3DKlocal(15, 1) = E_SFI_MVLEM_3DKlocal(1, 15);
	E_SFI_MVLEM_3DKlocal(15, 2) = E_SFI_MVLEM_3DKlocal(2, 15);
	E_SFI_MVLEM_3DKlocal(15, 3) = E_SFI_MVLEM_3DKlocal(3, 15);
	E_SFI_MVLEM_3DKlocal(15, 4) = E_SFI_MVLEM_3DKlocal(4, 15);
	E_SFI_MVLEM_3DKlocal(15, 5) = E_SFI_MVLEM_3DKlocal(5, 15);
	E_SFI_MVLEM_3DKlocal(15, 6) = E_SFI_MVLEM_3DKlocal(6, 15);
	E_SFI_MVLEM_3DKlocal(15, 7) = E_SFI_MVLEM_3DKlocal(7, 15);
	E_SFI_MVLEM_3DKlocal(15, 8) = E_SFI_MVLEM_3DKlocal(8, 15);
	E_SFI_MVLEM_3DKlocal(15, 9) = E_SFI_MVLEM_3DKlocal(9, 15);
	E_SFI_MVLEM_3DKlocal(15, 10) = E_SFI_MVLEM_3DKlocal(10, 15);
	E_SFI_MVLEM_3DKlocal(15, 11) = E_SFI_MVLEM_3DKlocal(11, 15);
	E_SFI_MVLEM_3DKlocal(15, 12) = E_SFI_MVLEM_3DKlocal(12, 15);
	E_SFI_MVLEM_3DKlocal(15, 13) = E_SFI_MVLEM_3DKlocal(13, 15);
	E_SFI_MVLEM_3DKlocal(15, 14) = E_SFI_MVLEM_3DKlocal(14, 15);
	E_SFI_MVLEM_3DKlocal(15, 15) = K13;
	E_SFI_MVLEM_3DKlocal(15, 16) = -K14;
	E_SFI_MVLEM_3DKlocal(15, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(15, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(15, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(15, 20) = -K5;
	E_SFI_MVLEM_3DKlocal(15, 21) = K15;
	E_SFI_MVLEM_3DKlocal(15, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(15, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(16, 0) = E_SFI_MVLEM_3DKlocal(0, 16);
	E_SFI_MVLEM_3DKlocal(16, 1) = E_SFI_MVLEM_3DKlocal(1, 16);
	E_SFI_MVLEM_3DKlocal(16, 2) = E_SFI_MVLEM_3DKlocal(2, 16);
	E_SFI_MVLEM_3DKlocal(16, 3) = E_SFI_MVLEM_3DKlocal(3, 16);
	E_SFI_MVLEM_3DKlocal(16, 4) = E_SFI_MVLEM_3DKlocal(4, 16);
	E_SFI_MVLEM_3DKlocal(16, 5) = E_SFI_MVLEM_3DKlocal(5, 16);
	E_SFI_MVLEM_3DKlocal(16, 6) = E_SFI_MVLEM_3DKlocal(6, 16);
	E_SFI_MVLEM_3DKlocal(16, 7) = E_SFI_MVLEM_3DKlocal(7, 16);
	E_SFI_MVLEM_3DKlocal(16, 8) = E_SFI_MVLEM_3DKlocal(8, 16);
	E_SFI_MVLEM_3DKlocal(16, 9) = E_SFI_MVLEM_3DKlocal(9, 16);
	E_SFI_MVLEM_3DKlocal(16, 10) = E_SFI_MVLEM_3DKlocal(10, 16);
	E_SFI_MVLEM_3DKlocal(16, 11) = E_SFI_MVLEM_3DKlocal(11, 16);
	E_SFI_MVLEM_3DKlocal(16, 12) = E_SFI_MVLEM_3DKlocal(12, 16);
	E_SFI_MVLEM_3DKlocal(16, 13) = E_SFI_MVLEM_3DKlocal(13, 16);
	E_SFI_MVLEM_3DKlocal(16, 14) = E_SFI_MVLEM_3DKlocal(14, 16);
	E_SFI_MVLEM_3DKlocal(16, 15) = E_SFI_MVLEM_3DKlocal(15, 16);
	E_SFI_MVLEM_3DKlocal(16, 16) = K18;
	E_SFI_MVLEM_3DKlocal(16, 17) = 0.0;
	E_SFI_MVLEM_3DKlocal(16, 18) = 0.0;
	E_SFI_MVLEM_3DKlocal(16, 19) = 0.0;
	E_SFI_MVLEM_3DKlocal(16, 20) = -K6;
	E_SFI_MVLEM_3DKlocal(16, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(16, 22) = K19;
	E_SFI_MVLEM_3DKlocal(16, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(17, 0) = E_SFI_MVLEM_3DKlocal(0, 17);
	E_SFI_MVLEM_3DKlocal(17, 1) = E_SFI_MVLEM_3DKlocal(1, 17);
	E_SFI_MVLEM_3DKlocal(17, 2) = E_SFI_MVLEM_3DKlocal(2, 17);
	E_SFI_MVLEM_3DKlocal(17, 3) = E_SFI_MVLEM_3DKlocal(3, 17);
	E_SFI_MVLEM_3DKlocal(17, 4) = E_SFI_MVLEM_3DKlocal(4, 17);
	E_SFI_MVLEM_3DKlocal(17, 5) = E_SFI_MVLEM_3DKlocal(5, 17);
	E_SFI_MVLEM_3DKlocal(17, 6) = E_SFI_MVLEM_3DKlocal(6, 17);
	E_SFI_MVLEM_3DKlocal(17, 7) = E_SFI_MVLEM_3DKlocal(7, 17);
	E_SFI_MVLEM_3DKlocal(17, 8) = E_SFI_MVLEM_3DKlocal(8, 17);
	E_SFI_MVLEM_3DKlocal(17, 9) = E_SFI_MVLEM_3DKlocal(9, 17);
	E_SFI_MVLEM_3DKlocal(17, 10) = E_SFI_MVLEM_3DKlocal(10, 17);
	E_SFI_MVLEM_3DKlocal(17, 11) = E_SFI_MVLEM_3DKlocal(11, 17);
	E_SFI_MVLEM_3DKlocal(17, 12) = E_SFI_MVLEM_3DKlocal(12, 17);
	E_SFI_MVLEM_3DKlocal(17, 13) = E_SFI_MVLEM_3DKlocal(13, 17);
	E_SFI_MVLEM_3DKlocal(17, 14) = E_SFI_MVLEM_3DKlocal(14, 17);
	E_SFI_MVLEM_3DKlocal(17, 15) = E_SFI_MVLEM_3DKlocal(15, 17);
	E_SFI_MVLEM_3DKlocal(17, 16) = E_SFI_MVLEM_3DKlocal(16, 17);
	E_SFI_MVLEM_3DKlocal(17, 17) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	E_SFI_MVLEM_3DKlocal(17, 18) = -(Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(17, 19) = e / (2.0 * (2.0 * (d * d) + 2.0)) - (6.0 * Eib * Iib) / (Lw * Lw) + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(17, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(17, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(17, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(17, 23) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (2.0 * Eib * Iib) / Lw;

	E_SFI_MVLEM_3DKlocal(18, 0) = E_SFI_MVLEM_3DKlocal(0, 18);
	E_SFI_MVLEM_3DKlocal(18, 1) = E_SFI_MVLEM_3DKlocal(1, 18);
	E_SFI_MVLEM_3DKlocal(18, 2) = E_SFI_MVLEM_3DKlocal(2, 18);
	E_SFI_MVLEM_3DKlocal(18, 3) = E_SFI_MVLEM_3DKlocal(3, 18);
	E_SFI_MVLEM_3DKlocal(18, 4) = E_SFI_MVLEM_3DKlocal(4, 18);
	E_SFI_MVLEM_3DKlocal(18, 5) = E_SFI_MVLEM_3DKlocal(5, 18);
	E_SFI_MVLEM_3DKlocal(18, 6) = E_SFI_MVLEM_3DKlocal(6, 18);
	E_SFI_MVLEM_3DKlocal(18, 7) = E_SFI_MVLEM_3DKlocal(7, 18);
	E_SFI_MVLEM_3DKlocal(18, 8) = E_SFI_MVLEM_3DKlocal(8, 18);
	E_SFI_MVLEM_3DKlocal(18, 9) = E_SFI_MVLEM_3DKlocal(9, 18);
	E_SFI_MVLEM_3DKlocal(18, 10) = E_SFI_MVLEM_3DKlocal(10, 18);
	E_SFI_MVLEM_3DKlocal(18, 11) = E_SFI_MVLEM_3DKlocal(11, 18);
	E_SFI_MVLEM_3DKlocal(18, 12) = E_SFI_MVLEM_3DKlocal(12, 18);
	E_SFI_MVLEM_3DKlocal(18, 13) = E_SFI_MVLEM_3DKlocal(13, 18);
	E_SFI_MVLEM_3DKlocal(18, 14) = E_SFI_MVLEM_3DKlocal(14, 18);
	E_SFI_MVLEM_3DKlocal(18, 15) = E_SFI_MVLEM_3DKlocal(15, 18);
	E_SFI_MVLEM_3DKlocal(18, 16) = E_SFI_MVLEM_3DKlocal(16, 18);
	E_SFI_MVLEM_3DKlocal(18, 17) = E_SFI_MVLEM_3DKlocal(17, 18);
	E_SFI_MVLEM_3DKlocal(18, 18) = Kh / 4.0 + (Aib * Eib) / Lw;
	E_SFI_MVLEM_3DKlocal(18, 19) = -(Kh * d * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));
	E_SFI_MVLEM_3DKlocal(18, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(18, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(18, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(18, 23) = -(Kh * h * (c - 1.0)) / (2.0 * (2.0 * (d * d) + 2.0));

	E_SFI_MVLEM_3DKlocal(19, 0) = E_SFI_MVLEM_3DKlocal(0, 19);
	E_SFI_MVLEM_3DKlocal(19, 1) = E_SFI_MVLEM_3DKlocal(1, 19);
	E_SFI_MVLEM_3DKlocal(19, 2) = E_SFI_MVLEM_3DKlocal(2, 19);
	E_SFI_MVLEM_3DKlocal(19, 3) = E_SFI_MVLEM_3DKlocal(3, 19);
	E_SFI_MVLEM_3DKlocal(19, 4) = E_SFI_MVLEM_3DKlocal(4, 19);
	E_SFI_MVLEM_3DKlocal(19, 5) = E_SFI_MVLEM_3DKlocal(5, 19);
	E_SFI_MVLEM_3DKlocal(19, 6) = E_SFI_MVLEM_3DKlocal(6, 19);
	E_SFI_MVLEM_3DKlocal(19, 7) = E_SFI_MVLEM_3DKlocal(7, 19);
	E_SFI_MVLEM_3DKlocal(19, 8) = E_SFI_MVLEM_3DKlocal(8, 19);
	E_SFI_MVLEM_3DKlocal(19, 9) = E_SFI_MVLEM_3DKlocal(9, 19);
	E_SFI_MVLEM_3DKlocal(19, 10) = E_SFI_MVLEM_3DKlocal(10, 19);
	E_SFI_MVLEM_3DKlocal(19, 11) = E_SFI_MVLEM_3DKlocal(11, 19);
	E_SFI_MVLEM_3DKlocal(19, 12) = E_SFI_MVLEM_3DKlocal(12, 19);
	E_SFI_MVLEM_3DKlocal(19, 13) = E_SFI_MVLEM_3DKlocal(13, 19);
	E_SFI_MVLEM_3DKlocal(19, 14) = E_SFI_MVLEM_3DKlocal(14, 19);
	E_SFI_MVLEM_3DKlocal(19, 15) = E_SFI_MVLEM_3DKlocal(15, 19);
	E_SFI_MVLEM_3DKlocal(19, 16) = E_SFI_MVLEM_3DKlocal(16, 19);
	E_SFI_MVLEM_3DKlocal(19, 17) = E_SFI_MVLEM_3DKlocal(17, 19);
	E_SFI_MVLEM_3DKlocal(19, 18) = E_SFI_MVLEM_3DKlocal(18, 19);
	E_SFI_MVLEM_3DKlocal(19, 19) = Kv / 4.0 + (d * e) / (2.0 * (2.0 * (d * d) + 2.0)) + (12.0 * Eib * Iib) / (Lw * Lw * Lw) + (d * (e / 2.0 + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	E_SFI_MVLEM_3DKlocal(19, 20) = 0.0;
	E_SFI_MVLEM_3DKlocal(19, 21) = 0.0;
	E_SFI_MVLEM_3DKlocal(19, 22) = 0.0;
	E_SFI_MVLEM_3DKlocal(19, 23) = (e / 2.0 + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) - (6.0 * Eib * Iib) / (Lw * Lw);

	E_SFI_MVLEM_3DKlocal(20, 0) = E_SFI_MVLEM_3DKlocal(0, 20);
	E_SFI_MVLEM_3DKlocal(20, 1) = E_SFI_MVLEM_3DKlocal(1, 20);
	E_SFI_MVLEM_3DKlocal(20, 2) = E_SFI_MVLEM_3DKlocal(2, 20);
	E_SFI_MVLEM_3DKlocal(20, 3) = E_SFI_MVLEM_3DKlocal(3, 20);
	E_SFI_MVLEM_3DKlocal(20, 4) = E_SFI_MVLEM_3DKlocal(4, 20);
	E_SFI_MVLEM_3DKlocal(20, 5) = E_SFI_MVLEM_3DKlocal(5, 20);
	E_SFI_MVLEM_3DKlocal(20, 6) = E_SFI_MVLEM_3DKlocal(6, 20);
	E_SFI_MVLEM_3DKlocal(20, 7) = E_SFI_MVLEM_3DKlocal(7, 20);
	E_SFI_MVLEM_3DKlocal(20, 8) = E_SFI_MVLEM_3DKlocal(8, 20);
	E_SFI_MVLEM_3DKlocal(20, 9) = E_SFI_MVLEM_3DKlocal(9, 20);
	E_SFI_MVLEM_3DKlocal(20, 10) = E_SFI_MVLEM_3DKlocal(10, 20);
	E_SFI_MVLEM_3DKlocal(20, 11) = E_SFI_MVLEM_3DKlocal(11, 20);
	E_SFI_MVLEM_3DKlocal(20, 12) = E_SFI_MVLEM_3DKlocal(12, 20);
	E_SFI_MVLEM_3DKlocal(20, 13) = E_SFI_MVLEM_3DKlocal(13, 20);
	E_SFI_MVLEM_3DKlocal(20, 14) = E_SFI_MVLEM_3DKlocal(14, 20);
	E_SFI_MVLEM_3DKlocal(20, 15) = E_SFI_MVLEM_3DKlocal(15, 20);
	E_SFI_MVLEM_3DKlocal(20, 16) = E_SFI_MVLEM_3DKlocal(16, 20);
	E_SFI_MVLEM_3DKlocal(20, 17) = E_SFI_MVLEM_3DKlocal(17, 20);
	E_SFI_MVLEM_3DKlocal(20, 18) = E_SFI_MVLEM_3DKlocal(18, 20);
	E_SFI_MVLEM_3DKlocal(20, 19) = E_SFI_MVLEM_3DKlocal(19, 20);
	E_SFI_MVLEM_3DKlocal(20, 20) = K1;
	E_SFI_MVLEM_3DKlocal(20, 21) = K2;
	E_SFI_MVLEM_3DKlocal(20, 22) = -K3;
	E_SFI_MVLEM_3DKlocal(20, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(21, 0) = E_SFI_MVLEM_3DKlocal(0, 21);
	E_SFI_MVLEM_3DKlocal(21, 1) = E_SFI_MVLEM_3DKlocal(1, 21);
	E_SFI_MVLEM_3DKlocal(21, 2) = E_SFI_MVLEM_3DKlocal(2, 21);
	E_SFI_MVLEM_3DKlocal(21, 3) = E_SFI_MVLEM_3DKlocal(3, 21);
	E_SFI_MVLEM_3DKlocal(21, 4) = E_SFI_MVLEM_3DKlocal(4, 21);
	E_SFI_MVLEM_3DKlocal(21, 5) = E_SFI_MVLEM_3DKlocal(5, 21);
	E_SFI_MVLEM_3DKlocal(21, 6) = E_SFI_MVLEM_3DKlocal(6, 21);
	E_SFI_MVLEM_3DKlocal(21, 7) = E_SFI_MVLEM_3DKlocal(7, 21);
	E_SFI_MVLEM_3DKlocal(21, 8) = E_SFI_MVLEM_3DKlocal(8, 21);
	E_SFI_MVLEM_3DKlocal(21, 9) = E_SFI_MVLEM_3DKlocal(9, 21);
	E_SFI_MVLEM_3DKlocal(21, 10) = E_SFI_MVLEM_3DKlocal(10, 21);
	E_SFI_MVLEM_3DKlocal(21, 11) = E_SFI_MVLEM_3DKlocal(11, 21);
	E_SFI_MVLEM_3DKlocal(21, 12) = E_SFI_MVLEM_3DKlocal(12, 21);
	E_SFI_MVLEM_3DKlocal(21, 13) = E_SFI_MVLEM_3DKlocal(13, 21);
	E_SFI_MVLEM_3DKlocal(21, 14) = E_SFI_MVLEM_3DKlocal(14, 21);
	E_SFI_MVLEM_3DKlocal(21, 15) = E_SFI_MVLEM_3DKlocal(15, 21);
	E_SFI_MVLEM_3DKlocal(21, 16) = E_SFI_MVLEM_3DKlocal(16, 21);
	E_SFI_MVLEM_3DKlocal(21, 17) = E_SFI_MVLEM_3DKlocal(17, 21);
	E_SFI_MVLEM_3DKlocal(21, 18) = E_SFI_MVLEM_3DKlocal(18, 21);
	E_SFI_MVLEM_3DKlocal(21, 19) = E_SFI_MVLEM_3DKlocal(19, 21);
	E_SFI_MVLEM_3DKlocal(21, 20) = E_SFI_MVLEM_3DKlocal(20, 21);
	E_SFI_MVLEM_3DKlocal(21, 21) = K13;
	E_SFI_MVLEM_3DKlocal(21, 22) = K14;
	E_SFI_MVLEM_3DKlocal(21, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(22, 0) = E_SFI_MVLEM_3DKlocal(0, 22);
	E_SFI_MVLEM_3DKlocal(22, 1) = E_SFI_MVLEM_3DKlocal(1, 22);
	E_SFI_MVLEM_3DKlocal(22, 2) = E_SFI_MVLEM_3DKlocal(2, 22);
	E_SFI_MVLEM_3DKlocal(22, 3) = E_SFI_MVLEM_3DKlocal(3, 22);
	E_SFI_MVLEM_3DKlocal(22, 4) = E_SFI_MVLEM_3DKlocal(4, 22);
	E_SFI_MVLEM_3DKlocal(22, 5) = E_SFI_MVLEM_3DKlocal(5, 22);
	E_SFI_MVLEM_3DKlocal(22, 6) = E_SFI_MVLEM_3DKlocal(6, 22);
	E_SFI_MVLEM_3DKlocal(22, 7) = E_SFI_MVLEM_3DKlocal(7, 22);
	E_SFI_MVLEM_3DKlocal(22, 8) = E_SFI_MVLEM_3DKlocal(8, 22);
	E_SFI_MVLEM_3DKlocal(22, 9) = E_SFI_MVLEM_3DKlocal(9, 22);
	E_SFI_MVLEM_3DKlocal(22, 10) = E_SFI_MVLEM_3DKlocal(10, 22);
	E_SFI_MVLEM_3DKlocal(22, 11) = E_SFI_MVLEM_3DKlocal(11, 22);
	E_SFI_MVLEM_3DKlocal(22, 12) = E_SFI_MVLEM_3DKlocal(12, 22);
	E_SFI_MVLEM_3DKlocal(22, 13) = E_SFI_MVLEM_3DKlocal(13, 22);
	E_SFI_MVLEM_3DKlocal(22, 14) = E_SFI_MVLEM_3DKlocal(14, 22);
	E_SFI_MVLEM_3DKlocal(22, 15) = E_SFI_MVLEM_3DKlocal(15, 22);
	E_SFI_MVLEM_3DKlocal(22, 16) = E_SFI_MVLEM_3DKlocal(16, 22);
	E_SFI_MVLEM_3DKlocal(22, 17) = E_SFI_MVLEM_3DKlocal(17, 22);
	E_SFI_MVLEM_3DKlocal(22, 18) = E_SFI_MVLEM_3DKlocal(18, 22);
	E_SFI_MVLEM_3DKlocal(22, 19) = E_SFI_MVLEM_3DKlocal(19, 22);
	E_SFI_MVLEM_3DKlocal(22, 20) = E_SFI_MVLEM_3DKlocal(20, 22);
	E_SFI_MVLEM_3DKlocal(22, 21) = E_SFI_MVLEM_3DKlocal(21, 22);
	E_SFI_MVLEM_3DKlocal(22, 22) = K18;
	E_SFI_MVLEM_3DKlocal(22, 23) = 0.0;

	E_SFI_MVLEM_3DKlocal(23, 0) = E_SFI_MVLEM_3DKlocal(0, 23);
	E_SFI_MVLEM_3DKlocal(23, 1) = E_SFI_MVLEM_3DKlocal(1, 23);
	E_SFI_MVLEM_3DKlocal(23, 2) = E_SFI_MVLEM_3DKlocal(2, 23);
	E_SFI_MVLEM_3DKlocal(23, 3) = E_SFI_MVLEM_3DKlocal(3, 23);
	E_SFI_MVLEM_3DKlocal(23, 4) = E_SFI_MVLEM_3DKlocal(4, 23);
	E_SFI_MVLEM_3DKlocal(23, 5) = E_SFI_MVLEM_3DKlocal(5, 23);
	E_SFI_MVLEM_3DKlocal(23, 6) = E_SFI_MVLEM_3DKlocal(6, 23);
	E_SFI_MVLEM_3DKlocal(23, 7) = E_SFI_MVLEM_3DKlocal(7, 23);
	E_SFI_MVLEM_3DKlocal(23, 8) = E_SFI_MVLEM_3DKlocal(8, 23);
	E_SFI_MVLEM_3DKlocal(23, 9) = E_SFI_MVLEM_3DKlocal(9, 23);
	E_SFI_MVLEM_3DKlocal(23, 10) = E_SFI_MVLEM_3DKlocal(10, 23);
	E_SFI_MVLEM_3DKlocal(23, 11) = E_SFI_MVLEM_3DKlocal(11, 23);
	E_SFI_MVLEM_3DKlocal(23, 12) = E_SFI_MVLEM_3DKlocal(12, 23);
	E_SFI_MVLEM_3DKlocal(23, 13) = E_SFI_MVLEM_3DKlocal(13, 23);
	E_SFI_MVLEM_3DKlocal(23, 14) = E_SFI_MVLEM_3DKlocal(14, 23);
	E_SFI_MVLEM_3DKlocal(23, 15) = E_SFI_MVLEM_3DKlocal(15, 23);
	E_SFI_MVLEM_3DKlocal(23, 16) = E_SFI_MVLEM_3DKlocal(16, 23);
	E_SFI_MVLEM_3DKlocal(23, 17) = E_SFI_MVLEM_3DKlocal(17, 23);
	E_SFI_MVLEM_3DKlocal(23, 18) = E_SFI_MVLEM_3DKlocal(18, 23);
	E_SFI_MVLEM_3DKlocal(23, 19) = E_SFI_MVLEM_3DKlocal(19, 23);
	E_SFI_MVLEM_3DKlocal(23, 20) = E_SFI_MVLEM_3DKlocal(20, 23);
	E_SFI_MVLEM_3DKlocal(23, 21) = E_SFI_MVLEM_3DKlocal(21, 23);
	E_SFI_MVLEM_3DKlocal(23, 22) = E_SFI_MVLEM_3DKlocal(22, 23);
	E_SFI_MVLEM_3DKlocal(23, 23) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;

	E_SFI_MVLEM_3DK.addMatrixTripleProduct(0.0, T, E_SFI_MVLEM_3DKlocal, 1.0); // Convert matrix from local to global cs

	// Return element Global stiffness matrix
	return E_SFI_MVLEM_3DK;

}

// Get element mass matrix assuming lumped mass
const Matrix& E_SFI_MVLEM_3D::getMass(void)
{

	E_SFI_MVLEM_3DM.Zero();
	E_SFI_MVLEM_3DMlocal.Zero();

	// No rotational mass
	E_SFI_MVLEM_3DMlocal(0, 0) = NodeMass;
	E_SFI_MVLEM_3DMlocal(1, 1) = NodeMass;
	E_SFI_MVLEM_3DMlocal(2, 2) = NodeMass;

	E_SFI_MVLEM_3DMlocal(6, 6) = NodeMass;
	E_SFI_MVLEM_3DMlocal(7, 7) = NodeMass;
	E_SFI_MVLEM_3DMlocal(8, 8) = NodeMass;

	E_SFI_MVLEM_3DMlocal(12, 12) = NodeMass;
	E_SFI_MVLEM_3DMlocal(13, 13) = NodeMass;
	E_SFI_MVLEM_3DMlocal(14, 14) = NodeMass;

	E_SFI_MVLEM_3DMlocal(18, 18) = NodeMass;
	E_SFI_MVLEM_3DMlocal(19, 19) = NodeMass;
	E_SFI_MVLEM_3DMlocal(20, 20) = NodeMass;

	// Convert matrix from local to global cs
	E_SFI_MVLEM_3DM.addMatrixTripleProduct(0.0, T, E_SFI_MVLEM_3DMlocal, 1.0);

	// Return element mass matrix
	return E_SFI_MVLEM_3DM;
}

// Get element damping matrix
const Matrix& E_SFI_MVLEM_3D::getDamp(void)
{
	E_SFI_MVLEM_3DD.Zero();

	E_SFI_MVLEM_3DD = this->Element::getDamp();

	// Return element damping matrix
	return E_SFI_MVLEM_3DD;
}

// N/A to this model - no element loads
void E_SFI_MVLEM_3D::zeroLoad(void)
{
	// does nothing - no elemental loads
}

// N/A to this model - no element loads
int E_SFI_MVLEM_3D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return 0;
}


int E_SFI_MVLEM_3D::addInertiaLoadToUnbalance(const Vector& accel)
{
	if (density == 0.0)
		return 0;

	// Get R * accel from the nodes
	const Vector& Raccel1 = theNodes[0]->getRV(accel);
	const Vector& Raccel2 = theNodes[1]->getRV(accel);
	const Vector& Raccel3 = theNodes[2]->getRV(accel);
	const Vector& Raccel4 = theNodes[3]->getRV(accel);

	if (6 != Raccel1.Size() || 6 != Raccel2.Size() || 6 != Raccel3.Size() || 6 != Raccel4.Size()) {
		opserr << "FourNodeQuad::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
		return -1;
	}

	Vector RaccelG(24);
	RaccelG.Zero();
	Vector RaccelL(24);
	RaccelL.Zero();

	// Assign nodal accelerations in global cs into a vector
	for (int i = 0; i < 6; i++) {
		RaccelG(i) = Raccel1(i);
		RaccelG(i + 6) = Raccel2(i);
		RaccelG(i + 12) = Raccel3(i);
		RaccelG(i + 18) = Raccel4(i);
	}

	// Transform accelerations from global to local cs
	RaccelL.addMatrixVector(0.0, T, RaccelG, 1.0);

	// Compute mass matrix
	this->getMass();

	// Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			E_SFI_MVLEM_3DRlocal(6 * i + j) += -E_SFI_MVLEM_3DMlocal(6 * i + j, 6 * i + j) * RaccelL(6 * i + j);
		}
	}

	// Transform forces from local to global cs
	E_SFI_MVLEM_3DR.addMatrixTransposeVector(1.0, T, E_SFI_MVLEM_3DRlocal, 1.0);

	return 0;
}

// Get element force vector
const Vector& E_SFI_MVLEM_3D::getResistingForce()
{

	E_SFI_MVLEM_3DR.Zero();
	E_SFI_MVLEM_3DRlocal.Zero();

	// Get Trial Displacements
	const Vector& disp1 = theNodes[0]->getTrialDisp();
	const Vector& disp2 = theNodes[1]->getTrialDisp();
	const Vector& disp3 = theNodes[2]->getTrialDisp();
	const Vector& disp4 = theNodes[3]->getTrialDisp();

	Vector dispG(24); // Vector of total 24 displacemets in global coordinates
	dispG.Zero();
	Vector dispL(24); // Vector of total 24 displacemets in local coordinates
	dispL.Zero();

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
		const Vector& Stress = theMaterial[i]->getStress();

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
		Fh += -1.0 * Fxy[i];
		Fysum += Fy[i];
	}

	R1 = Fh;
	R2 = -Fysum;
	R3 = -Fh * c * h;
	R4 = -Fh;
	R5 = Fysum;
	R6 = -Fh * (1.0 - c) * h;

	for (int i = 0; i < m; i++) {
		R3 -= Fy[i] * x[i];
		R6 += +Fy[i] * x[i];
	}

	// Calculate force vector in local cs
	E_SFI_MVLEM_3DRlocal(0) = R1 / 2.0 + (Aib * Eib * dispL(0)) / Lw - (Aib * Eib * dispL(6)) / Lw;
	E_SFI_MVLEM_3DRlocal(1) = R2 / 2.0 - (R3 * d) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib * dispL(1)) / (Lw * Lw * Lw) + (6.0 * Eib * Iib * dispL(5)) / (Lw * Lw) - (12.0 * Eib * Iib * dispL(7)) / (Lw * Lw * Lw) + (6.0 * Eib * Iib * dispL(11)) / (Lw * Lw);
	E_SFI_MVLEM_3DRlocal(2) = (Eave * (Tave * Tave * Tave) * (dispL(9)) * (4.0 * (h * h) * NUelastic + (h * h) - 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(3)) * (4.0 * (h * h) * NUelastic + (h * h) + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(15)) * ((h * h) - (h * h) * NUelastic + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(4)) * (4.0 * NUelastic * (Lw * Lw) + 10.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(10)) * (10.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(16)) * (4.0 * NUelastic * (Lw * Lw) - 5.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(2)) * (10.0 * (h * h * h * h) + 10.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(8)) * (10.0 * (h * h * h * h) - 5.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(14)) * (5.0 * (h * h * h * h) - 10.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(20)) * (5.0 * (h * h * h * h) + 5.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(21)) * ((h * h) * NUelastic - (h * h) + 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(22)) * (NUelastic * (Lw * Lw) + 5.0 * (h * h) - (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(3) = (Eave * NUelastic * (Tave * Tave * Tave) * (dispL(4))) / (12.0 * (NUelastic * NUelastic) - 12.0) - (Eave * (Tave * Tave * Tave) * (dispL(3)) * ((h * h) - (h * h) * NUelastic + 5 * (Lw * Lw))) / (45.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(2)) * (4.0 * (h * h) * NUelastic + (h * h) + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(8)) * (4.0 * (h * h) * NUelastic + (h * h) - 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1)) + (Eave * (Tave * Tave * Tave) * (dispL(14)) * ((h * h) - (h * h) * NUelastic + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(21)) * ((h * h) - (h * h) * NUelastic + 5.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(9)) * (2.0 * (h * h) * NUelastic - 2.0 * (h * h) + 5.0 * (Lw * Lw))) / (90.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(15)) * ((h * h) * NUelastic - (h * h) + 10.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(20)) * ((h * h) * NUelastic - (h * h) + 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(4) = (Eave * NUelastic * (Tave * Tave * Tave) * (dispL(3))) / (12.0 * (NUelastic * NUelastic) - 12.0) - (dispL(4)) * ((Eave * h * (Tave * Tave * Tave)) / (9.0 * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (45.0 * h * ((NUelastic * NUelastic) - 1.0))) - (dispL(10)) * ((Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0))) - (dispL(22)) * ((Eave * h * (Tave * Tave * Tave)) / (36.0 * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0))) - (dispL(16)) * ((Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * Lw * (2.0 * NUelastic - 2.0)) / (90.0 * h * ((NUelastic * NUelastic) - 1.0))) + (Eave * (Tave * Tave * Tave) * (dispL(2)) * (4.0 * NUelastic * (Lw * Lw) + 10.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(8)) * (10.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(14)) * (4.0 * NUelastic * (Lw * Lw) - 5.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(20)) * (NUelastic * (Lw * Lw) + 5.0 * (h * h) - (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(5) = R3 / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib * dispL(1)) / (Lw * Lw) + (4.0 * Eib * Iib * dispL(5)) / Lw - (6.0 * Eib * Iib * dispL(7)) / (Lw * Lw) + (2.0 * Eib * Iib * dispL(11)) / Lw;
	E_SFI_MVLEM_3DRlocal(6) = R1 / 2.0 - (Aib * Eib * dispL(0)) / Lw + (Aib * Eib * dispL(6)) / Lw;
	E_SFI_MVLEM_3DRlocal(7) = R2 / 2.0 + (R3 * d) / (2.0 * (d * d) + 2.0) - (12.0 * Eib * Iib * dispL(1)) / (Lw * Lw * Lw) - (6.0 * Eib * Iib * dispL(5)) / (Lw * Lw) + (12.0 * Eib * Iib * dispL(7)) / (Lw * Lw * Lw) - (6.0 * Eib * Iib * dispL(11)) / (Lw * Lw);
	E_SFI_MVLEM_3DRlocal(8) = (Eave * (Tave * Tave * Tave) * (dispL(3)) * (4.0 * (h * h) * NUelastic + (h * h) - 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(9)) * (4.0 * (h * h) * NUelastic + (h * h) + 10 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(21)) * ((h * h) - (h * h) * NUelastic + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(4)) * (10.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(10)) * (4.0 * NUelastic * (Lw * Lw) + 10.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(22)) * (4.0 * NUelastic * (Lw * Lw) - 5.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(2)) * (10.0 * (h * h * h * h) - 5.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(14)) * (5.0 * (h * h * h * h) + 5.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(8)) * (10.0 * (h * h * h * h) + 10.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(20)) * (5.0 * (h * h * h * h) - 10.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(15)) * ((h * h) * NUelastic - (h * h) + 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(16)) * (NUelastic * (Lw * Lw) + 5.0 * (h * h) - (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(9) = (Eave * (Tave * Tave * Tave) * (dispL(2)) * (4.0 * (h * h) * NUelastic + (h * h) - 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * NUelastic * (Tave * Tave * Tave) * (dispL(10))) / (12.0 * (NUelastic * NUelastic) - 12.0) - (Eave * (Tave * Tave * Tave) * (dispL(9)) * ((h * h) - (h * h) * NUelastic + 5.0 * (Lw * Lw))) / (45.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(15)) * ((h * h) - (h * h) * NUelastic + 5.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(8)) * (4.0 * (h * h) * NUelastic + (h * h) + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(20)) * ((h * h) - (h * h) * NUelastic + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(3)) * (2.0 * (h * h) * NUelastic - 2.0 * (h * h) + 5.0 * (Lw * Lw))) / (90.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(14)) * ((h * h) * NUelastic - (h * h) + 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(21)) * ((h * h) * NUelastic - (h * h) + 10.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(10) = (Eave * (Tave * Tave * Tave) * (dispL(2)) * (10.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (dispL(4)) * ((Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0))) - (dispL(16)) * ((Eave * h * (Tave * Tave * Tave)) / (36.0 * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0))) - (Eave * NUelastic * (Tave * Tave * Tave) * (dispL(9))) / (12.0 * (NUelastic * NUelastic) - 12.0) - (dispL(22)) * ((Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * Lw * (2.0 * NUelastic - 2.0)) / (90.0 * h * ((NUelastic * NUelastic) - 1.0))) - (Eave * (Tave * Tave * Tave) * (dispL(10)) * (5.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (45.0 * h * Lw * ((NUelastic * NUelastic) - 1)) - (Eave * (Tave * Tave * Tave) * (dispL(8)) * (4.0 * NUelastic * (Lw * Lw) + 10.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(20)) * (4.0 * NUelastic * (Lw * Lw) - 5.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(14)) * (NUelastic * (Lw * Lw) + 5.0 * (h * h) - (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(11) = R3 / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib * dispL(1)) / (Lw * Lw) + (2.0 * Eib * Iib * dispL(5)) / Lw - (6.0 * Eib * Iib * dispL(7)) / (Lw * Lw) + (4.0 * Eib * Iib * dispL(11)) / Lw;
	E_SFI_MVLEM_3DRlocal(12) = R4 / 2.0 + (Aib * Eib * dispL(12)) / Lw - (Aib * Eib * dispL(18)) / Lw;
	E_SFI_MVLEM_3DRlocal(13) = R5 / 2.0 - (R6 * d) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib * dispL(13)) / (Lw * Lw * Lw) + (6.0 * Eib * Iib * dispL(17)) / (Lw * Lw) - (12.0 * Eib * Iib * dispL(19)) / (Lw * Lw * Lw) + (6.0 * Eib * Iib * dispL(23)) / (Lw * Lw);
	E_SFI_MVLEM_3DRlocal(14) = (Eave * (Tave * Tave * Tave) * (dispL(3)) * ((h * h) - (h * h) * NUelastic + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(15)) * (4.0 * (h * h) * NUelastic + (h * h) + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(21)) * (4.0 * (h * h) * NUelastic + (h * h) - 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(4)) * (4.0 * NUelastic * (Lw * Lw) - 5.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(16)) * (4.0 * NUelastic * (Lw * Lw) + 10.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(22)) * (10.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(2)) * (5.0 * (h * h * h * h) - 10.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(8)) * (5.0 * (h * h * h * h) + 5.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(14)) * (10.0 * (h * h * h * h) + 10.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(20)) * (10.0 * (h * h * h * h) - 5.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(9)) * ((h * h) * NUelastic - (h * h) + 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(10)) * (NUelastic * (Lw * Lw) + 5.0 * (h * h) - (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(15) = (Eave * (Tave * Tave * Tave) * (dispL(14)) * (4.0 * (h * h) * NUelastic + (h * h) + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(2)) * ((h * h) - (h * h) * NUelastic + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(9)) * ((h * h) - (h * h) * NUelastic + 5.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(15)) * ((h * h) - (h * h) * NUelastic + 5.0 * (Lw * Lw))) / (45.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * NUelastic * (Tave * Tave * Tave) * (dispL(16))) / (12.0 * (NUelastic * NUelastic) - 12.0) - (Eave * (Tave * Tave * Tave) * (dispL(20)) * (4.0 * (h * h) * NUelastic + (h * h) - 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(3)) * ((h * h) * NUelastic - (h * h) + 10.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(8)) * ((h * h) * NUelastic - (h * h) + 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(21)) * (2.0 * (h * h) * NUelastic - 2.0 * (h * h) + 5.0 * (Lw * Lw))) / (90.0 * h * Lw * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(16) = (Eave * (Tave * Tave * Tave) * (dispL(14)) * (4.0 * NUelastic * (Lw * Lw) + 10.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (dispL(16)) * ((Eave * h * (Tave * Tave * Tave)) / (9.0 * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (45.0 * h * ((NUelastic * NUelastic) - 1.0))) - (dispL(22)) * ((Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0))) - (dispL(10)) * ((Eave * h * (Tave * Tave * Tave)) / (36.0 * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0))) - (Eave * NUelastic * (Tave * Tave * Tave) * (dispL(15))) / (12.0 * (NUelastic * NUelastic) - 12.0) - (Eave * (Tave * Tave * Tave) * (dispL(2)) * (4.0 * NUelastic * (Lw * Lw) - 5.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (dispL(4)) * ((Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * Lw * (2.0 * NUelastic - 2.0)) / (90.0 * h * ((NUelastic * NUelastic) - 1.0))) - (Eave * (Tave * Tave * Tave) * (dispL(20)) * (10.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(8)) * (NUelastic * (Lw * Lw) + 5.0 * (h * h) - (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(17) = R6 / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib * dispL(13)) / (Lw * Lw) + (4.0 * Eib * Iib * dispL(17)) / Lw - (6.0 * Eib * Iib * dispL(19)) / (Lw * Lw) + (2.0 * Eib * Iib * dispL(23)) / Lw;
	E_SFI_MVLEM_3DRlocal(18) = R4 / 2.0 - (Aib * Eib * dispL(12)) / Lw + (Aib * Eib * dispL(18)) / Lw;
	E_SFI_MVLEM_3DRlocal(19) = R5 / 2.0 + (R6 * d) / (2.0 * (d * d) + 2.0) - (12.0 * Eib * Iib * dispL(13)) / (Lw * Lw * Lw) - (6.0 * Eib * Iib * dispL(17)) / (Lw * Lw) + (12.0 * Eib * Iib * dispL(19)) / (Lw * Lw * Lw) - (6.0 * Eib * Iib * dispL(23)) / (Lw * Lw);
	E_SFI_MVLEM_3DRlocal(20) = (Eave * (Tave * Tave * Tave) * (dispL(9)) * ((h * h) - (h * h) * NUelastic + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(15)) * (4.0 * (h * h) * NUelastic + (h * h) - 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(21)) * (4.0 * (h * h) * NUelastic + (h * h) + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(10)) * (4.0 * NUelastic * (Lw * Lw) - 5.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(16)) * (10.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(22)) * (4.0 * NUelastic * (Lw * Lw) + 10.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(2)) * (5.0 * (h * h * h * h) + 5.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(8)) * (5.0 * (h * h * h * h) - 10.0 * (Lw * Lw * Lw * Lw) - 7.0 * (h * h) * (Lw * Lw) + 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(14)) * (10.0 * (h * h * h * h) - 5.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(20)) * (10.0 * (h * h * h * h) + 10.0 * (Lw * Lw * Lw * Lw) + 7.0 * (h * h) * (Lw * Lw) - 2.0 * (h * h) * NUelastic * (Lw * Lw))) / (30.0 * (h * h * h) * (Lw * Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(3)) * ((h * h) * NUelastic - (h * h) + 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(4)) * (NUelastic * (Lw * Lw) + 5.0 * (h * h) - (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(21) = (Eave * NUelastic * (Tave * Tave * Tave) * (dispL(22))) / (12.0 * (NUelastic * NUelastic) - 12.0) - (Eave * (Tave * Tave * Tave) * (dispL(3)) * ((h * h) - (h * h) * NUelastic + 5.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(8)) * ((h * h) - (h * h) * NUelastic + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(14)) * (4.0 * (h * h) * NUelastic + (h * h) - 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(21)) * ((h * h) - (h * h) * NUelastic + 5.0 * (Lw * Lw))) / (45.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(20)) * (4.0 * (h * h) * NUelastic + (h * h) + 10.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(2)) * ((h * h) * NUelastic - (h * h) + 5.0 * (Lw * Lw))) / (60.0 * (h * h) * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(9)) * ((h * h) * NUelastic - (h * h) + 10.0 * (Lw * Lw))) / (180.0 * h * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(15)) * (2.0 * (h * h) * NUelastic - 2.0 * (h * h) + 5.0 * (Lw * Lw))) / (90.0 * h * Lw * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(22) = (Eave * NUelastic * (Tave * Tave * Tave) * (dispL(21))) / (12.0 * (NUelastic * NUelastic) - 12.0) - (dispL(22)) * ((Eave * h * (Tave * Tave * Tave)) / (9.0 * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (45.0 * h * ((NUelastic * NUelastic) - 1.0))) - (dispL(16)) * ((Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0))) - (dispL(4)) * ((Eave * h * (Tave * Tave * Tave)) / (36.0 * Lw * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * Lw * (NUelastic - 1.0)) / (180.0 * h * ((NUelastic * NUelastic) - 1.0))) - (dispL(10)) * ((Eave * h * (Tave * Tave * Tave)) / (18.0 * Lw * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * Lw * (2.0 * NUelastic - 2.0)) / (90.0 * h * ((NUelastic * NUelastic) - 1.0))) + (Eave * (Tave * Tave * Tave) * (dispL(8)) * (4.0 * NUelastic * (Lw * Lw) - 5.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(14)) * (10.0 * (h * h) - NUelastic * (Lw * Lw) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) - (Eave * (Tave * Tave * Tave) * (dispL(20)) * (4.0 * NUelastic * (Lw * Lw) + 10.0 * (h * h) + (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0)) + (Eave * (Tave * Tave * Tave) * (dispL(2)) * (NUelastic * (Lw * Lw) + 5.0 * (h * h) - (Lw * Lw))) / (60.0 * h * (Lw * Lw) * ((NUelastic * NUelastic) - 1.0));
	E_SFI_MVLEM_3DRlocal(23) = R6 / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib * dispL(13)) / (Lw * Lw) + (2.0 * Eib * Iib * dispL(17)) / Lw - (6.0 * Eib * Iib * dispL(19)) / (Lw * Lw) + (4.0 * Eib * Iib * dispL(23)) / Lw;

	// Convert force vector from local to global cs
	E_SFI_MVLEM_3DR.addMatrixTransposeVector(0.0, T, E_SFI_MVLEM_3DRlocal, 1.0);

	// Return element force vector
	return E_SFI_MVLEM_3DR;
}

// Get resisting force increment from inertial forces
const Vector& E_SFI_MVLEM_3D::getResistingForceIncInertia()
{
	// if no mass terms .. just add damping terms
	if (density == 0.0) {
		// Compute the current resisting force
		this->getResistingForce();

		// Add the damping forces if rayleigh damping
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			E_SFI_MVLEM_3DR += this->getRayleighDampingForces();

		return E_SFI_MVLEM_3DR;
	}

	// Get nodal accelerations in global cs
	const Vector& accel1 = theNodes[0]->getTrialAccel();
	const Vector& accel2 = theNodes[1]->getTrialAccel();
	const Vector& accel3 = theNodes[2]->getTrialAccel();
	const Vector& accel4 = theNodes[3]->getTrialAccel();

	Vector accelG(24);
	accelG.Zero();
	Vector accelL(24);
	accelL.Zero();

	// Assign nodal accelerations in global cs into a vector
	for (int i = 0; i < 6; i++) {
		accelG(i) = accel1(i);
		accelG(i + 6) = accel2(i);
		accelG(i + 12) = accel3(i);
		accelG(i + 18) = accel4(i);
	}

	// Transform accelerations from global to local cs
	accelL.addMatrixVector(0.0, T, accelG, 1.0);

	// Compute the current resisting force
	this->getResistingForce();

	// Compute the mass matrix
	this->getMass();

	// Add inertia forces to force vector in local cs
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			E_SFI_MVLEM_3DRlocal(6 * i + j) += E_SFI_MVLEM_3DMlocal(6 * i + j, 6 * i + j) * accelL(6 * i + j);
		}
	}

	// Transform forces from local to global cs
	E_SFI_MVLEM_3DR.addMatrixTransposeVector(1.0, T, E_SFI_MVLEM_3DRlocal, 1.0);

	// Add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
		E_SFI_MVLEM_3DR += this->getRayleighDampingForces();

	return E_SFI_MVLEM_3DR;
}

// Send Self
int E_SFI_MVLEM_3D::sendSelf(int commitTag, Channel& theChannel)
{
	int res;
	int dataTag = this->getDbTag();

	static Vector data(6);

	data(0) = this->getTag();
	data(1) = density;
	data(2) = m;
	data(3) = c;
	data(4) = NUelastic;
	data(5) = Tfactor;

	// E_SFI_MVLEM_3D then sends the tags of it's nodes
	res = theChannel.sendID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING E_SFI_MVLEM_3D::sendSelf() - failed to send ID\n";
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

	return 0;
}

// Receive Self
int E_SFI_MVLEM_3D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res;
	int dataTag = this->getDbTag();

	// E_SFI_MVLEM_3D creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	// delete dynamic memory
	if (theMaterial != 0) {
		for (int i = 0; i < m; i++)
			if (theMaterial[i] != 0)
				delete theMaterial[i];
		delete[] theMaterial;
	}

	Vector data(6);
	res = theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING E_SFI_MVLEM_3D::recvSelf() - failed to receive Vector\n";
		return -1;
	}

	this->setTag((int)data(0));
	density = data(1);
	m = data(2);
	c = data(3);
	NUelastic = data(4);
	Tfactor = data(5);

	// E_SFI_MVLEM_3D now receives the tags of it's four external nodes
	res = theChannel.recvID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING E_SFI_MVLEM_3D::recvSelf() - failed to receive ID\n";
		return -2;
	}

	// Receive the material class tags
	ID matClassTags(m);
	res = theChannel.recvID(0, commitTag, matClassTags);

	// Allocate memory for the uniaxial materials
	theMaterial = new NDMaterial * [m];
	if (theMaterial == 0) {
		opserr << "E_SFI_MVLEM_3D::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Receive the material models
	for (int i = 0; i < m; i++) {
		theMaterial[i] = theBroker.getNewNDMaterial(matClassTags(i));
		if (theMaterial[i] == 0) {
			opserr << "E_SFI_MVLEM_3D::recvSelf() - "
				<< "failed to get blank uniaxial material.\n";
			return -3;
		}
		theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
	}

	return 0;
}

// Display model
int E_SFI_MVLEM_3D::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
	// First get the end points of the beam based on
	// the display factor (a measure of the distorted image)
	// get location of nodes

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

		const Vector& end1Disp = theNodes[0]->getDisp();
		const Vector& end2Disp = theNodes[1]->getDisp();
		const Vector& end3Disp = theNodes[2]->getDisp();
		const Vector& end4Disp = theNodes[3]->getDisp();

		for (int i = 0; i < 3; i++) { // loop over coordinates (3 for 3D elements)

	// add displacement (multiplied with the displacement facotr) to the original node location to obtain current node location
			Gv1(i) = nd1Crds(i) + end1Disp(i) * fact;
			Gv2(i) = nd2Crds(i) + end2Disp(i) * fact;
			Gv3(i) = nd3Crds(i) + end3Disp(i) * fact;
			Gv4(i) = nd4Crds(i) + end4Disp(i) * fact;

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

				Gv1(i) = nd1Crds(i) + eigen1(i, mode - 1) * fact;
				Gv2(i) = nd2Crds(i) + eigen2(i, mode - 1) * fact;
				Gv3(i) = nd3Crds(i) + eigen3(i, mode - 1) * fact;
				Gv4(i) = nd4Crds(i) + eigen4(i, mode - 1) * fact;

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
			const Vector& stress = theMaterial[panel]->getStrain();
			values(0) = stress(displayMode - 1);
		}

		// Determine the deformation - rotation - other is taken from v1, v2
		const Vector& end1Disp4G = theNodes[0]->getDisp();
		const Vector& end2Disp4G = theNodes[1]->getDisp();
		const Vector& end3Disp4G = theNodes[2]->getDisp();
		const Vector& end4Disp4G = theNodes[3]->getDisp();

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
		end1Disp(4) = end1Disp4(4) / (2.0 * (d * d) + 2.0) + end2Disp4(4) / (2.0 * (d * d) + 2.0) + (end1Disp4(2) * d) / (2.0 * (d * d) + 2.0) - (end2Disp4(2) * d) / (2.0 * (d * d) + 2.0);
		end1Disp(5) = end1Disp4(5) / (2.0 * (d * d) + 2.0) + end2Disp4(5) / (2.0 * (d * d) + 2.0) - (end1Disp4(1) * d) / (2.0 * (d * d) + 2.0) + (end2Disp4(1) * d) / (2.0 * (d * d) + 2.0);

		end2Disp(4) = end3Disp4(4) / (2.0 * (d * d) + 2.0) + end4Disp4(4) / (2.0 * (d * d) + 2.0) + (end3Disp4(2) * d) / (2.0 * (d * d) + 2.0) - (end4Disp4(2) * d) / (2.0 * (d * d) + 2.0);
		end2Disp(5) = end3Disp4(5) / (2.0 * (d * d) + 2.0) + end4Disp4(5) / (2.0 * (d * d) + 2.0) - (end3Disp4(1) * d) / (2.0 * (d * d) + 2.0) + (end4Disp4(1) * d) / (2.0 * (d * d) + 2.0);

		// Fiber nodes
		NodePLotCrds(panel, 0) = panel + 1; // panel id

		Vector LocCoord(3); LocCoord.Zero();
		Vector GlCoord(3); GlCoord.Zero();
		// Local node 1 - bottom left
		LocCoord(0) = Lv1_(0) + x[panel] - b[panel] / 2.0; // x 
		LocCoord(1) = Lv1_(1) + (x[panel] - b[panel] / 2.0) * end1Disp(5) * fact; // y
		LocCoord(2) = Lv1_(2) - (x[panel] - b[panel] / 2.0) * end1Disp(4) * fact; // z
		GlCoord.addMatrixTransposeVector(0.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 1) = GlCoord(0);
		NodePLotCrds(panel, 2) = GlCoord(1);
		NodePLotCrds(panel, 3) = GlCoord(2);
		LocCoord.Zero();
		GlCoord.Zero();
		// Local node 2 - bottom right
		LocCoord(0) = Lv1_(0) + x[panel] + b[panel] / 2.0; // x
		LocCoord(1) = Lv1_(1) + (x[panel] + b[panel] / 2.0) * end1Disp(5) * fact; // y
		LocCoord(2) = Lv1_(2) - (x[panel] + b[panel] / 2.0) * end1Disp(4) * fact; // z
		GlCoord.addMatrixTransposeVector(0.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 4) = GlCoord(0);
		NodePLotCrds(panel, 5) = GlCoord(1);
		NodePLotCrds(panel, 6) = GlCoord(2);
		LocCoord.Zero();
		GlCoord.Zero();
		// Local node 3 - top left
		LocCoord(0) = Lv2_(0) + x[panel] + b[panel] / 2.0; // x
		LocCoord(1) = Lv2_(1) + (x[panel] + b[panel] / 2.0) * end2Disp(5) * fact; // y
		LocCoord(2) = Lv2_(2) - (x[panel] + b[panel] / 2.0) * end2Disp(4) * fact; // z
		GlCoord.addMatrixTransposeVector(0.0, Tt, LocCoord, 1.0);
		NodePLotCrds(panel, 7) = GlCoord(0);
		NodePLotCrds(panel, 8) = GlCoord(1);
		NodePLotCrds(panel, 9) = GlCoord(2);
		LocCoord.Zero();
		GlCoord.Zero();
		// Local node 4 - top right
		LocCoord(0) = Lv2_(0) + x[panel] - b[panel] / 2.0; // x
		LocCoord(1) = Lv2_(1) + (x[panel] - b[panel] / 2.0) * end2Disp(5) * fact; // y
		LocCoord(2) = Lv2_(2) - (x[panel] - b[panel] / 2.0) * end2Disp(4) * fact; // z
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

	return error;

}

// Print Element Information
void E_SFI_MVLEM_3D::Print(OPS_Stream& s, int flag)
{
	if (flag == 0) {
		s << "E_SFI_MVLEM_3D Element tag: " << this->getTag() << endln;
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

}

// Set element responses
Response* E_SFI_MVLEM_3D::setResponse(const char** argv, int argc, OPS_Stream& s)
{

	Response* theResponse = 0;

	s.tag("ElementOutput");
	s.attr("eleType", "E_SFI_MVLEM_3D");
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

		theResponse = new ElementResponse(this, 1, Vector(24));

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

		theResponse = new ElementResponse(this, 2, Vector(24));

	}

	// Shear deformation
	else if (strcmp(argv[0], "ShearDef") == 0 || strcmp(argv[0], "sheardef") == 0) {

		s.tag("ResponseType", "Dsh");

		theResponse = new ElementResponse(this, 3, 0.0);

	}

	// Element curvature
	else if (strcmp(argv[0], "Curvature") == 0 || strcmp(argv[0], "curvature") == 0) {

		s.tag("ResponseType", "fi");

		theResponse = new ElementResponse(this, 4, 0.0);
	}

	// Material output
	else if (strcmp(argv[0], "RCpanel") == 0 || strcmp(argv[0], "RCPanel") == 0
		|| strcmp(argv[0], "RC_panel") == 0 || strcmp(argv[0], "RC_Panel") == 0
		|| strcmp(argv[0], "material") == 0)
	{

		// Check if correct # of arguments passed
		if (argc != 3) {
			opserr << "WARNING: Number of recorder input for RC Panel is: " << argc - 1 << "; should be 2: panTag (one panel only: 1 to m) and $Response_Type.\n";
			return 0;
		}

		int matNum = atoi(argv[1]);
		if (matNum > 0 && matNum <= m) {

			s.tag("GaussPointOutput");
			s.attr("number", matNum);
			s.attr("eta", x[matNum - 1] / Lw * 2.0);
			s.attr("weight", b[matNum - 1] / Lw * 2.0);

			theResponse = theMaterial[matNum - 1]->setResponse(&argv[argc - 1], argc - 2, s);
		}

	}

	s.endTag();

	return theResponse;
}

// Get shear deformation
double E_SFI_MVLEM_3D::getShearDef(void)
{
	return Dsh;
}

// Get curvature (from vertical strains)
double E_SFI_MVLEM_3D::getCurvature(void)
{
	double Curv;

	Curv = (E_SFI_MVLEM_3DStrainY[0] - E_SFI_MVLEM_3DStrainY[m - 1]) / (x[0] - x[m - 1]);

	return Curv;
}

// Get global forces at 24 DOFs (top and bottom node)
Vector E_SFI_MVLEM_3D::getResistingForce_24DOF(void)
{

	for (int i = 0; i < 24; i++) {
		P_24DOF(i) = E_SFI_MVLEM_3DR(i);
	}

	return P_24DOF;
}

// Obtain element responses
int E_SFI_MVLEM_3D::getResponse(int responseID, Information& eleInfo)
{

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

}

// Return element local forces 
Vector E_SFI_MVLEM_3D::getResistingForce_24DOF_local(void)
{
	for (int i = 0; i < 24; i++) {
		P_24DOF_local(i) = E_SFI_MVLEM_3DRlocal(i);
	}

	return P_24DOF_local;
}

// Compute element transformation matrix
void  E_SFI_MVLEM_3D::setTransformationMatrix(void) {

	T.Zero(); // element transformation matrix
	Tt.Zero(); // 3 x 3 - basic transformation matrix
	T6.Zero(); // 6 x 6 - nodal transformation matrix 

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
	X_ = pow(pow(Xx, 2) + pow(Xy, 2) + pow(Xz, 2), 0.5);

	// unit x components
	Xex = Xx / X_;
	Xey = Xy / X_;
	Xez = Xz / X_;

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
	Zey = -(Xex * Yez - Xez * Yex);
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

};
