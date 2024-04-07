// Code written/implemented by:	Carlos López Olea (carlos.lopez.o@ug.uchile.cl)
//								Leonardo M. Massone
//								Kristijan Kolozvari
//
// User documentation available at: https://github.com/carloslopezolea/E-SFI_Documentation
//
// Created: 06/2022
//
// Description: The Efficient Shear-Flexure Interaction (E-SFI) element was developed based on the SFI-MVLEM formulation. 
// The E-SFI element incorporates the shear-flexure interaction phenomenon by replacing the m number of uniaxial fibers of 
// the MVLEM, by two-dimensional RC panel elements subjected to membrane actions (FSAM). An E-SFI element is described by 
// six degrees of freedom, and therefore no additional degrees of freedom are incorporated into the original MVLEM formulation, 
// as in the SFI-MVLEM. The curvature of an E-SFI element is assumed to be uniform, and the resultant rotation is concentrated
// at height ch. The kinematic assumption of plane sections remain plane, as well as the assumption of constant shear strain 
// along the element length, are considered for computing the axial and shear strains for each panel over the entire section. 
// To complete the strain field of a panel element, a calibrated expression for the horizontal normal strain is implemented to 
// obtain accurate predictions from squat to slender RC walls.
//
// References:
// 1.- Massone, L. M., López, C. N., & Kolozvari, K. (2021). Formulation of an efficient shear-flexure interaction model for planar reinforced concrete walls. Engineering Structures, 243, 112680.
// 2.- López, C. N., Massone, L. M., & Kolozvari, K. (2022). Validation of an efficient shear-flexure interaction model for planar reinforced concrete walls. Engineering Structures, 252, 113590.
// 3.- López C. N. Efficient shear-flexure interaction model for nonlinear analysis of reinforced concrete structural walls. MS Dissertation. Santiago, Chile: University of Chile; 2021.
// 
// Source: /usr/local/cvs/OpenSees/SRC/element/mvlem/E_SFI.cpp
//
// Rev: 1.0

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>
#include <E_SFI.h>
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
#include <elementAPI.h>

// Read input parameters and build the material
void* OPS_E_SFI()
{
	// Pointer to a uniaxial material that will be returned                       
	Element* theElement = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numArgs < 7) {
		opserr << "Want: E_SFI eleTag iNode jNode m c -thick -width -mat\n";
		return 0;
	}

	int iData[4];
	double dData[1];

	int numData = 4;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid int data for element E_SFI" << endln;
		return 0;
	}

	numData = 1;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid c for element E_SFI " << iData[0] << endln;
		return 0;
	}

	int m = iData[3];
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
	while (numArgs >= (m + 1)) {
		//OPS_GetStringCopy(&str);
		str = OPS_GetString();
		if (strcmp(str, "-thick") == 0) {
			numData = m;
			if (OPS_GetDoubleInput(&numData, theThickness) != 0) {
				opserr << "Invalid thick parameter for E_SFI " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-width") == 0) {
			numData = m;
			if (OPS_GetDoubleInput(&numData, theWidth) != 0) {
				opserr << "Invalid width value for E_SFI " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-mat") == 0) {
			numData = m;
			if (OPS_GetIntInput(&numData, matTags) != 0) {
				opserr << "Invalid mat tags for E_SFI " << iData[0] << endln;
				return 0;
			}
			for (int i = 0; i < m; i++) {
				theMaterials[i] = 0;
				theMaterials[i] = OPS_getNDMaterial(matTags[i]);
				if (theMaterials[i] == 0) {
					opserr << "Invalid material tag " << matTags[i] << " for E_SFI " << iData[0] << endln;
					return 0;
				}
			}
		}

		numArgs = OPS_GetNumRemainingInputArgs();

	}

	theElement = new E_SFI(iData[0], iData[1], iData[2], theMaterials, theThickness, theWidth, iData[3], dData[0]);

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
E_SFI::E_SFI(int tag,
	int Nd1, int Nd2,
	NDMaterial** materials,
	double* thickness,
	double* width,
	int mm,
	double cc)
	:Element(tag, ELE_TAG_E_SFI),
	externalNodes(2),
	theNd1(0),
	theNd2(0),
	theNodesALL(0),
	theMaterial(0), theLoad(0),
	E_SFIStrainX(0), E_SFIStrainY(0), E_SFIStrainXY(0), E_SFIStrain(0),
	x(0), b(0), AcX(0), AcY(0), kx(0), ky(0), kh(0), Fx(0), Fy(0), Fxy(0), Dens(0), Dx(0), Dy(0), Dxy(0),
	E_SFIK(6, 6), E_SFIR(6), E_SFID(6, 6), E_SFIM(6, 6),
	P_6DOF(6), Dsh(0),
	m(mm), c(cc)

{
	// Fill with ZEROs all element matrices
	E_SFIK.Zero();		    
	E_SFIR.Zero();		  
	P_6DOF.Zero();			

	TotalMass = 0.0;
	NodeMass = 0.0;
	h = 0.0;

	externalNodes(0) = Nd1;
	externalNodes(1) = Nd2;

	// Set external node pointers to NULL - external nodes
	theNodes[0] = 0;
	theNodes[1] = 0;

	theNodesALL = new Node * [2];

	// Set theNodesALL pointers to NULL
	for (int i = 0; i < 2; i++) {
		theNodesALL[i] = 0;
	}

	// Check thickness and width input
	if (thickness == 0) {
		opserr << "E_SFI::E_SFI() - Null thickness array passed.\n";
		exit(-1);
	}

	if (width == 0) {
		opserr << "E_SFI::E_SFI() - Null width array passed.\n";
		exit(-1);
	}

	// Allocate memory for the thickness and width
	t = new double[m];
	b = new double[m];
	Lw = 0.0;

	for (int i = 0; i < m; i++) {
		t[i] = thickness[i];
		b[i] = width[i];
		Lw += b[i]; // Total length of the wall
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
		opserr << "E_SFI::E_SFI() - Null material array passed.\n";
		exit(-1);
	}

	// Allocate memory for the ND materials
	theMaterial = new NDMaterial * [m];

	if (theMaterial == 0) {
		opserr << "E_SFI::E_SFI() - Failed to allocate pointers for uniaxial materials.\n";
		exit(-1);
	}

	// Get copies of the ND materials
	for (int i = 0; i < m; i++) {
		if (materials[i] == 0) {
			opserr << "E_SFI::E_SFI() - Null ND material pointer passed.\n";
			exit(-1);
		}

		theMaterial[i] = materials[i]->getCopy("PlaneStress2D");

		if (theMaterial[i] == 0) {
			opserr << "E_SFI::E_SFI() - Failed to copy ND material.\n";
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
	kh = new double[1];

	// Panel force (trial)
	Fx = new double[m];
	Fy = new double[m];
	Fxy = new double[m];

	// Panel stiffness (trial)
	Dx = new double[m];
	Dy = new double[m];
	Dxy = new double[m];

	// Panel strains
	E_SFIStrainX = new double[m];
	E_SFIStrainY = new double[m];
	E_SFIStrainXY = new double[m];
	E_SFIStrain = new double[3 * m];

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

		E_SFIStrainX[i] = 0.0;
		E_SFIStrainY[i] = 0.0;
		E_SFIStrainXY[i] = 0.0;

		E_SFIStrain[i] = 0.0;
		E_SFIStrain[i + m] = 0.0;
		E_SFIStrain[i + 2 * m] = 0.0;

		Dens[i] = 0.0;
	}

	kh[0] = 0.0;

	// Calculate concrete areas in X and Y directions
	for (int i = 0; i < m; i++) {
		AcX[i] = h * t[i];
		AcY[i] = b[i] * t[i];
	}

	// Get panel density from 2-D materials
	for (int i = 0; i < m; i++) {
		Dens[i] = theMaterial[i]->getRho();
	}

}

// Constructor which should be invoked by an FE_ObjectBroker only
E_SFI::E_SFI()
	:Element(0, ELE_TAG_E_SFI),
	externalNodes(2),
	theNd1(0),
	theNd2(0),
	theNodesALL(0),
	theMaterial(0), theLoad(0),
	E_SFIStrainX(0), E_SFIStrainY(0), E_SFIStrainXY(0), E_SFIStrain(0),
	x(0), b(0), AcX(0), AcY(0), kx(0), ky(0), kh(0), Fx(0), Fy(0), Fxy(0), Dens(0), Dx(0), Dy(0), Dxy(0),
	E_SFIK(6, 6), E_SFIR(6), E_SFID(6, 6), E_SFIM(6, 6),
	P_6DOF(6), Dsh(0),
	m(0), c(0)
{

	theNodes[0] = 0;
	theNodes[1] = 0;

	theNodesALL = new Node * [2];	

	// Set theNodesALL pointers to zero
	for (int i = 0; i < 2; i++)
	{
		theNodesALL[i] = 0;
	}

	E_SFIK.Zero();		// element stiffness matrix
	E_SFIR.Zero();		// element force vector (6)
	P_6DOF.Zero();		// element force vector (6)
	E_SFID.Zero();      // element damping matrix
	E_SFIM.Zero();		// element mass matrix

}

//  Destructor - provided to clean up any memory 
E_SFI::~E_SFI()
{
	// clean up the memory associated with the element, this is
	// memory the E_SFI objects allocates and memory allocated 
	// by other objects that the E_SFI object is responsible for 
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
	if (t != 0)
		delete[]t;
	if (AcX != 0)
		delete[]AcX;
	if (AcY != 0)
		delete[]AcY;
	if (kx != 0)
		delete[]kx;
	if (ky != 0)
		delete[]ky;
	if (kh != 0)
		delete[]kh;
	if (Fx != 0)
		delete[]Fx;
	if (Fy != 0)
		delete[]Fy;
	if (Fxy != 0)
		delete[]Fxy;
	if (Dens != 0)
		delete[]Dens;
	if (Dx != 0)
		delete[]Dx;
	if (Dy != 0)
		delete[]Dy;
	if (Dxy != 0)
		delete[]Dxy;
	if (E_SFIStrainX != 0)
		delete[]E_SFIStrainX;
	if (E_SFIStrainY != 0)
		delete[]E_SFIStrainY;
	if (E_SFIStrainXY != 0)
		delete[]E_SFIStrainXY;
	if (E_SFIStrain != 0)
		delete[]E_SFIStrain;
	if (theNodesALL != 0)
		delete[] theNodesALL;
}

// Get number of nodes
int E_SFI::getNumExternalNodes(void) const
{
	return 2;
}

// Get node tags
const ID& E_SFI::getExternalNodes(void)
{
	return externalNodes;
}

// Get shear deformation
double E_SFI::getShearDef(void)
{
	return Dsh;
}

// Get curvature (from vertical strains)
double E_SFI::getCurvature(void)
{
	double Curv;

	Curv = (E_SFIStrainY[0] - E_SFIStrainY[m - 1]) / (x[0] - x[m - 1]);

	return Curv;
}

// Get global forces at 6 DOFs (top and bottom node)
Vector E_SFI::getResistingForce_6DOF(void)
{
	for (int i = 0; i < 6; i++) {
		P_6DOF(i) = E_SFIR(i);
	}

	return P_6DOF;
}

// Get node pointers
Node** E_SFI::getNodePtrs(void)
{

	return theNodes;
}

// Get number of DOFs 
int E_SFI::getNumDOF(void) {

	int NumDOF = 6; 

	return NumDOF;
}

// Set Domain
void E_SFI::setDomain(Domain* theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
	if (theDomain == 0)
	{
		return;
	}

	// Set node pointers to NULL
	theNodes[0] = 0;
	theNodes[1] = 0;

	// First ensure nodes (external) exist in Domain and set the node pointers
	int Nd1 = externalNodes(0);
	int Nd2 = externalNodes(1);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);

	// Get coordinates of end nodes
	const Vector& end1Crd = theNodes[0]->getCrds();
	const Vector& end2Crd = theNodes[1]->getCrds();

	// Calculate the element height and perform checks
	h = end2Crd(1) - end1Crd(1);

	if (h < 0.0) {
		opserr << "WARNING: Element height is negative. Define Nodes from bottom to top!";
		return;
	}

	if (h == 0.0) {
		opserr << "WARNING: Element height is ZERO!";
		return;
	}

	// Calculate concrete areas in X and Y directions
	for (int i = 0; i < m; i++) {
		AcX[i] = h * t[i];
	}

	// Currently element can be only vertical
	if (end2Crd(0) != end1Crd(0)) {
		opserr << "WARNING: Element is NOT vertical!";
	}	

	if (theNodes[0] == 0)
	{
		opserr << "WARNING E_SFI::setDomain() - at E_SFI " << this->getTag() << " node " << Nd1 << " does not exist in domain\n";
		return;  // Don't go any further - otherwise segemntation fault
	}
	if (theNodes[1] == 0)
	{
		opserr << "WARNING E_SFI::setDomain() - at E_SFI " << this->getTag() << " node " << Nd2 << " does not exist in domain\n";
		return;
	}

	// Call the DomainComponent class method 
	this->DomainComponent::setDomain(theDomain);

	// Ensure connected nodes have correct number of dof's
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();
	if ((dofNd1 != 3) || (dofNd2 != 3)) // 3 DOFS at nodes
	{
		opserr << "E_SFI::setDomain(): 3 dof required at nodes, " << dofNd1 << " and " << dofNd2 << " provided\n";
	}

	// Calculate the nodal mass (external nodes only) for lumped mass approach
	for (int i = 0; i < m; i++) {
		TotalMass += Dens[i] * AcY[i] * h;
	}

	NodeMass = TotalMass / 2.0;

	// Create a vector to hop applied loads - NOT used in the current model formulation (no element loads)
	if (theLoad == 0) { theLoad = new Vector(6); }
		
	if (theLoad == 0) {
		opserr << "E_SFI::setDomain() - element: " << this->getTag() << " out of memory creating vector of size: " << 6 << endln;
		return;
	}

}

// Commit state
int E_SFI::commitState()
{
	int errCode = 0;

	// Commit material models
	for (int i = 0; i < m; i++)

		errCode += theMaterial[i]->commitState();

	return errCode;
}

// Revert to last committed state (if convergence is not achieved)
int E_SFI::revertToLastCommit()
{
	int errCode = 0;

	// Revert material models
	for (int i = 0; i < m; i++)

		errCode += theMaterial[i]->revertToLastCommit();

	return errCode;
}

// Revert to start
int E_SFI::revertToStart()
{

	int errCode = 0;

	// Revert material models
	for (int i = 0; i < m; i++)
		errCode += theMaterial[i]->revertToStart();

	// Compute initial stiffness
	this->getInitialStiff();

	return errCode;

}

// Update state
int E_SFI::update()
{

	// Get the current strain given trial displacements at nodes
	this->computeCurrentStrain();

	// Set the strain in the materials
	int errCode1 = 0;

	for (int i = 0; i < m; i++) {

		Vector strain(3);

		strain(0) = E_SFIStrain[i];
		strain(1) = E_SFIStrain[i + m];
		strain(2) = E_SFIStrain[i + 2 * m];

		// Set trial response for material models
		errCode1 += theMaterial[i]->setTrialStrain(strain);
	}
	return errCode1;
}

// Send Self
int E_SFI::sendSelf(int commitTag, Channel& theChannel)
{
	int res;
	int dataTag = this->getDbTag();

	static Vector data(3);  // One bigger than needed so no clash later

	data(0) = this->getTag();
	data(1) = m;
	data(2) = c;

	// E_SFI then sends the tags of it's nodes
	res = theChannel.sendID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING E_SFI::sendSelf() - failed to send ID\n";
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
int E_SFI::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
	int res;
	int dataTag = this->getDbTag();

	// E_SFI creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	// delete dynamic memory
	if (theMaterial != 0) {
		for (int i = 0; i < m; i++)
			if (theMaterial[i] != 0)
				delete theMaterial[i];
		delete[] theMaterial;
	}

	Vector data(3); // One bigger than needed so no clash later
	res = theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING E_SFI::recvSelf() - failed to receive Vector\n";
		return -1;
	}

	this->setTag((int)data(0));

	data(0) = this->getTag();
	data(1) = m;
	data(2) = c;

	// E_SFI now receives the tags of it's two external nodes
	res = theChannel.recvID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING E_SFI::recvSelf() - failed to receive ID\n";
		return -2;
	}

	// Receive the material class tags
	ID matClassTags(m);
	res = theChannel.recvID(0, commitTag, matClassTags);

	// Allocate memory for the uniaxial materials
	theMaterial = new NDMaterial * [m];
	if (theMaterial == 0) {
		opserr << "E_SFI::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Receive the material models
	for (int i = 0; i < m; i++) {
		theMaterial[i] = theBroker.getNewNDMaterial(matClassTags(i));
		if (theMaterial[i] == 0) {
			opserr << "E_SFI::recvSelf() - "
				<< "failed to get blank uniaxial material.\n";
			return -3;
		}
		theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
	}
	return 0;
}

// Display model
int E_SFI::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
	// Get the end points of the beam for the display factor
	static Vector v1(3);
	static Vector v2(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);

	// determine the deformation - rotation - other is taken from v1, v2
	static Vector r1(1);
	theNodes[0]->getDisplayRots(r1, fact, displayMode);

	// Displaying wall axis
	int error = 0;
	Vector RGB(3);
	RGB(0) = 0.0;
	RGB(1) = 1.0;
	RGB(2) = 0.0;
	error += theViewer.drawLine(v1, v2, RGB, RGB, 1, 1);

	// Displaying Panels
	for (int panel = 0; panel < m; panel++) // loop over m panels
	{
		Matrix NodePLotCrds(m, 13); // (panel id, x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)

		// first set the quantity to be displayed at the nodes;
		// if displayMode is 1 through 3 we will plot material stresses otherwise 0.0
		static Vector values(1); // values to be plotted (either epsX, epsY, gammaXY)
		if (displayMode < 4 && displayMode > 0) {
			const Vector& stress = theMaterial[panel]->getStrain();
			values(0) = stress(displayMode - 1);
		}
		else {
			values(0) = 0.0;
		}

		// Panel nodes
		NodePLotCrds(panel, 0) = panel + 1; // panel id
		// Local node 1 - bottom left
		NodePLotCrds(panel, 1) = v1(0) + x[panel] - b[panel] / 2.0; // x 
		NodePLotCrds(panel, 2) = v1(1) + (x[panel] - b[panel] / 2.0) * r1(0); // y
		NodePLotCrds(panel, 3) = v1(2); // z
		// Local node 2 - bottom right
		NodePLotCrds(panel, 4) = v1(0) + x[panel] + b[panel] / 2.0; // x
		NodePLotCrds(panel, 5) = v1(1) + (x[panel] + b[panel] / 2.0) * r1(0); // y
		NodePLotCrds(panel, 6) = v1(2); // z
		// Local node 3 - top left
		NodePLotCrds(panel, 7) = v2(0) + x[panel] + b[panel] / 2.0; // x
		NodePLotCrds(panel, 8) = v2(1) + (x[panel] + b[panel] / 2.0) * r1(0); // y
		NodePLotCrds(panel, 9) = v2(2); // z
		// Local node 4 - top right
		NodePLotCrds(panel, 10) = v2(0) + x[panel] - b[panel] / 2.0; // x
		NodePLotCrds(panel, 11) = v2(1) + (x[panel] - b[panel] / 2.0) * r1(0); // y
		NodePLotCrds(panel, 12) = v2(2); // z

		Matrix coords(4, 3); // Temporary coordinates for plotting

		coords(0, 0) = NodePLotCrds(panel, 1); // node 1 x
		coords(1, 0) = NodePLotCrds(panel, 4); // node 2 x
		coords(2, 0) = NodePLotCrds(panel, 7); // node 3 x
		coords(3, 0) = NodePLotCrds(panel, 10); // node 4 x

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
void E_SFI::Print(OPS_Stream& s, int flag)
{
	if (flag == 0) {
		s << "E_SFI Element tag: " << this->getTag() << endln;
		s << "iNode: " << externalNodes(0) << ", jNode: " << externalNodes(1) << endln;
		s << "Element height: " << h << endln;
		s << "Number of RC panel elements: " << m << endln;

		// get resisting forces in global system
		s << "Global resisting forces: " << this->getResistingForce_6DOF();

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
Response* E_SFI::setResponse(const char** argv, int argc, OPS_Stream& s)

{
	Response* theResponse = 0;

	s.tag("ElementOutput");
	s.attr("eleType", "E_SFI");
	s.attr("eleTag", this->getTag());
	s.attr("node1", externalNodes[0]);
	s.attr("node2", externalNodes[1]);

	// Global forces
	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

		s.tag("ResponseType", "Fx_i");
		s.tag("ResponseType", "Fy_i");
		s.tag("ResponseType", "Mz_i");
		s.tag("ResponseType", "Fx_j");
		s.tag("ResponseType", "Fy_j");
		s.tag("ResponseType", "Mz_j");

		return theResponse = new ElementResponse(this, 1, Vector(6));

	}

	// Shear deformation
	else if (strcmp(argv[0], "ShearDef") == 0 || strcmp(argv[0], "sheardef") == 0) {

		s.tag("ResponseType", "Dsh");

		return theResponse = new ElementResponse(this, 2, 0.0);

	}

	// Element curvature
	else if (strcmp(argv[0], "Curvature") == 0 || strcmp(argv[0], "curvature") == 0) {

		s.tag("ResponseType", "fi");

		return theResponse = new ElementResponse(this, 3, 0.0);
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

	return 0;
}

// Obtain element responses
int E_SFI::getResponse(int responseID, Information& eleInfo)
{
	switch (responseID)
	{
	case 1:  // Global forces
		return eleInfo.setVector(this->getResistingForce_6DOF());

	case 2:  // Shear deformation
		return eleInfo.setDouble(this->getShearDef());

	case 3:  // Curvature
		return eleInfo.setDouble(this->getCurvature());

	default:

		return 0;

	}
}

// Get the element initial element tangent matrix
const Matrix& E_SFI::getInitialStiff(void)
{
	double Kh = 0.0;

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
	double Kv = 0.0; double Km = 0.0; double e = 0.0; double ex = 0.0;

	for (int i = 0; i < m; ++i)
	{
		Kv += ky[i];
		Km += ky[i] * x[i] * x[i];
		e += ky[i] * x[i];
	}

	E_SFIK(0, 0) = Kh;
	E_SFIK(0, 1) = 0.0;
	E_SFIK(0, 2) = -Kh * c * h;
	E_SFIK(0, 3) = -Kh;
	E_SFIK(0, 4) = 0.0;
	E_SFIK(0, 5) = -Kh * (1 - c) * h;

	E_SFIK(1, 0) = E_SFIK(0, 1);
	E_SFIK(1, 1) = Kv;
	E_SFIK(1, 2) = e;
	E_SFIK(1, 3) = 0.0;
	E_SFIK(1, 4) = -Kv;
	E_SFIK(1, 5) = -e;

	E_SFIK(2, 0) = E_SFIK(0, 2);
	E_SFIK(2, 1) = E_SFIK(1, 2);
	E_SFIK(2, 2) = h * h * c * c * Kh + Km;
	E_SFIK(2, 3) = h * c * Kh;
	E_SFIK(2, 4) = -e;
	E_SFIK(2, 5) = (1 - c) * c * h * h * Kh - Km;

	E_SFIK(3, 0) = E_SFIK(0, 3);
	E_SFIK(3, 1) = E_SFIK(1, 3);
	E_SFIK(3, 2) = E_SFIK(2, 3);
	E_SFIK(3, 3) = Kh;
	E_SFIK(3, 4) = 0.0;
	E_SFIK(3, 5) = Kh * (1 - c) * h;

	E_SFIK(4, 0) = E_SFIK(0, 4);
	E_SFIK(4, 1) = E_SFIK(1, 4);
	E_SFIK(4, 2) = E_SFIK(2, 4);
	E_SFIK(4, 3) = E_SFIK(3, 4);
	E_SFIK(4, 4) = Kv;
	E_SFIK(4, 5) = e;

	E_SFIK(5, 0) = E_SFIK(0, 5);
	E_SFIK(5, 1) = E_SFIK(1, 5);
	E_SFIK(5, 2) = E_SFIK(2, 5);
	E_SFIK(5, 3) = E_SFIK(3, 5);
	E_SFIK(5, 4) = E_SFIK(4, 5);
	E_SFIK(5, 5) = (1 - c) * (1 - c) * h * h * Kh + Km;

	// Return element stiffness matrix
	return E_SFIK;

}

// Get current element tangent stiffness matrix from the material for the last updated strain
const Matrix& E_SFI::getTangentStiff(void)
{

	double Kh = 0.0;

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
	double Kv = 0.0; double Km = 0.0; double e = 0.0; double ex = 0.0;

	for (int i = 0; i < m; ++i)
	{
		Kv += ky[i];
		Km += ky[i] * x[i] * x[i];
		e += ky[i] * x[i];
	}

	E_SFIK(0, 0) = Kh;
	E_SFIK(0, 1) = 0.0;
	E_SFIK(0, 2) = -Kh * c * h;
	E_SFIK(0, 3) = -Kh;
	E_SFIK(0, 4) = 0.0;
	E_SFIK(0, 5) = -Kh * (1 - c) * h;

	E_SFIK(1, 0) = E_SFIK(0, 1);
	E_SFIK(1, 1) = Kv;
	E_SFIK(1, 2) = e;
	E_SFIK(1, 3) = 0.0;
	E_SFIK(1, 4) = -Kv;
	E_SFIK(1, 5) = -e;

	E_SFIK(2, 0) = E_SFIK(0, 2);
	E_SFIK(2, 1) = E_SFIK(1, 2);
	E_SFIK(2, 2) = h * h * c * c * Kh + Km;
	E_SFIK(2, 3) = h * c * Kh;
	E_SFIK(2, 4) = -e;
	E_SFIK(2, 5) = (1 - c) * c * h * h * Kh - Km;

	E_SFIK(3, 0) = E_SFIK(0, 3);
	E_SFIK(3, 1) = E_SFIK(1, 3);
	E_SFIK(3, 2) = E_SFIK(2, 3);
	E_SFIK(3, 3) = Kh;
	E_SFIK(3, 4) = 0.0;
	E_SFIK(3, 5) = Kh * (1 - c) * h;

	E_SFIK(4, 0) = E_SFIK(0, 4);
	E_SFIK(4, 1) = E_SFIK(1, 4);
	E_SFIK(4, 2) = E_SFIK(2, 4);
	E_SFIK(4, 3) = E_SFIK(3, 4);
	E_SFIK(4, 4) = Kv;
	E_SFIK(4, 5) = e;

	E_SFIK(5, 0) = E_SFIK(0, 5);
	E_SFIK(5, 1) = E_SFIK(1, 5);
	E_SFIK(5, 2) = E_SFIK(2, 5);
	E_SFIK(5, 3) = E_SFIK(3, 5);
	E_SFIK(5, 4) = E_SFIK(4, 5);
	E_SFIK(5, 5) = (1 - c) * (1 - c) * h * h * Kh + Km;

	// Return element stiffness matrix
	return E_SFIK;

}

// Get element mass matrix assuming lumped mass
const Matrix& E_SFI::getMass(void)
{
	E_SFIM.Zero();

	// No rotational mass
	E_SFIM(0, 0) = NodeMass;
	E_SFIM(1, 1) = NodeMass;
	E_SFIM(3, 3) = NodeMass;
	E_SFIM(4, 4) = NodeMass;

	// Return element mass matrix
	return E_SFIM;
}

// Get element damping matrix
const Matrix& E_SFI::getDamp(void)
{
	E_SFID.Zero();

	E_SFID = this->Element::getDamp();

	// Return element damping matrix
	return E_SFID;
}

// N/A to this model - no element loads
void E_SFI::zeroLoad(void)
{
	// does nothing 
}

// N/A to this model - no element loads
int E_SFI::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return 0;
}

// N/A to this model - no element loads
int E_SFI::addInertiaLoadToUnbalance(const Vector& accel)
{
	return 0;
}

// Get element force vector
const Vector& E_SFI::getResistingForce()
{
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

	E_SFIR(0) = Fh;
	E_SFIR(1) = -Fysum;
	E_SFIR(2) = -Fh * c * h;
	E_SFIR(3) = -Fh;
	E_SFIR(4) = Fysum;
	E_SFIR(5) = -Fh * (1 - c) * h;

	for (int i = 0; i < m; i++) {
		E_SFIR(2) -= Fy[i] * x[i];
		E_SFIR(5) += +Fy[i] * x[i];
	}

	// Return element force vector
	return E_SFIR;
}

// Get current strains at RC panels (macro-fibers)
void E_SFI::computeCurrentStrain(void)
{
	const Vector& disp1 = theNodes[0]->getTrialDisp(); // DOFs D1,D2,D3
	const Vector& disp2 = theNodes[1]->getTrialDisp(); // DOFs D4,D5,D6

	// Deformations at each RC panel (macro-fiber)
	for (int i = 0; i < m; i++) {
		Dy[i] = -disp1(1) - x[i] * disp1(2) + disp2(1) + x[i] * disp2(2);
		Dxy[i] = disp1(0) - disp2(0) - c * h * disp1(2) - (1 - c) * h * disp2(2);
	}

	Dsh = -Dxy[0]; // Store shear deformations for the recorder

	// Strains at each RC panel (macro-fiber)
	for (int i = 0; i < m; i++) {
		E_SFIStrainY[i] = Dy[i] / h;
		E_SFIStrainXY[i] = -Dxy[i] / h;
		E_SFIStrainX[i] = 0.55 * 1.0 * (1.0 - pow(3.0, -800.0 * abs(E_SFIStrainXY[i]))) * abs(E_SFIStrainXY[i]);
	}

	// Store strains into a single vector
	for (int i = 0; i < m; i++) {
		E_SFIStrain[i] = E_SFIStrainX[i];
		E_SFIStrain[i + m] = E_SFIStrainY[i];
		E_SFIStrain[i + 2 * m] = E_SFIStrainXY[i];
	}
}

// Get resisting force increment from inertial forces
const Vector& E_SFI::getResistingForceIncInertia()
{
	// compute the current resisting force
	this->getResistingForce();

	if (TotalMass != 0.0) {

		// Get nodal accelerations
		const Vector& accel1 = theNodes[0]->getTrialAccel();
		const Vector& accel2 = theNodes[1]->getTrialAccel();

		E_SFIR(0) += NodeMass * accel1(0);
		E_SFIR(1) += NodeMass * accel1(1);
		E_SFIR(3) += NodeMass * accel2(0);
		E_SFIR(4) += NodeMass * accel2(1);

		// Add the damping forces if rayleigh damping
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			E_SFIR += this->getRayleighDampingForces();

	}
	else {

		// Add the damping forces if rayleigh damping
		if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			E_SFIR += this->getRayleighDampingForces();
	}

	return E_SFIR;
}



