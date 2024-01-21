// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								Kamiar Kalbasi
//								Kutay Orakcal
//								John Wallace
//								 
//
// Created: 03/2021
//
// Description: The MVLEM-3D model is a three-dimenaional four-node element with 24 DOFs for nonlinear analysis of 
// flexure-controlled non-rectangular reinforced concrete walls subjected to multidirectional loading. The model is 
// an extension of the two-dimensional, two-node Multiple-Vertical-Line-Element-Model (MVLEM). The baseline MVLEM, 
// which is essentially a line element for rectangular walls subjected to in-plane loading, is extended to a 
// three-dimensional model formulation by: 1) applying geometric transformation of the element in-plane degrees of 
// freedom that convert it into a four-node element formulation and by 2) incorporating linear 
// elastic out-of-plane behavior based on the Kirchhoff plate theory. The in-plane and the out-of-plane 
// element behaviors are uncoupled in the present model formulation.
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
// K. Kolozvari, K. Kalbasi, K. Orakcal & J. W. Wallace (2021), "Three-dimensional model for nonlinear analysis of slender flanged reinforced concrete walls", Engineering Structures.
//
// Source: /usr/local/cvs/OpenSees/SRC/element/mvlem/MVLEM_3D.cpp
//
// Rev: 1.0

#include "MVLEM_3D.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Globals.h>

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
Matrix MVLEM_3D::MVLEM_3DK(24, 24);
Matrix MVLEM_3D::MVLEM_3DM(24, 24);
Matrix MVLEM_3D::MVLEM_3DD(24, 24);
Vector MVLEM_3D::MVLEM_3DR(24);

// local coordinate system
Matrix MVLEM_3D::MVLEM_3DKlocal(24, 24);
Matrix MVLEM_3D::MVLEM_3DDlocal(24, 24);
Matrix MVLEM_3D::MVLEM_3DMlocal(24, 24);
Vector MVLEM_3D::MVLEM_3DRlocal(24);

#include <elementAPI.h>

// Read input parameters and build the material
void* OPS_MVLEM_3D(void)
{
	// Pointer to a uniaxial material that will be returned                       
	Element* theElement = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numArgs < 18) { 
		opserr << "Want: MVLEM_3D eleTag iNode jNode kNode lNode m -thick {fiberThick} -width {fiberWidth} -rho {Rho} -matConcrete {matTagsConcrete} -matSteel {matTagsSteel} -matShear {matTagShear} <-CoR c> <-ThickMod tMod> <-Poisson nu> <-Density Dens>\n";
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
		opserr << "WARNING invalid tag for element MVLEM_3D" << endln;
		return 0;
	}

	numData = 5;
	if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
		opserr << "WARNING iNode jNode kNode lNode or m for element MVLEM_3D" << iData[0] << endln;
		return 0;
	}

	int m = iData[5];
	const char* str = 0;

	double* theThickness = new double[m];
	double* theWidth = new double[m];
	double* theRho = new double[m];
	int* matTags = new int[m];

	UniaxialMaterial** theMaterialsConcrete = new UniaxialMaterial * [m];
	UniaxialMaterial** theMaterialsSteel = new UniaxialMaterial * [m];
	UniaxialMaterial** theMaterialsShear = new UniaxialMaterial * [1];
	for (int i = 0; i < m; i++) {
		theMaterialsConcrete[i] = 0;
		theMaterialsSteel[i] = 0;
	}
	theMaterialsShear[0] = 0;

	numArgs = OPS_GetNumRemainingInputArgs();
	while (numArgs > 0) {
		//OPS_GetStringCopy(&str);
		str = OPS_GetString();
		if (strcmp(str, "-thick") == 0) {
			numData = m;
			if (OPS_GetDoubleInput(&numData, theThickness) != 0) {
				opserr << "Invalid thick parameter for MVLEM   " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-width") == 0) {
			numData = m;
			if (OPS_GetDoubleInput(&numData, theWidth) != 0) {
				opserr << "Invalid width value for MVLEM  " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-rho") == 0) {
			numData = m;
			if (OPS_GetDoubleInput(&numData, theRho) != 0) {
				opserr << "Invalid rho value for MVLEM  " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-matConcrete") == 0) {
			numData = m;
			if (OPS_GetIntInput(&numData, matTags) != 0) {
				opserr << "Invalid concrete tags for MVLEM  " << iData[0] << endln;
				return 0;
			}
			for (int i = 0; i < m; i++) {
				theMaterialsConcrete[i] = 0;
				theMaterialsConcrete[i] = OPS_getUniaxialMaterial(matTags[i]);
				if (theMaterialsConcrete[i] == 0) {
					opserr << "Invalid material tag " << matTags[i] << "  for MVLEM  " << iData[0] << endln;
					return 0;
				}
			}
		}
		else if (strcmp(str, "-matSteel") == 0) {
			numData = m;
			if (OPS_GetIntInput(&numData, matTags) != 0) {
				opserr << "Invalid steel tags for MVLEM  " << iData[0] << endln;
				return 0;
			}
			for (int i = 0; i < m; i++) {
				theMaterialsSteel[i] = 0;
				theMaterialsSteel[i] = OPS_getUniaxialMaterial(matTags[i]);
				if (theMaterialsSteel[i] == 0) {
					opserr << "Invalid material tag " << matTags[i] << "  for MVLEM  " << iData[0] << endln;
					return 0;
				}
			}
		}
		else if (strcmp(str, "-matShear") == 0) {
			numData = 1;
			if (OPS_GetIntInput(&numData, matTags) != 0) {
				opserr << "Invalid shear tags for MVLEM  " << iData[0] << endln;
				return 0;
			}
			for (int i = 0; i < 1; i++) {
				theMaterialsShear[i] = 0;
				theMaterialsShear[i] = OPS_getUniaxialMaterial(matTags[i]);
				if (theMaterialsShear[i] == 0) {
					opserr << "Invalid material tag " << matTags[i] << "  for MVLEM  " << iData[0] << endln;
					return 0;
				}
			}
		}
		// optional parameters
		else if (strcmp(str, "-CoR") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
				opserr << "Invalid CoR parameter for MVLEM   " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-ThickMod") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
				opserr << "Invalid thickMod parameter for MVLEM   " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-Poisson") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[2]) != 0) {
				opserr << "Invalid Poisson parameter for MVLEM   " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-Density") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[3]) != 0) {
				opserr << "Invalid Dens parameter for MVLEM   " << iData[0] << endln;
				return 0;
			}
		}
		numArgs = OPS_GetNumRemainingInputArgs();

	}

	theElement = new MVLEM_3D(iData[0], dData[3], 
		iData[1], iData[2], iData[3], iData[4],
		theMaterialsConcrete, theMaterialsSteel, theMaterialsShear, theRho, theThickness, 
		theWidth, iData[5], dData[0], dData[2], dData[1]);

	// Cleanup dynamic memory
	if (theThickness != 0)
		delete[] theThickness;
	if (theWidth != 0)
		delete[] theWidth;
	if (theRho != 0)
		delete[] theRho;
	if (matTags != 0)
		delete[] matTags;

	if (theMaterialsConcrete != 0)
		delete[] theMaterialsConcrete;
	if (theMaterialsSteel != 0)
		delete[] theMaterialsSteel;
	if (theMaterialsShear != 0)
		delete[] theMaterialsShear;

	return theElement;
}


// typical constructor
MVLEM_3D::MVLEM_3D(int tag,
	double Dens,
	int Nd1, int Nd2, int Nd3, int Nd4,
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

	:Element(tag, ELE_TAG_MVLEM_3D),
	density(Dens),
	externalNodes(4),
	theNd1(0),
	theNd2(0),
	theNd3(0),
	theNd4(0),
	theMaterialsConcrete(0), theMaterialsSteel(0), theMaterialsShear(0),
	theLoad(0), MVLEM_3DStrain(0),
	c(cc), m(mm), NUelastic(nn), Tfactor(tf), Eave(0.0),
	T(24, 24), T6(6, 6), Tt(3, 3)
{
	// Fill with ZEROs all element matrices
	MVLEM_3DK.Zero();
	MVLEM_3DR.Zero();
	MVLEM_3DD.Zero();
	MVLEM_3DM.Zero();

	MVLEM_3DKlocal.Zero();
	MVLEM_3DRlocal.Zero();
	MVLEM_3DDlocal.Zero();
	MVLEM_3DMlocal.Zero();

	NodeMass = 0.0;
	h = 0.0;
	d = 0.0;
	Lw = 0.0;
	Tave = 0.0;
	Eave = 0.0;

	// Fill in the ID containing external node info with node id's    
	if (externalNodes.Size() != 4)
		opserr << "FATAL MVLEM_3D::MVLEM_3D() - out of memory, could not create an ID of size 4\n";

	// Assign node tags to external nodes - Node ordering switched on purpose to match the theoretical derivation
	externalNodes(0) = Nd1;
	externalNodes(1) = Nd2;
	externalNodes(3) = Nd3;	
	externalNodes(2) = Nd4; 

	//Set node pointers to NULL
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;

	// Check thickness and width input
	if (thickness == 0) {
		opserr << "MVLEM_3D::MVLEM_3D() - "
			<< "Null thickness array passed.\n";
		exit(-1);
	}

	if (width == 0) {
		opserr << "MVLEM_3D::MVLEM_3D() - "
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
	MVLEM_3DStrain = new double[m + 1];

	// Assign zero to element arrays
	for (int i = 0; i < m; i++) {

		Ac[i] = 0.0;
		As[i] = 0.0;

		ky[i] = 0.0;

		stressC[i] = 0.0;
		stressS[i] = 0.0;

		Ec[i] = 0.0;
		Es[i] = 0.0;

		MVLEM_3DStrain[i] = 0.0;
	}

	MVLEM_3DStrain[m] = 0.0;

	kh[0] = 0.0;

	// Check Concrete material input
	if (materialsConcrete == 0) {
		opserr << "MVLEM_3D::MVLEM_3D() - "
			<< "null Concrete material array passed.\n";
		exit(-1);
	}

	// Check Steel material input
	if (materialsSteel == 0) {
		opserr << "MVLEM_3D::MVLEM_3D() - "
			<< "null Steel material array passed.\n";
		exit(-1);
	}

	// Check Shear material input
	if (materialsShear == 0) {
		opserr << "MVLEM_3D::MVLEM_3D() - "
			<< "null Shear material passed.\n";
		exit(-1);
	}

	// Allocate memory for the Concrete uniaxial materials
	theMaterialsConcrete = new UniaxialMaterial*[m];
	if (theMaterialsConcrete == 0) {
		opserr << "MVLEM_3D::MVLEM_3D() - "
			<< "failed to allocate pointers for Concrete uniaxial materials.\n";
		exit(-1);
	}

	// Get copies of the Concrete uniaxial materials
	for (int i = 0; i < m; i++) {
		if (materialsConcrete[i] == 0) {
			opserr << "MVLEM_3D::MVLEM_3D() - "
				"null uniaxial Concrete material pointer passed.\n";
			exit(-1);
		}

		theMaterialsConcrete[i] = materialsConcrete[i]->getCopy();

		if (theMaterialsConcrete[i] == 0) {
			opserr << "MVLEM_3D::MVLEM_3D() - "
				<< "failed to copy Concrete uniaxial material.\n";
			exit(-1);
		}
	}

	// Allocate memory for the Steel uniaxial materials
	theMaterialsSteel = new UniaxialMaterial*[m];
	if (theMaterialsSteel == 0) {
		opserr << "MVLEM_3D::MVLEM_3D() - "
			<< "failed to allocate pointers for Steel uniaxial materials.\n";
		exit(-1);
	}

	// Get copies of the uniaxial materials
	for (int i = 0; i < m; i++) {
		if (materialsSteel[i] == 0) {
			opserr << "MVLEM_3D::MVLEM_3D() - "
				"null uniaxial Steel material pointer passed.\n";
			exit(-1);
		}

		theMaterialsSteel[i] = materialsSteel[i]->getCopy();

		if (theMaterialsSteel[i] == 0) {
			opserr << "MVLEM_3D::MVLEM_3D() - "
				<< "failed to copy Steel uniaxial material.\n";
			exit(-1);
		}
	}

	// Allocate memory for the Shear uniaxial materials
	theMaterialsShear = new UniaxialMaterial*[1];
	if (theMaterialsShear == 0) {
		opserr << "MVLEM_3D::MVLEM_3D() - "
			<< "failed to allocate pointers for Shear uniaxial materials.\n";
		exit(-1);
	}

	// Get copies of the uniaxial materials
	for (int i = 0; i < 1; i++) {
		if (materialsShear[i] == 0) {
			opserr << "MVLEM_3D::MVLEM_3D() - "
				"null uniaxial Shear material pointer passed.\n";
			exit(-1);
		}

		theMaterialsShear[i] = materialsShear[i]->getCopy();

		if (theMaterialsShear[i] == 0) {
			opserr << "MVLEM_3D::MVLEM_3D() - "
				<< "failed to copy Shear uniaxial material.\n";
			exit(-1);
		}
	}

	// Revert to start
	this->revertToStart();
}

// Constructor which should be invoked by an FE_ObjectBroker only
MVLEM_3D::MVLEM_3D()
	:Element(0, ELE_TAG_MVLEM_3D),
	density(0.0),
	externalNodes(4),
	theNd1(0),
	theNd2(0),
	theNd3(0),
	theNd4(0),
	theMaterialsConcrete(0), theMaterialsSteel(0), theMaterialsShear(0),
	theLoad(0), MVLEM_3DStrain(0),
	c(0.0), m(0.0), NUelastic(0.0), Tfactor(0.0), Eave(0.0),
	T(24, 24), T6(6, 6), Tt(3, 3)

{
	if (externalNodes.Size() != 4)
		opserr << "FATAL MVLEM_3D::MVLEM_3D() - out of memory, could not create an ID of size 2\n";
	theNodes[0] = 0;
	theNodes[1] = 0;
	theNodes[2] = 0;
	theNodes[3] = 0;
}

//  Destructor - provided to clean up any memory
MVLEM_3D::~MVLEM_3D()
{
	// clean up the memory associated with the element, this is
	// memory the MVLEM_3D objects allocates and memory allocated 
	// by other objects that the MVLEM_3D object is responsible for 
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
		delete []x;
	if (t != 0)
		delete []t;
	if (b != 0)
		delete []b;
	if (rho != 0)
		delete []rho;
	if (Ac != 0)
		delete []Ac;
	if (As != 0)
		delete []As;
	if (ky != 0)
		delete []ky;
	if (kh != 0)
		delete []kh;
	if (Ec != 0)
		delete []Ec;
	if (Es != 0)
		delete []Es;
	if (stressC != 0)
		delete []stressC;
	if (stressS != 0)
		delete []stressS;
	if (MVLEM_3DStrain != 0)
		delete []MVLEM_3DStrain;
}

int MVLEM_3D::getNumExternalNodes(void) const
{
	return 4;
}

const ID &MVLEM_3D::getExternalNodes(void)
{
	return externalNodes;
}

Node **MVLEM_3D::getNodePtrs(void)
{
	return theNodes;
}

int MVLEM_3D::getNumDOF(void) {
	return 24;
}

void MVLEM_3D::setDomain(Domain *theDomain)
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

	// get node pointers
	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);
	theNodes[2] = theDomain->getNode(Nd3);
	theNodes[3] = theDomain->getNode(Nd4);

	// Get coordinates of end nodes 
	nd1Crds = theNodes[0]->getCrds();
	nd2Crds = theNodes[1]->getCrds();
	nd3Crds = theNodes[2]->getCrds();
	nd4Crds = theNodes[3]->getCrds();

	if (theNodes[0] == 0) {
		opserr << "WARNING MVLEM_3D::setDomain() - at MVLEM_3D " << this->getTag() << " node " <<
			Nd1 << " does not exist in domain\n";
		return;  // don't go any further - otherwise segmentation fault
	}

	if (theNodes[1] == 0) {
		opserr << "WARNING MVLEM_3D::setDomain() - at MVLEM_3D " << this->getTag() << " node " <<
			Nd2 << " does not exist in domain\n";
		return;
	}

	if (theNodes[2] == 0) {
		opserr << "WARNING MVLEM_3D::setDomain() - at MVLEM_3D " << this->getTag() << " node " <<
			Nd3 << " does not exist in domain\n";
		return;
	}

	if (theNodes[3] == 0) {
		opserr << "WARNING MVLEM_3D::setDomain() - at MVLEM_3D " << this->getTag() << " node " <<
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
		opserr << "MVLEM_3D::setDomain(): 6 dof required at all nodes. " << dofNd1 << " provided at node 1, " << dofNd2 << " provided at node 2, "
			<< dofNd3 << " provided at node 3, " << dofNd4 << " provided at node 4";

	}

	// Calculate the element height based on distance between top and bottom nodes and perform checks
	double h1 = pow(pow(nd3Crds(0) - nd1Crds(0), 2.0) + pow(nd3Crds(1) - nd1Crds(1), 2.0) + pow(nd3Crds(2) - nd1Crds(2), 2.0), 0.5);
	double h2 = pow(pow(nd4Crds(0) - nd2Crds(0), 2.0) + pow(nd4Crds(1) - nd2Crds(1), 2.0) + pow(nd4Crds(2) - nd2Crds(2), 2.0), 0.5);

	// Check if element height is zero
	if ((h1 == 0.0) || (h2 == 0.0)) {
		opserr << "WARNING: MVLEM_3D element with tag " << this->getTag() <<
			" has ZERO height. Check geometry.";
		exit(-1);
	}

	// Check if element has constant height
	if ((h1 / h2 > 1.01) || (h1 / h2 < 0.99)) {
		opserr << "WARNING: MVLEM_3D element with tag " << this->getTag() << 
			" does not have constant height. Heights of the element are " << h1 << " and " << h2 << ". Check geometry.";
		exit(-1);
	}

	// Element height
	h = (h1 + h2) / 2.0;

	// Calculate average wall thickness (for out-of-plane behavior)
	for (int i = 0; i<m; i++) {
		Tave += t[i] * b[i] / Lw;
	}

	// Calculate the element length based on distance between left and right nodes and perform checks
	double L1 = pow(pow(nd1Crds(0) - nd2Crds(0), 2.0) + pow(nd1Crds(1) - nd2Crds(1), 2.0) + pow(nd1Crds(2) - nd2Crds(2), 2.0), 0.5);
	double L2 = pow(pow(nd4Crds(0) - nd3Crds(0), 2.0) + pow(nd4Crds(1) - nd3Crds(1), 2.0) + pow(nd4Crds(2) - nd3Crds(2), 2.0), 0.5);

	// Check the element width
	if ((L1 / L2 > 1.01) || (L1 / L2 < 0.99)) {
		opserr << "WARNING: MVLEM_3D element with tag " << this->getTag() <<
			" does not have constant length. Top and bottom lengths of the element are " << L1 << " and " << L2 << ". Check geometry.";
		exit(-1);
	}
	
	if ((Lw / L1 > 1.01) || (Lw / L1 < 0.99)) {
		opserr << "WARNING: Node coordinates do not match sum of fiber widths for MVLEM_3D element with tag " << this->getTag() <<
			". Element width based on model geometry is " << L1 << " and sum of fiber widths is " << Lw << ". Check input and geometry.";
		exit(-1);
	}

	if ((Lw / L2 > 1.01) || (Lw / L2 < 0.99)) {
		opserr << "WARNING: Node coordinates do not match sum of fiber widths for MVLEM_3D element with tag " << this->getTag() <<
			". Element width based on model geometry is " << L2 << " and sum of fiber widths is " << Lw << ". Check input and geometry.";
		exit(-1);
	} //*/

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

	// Internal beam parameters
	Eib = Eave;
	Hib = h;		
	Aib = Tave * Hib;
	Iib = 0.5 * (Tave * Hib*Hib*Hib / 12.0);
	
	Tave *= Tfactor; // multiply Tave with the modification factor

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

	// Determine the transformation matrix
	setTransformationMatrix();

	// Create a vector to hop applied loads
	if (theLoad == 0)
		theLoad = new Vector(24);
	if (theLoad == 0) {
		opserr << "MVLEM_3D::setDomain() - element: " << this->getTag()
			<< " out of memory creating vector of size: " << 24 << endln;
		return;
	}
}

// Commit state of the materials
int MVLEM_3D::commitState() {

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
int MVLEM_3D::revertToLastCommit() {

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
int MVLEM_3D::revertToStart() {

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
int MVLEM_3D::update() {

	// Determine the current strain given trial displacements at nodes
	MVLEM_3DStrain = this->computeCurrentStrain();

	// Set the strain in the materials
	int errCode1 = 0;

	// Set trial response for Concrete material models
	for (int i = 0; i < m; i++)
		errCode1 += theMaterialsConcrete[i]->setTrialStrain(MVLEM_3DStrain[i]);

	// Set trial response for Steel material models
	for (int i = 0; i < m; i++)
		errCode1 += theMaterialsSteel[i]->setTrialStrain(MVLEM_3DStrain[i]);

	// Set trial response for Shear material model
	errCode1 += theMaterialsShear[0]->setTrialStrain(MVLEM_3DStrain[m]);

	return errCode1;
}

double * MVLEM_3D::computeCurrentStrain(void)
{

	// get nodal displacements in global cs
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

	// transform nodal displacements from global to local cs
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
		MVLEM_3DStrain[i] = (-dispL_inPlan2N(1) - x[i] * dispL_inPlan2N(2) + dispL_inPlan2N(4) + x[i] * dispL_inPlan2N(5)) / h;	
	}

	// Shear deformation
	MVLEM_3DStrain[m] = dispL_inPlan2N(0) - dispL_inPlan2N(3) - c * h*dispL_inPlan2N(2) - (1.0 - c)*h*dispL_inPlan2N(5);

	return MVLEM_3DStrain;

}

// Get initial stiffness matrix
const Matrix& MVLEM_3D::getInitialStiff(void)
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

	// Assemble element stiffness matrix 
	MVLEM_3DKlocal(0, 0) = Kh / 4.0 + (Aib * Eib) / Lw;
	MVLEM_3DKlocal(0, 1) = (Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 2) = 0.0;
	MVLEM_3DKlocal(0, 3) = 0.0;
	MVLEM_3DKlocal(0, 4) = 0.0;
	MVLEM_3DKlocal(0, 5) = -(Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 6) = Kh / 4.0 - (Aib * Eib) / Lw;
	MVLEM_3DKlocal(0, 7) = -(Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 8) = 0.0;
	MVLEM_3DKlocal(0, 9) = 0.0;
	MVLEM_3DKlocal(0, 10) = 0.0;
	MVLEM_3DKlocal(0, 11) = -(Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 12) = -Kh / 4.0;
	MVLEM_3DKlocal(0, 13) = -(Kh * d * h * (c - 1)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 14) = 0.0;
	MVLEM_3DKlocal(0, 15) = 0.0;
	MVLEM_3DKlocal(0, 16) = 0.0;
	MVLEM_3DKlocal(0, 17) = (Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 18) = -Kh / 4.0;
	MVLEM_3DKlocal(0, 19) = (Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 20) = 0.0;
	MVLEM_3DKlocal(0, 21) = 0.0;
	MVLEM_3DKlocal(0, 22) = 0.0;
	MVLEM_3DKlocal(0, 23) = (Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);

	MVLEM_3DKlocal(1, 0) = MVLEM_3DKlocal(0, 1);
	MVLEM_3DKlocal(1, 1) = Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) - (d * (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	MVLEM_3DKlocal(1, 2) = 0.0;
	MVLEM_3DKlocal(1, 3) = 0.0;
	MVLEM_3DKlocal(1, 4) = 0.0;
	MVLEM_3DKlocal(1, 5) = (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(1, 6) = (Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(1, 7) = Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) + (d * (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) - (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	MVLEM_3DKlocal(1, 8) = 0.0;
	MVLEM_3DKlocal(1, 9) = 0.0;
	MVLEM_3DKlocal(1, 10) = 0.0;
	MVLEM_3DKlocal(1, 11) = (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(1, 12) = -(Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(1, 13) = (d * e) / (4.0 * (d * d) + 4.0) - Kv / 4.0 + (d * (e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(1, 14) = 0.0;
	MVLEM_3DKlocal(1, 15) = 0.0;
	MVLEM_3DKlocal(1, 16) = 0.0;
	MVLEM_3DKlocal(1, 17) = -(e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(1, 18) = -(Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(1, 19) = (d * e) / (4.0 * (d * d) + 4.0) - Kv / 4.0 - (d * (e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(1, 20) = 0.0;
	MVLEM_3DKlocal(1, 21) = 0.0;
	MVLEM_3DKlocal(1, 22) = 0.0;
	MVLEM_3DKlocal(1, 23) = -(e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);

	MVLEM_3DKlocal(2, 0) = MVLEM_3DKlocal(0, 2);
	MVLEM_3DKlocal(2, 1) = MVLEM_3DKlocal(1, 2);
	MVLEM_3DKlocal(2, 2) = K1;
	MVLEM_3DKlocal(2, 3) = -K2;
	MVLEM_3DKlocal(2, 4) = K3;
	MVLEM_3DKlocal(2, 5) = 0.0;
	MVLEM_3DKlocal(2, 6) = 0.0;
	MVLEM_3DKlocal(2, 7) = 0.0;
	MVLEM_3DKlocal(2, 8) = K4;
	MVLEM_3DKlocal(2, 9) = K5;
	MVLEM_3DKlocal(2, 10) = K6;
	MVLEM_3DKlocal(2, 11) = 0.0;
	MVLEM_3DKlocal(2, 12) = 0.0;
	MVLEM_3DKlocal(2, 13) = 0.0;
	MVLEM_3DKlocal(2, 14) = K7;
	MVLEM_3DKlocal(2, 15) = -K8;
	MVLEM_3DKlocal(2, 16) = -K9;
	MVLEM_3DKlocal(2, 17) = 0.0;
	MVLEM_3DKlocal(2, 18) = 0.0;
	MVLEM_3DKlocal(2, 19) = 0.0;
	MVLEM_3DKlocal(2, 20) = K10;
	MVLEM_3DKlocal(2, 21) = -K11;
	MVLEM_3DKlocal(2, 22) = K12;
	MVLEM_3DKlocal(2, 23) = 0.0;

	MVLEM_3DKlocal(3, 0) = MVLEM_3DKlocal(0, 3);
	MVLEM_3DKlocal(3, 1) = MVLEM_3DKlocal(1, 3);
	MVLEM_3DKlocal(3, 2) = MVLEM_3DKlocal(2, 3);
	MVLEM_3DKlocal(3, 3) = K13;
	MVLEM_3DKlocal(3, 4) = K14;
	MVLEM_3DKlocal(3, 5) = 0.0;
	MVLEM_3DKlocal(3, 6) = 0.0;
	MVLEM_3DKlocal(3, 7) = 0.0;
	MVLEM_3DKlocal(3, 8) = K5;
	MVLEM_3DKlocal(3, 9) = K15;
	MVLEM_3DKlocal(3, 10) = 0.0;
	MVLEM_3DKlocal(3, 11) = 0.0;
	MVLEM_3DKlocal(3, 12) = 0.0;
	MVLEM_3DKlocal(3, 13) = 0.0;
	MVLEM_3DKlocal(3, 14) = K8;
	MVLEM_3DKlocal(3, 15) = K16;
	MVLEM_3DKlocal(3, 16) = 0.0;
	MVLEM_3DKlocal(3, 17) = 0.0;
	MVLEM_3DKlocal(3, 18) = 0.0;
	MVLEM_3DKlocal(3, 19) = 0.0;
	MVLEM_3DKlocal(3, 20) = K11;
	MVLEM_3DKlocal(3, 21) = K17;
	MVLEM_3DKlocal(3, 22) = 0.0;
	MVLEM_3DKlocal(3, 23) = 0.0;

	MVLEM_3DKlocal(4, 0) = MVLEM_3DKlocal(0, 4);
	MVLEM_3DKlocal(4, 1) = MVLEM_3DKlocal(1, 4);
	MVLEM_3DKlocal(4, 2) = MVLEM_3DKlocal(2, 4);
	MVLEM_3DKlocal(4, 3) = MVLEM_3DKlocal(3, 4);
	MVLEM_3DKlocal(4, 4) = K18;
	MVLEM_3DKlocal(4, 5) = 0.0;
	MVLEM_3DKlocal(4, 6) = 0.0;
	MVLEM_3DKlocal(4, 7) = 0.0;
	MVLEM_3DKlocal(4, 8) = -K6;
	MVLEM_3DKlocal(4, 9) = 0.0;
	MVLEM_3DKlocal(4, 10) = K19;
	MVLEM_3DKlocal(4, 11) = 0.0;
	MVLEM_3DKlocal(4, 12) = 0.0;
	MVLEM_3DKlocal(4, 13) = 0.0;
	MVLEM_3DKlocal(4, 14) = -K9;
	MVLEM_3DKlocal(4, 15) = 0.0;
	MVLEM_3DKlocal(4, 16) = K20;
	MVLEM_3DKlocal(4, 17) = 0.0;
	MVLEM_3DKlocal(4, 18) = 0.0;
	MVLEM_3DKlocal(4, 19) = 0.0;
	MVLEM_3DKlocal(4, 20) = -K12;
	MVLEM_3DKlocal(4, 21) = 0.0;
	MVLEM_3DKlocal(4, 22) = K21;
	MVLEM_3DKlocal(4, 23) = 0.0;

	MVLEM_3DKlocal(5, 0) = MVLEM_3DKlocal(0, 5);
	MVLEM_3DKlocal(5, 1) = MVLEM_3DKlocal(1, 5);
	MVLEM_3DKlocal(5, 2) = MVLEM_3DKlocal(2, 5);
	MVLEM_3DKlocal(5, 3) = MVLEM_3DKlocal(3, 5);
	MVLEM_3DKlocal(5, 4) = MVLEM_3DKlocal(4, 5);
	MVLEM_3DKlocal(5, 5) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	MVLEM_3DKlocal(5, 6) = -(Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(5, 7) = e / (4.0 * (d * d) + 4.0) + (d * (Km + Kh * (c * c) * (h * h))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(5, 8) = 0.0;
	MVLEM_3DKlocal(5, 9) = 0.0;
	MVLEM_3DKlocal(5, 10) = 0.0;
	MVLEM_3DKlocal(5, 11) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (2.0 * Eib * Iib) / Lw;
	MVLEM_3DKlocal(5, 12) = (Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(5, 13) = (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - e / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(5, 14) = 0.0;
	MVLEM_3DKlocal(5, 15) = 0.0;
	MVLEM_3DKlocal(5, 16) = 0.0;
	MVLEM_3DKlocal(5, 17) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(5, 18) = (Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(5, 19) = -e / (4.0 * (d * d) + 4.0) - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(5, 20) = 0.0;
	MVLEM_3DKlocal(5, 21) = 0.0;
	MVLEM_3DKlocal(5, 22) = 0.0;
	MVLEM_3DKlocal(5, 23) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));

	MVLEM_3DKlocal(6, 0) = MVLEM_3DKlocal(0, 6);
	MVLEM_3DKlocal(6, 1) = MVLEM_3DKlocal(1, 6);
	MVLEM_3DKlocal(6, 2) = MVLEM_3DKlocal(2, 6);
	MVLEM_3DKlocal(6, 3) = MVLEM_3DKlocal(3, 6);
	MVLEM_3DKlocal(6, 4) = MVLEM_3DKlocal(4, 6);
	MVLEM_3DKlocal(6, 5) = MVLEM_3DKlocal(5, 6);
	MVLEM_3DKlocal(6, 6) = Kh / 4.0 + (Aib * Eib) / Lw;
	MVLEM_3DKlocal(6, 7) = -(Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 8) = 0.0;
	MVLEM_3DKlocal(6, 9) = 0.0;
	MVLEM_3DKlocal(6, 10) = 0.0;
	MVLEM_3DKlocal(6, 11) = -(Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 12) = -Kh / 4.0;
	MVLEM_3DKlocal(6, 13) = -(Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 14) = 0.0;
	MVLEM_3DKlocal(6, 15) = 0.0;
	MVLEM_3DKlocal(6, 16) = 0.0;
	MVLEM_3DKlocal(6, 17) = (Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 18) = -Kh / 4;
	MVLEM_3DKlocal(6, 19) = (Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 20) = 0.0;
	MVLEM_3DKlocal(6, 21) = 0.0;
	MVLEM_3DKlocal(6, 22) = 0.0;
	MVLEM_3DKlocal(6, 23) = (Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);

	MVLEM_3DKlocal(7, 0) = MVLEM_3DKlocal(0, 7);
	MVLEM_3DKlocal(7, 1) = MVLEM_3DKlocal(1, 7);
	MVLEM_3DKlocal(7, 2) = MVLEM_3DKlocal(2, 7);
	MVLEM_3DKlocal(7, 3) = MVLEM_3DKlocal(3, 7);
	MVLEM_3DKlocal(7, 4) = MVLEM_3DKlocal(4, 7);
	MVLEM_3DKlocal(7, 5) = MVLEM_3DKlocal(5, 7);
	MVLEM_3DKlocal(7, 6) = MVLEM_3DKlocal(6, 7);
	MVLEM_3DKlocal(7, 7) = Kv / 4.0 + (d * e) / (4.0 * (d * d) + 4.0) + (d * (e / 2.0 + (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	MVLEM_3DKlocal(7, 8) = 0.0;
	MVLEM_3DKlocal(7, 9) = 0.0;
	MVLEM_3DKlocal(7, 10) = 0.0;
	MVLEM_3DKlocal(7, 11) = (e / 2.0 + (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) - (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(7, 12) = (Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(7, 13) = (d * (e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) - (d * e) / (4.0 * (d * d) + 4.0) - Kv / 4.0;
	MVLEM_3DKlocal(7, 14) = 0.0;
	MVLEM_3DKlocal(7, 15) = 0.0;
	MVLEM_3DKlocal(7, 16) = 0.0;
	MVLEM_3DKlocal(7, 17) = -(e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(7, 18) = (Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(7, 19) = -Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) - (d * (e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(7, 20) = 0.0;
	MVLEM_3DKlocal(7, 21) = 0.0;
	MVLEM_3DKlocal(7, 22) = 0.0;
	MVLEM_3DKlocal(7, 23) = -(e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);

	MVLEM_3DKlocal(8, 0) = MVLEM_3DKlocal(0, 8);
	MVLEM_3DKlocal(8, 1) = MVLEM_3DKlocal(1, 8);
	MVLEM_3DKlocal(8, 2) = MVLEM_3DKlocal(2, 8);
	MVLEM_3DKlocal(8, 3) = MVLEM_3DKlocal(3, 8);
	MVLEM_3DKlocal(8, 4) = MVLEM_3DKlocal(4, 8);
	MVLEM_3DKlocal(8, 5) = MVLEM_3DKlocal(5, 8);
	MVLEM_3DKlocal(8, 6) = MVLEM_3DKlocal(6, 8);
	MVLEM_3DKlocal(8, 7) = MVLEM_3DKlocal(7, 8);
	MVLEM_3DKlocal(8, 8) = K1;
	MVLEM_3DKlocal(8, 9) = -K2;
	MVLEM_3DKlocal(8, 10) = -K3;
	MVLEM_3DKlocal(8, 11) = 0.0;
	MVLEM_3DKlocal(8, 12) = 0.0;
	MVLEM_3DKlocal(8, 13) = 0.0;
	MVLEM_3DKlocal(8, 14) = K10;
	MVLEM_3DKlocal(8, 15) = -K11;
	MVLEM_3DKlocal(8, 16) = -K12;
	MVLEM_3DKlocal(8, 17) = 0.0;
	MVLEM_3DKlocal(8, 18) = 0.0;
	MVLEM_3DKlocal(8, 19) = 0.0;
	MVLEM_3DKlocal(8, 20) = K7;
	MVLEM_3DKlocal(8, 21) = -K8;
	MVLEM_3DKlocal(8, 22) = K9;
	MVLEM_3DKlocal(8, 23) = 0.0;

	MVLEM_3DKlocal(9, 0) = MVLEM_3DKlocal(0, 9);
	MVLEM_3DKlocal(9, 1) = MVLEM_3DKlocal(1, 9);
	MVLEM_3DKlocal(9, 2) = MVLEM_3DKlocal(2, 9);
	MVLEM_3DKlocal(9, 3) = MVLEM_3DKlocal(3, 9);
	MVLEM_3DKlocal(9, 4) = MVLEM_3DKlocal(4, 9);
	MVLEM_3DKlocal(9, 5) = MVLEM_3DKlocal(5, 9);
	MVLEM_3DKlocal(9, 6) = MVLEM_3DKlocal(6, 9);
	MVLEM_3DKlocal(9, 7) = MVLEM_3DKlocal(7, 9);
	MVLEM_3DKlocal(9, 8) = MVLEM_3DKlocal(8, 9);
	MVLEM_3DKlocal(9, 9) = K13;
	MVLEM_3DKlocal(9, 10) = -K14;
	MVLEM_3DKlocal(9, 11) = 0.0;
	MVLEM_3DKlocal(9, 12) = 0.0;
	MVLEM_3DKlocal(9, 13) = 0.0;
	MVLEM_3DKlocal(9, 14) = K11;
	MVLEM_3DKlocal(9, 15) = K17;
	MVLEM_3DKlocal(9, 16) = 0.0;
	MVLEM_3DKlocal(9, 17) = 0.0;
	MVLEM_3DKlocal(9, 18) = 0.0;
	MVLEM_3DKlocal(9, 19) = 0.0;
	MVLEM_3DKlocal(9, 20) = K8;
	MVLEM_3DKlocal(9, 21) = K16;
	MVLEM_3DKlocal(9, 22) = 0.0;
	MVLEM_3DKlocal(9, 23) = 0.0;

	MVLEM_3DKlocal(10, 0) = MVLEM_3DKlocal(0, 10);
	MVLEM_3DKlocal(10, 1) = MVLEM_3DKlocal(1, 10);
	MVLEM_3DKlocal(10, 2) = MVLEM_3DKlocal(2, 10);
	MVLEM_3DKlocal(10, 3) = MVLEM_3DKlocal(3, 10);
	MVLEM_3DKlocal(10, 4) = MVLEM_3DKlocal(4, 10);
	MVLEM_3DKlocal(10, 5) = MVLEM_3DKlocal(5, 10);
	MVLEM_3DKlocal(10, 6) = MVLEM_3DKlocal(6, 10);
	MVLEM_3DKlocal(10, 7) = MVLEM_3DKlocal(7, 10);
	MVLEM_3DKlocal(10, 8) = MVLEM_3DKlocal(8, 10);
	MVLEM_3DKlocal(10, 9) = MVLEM_3DKlocal(9, 10);
	MVLEM_3DKlocal(10, 10) = K22;
	MVLEM_3DKlocal(10, 11) = 0.0;
	MVLEM_3DKlocal(10, 12) = 0.0;
	MVLEM_3DKlocal(10, 13) = 0.0;
	MVLEM_3DKlocal(10, 14) = K12;
	MVLEM_3DKlocal(10, 15) = 0.0;
	MVLEM_3DKlocal(10, 16) = K21;
	MVLEM_3DKlocal(10, 17) = 0.0;
	MVLEM_3DKlocal(10, 18) = 0.0;
	MVLEM_3DKlocal(10, 19) = 0.0;
	MVLEM_3DKlocal(10, 20) = K9;
	MVLEM_3DKlocal(10, 21) = 0.0;
	MVLEM_3DKlocal(10, 22) = K20;
	MVLEM_3DKlocal(10, 23) = 0.0;

	MVLEM_3DKlocal(11, 0) = MVLEM_3DKlocal(0, 11);
	MVLEM_3DKlocal(11, 1) = MVLEM_3DKlocal(1, 11);
	MVLEM_3DKlocal(11, 2) = MVLEM_3DKlocal(2, 11);
	MVLEM_3DKlocal(11, 3) = MVLEM_3DKlocal(3, 11);
	MVLEM_3DKlocal(11, 4) = MVLEM_3DKlocal(4, 11);
	MVLEM_3DKlocal(11, 5) = MVLEM_3DKlocal(5, 11);
	MVLEM_3DKlocal(11, 6) = MVLEM_3DKlocal(6, 11);
	MVLEM_3DKlocal(11, 7) = MVLEM_3DKlocal(7, 11);
	MVLEM_3DKlocal(11, 8) = MVLEM_3DKlocal(8, 11);
	MVLEM_3DKlocal(11, 9) = MVLEM_3DKlocal(9, 11);
	MVLEM_3DKlocal(11, 10) = MVLEM_3DKlocal(10, 11);
	MVLEM_3DKlocal(11, 11) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	MVLEM_3DKlocal(11, 12) = (Kh * c * h) / (4 * (d * d) + 4);
	MVLEM_3DKlocal(11, 13) = (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - e / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(11, 14) = 0.0;
	MVLEM_3DKlocal(11, 15) = 0.0;
	MVLEM_3DKlocal(11, 16) = 0.0;
	MVLEM_3DKlocal(11, 17) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(11, 18) = (Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(11, 19) = -e / (4.0 * (d * d) + 4.0) - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(11, 20) = 0.0;
	MVLEM_3DKlocal(11, 21) = 0.0;
	MVLEM_3DKlocal(11, 22) = 0.0;
	MVLEM_3DKlocal(11, 23) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));

	MVLEM_3DKlocal(12, 0) = MVLEM_3DKlocal(0, 12);
	MVLEM_3DKlocal(12, 1) = MVLEM_3DKlocal(1, 12);
	MVLEM_3DKlocal(12, 2) = MVLEM_3DKlocal(2, 12);
	MVLEM_3DKlocal(12, 3) = MVLEM_3DKlocal(3, 12);
	MVLEM_3DKlocal(12, 4) = MVLEM_3DKlocal(4, 12);
	MVLEM_3DKlocal(12, 5) = MVLEM_3DKlocal(5, 12);
	MVLEM_3DKlocal(12, 6) = MVLEM_3DKlocal(6, 12);
	MVLEM_3DKlocal(12, 7) = MVLEM_3DKlocal(7, 12);
	MVLEM_3DKlocal(12, 8) = MVLEM_3DKlocal(8, 12);
	MVLEM_3DKlocal(12, 9) = MVLEM_3DKlocal(9, 12);
	MVLEM_3DKlocal(12, 10) = MVLEM_3DKlocal(10, 12);
	MVLEM_3DKlocal(12, 11) = MVLEM_3DKlocal(11, 12);
	MVLEM_3DKlocal(12, 12) = Kh / 4.0 + (Aib * Eib) / Lw;
	MVLEM_3DKlocal(12, 13) = (Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(12, 14) = 0.0;
	MVLEM_3DKlocal(12, 15) = 0.0;
	MVLEM_3DKlocal(12, 16) = 0.0;
	MVLEM_3DKlocal(12, 17) = -(Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(12, 18) = Kh / 4.0 - (Aib * Eib) / Lw;
	MVLEM_3DKlocal(12, 19) = -(Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(12, 20) = 0.0;
	MVLEM_3DKlocal(12, 21) = 0.0;
	MVLEM_3DKlocal(12, 22) = 0.0;
	MVLEM_3DKlocal(12, 23) = -(Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);

	MVLEM_3DKlocal(13, 0) = MVLEM_3DKlocal(0, 13);
	MVLEM_3DKlocal(13, 1) = MVLEM_3DKlocal(1, 13);
	MVLEM_3DKlocal(13, 2) = MVLEM_3DKlocal(2, 13);
	MVLEM_3DKlocal(13, 3) = MVLEM_3DKlocal(3, 13);
	MVLEM_3DKlocal(13, 4) = MVLEM_3DKlocal(4, 13);
	MVLEM_3DKlocal(13, 5) = MVLEM_3DKlocal(5, 13);
	MVLEM_3DKlocal(13, 6) = MVLEM_3DKlocal(6, 13);
	MVLEM_3DKlocal(13, 7) = MVLEM_3DKlocal(7, 13);
	MVLEM_3DKlocal(13, 8) = MVLEM_3DKlocal(8, 13);
	MVLEM_3DKlocal(13, 9) = MVLEM_3DKlocal(9, 13);
	MVLEM_3DKlocal(13, 10) = MVLEM_3DKlocal(10, 13);
	MVLEM_3DKlocal(13, 11) = MVLEM_3DKlocal(11, 13);
	MVLEM_3DKlocal(13, 12) = MVLEM_3DKlocal(12, 13);
	MVLEM_3DKlocal(13, 13) = Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw) - (d * (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(13, 14) = 0.0;
	MVLEM_3DKlocal(13, 15) = 0.0;
	MVLEM_3DKlocal(13, 16) = 0.0;
	MVLEM_3DKlocal(13, 17) = (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(13, 18) = (Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(13, 19) = Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) - (12.0 * Eib * Iib) / (Lw * Lw * Lw) + (d * (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(13, 20) = 0.0;
	MVLEM_3DKlocal(13, 21) = 0.0;
	MVLEM_3DKlocal(13, 22) = 0.0;
	MVLEM_3DKlocal(13, 23) = (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);

	MVLEM_3DKlocal(14, 0) = MVLEM_3DKlocal(0, 14);
	MVLEM_3DKlocal(14, 1) = MVLEM_3DKlocal(1, 14);
	MVLEM_3DKlocal(14, 2) = MVLEM_3DKlocal(2, 14);
	MVLEM_3DKlocal(14, 3) = MVLEM_3DKlocal(3, 14);
	MVLEM_3DKlocal(14, 4) = MVLEM_3DKlocal(4, 14);
	MVLEM_3DKlocal(14, 5) = MVLEM_3DKlocal(5, 14);
	MVLEM_3DKlocal(14, 6) = MVLEM_3DKlocal(6, 14);
	MVLEM_3DKlocal(14, 7) = MVLEM_3DKlocal(7, 14);
	MVLEM_3DKlocal(14, 8) = MVLEM_3DKlocal(8, 14);
	MVLEM_3DKlocal(14, 9) = MVLEM_3DKlocal(9, 14);
	MVLEM_3DKlocal(14, 10) = MVLEM_3DKlocal(10, 14);
	MVLEM_3DKlocal(14, 11) = MVLEM_3DKlocal(11, 14);
	MVLEM_3DKlocal(14, 12) = MVLEM_3DKlocal(12, 14);
	MVLEM_3DKlocal(14, 13) = MVLEM_3DKlocal(13, 14);
	MVLEM_3DKlocal(14, 14) = K1;
	MVLEM_3DKlocal(14, 15) = K2;
	MVLEM_3DKlocal(14, 16) = K3;
	MVLEM_3DKlocal(14, 17) = 0.0;
	MVLEM_3DKlocal(14, 18) = 0.0;
	MVLEM_3DKlocal(14, 19) = 0.0;
	MVLEM_3DKlocal(14, 20) = K4;
	MVLEM_3DKlocal(14, 21) = -K5;
	MVLEM_3DKlocal(14, 22) = K6;
	MVLEM_3DKlocal(14, 23) = 0.0;

	MVLEM_3DKlocal(15, 0) = MVLEM_3DKlocal(0, 15);
	MVLEM_3DKlocal(15, 1) = MVLEM_3DKlocal(1, 15);
	MVLEM_3DKlocal(15, 2) = MVLEM_3DKlocal(2, 15);
	MVLEM_3DKlocal(15, 3) = MVLEM_3DKlocal(3, 15);
	MVLEM_3DKlocal(15, 4) = MVLEM_3DKlocal(4, 15);
	MVLEM_3DKlocal(15, 5) = MVLEM_3DKlocal(5, 15);
	MVLEM_3DKlocal(15, 6) = MVLEM_3DKlocal(6, 15);
	MVLEM_3DKlocal(15, 7) = MVLEM_3DKlocal(7, 15);
	MVLEM_3DKlocal(15, 8) = MVLEM_3DKlocal(8, 15);
	MVLEM_3DKlocal(15, 9) = MVLEM_3DKlocal(9, 15);
	MVLEM_3DKlocal(15, 10) = MVLEM_3DKlocal(10, 15);
	MVLEM_3DKlocal(15, 11) = MVLEM_3DKlocal(11, 15);
	MVLEM_3DKlocal(15, 12) = MVLEM_3DKlocal(12, 15);
	MVLEM_3DKlocal(15, 13) = MVLEM_3DKlocal(13, 15);
	MVLEM_3DKlocal(15, 14) = MVLEM_3DKlocal(14, 15);
	MVLEM_3DKlocal(15, 15) = K13;
	MVLEM_3DKlocal(15, 16) = -K14;
	MVLEM_3DKlocal(15, 17) = 0.0;
	MVLEM_3DKlocal(15, 18) = 0.0;
	MVLEM_3DKlocal(15, 19) = 0.0;
	MVLEM_3DKlocal(15, 20) = -K5;
	MVLEM_3DKlocal(15, 21) = K15;
	MVLEM_3DKlocal(15, 22) = 0.0;
	MVLEM_3DKlocal(15, 23) = 0.0;

	MVLEM_3DKlocal(16, 0) = MVLEM_3DKlocal(0, 16);
	MVLEM_3DKlocal(16, 1) = MVLEM_3DKlocal(1, 16);
	MVLEM_3DKlocal(16, 2) = MVLEM_3DKlocal(2, 16);
	MVLEM_3DKlocal(16, 3) = MVLEM_3DKlocal(3, 16);
	MVLEM_3DKlocal(16, 4) = MVLEM_3DKlocal(4, 16);
	MVLEM_3DKlocal(16, 5) = MVLEM_3DKlocal(5, 16);
	MVLEM_3DKlocal(16, 6) = MVLEM_3DKlocal(6, 16);
	MVLEM_3DKlocal(16, 7) = MVLEM_3DKlocal(7, 16);
	MVLEM_3DKlocal(16, 8) = MVLEM_3DKlocal(8, 16);
	MVLEM_3DKlocal(16, 9) = MVLEM_3DKlocal(9, 16);
	MVLEM_3DKlocal(16, 10) = MVLEM_3DKlocal(10, 16);
	MVLEM_3DKlocal(16, 11) = MVLEM_3DKlocal(11, 16);
	MVLEM_3DKlocal(16, 12) = MVLEM_3DKlocal(12, 16);
	MVLEM_3DKlocal(16, 13) = MVLEM_3DKlocal(13, 16);
	MVLEM_3DKlocal(16, 14) = MVLEM_3DKlocal(14, 16);
	MVLEM_3DKlocal(16, 15) = MVLEM_3DKlocal(15, 16);
	MVLEM_3DKlocal(16, 16) = K18;
	MVLEM_3DKlocal(16, 17) = 0.0;
	MVLEM_3DKlocal(16, 18) = 0.0;
	MVLEM_3DKlocal(16, 19) = 0.0;
	MVLEM_3DKlocal(16, 20) = -K6;
	MVLEM_3DKlocal(16, 21) = 0.0;
	MVLEM_3DKlocal(16, 22) = K19;
	MVLEM_3DKlocal(16, 23) = 0.0;

	MVLEM_3DKlocal(17, 0) = MVLEM_3DKlocal(0, 17);
	MVLEM_3DKlocal(17, 1) = MVLEM_3DKlocal(1, 17);
	MVLEM_3DKlocal(17, 2) = MVLEM_3DKlocal(2, 17);
	MVLEM_3DKlocal(17, 3) = MVLEM_3DKlocal(3, 17);
	MVLEM_3DKlocal(17, 4) = MVLEM_3DKlocal(4, 17);
	MVLEM_3DKlocal(17, 5) = MVLEM_3DKlocal(5, 17);
	MVLEM_3DKlocal(17, 6) = MVLEM_3DKlocal(6, 17);
	MVLEM_3DKlocal(17, 7) = MVLEM_3DKlocal(7, 17);
	MVLEM_3DKlocal(17, 8) = MVLEM_3DKlocal(8, 17);
	MVLEM_3DKlocal(17, 9) = MVLEM_3DKlocal(9, 17);
	MVLEM_3DKlocal(17, 10) = MVLEM_3DKlocal(10, 17);
	MVLEM_3DKlocal(17, 11) = MVLEM_3DKlocal(11, 17);
	MVLEM_3DKlocal(17, 12) = MVLEM_3DKlocal(12, 17);
	MVLEM_3DKlocal(17, 13) = MVLEM_3DKlocal(13, 17);
	MVLEM_3DKlocal(17, 14) = MVLEM_3DKlocal(14, 17);
	MVLEM_3DKlocal(17, 15) = MVLEM_3DKlocal(15, 17);
	MVLEM_3DKlocal(17, 16) = MVLEM_3DKlocal(16, 17);
	MVLEM_3DKlocal(17, 17) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	MVLEM_3DKlocal(17, 18) = -(Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(17, 19) = e / (4.0 * (d * d) + 4.0) - (6.0 * Eib * Iib) / (Lw * Lw) + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(17, 20) = 0.0;
	MVLEM_3DKlocal(17, 21) = 0.0;
	MVLEM_3DKlocal(17, 22) = 0.0;
	MVLEM_3DKlocal(17, 23) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (2.0 * Eib * Iib) / Lw;

	MVLEM_3DKlocal(18, 0) = MVLEM_3DKlocal(0, 18);
	MVLEM_3DKlocal(18, 1) = MVLEM_3DKlocal(1, 18);
	MVLEM_3DKlocal(18, 2) = MVLEM_3DKlocal(2, 18);
	MVLEM_3DKlocal(18, 3) = MVLEM_3DKlocal(3, 18);
	MVLEM_3DKlocal(18, 4) = MVLEM_3DKlocal(4, 18);
	MVLEM_3DKlocal(18, 5) = MVLEM_3DKlocal(5, 18);
	MVLEM_3DKlocal(18, 6) = MVLEM_3DKlocal(6, 18);
	MVLEM_3DKlocal(18, 7) = MVLEM_3DKlocal(7, 18);
	MVLEM_3DKlocal(18, 8) = MVLEM_3DKlocal(8, 18);
	MVLEM_3DKlocal(18, 9) = MVLEM_3DKlocal(9, 18);
	MVLEM_3DKlocal(18, 10) = MVLEM_3DKlocal(10, 18);
	MVLEM_3DKlocal(18, 11) = MVLEM_3DKlocal(11, 18);
	MVLEM_3DKlocal(18, 12) = MVLEM_3DKlocal(12, 18);
	MVLEM_3DKlocal(18, 13) = MVLEM_3DKlocal(13, 18);
	MVLEM_3DKlocal(18, 14) = MVLEM_3DKlocal(14, 18);
	MVLEM_3DKlocal(18, 15) = MVLEM_3DKlocal(15, 18);
	MVLEM_3DKlocal(18, 16) = MVLEM_3DKlocal(16, 18);
	MVLEM_3DKlocal(18, 17) = MVLEM_3DKlocal(17, 18);
	MVLEM_3DKlocal(18, 18) = Kh / 4.0 + (Aib * Eib) / Lw;
	MVLEM_3DKlocal(18, 19) = -(Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(18, 20) = 0.0;
	MVLEM_3DKlocal(18, 21) = 0.0;
	MVLEM_3DKlocal(18, 22) = 0.0;
	MVLEM_3DKlocal(18, 23) = -(Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);

	MVLEM_3DKlocal(19, 0) = MVLEM_3DKlocal(0, 19);
	MVLEM_3DKlocal(19, 1) = MVLEM_3DKlocal(1, 19);
	MVLEM_3DKlocal(19, 2) = MVLEM_3DKlocal(2, 19);
	MVLEM_3DKlocal(19, 3) = MVLEM_3DKlocal(3, 19);
	MVLEM_3DKlocal(19, 4) = MVLEM_3DKlocal(4, 19);
	MVLEM_3DKlocal(19, 5) = MVLEM_3DKlocal(5, 19);
	MVLEM_3DKlocal(19, 6) = MVLEM_3DKlocal(6, 19);
	MVLEM_3DKlocal(19, 7) = MVLEM_3DKlocal(7, 19);
	MVLEM_3DKlocal(19, 8) = MVLEM_3DKlocal(8, 19);
	MVLEM_3DKlocal(19, 9) = MVLEM_3DKlocal(9, 19);
	MVLEM_3DKlocal(19, 10) = MVLEM_3DKlocal(10, 19);
	MVLEM_3DKlocal(19, 11) = MVLEM_3DKlocal(11, 19);
	MVLEM_3DKlocal(19, 12) = MVLEM_3DKlocal(12, 19);
	MVLEM_3DKlocal(19, 13) = MVLEM_3DKlocal(13, 19);
	MVLEM_3DKlocal(19, 14) = MVLEM_3DKlocal(14, 19);
	MVLEM_3DKlocal(19, 15) = MVLEM_3DKlocal(15, 19);
	MVLEM_3DKlocal(19, 16) = MVLEM_3DKlocal(16, 19);
	MVLEM_3DKlocal(19, 17) = MVLEM_3DKlocal(17, 19);
	MVLEM_3DKlocal(19, 18) = MVLEM_3DKlocal(18, 19);
	MVLEM_3DKlocal(19, 19) = Kv / 4.0 + (d * e) / (4.0 * (d * d) + 4.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw) + (d * (e / 2.0 + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(19, 20) = 0.0;
	MVLEM_3DKlocal(19, 21) = 0.0;
	MVLEM_3DKlocal(19, 22) = 0.0;
	MVLEM_3DKlocal(19, 23) = (e / 2.0 + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) - (6.0 * Eib * Iib) / (Lw * Lw);

	MVLEM_3DKlocal(20, 0) = MVLEM_3DKlocal(0, 20);
	MVLEM_3DKlocal(20, 1) = MVLEM_3DKlocal(1, 20);
	MVLEM_3DKlocal(20, 2) = MVLEM_3DKlocal(2, 20);
	MVLEM_3DKlocal(20, 3) = MVLEM_3DKlocal(3, 20);
	MVLEM_3DKlocal(20, 4) = MVLEM_3DKlocal(4, 20);
	MVLEM_3DKlocal(20, 5) = MVLEM_3DKlocal(5, 20);
	MVLEM_3DKlocal(20, 6) = MVLEM_3DKlocal(6, 20);
	MVLEM_3DKlocal(20, 7) = MVLEM_3DKlocal(7, 20);
	MVLEM_3DKlocal(20, 8) = MVLEM_3DKlocal(8, 20);
	MVLEM_3DKlocal(20, 9) = MVLEM_3DKlocal(9, 20);
	MVLEM_3DKlocal(20, 10) = MVLEM_3DKlocal(10, 20);
	MVLEM_3DKlocal(20, 11) = MVLEM_3DKlocal(11, 20);
	MVLEM_3DKlocal(20, 12) = MVLEM_3DKlocal(12, 20);
	MVLEM_3DKlocal(20, 13) = MVLEM_3DKlocal(13, 20);
	MVLEM_3DKlocal(20, 14) = MVLEM_3DKlocal(14, 20);
	MVLEM_3DKlocal(20, 15) = MVLEM_3DKlocal(15, 20);
	MVLEM_3DKlocal(20, 16) = MVLEM_3DKlocal(16, 20);
	MVLEM_3DKlocal(20, 17) = MVLEM_3DKlocal(17, 20);
	MVLEM_3DKlocal(20, 18) = MVLEM_3DKlocal(18, 20);
	MVLEM_3DKlocal(20, 19) = MVLEM_3DKlocal(19, 20);
	MVLEM_3DKlocal(20, 20) = K1;
	MVLEM_3DKlocal(20, 21) = K2;
	MVLEM_3DKlocal(20, 22) = -K3;
	MVLEM_3DKlocal(20, 23) = 0.0;

	MVLEM_3DKlocal(21, 0) = MVLEM_3DKlocal(0, 21);
	MVLEM_3DKlocal(21, 1) = MVLEM_3DKlocal(1, 21);
	MVLEM_3DKlocal(21, 2) = MVLEM_3DKlocal(2, 21);
	MVLEM_3DKlocal(21, 3) = MVLEM_3DKlocal(3, 21);
	MVLEM_3DKlocal(21, 4) = MVLEM_3DKlocal(4, 21);
	MVLEM_3DKlocal(21, 5) = MVLEM_3DKlocal(5, 21);
	MVLEM_3DKlocal(21, 6) = MVLEM_3DKlocal(6, 21);
	MVLEM_3DKlocal(21, 7) = MVLEM_3DKlocal(7, 21);
	MVLEM_3DKlocal(21, 8) = MVLEM_3DKlocal(8, 21);
	MVLEM_3DKlocal(21, 9) = MVLEM_3DKlocal(9, 21);
	MVLEM_3DKlocal(21, 10) = MVLEM_3DKlocal(10, 21);
	MVLEM_3DKlocal(21, 11) = MVLEM_3DKlocal(11, 21);
	MVLEM_3DKlocal(21, 12) = MVLEM_3DKlocal(12, 21);
	MVLEM_3DKlocal(21, 13) = MVLEM_3DKlocal(13, 21);
	MVLEM_3DKlocal(21, 14) = MVLEM_3DKlocal(14, 21);
	MVLEM_3DKlocal(21, 15) = MVLEM_3DKlocal(15, 21);
	MVLEM_3DKlocal(21, 16) = MVLEM_3DKlocal(16, 21);
	MVLEM_3DKlocal(21, 17) = MVLEM_3DKlocal(17, 21);
	MVLEM_3DKlocal(21, 18) = MVLEM_3DKlocal(18, 21);
	MVLEM_3DKlocal(21, 19) = MVLEM_3DKlocal(19, 21);
	MVLEM_3DKlocal(21, 20) = MVLEM_3DKlocal(20, 21);
	MVLEM_3DKlocal(21, 21) = K13;
	MVLEM_3DKlocal(21, 22) = K14;
	MVLEM_3DKlocal(21, 23) = 0.0;

	MVLEM_3DKlocal(22, 0) = MVLEM_3DKlocal(0, 22);
	MVLEM_3DKlocal(22, 1) = MVLEM_3DKlocal(1, 22);
	MVLEM_3DKlocal(22, 2) = MVLEM_3DKlocal(2, 22);
	MVLEM_3DKlocal(22, 3) = MVLEM_3DKlocal(3, 22);
	MVLEM_3DKlocal(22, 4) = MVLEM_3DKlocal(4, 22);
	MVLEM_3DKlocal(22, 5) = MVLEM_3DKlocal(5, 22);
	MVLEM_3DKlocal(22, 6) = MVLEM_3DKlocal(6, 22);
	MVLEM_3DKlocal(22, 7) = MVLEM_3DKlocal(7, 22);
	MVLEM_3DKlocal(22, 8) = MVLEM_3DKlocal(8, 22);
	MVLEM_3DKlocal(22, 9) = MVLEM_3DKlocal(9, 22);
	MVLEM_3DKlocal(22, 10) = MVLEM_3DKlocal(10, 22);
	MVLEM_3DKlocal(22, 11) = MVLEM_3DKlocal(11, 22);
	MVLEM_3DKlocal(22, 12) = MVLEM_3DKlocal(12, 22);
	MVLEM_3DKlocal(22, 13) = MVLEM_3DKlocal(13, 22);
	MVLEM_3DKlocal(22, 14) = MVLEM_3DKlocal(14, 22);
	MVLEM_3DKlocal(22, 15) = MVLEM_3DKlocal(15, 22);
	MVLEM_3DKlocal(22, 16) = MVLEM_3DKlocal(16, 22);
	MVLEM_3DKlocal(22, 17) = MVLEM_3DKlocal(17, 22);
	MVLEM_3DKlocal(22, 18) = MVLEM_3DKlocal(18, 22);
	MVLEM_3DKlocal(22, 19) = MVLEM_3DKlocal(19, 22);
	MVLEM_3DKlocal(22, 20) = MVLEM_3DKlocal(20, 22);
	MVLEM_3DKlocal(22, 21) = MVLEM_3DKlocal(21, 22);
	MVLEM_3DKlocal(22, 22) = K18;
	MVLEM_3DKlocal(22, 23) = 0.0;

	MVLEM_3DKlocal(23, 0) = MVLEM_3DKlocal(0, 23);
	MVLEM_3DKlocal(23, 1) = MVLEM_3DKlocal(1, 23);
	MVLEM_3DKlocal(23, 2) = MVLEM_3DKlocal(2, 23);
	MVLEM_3DKlocal(23, 3) = MVLEM_3DKlocal(3, 23);
	MVLEM_3DKlocal(23, 4) = MVLEM_3DKlocal(4, 23);
	MVLEM_3DKlocal(23, 5) = MVLEM_3DKlocal(5, 23);
	MVLEM_3DKlocal(23, 6) = MVLEM_3DKlocal(6, 23);
	MVLEM_3DKlocal(23, 7) = MVLEM_3DKlocal(7, 23);
	MVLEM_3DKlocal(23, 8) = MVLEM_3DKlocal(8, 23);
	MVLEM_3DKlocal(23, 9) = MVLEM_3DKlocal(9, 23);
	MVLEM_3DKlocal(23, 10) = MVLEM_3DKlocal(10, 23);
	MVLEM_3DKlocal(23, 11) = MVLEM_3DKlocal(11, 23);
	MVLEM_3DKlocal(23, 12) = MVLEM_3DKlocal(12, 23);
	MVLEM_3DKlocal(23, 13) = MVLEM_3DKlocal(13, 23);
	MVLEM_3DKlocal(23, 14) = MVLEM_3DKlocal(14, 23);
	MVLEM_3DKlocal(23, 15) = MVLEM_3DKlocal(15, 23);
	MVLEM_3DKlocal(23, 16) = MVLEM_3DKlocal(16, 23);
	MVLEM_3DKlocal(23, 17) = MVLEM_3DKlocal(17, 23);
	MVLEM_3DKlocal(23, 18) = MVLEM_3DKlocal(18, 23);
	MVLEM_3DKlocal(23, 19) = MVLEM_3DKlocal(19, 23);
	MVLEM_3DKlocal(23, 20) = MVLEM_3DKlocal(20, 23);
	MVLEM_3DKlocal(23, 21) = MVLEM_3DKlocal(21, 23);
	MVLEM_3DKlocal(23, 22) = MVLEM_3DKlocal(22, 23);
	MVLEM_3DKlocal(23, 23) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;

	// Convert matrix from local to global cs 
	MVLEM_3DK.addMatrixTripleProduct(0.0, T, MVLEM_3DKlocal, 1.0);

	// Return element stiffness matrix
	return MVLEM_3DK;

}


// Get current stiffness matrix
const Matrix& MVLEM_3D::getTangentStiff(void)
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
	MVLEM_3DKlocal(0, 0) = Kh / 4.0 + (Aib * Eib) / Lw;
	MVLEM_3DKlocal(0, 1) = (Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 2) = 0.0;
	MVLEM_3DKlocal(0, 3) = 0.0;
	MVLEM_3DKlocal(0, 4) = 0.0;
	MVLEM_3DKlocal(0, 5) = -(Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 6) = Kh / 4.0 - (Aib * Eib) / Lw;
	MVLEM_3DKlocal(0, 7) = -(Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 8) = 0.0;
	MVLEM_3DKlocal(0, 9) = 0.0;
	MVLEM_3DKlocal(0, 10) = 0.0;
	MVLEM_3DKlocal(0, 11) = -(Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 12) = -Kh / 4.0;
	MVLEM_3DKlocal(0, 13) = -(Kh * d * h * (c - 1)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 14) = 0.0;
	MVLEM_3DKlocal(0, 15) = 0.0;
	MVLEM_3DKlocal(0, 16) = 0.0;
	MVLEM_3DKlocal(0, 17) = (Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 18) = -Kh / 4.0;
	MVLEM_3DKlocal(0, 19) = (Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(0, 20) = 0.0;
	MVLEM_3DKlocal(0, 21) = 0.0;
	MVLEM_3DKlocal(0, 22) = 0.0;
	MVLEM_3DKlocal(0, 23) = (Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);

	MVLEM_3DKlocal(1, 0) = MVLEM_3DKlocal(0, 1);
	MVLEM_3DKlocal(1, 1) = Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) - (d * (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	MVLEM_3DKlocal(1, 2) = 0.0;
	MVLEM_3DKlocal(1, 3) = 0.0;
	MVLEM_3DKlocal(1, 4) = 0.0;
	MVLEM_3DKlocal(1, 5) = (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(1, 6) = (Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(1, 7) = Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) + (d * (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) - (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	MVLEM_3DKlocal(1, 8) = 0.0;
	MVLEM_3DKlocal(1, 9) = 0.0;
	MVLEM_3DKlocal(1, 10) = 0.0;
	MVLEM_3DKlocal(1, 11) = (e / 2.0 - (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(1, 12) = -(Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(1, 13) = (d * e) / (4.0 * (d * d) + 4.0) - Kv / 4.0 + (d * (e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(1, 14) = 0.0;
	MVLEM_3DKlocal(1, 15) = 0.0;
	MVLEM_3DKlocal(1, 16) = 0.0;
	MVLEM_3DKlocal(1, 17) = -(e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(1, 18) = -(Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(1, 19) = (d * e) / (4.0 * (d * d) + 4.0) - Kv / 4.0 - (d * (e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(1, 20) = 0.0;
	MVLEM_3DKlocal(1, 21) = 0.0;
	MVLEM_3DKlocal(1, 22) = 0.0;
	MVLEM_3DKlocal(1, 23) = -(e / 2.0 - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);

	MVLEM_3DKlocal(2, 0) = MVLEM_3DKlocal(0, 2);
	MVLEM_3DKlocal(2, 1) = MVLEM_3DKlocal(1, 2);
	MVLEM_3DKlocal(2, 2) = K1;
	MVLEM_3DKlocal(2, 3) = -K2;
	MVLEM_3DKlocal(2, 4) = K3;
	MVLEM_3DKlocal(2, 5) = 0.0;
	MVLEM_3DKlocal(2, 6) = 0.0;
	MVLEM_3DKlocal(2, 7) = 0.0;
	MVLEM_3DKlocal(2, 8) = K4;
	MVLEM_3DKlocal(2, 9) = K5;
	MVLEM_3DKlocal(2, 10) = K6;
	MVLEM_3DKlocal(2, 11) = 0.0;
	MVLEM_3DKlocal(2, 12) = 0.0;
	MVLEM_3DKlocal(2, 13) = 0.0;
	MVLEM_3DKlocal(2, 14) = K7;
	MVLEM_3DKlocal(2, 15) = -K8;
	MVLEM_3DKlocal(2, 16) = -K9;
	MVLEM_3DKlocal(2, 17) = 0.0;
	MVLEM_3DKlocal(2, 18) = 0.0;
	MVLEM_3DKlocal(2, 19) = 0.0;
	MVLEM_3DKlocal(2, 20) = K10;
	MVLEM_3DKlocal(2, 21) = -K11;
	MVLEM_3DKlocal(2, 22) = K12;
	MVLEM_3DKlocal(2, 23) = 0.0;

	MVLEM_3DKlocal(3, 0) = MVLEM_3DKlocal(0, 3);
	MVLEM_3DKlocal(3, 1) = MVLEM_3DKlocal(1, 3);
	MVLEM_3DKlocal(3, 2) = MVLEM_3DKlocal(2, 3);
	MVLEM_3DKlocal(3, 3) = K13;
	MVLEM_3DKlocal(3, 4) = K14;
	MVLEM_3DKlocal(3, 5) = 0.0;
	MVLEM_3DKlocal(3, 6) = 0.0;
	MVLEM_3DKlocal(3, 7) = 0.0;
	MVLEM_3DKlocal(3, 8) = K5;
	MVLEM_3DKlocal(3, 9) = K15;
	MVLEM_3DKlocal(3, 10) = 0.0;
	MVLEM_3DKlocal(3, 11) = 0.0;
	MVLEM_3DKlocal(3, 12) = 0.0;
	MVLEM_3DKlocal(3, 13) = 0.0;
	MVLEM_3DKlocal(3, 14) = K8;
	MVLEM_3DKlocal(3, 15) = K16;
	MVLEM_3DKlocal(3, 16) = 0.0;
	MVLEM_3DKlocal(3, 17) = 0.0;
	MVLEM_3DKlocal(3, 18) = 0.0;
	MVLEM_3DKlocal(3, 19) = 0.0;
	MVLEM_3DKlocal(3, 20) = K11;
	MVLEM_3DKlocal(3, 21) = K17;
	MVLEM_3DKlocal(3, 22) = 0.0;
	MVLEM_3DKlocal(3, 23) = 0.0;

	MVLEM_3DKlocal(4, 0) = MVLEM_3DKlocal(0, 4);
	MVLEM_3DKlocal(4, 1) = MVLEM_3DKlocal(1, 4);
	MVLEM_3DKlocal(4, 2) = MVLEM_3DKlocal(2, 4);
	MVLEM_3DKlocal(4, 3) = MVLEM_3DKlocal(3, 4);
	MVLEM_3DKlocal(4, 4) = K18;
	MVLEM_3DKlocal(4, 5) = 0.0;
	MVLEM_3DKlocal(4, 6) = 0.0;
	MVLEM_3DKlocal(4, 7) = 0.0;
	MVLEM_3DKlocal(4, 8) = -K6;
	MVLEM_3DKlocal(4, 9) = 0.0;
	MVLEM_3DKlocal(4, 10) = K19;
	MVLEM_3DKlocal(4, 11) = 0.0;
	MVLEM_3DKlocal(4, 12) = 0.0;
	MVLEM_3DKlocal(4, 13) = 0.0;
	MVLEM_3DKlocal(4, 14) = -K9;
	MVLEM_3DKlocal(4, 15) = 0.0;
	MVLEM_3DKlocal(4, 16) = K20;
	MVLEM_3DKlocal(4, 17) = 0.0;
	MVLEM_3DKlocal(4, 18) = 0.0;
	MVLEM_3DKlocal(4, 19) = 0.0;
	MVLEM_3DKlocal(4, 20) = -K12;
	MVLEM_3DKlocal(4, 21) = 0.0;
	MVLEM_3DKlocal(4, 22) = K21;
	MVLEM_3DKlocal(4, 23) = 0.0;

	MVLEM_3DKlocal(5, 0) = MVLEM_3DKlocal(0, 5);
	MVLEM_3DKlocal(5, 1) = MVLEM_3DKlocal(1, 5);
	MVLEM_3DKlocal(5, 2) = MVLEM_3DKlocal(2, 5);
	MVLEM_3DKlocal(5, 3) = MVLEM_3DKlocal(3, 5);
	MVLEM_3DKlocal(5, 4) = MVLEM_3DKlocal(4, 5);
	MVLEM_3DKlocal(5, 5) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	MVLEM_3DKlocal(5, 6) = -(Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(5, 7) = e / (4.0 * (d * d) + 4.0) + (d * (Km + Kh * (c * c) * (h * h))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(5, 8) = 0.0;
	MVLEM_3DKlocal(5, 9) = 0.0;
	MVLEM_3DKlocal(5, 10) = 0.0;
	MVLEM_3DKlocal(5, 11) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (2.0 * Eib * Iib) / Lw;
	MVLEM_3DKlocal(5, 12) = (Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(5, 13) = (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - e / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(5, 14) = 0.0;
	MVLEM_3DKlocal(5, 15) = 0.0;
	MVLEM_3DKlocal(5, 16) = 0.0;
	MVLEM_3DKlocal(5, 17) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(5, 18) = (Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(5, 19) = -e / (4.0 * (d * d) + 4.0) - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(5, 20) = 0.0;
	MVLEM_3DKlocal(5, 21) = 0.0;
	MVLEM_3DKlocal(5, 22) = 0.0;
	MVLEM_3DKlocal(5, 23) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));

	MVLEM_3DKlocal(6, 0) = MVLEM_3DKlocal(0, 6);
	MVLEM_3DKlocal(6, 1) = MVLEM_3DKlocal(1, 6);
	MVLEM_3DKlocal(6, 2) = MVLEM_3DKlocal(2, 6);
	MVLEM_3DKlocal(6, 3) = MVLEM_3DKlocal(3, 6);
	MVLEM_3DKlocal(6, 4) = MVLEM_3DKlocal(4, 6);
	MVLEM_3DKlocal(6, 5) = MVLEM_3DKlocal(5, 6);
	MVLEM_3DKlocal(6, 6) = Kh / 4.0 + (Aib * Eib) / Lw;
	MVLEM_3DKlocal(6, 7) = -(Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 8) = 0.0;
	MVLEM_3DKlocal(6, 9) = 0.0;
	MVLEM_3DKlocal(6, 10) = 0.0;
	MVLEM_3DKlocal(6, 11) = -(Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 12) = -Kh / 4.0;
	MVLEM_3DKlocal(6, 13) = -(Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 14) = 0.0;
	MVLEM_3DKlocal(6, 15) = 0.0;
	MVLEM_3DKlocal(6, 16) = 0.0;
	MVLEM_3DKlocal(6, 17) = (Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 18) = -Kh / 4;
	MVLEM_3DKlocal(6, 19) = (Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(6, 20) = 0.0;
	MVLEM_3DKlocal(6, 21) = 0.0;
	MVLEM_3DKlocal(6, 22) = 0.0;
	MVLEM_3DKlocal(6, 23) = (Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);

	MVLEM_3DKlocal(7, 0) = MVLEM_3DKlocal(0, 7);
	MVLEM_3DKlocal(7, 1) = MVLEM_3DKlocal(1, 7);
	MVLEM_3DKlocal(7, 2) = MVLEM_3DKlocal(2, 7);
	MVLEM_3DKlocal(7, 3) = MVLEM_3DKlocal(3, 7);
	MVLEM_3DKlocal(7, 4) = MVLEM_3DKlocal(4, 7);
	MVLEM_3DKlocal(7, 5) = MVLEM_3DKlocal(5, 7);
	MVLEM_3DKlocal(7, 6) = MVLEM_3DKlocal(6, 7);
	MVLEM_3DKlocal(7, 7) = Kv / 4.0 + (d * e) / (4.0 * (d * d) + 4.0) + (d * (e / 2.0 + (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw);
	MVLEM_3DKlocal(7, 8) = 0.0;
	MVLEM_3DKlocal(7, 9) = 0.0;
	MVLEM_3DKlocal(7, 10) = 0.0;
	MVLEM_3DKlocal(7, 11) = (e / 2.0 + (d * (Km + Kh * (c * c) * (h * h))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) - (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(7, 12) = (Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(7, 13) = (d * (e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0) - (d * e) / (4.0 * (d * d) + 4.0) - Kv / 4.0;
	MVLEM_3DKlocal(7, 14) = 0.0;
	MVLEM_3DKlocal(7, 15) = 0.0;
	MVLEM_3DKlocal(7, 16) = 0.0;
	MVLEM_3DKlocal(7, 17) = -(e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(7, 18) = (Kh * c * d * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(7, 19) = -Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) - (d * (e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(7, 20) = 0.0;
	MVLEM_3DKlocal(7, 21) = 0.0;
	MVLEM_3DKlocal(7, 22) = 0.0;
	MVLEM_3DKlocal(7, 23) = -(e / 2.0 + (d * (Km + Kh * c * (h * h) * (c - 1.0))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0);

	MVLEM_3DKlocal(8, 0) = MVLEM_3DKlocal(0, 8);
	MVLEM_3DKlocal(8, 1) = MVLEM_3DKlocal(1, 8);
	MVLEM_3DKlocal(8, 2) = MVLEM_3DKlocal(2, 8);
	MVLEM_3DKlocal(8, 3) = MVLEM_3DKlocal(3, 8);
	MVLEM_3DKlocal(8, 4) = MVLEM_3DKlocal(4, 8);
	MVLEM_3DKlocal(8, 5) = MVLEM_3DKlocal(5, 8);
	MVLEM_3DKlocal(8, 6) = MVLEM_3DKlocal(6, 8);
	MVLEM_3DKlocal(8, 7) = MVLEM_3DKlocal(7, 8);
	MVLEM_3DKlocal(8, 8) = K1;
	MVLEM_3DKlocal(8, 9) = -K2;
	MVLEM_3DKlocal(8, 10) = -K3;
	MVLEM_3DKlocal(8, 11) = 0.0;
	MVLEM_3DKlocal(8, 12) = 0.0;
	MVLEM_3DKlocal(8, 13) = 0.0;
	MVLEM_3DKlocal(8, 14) = K10;
	MVLEM_3DKlocal(8, 15) = -K11;
	MVLEM_3DKlocal(8, 16) = -K12;
	MVLEM_3DKlocal(8, 17) = 0.0;
	MVLEM_3DKlocal(8, 18) = 0.0;
	MVLEM_3DKlocal(8, 19) = 0.0;
	MVLEM_3DKlocal(8, 20) = K7;
	MVLEM_3DKlocal(8, 21) = -K8;
	MVLEM_3DKlocal(8, 22) = K9;
	MVLEM_3DKlocal(8, 23) = 0.0;

	MVLEM_3DKlocal(9, 0) = MVLEM_3DKlocal(0, 9);
	MVLEM_3DKlocal(9, 1) = MVLEM_3DKlocal(1, 9);
	MVLEM_3DKlocal(9, 2) = MVLEM_3DKlocal(2, 9);
	MVLEM_3DKlocal(9, 3) = MVLEM_3DKlocal(3, 9);
	MVLEM_3DKlocal(9, 4) = MVLEM_3DKlocal(4, 9);
	MVLEM_3DKlocal(9, 5) = MVLEM_3DKlocal(5, 9);
	MVLEM_3DKlocal(9, 6) = MVLEM_3DKlocal(6, 9);
	MVLEM_3DKlocal(9, 7) = MVLEM_3DKlocal(7, 9);
	MVLEM_3DKlocal(9, 8) = MVLEM_3DKlocal(8, 9);
	MVLEM_3DKlocal(9, 9) = K13;
	MVLEM_3DKlocal(9, 10) = -K14;
	MVLEM_3DKlocal(9, 11) = 0.0;
	MVLEM_3DKlocal(9, 12) = 0.0;
	MVLEM_3DKlocal(9, 13) = 0.0;
	MVLEM_3DKlocal(9, 14) = K11;
	MVLEM_3DKlocal(9, 15) = K17;
	MVLEM_3DKlocal(9, 16) = 0.0;
	MVLEM_3DKlocal(9, 17) = 0.0;
	MVLEM_3DKlocal(9, 18) = 0.0;
	MVLEM_3DKlocal(9, 19) = 0.0;
	MVLEM_3DKlocal(9, 20) = K8;
	MVLEM_3DKlocal(9, 21) = K16;
	MVLEM_3DKlocal(9, 22) = 0.0;
	MVLEM_3DKlocal(9, 23) = 0.0;

	MVLEM_3DKlocal(10, 0) = MVLEM_3DKlocal(0, 10);
	MVLEM_3DKlocal(10, 1) = MVLEM_3DKlocal(1, 10);
	MVLEM_3DKlocal(10, 2) = MVLEM_3DKlocal(2, 10);
	MVLEM_3DKlocal(10, 3) = MVLEM_3DKlocal(3, 10);
	MVLEM_3DKlocal(10, 4) = MVLEM_3DKlocal(4, 10);
	MVLEM_3DKlocal(10, 5) = MVLEM_3DKlocal(5, 10);
	MVLEM_3DKlocal(10, 6) = MVLEM_3DKlocal(6, 10);
	MVLEM_3DKlocal(10, 7) = MVLEM_3DKlocal(7, 10);
	MVLEM_3DKlocal(10, 8) = MVLEM_3DKlocal(8, 10);
	MVLEM_3DKlocal(10, 9) = MVLEM_3DKlocal(9, 10);
	MVLEM_3DKlocal(10, 10) = K22;
	MVLEM_3DKlocal(10, 11) = 0.0;
	MVLEM_3DKlocal(10, 12) = 0.0;
	MVLEM_3DKlocal(10, 13) = 0.0;
	MVLEM_3DKlocal(10, 14) = K12;
	MVLEM_3DKlocal(10, 15) = 0.0;
	MVLEM_3DKlocal(10, 16) = K21;
	MVLEM_3DKlocal(10, 17) = 0.0;
	MVLEM_3DKlocal(10, 18) = 0.0;
	MVLEM_3DKlocal(10, 19) = 0.0;
	MVLEM_3DKlocal(10, 20) = K9;
	MVLEM_3DKlocal(10, 21) = 0.0;
	MVLEM_3DKlocal(10, 22) = K20;
	MVLEM_3DKlocal(10, 23) = 0.0;

	MVLEM_3DKlocal(11, 0) = MVLEM_3DKlocal(0, 11);
	MVLEM_3DKlocal(11, 1) = MVLEM_3DKlocal(1, 11);
	MVLEM_3DKlocal(11, 2) = MVLEM_3DKlocal(2, 11);
	MVLEM_3DKlocal(11, 3) = MVLEM_3DKlocal(3, 11);
	MVLEM_3DKlocal(11, 4) = MVLEM_3DKlocal(4, 11);
	MVLEM_3DKlocal(11, 5) = MVLEM_3DKlocal(5, 11);
	MVLEM_3DKlocal(11, 6) = MVLEM_3DKlocal(6, 11);
	MVLEM_3DKlocal(11, 7) = MVLEM_3DKlocal(7, 11);
	MVLEM_3DKlocal(11, 8) = MVLEM_3DKlocal(8, 11);
	MVLEM_3DKlocal(11, 9) = MVLEM_3DKlocal(9, 11);
	MVLEM_3DKlocal(11, 10) = MVLEM_3DKlocal(10, 11);
	MVLEM_3DKlocal(11, 11) = (Km + Kh * (c * c) * (h * h)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	MVLEM_3DKlocal(11, 12) = (Kh * c * h) / (4 * (d * d) + 4);
	MVLEM_3DKlocal(11, 13) = (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) - e / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(11, 14) = 0.0;
	MVLEM_3DKlocal(11, 15) = 0.0;
	MVLEM_3DKlocal(11, 16) = 0.0;
	MVLEM_3DKlocal(11, 17) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(11, 18) = (Kh * c * h) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(11, 19) = -e / (4.0 * (d * d) + 4.0) - (d * (Km + Kh * c * (h * h) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(11, 20) = 0.0;
	MVLEM_3DKlocal(11, 21) = 0.0;
	MVLEM_3DKlocal(11, 22) = 0.0;
	MVLEM_3DKlocal(11, 23) = -(Km + Kh * c * (h * h) * (c - 1.0)) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));

	MVLEM_3DKlocal(12, 0) = MVLEM_3DKlocal(0, 12);
	MVLEM_3DKlocal(12, 1) = MVLEM_3DKlocal(1, 12);
	MVLEM_3DKlocal(12, 2) = MVLEM_3DKlocal(2, 12);
	MVLEM_3DKlocal(12, 3) = MVLEM_3DKlocal(3, 12);
	MVLEM_3DKlocal(12, 4) = MVLEM_3DKlocal(4, 12);
	MVLEM_3DKlocal(12, 5) = MVLEM_3DKlocal(5, 12);
	MVLEM_3DKlocal(12, 6) = MVLEM_3DKlocal(6, 12);
	MVLEM_3DKlocal(12, 7) = MVLEM_3DKlocal(7, 12);
	MVLEM_3DKlocal(12, 8) = MVLEM_3DKlocal(8, 12);
	MVLEM_3DKlocal(12, 9) = MVLEM_3DKlocal(9, 12);
	MVLEM_3DKlocal(12, 10) = MVLEM_3DKlocal(10, 12);
	MVLEM_3DKlocal(12, 11) = MVLEM_3DKlocal(11, 12);
	MVLEM_3DKlocal(12, 12) = Kh / 4.0 + (Aib * Eib) / Lw;
	MVLEM_3DKlocal(12, 13) = (Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(12, 14) = 0.0;
	MVLEM_3DKlocal(12, 15) = 0.0;
	MVLEM_3DKlocal(12, 16) = 0.0;
	MVLEM_3DKlocal(12, 17) = -(Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(12, 18) = Kh / 4.0 - (Aib * Eib) / Lw;
	MVLEM_3DKlocal(12, 19) = -(Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(12, 20) = 0.0;
	MVLEM_3DKlocal(12, 21) = 0.0;
	MVLEM_3DKlocal(12, 22) = 0.0;
	MVLEM_3DKlocal(12, 23) = -(Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);

	MVLEM_3DKlocal(13, 0) = MVLEM_3DKlocal(0, 13);
	MVLEM_3DKlocal(13, 1) = MVLEM_3DKlocal(1, 13);
	MVLEM_3DKlocal(13, 2) = MVLEM_3DKlocal(2, 13);
	MVLEM_3DKlocal(13, 3) = MVLEM_3DKlocal(3, 13);
	MVLEM_3DKlocal(13, 4) = MVLEM_3DKlocal(4, 13);
	MVLEM_3DKlocal(13, 5) = MVLEM_3DKlocal(5, 13);
	MVLEM_3DKlocal(13, 6) = MVLEM_3DKlocal(6, 13);
	MVLEM_3DKlocal(13, 7) = MVLEM_3DKlocal(7, 13);
	MVLEM_3DKlocal(13, 8) = MVLEM_3DKlocal(8, 13);
	MVLEM_3DKlocal(13, 9) = MVLEM_3DKlocal(9, 13);
	MVLEM_3DKlocal(13, 10) = MVLEM_3DKlocal(10, 13);
	MVLEM_3DKlocal(13, 11) = MVLEM_3DKlocal(11, 13);
	MVLEM_3DKlocal(13, 12) = MVLEM_3DKlocal(12, 13);
	MVLEM_3DKlocal(13, 13) = Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw) - (d * (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(13, 14) = 0.0;
	MVLEM_3DKlocal(13, 15) = 0.0;
	MVLEM_3DKlocal(13, 16) = 0.0;
	MVLEM_3DKlocal(13, 17) = (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);
	MVLEM_3DKlocal(13, 18) = (Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(13, 19) = Kv / 4.0 - (d * e) / (4.0 * (d * d) + 4.0) - (12.0 * Eib * Iib) / (Lw * Lw * Lw) + (d * (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(13, 20) = 0.0;
	MVLEM_3DKlocal(13, 21) = 0.0;
	MVLEM_3DKlocal(13, 22) = 0.0;
	MVLEM_3DKlocal(13, 23) = (e / 2.0 - (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) + (6.0 * Eib * Iib) / (Lw * Lw);

	MVLEM_3DKlocal(14, 0) = MVLEM_3DKlocal(0, 14);
	MVLEM_3DKlocal(14, 1) = MVLEM_3DKlocal(1, 14);
	MVLEM_3DKlocal(14, 2) = MVLEM_3DKlocal(2, 14);
	MVLEM_3DKlocal(14, 3) = MVLEM_3DKlocal(3, 14);
	MVLEM_3DKlocal(14, 4) = MVLEM_3DKlocal(4, 14);
	MVLEM_3DKlocal(14, 5) = MVLEM_3DKlocal(5, 14);
	MVLEM_3DKlocal(14, 6) = MVLEM_3DKlocal(6, 14);
	MVLEM_3DKlocal(14, 7) = MVLEM_3DKlocal(7, 14);
	MVLEM_3DKlocal(14, 8) = MVLEM_3DKlocal(8, 14);
	MVLEM_3DKlocal(14, 9) = MVLEM_3DKlocal(9, 14);
	MVLEM_3DKlocal(14, 10) = MVLEM_3DKlocal(10, 14);
	MVLEM_3DKlocal(14, 11) = MVLEM_3DKlocal(11, 14);
	MVLEM_3DKlocal(14, 12) = MVLEM_3DKlocal(12, 14);
	MVLEM_3DKlocal(14, 13) = MVLEM_3DKlocal(13, 14);
	MVLEM_3DKlocal(14, 14) = K1;
	MVLEM_3DKlocal(14, 15) = K2;
	MVLEM_3DKlocal(14, 16) = K3;
	MVLEM_3DKlocal(14, 17) = 0.0;
	MVLEM_3DKlocal(14, 18) = 0.0;
	MVLEM_3DKlocal(14, 19) = 0.0;
	MVLEM_3DKlocal(14, 20) = K4;
	MVLEM_3DKlocal(14, 21) = -K5;
	MVLEM_3DKlocal(14, 22) = K6;
	MVLEM_3DKlocal(14, 23) = 0.0;

	MVLEM_3DKlocal(15, 0) = MVLEM_3DKlocal(0, 15);
	MVLEM_3DKlocal(15, 1) = MVLEM_3DKlocal(1, 15);
	MVLEM_3DKlocal(15, 2) = MVLEM_3DKlocal(2, 15);
	MVLEM_3DKlocal(15, 3) = MVLEM_3DKlocal(3, 15);
	MVLEM_3DKlocal(15, 4) = MVLEM_3DKlocal(4, 15);
	MVLEM_3DKlocal(15, 5) = MVLEM_3DKlocal(5, 15);
	MVLEM_3DKlocal(15, 6) = MVLEM_3DKlocal(6, 15);
	MVLEM_3DKlocal(15, 7) = MVLEM_3DKlocal(7, 15);
	MVLEM_3DKlocal(15, 8) = MVLEM_3DKlocal(8, 15);
	MVLEM_3DKlocal(15, 9) = MVLEM_3DKlocal(9, 15);
	MVLEM_3DKlocal(15, 10) = MVLEM_3DKlocal(10, 15);
	MVLEM_3DKlocal(15, 11) = MVLEM_3DKlocal(11, 15);
	MVLEM_3DKlocal(15, 12) = MVLEM_3DKlocal(12, 15);
	MVLEM_3DKlocal(15, 13) = MVLEM_3DKlocal(13, 15);
	MVLEM_3DKlocal(15, 14) = MVLEM_3DKlocal(14, 15);
	MVLEM_3DKlocal(15, 15) = K13;
	MVLEM_3DKlocal(15, 16) = -K14;
	MVLEM_3DKlocal(15, 17) = 0.0;
	MVLEM_3DKlocal(15, 18) = 0.0;
	MVLEM_3DKlocal(15, 19) = 0.0;
	MVLEM_3DKlocal(15, 20) = -K5;
	MVLEM_3DKlocal(15, 21) = K15;
	MVLEM_3DKlocal(15, 22) = 0.0;
	MVLEM_3DKlocal(15, 23) = 0.0;

	MVLEM_3DKlocal(16, 0) = MVLEM_3DKlocal(0, 16);
	MVLEM_3DKlocal(16, 1) = MVLEM_3DKlocal(1, 16);
	MVLEM_3DKlocal(16, 2) = MVLEM_3DKlocal(2, 16);
	MVLEM_3DKlocal(16, 3) = MVLEM_3DKlocal(3, 16);
	MVLEM_3DKlocal(16, 4) = MVLEM_3DKlocal(4, 16);
	MVLEM_3DKlocal(16, 5) = MVLEM_3DKlocal(5, 16);
	MVLEM_3DKlocal(16, 6) = MVLEM_3DKlocal(6, 16);
	MVLEM_3DKlocal(16, 7) = MVLEM_3DKlocal(7, 16);
	MVLEM_3DKlocal(16, 8) = MVLEM_3DKlocal(8, 16);
	MVLEM_3DKlocal(16, 9) = MVLEM_3DKlocal(9, 16);
	MVLEM_3DKlocal(16, 10) = MVLEM_3DKlocal(10, 16);
	MVLEM_3DKlocal(16, 11) = MVLEM_3DKlocal(11, 16);
	MVLEM_3DKlocal(16, 12) = MVLEM_3DKlocal(12, 16);
	MVLEM_3DKlocal(16, 13) = MVLEM_3DKlocal(13, 16);
	MVLEM_3DKlocal(16, 14) = MVLEM_3DKlocal(14, 16);
	MVLEM_3DKlocal(16, 15) = MVLEM_3DKlocal(15, 16);
	MVLEM_3DKlocal(16, 16) = K18;
	MVLEM_3DKlocal(16, 17) = 0.0;
	MVLEM_3DKlocal(16, 18) = 0.0;
	MVLEM_3DKlocal(16, 19) = 0.0;
	MVLEM_3DKlocal(16, 20) = -K6;
	MVLEM_3DKlocal(16, 21) = 0.0;
	MVLEM_3DKlocal(16, 22) = K19;
	MVLEM_3DKlocal(16, 23) = 0.0;

	MVLEM_3DKlocal(17, 0) = MVLEM_3DKlocal(0, 17);
	MVLEM_3DKlocal(17, 1) = MVLEM_3DKlocal(1, 17);
	MVLEM_3DKlocal(17, 2) = MVLEM_3DKlocal(2, 17);
	MVLEM_3DKlocal(17, 3) = MVLEM_3DKlocal(3, 17);
	MVLEM_3DKlocal(17, 4) = MVLEM_3DKlocal(4, 17);
	MVLEM_3DKlocal(17, 5) = MVLEM_3DKlocal(5, 17);
	MVLEM_3DKlocal(17, 6) = MVLEM_3DKlocal(6, 17);
	MVLEM_3DKlocal(17, 7) = MVLEM_3DKlocal(7, 17);
	MVLEM_3DKlocal(17, 8) = MVLEM_3DKlocal(8, 17);
	MVLEM_3DKlocal(17, 9) = MVLEM_3DKlocal(9, 17);
	MVLEM_3DKlocal(17, 10) = MVLEM_3DKlocal(10, 17);
	MVLEM_3DKlocal(17, 11) = MVLEM_3DKlocal(11, 17);
	MVLEM_3DKlocal(17, 12) = MVLEM_3DKlocal(12, 17);
	MVLEM_3DKlocal(17, 13) = MVLEM_3DKlocal(13, 17);
	MVLEM_3DKlocal(17, 14) = MVLEM_3DKlocal(14, 17);
	MVLEM_3DKlocal(17, 15) = MVLEM_3DKlocal(15, 17);
	MVLEM_3DKlocal(17, 16) = MVLEM_3DKlocal(16, 17);
	MVLEM_3DKlocal(17, 17) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;
	MVLEM_3DKlocal(17, 18) = -(Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(17, 19) = e / (4.0 * (d * d) + 4.0) - (6.0 * Eib * Iib) / (Lw * Lw) + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0));
	MVLEM_3DKlocal(17, 20) = 0.0;
	MVLEM_3DKlocal(17, 21) = 0.0;
	MVLEM_3DKlocal(17, 22) = 0.0;
	MVLEM_3DKlocal(17, 23) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (2.0 * Eib * Iib) / Lw;

	MVLEM_3DKlocal(18, 0) = MVLEM_3DKlocal(0, 18);
	MVLEM_3DKlocal(18, 1) = MVLEM_3DKlocal(1, 18);
	MVLEM_3DKlocal(18, 2) = MVLEM_3DKlocal(2, 18);
	MVLEM_3DKlocal(18, 3) = MVLEM_3DKlocal(3, 18);
	MVLEM_3DKlocal(18, 4) = MVLEM_3DKlocal(4, 18);
	MVLEM_3DKlocal(18, 5) = MVLEM_3DKlocal(5, 18);
	MVLEM_3DKlocal(18, 6) = MVLEM_3DKlocal(6, 18);
	MVLEM_3DKlocal(18, 7) = MVLEM_3DKlocal(7, 18);
	MVLEM_3DKlocal(18, 8) = MVLEM_3DKlocal(8, 18);
	MVLEM_3DKlocal(18, 9) = MVLEM_3DKlocal(9, 18);
	MVLEM_3DKlocal(18, 10) = MVLEM_3DKlocal(10, 18);
	MVLEM_3DKlocal(18, 11) = MVLEM_3DKlocal(11, 18);
	MVLEM_3DKlocal(18, 12) = MVLEM_3DKlocal(12, 18);
	MVLEM_3DKlocal(18, 13) = MVLEM_3DKlocal(13, 18);
	MVLEM_3DKlocal(18, 14) = MVLEM_3DKlocal(14, 18);
	MVLEM_3DKlocal(18, 15) = MVLEM_3DKlocal(15, 18);
	MVLEM_3DKlocal(18, 16) = MVLEM_3DKlocal(16, 18);
	MVLEM_3DKlocal(18, 17) = MVLEM_3DKlocal(17, 18);
	MVLEM_3DKlocal(18, 18) = Kh / 4.0 + (Aib * Eib) / Lw;
	MVLEM_3DKlocal(18, 19) = -(Kh * d * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);
	MVLEM_3DKlocal(18, 20) = 0.0;
	MVLEM_3DKlocal(18, 21) = 0.0;
	MVLEM_3DKlocal(18, 22) = 0.0;
	MVLEM_3DKlocal(18, 23) = -(Kh * h * (c - 1.0)) / (4.0 * (d * d) + 4.0);

	MVLEM_3DKlocal(19, 0) = MVLEM_3DKlocal(0, 19);
	MVLEM_3DKlocal(19, 1) = MVLEM_3DKlocal(1, 19);
	MVLEM_3DKlocal(19, 2) = MVLEM_3DKlocal(2, 19);
	MVLEM_3DKlocal(19, 3) = MVLEM_3DKlocal(3, 19);
	MVLEM_3DKlocal(19, 4) = MVLEM_3DKlocal(4, 19);
	MVLEM_3DKlocal(19, 5) = MVLEM_3DKlocal(5, 19);
	MVLEM_3DKlocal(19, 6) = MVLEM_3DKlocal(6, 19);
	MVLEM_3DKlocal(19, 7) = MVLEM_3DKlocal(7, 19);
	MVLEM_3DKlocal(19, 8) = MVLEM_3DKlocal(8, 19);
	MVLEM_3DKlocal(19, 9) = MVLEM_3DKlocal(9, 19);
	MVLEM_3DKlocal(19, 10) = MVLEM_3DKlocal(10, 19);
	MVLEM_3DKlocal(19, 11) = MVLEM_3DKlocal(11, 19);
	MVLEM_3DKlocal(19, 12) = MVLEM_3DKlocal(12, 19);
	MVLEM_3DKlocal(19, 13) = MVLEM_3DKlocal(13, 19);
	MVLEM_3DKlocal(19, 14) = MVLEM_3DKlocal(14, 19);
	MVLEM_3DKlocal(19, 15) = MVLEM_3DKlocal(15, 19);
	MVLEM_3DKlocal(19, 16) = MVLEM_3DKlocal(16, 19);
	MVLEM_3DKlocal(19, 17) = MVLEM_3DKlocal(17, 19);
	MVLEM_3DKlocal(19, 18) = MVLEM_3DKlocal(18, 19);
	MVLEM_3DKlocal(19, 19) = Kv / 4.0 + (d * e) / (4.0 * (d * d) + 4.0) + (12.0 * Eib * Iib) / (Lw * Lw * Lw) + (d * (e / 2.0 + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0))) / (2.0 * (d * d) + 2.0);
	MVLEM_3DKlocal(19, 20) = 0.0;
	MVLEM_3DKlocal(19, 21) = 0.0;
	MVLEM_3DKlocal(19, 22) = 0.0;
	MVLEM_3DKlocal(19, 23) = (e / 2.0 + (d * (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0)))) / (2.0 * (d * d) + 2.0)) / (2.0 * (d * d) + 2.0) - (6.0 * Eib * Iib) / (Lw * Lw);

	MVLEM_3DKlocal(20, 0) = MVLEM_3DKlocal(0, 20);
	MVLEM_3DKlocal(20, 1) = MVLEM_3DKlocal(1, 20);
	MVLEM_3DKlocal(20, 2) = MVLEM_3DKlocal(2, 20);
	MVLEM_3DKlocal(20, 3) = MVLEM_3DKlocal(3, 20);
	MVLEM_3DKlocal(20, 4) = MVLEM_3DKlocal(4, 20);
	MVLEM_3DKlocal(20, 5) = MVLEM_3DKlocal(5, 20);
	MVLEM_3DKlocal(20, 6) = MVLEM_3DKlocal(6, 20);
	MVLEM_3DKlocal(20, 7) = MVLEM_3DKlocal(7, 20);
	MVLEM_3DKlocal(20, 8) = MVLEM_3DKlocal(8, 20);
	MVLEM_3DKlocal(20, 9) = MVLEM_3DKlocal(9, 20);
	MVLEM_3DKlocal(20, 10) = MVLEM_3DKlocal(10, 20);
	MVLEM_3DKlocal(20, 11) = MVLEM_3DKlocal(11, 20);
	MVLEM_3DKlocal(20, 12) = MVLEM_3DKlocal(12, 20);
	MVLEM_3DKlocal(20, 13) = MVLEM_3DKlocal(13, 20);
	MVLEM_3DKlocal(20, 14) = MVLEM_3DKlocal(14, 20);
	MVLEM_3DKlocal(20, 15) = MVLEM_3DKlocal(15, 20);
	MVLEM_3DKlocal(20, 16) = MVLEM_3DKlocal(16, 20);
	MVLEM_3DKlocal(20, 17) = MVLEM_3DKlocal(17, 20);
	MVLEM_3DKlocal(20, 18) = MVLEM_3DKlocal(18, 20);
	MVLEM_3DKlocal(20, 19) = MVLEM_3DKlocal(19, 20);
	MVLEM_3DKlocal(20, 20) = K1;
	MVLEM_3DKlocal(20, 21) = K2;
	MVLEM_3DKlocal(20, 22) = -K3;
	MVLEM_3DKlocal(20, 23) = 0.0;

	MVLEM_3DKlocal(21, 0) = MVLEM_3DKlocal(0, 21);
	MVLEM_3DKlocal(21, 1) = MVLEM_3DKlocal(1, 21);
	MVLEM_3DKlocal(21, 2) = MVLEM_3DKlocal(2, 21);
	MVLEM_3DKlocal(21, 3) = MVLEM_3DKlocal(3, 21);
	MVLEM_3DKlocal(21, 4) = MVLEM_3DKlocal(4, 21);
	MVLEM_3DKlocal(21, 5) = MVLEM_3DKlocal(5, 21);
	MVLEM_3DKlocal(21, 6) = MVLEM_3DKlocal(6, 21);
	MVLEM_3DKlocal(21, 7) = MVLEM_3DKlocal(7, 21);
	MVLEM_3DKlocal(21, 8) = MVLEM_3DKlocal(8, 21);
	MVLEM_3DKlocal(21, 9) = MVLEM_3DKlocal(9, 21);
	MVLEM_3DKlocal(21, 10) = MVLEM_3DKlocal(10, 21);
	MVLEM_3DKlocal(21, 11) = MVLEM_3DKlocal(11, 21);
	MVLEM_3DKlocal(21, 12) = MVLEM_3DKlocal(12, 21);
	MVLEM_3DKlocal(21, 13) = MVLEM_3DKlocal(13, 21);
	MVLEM_3DKlocal(21, 14) = MVLEM_3DKlocal(14, 21);
	MVLEM_3DKlocal(21, 15) = MVLEM_3DKlocal(15, 21);
	MVLEM_3DKlocal(21, 16) = MVLEM_3DKlocal(16, 21);
	MVLEM_3DKlocal(21, 17) = MVLEM_3DKlocal(17, 21);
	MVLEM_3DKlocal(21, 18) = MVLEM_3DKlocal(18, 21);
	MVLEM_3DKlocal(21, 19) = MVLEM_3DKlocal(19, 21);
	MVLEM_3DKlocal(21, 20) = MVLEM_3DKlocal(20, 21);
	MVLEM_3DKlocal(21, 21) = K13;
	MVLEM_3DKlocal(21, 22) = K14;
	MVLEM_3DKlocal(21, 23) = 0.0;

	MVLEM_3DKlocal(22, 0) = MVLEM_3DKlocal(0, 22);
	MVLEM_3DKlocal(22, 1) = MVLEM_3DKlocal(1, 22);
	MVLEM_3DKlocal(22, 2) = MVLEM_3DKlocal(2, 22);
	MVLEM_3DKlocal(22, 3) = MVLEM_3DKlocal(3, 22);
	MVLEM_3DKlocal(22, 4) = MVLEM_3DKlocal(4, 22);
	MVLEM_3DKlocal(22, 5) = MVLEM_3DKlocal(5, 22);
	MVLEM_3DKlocal(22, 6) = MVLEM_3DKlocal(6, 22);
	MVLEM_3DKlocal(22, 7) = MVLEM_3DKlocal(7, 22);
	MVLEM_3DKlocal(22, 8) = MVLEM_3DKlocal(8, 22);
	MVLEM_3DKlocal(22, 9) = MVLEM_3DKlocal(9, 22);
	MVLEM_3DKlocal(22, 10) = MVLEM_3DKlocal(10, 22);
	MVLEM_3DKlocal(22, 11) = MVLEM_3DKlocal(11, 22);
	MVLEM_3DKlocal(22, 12) = MVLEM_3DKlocal(12, 22);
	MVLEM_3DKlocal(22, 13) = MVLEM_3DKlocal(13, 22);
	MVLEM_3DKlocal(22, 14) = MVLEM_3DKlocal(14, 22);
	MVLEM_3DKlocal(22, 15) = MVLEM_3DKlocal(15, 22);
	MVLEM_3DKlocal(22, 16) = MVLEM_3DKlocal(16, 22);
	MVLEM_3DKlocal(22, 17) = MVLEM_3DKlocal(17, 22);
	MVLEM_3DKlocal(22, 18) = MVLEM_3DKlocal(18, 22);
	MVLEM_3DKlocal(22, 19) = MVLEM_3DKlocal(19, 22);
	MVLEM_3DKlocal(22, 20) = MVLEM_3DKlocal(20, 22);
	MVLEM_3DKlocal(22, 21) = MVLEM_3DKlocal(21, 22);
	MVLEM_3DKlocal(22, 22) = K18;
	MVLEM_3DKlocal(22, 23) = 0.0;

	MVLEM_3DKlocal(23, 0) = MVLEM_3DKlocal(0, 23);
	MVLEM_3DKlocal(23, 1) = MVLEM_3DKlocal(1, 23);
	MVLEM_3DKlocal(23, 2) = MVLEM_3DKlocal(2, 23);
	MVLEM_3DKlocal(23, 3) = MVLEM_3DKlocal(3, 23);
	MVLEM_3DKlocal(23, 4) = MVLEM_3DKlocal(4, 23);
	MVLEM_3DKlocal(23, 5) = MVLEM_3DKlocal(5, 23);
	MVLEM_3DKlocal(23, 6) = MVLEM_3DKlocal(6, 23);
	MVLEM_3DKlocal(23, 7) = MVLEM_3DKlocal(7, 23);
	MVLEM_3DKlocal(23, 8) = MVLEM_3DKlocal(8, 23);
	MVLEM_3DKlocal(23, 9) = MVLEM_3DKlocal(9, 23);
	MVLEM_3DKlocal(23, 10) = MVLEM_3DKlocal(10, 23);
	MVLEM_3DKlocal(23, 11) = MVLEM_3DKlocal(11, 23);
	MVLEM_3DKlocal(23, 12) = MVLEM_3DKlocal(12, 23);
	MVLEM_3DKlocal(23, 13) = MVLEM_3DKlocal(13, 23);
	MVLEM_3DKlocal(23, 14) = MVLEM_3DKlocal(14, 23);
	MVLEM_3DKlocal(23, 15) = MVLEM_3DKlocal(15, 23);
	MVLEM_3DKlocal(23, 16) = MVLEM_3DKlocal(16, 23);
	MVLEM_3DKlocal(23, 17) = MVLEM_3DKlocal(17, 23);
	MVLEM_3DKlocal(23, 18) = MVLEM_3DKlocal(18, 23);
	MVLEM_3DKlocal(23, 19) = MVLEM_3DKlocal(19, 23);
	MVLEM_3DKlocal(23, 20) = MVLEM_3DKlocal(20, 23);
	MVLEM_3DKlocal(23, 21) = MVLEM_3DKlocal(21, 23);
	MVLEM_3DKlocal(23, 22) = MVLEM_3DKlocal(22, 23);
	MVLEM_3DKlocal(23, 23) = (Km + Kh * (h * h) * ((c - 1.0) * (c - 1.0))) / ((2.0 * (d * d) + 2.0) * (2.0 * (d * d) + 2.0)) + (4.0 * Eib * Iib) / Lw;

	// Convert matrix from local to global cs
	MVLEM_3DK.addMatrixTripleProduct(0.0, T, MVLEM_3DKlocal, 1.0);

	// Return element stiffness matrix
	return MVLEM_3DK;
}

// Get element mass matrix assuming lumped mass
const Matrix & MVLEM_3D::getMass(void)
{

	// No rotational mass
	MVLEM_3DMlocal(0, 0) = NodeMass;
	MVLEM_3DMlocal(1, 1) = NodeMass;
	MVLEM_3DMlocal(2, 2) = NodeMass;

	MVLEM_3DMlocal(6, 6) = NodeMass;
	MVLEM_3DMlocal(7, 7) = NodeMass;
	MVLEM_3DMlocal(8, 8) = NodeMass;

	MVLEM_3DMlocal(12, 12) = NodeMass;
	MVLEM_3DMlocal(13, 13) = NodeMass;
	MVLEM_3DMlocal(14, 14) = NodeMass;

	MVLEM_3DMlocal(18, 18) = NodeMass;
	MVLEM_3DMlocal(19, 19) = NodeMass;
	MVLEM_3DMlocal(20, 20) = NodeMass;

	// Convert matrix from local to global cs
	MVLEM_3DM.addMatrixTripleProduct(0.0, T, MVLEM_3DMlocal, 1.0); 

	// Return element mass matrix
	return MVLEM_3DM;
}

// Get element damping matrix
const Matrix & MVLEM_3D::getDamp(void)
{
	MVLEM_3DD.Zero();

	MVLEM_3DD = this->Element::getDamp();

	// Return element damping matrix
	return MVLEM_3DD;
}

// zeroLoad
void MVLEM_3D::zeroLoad(void)
{
	// does nothing - no elemental loads
}

// addLoad
int MVLEM_3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	return 0;
}

int MVLEM_3D::addInertiaLoadToUnbalance(const Vector &accel)
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
			MVLEM_3DRlocal(6 * i + j) += - MVLEM_3DMlocal(6 * i + j, 6 * i + j) * RaccelL(6 * i + j);
		}
	}

	// Transform forces from local to global cs
	MVLEM_3DR.addMatrixTransposeVector(1.0, T, MVLEM_3DRlocal, 1.0);

	return 0;
}

// Get element force vector
const Vector & MVLEM_3D::getResistingForce()
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

	// Assemble force vector in local cs
	MVLEM_3DRlocal(0) = R1 / 2.0 + (Aib*Eib*dispL(0)) / Lw - (Aib*Eib*dispL(6)) / Lw;
	MVLEM_3DRlocal(1) = R2 / 2.0 - (R3*d) / (2.0 * (d*d) + 2.0) + (12.0 * Eib*Iib*dispL(1)) / (Lw*Lw*Lw) + (6.0 * Eib*Iib*dispL(5)) / (Lw*Lw) - (12.0 * Eib*Iib*dispL(7)) / (Lw*Lw*Lw) + (6.0 * Eib*Iib*dispL(11)) / (Lw*Lw);
	MVLEM_3DRlocal(2) = (Eave*(Tave*Tave*Tave)*(dispL(9))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(3))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(15))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(4))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(10))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(16))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(2))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(8))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(14))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(20))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(21))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(22))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(3) = (Eave*NUelastic*(Tave*Tave*Tave)*(dispL(4))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave)*(dispL(3))*((h*h) - (h*h)*NUelastic + 5 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(2))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(8))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1)) + (Eave*(Tave*Tave*Tave)*(dispL(14))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(21))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(9))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(15))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(20))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(4) = (Eave*NUelastic*(Tave*Tave*Tave)*(dispL(3))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (dispL(4))*((Eave*h*(Tave*Tave*Tave)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(10))*((Eave*h*(Tave*Tave*Tave)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(22))*((Eave*h*(Tave*Tave*Tave)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(16))*((Eave*h*(Tave*Tave*Tave)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) + (Eave*(Tave*Tave*Tave)*(dispL(2))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(8))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(14))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(20))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(5) = R3 / (2.0 * (d*d) + 2.0) + (6.0 * Eib*Iib*dispL(1)) / (Lw*Lw) + (4.0 * Eib*Iib*dispL(5)) / Lw - (6.0 * Eib*Iib*dispL(7)) / (Lw*Lw) + (2.0 * Eib*Iib*dispL(11)) / Lw;
	MVLEM_3DRlocal(6) = R1 / 2.0 - (Aib*Eib*dispL(0)) / Lw + (Aib*Eib*dispL(6)) / Lw;
	MVLEM_3DRlocal(7) = R2 / 2.0 + (R3*d) / (2.0 * (d*d) + 2.0) - (12.0 * Eib*Iib*dispL(1)) / (Lw*Lw*Lw) - (6.0 * Eib*Iib*dispL(5)) / (Lw*Lw) + (12.0 * Eib*Iib*dispL(7)) / (Lw*Lw*Lw) - (6.0 * Eib*Iib*dispL(11)) / (Lw*Lw);
	MVLEM_3DRlocal(8) = (Eave*(Tave*Tave*Tave)*(dispL(3))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(9))*(4.0 * (h*h)*NUelastic + (h*h) + 10 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(21))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(4))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(10))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(22))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(2))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(14))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(8))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(20))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(15))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(16))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(9) = (Eave*(Tave*Tave*Tave)*(dispL(2))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*NUelastic*(Tave*Tave*Tave)*(dispL(10))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave)*(dispL(9))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(15))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(8))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(20))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(3))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(14))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(21))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(10) = (Eave*(Tave*Tave*Tave)*(dispL(2))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (dispL(4))*((Eave*h*(Tave*Tave*Tave)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(16))*((Eave*h*(Tave*Tave*Tave)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*NUelastic*(Tave*Tave*Tave)*(dispL(9))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (dispL(22))*((Eave*h*(Tave*Tave*Tave)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*(Tave*Tave*Tave)*(dispL(10))*(5.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1)) - (Eave*(Tave*Tave*Tave)*(dispL(8))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(20))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(14))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(11) = R3 / (2.0 * (d*d) + 2.0) + (6.0 * Eib*Iib*dispL(1)) / (Lw*Lw) + (2.0 * Eib*Iib*dispL(5)) / Lw - (6.0 * Eib*Iib*dispL(7)) / (Lw*Lw) + (4.0 * Eib*Iib*dispL(11)) / Lw;
	MVLEM_3DRlocal(12) = R4 / 2.0 + (Aib*Eib*dispL(12)) / Lw - (Aib*Eib*dispL(18)) / Lw;
	MVLEM_3DRlocal(13) = R5 / 2.0 - (R6*d) / (2.0 * (d*d) + 2.0) + (12.0 * Eib*Iib*dispL(13)) / (Lw*Lw*Lw) + (6.0 * Eib*Iib*dispL(17)) / (Lw*Lw) - (12.0 * Eib*Iib*dispL(19)) / (Lw*Lw*Lw) + (6.0 * Eib*Iib*dispL(23)) / (Lw*Lw);
	MVLEM_3DRlocal(14) = (Eave*(Tave*Tave*Tave)*(dispL(3))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(15))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(21))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(4))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(16))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(22))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(2))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(8))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(14))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(20))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(9))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(10))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(15) = (Eave*(Tave*Tave*Tave)*(dispL(14))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(2))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(9))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(15))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*NUelastic*(Tave*Tave*Tave)*(dispL(16))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave)*(dispL(20))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(3))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(8))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(21))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(16) = (Eave*(Tave*Tave*Tave)*(dispL(14))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (dispL(16))*((Eave*h*(Tave*Tave*Tave)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(22))*((Eave*h*(Tave*Tave*Tave)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(10))*((Eave*h*(Tave*Tave*Tave)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*NUelastic*(Tave*Tave*Tave)*(dispL(15))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave)*(dispL(2))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (dispL(4))*((Eave*h*(Tave*Tave*Tave)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) - (Eave*(Tave*Tave*Tave)*(dispL(20))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(8))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(17) = R6 / (2.0 * (d*d) + 2.0) + (6.0 * Eib*Iib*dispL(13)) / (Lw*Lw) + (4.0 * Eib*Iib*dispL(17)) / Lw - (6.0 * Eib*Iib*dispL(19)) / (Lw*Lw) + (2.0 * Eib*Iib*dispL(23)) / Lw;
	MVLEM_3DRlocal(18) = R4 / 2.0 - (Aib*Eib*dispL(12)) / Lw + (Aib*Eib*dispL(18)) / Lw;
	MVLEM_3DRlocal(19) = R5 / 2.0 + (R6*d) / (2.0 * (d*d) + 2.0) - (12.0 * Eib*Iib*dispL(13)) / (Lw*Lw*Lw) - (6.0 * Eib*Iib*dispL(17)) / (Lw*Lw) + (12.0 * Eib*Iib*dispL(19)) / (Lw*Lw*Lw) - (6.0 * Eib*Iib*dispL(23)) / (Lw*Lw);
	MVLEM_3DRlocal(20) = (Eave*(Tave*Tave*Tave)*(dispL(9))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(15))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(21))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(10))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(16))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(22))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(2))*(5.0 * (h*h*h*h) + 5.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(8))*(5.0 * (h*h*h*h) - 10.0 * (Lw*Lw*Lw*Lw) - 7.0 * (h*h)*(Lw*Lw) + 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(14))*(10.0 * (h*h*h*h) - 5.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(20))*(10.0 * (h*h*h*h) + 10.0 * (Lw*Lw*Lw*Lw) + 7.0 * (h*h)*(Lw*Lw) - 2.0 * (h*h)*NUelastic*(Lw*Lw))) / (30.0 * (h*h*h)*(Lw*Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(3))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(4))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(21) = (Eave*NUelastic*(Tave*Tave*Tave)*(dispL(22))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (Eave*(Tave*Tave*Tave)*(dispL(3))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(8))*((h*h) - (h*h)*NUelastic + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(14))*(4.0 * (h*h)*NUelastic + (h*h) - 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(21))*((h*h) - (h*h)*NUelastic + 5.0 * (Lw*Lw))) / (45.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(20))*(4.0 * (h*h)*NUelastic + (h*h) + 10.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(2))*((h*h)*NUelastic - (h*h) + 5.0 * (Lw*Lw))) / (60.0 * (h*h)*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(9))*((h*h)*NUelastic - (h*h) + 10.0 * (Lw*Lw))) / (180.0 * h*Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(15))*(2.0 * (h*h)*NUelastic - 2.0 * (h*h) + 5.0 * (Lw*Lw))) / (90.0 * h*Lw*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(22) = (Eave*NUelastic*(Tave*Tave*Tave)*(dispL(21))) / (12.0 * (NUelastic*NUelastic) - 12.0) - (dispL(22))*((Eave*h*(Tave*Tave*Tave)) / (9.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (45.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(16))*((Eave*h*(Tave*Tave*Tave)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(4))*((Eave*h*(Tave*Tave*Tave)) / (36.0 * Lw*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*Lw*(NUelastic - 1.0)) / (180.0 * h*((NUelastic*NUelastic) - 1.0))) - (dispL(10))*((Eave*h*(Tave*Tave*Tave)) / (18.0 * Lw*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*Lw*(2.0 * NUelastic - 2.0)) / (90.0 * h*((NUelastic*NUelastic) - 1.0))) + (Eave*(Tave*Tave*Tave)*(dispL(8))*(4.0 * NUelastic*(Lw*Lw) - 5.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(14))*(10.0 * (h*h) - NUelastic*(Lw*Lw) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) - (Eave*(Tave*Tave*Tave)*(dispL(20))*(4.0 * NUelastic*(Lw*Lw) + 10.0 * (h*h) + (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0)) + (Eave*(Tave*Tave*Tave)*(dispL(2))*(NUelastic*(Lw*Lw) + 5.0 * (h*h) - (Lw*Lw))) / (60.0 * h*(Lw*Lw)*((NUelastic*NUelastic) - 1.0));
	MVLEM_3DRlocal(23) = R6 / (2.0 * (d*d) + 2.0) + (6.0 * Eib*Iib*dispL(13)) / (Lw*Lw) + (2.0 * Eib*Iib*dispL(17)) / Lw - (6.0 * Eib*Iib*dispL(19)) / (Lw*Lw) + (4.0 * Eib*Iib*dispL(23)) / Lw;

	// Convert force vector from local to global cs
	MVLEM_3DR.addMatrixTransposeVector(0.0, T, MVLEM_3DRlocal, 1.0);

	// Return element force vector
	return MVLEM_3DR;
}

// getResistingForceIncInertia
const Vector & MVLEM_3D::getResistingForceIncInertia()
{
	// if no mass terms .. just add damping terms
	if (density == 0.0) {
		// Compute the current resisting force
		this->getResistingForce();

		// Add the damping forces if rayleigh damping
		if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
			MVLEM_3DR += this->getRayleighDampingForces();	

		return MVLEM_3DR;
	}

	// Get nodal accelerations in global cs
	const Vector &accel1 = theNodes[0]->getTrialAccel();
	const Vector &accel2 = theNodes[1]->getTrialAccel();
	const Vector &accel3 = theNodes[2]->getTrialAccel();
	const Vector &accel4 = theNodes[3]->getTrialAccel();

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
			MVLEM_3DRlocal(6*i+j) += MVLEM_3DMlocal(6*i+j, 6*i+j) * accelL(6*i+j);
		}
	}

	// Transform forces from local to global cs
	MVLEM_3DR.addMatrixTransposeVector(1.0, T, MVLEM_3DRlocal, 1.0);

	// Add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
		MVLEM_3DR += this->getRayleighDampingForces();

	return MVLEM_3DR;
}

// sendSelf 
int MVLEM_3D::sendSelf(int commitTag, Channel &theChannel)
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

	// MVLEM_3D then sends the tags of it's four end nodes
	res = theChannel.sendID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING MVLEM_3D::sendSelf() - failed to send ID\n";
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
int MVLEM_3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res;
	int dataTag = this->getDbTag();

	// MVLEM_3D creates a Vector, receives the Vector and then sets the 
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
		opserr << "WARNING MVLEM_3D::recvSelf() - failed to receive Vector\n";
		return -1;
	}

	this->setTag((int)data(0));
	density = data(1);
	m = data(2);
	c = data(3);
	NUelastic = data(4);
	Tfactor = data(5);

	// MVLEM_3D now receives the tags of it's four external nodes
	res = theChannel.recvID(dataTag, commitTag, externalNodes);
	if (res < 0) {
		opserr << "WARNING MVLEM_3D::recvSelf() - failed to receive ID\n";
		return -2;
	}

	// Receive the material class tags
	ID matClassTags(2 * m + 1);
	res = theChannel.recvID(0, commitTag, matClassTags);

	// Allocate memory for the Concrete uniaxial materials
	theMaterialsConcrete = new UniaxialMaterial*[m];
	if (theMaterialsConcrete == 0) {
		opserr << "MVLEM_3D::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Allocate memory for the Steel uniaxial materials
	theMaterialsSteel = new UniaxialMaterial*[m];
	if (theMaterialsSteel == 0) {
		opserr << "MVLEM_3D::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Allocate memory for the Shear uniaxial material
	theMaterialsShear = new UniaxialMaterial*[1];
	if (theMaterialsShear == 0) {
		opserr << "MVLEM_3D::recvSelf() - "
			<< "failed to allocate pointers for uniaxial materials.\n";
		return -2;
	}

	// Receive the Concrete material models
	for (int i = 0; i < m; i++) {
		theMaterialsConcrete[i] = theBroker.getNewUniaxialMaterial(matClassTags(i));
		if (theMaterialsConcrete[i] == 0) {
			opserr << "MVLEM_3D::recvSelf() - "
				<< "failed to get blank uniaxial material.\n";
			return -3;
		}
		theMaterialsConcrete[i]->recvSelf(commitTag, theChannel, theBroker);
	}

	// Receive the Steel material models
	for (int i = 0; i < m; i++) {
		theMaterialsSteel[i] = theBroker.getNewUniaxialMaterial(matClassTags(i + m));
		if (theMaterialsSteel[i] == 0) {
			opserr << "MVLEM_3D::recvSelf() - "
				<< "failed to get blank uniaxial material.\n";
			return -3;
		}
		theMaterialsSteel[i]->recvSelf(commitTag, theChannel, theBroker);
	}

	// Receive the Shear material model
	theMaterialsShear[0] = theBroker.getNewUniaxialMaterial(matClassTags(2 * m));
	if (theMaterialsShear[0] == 0) {
		opserr << "MVLEM_3D::recvSelf() - "
			<< "failed to get blank uniaxial material.\n";
		return -3;
	}
	theMaterialsShear[0]->recvSelf(commitTag, theChannel, theBroker);

	return 0;
}

// Display model
int MVLEM_3D::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
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

		// add displacement (multiplied with the displacement factor) to the original node location to obtain current node location
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
		// Local node 4 - top right
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


void MVLEM_3D::Print(OPS_Stream &s, int flag)
{
	if (flag == 0)
	{
		// Print out element properties
		s << "Element: " << this->getTag() << endln;
		s << "  type: MVLEM_3D" << endln;
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
Response *MVLEM_3D::setResponse(const char **argv, int argc, OPS_Stream &s)
{
	Response *theResponse = 0;

	s.tag("ElementOutput");
	s.attr("eleType", "MVLEM_3D");
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

		theResponse = new ElementResponse(this, 1, Vector(24));

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

		theResponse = new ElementResponse(this, 2, Vector(24));

	}

	// Element curvature
	else if (strcmp(argv[0], "Curvature") == 0 || strcmp(argv[0], "curvature") == 0) {

		s.tag("ResponseType", "fi");

		theResponse = new ElementResponse(this, 3, 0.0);
	}

	// Fiber strain
	else if (strcmp(argv[0], "Fiber_Strain") == 0 || strcmp(argv[0], "fiber_strain") == 0) {

		for (int pointNum = 1; pointNum <= m; ++pointNum) {
			s.tag("GaussPointOutput");
			s.attr("number", pointNum);
			s.attr("eta", x[pointNum - 1] / Lw * 2.0);
			s.attr("weight", b[pointNum - 1] / Lw * 2.0);
			s.tag("ResponseType", "epsy");
			s.endTag();
		}

		theResponse = new ElementResponse(this, 4, Vector(m));
	}

	// Fiber concrete stresses
	else if (strcmp(argv[0], "Fiber_Stress_Concrete") == 0 || strcmp(argv[0], "fiber_stress_concrete") == 0) {

		for (int pointNum = 1; pointNum <= m; ++pointNum) {
			s.tag("GaussPointOutput");
			s.attr("number", pointNum);
			s.attr("eta", x[pointNum - 1] / Lw * 2.0);
			s.attr("weight", b[pointNum - 1] / Lw * 2.0);
			s.tag("ResponseType", "sigmayc");
			s.endTag();
		}

		theResponse = new ElementResponse(this, 5, Vector(m));
	}

	// Fiber steel stresses
	else if (strcmp(argv[0], "Fiber_Stress_Steel") == 0 || strcmp(argv[0], "fiber_stress_steel") == 0) {

		for (int pointNum = 1; pointNum <= m; ++pointNum) {
			s.tag("GaussPointOutput");
			s.attr("number", pointNum);
			s.attr("eta", x[pointNum - 1] / Lw * 2.0);
			s.attr("weight", b[pointNum - 1] / Lw * 2.0);
			s.tag("ResponseType", "sigmays");
			s.endTag();
		}

		theResponse = new ElementResponse(this, 6, Vector(m));
	}

	// Shear force deformation
	else if (strcmp(argv[0], "Shear_Force_Deformation") == 0 || strcmp(argv[0], "shear_force_deformation") == 0) {

		s.tag("ResponseType", "shearDef");
		s.tag("ResponseType", "shearFrc");

		theResponse = new ElementResponse(this, 7, Vector(2));
	}

	// Shear Deformation
	else if (strcmp(argv[0], "ShearDef") == 0 || strcmp(argv[0], "sheardef") == 0) {

		s.tag("ResponseType", "shearDef");

		theResponse = new ElementResponse(this, 8, 0.0);
	}

	// Results at specified integration point
	else if (strcmp(argv[0], "material") == 0 && (argc > 2)) {
		int pointNum = atoi(argv[1]);
		if (pointNum > 0 && pointNum <= m) {
			s.tag("GaussPointOutput");
			s.attr("number", pointNum);
			s.attr("eta", x[pointNum - 1] / Lw * 2.0);
			s.attr("weight", b[pointNum - 1] / Lw * 2.0);
			if (argc > 3) {
				// specify material and forward result type
				if (strcmp(argv[2], "concrete") == 0 || strcmp(argv[2], "Concrete") == 0) {
					theResponse = theMaterialsConcrete[pointNum - 1]->setResponse(&argv[3], argc - 3, s);
				}
				else if (strcmp(argv[2], "steel") == 0 || strcmp(argv[2], "Steel") == 0) {
					theResponse = theMaterialsSteel[pointNum - 1]->setResponse(&argv[3], argc - 3, s);
				}
				else if (strcmp(argv[2], "shear") == 0 || strcmp(argv[2], "Shear") == 0) {
					theResponse = theMaterialsShear[0]->setResponse(&argv[3], argc - 3, s);
				}
			}
			s.endTag();
		}
	}

	s.endTag();

	return theResponse;
}

// get recorders
int MVLEM_3D::getResponse(int responseID, Information &eleInfo)
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
		return eleInfo.setDouble(this->getShearDef());

	default:

		return 0;

	}
}

// Return element local forces
Vector MVLEM_3D::getResistingForceLocal(void)
{
	return MVLEM_3DRlocal;
}

// Get curvature (from vertical strains)
double MVLEM_3D::getCurvature(void)
{
	double Curv;

	Curv = (MVLEM_3DStrain[0] - MVLEM_3DStrain[m - 1]) / (x[0] - x[m - 1]);

	return Curv;
}

// Get fiber strains
Vector MVLEM_3D::getStrain(void)
{
	Vector fiberStrain(m);

	for (int i = 0; i<m; i++) {
		fiberStrain(i) = MVLEM_3DStrain[i];
	}

	return fiberStrain;
}

// Get Concrete Stress 
Vector MVLEM_3D::getStressConcrete(void)
{
	Vector concreteStress(m);

	for (int i = 0; i<m; i++) {
		concreteStress(i) = theMaterialsConcrete[i]->getStress();
	}

	return concreteStress;
}

// Get Steel Stress 
Vector MVLEM_3D::getStressSteel(void)
{
	Vector steelStress(m);

	for (int i = 0; i<m; i++) {
		steelStress(i) = theMaterialsSteel[i]->getStress();
	}

	return steelStress;
}

// Get Shear Stress-Strain 
Vector MVLEM_3D::getShearFD(void)
{
	Vector shearStrainStress(2);

	shearStrainStress(0) = theMaterialsShear[0]->getStrain();
	shearStrainStress(1) = theMaterialsShear[0]->getStress();

	return shearStrainStress;
}

// Get Shear Stress-Strain 
double MVLEM_3D::getShearDef(void)
{
	double shearDef;

	shearDef = theMaterialsShear[0]->getStrain();

	return shearDef;
}

// Compute element transformation matrix
void  MVLEM_3D::setTransformationMatrix(void) {

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

}
