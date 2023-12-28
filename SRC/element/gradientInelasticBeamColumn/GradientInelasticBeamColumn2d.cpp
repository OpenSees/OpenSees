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

/* Written by: Mohammad Salehi (mohammad.salehi@tamu.edu)
** Created: 07/19
** Description: The source code for the 2D gradient inelastic (GI) force-based beam-column element formulation
**
**
** References:
**
** Mohammad Salehi and Petros Sideris (2017)
** “Refined Gradient Inelastic Flexibility-Based Formulation for Members Subjected to Arbitrary Loading”
** ASCE Journal of Engineering Mechanics, 143(9): 04017090
**
** Petros Sideris and Mohammad Salehi (2016)
** “A Gradient Inelastic Flexibility-Based Frame Element Formulation”
** ASCE Journal of Engineering Mechanics, 142(7): 04016039
*/

#include "GradientInelasticBeamColumn2d.h"

#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>

#include <NewtonCotesBeamIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <SimpsonBeamIntegration.h>
#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

// Method to Read Command Arguments
void* OPS_GradientInelasticBeamColumn2d()
{
	// Necessary Arguments
	if (OPS_GetNumRemainingInputArgs() < 6) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - insufficient arguments\n" <<
			"         Want: eleTag? iNode? jNode? transfTag? integrationTag? lc?\n" <<
			"         <-constH> <-iter maxIter? minTol? maxTol?> <-corControl maxEpsInc? maxPhiInc?>\n";
		return 0;
	}

	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();
	if (ndm != 2 || ndf != 3) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - ndm must be 2 and ndf must be 3\n";
		return 0;
	}

	// inputs: 
	int iData[5];
	int numData = 5;
	if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - invalid input tags\n";
		return 0;
	}

	int eleTag = iData[0];
	int nodeTagI = iData[1];
	int nodeTagJ = iData[2];
	int transfTag = iData[3];
	int integrTag = iData[4];

	double lc;
	numData = 1;
	if (OPS_GetDoubleInput(&numData, &lc) < 0) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - invalid double input\n";
		return 0;
	}
	
	double lam1 = 0.1, lam2 = 0.1;	// would not affect the results
	
	// Optional Arguments
	int maxIter = 50;
	double minTol = 1E-10, maxTol = 1E-8;
	bool correctionControl = false;
	bool constH = false;
	double maxEpsInc = 0.0, maxPhiInc = 0.0;

	numData = 1;
	while (OPS_GetNumRemainingInputArgs() > 0) {
		const char* word = OPS_GetString();

		if (strcmp(word, "-constH") == 0)
			constH = true;
		else if (strcmp(word, "-iter") == 0) {
			if (OPS_GetNumRemainingInputArgs() > 2) {
				if (OPS_GetIntInput(&numData, &maxIter) < 0) {
					opserr << "WARNING! gradientInelasticBeamColumn2d - invalid maxIter\n";
					return 0;
				}
				if (OPS_GetDoubleInput(&numData, &minTol) < 0) {
					opserr << "WARNING! gradientInelasticBeamColumn2d - invalid minTol\n";
					return 0;
				}
				if (OPS_GetDoubleInput(&numData, &maxTol) < 0) {
					opserr << "WARNING! gradientInelasticBeamColumn2d - invalid maxTol\n";
					return 0;
				}
			}
			else {
				opserr << "WARNING! gradientInelasticBeamColumn2d - need maxIter? minTol? maxTol? after -iter \n";
				return 0;
			}
		}
		else if (strcmp(word, "-corControl") == 0) {
			correctionControl = true;

			if (OPS_GetNumRemainingInputArgs() > 1) {
				if (OPS_GetDoubleInput(&numData, &maxEpsInc) < 0) {
					opserr << "WARNING! gradientInelasticBeamColumn2d - invalid maxEpsInc\n";
					return 0;
				}
				if (OPS_GetDoubleInput(&numData, &maxPhiInc) < 0) {
					opserr << "WARNING! gradientInelasticBeamColumn2d - invalid maxPhiInc\n";
					return 0;
				}
			}
			else
				opserr << "WARNING! gradientInelasticBeamColumn2d - no max. correction increments set\n" <<
				"         -> setting them automatically|\n";
		}
	}

	// check transf
	CrdTransf* theTransf = OPS_getCrdTransf(transfTag);
	if (theTransf == 0) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - CrdTransf with tag " << transfTag << " not found\n";
		return 0;
	}

	// check beam integrataion
	BeamIntegrationRule* theRule = OPS_getBeamIntegrationRule(integrTag);
	if (theRule == 0) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - BeamIntegrationRule with tag " << integrTag << " not found\n";
		return 0;
	}

	BeamIntegration* beamIntegr = theRule->getBeamIntegration();
	if (beamIntegr == 0) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - failed to create beam integration\n";
		return 0;
	}

	// check sections
	const ID& secTags = theRule->getSectionTags();
	int numIntegrPoints = secTags.Size();

	for (int i = 2; i < numIntegrPoints - 1; i++) {
		if (secTags(i) != secTags(i - 1)) {
			opserr << "WARNING! gradientInelasticBeamColumn2d - internal integration points should have identical tags\n"
				<< "continued using section tag of integration point 2 for all internal integration points\n";
			return 0;
		}
	}

	SectionForceDeformation* endSection1 = OPS_getSectionForceDeformation(secTags(0));
	if (!endSection1) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - section with tag " << secTags(0) << " not found\n";
		return 0;
	}

	SectionForceDeformation* intSection = OPS_getSectionForceDeformation(secTags(1));
	if (!intSection) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - section with tag " << secTags(1) << " not found\n";
		return 0;
	}

	SectionForceDeformation* endSection2 = OPS_getSectionForceDeformation(secTags(numIntegrPoints - 1));
	if (!endSection2) {
		opserr << "WARNING! gradientInelasticBeamColumn2d - section with tag " << secTags(numIntegrPoints - 1) << " not found\n";
		return 0;
	}

	Element* theEle = new GradientInelasticBeamColumn2d(eleTag, nodeTagI, nodeTagJ, numIntegrPoints, *endSection1, *intSection, *endSection2,
		lam1, lam2, *beamIntegr, *theTransf, lc, minTol, maxTol, maxIter, constH, correctionControl, maxEpsInc, maxPhiInc);

	return theEle;
}

// Initialize Class Wide Variables
Matrix GradientInelasticBeamColumn2d::theMatrix(6, 6);
Vector GradientInelasticBeamColumn2d::theVector(6);

// Constructor 1 (for normal processing)
GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d(int tag, int nodeI, int nodeJ,
	int numSec, SectionForceDeformation &endSec1, SectionForceDeformation &sec, SectionForceDeformation &endSec2, double R1, double R2,
	BeamIntegration &BI, CrdTransf &CT, double LC,
	double minTolerance, double maxTolerance, int maxNumIters,
	bool constH,
	bool corControl, double maxEps, double maxPhi)
	: Element(tag, ELE_TAG_GradientInelasticBeamColumn2d), connectedExternalNodes(2),
	numSections(numSec), sections(0),
	beamIntegr(0), crdTransf(0), lc(LC),
	maxIters(maxNumIters), minTol(minTolerance), maxTol(maxTolerance), 
	cnstH(constH), secLR1(R1), secLR2(R2),
	correctionControl(corControl), maxEpsInc(maxEps), maxPhiInc(maxPhi),
	L(0.0), F_tol_q(0.0), F_tol_f_ms(0.0), 
	hh(0), H(0), H_init(0), H_inv(0),
	B_q(0), B_Q(0), B_q_H_inv_init(0), K0(0),
	J(0), J_init(0), J_commit(0),
	k_init(3), flex_ms_init(0), trial_change(0), max_trial_change(0),
	initialFlag(0), Q(3), Q_commit(3),
	d_sec(0), d_sec_commit(0), d_tot(0), d_tot_commit(0), d_nl_tot(0), d_nl_tot_commit(0),
	F_ms(0), F_ms_commit(0),
	iterNo(0), strIterNo(0), totStrIterNo(0), commitNo(0), iters(3)
	// complete
{
	// Pointers to Nodes and Their IDs
	if (connectedExternalNodes.Size() != 2) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d() - element: " << this->getTag() << " - failed to create an ID of size 2\n";
		exit(-1);
	}

	connectedExternalNodes(0) = nodeI;
	connectedExternalNodes(1) = nodeJ;

	theNodes[0] = 0;
	theNodes[1] = 0;

	// Get Copy of Integration Method
	beamIntegr = BI.getCopy();

	if (!beamIntegr) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d() - element: " << this->getTag() << " - could not create copy of beam integration object" << endln;
		exit(-1);
	}

	sections = new SectionForceDeformation *[numSections];
	if (!sections) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d() - element: " << this->getTag() << " - could not allocate section pointers\n";
		exit(-1);
	}

	double *secX = new double[numSections];
	beamIntegr->getSectionLocations(numSections, L, secX);	// relative locations of sections (x/L)

	for (int i = 0; i < numSections; i++) {
		if (secX[i] >= 1.0 - secLR2)
			sections[i] = endSec2.getCopy();
		else if (secX[i] > secLR1)
			sections[i] = sec.getCopy();
		else
			sections[i] = endSec1.getCopy();

		if (!sections[i]) {
			opserr << "WARNING! GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d() - element: " << this->getTag() << " - could not create copy of section " << i + 1 << endln;
			exit(-1);
		}
	}

	if (secX)
		delete[] secX;

	// Check Sections Order
	secOrder = sec.getOrder();

	if (secOrder < 2) {
		opserr << "ERROR! GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d() - element: " << this->getTag() << " - section order must be at least 2" << endln;
		exit(-1);
	}

	// Initialize Matrices
	B_q = new Matrix(3, numSections * secOrder);

	B_Q = new Matrix(numSections * secOrder, 3);

	H = new Matrix(numSections * secOrder, numSections * secOrder);
	H_init = new Matrix(numSections * secOrder, numSections * secOrder);
	H_inv = new Matrix(numSections * secOrder, numSections * secOrder);
	hh = new Vector(numSections * secOrder);

	B_q_H_inv_init = new Matrix(3, numSections * secOrder);

	J = new Matrix(3 + numSections * secOrder, 3 + numSections * secOrder);
	J_init = new Matrix(3 + numSections * secOrder, 3 + numSections * secOrder);
	J_commit = new Matrix(3 + numSections * secOrder, 3 + numSections * secOrder);

	flex_ms_init = new Vector(numSections * secOrder);
	trial_change = new Vector(numSections * secOrder + 3);
	max_trial_change = new Vector(numSections * secOrder + 3);

	d_tot = new Vector(numSections * secOrder);
	d_tot_commit = new Vector(numSections * secOrder);
	d_nl_tot = new Vector(numSections * secOrder);
	d_nl_tot_commit = new Vector(numSections * secOrder);

	F_ms = new Vector(numSections * secOrder);
	F_ms_commit = new Vector(numSections * secOrder);

	// Get Copy of Coordinate Transformation Method 
	crdTransf = CT.getCopy2d();

	if (!crdTransf) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d() - element: " << this->getTag() << " - could not create copy of coordinate transformation object " << endln;
		exit(-1);
	}

	// Allocate Array of Pointers to Analysis State Variables
	d_sec = new Vector[numSections];
	if (!d_sec) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d() - element: " << this->getTag() << " - could not allocate section deformation pointers\n";
		exit(-1);
	}

	d_sec_commit = new Vector[numSections];
	if (!d_sec_commit) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d() - element: " << this->getTag() << " - could not allocate committed section deformation pointers\n";
		exit(-1);
	}
}

void
GradientInelasticBeamColumn2d::setSectionPointers(void)
{
  if (numSections < 1)
    return;

  // Initialize Matrices
  B_q = new Matrix(3, numSections * secOrder);
  
  B_Q = new Matrix(numSections * secOrder, 3);
  
  H = new Matrix(numSections * secOrder, numSections * secOrder);
  H_init = new Matrix(numSections * secOrder, numSections * secOrder);
  H_inv = new Matrix(numSections * secOrder, numSections * secOrder);
  hh = new Vector(numSections * secOrder);
  
  B_q_H_inv_init = new Matrix(3, numSections * secOrder);
  
  J = new Matrix(3 + numSections * secOrder, 3 + numSections * secOrder);
  J_init = new Matrix(3 + numSections * secOrder, 3 + numSections * secOrder);
  J_commit = new Matrix(3 + numSections * secOrder, 3 + numSections * secOrder);
  
  flex_ms_init = new Vector(numSections * secOrder);
  trial_change = new Vector(numSections * secOrder + 3);
  max_trial_change = new Vector(numSections * secOrder + 3);
  
  d_tot = new Vector(numSections * secOrder);
  d_tot_commit = new Vector(numSections * secOrder);
  d_nl_tot = new Vector(numSections * secOrder);
  d_nl_tot_commit = new Vector(numSections * secOrder);
  
  F_ms = new Vector(numSections * secOrder);
  F_ms_commit = new Vector(numSections * secOrder);  
}

// Constructor 2 (for parallel processing)
GradientInelasticBeamColumn2d::GradientInelasticBeamColumn2d()
  : Element(0, ELE_TAG_GradientInelasticBeamColumn2d), connectedExternalNodes(2),
    numSections(0), sections(0),
    beamIntegr(0), crdTransf(0), lc(0.0),
    maxIters(0), minTol(0.0), maxTol(0.0),
    cnstH(0), secLR1(0.0), secLR2(0.0),
    correctionControl(0), maxEpsInc(0.0), maxPhiInc(0.0),    
    L(0.0), F_tol_q(0.0), F_tol_f_ms(0.0),
    hh(0), H(0), H_init(0), H_inv(0),
    B_q(0), B_Q(0), B_q_H_inv_init(0), K0(0),
    J(0), J_init(0), J_commit(0),
    k_init(3), flex_ms_init(0), trial_change(0), max_trial_change(0),
    initialFlag(0), Q(3), Q_commit(3),
    d_sec(0), d_sec_commit(0), d_tot(0), d_tot_commit(0), d_nl_tot(0), d_nl_tot_commit(0),
    F_ms(0), F_ms_commit(0),
    iterNo(0), strIterNo(0), totStrIterNo(0), commitNo(0), iters(3)
    // complete
{
	// Set Node Pointers to 0
	theNodes[0] = 0;
	theNodes[1] = 0;
}

// Destructor
GradientInelasticBeamColumn2d::~GradientInelasticBeamColumn2d()
{
	// Delete Matrix Pointers
	if (B_q != 0)
		delete B_q;

	if (B_Q != 0)
		delete B_Q;

	if (H != 0)
		delete H;

	if (H_init != 0)
		delete H_init;

	if (H_inv != 0)
		delete H_inv;

	if (hh != 0)
		delete hh;

	if (B_q_H_inv_init != 0)
		delete B_q_H_inv_init;

	if (J != 0)
		delete J;

	if (J_init != 0)
		delete J_init;

	if (J_commit != 0)
		delete J_commit;

	if (flex_ms_init != 0)
		delete flex_ms_init;

	if (trial_change != 0)
		delete trial_change;

	if (max_trial_change != 0)
		delete max_trial_change;

	if (d_tot != 0)
		delete d_tot;

	if (d_tot_commit != 0)
		delete d_tot_commit;

	if (d_nl_tot != 0)
		delete d_nl_tot;

	if (d_nl_tot_commit != 0)
		delete d_nl_tot_commit;

	if (F_ms != 0)
		delete F_ms;

	if (F_ms_commit != 0)
		delete F_ms_commit;

	if (K0 != 0)
		delete K0;

	// Delete Section Pointers
	if (sections) {
		for (int i = 0; i < numSections; i++)
			if (sections[i])
				delete sections[i];
		delete[] sections;
	}

	// Delete Beam Integration Pointer and Coordinate Transformation Pointer
	if (beamIntegr != 0)
		delete beamIntegr;

	if (crdTransf != 0)
		delete crdTransf;

	// Delete Analysis Arrays
	if (d_sec != 0)
		delete[] d_sec;

	if (d_sec_commit != 0)
		delete[] d_sec_commit;
}

// Definition of setDomain()
void
GradientInelasticBeamColumn2d::setDomain(Domain *theDomain)
{
	// Check Domain is not Null
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;

		opserr << "ERROR! GradientInelasticBeamColumn2d::setDomain() - element: " << this->getTag() << " - the domain is null\n";
		exit(0);
	}

	// Get Node Pointers
	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);

	if (!theNodes[0]) {
		opserr << "ERROR! GradientInelasticBeamColumn2d::setDomain() - element: " << this->getTag() << " - node " << Nd1 << " does not exist in the domain\n";
		exit(0);
	}

	if (!theNodes[1]) {
		opserr << "ERROR! GradientInelasticBeamColumn2d::setDomain() - element: " << this->getTag() << " - node " << Nd2 << " does not exist in the domain\n";
		exit(0);
	}

	// Call DomainComponent Class Method
	this->DomainComponent::setDomain(theDomain);

	// Check Node DOFs are 3
	int dofNd1 = theNodes[0]->getNumberDOF();
	int dofNd2 = theNodes[1]->getNumberDOF();

	if (dofNd1 != 3) {
		opserr << "ERROR! GradientInelasticBeamColumn2d::setDomain() - element: " << this->getTag() << " - node " << Nd1 << " has incorrect number of DOFs (not 3)\n";
		exit(0);
	}

	if (dofNd2 != 3) {
		opserr << "ERROR! GradientInelasticBeamColumn2d::setDomain() - element: " << this->getTag() << " - node " << Nd2 << " has incorrect number of DOFs (not 3)\n";
		exit(0);
	}

	// Initialize Coordinate Transformation
	if (crdTransf->initialize(theNodes[0], theNodes[1])) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::setDomain() - element: " << this->getTag() << " - coordinate transformation object could not be initialized\n";
		exit(0);
	}

	// Determine Element Initial Length
	L = crdTransf->getInitialLength();

	if (L == 0.0) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::setDomain() - element: " << this->getTag() << " - element length is zero\n";
		exit(0);
	}

	// Form Total Force Interpolation and Total Displacement Integration Matrices
	double *secX = new double[numSections];
	beamIntegr->getSectionLocations(numSections, L, secX);	// relative locations of sections (x/L)

	double *secW = new double[numSections];
	beamIntegr->getSectionWeights(numSections, L, secW);	// relative weights of sections (w/L)

	B_Q->Zero();
	B_q->Zero();

	double w, x;
	Vector dx(numSections - 1);	// spaces between integration points

	for (int j = 0; j < numSections; j++) {
		const ID &code = sections[j]->getType();

		w = L * secW[j];
		x = L * secX[j];

		if (j < numSections - 1)
			dx(j) = L * (secX[j + 1] - secX[j]);

		for (int i = 0; i < secOrder; i++) {
			switch (code(i)) {
			case SECTION_RESPONSE_P:
				(*B_Q)(j * secOrder + i, 0) = 1.0;
				(*B_q)(0, j * secOrder + i) = w;
				break;
			case SECTION_RESPONSE_MZ:
				(*B_Q)(j * secOrder + i, 1) = x / L - 1.0;
				(*B_Q)(j * secOrder + i, 2) = x / L;
				(*B_q)(1, j * secOrder + i) = w * (x / L - 1.0);
				(*B_q)(2, j * secOrder + i) = w * x / L;
				break;
			case SECTION_RESPONSE_VY:
				(*B_Q)(j * secOrder + i, 1) = (*B_Q)(j * secOrder + i, 2) = -1.0 / L;
				(*B_q)(1, j * secOrder + i) = (*B_q)(2, j * secOrder + i) = -w / L;
				break;
			default:
				break;
			}
		}
	}

	if (secX)
		delete[] secX;

	// Form H Matrix
	H_init->Zero();

	/// 2nd Order PDE, Dirichlet BCs
	for (int i = 0; i < secOrder; i++) {
		(*H_init)(i, i) = 1.0;
		(*H_init)(secOrder * numSections - i - 1, secOrder * numSections - i - 1) = 1.0;
	}

	for (int j = 1; j < numSections - 1; j++) {
		for (int i = 0; i < secOrder; i++) {
			(*H_init)(j * secOrder + i, (j - 1) * secOrder + i) = -lc * lc / (dx(j - 1) * (dx(j - 1) + dx(j)));
			(*H_init)(j * secOrder + i, j * secOrder + i) = 1 + lc * lc / (dx(j - 1) * dx(j));
			(*H_init)(j * secOrder + i, (j + 1) * secOrder + i) = -lc * lc / (dx(j) * (dx(j - 1) + dx(j)));
		}
	}

	/// 4th Order PDE, Neumann BCs with Zero 1st Derivative of Nonlocal Strains at the Ends
	/*double Ac = 1.0 + pow(lc / dx(1), 2.0) + 0.75 * pow(lc / dx(1), 4.0);
	double Bc = -0.5 * (pow(lc / dx(1), 2.0) + pow(lc / dx(1), 4.0));
	double Cc = 0.125 * pow(lc / dx(1), 4.0);

	for (int i = 0; i < secOrder; i++) {
		(*H_init)(i, i) = 1.0;
		(*H_init)(secOrder * numSections - i - 1, secOrder * numSections - i - 1) = 1.0;

		(*H_init)(i + secOrder, i + secOrder) = Ac + Cc;
		(*H_init)(i + secOrder, i) = (*H_init)(i + secOrder, i + 2 * secOrder) = Bc;
		(*H_init)(i + secOrder, i + 3 * secOrder) = Cc;

		(*H_init)(secOrder * numSections - secOrder - i - 1, secOrder * numSections - secOrder - i - 1) = Ac + Cc;
		(*H_init)(secOrder * numSections - secOrder - i - 1, secOrder * numSections - i - 1) = (*H_init)(secOrder * numSections - secOrder - i - 1, secOrder * numSections - 2 * secOrder - i - 1) = Bc;
		(*H_init)(secOrder * numSections - secOrder - i - 1, secOrder * numSections - 3 * secOrder - i - 1) = Cc;
	}

	for (int i = 2 * secOrder; i < ((numSections - 2) * secOrder); i++) {
		(*H_init)(i, i - 2 * secOrder) = (*H_init)(i, i + 2 * secOrder) = Cc;
		(*H_init)(i, i - secOrder) = (*H_init)(i, i + secOrder) = Bc;
		(*H_init)(i, i) = Ac;
	}*/

	/// 4th Order PDE, Neumann BCs with Zero 2nd Derivative of Nonlocal Strains at the Ends
	/*double Ac = 1.0 + pow(lc / dx(1), 2.0) + 0.75 * pow(lc / dx(1), 4.0);
	double Bc = -0.5 * (pow(lc / dx(1), 2.0) + pow(lc / dx(1), 4.0));
	double Cc = 0.125 * pow(lc / dx(1), 4.0);

	for (int i = 0; i < secOrder; i++) {
		(*H_init)(i, i) = 1.0;
		(*H_init)(secOrder * numSections - i - 1, secOrder * numSections - i - 1) = 1.0;

		(*H_init)(i + secOrder, i) = Bc + 2.0 * Cc;
		(*H_init)(i + secOrder, i + secOrder) = Ac - Cc;
		(*H_init)(i + secOrder, i + 2 * secOrder) = Bc;
		(*H_init)(i + secOrder, i + 3 * secOrder) = Cc;

		(*H_init)(secOrder * numSections - secOrder - i - 1, secOrder * numSections - i - 1) = Bc + 2.0 * Cc;
		(*H_init)(secOrder * numSections - secOrder - i - 1, secOrder * numSections - secOrder - i - 1) = Ac - Cc;
		(*H_init)(secOrder * numSections - secOrder - i - 1, secOrder * numSections - 2 * secOrder - i - 1) = Bc;
		(*H_init)(secOrder * numSections - secOrder - i - 1, secOrder * numSections - 3 * secOrder - i - 1) = Cc;
	}

	for (int i = 2 * secOrder; i < ((numSections - 2) * secOrder); i++) {
		(*H_init)(i, i - 2 * secOrder) = (*H_init)(i, i + 2 * secOrder) = Cc;
		(*H_init)(i, i - secOrder) = (*H_init)(i, i + secOrder) = Bc;
		(*H_init)(i, i) = Ac;
	}*/

	*H = *H_init;

	// Initialize Section Deformations
	if (!initialFlag) {
		for (int i = 0; i < numSections; i++) {
			d_sec[i] = Vector(secOrder);
			d_sec_commit[i] = Vector(secOrder);

			d_sec[i].Zero();
			d_sec_commit[i].Zero();
		}

		d_tot->Zero();
		d_tot_commit->Zero();

		d_nl_tot->Zero();
		d_nl_tot_commit->Zero();

		F_ms->Zero();
		F_ms_commit->Zero();
	}
	
	// Compute H_inv and B_q_H_inv
	if (H_init->Invert(*H_inv) < 0) {
	  opserr << "WARNING! GradientInelasticBeamColumn2d::setDomain() - element: " << this->getTag() << " - could not invert H matrix\n";
	  exit(0);
	}
	
	*B_q_H_inv_init = (*B_q) * (*H_inv);
	//}
		
	// Form Initial Jacobian Matrix
	Matrix K_ms(numSections * secOrder, numSections * secOrder);
	this->getSectionsInitialStiff(K_ms);

	J_init->Zero();
	this->assembleMatrix(*J_init, *B_Q, 0, numSections * secOrder - 1, 0, 2, 1.0);
	this->assembleMatrix(*J_init, K_ms, 0, numSections * secOrder - 1, 3, numSections * secOrder + 2, -1.0);
	this->assembleMatrix(*J_init, *B_q_H_inv_init, numSections * secOrder, numSections * secOrder + 2, 3, numSections * secOrder + 2, -1.0);

	*J = *J_init;
	*J_commit = *J;

	// Determine Element's Initial Stiffness Matrix to Use as Weight Matrix for Norm Calculations
	Matrix k_el_init = this->getInitialBasicStiff();

	for (int i = 0; i < 3; i++)
		k_init(i) = k_el_init(i, i);

	// Determine Element's Initial Sections' Felxibility Matrix
	for (int i = 0; i < numSections; i++) {
		double W = secW[i] * L;

		for (int j = 0; j < secOrder; j++) {
			int k = i * secOrder + j;
			(*flex_ms_init)(k) = W / K_ms(k, k);
		}
	}

	if (secW)
		delete[] secW;

	// Determine Maximum Corrections (for d_tot and Q in iterations)
	if (correctionControl) {
		// Zero trial_change
		trial_change->Zero();

		// Decide on Maximum Trial Change
		if (maxEpsInc != 0.0) {
			const ID &code = sections[0]->getType();

			for (int i = 0; i < secOrder; i++) {
				for (int j = 0; j < numSections; j++) {
					int k = j * secOrder + i;

					switch (code(i)) {
					case SECTION_RESPONSE_P:
						(*max_trial_change)(3 + k) = maxEpsInc;
						break;
					case SECTION_RESPONSE_MZ:
						(*max_trial_change)(3 + k) = maxPhiInc;
						break;
					case SECTION_RESPONSE_VY:
						(*max_trial_change)(3 + k) = maxEpsInc;
						break;
					default:
						break;
					}
				}
			}

			Vector d_tot_max(numSections * secOrder);
			d_tot_max.Extract(*max_trial_change, 3);
			Vector Q_max = k_el_init * (*B_q) * d_tot_max;

			for (int i = 0; i < 3; i++)
				(*max_trial_change)(i) = fabs(Q_max(i));
		}
		else
			max_trial_change->Zero();
	}

	// Determine factors for tol_q and tol_f_ms
	F_tol_q = sqrt(k_init(0) + k_init(1) + k_init(2));	// for nodal displacements

	F_tol_f_ms = 0.0;
	for (int i = 0; i < numSections * secOrder; i++)
		F_tol_f_ms += (*flex_ms_init)(i);

	F_tol_f_ms = sqrt(F_tol_f_ms);	// for section forces
}

// Definition of Methods Dealing with Nodes Information
int
GradientInelasticBeamColumn2d::getNumExternalNodes(void) const
{
	return 2;
}

const ID &
GradientInelasticBeamColumn2d::getExternalNodes(void)
{
	return connectedExternalNodes;
}

Node **
GradientInelasticBeamColumn2d::getNodePtrs(void)
{
	return theNodes;
}

int
GradientInelasticBeamColumn2d::getNumDOF(void)
{
	return 6;
}

// Definition of Methods Dealing with Element's State
int
GradientInelasticBeamColumn2d::commitState(void)
{
	int err = 0;

	// Element commitState()
	if ((err = this->Element::commitState()))
		opserr << "WARNING! GradientInelasticBeamColumn2d::commitState() - element: " << this->getTag() << " - failed in committing base class\n";

	// Record [H] Diagonal Elements
	for (int i = 0; i < (numSections * secOrder); i++) {
		(*hh)(i) = (*H)(i, i);
	}

	// Commit Section State Variables
	for (int i = 0; i < numSections; i++) {
		err += sections[i]->commitState();
		d_sec_commit[i] = d_sec[i];
	}

	*d_tot_commit = *d_tot;
	*d_nl_tot_commit = *d_nl_tot;

	*F_ms_commit = *F_ms;

	// Commit Coordinate Transformation Object
	if ((err = crdTransf->commitState()))
		opserr << "WARNING! GradientInelasticBeamColumn2d::commitState() - element: " << this->getTag() << " - coordinate transformation object failed to commit\n";

	// Commit Jacobian Matrix
	*J_commit = *J;

	// Commit Element State Variables
	Q_commit = Q;

	// Set iters vector and zero iterNo and strIterNo
	totStrIterNo--;
	iters(0) = totStrIterNo;
	iters(1) = strIterNo;
	iters(2) = iterNo;

	iterNo = 0;
	strIterNo = 0;

	// Update {max_trial_change}
	commitNo++;

	if (correctionControl && (maxEpsInc == 0.0))
		for (int i = 0; i < numSections * secOrder + 3; i++)
			(*max_trial_change)(i) = ((commitNo - 1.0) * (*max_trial_change)(i)+fabs((*trial_change)(i))) / commitNo;

	// complete committing the variables

	return err;
}

int
GradientInelasticBeamColumn2d::revertToLastCommit(void)
{
	int err = 0;

	// Revert Section State Variables to Last Committed State
	for (int i = 0; i < numSections; i++) {
		err += sections[i]->revertToLastCommit();

		d_sec[i] = d_sec_commit[i];
		sections[i]->setTrialSectionDeformation(d_sec[i]);
	}

	*d_tot = *d_tot_commit;
	*d_nl_tot = *d_nl_tot_commit;

	// Revert Coordinate Transformation Object to Last Committed State
	if ((err = crdTransf->revertToLastCommit()))
		opserr << "WARNING! GradientInelasticBeamColumn2d::revertToLastCommit() - element: " << this->getTag() << " - coordinate transformation object failed to revert to last committed state\n";

	// Revert Element State Variables to Last Committed State
	Q = Q_commit;

	// Iteration Numbers
	initialFlag = 0;
	iterNo = 0;
	strIterNo = 0;
	iters.Zero();

	// complete reverting the variables

	return err;
}

int
GradientInelasticBeamColumn2d::revertToStart(void)
{
	int err = 0;

	// Revert Section State Variables to Start
	for (int i = 0; i < numSections; i++) {
		err += sections[i]->revertToStart();

		d_sec[i].Zero();
	}

	d_tot->Zero();
	d_tot_commit->Zero();

	d_nl_tot->Zero();
	d_nl_tot_commit->Zero();

	// Revert Coordinate Transformation Object to Start
	if ((err = crdTransf->revertToStart()))
		opserr << "WARNING! GradientInelasticBeamColumn2d::revertToStart() - element: " << this->getTag() << " - coordinate transformation object failed to revert to start\n";

	// Revert Element State Variables to Start
	Q.Zero();
	Q_commit.Zero();

	// Zero iteration numbers
	totStrIterNo = 0;
	iterNo = 0;
	strIterNo = 0;
	commitNo = 0;

	initialFlag = 0;
	return err;
}

void
GradientInelasticBeamColumn2d::assembleMatrix(Matrix& A, const Matrix& B, int rowStart, int rowEnd, int colStart, int colEnd, double fact)
{
	int rowsNo = rowEnd - rowStart + 1;
	int colsNo = colEnd - colStart + 1;

	if (B.noRows() != rowsNo)
		opserr << "ERROR! GradientInelasticBeamColumn2d::assembleMatrix() - element: " << this->getTag() << " - incompatible number of rows to assemble\n";

	if (B.noCols() != colsNo)
		opserr << "ERROR! GradientInelasticBeamColumn2d::assembleMatrix() - element: " << this->getTag() << " - incompatible number of columns to assemble\n";

	if ((A.noRows() - 1) < rowEnd)
		opserr << "ERROR! GradientInelasticBeamColumn2d::assembleMatrix() - element: " << this->getTag() << " - receiving matrix has less rows than needed\n";

	if ((A.noCols() - 1) < colEnd)
		opserr << "ERROR! GradientInelasticBeamColumn2d::assembleMatrix() - element: " << this->getTag() << " - receiving matrix has less columns than needed\n";

	int i = 0;
	int j = 0;

	for (int row = rowStart; row <= rowEnd; row++) {
		for (int col = colStart; col <= colEnd; col++) {
			A(row, col) = fact * B(i, j);
			j++;
		}

		j = 0;
		i++;
	}
}

void
GradientInelasticBeamColumn2d::assembleMatrix(Matrix &A, const Vector &B, int col, double fact)
{
	if (A.noRows() != B.Size())
		opserr << "ERROR! NonlocalBeamColumn2d::assembleMatrix - element: " << this->getTag() << " - incompatible matrix column number and vector size\n";

	for (int row = 0; row < B.Size(); row++)
		A(row, col) = fact * B(row);
}

void
GradientInelasticBeamColumn2d::assembleVector(Vector &A, const Vector &B, int rowStart, int rowEnd, double fact)
{
	int rowsNo = rowEnd - rowStart + 1;

	if (B.Size() != rowsNo)
		opserr << "ERROR! GradientInelasticBeamColumn2d::assembleVector() - element: " << this->getTag() << " - incompatible number of rows to assemble\n";

	if ((A.Size() - 1) < rowEnd)
		opserr << "ERROR! GradientInelasticBeamColumn2d::assembleVector() - element: " << this->getTag() << " - receiving matrix has less rows than needed\n";

	int i = 0;

	for (int row = rowStart; row <= rowEnd; row++) {
		A(row) = fact * B(i);

		i++;
	}
}

// Definition of Method to Solve for Nodal Forces
int
GradientInelasticBeamColumn2d::update(void)
{
	totStrIterNo++;

	// Get Section Behaviors
	const ID &code = sections[0]->getType();

	// Get Trial Nodal Basic Displacements
	crdTransf->update();
	const Vector &q_t = crdTransf->getBasicTrialDisp();		// target
	const Vector &q_inc = crdTransf->getBasicIncrDisp();		// increment from last committed step
	static Vector q(3);											// trial for each iteration (could be smaller than target)

	if (initialFlag != 0 && q_inc.Norm() <= DBL_EPSILON)
		return 0;

	// Declare Previous Step Variables for Subdivision
	Vector d_tot_prev(numSections * secOrder);
	Vector d_nl_tot_prev(numSections * secOrder);
	Vector F_ms_prev(numSections * secOrder);
	static Vector Q_prev(3);

	// Initialize Previous Sub-step Variables
	d_tot_prev = *d_tot_commit;
	d_nl_tot_prev = *d_nl_tot_commit;
	Q_prev = Q_commit;
	F_ms_prev = *F_ms_commit; // B_Q * Q_commit;

	// Declare Variables Requaired for Iterations
	static Vector dq(3);
	Vector dF_ms(numSections * secOrder);
	Vector d_inc_tot(numSections * secOrder);
	Vector trial_diff(3 + numSections * secOrder);
	Vector trial_old(3 + numSections * secOrder);
	Vector trial_new(3 + numSections * secOrder);
	Matrix K_ms(numSections * secOrder, numSections * secOrder);
	Matrix B_q_H_inv(3, numSections * secOrder);
	Matrix J_prev = *J_commit;

	bool bothConverged;
	bool converged = false;

	double dq_norm;
	double dF_ms_norm;

	// Subdivision of Basic Displacements in Case Convergence is not Achieved
	double r = 1.0;
	double r_fail;
	double r_commit = 0.0;
	double dr = 1.0;

	const double div = 10;
	const int maxDivNo = 7;

	for (int divNo = 0; divNo < maxDivNo; divNo++) {

		// Determine Current Trial Basic Displacements
		q = q_t - (1.0 - r) * q_inc;

		// Do Newton-Raphson Iteration to Achieve Compatible Forces
		for (int l = 0; l < 4; l++) {
			// Set Initial Solution Guesses
			*d_tot = d_tot_prev;
			Q = Q_prev;

			for (int iter = 1; iter <= maxIters; iter++) {

				// Set Trial Section Strains and Get New Section Forces
				for (int i = 0; i < numSections; i++) {
					d_sec[i].Extract(*d_tot, secOrder * i, 1.0);

					if (sections[i]->setTrialSectionDeformation(d_sec[i]) < 0) {
						opserr << "WARNING! GradientInelasticBeamColumn2d::update() - element: " << this->getTag() << " - section " << i << " failed in setTrialSectionDeformation\n";
						return -1;
					}

					this->assembleVector(*F_ms, sections[i]->getStressResultant(), i * secOrder, ((i + 1) * secOrder - 1), 1.0);
				}

				// Compute Material Section Strain Increments
				d_inc_tot = (*d_tot) - d_tot_prev;

				if (!cnstH) {
					//Update[H]
					*H = *H_init;

					double E_sec, eps_sec;

					for (int i = secOrder; i < (numSections - 1) * secOrder; i += secOrder) {
						E_sec = eps_sec = 0.0;

						for (int j = 0; j < secOrder; j++) {
							E_sec += (((*F_ms)(i + j) - F_ms_prev(i + j)) * d_inc_tot(i + j));
							eps_sec += fabs(d_inc_tot(i + j));
						}

						if (E_sec <= 0.0 && eps_sec > DBL_EPSILON) {
							for (int k = 0; k < secOrder; k++) {
								for (int j = 0; j < numSections * secOrder; j++)
									(*H)(i + k, j) = 0.0;

								(*H)(i + k, i + k) = 1.0;
							}
						}
					}

					// Invert [H]
					if (H->Invert(*H_inv) < 0) {
						opserr << "WARNING! GradientInelasticBeamColumn2d::update() - element: " << this->getTag() << " - could not invert [H]\n";
						return -1;
					}
				}

				// Compute Macroscopic Section Strains
				*d_nl_tot = d_nl_tot_prev + ((*H_inv) * d_inc_tot);

				// Check Convergence
				bothConverged = false;

				if ((this->qConvergence(iter, q, *d_nl_tot, dq, dq_norm) && this->fConvergence(iter, Q, dF_ms, dF_ms_norm)) || initialFlag == 0) {
					bothConverged = true;

					iterNo += iter;
					strIterNo++;

					iter = 10 * maxIters;
					l = 10;
				}
				else {
					// Assemble Jacobian Metrix
					switch (l) {
					case 0:
						if (!cnstH) {
							B_q_H_inv = (*B_q) * (*H_inv);
							this->assembleMatrix(*J, B_q_H_inv, numSections * secOrder, numSections * secOrder + 2, 3, numSections * secOrder + 2, -1.0);
						}
						else
							this->assembleMatrix(*J, *B_q_H_inv_init, numSections * secOrder, numSections * secOrder + 2, 3, numSections * secOrder + 2, -1.0);

						this->getSectionsTangentStiff(K_ms);
						this->assembleMatrix(*J, K_ms, 0, numSections * secOrder - 1, 3, numSections * secOrder + 2, -1.0);

						break;
					case 1:
						if (!cnstH) {
							B_q_H_inv = (*B_q) * (*H_inv);
							this->assembleMatrix(*J, B_q_H_inv, numSections * secOrder, numSections * secOrder + 2, 3, numSections * secOrder + 2, -1.0);
						}
						else
							this->assembleMatrix(*J, *B_q_H_inv_init, numSections * secOrder, numSections * secOrder + 2, 3, numSections * secOrder + 2, -1.0);

						this->getSectionsTangentStiff(K_ms);
						this->assembleMatrix(*J, K_ms, 0, numSections * secOrder - 1, 3, numSections * secOrder + 2, -1.0);

						break;
					case 2:
						*J = J_prev;

						break;
					case 3:
						*J = *J_init;

						break;
					}

					// Assemble Trial Variables
					this->assembleVector(trial_diff, dF_ms, 0, numSections * secOrder - 1, 1.0);
					this->assembleVector(trial_diff, dq, numSections * secOrder, numSections * secOrder + 2, 1.0);

					this->assembleVector(trial_old, Q, 0, 2, 1.0);
					this->assembleVector(trial_old, *d_tot, 3, numSections * secOrder + 2, 1.0);

					// Apply Newton-Raphson
					if (J->Solve(trial_diff, *trial_change) < 0) {
						opserr << "WARNING! GradientInelasticBeamColumn2d::update() - element: " << this->getTag() << " - could not invert Jacobian\n";
						return -1;
					}

					// Scale down trial_change if necessary
					if (correctionControl && !initialFlag && l > 0) {
						double gamma = 1.0;
						for (int i = 0; i < numSections * secOrder + 3; i++)
							if (fabs((*trial_change)(i)) > (*max_trial_change)(i))
								gamma = fmin((*max_trial_change)(i) / fabs((*trial_change)(i)), gamma);

						trial_new = trial_old - (gamma * (*trial_change));
					}
					else
						trial_new = trial_old - (*trial_change);

					Q.Extract(trial_new, 0, 1.0);
					d_tot->Extract(trial_new, 3, 1.0);
				} // else of (if (q_converged && F_ms_converged))
			} // for (int iter = 0; iter <= maxIters; iter++)
		} // for (int l = 0; l < 3; l++)

		if (bothConverged) {
			if (r == 1.0) {
				converged = true;
				break;
			}
			else {
				divNo--;
				r_commit = r;

				if (r > r_fail) {	// increase dr if passed previous failure point
					dr *= div;
					divNo--;
				}

				r += dr;

				if (r > 1.0)
					r = 1.0;

				d_tot_prev = *d_tot;
				d_nl_tot_prev = *d_nl_tot;
				F_ms_prev = *F_ms;
				Q_prev = Q;
				J_prev = (*J);
			}
		}
		else {
			r_fail = r;
			dr /= div;
			r = r_commit + dr;
		}
	} // for (int divNo = 0; divNo < maxDivNo; divNo++)

	// Check if Convergence has Occurred
	if (!converged) {
		opserr << "\nWARNING! GradientInelasticBeamColumn2d::update() - element: " << this->getTag() << " - failed to get compatible forces"
			<< "\ntarget basic displacements:    " << q_t(0) << ", " << q_t(1) << ", " << q_t(2)
			<< "\nbasic displacement increments: " << q_inc(0) << ", " << q_inc(1) << ", " << q_inc(2)
			<< "\ndq_norm: " << dq_norm << ", dF_ms_norm: " << dF_ms_norm << "\n\n";

		return -1;
	}

	// Set Convergence Indicators
	initialFlag = 1;
	return 0;
}

// Definition of Methods Dealing with Sections
void
GradientInelasticBeamColumn2d::getSectionsTangentStiff(Matrix &tStiff)
{
	tStiff.Zero();

	for (int i = 0; i < numSections; i++) {
		const Matrix &k_ms = sections[i]->getSectionTangent();
		//opserr << k_ms << endln;
		this->assembleMatrix(tStiff, k_ms, (i * secOrder), ((i + 1) * secOrder - 1), (i * secOrder), ((i + 1) * secOrder - 1), 1.0);
	}
}

void
GradientInelasticBeamColumn2d::getSectionsInitialStiff(Matrix &iStiff)
{
	iStiff.Zero();

	for (int i = 0; i < numSections; i++) {
		const Matrix &k_ms0 = sections[i]->getInitialTangent();
		this->assembleMatrix(iStiff, k_ms0, (i * secOrder), ((i + 1) * secOrder - 1), (i * secOrder), ((i + 1) * secOrder - 1), 1.0);
	}
}

// Definition of Convergence Test Methods
double
GradientInelasticBeamColumn2d::weightedNorm(const Vector &W, const Vector &V, bool sqRt)
{
	if (W.Size() != V.Size())
		opserr << "WARNING! GradientInelasticBeamColumnPF3d::weightedNorm() - element: " << this->getTag() << " - inequal number of elements in vectors\n";

	double sqSum = 0.0;
	for (int i = 0; i < V.Size(); i++)
		sqSum += V(i) * W(i) * V(i);

	if (sqRt)
		return sqrt(sqSum);
	else
		return sqSum;
}

bool
GradientInelasticBeamColumn2d::qConvergence(const int &iter, const Vector &qt, const Vector &dnl_tot, Vector &Dq, double &dqNorm)
{
	bool q_converged;

	Dq = qt - ((*B_q) * (*d_nl_tot));
	dqNorm = this->weightedNorm(k_init, Dq);

	if (iter < (maxIters / 3))
		q_converged = (dqNorm <= fmin(minTol * this->weightedNorm(k_init, qt), minTol * F_tol_q));
	else if (iter < (2 * maxIters / 3))
		q_converged = (dqNorm <= fmax(minTol * this->weightedNorm(k_init, qt), minTol * F_tol_q));
	else
		q_converged = (dqNorm <= fmax(maxTol * this->weightedNorm(k_init, qt), maxTol * F_tol_q));

	return q_converged;
}

bool
GradientInelasticBeamColumn2d::fConvergence(const int &iter, const Vector &Qt, Vector &DF_ms, double &dfNorm)
{
	bool F_ms_converged;

	Vector F_ms_trial = (*B_Q) * Q;
	DF_ms = F_ms_trial - (*F_ms);

	dfNorm = this->weightedNorm(*flex_ms_init, DF_ms);

	if (iter < (maxIters / 3))
		F_ms_converged = (dfNorm <= fmin(minTol * this->weightedNorm(*flex_ms_init, *F_ms), fmin(minTol * this->weightedNorm(*flex_ms_init, F_ms_trial), 100.0 * minTol * F_tol_f_ms)));
	else if (iter < (2 * maxIters / 3))
		F_ms_converged = (dfNorm <= fmax(minTol * this->weightedNorm(*flex_ms_init, *F_ms), fmax(minTol * this->weightedNorm(*flex_ms_init, F_ms_trial), 100.0 * minTol * F_tol_f_ms)));
	else
		F_ms_converged = (dfNorm <= fmax(maxTol * this->weightedNorm(*flex_ms_init, *F_ms), fmax(maxTol * this->weightedNorm(*flex_ms_init, F_ms_trial), 100.0 * maxTol * F_tol_f_ms)));

	return F_ms_converged;
}

// Definition of Methods Used to Determine Stiffness Matrices and Mass Matrix
const Matrix &
GradientInelasticBeamColumn2d::getBasicStiff(void)
{
	// Determine Element Stiffness Matrix in Basic System
	Matrix K_ms(numSections * secOrder, numSections * secOrder);
	Matrix K_ms_inv_BQ(numSections * secOrder, 3);
	static Matrix F(3, 3);       // flexibility matrix in the basic system
	static Matrix K(3, 3);       // stiffness matrix in the basic system

	this->getSectionsTangentStiff(K_ms);

	//opserr << K_ms << endln;
	//opserr << *B_Q << endln;
	
	if (K_ms.Solve(*B_Q, K_ms_inv_BQ) < 0) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::getBasicStiff() - element: " << this->getTag() << " - could not invert K_ms\n";
	}

	if (!cnstH)
		F = ((*B_q) * (*H_inv)) * K_ms_inv_BQ;
	else
		F = (*B_q_H_inv_init) * K_ms_inv_BQ;

	if (F.Invert(K) < 0) {
		opserr << "WARNING! GradientInelasticBeamColumn2d::getBasicStiff() - element: " << this->getTag() << " - could not invert element flexibility matrix\n";
	}

	return K;
}

const Matrix &
GradientInelasticBeamColumn2d::getTangentStiff(void)
{
	crdTransf->update();
	return crdTransf->getGlobalStiffMatrix(this->getBasicStiff(), Q);
}

const Matrix &
GradientInelasticBeamColumn2d::getInitialBasicStiff(void)
{
	// Determine Sections Initial Stiffness Matrices
	Matrix K_ms_init(numSections * secOrder, numSections * secOrder);

	this->getSectionsInitialStiff(K_ms_init);

	// Determine Element Stiffness Matrix in Basic System
	Matrix K_ms_inv_BQ(numSections * secOrder, 3);
	static Matrix F_init(3, 3);      // initial flexibility matrix in the basic system
	static Matrix K_init(3, 3);      // initial stiffness matrix in the basic system

	if (K_ms_init.Solve(*B_Q, K_ms_inv_BQ) < 0)
		opserr << "WARNING! GradientInelasticBeamColumn2d::getInitialBasicStiff() - element: " << this->getTag() << " - could not invert K_ms_init\n";

	F_init = (*B_q_H_inv_init) * K_ms_inv_BQ;

	if (F_init.Invert(K_init) < 0)
		opserr << "WARNING! GradientInelasticBeamColumn2d::getInitialBasicStiff() - element: " << this->getTag() << " - could not invert element initial flexibility matrix\n";

	return K_init;
}

const Matrix &
GradientInelasticBeamColumn2d::getInitialStiff(void)
{
	// Check for Quick Return
	if (K0 != 0)
		return *K0;

	// Transform to Global System and Return
	K0 = new Matrix(crdTransf->getInitialGlobalStiffMatrix(this->getInitialBasicStiff()));

	return *K0;
}

const Matrix &
GradientInelasticBeamColumn2d::getMass(void)
{
	theMatrix.Zero();

	return theMatrix;
}

// Definition of Methods to Get Forces
const Vector &
GradientInelasticBeamColumn2d::getResistingForce(void)
{
	double Q0[3];
	Vector Q0Vec(Q0, 3);
	Q0Vec.Zero();

	crdTransf->update();
	return crdTransf->getGlobalResistingForce(Q, Q0Vec);
}

const Vector &
GradientInelasticBeamColumn2d::getResistingForceIncInertia()
{
	// Compute the current resisting force
	theVector = this->getResistingForce();

	// add the damping forces if rayleigh damping
	if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
		theVector += this->getRayleighDampingForces();

	return theVector;
}

// Definition of Print Command
void
GradientInelasticBeamColumn2d::Print(OPS_Stream &s, int flag)
{
	s << "Element Tag: " << this->getTag() << endln;
	s << "Type: GradientInelasticBeamColumn2d" << endln;
	s << "Connected Node Tags: iNode " << connectedExternalNodes(0)
		<< ", jNode " << connectedExternalNodes(1) << endln;
	s << "Section Tag: " << sections[0]->getTag() << endln;
	s << "Number of Sections: " << numSections << endln;
	s << "Characteristic Length: " << lc << endln;
}

// Definition of Response Parameters
Response *
GradientInelasticBeamColumn2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	// Define and Initialize theResponse
	Response *theResponse = 0;

	output.tag("ElementOutput");
	output.attr("eleType", this->getClassType());
	output.attr("eleTag", this->getTag());
	output.attr("node1", connectedExternalNodes[0]);
	output.attr("node2", connectedExternalNodes[1]);

	// Global Forces
	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
		strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {
		output.tag("ResponseType", "Px_1");
		output.tag("ResponseType", "Py_1");
		output.tag("ResponseType", "Mz_1");
		output.tag("ResponseType", "Px_2");
		output.tag("ResponseType", "Py_2");
		output.tag("ResponseType", "Mz_2");

		theResponse = new ElementResponse(this, 1, theVector);
	}

	// Local Forces
	else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {
		output.tag("ResponseType", "N_ 1");
		output.tag("ResponseType", "My_1");
		output.tag("ResponseType", "Vy_1");
		output.tag("ResponseType", "N_2");
		output.tag("ResponseType", "Mz_2");
		output.tag("ResponseType", "Vy_2");

		theResponse = new ElementResponse(this, 2, theVector);
	}

	// Basic Forces
	else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {
		output.tag("ResponseType", "N_J");
		output.tag("ResponseType", "Mz_I");
		output.tag("ResponseType", "Mz_J");

		theResponse = new ElementResponse(this, 3, Vector(3));
	}

	// Nonlocal Strains
	else if (strcmp(argv[0], "nonlocalStrain") == 0 || strcmp(argv[0], "nonlocalStrains") == 0) {
		theResponse = new ElementResponse(this, 4, Vector(secOrder * numSections));
	}

	// Local Strains
	else if (strcmp(argv[0], "localStrain") == 0 || strcmp(argv[0], "localStrains") == 0) {
		theResponse = new ElementResponse(this, 5, Vector(secOrder * numSections));
	}

	// [H] Diagonal Elements
	else if (strcmp(argv[0], "Hdiagonal") == 0) {
		theResponse = new ElementResponse(this, 6, Vector(secOrder * numSections));
	}

	// Global Damping Forces
	else if (strcmp(argv[0], "dampingForce") == 0 || strcmp(argv[0], "dampingForces") == 0) {
		theResponse = new ElementResponse(this, 7, theVector);
	}

	// Iteration Number
	else if (strcmp(argv[0], "iterNo") == 0) {
		theResponse = new ElementResponse(this, 8, iters);
	}

	// Sections Responses
	else if (strstr(argv[0], "section") != 0) {

		if (argc > 1) {

			int sectionNum = atoi(argv[1]);

			if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

				double *secX = new double[numSections]; //double secX[50];
				beamIntegr->getSectionLocations(numSections, L, secX);

				output.tag("GaussPointOutput");
				output.attr("number", sectionNum);
				output.attr("eta", secX[sectionNum - 1] * L);

				if (strcmp(argv[2], "dsdh") != 0) {
					theResponse = sections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);
				}
				else {
					theResponse = new ElementResponse(this, 76, Vector(secOrder));
					Information &info = theResponse->getInformation();
					info.theInt = sectionNum;
				}

				output.endTag();

				if (secX)
					delete[] secX;
			}
		}
	}

	return theResponse;
}

// Definition of Method to Get Response Parameters
int
GradientInelasticBeamColumn2d::getResponse(int responseID, Information &eleInfo)
{
	switch (responseID) {
	case 1: // Global Forces
		return eleInfo.setVector(this->getResistingForce());

	case 2: // Local Forces
		theVector.Zero();

		theVector(0) = -Q(0);		// N_1
		theVector(3) = Q(0);		// N_2

		theVector(1) = (Q(1) + Q(2)) / L;		// V_1
		theVector(4) = -(Q(1) + Q(2)) / L;		// V_2

		theVector(2) = Q(1);		// M_1
		theVector(5) = Q(2);		// M_2

		return eleInfo.setVector(theVector);

	case 3: // Basic Forces
		return eleInfo.setVector(Q);

	case 4: // Nonlocal Strains
		return eleInfo.setVector(*d_nl_tot);

	case 5: // Local Strains
		return eleInfo.setVector(*d_tot);

	case 6:
		return eleInfo.setVector(*hh);

	case 7:
		return eleInfo.setVector(this->getRayleighDampingForces());

	case 8:
		return eleInfo.setVector(iters);

	default:
		return -1;
	}
}

// Definition of Method to Display Element
int
GradientInelasticBeamColumn2d::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** displayModes, int numModes)
{
    static Vector v1(3);
    static Vector v2(3);

    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);

    return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

// Definition of Methods Dealing with Parallel Processing
int
GradientInelasticBeamColumn2d::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID idData(17);
  idData(0) = this->getTag();
  idData(1) = connectedExternalNodes(0);
  idData(2) = connectedExternalNodes(1);  
  idData(3) = numSections;
  idData(4) = secOrder;
  idData(5) = maxIters;
  idData(6) = correctionControl;
  idData(7) = cnstH;
  idData(8) = crdTransf->getClassTag();
  int crdTransfDbTag  = crdTransf->getDbTag();
  if (crdTransfDbTag  == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag  != 0) 
      crdTransf->setDbTag(crdTransfDbTag);
  }
  idData(9) = crdTransfDbTag;
  idData(10) = beamIntegr->getClassTag();
  int beamIntDbTag  = beamIntegr->getDbTag();
  if (beamIntDbTag  == 0) {
    beamIntDbTag = theChannel.getDbTag();
    if (beamIntDbTag  != 0) 
      beamIntegr->setDbTag(beamIntDbTag);
  }
  idData(11) = beamIntDbTag;
  idData(12) = initialFlag;
  idData(13) = iterNo;
  idData(14) = strIterNo;
  idData(15) = totStrIterNo;  
  idData(16) = commitNo;
  
  if (theChannel.sendID(dbTag, commitTag, idData) < 0) {
    opserr << "GradientInelasticBeamColumn2d::sendSelf() - failed to send data ID" << endln;
    return -1;
  }
  
  static Vector data(10);
  data(0) = secLR1;
  data(1) = secLR2;
  data(2) = lc;
  data(3) = minTol;
  data(4) = maxTol;
  data(5) = F_tol_q; // Formed in setDomain
  data(6) = F_tol_f_ms; // Formed in setDomain
  data(7) = maxEpsInc;
  data(8) = maxPhiInc;
  data(9) = L;
  
  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "GradientInelasticBeamColumn2d::sendSelf() - failed to send data Vector" << endln;
    return -2;
  }

  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "GradientInelasticBeamColumn2d::sendSelf() - failed to send crdTranf" << endln;
    return -3;
  }      

  // send the beam integration
  if (beamIntegr->sendSelf(commitTag, theChannel) < 0) {
    opserr << "GradientInelasticBeamColumn2d::sendSelf() - failed to send beamInt" << endln;
    return -4;
  }


  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //
  ID idSections(2*numSections);
  int loc = 0;
  for (int i = 0; i<numSections; i++) {
    int sectClassTag = sections[i]->getClassTag();
    int sectDbTag = sections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      sections[i]->setDbTag(sectDbTag);
    }

    idSections(loc) = sectClassTag;
    idSections(loc+1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0)  {
    opserr << "GradientInelasticBeamColumn2d::sendSelf() - failed to send ID data" << endln;
    return -5;
  }    

  //
  // send the sections
  //
  
  for (int i = 0; i<numSections; i++) {
    if (sections[i]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "GradientInelasticBeamColumn2d::sendSelf() - section " << 
	i << " failed to send itself" << endln;
      return -6;
    }
  }

  /*
  Matrix *B_q;	3 x (Np*order)	// Formed in setDomain
  Matrix *B_Q;	(Np*order) x 3  // Formed in setDomain
  Matrix *H_init;   (Np*order) x (Np*order) // Formed in setDomain
  Matrix *H_inv;    (Np*order) x (Np*order) // Formed in setDomain
  Matrix *B_q_H_inv_init; 3 x (Np*order) // Formed in setDomain	
  Matrix *J_init;    (3+Np*order) x (3+Np*order) // Formed in setDomain
  Vector k_init;     3 // Formed in setDomain
  Vector *flex_ms_init;	Np*order // Formed in setDomain
  Vector *trial_change; 3+Np*order // Formed in setDomain
  Vector *max_trial_change; 3+Np*order // Formed in setDomain

  Vector Q_commit;      3
  Vector *d_tot_commit; Np*order
  Vector *d_nl_tot_commit; Np*order
  Vector *F_ms_commit;	Np*order
  Vector *hh;    Np*order
  Matrix *H;	(Np*order) x (Np*order)
  Matrix *J_commit; (3+Np*order) x (3+Np*order)

  Vector *d_sec_commit;  send the Np vectors
  */

  int totalSectionOrder = numSections*secOrder;
  Vector elementData(4*totalSectionOrder +
		     totalSectionOrder*totalSectionOrder +
		     (3 + totalSectionOrder)*(3 + totalSectionOrder) +
		     3);

  for (int i = 0; i < totalSectionOrder; i++) {
    elementData(i                      ) = (*d_tot_commit)(i);
    elementData(i +   totalSectionOrder) = (*d_nl_tot_commit)(i);
    elementData(i + 2*totalSectionOrder) = (*F_ms_commit)(i);
    elementData(i + 3*totalSectionOrder) = (*hh)(i);            
  }
  loc = 4*totalSectionOrder;
  for (int i = 0; i < totalSectionOrder; i++) {
    for (int j = 0; j < totalSectionOrder; j++) {
      elementData(loc++) = (*H)(i,j);
    }
  }
  for (int i = 0; i < 3+totalSectionOrder; i++) {
    for (int j = 0; j < 3+totalSectionOrder; j++) {
      elementData(loc++) = (*J_commit)(i,j);
    }
  }  
  elementData(loc++) = Q_commit(0);
  elementData(loc++) = Q_commit(1);
  elementData(loc++) = Q_commit(2);  

  if (theChannel.sendVector(dbTag, commitTag, elementData) < 0) {
    opserr << "GradientInelasaticBeamColumn2d::sendSelf() - failed to send elementData Vector" << endln;
    return -7;
  }
  
  return 0;
}

int
GradientInelasticBeamColumn2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID idData(17);
  if (theChannel.recvID(dbTag, commitTag, idData) < 0) {
    opserr << "GradientInelasticBeamColumn2d::recvSelf() - failed to receive data ID" << endln;
    return -1;
  }

  this->setTag(idData(0));
  connectedExternalNodes(0) = idData(1);
  connectedExternalNodes(1) = idData(2);  
  int nSect = idData(3);
  secOrder = idData(4);
  maxIters = idData(5);
  correctionControl = idData(6);
  cnstH = idData(7);
  int crdTransfClassTag = idData(8);
  int crdTransfDbTag = idData(9);
  int beamIntClassTag = idData(10);
  int beamIntDbTag = idData(11);  
  initialFlag = idData(12);
  iterNo = idData(13);
  strIterNo = idData(14);
  totStrIterNo = idData(15);
  commitNo = idData(16);
  
  static Vector data(10);
  if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
    opserr << "GradientInelasticBeamColumn2d::recvSelf() - failed to receive data Vector" << endln;
    return -2;
  }

  secLR1 = data(0);
  secLR2 = data(1);
  lc = data(2);
  minTol = data(3);
  maxTol = data(4);
  F_tol_q = data(5);
  F_tol_f_ms = data(6);
  maxEpsInc = data(7);
  maxPhiInc = data(8);  
  L = data(9);

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
    if (crdTransf != 0)
      delete crdTransf;
    
    crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);
    
    if (crdTransf == 0) {
      opserr << "GradientInelasticBeamColumn2d::recvSelf() - failed to obtain a CrdTransf object with classTag " <<
	crdTransfClassTag << endln;
      exit(-1);
    }
  }
  crdTransf->setDbTag(crdTransfDbTag);
  
  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "GradientInelasticBeamColumn2d::sendSelf() - failed to recv crdTranf" << endln;
    return -3;
  }      
  
  // create a new beamInt object if one needed
  if (beamIntegr == 0 || beamIntegr->getClassTag() != beamIntClassTag) {
    if (beamIntegr != 0)
      delete beamIntegr;
    
    beamIntegr = theBroker.getNewBeamIntegration(beamIntClassTag);
    
    if (beamIntegr == 0) {
      opserr << "GradientInelasticBeamColumn2d::recvSelf() - failed to obtain the beam integration object with classTag" <<
	beamIntClassTag << endln;
      exit(-1);
    }
  }
  
  beamIntegr->setDbTag(beamIntDbTag);
  
  // invoke recvSelf on the beamInt object
  if (beamIntegr->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "GradientInelasticBeamColumn2d::sendSelf() - failed to recv beam integration" << endln;
    return -4;
  }      


  ID idSections(2*nSect);

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0)  {
    opserr << "DispBeamColumn2d::recvSelf() - failed to recv ID data\n";
    return -5;
  }    

  //
  // now receive the sections
  //
  
  if (numSections != nSect) {

    //
    // we do not have correct number of sections, must delete the old and create
    // new ones before can recvSelf on the sections
    //

    // delete the old
    if (numSections != 0) {
      for (int i=0; i<numSections; i++)
	delete sections[i];
      delete [] sections;
    }

    // create a new array to hold pointers
    sections = new SectionForceDeformation *[nSect];
    if (sections == 0) {
      opserr << "GradientInelasticBeamColumn2d::recvSelf() - out of memory creating sections array of size " <<
	nSect << endln;
      return -5;
    }    

    // create a section and recvSelf on it
    numSections = nSect;
    int loc = 0;
    for (int i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;
      sections[i] = theBroker.getNewSection(sectClassTag);
      if (sections[i] == 0) {
	opserr << "GradientInelasticBeamColumn2d::recvSelf() - Broker could not create Section of class type " <<
	  sectClassTag << endln;
	exit(-1);
      }
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "GradientInelasticBeamColumn2d::recvSelf() - section " << i << " failed to recv itself" << endln;
	return -5;
      }     
    }

  } else {

    // 
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //
    
    int loc = 0;
    for (int i=0; i<numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag = idSections(loc+1);
      loc += 2;

      // check of correct type
      if (sections[i]->getClassTag() !=  sectClassTag) {
	// delete the old section[i] and create a new one
	delete sections[i];
	sections[i] = theBroker.getNewSection(sectClassTag);
	if (sections[i] == 0) {
	  opserr << "GradientInelasticBeamColumn2d::recvSelf() - Broker could not create Section of class type " <<
	    sectClassTag << endln;
	  exit(-1);
	}
      }

      // recvSelf on it
      sections[i]->setDbTag(sectDbTag);
      if (sections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
	opserr << "GradientInelasticBeamColumn2d::recvSelf() - section " << i << " failed to recv itself" << endln;
	return -5;
      }     
    }
  }

  //opserr << "recvSelf " << numSections << endln;
  this->setSectionPointers();
  
  int totalSectionOrder = numSections*secOrder;
  Vector elementData(4*totalSectionOrder +
		     totalSectionOrder*totalSectionOrder +
		     (3 + totalSectionOrder)*(3 + totalSectionOrder) +
		     3);

  if (theChannel.recvVector(dbTag, commitTag, elementData) < 0) {
    opserr << "GradientInelasticBeamColumn2d::recvSelf() - failed to receive elementData Vector" << endln;
    return -6;
  }

  for (int i = 0; i < totalSectionOrder; i++) {
    (*d_tot_commit)(i)    = elementData(i);
    (*d_nl_tot_commit)(i) = elementData(i +   totalSectionOrder);
    (*F_ms_commit)(i)     = elementData(i + 2*totalSectionOrder);
    (*hh)(i)              = elementData(i + 3*totalSectionOrder);
  }
  int loc = 4*totalSectionOrder;
  for (int i = 0; i < totalSectionOrder; i++) {
    for (int j = 0; j < totalSectionOrder; j++) {
      (*H)(i,j) = elementData(loc++);
    }
  }
  for (int i = 0; i < 3+totalSectionOrder; i++) {
    for (int j = 0; j < 3+totalSectionOrder; j++) {
      (*J_commit)(i,j) = elementData(loc++);
    }
  }  
  Q_commit(0) = elementData(loc++);
  Q_commit(1) = elementData(loc++);
  Q_commit(2) = elementData(loc++);
  
  return 0;
}
