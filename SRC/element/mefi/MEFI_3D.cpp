// Code written/implemented by:	Carlos Lopez Olea (carlos.lopez.o@ug.uchile.cl)
//
// User documentation available at: https://github.com/carloslopezolea/MEFI
//
// Created: 01/2026
//
// Description: The three-dimensional Membrane Fiber element (MEFI_3D) is a four-node element with six degrees of freedom (DOFs) per node: three translational DOFs and three rotational DOFs.
// The in-plane response is based on the MEFI_2D formulation, whereas the out-of-plane response follows the Kirchhoff plate formulation with four integration points per element. Both behaviors are formulated independently, 
// providing an uncoupled representation of membrane and bending actions in reinforced concrete walls.
//
// Reference:
// 1.- Lopez, C. N., Rojas, F., & Massone, L. M. (2022). Membrane fiber element for reinforced concrete walls - the benefits of macro and micro modeling approaches. Engineering Structures, 254, 113819.
// 2.- Suquillo, B., Rojas, F., López, C. et al. MEFI-3D: A membrane fiber element for non-planar reinforced concrete structural walls. Bull Earthquake Eng 24, 211–238 (2026).
//
// Source: /usr/local/cvs/OpenSees/SRC/element/MEFI/MEFI_3D.cpp
//
// Rev: 1.0                                                       

#include <MEFI_3D.h>
#include <Node.h>
#include <SectionForceDeformation.h>
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
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <elementAPI.h>
#include <DummyStream.h>

void* OPS_MEFI_3D()
{
    //pointer to a element that will be returned   
	Element* theEle = 0;

	//check model dimensions and nodal dofs 
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
	if (ndm != 3 || ndf != 6) {
		opserr << "WARNING element MEFI_3D: invalid model dimensions and/or nodal DOFs" << endln;
		return 0;
    }
	
	//check number of arguments provided
    if (OPS_GetNumRemainingInputArgs() < 8) {
		opserr << "WARNING element MEFI_3D: not enough args provided, want: element MEFI_3D eleTag iNode jNode kNode lNode numFibers -thick -sec" << endln;
		return 0;
    } 
	
	//create array to store integer data
	int iData[6];

	//check element tag and add to integer data
	int numData = 1;
	if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
		opserr << "WARNING element MEFI_3D: invalid integer element tag" << endln;
		return 0;
	}

	//check element nodes and add to integer data
	numData = 4;
	if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
		opserr << "WARNING element MEFI_3D: invalid integer node tag for element with tag " << iData[0] << endln;
		return 0;
	}

	//check element number of fibers and add to integer data
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[5]) != 0) {
		opserr << "WARNING element MEFI_3D: invalid integer numFibers for element with tag " << iData[0] << endln;
		return 0;
	}
	
	//create array to store double data
	double dData[2];
	dData[0] = 0.63;
	dData[1] = 0.25;

	//create pointers to store array data
	const char* str = 0;
	double* theWidth = new double[iData[5]];
	int* secTags = new int[iData[5]];
	SectionForceDeformation** theSec = new SectionForceDeformation * [iData[5]];
	
	int numArgs = OPS_GetNumRemainingInputArgs();
	while (numArgs > 0) {
		str = OPS_GetString();
		if (strcmp(str, "-width") == 0) {
			numData = iData[5];
			if (OPS_GetDoubleInput(&numData, theWidth) != 0) {
				opserr << "WARNING element MEFI_3D: invalid width value for element with tag " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-sec") == 0) {
			numData = iData[5];
			if (OPS_GetIntInput(&numData, secTags) != 0) {
				opserr << "WARNING element MEFI_3D: invalid section tag for element with tag " << iData[0] << endln;
				return 0;
			}
			for (int i = 0; i < iData[5]; i++) {
				theSec[i] = 0;
				theSec[i] = OPS_getSectionForceDeformation(secTags[i]);
				if (theSec[i] == 0) {
					opserr << "WARNING element MEFI_3D: invalid section tag " << secTags[i] << " for element with tag " << iData[0] << endln;
					return 0;
				}
			}
		}
		else if (strcmp(str, "-thickMod") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
				opserr << "WARNING element MEFI_3D: invalid thickness modifier for element with tag " << iData[0] << endln;
				return 0;
			}
			else if (dData[0] < 0.01 || dData[0] > 1.0) {
				opserr << "WARNING element MEFI_3D: invalid thickness modifier for element with tag " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-poisson") == 0) {
			numData = 1;
			if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
				opserr << "WARNING element MEFI_3D: invalid poisson ratio for element with tag " << iData[0] << endln;
				return 0;
			}
			else if (dData[1] < 0.01 || dData[1] > 1.0) {
				opserr << "WARNING element MEFI_3D: invalid poisson ratio for element with tag " << iData[0] << endln;
				return 0;
			}
		}
		numArgs = OPS_GetNumRemainingInputArgs();

	}//end while

    //now create the element and add it to the Domain
	theEle = new MEFI_3D(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], theSec, theWidth, dData[0], dData[1]);

	//cleanup dynamic memory
	if (theWidth != 0) { delete[] theWidth; }
	if (secTags != 0)  { delete[] secTags; }
	if (theSec != 0) { delete[] theSec; }

	//check if the element is not null
	if (theEle == 0) {
		opserr << "WARNING element MEFI_3D: ran out of memory creating element with tag " << iData[0] << endln;
		return 0;
	}
	
	else { return theEle; }

}//end OPS_MEFI_3D

Matrix MEFI_3D::K(24, 24);
Matrix MEFI_3D::KL(24, 24);
Matrix MEFI_3D::K12(12, 12);
Matrix MEFI_3D::K24(24, 24);
Vector MEFI_3D::P(24);
Vector MEFI_3D::PL(24);
Vector MEFI_3D::P12(12);
Matrix MEFI_3D::BSD(3, 12);
Matrix MEFI_3D::DB(3, 12);
double MEFI_3D::detJac;
int MEFI_3D::iMemb[12];
int MEFI_3D::iBend[12];
double MEFI_3D::qdtLocationsB[4][2];
double MEFI_3D::qdtWeightsB[4];

MEFI_3D::MEFI_3D(int tag, 
	int nd1, int nd2, int nd3, int nd4, 
	int numFibers,
	SectionForceDeformation** sec,
	double* width,
	double thickMod,
	double nu)
	:Element (tag, ELE_TAG_MEFI_3D), 
	connectedExternalNodes(4), nd1Crds(3), nd2Crds(3), nd3Crds(3), nd4Crds(3), Q(24),
	nd1CrdsL(3), nd2CrdsL(3), nd3CrdsL(3), nd4CrdsL(3), T(24, 24), Tt(3, 3), T6(6, 6), KP(12, 12),
	theSection(0), nip(numFibers), qdtLocationsM(numFibers, 2), qdtWeightsM(numFibers), lw(0), Ki(0),
	plateTangentNDM(3, 3), Eave(0), Tave(0), Nu(nu), tMod(thickMod), dispB(12), 
	BM(3, 12* numFibers), detJacM(numFibers), BP(3, 48), detJacP(4)
	
{
	//check width input
	if (width == 0) {
		opserr << "MEFI_3D::MEFI_3D() - Null width array passed.\n";
		exit(-1);
	}

	//create variables to allocate temporal data
	Vector x(nip);
	Vector b(nip);
	
	//assign some variables
	lw = 0;
	b.Zero();
	for (int i = 0; i < nip; i++) {
		b(i) = width[i];
		lw += b(i);
	}

	//calculate the location of the IPs along the element length
	x.Zero();
	for (int i = 0; i < nip; i++) {
		double sumb_i = 0.0;
		for (int j = 0; j < i + 1; j++) {
			sumb_i += b(j);
		}
		x(i) = (sumb_i - b(i) / 2.0) - lw / 2.0;
	}

	//calculate locations of quadrature points and weights
	for (int i = 0; i < nip; i++) {
		qdtLocationsM(i, 0) = x(i) * (2 / lw);
		qdtLocationsM(i, 1) = 0.001;
		qdtWeightsM(i) = b(i) * (2 / lw) * 2;
	}

	// we calculate locations of quadrature points and weights for bending behavior
	qdtLocationsB[0][0] = -1 / sqrt(3);	qdtLocationsB[0][1] = -1 / sqrt(3);
	qdtLocationsB[1][0] = +1 / sqrt(3);	qdtLocationsB[1][1] = -1 / sqrt(3);
	qdtLocationsB[2][0] = +1 / sqrt(3);	qdtLocationsB[2][1] = +1 / sqrt(3);
	qdtLocationsB[3][0] = -1 / sqrt(3);	qdtLocationsB[3][1] = +1 / sqrt(3);

	
	for (int i = 0; i < 4; i++) {
		qdtWeightsB[i] = 1;
	} // end for i

	for (int i = 0; i < 4; i++) {
		iMemb[3 * i + 0] = 6 * i + 0;
		iMemb[3 * i + 1] = 6 * i + 1;
		iMemb[3 * i + 2] = 6 * i + 5;
		iBend[3 * i + 0] = 6 * i + 2;
		iBend[3 * i + 1] = 6 * i + 3;
		iBend[3 * i + 2] = 6 * i + 4;
	}

	//allocate arrays of pointers to sections
	theSection = new SectionForceDeformation * [nip];
	if (theSection == 0) {
		opserr << "MEFI_3D::MEFI_3D() - failed allocate section model pointer" << endln;
		exit(-1);
	}

	//get copies of the sections
	for (int i = 0; i < nip; i++) {
		if (sec[i] == 0) {
			opserr << "MEFI_3D::MEFI_3D() - null section pointer passed" << endln;
			exit(-1);
		}
		theSection[i] = sec[i]->getCopy();
		if (theSection[i] == 0) {
			opserr << "MEFI_3D::MEFI_3D() - failed to copy section" << endln;
			exit(-1);
		}
	}

	//calculate bending parameters
	Eave = 0.0;
	Tave = 0.0;
	DummyStream theDummyStream;
	char aa[80] = "getBendingParameters";
	const char* argv[1];
	argv[0] = aa;

	for (int i = 0; i < nip; i++) {
		Response* theResponse = theSection[i]->setResponse(argv, 1, theDummyStream);
		if (theResponse == 0) {
			opserr << "MEFI_3D::MEFI_3D() - failled to get out of plane parameters for element with tag " << this->getTag() << endln;
			exit(-1);
		}
		theResponse->getResponse();
		Information& theInfoPar = theResponse->getInformation();
		const Vector& BendPar = theInfoPar.getData();
		Eave += (BendPar[0] * b(i)) / lw;
		Tave += (BendPar[1] * b(i)) / lw * tMod;
		delete theResponse;
	}

    //set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;

	//set node pointer to null
	for (int i = 0; i < 4; i++) { theNodes[i] = 0; }

}//end MEFI_3D

MEFI_3D::MEFI_3D()
	:Element(0, ELE_TAG_MEFI_3D),
	connectedExternalNodes(4), nd1Crds(3), nd2Crds(3), nd3Crds(3), nd4Crds(3), Q(24),
	nd1CrdsL(3), nd2CrdsL(3), nd3CrdsL(3), nd4CrdsL(3), T(24, 24), Tt(3, 3), T6(6, 6), KP(12, 12),
	theSection(0), nip(0), qdtLocationsM(1, 2), qdtWeightsM(0), lw(0), Ki(0),
	plateTangentNDM(3, 3), Eave(0), Tave(0), Nu(0), tMod(0), dispB(12),
	BM(3, 12), detJacM(0), BP(3, 12), detJacP(0)
{
	for (int i = 0; i < 4; i++) { theNodes[i] = 0; }

}//end MEFI_3D

MEFI_3D::~MEFI_3D()
{    
	for (int i = 0; i < nip; i++) { delete theSection[i]; }
	if (theSection != 0) { delete[] theSection; }
	if (Ki != 0) { delete Ki; }

}//end ~MEFI_3D

int
MEFI_3D::getNumExternalNodes() const
{
    return 4;

}//end getNumExternalNodes

const ID&
MEFI_3D::getExternalNodes()
{
    return connectedExternalNodes;

}//end getExternalNodes

Node **
MEFI_3D::getNodePtrs(void)
{
  return theNodes;

}//end getNodePtrs

int
MEFI_3D::getNumDOF()
{
    return 24;

}//end getNumDOF

void
MEFI_3D::setDomain(Domain *theDomain)
{
	//check domain is not null - invoked when object removed from a domain
	if (theDomain == 0) {
		theNodes[0] = 0;
		theNodes[1] = 0;
		theNodes[2] = 0;
		theNodes[3] = 0;
		return;
	}
	//first ensure external nodes exist in domain and set the node pointers
	int Nd1 = connectedExternalNodes(0);
	int Nd2 = connectedExternalNodes(1);
	int Nd3 = connectedExternalNodes(2);
	int Nd4 = connectedExternalNodes(3);

	theNodes[0] = theDomain->getNode(Nd1);
	theNodes[1] = theDomain->getNode(Nd2);
	theNodes[2] = theDomain->getNode(Nd3);
	theNodes[3] = theDomain->getNode(Nd4);

	//check the node pointer is not null
	if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 || theNodes[3] == 0) {
		opserr << "MEFI_3D::setDomain(): node not found in domain for element with tag " << this->getTag() << endln;
		return;
	}

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
    
    //check dofs compatibility
    if (dofNd1 != 6 || dofNd2 != 6 || dofNd3 != 6 || dofNd4 != 6) {
		opserr << "MEFI_3D::setDomain(): 6 dofs required at all nodes for element with tag " << this->getTag() << endln;
		return;
    }

	//get coordinates of end nodes
	nd1Crds = theNodes[0]->getCrds();
	nd2Crds = theNodes[1]->getCrds();
	nd3Crds = theNodes[2]->getCrds();
	nd4Crds = theNodes[3]->getCrds();

	//compute coordinate transformation matrix
	this->setTransformationMatrix();
	nd1CrdsL.addMatrixVector(0.0, Tt, nd1Crds, 1.0);
	nd2CrdsL.addMatrixVector(0.0, Tt, nd2Crds, 1.0);
	nd3CrdsL.addMatrixVector(0.0, Tt, nd3Crds, 1.0);
	nd4CrdsL.addMatrixVector(0.0, Tt, nd4Crds, 1.0);
	
	//calculate the element height and perform checks
	double h1 = pow(pow(nd4CrdsL(0) - nd1CrdsL(0), 2.0) + pow(nd4CrdsL(1) - nd1CrdsL(1), 2.0), 0.5);
	double h2 = pow(pow(nd3CrdsL(0) - nd2CrdsL(0), 2.0) + pow(nd3CrdsL(1) - nd2CrdsL(1), 2.0), 0.5);

	//check if element height is zero
	if ((h1 == 0.0) || (h2 == 0.0)) {
		opserr << "MEFI_3D::setDomain(): one of the sides is zero for element with tag " << this->getTag() << endln;
		exit(-1);
	}

	//check if element has constant height
	if ((h1 / h2 > 1.01) || (h1 / h2 < 0.99)) {
		opserr << "MEFI_3D::setDomain(): not constant height for element with tag " << this->getTag() << endln;
		exit(-1);
	}

	//calculate the element width and perform checks
	double b1 = pow(pow(nd2CrdsL(0) - nd1CrdsL(0), 2.0) + pow(nd2CrdsL(1) - nd1CrdsL(1), 2.0), 0.5);
	double b2 = pow(pow(nd3CrdsL(0) - nd4CrdsL(0), 2.0) + pow(nd3CrdsL(1) - nd4CrdsL(1), 2.0), 0.5);

	//check width of element
	if ((lw / b1 > 1.01) || (lw / b1 < 0.99) || (lw / b2 > 1.01) || (lw / b2 < 0.99)) {
		opserr << "MEFI_3D::setDomain(): nodes coordinates are not matched with fibers width for element with tag " << this->getTag() << endln;
		exit(-1);
	}

	//compute plate tangent
	this->setPlateTangent();

	//loop over the integration points
	for (int i = 0; i < nip; i++) {
		this->membraneFieldInterpolation(qdtLocationsM(i, 0), qdtLocationsM(i, 1));
		detJacM(i) = detJac;
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 12; k++) {
				BM(j , 12 * i + k) = BSD(j,k);
			}	
		}
	}

	//loop over the integration points
	for (int i = 0; i < 4; i++) {
		this->plateFieldInterpolation(qdtLocationsB[i][0], qdtLocationsB[i][1]);
		detJacP(i) = detJac;
		DB.addMatrixProduct(0.0, plateTangentNDM, BSD, 1.0);
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 12; k++) {
				BP(j, 12 * i + k) = BSD(j, k);
			}
		}
	}

    this->DomainComponent::setDomain(theDomain);

}//end setDomain

int
MEFI_3D::commitState()
{
	int retVal = 0;
	//call element commitState to do any base class stuff
	if ((retVal = this->Element::commitState()) != 0) {
		opserr << "MEFI_3D::commitState(): failed in base class for element with tag " << this->getTag() << endln;
	}
	//loop over the integration points and commit the material states
	for (int i = 0; i < nip; i++) {
		retVal += theSection[i]->commitState();
	}
	return retVal;

}//end commitState

int
MEFI_3D::revertToLastCommit()
{
	int retVal = 0;
	//loop over the integration points and revert to last committed state
	for (int i = 0; i < nip; i++) {
		retVal += theSection[i]->revertToLastCommit();
	}
	return retVal;

}//end revertToLastCommit

int
MEFI_3D::revertToStart()
{
	int retVal = 0;
	//loop over the integration points and revert states to start
	for (int i = 0; i < nip; i++) {
		retVal += theSection[i]->revertToStart();
	}
	return retVal;

}//end revertToStart


int
MEFI_3D::update()
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();

	// Store nodal displacements in global coordinates
	P.Zero();
	for (int i = 0; i < 6; i++) {
		P(i + 0) = disp1(i);
		P(i + 6) = disp2(i);
		P(i + 12) = disp3(i);
		P(i + 18) = disp4(i);
	}

	// Convert nodal displacements from global to local coordinate system
	PL.addMatrixVector(0.0, T, P, 1.0);

	// Get nodal displacements for membrane behavior
	for (int i = 0; i < 12; i++) {
		P12(i) = PL(iMemb[i]);
		dispB(i) = PL(iBend[i]);
	}

	// Strains at each quadrature point for membrane behavior
	Vector strainAtQuadraturePointM(3);
	strainAtQuadraturePointM.Zero();

	int ret = 0;
	//loop over the integration points
	for (int i = 0; i < nip; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 12; k++) {
				BSD(j, k) = BM(j, 12 * i + k);
			}
		}
		strainAtQuadraturePointM = BSD * P12;
		ret += theSection[i]->setTrialSectionDeformation(strainAtQuadraturePointM);
	}
	
	return ret;

}//end update


const Matrix&
MEFI_3D::getTangentStiff()
{

	K12.Zero();
	// loop over the integration points for membrane behavior
	for (int i = 0; i < nip; i++) {
		const Matrix& D = theSection[i]->getSectionTangent();
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 12; k++) {
				BSD(j, k) = BM(j, 12 * i + k);
			}	
		}
		K12.addMatrixTripleProduct(1.0, BSD, D, qdtWeightsM(i) * detJacM(i));
	}

	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			KL(iMemb[i], iMemb[j]) = K12(i, j);
			KL(iBend[i], iBend[j]) = KP(i, j);
		}
	}

	K.addMatrixTripleProduct(0.0, T, KL, 1.0); // Convert matrix from local to global cs

	return K;

}//end getTangentStiff


const Matrix&
MEFI_3D::getInitialStiff()
{
	if (Ki != 0) { return *Ki; }

	K12.Zero();
	// loop over the integration points for membrane behavior
	for (int i = 0; i < nip; i++) {
		const Matrix& D = theSection[i]->getInitialTangent();
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 12; k++) {
				BSD(j, k) = BM(j, 12 * i + k);
			}
		}
		K12.addMatrixTripleProduct(1.0, BSD, D, qdtWeightsM(i) * detJacM(i));
	}

	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 12; j++) {
			KL(iMemb[i], iMemb[j]) = K12(i, j);
			KL(iBend[i], iBend[j]) = KP(i, j);
		}
	}

	K.addMatrixTripleProduct(0.0, T, KL, 1.0); // Convert matrix from local to global cs

	Ki = new Matrix(K);
	return K;

}//end getInitialStiff

const Matrix&
MEFI_3D::getMass()
{
	K.Zero();
	double rhoi = 0;
	double massElement = 0;
	//we calculate the mass of the Element 
	for (int i = 0; i < nip; i++) {
		rhoi = theSection[i]->getRho();
		massElement += rhoi * qdtWeightsM(i) * detJacM(i);
	}

	//we assume lump mass for each node, only in the 1,2,3 direction
	for (int i = 0; i < 4; i++) {
		K(6 * i, 6 * i) = 0.25 * massElement;
		K(6 * i + 1, 6 * i + 1) = 0.25 * massElement;
		K(6 * i + 2, 6 * i + 2) = 0.25 * massElement;
	}

	return K;

}//end getMass

// N/A to this model - no element loads
void 
MEFI_3D::zeroLoad(void)
{
  	return;

}//end zeroLoad

// N/A to this model - no element loads
int 
MEFI_3D::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	return 0;

}//end addLoad

int 
MEFI_3D::addInertiaLoadToUnbalance(const Vector &accel)
{
	Q.Zero();
	K.Zero();
	//compute mass matrix
	K = this->getMass();

	//want to add ( - fact * M R * accel ) to unbalance
	for (int i = 0; i < 4; i++) {
		const Vector& Raccel = theNodes[i]->getRV(accel);
		Q(6 * i + 0) += -K(6 * i + 0, 6 * i + 0) * Raccel(0);
		Q(6 * i + 1) += -K(6 * i + 1, 6 * i + 1) * Raccel(1);
		Q(6 * i + 2) += -K(6 * i + 2, 6 * i + 2) * Raccel(2);
	}

	return 0;

}//end addInertiaLoadToUnbalance

const Vector&
MEFI_3D::getResistingForce()
{
	
	// loop over the integration points for membrane behavior
	P12.Zero();
	for (int i = 0; i < nip; i++) {
		const Vector& Stress = theSection[i]->getStressResultant();
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 12; k++) {
				BSD(j, k) = BM(j, 12 * i + k);
			}
		}

		P12.addMatrixTransposeVector(1.0, BSD, Stress, qdtWeightsM(i) * detJacM(i));
	}

	for (int i = 0; i < 12; i++) {
		PL(iMemb[i]) = P12(i);
	}

	// loop over the integration points for bending behavior
	P12.Zero();
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 12; k++) {
				BSD(j, k) = BP(j, 12 * i + k);
			}
		}

		DB.addMatrixProduct(0.0, plateTangentNDM, BSD, 1.0);
		Vector Stress(3);
		Stress.addMatrixVector(0.0, DB, dispB, 1.0);
		
		P12.addMatrixTransposeVector(1.0, BSD, Stress, qdtWeightsB[i] * detJacP(i));
	}

	for (int i = 0; i < 12; i++) {
		PL(iBend[i]) = P12(i);
	}

	P.addMatrixTransposeVector(0.0, T, PL, 1.0); // Convert matrix from local to global cs

	// Subtract other external nodal loads ... P_res = P_int - P_ext
	//P = P - Q;
	if (Q != 0) { P.addVector(1.0, Q, -1.0); }

	return P;

}//end getResistingForce

const Vector&
MEFI_3D::getResistingForceIncInertia()
{
	P.Zero();
	K.Zero();
	// compute the current resisting force
	P = this->getResistingForce();

	// add the damping forces if rayleigh damping
	if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
		P += this->getRayleighDampingForces();
	}

	// compute mass matrix
	K = this->getMass();

	// add inertia load
	for (int i = 0; i < 4; i++) {
		const Vector& accel = theNodes[i]->getTrialAccel();
		P(6 * i + 0) += K(6 * i + 0, 6 * i + 0) * accel(0);
		P(6 * i + 1) += K(6 * i + 1, 6 * i + 1) * accel(1);
		P(6 * i + 2) += K(6 * i + 2, 6 * i + 2) * accel(2);
	}

	return P;

}//end getResistingForceIncInertia

int MEFI_3D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;
	int dataTag = this->getDbTag();
	//MEFI_3D packs its data into a Vector and sends this to theChannel
	//along with its dbTag and the commitTag passed in the arguments
	static Vector data(6);
	data(0) = this->getTag();
	data(1) = nip;
	data(2) = alphaM;
	data(3) = betaK;
	data(4) = betaK0;
	data(5) = betaKc;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING MEFI_3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}
	//MEFI_3D sends the ids of its sections and nodes
	int secDbTag;
	static ID idData(2 * nip + 4);

	for (int i = 0; i < nip; i++) {
		idData(i) = theSection[i]->getClassTag();
		secDbTag = theSection[i]->getDbTag();
		//NOTE: we do have to ensure that the material has a database
		//tag if we are sending to a database channel.
		if (secDbTag == 0) {
			secDbTag = theChannel.getDbTag();
			if (secDbTag != 0)
				theSection[i]->setDbTag(secDbTag);
		}
		idData(i + nip) = secDbTag;
	}
	idData(2 * nip + 0) = connectedExternalNodes(0);
	idData(2 * nip + 1) = connectedExternalNodes(1);
	idData(2 * nip + 2) = connectedExternalNodes(2);
	idData(2 * nip + 3) = connectedExternalNodes(3);

	res += theChannel.sendID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING MEFI_3D::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}
	//MEFI_3D asks its material objects to send themselves
	for (int i = 0; i < nip; i++) {
		res += theSection[i]->sendSelf(commitTag, theChannel);
		if (res < 0) {
			opserr << "WARNING MEFI_3D::sendSelf() - " << this->getTag() << " failed to send its Material\n";
			return res;
		}
	}
	return res;
}

int MEFI_3D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	int dataTag = this->getDbTag();
	//MEFI_3D creates a Vector, receives the Vector and then sets the 
	//internal data with the data in the Vector
	static Vector data(6);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING MEFI_3D::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	nip = (int)data(1);
	alphaM = data(2);
	betaK = data(3);
	betaK0 = data(4);
	betaKc = data(5);

	static ID idData(2 * nip + 4);
	//MEFI_3D receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING MEFI_3D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	connectedExternalNodes(0) = idData(2 * nip + 0);
	connectedExternalNodes(1) = idData(2 * nip + 1);
	connectedExternalNodes(2) = idData(2 * nip + 2);
	connectedExternalNodes(3) = idData(2 * nip + 3);

	if (theSection == 0) {
		//allocate new materials
		theSection = new SectionForceDeformation *[nip];
		if (theSection == 0) {
			opserr << "MEFI_3D::recvSelf() - Could not allocate section array\n";
			return -1;
		}
		for (int i = 0; i < nip; i++) {
			int matClassTag = idData(i);
			int matDbTag = idData(i + nip);
			//allocate new material with the sent class tag
			theSection[i] = theBroker.getNewSection(matClassTag);
			if (theSection[i] == 0) {
				opserr << "MEFI_3D::recvSelf() - Broker could not create section of class type " << matClassTag << endln;
				return -1;
			}
			//now receive materials into the newly allocated space
			theSection[i]->setDbTag(matDbTag);
			res += theSection[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "MEFI_3D::recvSelf() - material " << i << "failed to recv itself\n";
				return res;
			}
		}
	}

	//materials exist, ensure materials of correct type and recvSelf on them
	else {
		for (int i = 0; i < nip; i++) {
			int matClassTag = idData(i);
			int matDbTag = idData(i + 4);
			//check that material is of the right type; if not,
			//delete it and create a new one of the right type
			if (theSection[i]->getClassTag() != matClassTag) {
				delete theSection[i];
				theSection[i] = theBroker.getNewSection(matClassTag);
				if (theSection[i] == 0) {
					opserr << "MEFI_3D::recvSelf() - material " << i << "failed to create\n";
					return -1;
				}
			}
			//receive the material
			theSection[i]->setDbTag(matDbTag);
			res += theSection[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "MEFI_3D::recvSelf() - material " << i << "failed to recv itself\n";
				return res;
			}
		}
	}
	
	return res;
}

void MEFI_3D::Print(OPS_Stream &s, int flag)
{
	if (flag == 0) {
		s << "MEFI_3D Element tag: " << this->getTag() << endln;
		s << "connected external nodes: " << connectedExternalNodes << endln;
		s << "number of fibers: " << nip << endln;
		s << "global resisting forces: " << this->getResistingForce() << endln;
		s << "Fiber responses: " << endln;
		for (int i = 0; i < nip; i++) {
			s << "Panel #: " << i + 1 << endln;
			s << "Section with tag: " << theSection[i]->getTag() << endln;
			theSection[i]->Print(s, flag);
		}
	}
	else if (flag == 1) {
		// does nothing
	}
}

int MEFI_3D::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	// get the end point display coords
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);
	theNodes[2]->getDisplayCrds(v3, fact, displayMode);
	theNodes[3]->getDisplayCrds(v4, fact, displayMode);

	// place values in coords matrix
	static Matrix coords(4, 3);
	for (int i = 0; i < 3; i++) {
		coords(0, i) = v1(i);
		coords(1, i) = v2(i);
		coords(2, i) = v3(i);
		coords(3, i) = v4(i);
	}

	// Display mode is positive:
	// display mode = 0 -> plot no contour
	static Vector values(4);
	for (int i = 0; i < 4; i++)
		values(i) = 0.0;

	// draw the polygon
	return theViewer.drawPolygon(coords, values, this->getTag());

}

Response* MEFI_3D::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "MEFI_3D");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes(0));
  output.attr("node2", connectedExternalNodes(1));
  output.attr("node3", connectedExternalNodes(2));
  output.attr("node4", connectedExternalNodes(3));

  char dataOut[40];

  //global forces
  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
	  strcmp(argv[0], "Force") == 0 || strcmp(argv[0], "Forces") == 0 ||
	  strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0 ||
	  strcmp(argv[0], "GlobalForce") == 0 || strcmp(argv[0], "GlobalForces") == 0) {

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

  // stresses at integration points
  else if ((strcmp(argv[0], "stresses") == 0) || (strcmp(argv[0], "stress") == 0) || 
			(strcmp(argv[0], "Stresses") == 0) || (strcmp(argv[0], "Stress") == 0)) {

	  for (int i = 1; i <= nip; i++) {
		  sprintf(dataOut, "sigma11_%d", i);
		  output.tag("ResponseType", dataOut);
		  sprintf(dataOut, "sigma22_%d", i);
		  output.tag("ResponseType", dataOut);
		  sprintf(dataOut, "sigma12_%d", i);
		  output.tag("ResponseType", dataOut);
	  } //end for i

	  theResponse = new ElementResponse(this, 2, Vector(3*nip));
  }

  // strains at integration points
  else if ((strcmp(argv[0], "strains") == 0) || (strcmp(argv[0], "strain") == 0) ||
			(strcmp(argv[0], "Strains") == 0) || (strcmp(argv[0], "Strain") == 0)) {

	  for (int i = 1; i <= nip; i++) {
		  sprintf(dataOut, "eps11_%d", i);
		  output.tag("ResponseType", dataOut);
		  sprintf(dataOut, "eps22_%d", i);
		  output.tag("ResponseType", dataOut);
		  sprintf(dataOut, "eps12_%d", i);
		  output.tag("ResponseType", dataOut);
	  } //end for i

	  theResponse = new ElementResponse(this, 3, Vector(3 * nip));
  }

 
  //material output
  else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0 ||
			strcmp(argv[0], "Material") == 0 || strcmp(argv[0], "IntegrPoint") == 0 ||
			strcmp(argv[0], "RCpanel") == 0 || strcmp(argv[0], "RCPanel") || 
			strcmp(argv[0], "RC_panel") || strcmp(argv[0], "RC_Panel") == 0) {

	  //check if correct # of arguments passed
	  if (argc < 3) {
		  opserr << "WARNING: Number of recorder input for section is: " << argc - 1 << "; should be at least 3.\n";
		  return 0;
	  }

	  int secNum = atoi(argv[1]);

	  output.tag("Material");
	  output.attr("number", secNum);

	  theResponse = theSection[secNum - 1]->setResponse(&argv[2], argc - 2, output);

  }

  output.endTag(); //elementOutput

  return theResponse;
}

int  MEFI_3D::getResponse(int responseID, Information &eleInfo)
{
	//global forces
	if (responseID == 1) {
		return eleInfo.setVector(this->getResistingForce());
	}

	//stresses at integration points
	else if (responseID == 2) {

		//loop over the integration points
		static Vector stresses(3*nip);
		int cnt = 0;
		for (int i = 0; i < nip; i++) {

			//get section stress response
			const Vector& sigma = theSection[i]->getStressResultant();
			stresses(cnt) = sigma(0);
			stresses(cnt + 1) = sigma(1);
			stresses(cnt + 2) = sigma(2);
			cnt += 3;
		}

		return eleInfo.setVector(stresses);

	}

	//strains at integration points
	else if (responseID == 3) {

		//loop over the integration points
		static Vector strains(3 * nip);
		int cnt = 0;
		for (int i = 0; i < nip; i++) {

			//get section stress response
			const Vector& eps = theSection[i]->getSectionDeformation();
			strains(cnt) = eps(0);
			strains(cnt + 1) = eps(1);
			strains(cnt + 2) = eps(2);
			cnt += 3;
		}

		return eleInfo.setVector(strains);

	}
	
	//nothing to do
	else { return -1; }
}

void MEFI_3D::membraneFieldInterpolation(double xi, double eta)
{
	// we calculate the Derivatives of the Shape Function at the natural coordinate (xi, eta)
	// we assume the Serendipity 8 Node Interpolation to calculate the Jacobian
	// we calculate the Shape function Derivatives Matrix w.r.t.the natural Coord System
	Matrix dShapeFunction(2, 8); 
	dShapeFunction(0, 0) = 0.25 * (1.0 - eta) * (2.0 * xi + eta);
	dShapeFunction(0, 1) = 0.25 * (1.0 - eta) * (2.0 * xi - eta);
	dShapeFunction(0, 2) = 0.25 * (1.0 + eta) * (2.0 * xi + eta);
	dShapeFunction(0, 3) = 0.25 * (1.0 + eta) * (2.0 * xi - eta);
	dShapeFunction(0, 4) = 0.5 * (-2.0 * xi) * (1.0 - eta);
	dShapeFunction(0, 5) = 0.5 * (1.0 - eta * eta);
	dShapeFunction(0, 6) = 0.5 * (-2.0 * xi) * (1.0 + eta);
	dShapeFunction(0, 7) = 0.5 * (-1.0) * (1.0 - eta * eta);
	dShapeFunction(1, 0) = 0.25 * (1.0 - xi) * (2.0 * eta + xi);
	dShapeFunction(1, 1) = 0.25 * (1.0 + xi) * (2.0 * eta - xi);
	dShapeFunction(1, 2) = 0.25 * (1.0 + xi) * (2.0 * eta + xi);
	dShapeFunction(1, 3) = 0.25 * (1.0 - xi) * (2.0 * eta - xi);
	dShapeFunction(1, 4) = 0.5 * (-1.0) * (1.0 - xi * xi);
	dShapeFunction(1, 5) = 0.5 * (-2.0 * eta) * (1.0 + xi);
	dShapeFunction(1, 6) = 0.5 * (1.0 - xi * xi);
	dShapeFunction(1, 7) = 0.5 * (-2.0 * eta) * (1.0 - xi);

	// we calculate the intermedial node coord for a 8 Node Serendipity Field Interpolation
	Matrix localCoord8Nodes(8, 2); 
	localCoord8Nodes(0, 0) = nd1CrdsL(0);						localCoord8Nodes(0, 1) = nd1CrdsL(1);
	localCoord8Nodes(1, 0) = nd2CrdsL(0);						localCoord8Nodes(1, 1) = nd2CrdsL(1);
	localCoord8Nodes(2, 0) = nd3CrdsL(0);						localCoord8Nodes(2, 1) = nd3CrdsL(1);
	localCoord8Nodes(3, 0) = nd4CrdsL(0);						localCoord8Nodes(3, 1) = nd4CrdsL(1);
	localCoord8Nodes(4, 0) = 0.5 * (nd1CrdsL(0) + nd2CrdsL(0)); localCoord8Nodes(4, 1) = 0.5 * (nd1CrdsL(1) + nd2CrdsL(1));
	localCoord8Nodes(5, 0) = 0.5 * (nd2CrdsL(0) + nd3CrdsL(0)); localCoord8Nodes(5, 1) = 0.5 * (nd2CrdsL(1) + nd3CrdsL(1));
	localCoord8Nodes(6, 0) = 0.5 * (nd3CrdsL(0) + nd4CrdsL(0)); localCoord8Nodes(6, 1) = 0.5 * (nd3CrdsL(1) + nd4CrdsL(1));
	localCoord8Nodes(7, 0) = 0.5 * (nd4CrdsL(0) + nd1CrdsL(0)); localCoord8Nodes(7, 1) = 0.5 * (nd4CrdsL(1) + nd1CrdsL(1));

	// we calculate the Jacobian Transformation Matrix
	Matrix JacobianMatrix(2, 2);
	JacobianMatrix.addMatrixProduct(0.0, dShapeFunction, localCoord8Nodes, 1.0);

	// we calculate the Jacobian inverse and determinant 
	detJac = JacobianMatrix(0, 0) * JacobianMatrix(1, 1) - JacobianMatrix(0, 1) * JacobianMatrix(1, 0);
	Matrix inverseJacobian(2, 2);
	inverseJacobian(0, 0) = +JacobianMatrix(1, 1) / detJac;
	inverseJacobian(0, 1) = -JacobianMatrix(0, 1) / detJac;
	inverseJacobian(1, 0) = -JacobianMatrix(1, 0) / detJac;
	inverseJacobian(1, 1) = +JacobianMatrix(0, 0) / detJac;

	// we calculate the J Matrix = [inverseJacobian zeros(2,2) ; zeros(2,2) inverseJacobian]
	Matrix jMatrix(4, 4); 
	jMatrix(0, 0) = inverseJacobian(0, 0);
	jMatrix(0, 1) = inverseJacobian(0, 1);
	jMatrix(1, 0) = inverseJacobian(1, 0);
	jMatrix(1, 1) = inverseJacobian(1, 1);
	jMatrix(2, 2) = inverseJacobian(0, 0);
	jMatrix(2, 3) = inverseJacobian(0, 1);
	jMatrix(3, 2) = inverseJacobian(1, 0);
	jMatrix(3, 3) = inverseJacobian(1, 1);

	// we calculate the C matrix = ShapeFunction Derivatives w.r.t. the natural Coord System (xi and eta) for the Field Interpolation
	Matrix cMatrix(4, 16); 
	cMatrix(0, 0)  = (0.5) * (-(0.5) + (0.75) * (eta) - 0.25 * (eta * eta * eta));
	cMatrix(0, 2)  = (0.5) * (0.25 - 0.25 * eta - 0.25 * (eta * eta) + 0.25 * (eta * eta * eta));
	cMatrix(0, 4)  = (0.5) * ((0.5) - (0.75) * (eta) +0.25 * (eta * eta * eta));
	cMatrix(0, 6)  = (0.5) * (-(0.25) + 0.25 * eta + 0.25 * (eta * eta) - 0.25 * (eta * eta * eta));
	cMatrix(0, 8)  = (0.5) * ((0.5) + (0.75) * (eta) - 0.25 * (eta * eta * eta));
	cMatrix(0, 10) = (0.5) * (0.25 + 0.25 * eta - 0.25 * (eta * eta) - 0.25 * (eta * eta * eta));
	cMatrix(0, 12) = (0.5) * (-(0.5) - (0.75) * (eta) + 0.25 * (eta * eta * eta));
	cMatrix(0, 14) = (0.5) * (-(0.25) - 0.25 * eta + 0.25 * (eta * eta) + 0.25 * (eta * eta * eta));

	cMatrix(1, 0)  = (0.5) * (-(0.75) + (0.75) * (eta * eta)) * (1.0 - xi);
	cMatrix(1, 2)  = (0.5) * (-(0.25) - (0.5) * eta + (0.75) * (eta * eta)) * (-1.0 + xi);
	cMatrix(1, 4)  = (0.5) * (-(0.75) + (0.75) * (eta * eta)) * (1.0 + xi);
	cMatrix(1, 6)  = (0.5) * (-(0.25) - (0.5) * eta + (0.75) * (eta * eta)) * (-1.0 - xi);
	cMatrix(1, 8)  = (0.5) * ((0.75) - (0.75) * (eta * eta)) * (1.0 + xi);
	cMatrix(1, 10) = (0.5) * (-(0.25) + (0.5) * eta + (0.75) * (eta * eta)) * (-1.0 - xi);
	cMatrix(1, 12) = (0.5) * ((0.75) - (0.75) * (eta * eta)) * (1.0 - xi);
	cMatrix(1, 14) = (0.5) * (-(0.25) + (0.5) * eta + (0.75) * (eta * eta)) * (-1.0 + xi);

	cMatrix(2, 1)  = (0.5) * (1.0 - eta) * (-(0.75) + (0.75) * (xi * xi));
	cMatrix(2, 3)  = (0.5) * (1.0 - eta) * (-(0.25) - 0.5 * xi + (0.75) * (xi * xi));
	cMatrix(2, 5)  = (0.5) * (1.0 - eta) * ((0.75) - (0.75) * (xi * xi));
	cMatrix(2, 7)  = (0.5) * (1.0 - eta) * (-(0.25) + 0.5 * xi + (0.75) * (xi * xi));
	cMatrix(2, 9)  = (0.5) * (1.0 + eta) * ((0.75) - (0.75) * (xi * xi));
	cMatrix(2, 11) = (0.5) * (1.0 + eta) * (-(0.25) + 0.5 * xi + (0.75) * (xi * xi));
	cMatrix(2, 13) = (0.5) * (1.0 + eta) * (-(0.75) + (0.75) * (xi * xi));
	cMatrix(2, 15) = (0.5) * (1.0 + eta) * (-(0.25) - 0.5 * xi + (0.75) * (xi * xi));

	cMatrix(3, 1)  = (0.5) * (-(0.5) + (0.75) * (xi) - (0.25) * (xi * xi * xi));
	cMatrix(3, 3)  = (0.5) * (-(0.25) + (0.25) * xi + (0.25) * (xi * xi) - (0.25) * (xi * xi * xi));
	cMatrix(3, 5)  = (0.5) * (-(0.5) - (0.75) * (xi) + (0.25) * (xi * xi * xi));
	cMatrix(3, 7)  = (0.5) * ((0.25) + (0.25) * (xi) - (0.25) * (xi * xi) - (0.25) * (xi * xi * xi));
	cMatrix(3, 9)  = (0.5) * ((0.5) + (0.75) * (xi) - (0.25) * (xi * xi * xi));
	cMatrix(3, 11) = (0.5) * (-(0.25) - (0.25) * xi + (0.25) * (xi * xi) + (0.25) * (xi * xi * xi));
	cMatrix(3, 13) = (0.5) * ((0.5) - (0.75) * (xi) + (0.25) * (xi * xi * xi));
	cMatrix(3, 15) = (0.5) * ((0.25) - (0.25) * xi - (0.25) * (xi * xi) + (0.25) * (xi * xi * xi));
	
	// we calculate the Transformation matrix from a 4 Node Blended interpolation with 4DOf per Node to a 4 Node with Rotational DOF
	double x21 = 0.5 * (nd2CrdsL(0) - nd1CrdsL(0));
	double x34 = 0.5 * (nd3CrdsL(0) - nd4CrdsL(0));
	double y41 = 0.5 * (nd4CrdsL(1) - nd1CrdsL(1));
	double y32 = 0.5 * (nd3CrdsL(1) - nd2CrdsL(1));

	Matrix trMatrix(16, 12); 
	trMatrix(0, 0) = 1.0;
	trMatrix(1, 1) = 1.0;
	trMatrix(2, 2) = y41;
	trMatrix(3, 2) = x21;
	trMatrix(4, 3) = 1.0;
	trMatrix(5, 4) = 1.0;
	trMatrix(6, 5) = y32;
	trMatrix(7, 5) = x21;
	trMatrix(8, 6) = 1.0;
	trMatrix(9, 7) = 1.0;
	trMatrix(10, 8) = y32;
	trMatrix(11, 8) = x34;
	trMatrix(12, 9) = 1.0;
	trMatrix(13, 10) = 1.0;
	trMatrix(14, 11) = y41;
	trMatrix(15, 11) = x34;

	// we create the Kinematic matrix between the Curvature and the Displacement derivatives
	Matrix aMatrix(3, 4); 
	aMatrix(0, 0) = 1.0;
	aMatrix(1, 3) = 1.0;
	aMatrix(2, 1) = 1.0;
	aMatrix(2, 2) = 1.0;

	// we calculate B = strain-displacement relationship matrix , B = (A*J*C*T)
	Matrix M_3_4(3, 4);
	M_3_4.addMatrixProduct(0.0, aMatrix, jMatrix, 1.0);
	Matrix M_3_16(3, 16);
	M_3_16.addMatrixProduct(0.0, M_3_4, cMatrix, 1.0);
	BSD.addMatrixProduct(0.0, M_3_16, trMatrix, 1.0);
	
} // end membraneFieldInterpolation

void MEFI_3D::plateFieldInterpolation(double xi, double eta)
{
	// we calculate the Jacobian using a Bilinear Interpolation
	// we calculate the Shape function Derivatives Matrix w.r.t. the natural Coord System
	Matrix dShapeFunction(2, 4);
	dShapeFunction(0, 0) = -0.25 * (1.0 - eta);
	dShapeFunction(0, 1) = +0.25 * (1.0 - eta);
	dShapeFunction(0, 2) = +0.25 * (1.0 + eta);
	dShapeFunction(0, 3) = -0.25 * (1.0 + eta);
	dShapeFunction(1, 0) = -0.25 * (1.0 - xi);
	dShapeFunction(1, 1) = -0.25 * (1.0 + xi);
	dShapeFunction(1, 2) = +0.25 * (1.0 + xi);
	dShapeFunction(1, 3) = +0.25 * (1.0 - xi);

	// we set the local Coord 4 Node Matrix
	Matrix localCoord4Nodes(4, 2);
	localCoord4Nodes(0, 0) = nd1CrdsL(0);	localCoord4Nodes(0, 1) = nd1CrdsL(1);
	localCoord4Nodes(1, 0) = nd2CrdsL(0);	localCoord4Nodes(1, 1) = nd2CrdsL(1);
	localCoord4Nodes(2, 0) = nd3CrdsL(0);	localCoord4Nodes(2, 1) = nd3CrdsL(1);
	localCoord4Nodes(3, 0) = nd4CrdsL(0);	localCoord4Nodes(3, 1) = nd4CrdsL(1);

	// we calculate the Jacobian Transformation Matrix
	Matrix JacobianMatrix(2, 2);
	JacobianMatrix.addMatrixProduct(0.0, dShapeFunction, localCoord4Nodes, 1.0);

	// we calculate the Jacobian inverse and determinant 
	detJac = JacobianMatrix(0, 0) * JacobianMatrix(1, 1) - JacobianMatrix(0, 1) * JacobianMatrix(1, 0);
	Matrix inverseJacobian(2, 2);
	inverseJacobian(0, 0) = +JacobianMatrix(1, 1) / detJac;
	inverseJacobian(0, 1) = -JacobianMatrix(0, 1) / detJac;
	inverseJacobian(1, 0) = -JacobianMatrix(1, 0) / detJac;
	inverseJacobian(1, 1) = +JacobianMatrix(0, 0) / detJac;

	// we calculate the Shape function Derivatives Matrix w.r.t. the natural Coord System for a Serendipity 8 Node Element
	Matrix dPsi(2, 8);
	dPsi(0, 0) = 0.25 * (1.0 - eta) * (2.0 * xi + eta);
	dPsi(0, 1) = 0.25 * (1.0 - eta) * (2.0 * xi - eta);
	dPsi(0, 2) = 0.25 * (1.0 + eta) * (2.0 * xi + eta);
	dPsi(0, 3) = 0.25 * (1.0 + eta) * (2.0 * xi - eta);
	dPsi(0, 4) = 0.5 * (-2.0 * xi) * (1.0 - eta);
	dPsi(0, 5) = 0.5 * (1.0 - eta * eta);
	dPsi(0, 6) = 0.5 * (-2.0 * xi) * (1.0 + eta);
	dPsi(0, 7) = 0.5 * (-1.0) * (1.0 - eta * eta);
	dPsi(1, 0) = 0.25 * (1.0 - xi) * (2.0 * eta + xi);
	dPsi(1, 1) = 0.25 * (1.0 + xi) * (2.0 * eta - xi);
	dPsi(1, 2) = 0.25 * (1.0 + xi) * (2.0 * eta + xi);
	dPsi(1, 3) = 0.25 * (1.0 - xi) * (2.0 * eta - xi);
	dPsi(1, 4) = 0.5 * (-1.0) * (1.0 - xi * xi);
	dPsi(1, 5) = 0.5 * (-2.0 * eta) * (1.0 + xi);
	dPsi(1, 6) = 0.5 * (1.0 - xi * xi);
	dPsi(1, 7) = 0.5 * (-2.0 * eta) * (1.0 - xi);

	// we set the dH_X and dH_Y matrix = ShapeFunction Derivatives w.r.t. the natural Coord System (xi and eta) for the Field Interpolation
	Matrix dH_X(2, 12);
	Matrix dH_Y(2, 12);

	// we calculate the a, b, c, d, e factors used to calculate the transformation from a serendipity 8 node
	int ij[4][2] = { 0 };
	ij[0][0] = 0;  ij[0][1] = 1;
	ij[1][0] = 1;  ij[1][1] = 2;
	ij[2][0] = 2;  ij[2][1] = 3;
	ij[3][0] = 3;  ij[3][1] = 0;

	Vector a(4);
	Vector b(4);
	Vector c(4);
	Vector d(4);
	Vector e(4);

	for (int k = 0; k < 4; k++) {
		double x_ij = localCoord4Nodes(ij[k][0], 0) - localCoord4Nodes(ij[k][1], 0);
		double y_ij = localCoord4Nodes(ij[k][0], 1) - localCoord4Nodes(ij[k][1], 1);
		double L_ij2 = (x_ij * x_ij) + (y_ij * y_ij);
		a(k) = -x_ij / L_ij2;
		b(k) = 0.75 * (x_ij * y_ij) / L_ij2;
		c(k) = (0.25 * (x_ij * x_ij) - 0.5 * (y_ij * y_ij)) / L_ij2;
		d(k) = -y_ij / L_ij2;
		e(k) = (-0.5 * (x_ij * x_ij) + 0.25 * (y_ij * y_ij)) / L_ij2;
	}

	int ik[4][2] = { 0 };
	ik[0][0] = 0;  ik[0][1] = 4;
	ik[1][0] = 1;  ik[1][1] = 5;
	ik[2][0] = 2;  ik[2][1] = 6;
	ik[3][0] = 3;  ik[3][1] = 7;

	int jk[4][2] = { 0 };
	jk[0][0] = 3;  jk[0][1] = 7;
	jk[1][0] = 0;  jk[1][1] = 4;
	jk[2][0] = 1;  jk[2][1] = 5;
	jk[3][0] = 2;  jk[3][1] = 6;

	// we proceed to calculate the derivatives of the FieldInterpolation w.r.t.the Natural Coordinate System
	for (int iNode = 0; iNode < 4; iNode++) {
		for (int iDir = 0; iDir < 2; iDir++) {
			// we calculate the derivatives
			dH_X(iDir, 3 * iNode + 0) = 1.5 * (a(ik[iNode][0]) * dPsi(iDir, ik[iNode][1]) - a(jk[iNode][0]) * dPsi(iDir, jk[iNode][1]));
			dH_X(iDir, 3 * iNode + 1) = b(ik[iNode][0]) * dPsi(iDir, ik[iNode][1]) + b(jk[iNode][0]) * dPsi(iDir, jk[iNode][1]);
			dH_X(iDir, 3 * iNode + 2) = dPsi(iDir, iNode) - c(ik[iNode][0]) * dPsi(iDir, ik[iNode][1]) - c(jk[iNode][0]) * dPsi(iDir, jk[iNode][1]);
			dH_Y(iDir, 3 * iNode + 0) = 1.5 * (d(ik[iNode][0]) * dPsi(iDir, ik[iNode][1]) - d(jk[iNode][0]) * dPsi(iDir, jk[iNode][1]));
			dH_Y(iDir, 3 * iNode + 1) = -dPsi(iDir, iNode) + e(ik[iNode][0]) * dPsi(iDir, ik[iNode][1]) + e(jk[iNode][0]) * dPsi(iDir, jk[iNode][1]);
			dH_Y(iDir, 3 * iNode + 2) = -b(ik[iNode][0]) * dPsi(iDir, ik[iNode][1]) - b(jk[iNode][0]) * dPsi(iDir, jk[iNode][1]);
		} // end for iDir
	} // end for iNode

	// we calculate the Matrix that Relate the Derivative of the Field Interpolation and the value at the nodes
	Matrix dH(4, 12); 
	for (int i = 0; i < 12; i++) {
		dH(0, i) = inverseJacobian(0, 0) * dH_X(0, i) + inverseJacobian(0, 1) * dH_X(1, i);
		dH(1, i) = inverseJacobian(1, 0) * dH_X(0, i) + inverseJacobian(1, 1) * dH_X(1, i);
		dH(2, i) = inverseJacobian(0, 0) * dH_Y(0, i) + inverseJacobian(0, 1) * dH_Y(1, i);
		dH(3, i) = inverseJacobian(1, 0) * dH_Y(0, i) + inverseJacobian(1, 1) * dH_Y(1, i);
	}
	
	// we create the Kinematic matrix between the Curvature and the Displacement derivatives
	Matrix aMatrix(3, 4);
	aMatrix(0, 0) = 1.0;
	aMatrix(1, 3) = 1.0;
	aMatrix(2, 1) = 1.0;
	aMatrix(2, 2) = 1.0;

	// we calculate B = strain-displacement relationship matrix , B = (A*dH)
	BSD.addMatrixProduct(0.0, aMatrix, dH, 1.0);

} // end plateFieldInterpolation

void 
MEFI_3D::setPlateTangent()
{
	//compute material tangent 
	plateTangentNDM.Zero();
	plateTangentNDM(0, 0) = 1.0;
	plateTangentNDM(0, 1) = Nu;
	plateTangentNDM(1, 0) = plateTangentNDM(0, 1);
	plateTangentNDM(1, 1) = plateTangentNDM(0, 0);
	plateTangentNDM(2, 2) = (1.0 - Nu) / 2.0;
	plateTangentNDM *= (Eave * Tave * Tave * Tave / 12.0) / (1.0 - Nu * Nu);

	//compute plate tangent 
	KP.Zero();
	for (int i = 0; i < 4; i++) {
		this->plateFieldInterpolation(qdtLocationsB[i][0], qdtLocationsB[i][1]);
		KP.addMatrixTripleProduct(1.0, BSD, plateTangentNDM, qdtWeightsB[i] * detJac);
	}
}

const Matrix&
MEFI_3D::transpose(const Matrix& M)
{
	//we're always transposing 3x12 matrices for this element,
	//so always return a 12x3 .

	static int dim1 = 12;
	static int dim2 = 3;
	static Matrix Mtran(dim1, dim2);

	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim2; j++)
			Mtran(i, j) = M(j, i);
	} // end for i

	return Mtran;
}

void  
MEFI_3D::setTransformationMatrix(void) {

	T.Zero(); // element transformation matrix
	Tt.Zero(); // 3 x 3 - basic transformation matrix
	T6.Zero(); // 6 x 6 - nodal transformation matrix 

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
	Yx = nd4Crds(0) - nd1Crds(0);
	Yy = nd4Crds(1) - nd1Crds(1);
	Yz = nd4Crds(2) - nd1Crds(2);

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