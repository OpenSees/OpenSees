// Code written/implemented by:	Carlos Lopez Olea (carlos.lopez.o@ug.uchile.cl)
//
// User documentation available at: https://github.com/carloslopezolea/MEFI
//
// Created: 06/2022
//
// Description: The Membrane Fiber (MEFI) element, is described by four nodes, each containing three degrees of freedom (DOFs), two translations, and one in-plane rotation 
// (drilling) DOF, which incorporates a blended interpolation function for the displacements over the element. The element formulation accommodates the quadrature 
// points and weights of the classical finite element formulation of membrane elements to resemble strips (fibers), similarly to macroscopic elements.
//
// Reference:
// 1.- López, C. N., Rojas, F., & Massone, L. M. (2022). Membrane fiber element for reinforced concrete walls – the benefits of macro and micro modeling approaches. Engineering Structures, 254, 113819.
//
// Source: /usr/local/cvs/OpenSees/SRC/element/mefi/MEFI.cpp
//
// Rev: 1.0                                                       

#include <MEFI.h>
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

void* OPS_MEFI() 
{
	//pointer to a element that will be returned   
	Element* theEle = 0;

	//check model dimensions and nodal dofs 
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
	if (ndm != 2 || ndf != 3) {
		opserr << "WARNING element MEFI: invalid model dimensions and/or nodal DOFs" << endln;
		return 0;
    }

	//check number of arguments provided
    if (OPS_GetNumRemainingInputArgs() < 8) {
		opserr << "WARNING element MEFI: not enough args provided, want: element MEFI eleTag iNode jNode kNode lNode numFibers -thick -sec" << endln;
		return 0;
    } 

	//create array to store integer data
	int iData[6];
	
	//check element tag and add to integer data
	int numData = 1;
	if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
		opserr << "WARNING element MEFI: invalid integer element tag" << endln;
		return 0;
	}

	//check element nodes and add to integer data
	numData = 4;
	if (OPS_GetIntInput(&numData, &iData[1]) != 0) {
		opserr << "WARNING element MEFI: invalid integer node tag for element with tag " << iData[0] << endln;
		return 0;
	}

	//check element number of fibers and add to integer data
	numData = 1;
	if (OPS_GetIntInput(&numData, &iData[5]) != 0) {
		opserr << "WARNING element MEFI: invalid integer numFibers for element with tag " << iData[0] << endln;
		return 0;
	}

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
				opserr << "WARNING element MEFI: invalid width value for element with tag " << iData[0] << endln;
				return 0;
			}
		}
		else if (strcmp(str, "-sec") == 0) {
			numData = iData[5];
			if (OPS_GetIntInput(&numData, secTags) != 0) {
				opserr << "WARNING element MEFI: invalid section tag for element with tag " << iData[0] << endln;
				return 0;
			}
			for (int i = 0; i < iData[5]; i++) {
				theSec[i] = 0;
				theSec[i] = OPS_getSectionForceDeformation(secTags[i]);
				if (theSec[i] == 0) {
					opserr << "WARNING element MEFI: invalid section tag " << secTags[i] << " for element with tag " << iData[0] << endln;
					return 0;
				}
			}
		}
		numArgs = OPS_GetNumRemainingInputArgs();

	}//end while

	//now create the element and add it to the Domain
	theEle = new MEFI(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], theSec, theWidth);

	//cleanup dynamic memory
	if (theWidth != 0) { delete[] theWidth; }
	if (secTags != 0)  { delete[] secTags; }
	if (theSec != 0) { delete[] theSec; }

	//check if the element is not null
	if (theEle == 0) {
		opserr << "WARNING element MEFI: ran out of memory creating element with tag " << iData[0] << endln;
		return 0;
	}
	
	else { return theEle; }
	
}//end OPS_MEFI

Matrix MEFI::K(12, 12);
Vector MEFI::P(12);
double MEFI::detJac;
Matrix MEFI::BSD(3, 12);

MEFI::MEFI(int tag, 
	int nd1, int nd2, int nd3, int nd4, 
	int numFibers,
	SectionForceDeformation** sec,
	double* width)
	:Element (tag, ELE_TAG_MEFI), 
	connectedExternalNodes(4), nd1Crds(2), nd2Crds(2), nd3Crds(2), nd4Crds(2), Q(12),
	theSection(0), nip(numFibers), qdtLocations(numFibers, 2), qdtWeights(numFibers), lw(0), Ki(0)
{
	//check width input
	if (width == 0) {
		opserr << "MEFI::MEFI() - null width array passed" << endln;
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
		qdtLocations(i, 0) = x(i) * (2 / lw);
		qdtLocations(i, 1) = 0.001;
		qdtWeights(i) = b(i) * (2 / lw) * 2;
	}

    //allocate arrays of pointers to sections
	theSection = new SectionForceDeformation *[nip];
    if (theSection == 0) {
      opserr << "MEFI::MEFI() - failed allocate section model pointer" << endln;
      exit(-1);
    }

	//get copies of the sections
	for (int i = 0; i < nip; i++) {
		if (sec[i] == 0) {
			opserr << "MEFI::MEFI() - null section pointer passed" << endln;
			exit(-1);
		}
		theSection[i] = sec[i]->getCopy();
		if (theSection[i] == 0) {
			opserr << "MEFI::MEFI() - failed to copy section" << endln;
			exit(-1);
		}
	}

    //set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;

	//set node pointer to null
    for (int i=0; i<4; i++) { theNodes[i] = 0; }

}//end MEFI


MEFI::MEFI()
	
	:Element(0, ELE_TAG_MEFI),
	connectedExternalNodes(4), nd1Crds(2), nd2Crds(2), nd3Crds(2), nd4Crds(2), Q(12),
	theSection(0), nip(0), qdtLocations(0, 2), qdtWeights(0), lw(0), Ki(0)
{

	for (int i = 0; i < 4; i++) { theNodes[i] = 0; }

}//end MEFI


MEFI::~MEFI()
{    
	for (int i = 0; i < nip; i++) { delete theSection[i]; }
	if (theSection != 0) { delete [] theSection; }
	if (Ki != 0) { delete Ki; }

}//end ~MEFI

int 
MEFI::getNumExternalNodes() const
{
    return 4;

}//end getNumExternalNodes

const ID&
MEFI::getExternalNodes()
{
    return connectedExternalNodes;

}//end getExternalNodes

Node **
MEFI::getNodePtrs(void)
{
  return theNodes;

}//end getNodePtrs

int
MEFI::getNumDOF()
{
    return 12;

}//end getNumDOF

void
MEFI::setDomain(Domain *theDomain)
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
		opserr << "MEFI::setDomain(): node not found in domain for element with tag " << this->getTag() << endln;
		return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    int dofNd3 = theNodes[2]->getNumberDOF();
    int dofNd4 = theNodes[3]->getNumberDOF();
    
	//check dofs compatibility
    if (dofNd1 != 3 || dofNd2 != 3 || dofNd3 != 3 || dofNd4 != 3) {
		opserr << "MEFI::setDomain(): 3 dofs required at all nodes for element with tag " << this->getTag() << endln;
		return;
    }

	//get coordinates of end nodes
	nd1Crds = theNodes[0]->getCrds();
	nd2Crds = theNodes[1]->getCrds();
	nd3Crds = theNodes[2]->getCrds();
	nd4Crds = theNodes[3]->getCrds();

	//calculate the element height and perform checks
	double h1 = pow(pow(nd4Crds(0) - nd1Crds(0), 2.0) + pow(nd4Crds(1) - nd1Crds(1), 2.0), 0.5);
	double h2 = pow(pow(nd3Crds(0) - nd2Crds(0), 2.0) + pow(nd3Crds(1) - nd2Crds(1), 2.0), 0.5);

	//check if element height is zero
	if ((h1 == 0.0) || (h2 == 0.0)) {
		opserr << "MEFI::setDomain(): one of the sides is zero for element with tag " << this->getTag() << endln;
		exit(-1);
	}

	//check if element has constant height
	if ((h1 / h2 > 1.01) || (h1 / h2 < 0.99)) {
		opserr << "MEFI::setDomain(): not constant height for element with tag " << this->getTag() << endln;
		exit(-1);
	}

	//calculate the element width and perform checks
	double b1 = pow(pow(nd2Crds(0) - nd1Crds(0), 2.0) + pow(nd2Crds(1) - nd1Crds(1), 2.0), 0.5);
	double b2 = pow(pow(nd3Crds(0) - nd4Crds(0), 2.0) + pow(nd3Crds(1) - nd4Crds(1), 2.0), 0.5);

	//check width of element
	if ((lw / b1 > 1.01) || (lw / b1 < 0.99) || (lw / b2 > 1.01) || (lw / b2 < 0.99)) {
		opserr << "MEFI::setDomain(): nodes coordinates are not matched with fibers width for element with tag " << this->getTag() << endln;
		exit(-1);
	}

    this->DomainComponent::setDomain(theDomain);

}//end setDomain

int
MEFI::commitState()
{
    int retVal = 0;
    //call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "MEFI::commitState(): failed in base class for element with tag " << this->getTag() << endln;
    }    
    //loop over the integration points and commit the material states
    for (int i = 0; i < nip; i++) {
      retVal += theSection[i]->commitState();
	}
    return retVal;

}//end commitState

int 
MEFI::revertToLastCommit()
{
    int retVal = 0;
    //loop over the integration points and revert to last committed state
    for (int i = 0; i < nip; i++) {
		retVal += theSection[i]->revertToLastCommit();
	}
    return retVal;

}//end revertToLastCommit

int 
MEFI::revertToStart()
{
    int retVal = 0;
    //loop over the integration points and revert states to start
    for (int i = 0; i < nip; i++) {
		retVal += theSection[i]->revertToStart();
	}
    return retVal;

}//end revertToStart

int
MEFI::update()
{
	const Vector &disp1 = theNodes[0]->getTrialDisp();
	const Vector &disp2 = theNodes[1]->getTrialDisp();
	const Vector &disp3 = theNodes[2]->getTrialDisp();
	const Vector &disp4 = theNodes[3]->getTrialDisp();

	//in-plane displacements for membrane behavior
	P.Zero();
	for (int i = 0; i < 3; i++) {
		P(i) = disp1(i);
		P(i + 3) = disp2(i);
		P(i + 6) = disp3(i);
		P(i + 9) = disp4(i);
	}
	//strains at each fiber
	Vector strainAtQuadraturePoint(3);
	strainAtQuadraturePoint.Zero();

	int ret = 0;

	//loop over the integration points
	for (int i = 0; i < nip; i++) {
		this->membraneFieldInterpolation(qdtLocations(i,0), qdtLocations(i, 1));
		strainAtQuadraturePoint = BSD * P;
		ret += theSection[i]->setTrialSectionDeformation(strainAtQuadraturePoint);
	}
	return ret;

}//end update

const Matrix&
MEFI::getTangentStiff()
{
	K.Zero();
	//loop over the integration points
	for (int i = 0; i < nip; i++) {
		const Matrix& D = theSection[i]->getSectionTangent();
		this->membraneFieldInterpolation(qdtLocations(i, 0), qdtLocations(i, 1));
		K.addMatrixTripleProduct(1.0, BSD, D, qdtWeights(i) * detJac);
	}
	
	return K;

}//end getTangentStiff

const Matrix&
MEFI::getInitialStiff()
{
	if (Ki != 0) { return *Ki; }
	K.Zero();
	// Loop over the integration points
	for (int i = 0; i < nip; i++) {
		const Matrix& D = theSection[i]->getInitialTangent();
		this->membraneFieldInterpolation(qdtLocations(i, 0), qdtLocations(i, 1));
		K.addMatrixTripleProduct(1.0, BSD, D, qdtWeights(i) * detJac);
	}

	Ki = new Matrix(K);
	return K;

}//end getInitialStiff

const Matrix&
MEFI::getMass()
{
	K.Zero();
	double rhoi = 0;
	double massElement = 0;
	//we calculate the mass of the Element 
	for (int i = 0; i < nip; i++) {
		rhoi = theSection[i]->getRho();
		this->membraneFieldInterpolation(qdtLocations(i, 0), qdtLocations(i, 1));
		massElement += rhoi * qdtWeights(i) * detJac;
	}

	//we assume lump mass for each node, only in the 1,2 direction
	for (int i = 0; i < 4; i++) {
		K(3 * i, 3 * i) = 0.25 * massElement;
		K(3 * i + 1, 3 * i + 1) = 0.25 * massElement;
	}

	return K;

}//end getMass

//N/A to this model - no element loads
void 
MEFI::zeroLoad(void)
{
  	return;

}//end zeroLoad

// N/A to this model - no element loads
int 
MEFI::addLoad(ElementalLoad *theLoad, double loadFactor)
{
	return 0;

}//end addLoad

int 
MEFI::addInertiaLoadToUnbalance(const Vector &accel)
{
	Q.Zero();
	K.Zero();
	//compute mass matrix
	K = this->getMass();

	//want to add ( - fact * M R * accel ) to unbalance
	for (int i = 0; i < 4; i++) {
		const Vector& Raccel = theNodes[i]->getRV(accel);
		Q(3 * i + 0) += -K(3 * i + 0, 3 * i + 0) * Raccel(0);
		Q(3 * i + 1) += -K(3 * i + 1, 3 * i + 1) * Raccel(1);
	} 

	return 0;

}//end addInertiaLoadToUnbalance

const Vector&
MEFI::getResistingForce()
{
	P.Zero();
	//loop over the integration points
	for (int i = 0; i < nip; i++) {
		const Vector& Stress = theSection[i]->getStressResultant();
		this->membraneFieldInterpolation(qdtLocations(i, 0), qdtLocations(i, 1));
		P += transpose(BSD) * Stress * qdtWeights(i) * detJac;	
	}

	// Subtract other external nodal loads ... P_res = P_int - P_ext
	//P = P - Q;
	P.addVector(1.0, Q, -1.0);

	return P;

}//end getResistingForce

const Vector&
MEFI::getResistingForceIncInertia()
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
		P(3 * i + 0) += K(3 * i + 0, 3 * i + 0) * accel(0);
		P(3 * i + 1) += K(3 * i + 1, 3 * i + 1) * accel(1);
	}
	
	return P;

}//end getResistingForceIncInertia

int 
MEFI::sendSelf(int commitTag, Channel& theChannel)
{
	int res = 0;
	int dataTag = this->getDbTag();
	//MEFI packs its data into a Vector and sends this to theChannel
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
		opserr << "WARNING MEFI::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}
	//MEFI sends the ids of its sections and nodes
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
		opserr << "WARNING MEFI::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}
	//MEFI asks its material objects to send themselves
	for (int i = 0; i < nip; i++) {
		res += theSection[i]->sendSelf(commitTag, theChannel);
		if (res < 0) {
			opserr << "WARNING MEFI::sendSelf() - " << this->getTag() << " failed to send its Material\n";
			return res;
		}
	}
	return res;

}//end sendSelf

int 
MEFI::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	int dataTag = this->getDbTag();
	//MEFI creates a Vector, receives the Vector and then sets the 
	//internal data with the data in the Vector
	static Vector data(6);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING MEFI::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	nip = (int)data(1);
	alphaM = data(2);
	betaK = data(3);
	betaK0 = data(4);
	betaKc = data(5);

	static ID idData(2 * nip + 4);
	//MEFI receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, idData);
	if (res < 0) {
		opserr << "WARNING MEFI::recvSelf() - " << this->getTag() << " failed to receive ID\n";
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
			opserr << "MEFI::recvSelf() - Could not allocate Section array\n";
			return -1;
		}
		for (int i = 0; i < nip; i++) {
			int matClassTag = idData(i);
			int matDbTag = idData(i + nip);
			//allocate new material with the sent class tag
			theSection[i] = theBroker.getNewSection(matClassTag);
			if (theSection[i] == 0) {
				opserr << "MEFI::recvSelf() - Broker could not create Section of class type " << matClassTag << endln;
				return -1;
			}
			//now receive materials into the newly allocated space
			theSection[i]->setDbTag(matDbTag);
			res += theSection[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "MEFI::recvSelf() - material " << i << "failed to recv itself\n";
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
					opserr << "MEFI::recvSelf() - material " << i << "failed to create\n";
					return -1;
				}
			}
			//receive the material
			theSection[i]->setDbTag(matDbTag);
			res += theSection[i]->recvSelf(commitTag, theChannel, theBroker);
			if (res < 0) {
				opserr << "MEFI::recvSelf() - material " << i << "failed to recv itself\n";
				return res;
			}
		}
	}
	
	return res;

}//end recvSelf

void 
MEFI::Print(OPS_Stream &s, int flag)
{
	if (flag == 0) {
		s << "MEFI Element tag: " << this->getTag() << endln;
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

}//end Print

int 
MEFI::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
	//get the end point display coords
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	theNodes[0]->getDisplayCrds(v1, fact, displayMode);
	theNodes[1]->getDisplayCrds(v2, fact, displayMode);
	theNodes[2]->getDisplayCrds(v3, fact, displayMode);
	theNodes[3]->getDisplayCrds(v4, fact, displayMode);

	//place values in coords matrix
	static Matrix coords(4, 3);
	for (int i = 0; i < 3; i++) {
		coords(0, i) = v1(i);
		coords(1, i) = v2(i);
		coords(2, i) = v3(i);
		coords(3, i) = v4(i);
	}

	static Vector values(4);
	for (int i = 0; i < 4; i++) {
		values(i) = 0.0;
	}

	// draw the polygon
	return theViewer.drawPolygon(coords, values, this->getTag());

}//end displaySelf

Response* MEFI::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "MEFI");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes(0));
  output.attr("node2", connectedExternalNodes(1));
  output.attr("node3", connectedExternalNodes(2));
  output.attr("node4", connectedExternalNodes(3));

  char dataOut[10];

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
		  sprintf(dataOut, "M3_%d", i);
		  output.tag("ResponseType", dataOut);
	  }

	  theResponse = new ElementResponse(this, 1, Vector(12));
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
	  if (argc != 3) {
		  opserr << "WARNING: Number of recorder input for section is: " << argc - 1 << "; should be 2: secTag and $Response_Type.\n";
		  return 0;
	  }

	  int secNum = atoi(argv[1]);

	  output.tag("Material");
	  output.attr("number", secNum);

	  theResponse = theSection[secNum - 1]->setResponse(&argv[argc - 1], argc - 2, output);

  }

  output.endTag(); //elementOutput

  return theResponse;

}//end setResponse

int  MEFI::getResponse(int responseID, Information &eleInfo)
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

}//end getResponse

void 
MEFI::membraneFieldInterpolation(double xi, double eta)
{
	//get nodal coordinates
	double x1 = nd1Crds(0); double y1 = nd1Crds(1);
	double x2 = nd2Crds(0); double y2 = nd2Crds(1);
	double x3 = nd3Crds(0); double y3 = nd3Crds(1); 
	double x4 = nd4Crds(0); double y4 = nd4Crds(1);

	//calculate the shape function derivatives matrix w.r.t. the natural coord system
	Matrix dShapeFunction(2, 8); 
	dShapeFunction.Zero();
	dShapeFunction(0, 0) = 0.25 * (1.0 - eta) * (2.0 * xi + eta);
	dShapeFunction(0, 1) = 0.25 * (1.0 - eta) * (2.0 * xi - eta);
	dShapeFunction(0, 2) = 0.25 * (1.0 + eta) * (2.0 * xi + eta);
	dShapeFunction(0, 3) = 0.25 * (1.0 + eta) * (2.0 * xi - eta);
	dShapeFunction(0, 4) = 0.5 * (-2.0 * xi) * (1.0 - eta);
	dShapeFunction(0, 5) = 0.5 * (1.0 - pow(eta, 2));
	dShapeFunction(0, 6) = 0.5 * (-2.0 * xi) * (1.0 + eta);
	dShapeFunction(0, 7) = 0.5 * (-1.0) * (1.0 - pow(eta, 2));
	dShapeFunction(1, 0) = 0.25 * (1.0 - xi) * (2.0 * eta + xi);
	dShapeFunction(1, 1) = 0.25 * (1.0 + xi) * (2.0 * eta - xi);
	dShapeFunction(1, 2) = 0.25 * (1.0 + xi) * (2.0 * eta + xi);
	dShapeFunction(1, 3) = 0.25 * (1.0 - xi) * (2.0 * eta - xi);
	dShapeFunction(1, 4) = 0.5 * (-1.0) * (1.0 - pow(xi, 2));
	dShapeFunction(1, 5) = 0.5 * (-2.0 * eta) * (1.0 + xi);
	dShapeFunction(1, 6) = 0.5 * (1.0 - pow(xi, 2));
	dShapeFunction(1, 7) = 0.5 * (-2.0 * eta) * (1.0 - xi);

	//calculate the intermedial node coord for a 8 node serendipity field interpolation
	Matrix localCoord8Nodes(8, 2); 
	localCoord8Nodes.Zero();
	localCoord8Nodes(0, 0) = x1;              localCoord8Nodes(0, 1) = y1;
	localCoord8Nodes(1, 0) = x2;              localCoord8Nodes(1, 1) = y2;
	localCoord8Nodes(2, 0) = x3;              localCoord8Nodes(2, 1) = y3;
	localCoord8Nodes(3, 0) = x4;              localCoord8Nodes(3, 1) = y4;
	localCoord8Nodes(4, 0) = 0.5 * (x1 + x2); localCoord8Nodes(4, 1) = 0.5 * (y1 + y2);
	localCoord8Nodes(5, 0) = 0.5 * (x2 + x3); localCoord8Nodes(5, 1) = 0.5 * (y2 + y3);
	localCoord8Nodes(6, 0) = 0.5 * (x3 + x4); localCoord8Nodes(6, 1) = 0.5 * (y3 + y4);
	localCoord8Nodes(7, 0) = 0.5 * (x4 + x1); localCoord8Nodes(7, 1) = 0.5 * (y4 + y1);

	//calculate jacobian matrix
	Matrix jacMatrix(2, 2); 
	jacMatrix.Zero();
	jacMatrix = dShapeFunction * localCoord8Nodes;

	//calculate jacobian determinant
	detJac = jacMatrix(0, 0) * jacMatrix(1, 1) - jacMatrix(0, 1) * jacMatrix(1, 0);

	//calculate jacobian matrix inverse
	Matrix invJacMatrix(2, 2); 
	invJacMatrix.Zero();
	invJacMatrix(0, 0) = jacMatrix(1, 1) / detJac;
	invJacMatrix(1, 0) = -jacMatrix(0, 1) / detJac;
	invJacMatrix(0, 1) = -jacMatrix(1, 0) / detJac;
	invJacMatrix(1, 1) = jacMatrix(0, 0) / detJac;

	//calculate jMatrix = [inverseJacobian zeros(2,2) ; zeros(2,2) inverseJacobian]
	Matrix jMatrix(4, 4); 
	jMatrix.Zero();
	jMatrix(0, 0) = invJacMatrix(0, 0);
	jMatrix(0, 1) = invJacMatrix(0, 1);
	jMatrix(1, 0) = invJacMatrix(1, 0);
	jMatrix(1, 1) = invJacMatrix(1, 1);
	jMatrix(2, 2) = invJacMatrix(0, 0);
	jMatrix(2, 3) = invJacMatrix(0, 1);
	jMatrix(3, 2) = invJacMatrix(1, 0);
	jMatrix(3, 3) = invJacMatrix(1, 1);

	//calculate cMatrix = shape function derivatives w.r.t. the natural coord system (xi and eta) for the field interpolation
	Matrix cMatrix(4, 16); 
	cMatrix.Zero();
	cMatrix(0, 0) = (0.5) * (-(0.5) + (0.75) * (eta)-0.25 * (pow(eta, 3)));
	cMatrix(0, 1) = 0.0;
	cMatrix(0, 2) = (0.5) * (0.25 - 0.25 * eta - 0.25 * (pow(eta, 2)) + 0.25 * (pow(eta, 3)));
	cMatrix(0, 3) = 0.0;
	cMatrix(0, 4) = (0.5) * ((0.5) - (0.75) * (eta)+0.25 * (pow(eta, 3)));
	cMatrix(0, 5) = 0.0;
	cMatrix(0, 6) = (0.5) * (-(0.25) + 0.25 * eta + 0.25 * (pow(eta, 2)) - 0.25 * (pow(eta, 3)));
	cMatrix(0, 7) = 0.0;
	cMatrix(0, 8) = (0.5) * ((0.5) + (0.75) * (eta)-0.25 * (pow(eta, 3)));
	cMatrix(0, 9) = 0.0;
	cMatrix(0, 10) = (0.5) * (0.25 + 0.25 * eta - 0.25 * (pow(eta, 2)) - 0.25 * (pow(eta, 3)));
	cMatrix(0, 11) = 0.0;
	cMatrix(0, 12) = (0.5) * (-(0.5) - (0.75) * (eta)+0.25 * (pow(eta, 3)));
	cMatrix(0, 13) = 0.0;
	cMatrix(0, 14) = (0.5) * (-(0.25) - 0.25 * eta + 0.25 * (pow(eta, 2)) + 0.25 * (pow(eta, 3)));
	cMatrix(0, 15) = 0.0;

	cMatrix(1, 0) = (0.5) * (-(0.75) + (0.75) * (pow(eta, 2))) * (1.0 - xi);
	cMatrix(1, 1) = 0.0;
	cMatrix(1, 2) = (0.5) * (-(0.25) - (0.5) * eta + (0.75) * (pow(eta, 2))) * (-1.0 + xi);
	cMatrix(1, 3) = 0.0;
	cMatrix(1, 4) = (0.5) * (-(0.75) + (0.75) * (pow(eta, 2))) * (1.0 + xi);
	cMatrix(1, 5) = 0.0;
	cMatrix(1, 6) = (0.5) * (-(0.25) - (0.5) * eta + (0.75) * (pow(eta, 2))) * (-1.0 - xi);
	cMatrix(1, 7) = 0.0;
	cMatrix(1, 8) = (0.5) * ((0.75) - (0.75) * (pow(eta, 2))) * (1 + xi);
	cMatrix(1, 9) = 0.0;
	cMatrix(1, 10) = (0.5) * (-(0.25) + (0.5) * eta + (0.75) * (pow(eta, 2))) * (-1 - xi);
	cMatrix(1, 11) = 0.0;
	cMatrix(1, 12) = (0.5) * ((0.75) - (0.75) * (pow(eta, 2))) * (1 - xi);
	cMatrix(1, 13) = 0.0;
	cMatrix(1, 14) = (0.5) * (-(0.25) + (0.5) * eta + (0.75) * (pow(eta, 2))) * (-1 + xi);
	cMatrix(1, 15) = 0.0;

	cMatrix(2, 0) = 0.0;
	cMatrix(2, 1) = (0.5) * (1 - eta) * (-(0.75) + (0.75) * (pow(xi, 2)));
	cMatrix(2, 2) = 0.0;
	cMatrix(2, 3) = (0.5) * (1 - eta) * (-(0.25) - 0.5 * xi + (0.75) * (pow(xi, 2)));
	cMatrix(2, 4) = 0.0;
	cMatrix(2, 5) = (0.5) * (1 - eta) * ((0.75) - (0.75) * (pow(xi, 2)));
	cMatrix(2, 6) = 0.0;
	cMatrix(2, 7) = (0.5) * (1 - eta) * (-(0.25) + 0.5 * xi + (0.75) * (pow(xi, 2)));
	cMatrix(2, 8) = 0.0;
	cMatrix(2, 9) = (0.5) * (1 + eta) * ((0.75) - (0.75) * (pow(xi, 2)));
	cMatrix(2, 10) = 0.0;
	cMatrix(2, 11) = (0.5) * (1 + eta) * (-(0.25) + 0.5 * xi + (0.75) * (pow(xi, 2)));
	cMatrix(2, 12) = 0.0;
	cMatrix(2, 13) = (0.5) * (1 + eta) * (-(0.75) + (0.75) * (pow(xi, 2)));
	cMatrix(2, 14) = 0.0;
	cMatrix(2, 15) = (0.5) * (1 + eta) * (-0.25 - 0.5 * xi + (0.75) * (pow(xi, 2)));

	cMatrix(3, 0) = 0.0;
	cMatrix(3, 1) = (0.5) * (-0.5 + (0.75) * (xi)-(0.25) * (pow(xi, 3)));
	cMatrix(3, 2) = 0.0;
	cMatrix(3, 3) = (0.5) * (-(0.25) + (0.25) * xi + (0.25) * (pow(xi, 2)) - (0.25) * (pow(xi, 3)));
	cMatrix(3, 4) = 0.0;
	cMatrix(3, 5) = (0.5) * (-(0.5) - (0.75) * (xi)+(0.25) * (pow(xi, 3)));
	cMatrix(3, 6) = 0.0;
	cMatrix(3, 7) = (0.5) * (0.25 + 0.25 * xi - (0.25) * (pow(xi, 2)) - (0.25) * (pow(xi, 3)));
	cMatrix(3, 8) = 0.0;
	cMatrix(3, 9) = (0.5) * ((0.5) + (0.75) * (xi)-(0.25) * (pow(xi, 3)));
	cMatrix(3, 10) = 0.0;
	cMatrix(3, 11) = (0.5) * (-(0.25) - (0.25) * xi + (0.25) * (pow(xi, 2)) + (0.25) * (pow(xi, 3)));
	cMatrix(3, 12) = 0.0;
	cMatrix(3, 13) = (0.5) * ((0.5) - (0.75) * (xi)+(0.25) * (pow(xi, 3)));
	cMatrix(3, 14) = 0.0;
	cMatrix(3, 15) = (0.5) * ((0.25) - (0.25) * xi - (0.25) * (pow(xi, 2)) + (0.25) * (pow(xi, 3)));
	
	//calculate the transformation matrix from a 4 node blended interpolation with 4DOf per node to a 4 node with rotational DOF
	double x21 = 0.5 * (x2 - x1);
	double x34 = 0.5 * (x3 - x4);
	double y41 = 0.5 * (y4 - y1);
	double y32 = 0.5 * (y3 - y2);

	Matrix trMatrix(16, 12); 
	trMatrix.Zero();
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

	//create the kinematic matrix between the strain and the displacement derivatives
	Matrix aMatrix(3, 4); 
	aMatrix.Zero();
	aMatrix(0, 0) = 1.0;
	aMatrix(1, 3) = 1.0;
	aMatrix(2, 1) = 1.0;
	aMatrix(2, 2) = 1.0;

	//calculate B = strain-displacement relationship matrix , B = (A*J*C*T)
	BSD.Zero();
	BSD = aMatrix * jMatrix * cMatrix * trMatrix;

} //end membraneFieldInterpolation

const Matrix&
MEFI::transpose(const Matrix& M)
{
	static Matrix Mtran(12, 3);

	for (int i = 0; i < 12; i++) {
		for (int j = 0; j < 3; j++)
			Mtran(i, j) = M(j, i);
	}

	return Mtran;

}//end transpose

