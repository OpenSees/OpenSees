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

// Written by: Long Chen, Pedro Arduino (parduino@uw.edu), Wenyang Zhang and fmk
//
// Eight node PML3D element .. a c++ wrapper to fortran routine 
// providewd by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2020


#include "PML3D.h"

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <OPS_Globals.h>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Parameter.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

void* OPS_PML3D()
{
	if (OPS_GetNumRemainingInputArgs() < (9 + PML3D_NUM_PROPS)) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element PML3D eleTag? [8 integer nodeTags] [PML3D_NUM_PARAMS material properties]\n";
		return 0;
	}

	int idata[9];
	int num = 9;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer data\n";
		return 0;
	}

	double dData[PML3D_NUM_PROPS]; num = PML3D_NUM_PROPS;
	if (OPS_GetDoubleInput(&num, dData) < 0) {
		opserr << "WARNING: invalid double data\n";
		return 0;
	}

	return new PML3D(idata[0], &idata[1], dData);
}

//static data
Matrix  PML3D::tangent(PML3D_NUM_DOF, PML3D_NUM_DOF);
Vector  PML3D::resid(PML3D_NUM_DOF);
Matrix  PML3D::mass(PML3D_NUM_DOF, PML3D_NUM_DOF);
Matrix  PML3D::damping(PML3D_NUM_DOF, PML3D_NUM_DOF);


//null constructor
PML3D::PML3D()
	:Element(0, ELE_TAG_PML3D),
	connectedExternalNodes(PML3D_NUM_NODES)
{
	for (int i = 0; i < PML3D_NUM_NODES; i++) {
		nodePointers[i] = 0;
	}
}


//*********************************************************************
//full constructor
PML3D::PML3D(int tag,
	int* nodeTags,
	double* eleData)
	:Element(tag, ELE_TAG_PML3D),
	connectedExternalNodes(PML3D_NUM_NODES)
{
	for (int i = 0; i < PML3D_NUM_NODES; i++) {
		connectedExternalNodes(i) = nodeTags[i];
		nodePointers[i] = 0;
	}
	for (int i = 0; i < PML3D_NUM_PROPS; i++)
		props[i] = eleData[i];
	for (int i = 0; i < PML3D_NUM_PROPS; i++)
		opserr << props[i] << "\n";
}

//destructor 
PML3D::~PML3D()
{

}


//set domain
void  PML3D::setDomain(Domain* theDomain)
{

	int i;

	//node pointers
	for (i = 0; i < PML3D_NUM_NODES; i++)
		nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

	this->DomainComponent::setDomain(theDomain);

}


//get the number of external nodes
int  PML3D::getNumExternalNodes() const
{
	return PML3D_NUM_NODES;
}


//return connected external nodes
const ID& PML3D::getExternalNodes()
{
	return connectedExternalNodes;
}

Node**
PML3D::getNodePtrs(void)
{
	return nodePointers;
}

//return number of dofs
int  PML3D::getNumDOF()
{
	return PML3D_NUM_DOF;
}


//commit state
int  PML3D::commitState()
{
	// call element commitState to do any base class stuff
	int success = 0;
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PML3D::commitState () - failed in base class";
	}

	return success;
}



//revert to last commit 
int  PML3D::revertToLastCommit()
{
	int success = 0;
	return success;
}


//revert to start 
int  PML3D::revertToStart()
{
	int success = 0;
	return success;
}

//print out element data
void  PML3D::Print(OPS_Stream& s, int flag)
{
	if (flag == 2) {

		s << "#PML3D\n";

		int i;
		const int numNodes = 8;

		for (i = 0; i < numNodes; i++) {
			const Vector& nodeCrd = nodePointers[i]->getCrds();
			const Vector& nodeDisp = nodePointers[i]->getDisp();
			s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << nodeCrd(2)
				<< " " << nodeDisp(0) << " " << nodeDisp(1) << " " << nodeDisp(2) << endln;
		}
	}
	if (flag == OPS_PRINT_CURRENTSTATE) {

		s << "Standard Eight Node PML3D \n";
		s << "Element Number: " << this->getTag() << endln;
		s << "Nodes: " << connectedExternalNodes;

		s << endln;
		s << this->getTag() << " " << connectedExternalNodes(0)
			<< " " << connectedExternalNodes(1)
			<< " " << connectedExternalNodes(2)
			<< " " << connectedExternalNodes(3)
			<< " " << connectedExternalNodes(4)
			<< " " << connectedExternalNodes(5)
			<< " " << connectedExternalNodes(6)
			<< " " << connectedExternalNodes(7)
			<< endln;

		s << "Resisting Force (no inertia): " << this->getResistingForce();
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"PML3D\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
		for (int i = 1; i < 6; i++)
			s << connectedExternalNodes(i) << ", ";
		s << connectedExternalNodes(7) << "], ";
	}
}


//return stiffness matrix 
const Matrix& PML3D::getTangentStiff()
{
	tangent.setData(K, PML3D_NUM_DOF, PML3D_NUM_DOF);

	return tangent;
}


const Matrix& PML3D::getInitialStiff()
{
	return this->getTangentStiff();
}


//return mass matrix
const Matrix& PML3D::getMass()
{
	mass.setData(M, PML3D_NUM_DOF, PML3D_NUM_DOF);

	return mass;
}

//return damping matrix
const Matrix& PML3D::getDamp()
{
	damping.setData(C, PML3D_NUM_DOF, PML3D_NUM_DOF);

	return damping;
}

void  PML3D::zeroLoad()
{
	return;
}


int
PML3D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return -1;
}

//get residual
const Vector&
PML3D::getResistingForce()
{
	int numNodalDOF = 18;
	static Vector theVector(PML3D_NUM_DOF);

	// get K into stiff
	tangent.setData(K, PML3D_NUM_DOF, PML3D_NUM_DOF);

	//
	// perform: R = K * u
	//

	int loc = 0;
	for (int i = 0; i < PML3D_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < numNodalDOF; j++)
			theVector(loc++) = uNode(j);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);

	return resid;
}


//get residual with inertia terms
const Vector&
PML3D::getResistingForceIncInertia()
{
	int numNodalDOF = 18;
	static Vector theVector(PML3D_NUM_DOF);
	static Matrix theMatrix(PML3D_NUM_DOF, PML3D_NUM_DOF);

	//
	// perform: R = P(U) - Pext(t)
	//

	this->getResistingForce();

	//
	// perform: R = R - M * a
	//

	int loc = 0;
	Node** theNodes = this->getNodePtrs();
	for (int i = 0; i < PML3D_NUM_NODES; i++) {
		const Vector& acc = theNodes[i]->getTrialAccel();
		for (int j = 0; j < numNodalDOF; j++) {
			theVector(loc++) = acc(j);
		}
	}
	resid.addMatrixVector(1.0, this->getMass(), theVector, 1.0);

	//
	// perform: R = R + (alphaM * M + betaK0 * K0 + betaK * K) * v
	//            = R + D * v
	//

	// determine the vel vector from ele nodes
	loc = 0;
	for (int i = 0; i < PML3D_NUM_NODES; i++) {
		const Vector& vel = theNodes[i]->getTrialVel();
		for (int j = 0; j < numNodalDOF; j++) {
			theVector(loc++) = vel[j];
		}
	}

	// finally the C * v
	resid.addMatrixVector(1.0, this->getDamp(), theVector, 1.0);

	return resid;
}

int
PML3D::update(void)
{
	static const int numberNodes = 8;

	static double coords[PML3D_NUM_DOF];
	static double U[PML3D_NUM_DOF];
	static double V[PML3D_NUM_DOF];
	static double A[PML3D_NUM_DOF];

	for (int i = 0; i < numberNodes; i++) {
		const Vector& loc = nodePointers[i]->getCrds();
		coords[i * 3] = loc(0);
		coords[i * 3 + 1] = loc(1);
		coords[i * 3 + 2] = loc(2);

		const Vector& u = nodePointers[i]->getTrialDisp();
		U[i * 3] = u(0);
		U[i * 3 + 1] = u(1);
		U[i * 3 + 2] = u(2);
		const Vector& v = nodePointers[i]->getTrialVel();
		V[i * 3] = v(0);
		V[i * 3 + 1] = v(1);
		V[i * 3 + 2] = v(2);
		const Vector& a = nodePointers[i]->getTrialAccel();
		A[i * 3] = a(0);
		A[i * 3 + 1] = a(1);
		A[i * 3 + 2] = a(2);
	}

	//  static double residualV[48];

	// double residualV[PML3D_NUM_DOF];

	int NDOFEL = PML3D_NUM_DOF;
	int NPROPS = 12;
	int MCRD = 3;
	int NNODE = 8;

	pml_(M,
		C,
		K,
		&NDOFEL,
		props,
		&NPROPS,
		coords,
		&MCRD,
		&NNODE);


	// opserr << this->getTag() << "\n";
	// opserr << K[0] << "\n";
	// opserr << C[0] << "\n";
	// opserr << M[0] << "\n";

	// opserr << "DONE CALLING PML\n";
	return 0;
}



int  PML3D::sendSelf(int commitTag,
	Channel& theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// PML3D packs its data into a Vector and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
	static Vector data(PML3D_NUM_PROPS + 1);
	data(0) = this->getTag();

	for (int ii = 1; ii <= PML3D_NUM_PROPS; ii++) {
		data(ii) = props[ii - 1];
	}

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}


	// PML3D then sends the tags of its four nodes
	res += theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML3D::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	return res;
}

int  PML3D::recvSelf(int commitTag,
	Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// PML3D creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(PML3D_NUM_PROPS + 1);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML3D::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));

	for (int ii = 1; ii <= PML3D_NUM_PROPS; ii++) {
		props[ii - 1] = data(ii);
	}

	// PML3D now receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML3D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	return res;
}


int
PML3D::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
	const Vector& end1Crd = nodePointers[0]->getCrds();
	const Vector& end2Crd = nodePointers[1]->getCrds();
	const Vector& end3Crd = nodePointers[2]->getCrds();
	const Vector& end4Crd = nodePointers[3]->getCrds();

	const Vector& end5Crd = nodePointers[4]->getCrds();
	const Vector& end6Crd = nodePointers[5]->getCrds();
	const Vector& end7Crd = nodePointers[6]->getCrds();
	const Vector& end8Crd = nodePointers[7]->getCrds();

	static Matrix coords(8, 3);
	static Vector values(8);
	static Vector P(PML3D_NUM_DOF);

	for (int i = 0; i < 8; i++)
		values(i) = 1.0;

	int error = 0;
	int i;

	if (displayMode >= 0) {

		const Vector& end1Disp = nodePointers[0]->getDisp();
		const Vector& end2Disp = nodePointers[1]->getDisp();
		const Vector& end3Disp = nodePointers[2]->getDisp();
		const Vector& end4Disp = nodePointers[3]->getDisp();
		const Vector& end5Disp = nodePointers[4]->getDisp();
		const Vector& end6Disp = nodePointers[5]->getDisp();
		const Vector& end7Disp = nodePointers[6]->getDisp();
		const Vector& end8Disp = nodePointers[7]->getDisp();

		// for each face of the brick we:
		//   1) determine the coordinates of the displaced point
		//   2) determine the value to be drawn, the stress at nearest gauss point in displayMode dirn
		//   3) get the renderer to draw the face

		for (i = 0; i < 3; i++) {
			coords(0, i) = end1Crd(i) + end1Disp(i) * fact;
			coords(1, i) = end2Crd(i) + end2Disp(i) * fact;
			coords(2, i) = end3Crd(i) + end3Disp(i) * fact;
			coords(3, i) = end4Crd(i) + end4Disp(i) * fact;
			coords(4, i) = end5Crd(i) + end5Disp(i) * fact;
			coords(5, i) = end6Crd(i) + end6Disp(i) * fact;
			coords(6, i) = end7Crd(i) + end7Disp(i) * fact;
			coords(7, i) = end8Crd(i) + end8Disp(i) * fact;
		}

		error = theViewer.drawCube(coords, values, this->getTag());

	}
	else {

		int mode = displayMode * -1;

		const Matrix& eigen1 = nodePointers[0]->getEigenvectors();
		const Matrix& eigen2 = nodePointers[1]->getEigenvectors();
		const Matrix& eigen3 = nodePointers[2]->getEigenvectors();
		const Matrix& eigen4 = nodePointers[3]->getEigenvectors();
		const Matrix& eigen5 = nodePointers[4]->getEigenvectors();
		const Matrix& eigen6 = nodePointers[5]->getEigenvectors();
		const Matrix& eigen7 = nodePointers[6]->getEigenvectors();
		const Matrix& eigen8 = nodePointers[7]->getEigenvectors();

		if (eigen1.noCols() >= mode) {

			for (i = 0; i < 3; i++) {
				coords(0, i) = end1Crd(i) + eigen1(i, mode - 1) * fact;
				coords(1, i) = end2Crd(i) + eigen2(i, mode - 1) * fact;
				coords(2, i) = end3Crd(i) + eigen3(i, mode - 1) * fact;
				coords(3, i) = end4Crd(i) + eigen4(i, mode - 1) * fact;
				coords(4, i) = end5Crd(i) + eigen5(i, mode - 1) * fact;
				coords(5, i) = end6Crd(i) + eigen6(i, mode - 1) * fact;
				coords(6, i) = end7Crd(i) + eigen7(i, mode - 1) * fact;
				coords(7, i) = end8Crd(i) + eigen8(i, mode - 1) * fact;
			}

			for (int i = 0; i < 8; i++)
				values(i) = 0.0;

			error = theViewer.drawCube(coords, values, this->getTag());
			opserr << "\n tag: " << this->getTag() << endln;
			opserr << coords;
			opserr << values;

		}
	}

	return error;
}

Response*
PML3D::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;

	char outputData[32];

	output.tag("ElementOutput");
	output.attr("eleType", "PML3D");
	output.attr("eleTag", this->getTag());
	for (int i = 1; i <= 8; i++) {
		sprintf(outputData, "node%d", i);
		output.attr(outputData, nodePointers[i - 1]->getTag());
	}

	if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {

		for (int i = 1; i <= 8; i++) {
			sprintf(outputData, "P1_%d", i);
			output.tag("ResponseType", outputData);
			sprintf(outputData, "P2_%d", i);
			output.tag("ResponseType", outputData);
			sprintf(outputData, "P3_%d", i);
			output.tag("ResponseType", outputData);
		}

		theResponse = new ElementResponse(this, 1, resid);
	}
	output.endTag(); // ElementOutput
	return theResponse;
}

int
PML3D::getResponse(int responseID, Information& eleInfo)
{
	static Vector stresses(48);

	if (responseID == 1)
		return eleInfo.setVector(this->getResistingForce());

	return -1;
}

int
PML3D::setParameter(const char** argv, int argc, Parameter& param)
{
	int res = -1;
	return res;
}

int
PML3D::updateParameter(int parameterID, Information& info)
{
	int res = -1;
	return res;
}

