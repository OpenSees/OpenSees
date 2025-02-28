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

// Written by: Amin Pakzad, Pedro Arduino (parduino@uw.edu)
//
// Eight node PML3DVISCOUS element .. a c++ wrapper to fortran routine 
// provided by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2020


#include "PML3DVISCOUS.h"

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
#include <fstream>


// =======================================================================
// PML3DVISCOUS element tcl command
// =======================================================================
void* OPS_PML3DVISCOUS()
{
	// check if the total number of arguments passed is correct
	if (OPS_GetNumRemainingInputArgs() < (9 + PML3DVISCOUS_NUM_PROPS + 3)) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element PML3DVISCOUS eleTag? [8 integer nodeTags] [PML3DVISCOUS_NUM_PARAMS material properties]\n";
		return 0;
	}

	// reading element tag and node numbers 
	int idata[9];
	int num = 9;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer data : could be the tag or the node numbers \n";
		return 0;
	}

	// reading Newmark parameters
	double Newmark[3];
	num = 3;
	if (OPS_GetDoubleInput(&num, Newmark) < 0) {
		opserr << "WARNING: invalid double data: could be Newmark parameters\n";
		return 0;
	}

	// reading material properties
	double dData[PML3DVISCOUS_NUM_PROPS]; num = PML3DVISCOUS_NUM_PROPS;
	if (OPS_GetDoubleInput(&num, dData) < 0) {
		opserr << "WARNING: invalid double data\n";
		return 0;
	}

	// create a new PML3DVISCOUS element and add it to the Domain
	return new PML3DVISCOUS(idata[0], &idata[1], Newmark, dData);
}

// =======================================================================
// static data
// =======================================================================
Matrix  PML3DVISCOUS::tangent(PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::mass(PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Matrix  PML3DVISCOUS::damping(PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
Vector  PML3DVISCOUS::resid(PML3DVISCOUS_NUM_DOF);
double  PML3DVISCOUS::eta = 0.;
double  PML3DVISCOUS::beta = 0.;
double  PML3DVISCOUS::gamma = 0.;
double  PML3DVISCOUS::dt = 0.;
int     PML3DVISCOUS::eleCount = 0;
// int     PML3DVISCOUS::numberOfElements = 0;


// =======================================================================
// null constructor
// =======================================================================
PML3DVISCOUS::PML3DVISCOUS()
	:Element(0, ELE_TAG_PML3DVISCOUS),
	connectedExternalNodes(PML3DVISCOUS_NUM_NODES),
	ubar(PML3DVISCOUS_NUM_DOF),
	ubart(PML3DVISCOUS_NUM_DOF),
    ubarbar(PML3DVISCOUS_NUM_DOF),
    ubarbart(PML3DVISCOUS_NUM_DOF)
{
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		nodePointers[i] = 0;
	}
	dt = 0;
	ubar.Zero();
	ubart.Zero();
    ubarbar.Zero();
    ubarbart.Zero();
	updateflag = 0;
	update_dt = 0;
	eta = 0;
	beta = 0;
	gamma = 0;
	
}

// =======================================================================
// Full constructor
// =======================================================================
PML3DVISCOUS::PML3DVISCOUS(int tag, int* nodeTags, double* nemwarks, double* eleData)
	:Element(tag, ELE_TAG_PML3DVISCOUS),
	connectedExternalNodes(PML3DVISCOUS_NUM_NODES),
	ubar(PML3DVISCOUS_NUM_DOF),
	ubart(PML3DVISCOUS_NUM_DOF),
    ubarbar(PML3DVISCOUS_NUM_DOF),
    ubarbart(PML3DVISCOUS_NUM_DOF)
{
	eleCount++;
	if (eleCount == 1) {
		opserr << "Perfectly Matched Layer 3D (PMLVISCOUS) element -  Written: W. Zhang, E. Taciroglu, A. Pakzad, P. Arduino, UCLA, UCLA, U.Washington, U.Washington\n ";
	}
	// initialize node pointers
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		connectedExternalNodes(i) = nodeTags[i];
		nodePointers[i] = 0;
	}

	// initialize Newmark parameters
	eta   = nemwarks[0];
	beta  = nemwarks[1];
	gamma = nemwarks[2];

	// initialize material properties
	for (int i = 0; i < PML3DVISCOUS_NUM_PROPS; i++)
		props[i] = eleData[i];

	// initialize the ubar and ubart vectors to zero
	ubart.Zero();
	ubar.Zero();
    ubarbar.Zero();
    ubarbart.Zero();
	updateflag = 0;
	update_dt = 0;
}

// =======================================================================
//  destructor
// ======================================================================= 
PML3DVISCOUS::~PML3DVISCOUS()
{

}

// =======================================================================
// Set Domain
// =======================================================================
void  PML3DVISCOUS::setDomain(Domain* theDomain)
{

	Domainptr = theDomain;

	// node pointers
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++)
		nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));

	this->DomainComponent::setDomain(theDomain);

	// create the coordinate vectors
	double coords[PML3DVISCOUS_NUM_DOF];
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& loc = nodePointers[i]->getCrds();
		coords[i * 3] = loc(0);
		coords[i * 3 + 1] = loc(1);
		coords[i * 3 + 2] = loc(2);
	}
	int NDOFEL = PML3DVISCOUS_NUM_DOF;
	int NPROPS = 13;
	int MCRD = 3;
	int NNODE = 8;
	int LFLAGS = 12;
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
		C[i] = 0.0;
		K[i] = 0.0;
		M[i] = 0.0;
		G[i] = 0.0;
		H[i] = 0.0;
	}
	pml3d_(M, C, K, G, H, &NDOFEL, props, coords, &MCRD, &NNODE, &LFLAGS);
	dt = theDomain->getDT();

	// // make C zero 
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
		// C[i] = 0.0;
		// K[i] = 0.0;
		// M[i] = 0.0;
		// G[i] = 0.0;
		// H[i] = 0.0;
	}


	std::ofstream myfile;
	std::string filename;
	int tag = this->getTag();

	// // save M matrix in a file
	// filename = "./PML/M" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << M[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

	// // save C matrix in a file
	// filename = "./PML/C" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << C[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

	// // save K matrix in a file
	// filename = "./PML/K" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << K[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

	// // save G matrix in a file
	// filename = "./PML/G" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << G[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

	// save H matrix in a file
	// filename = "./PML/H" + std::to_string(tag) + ".mat";
	// myfile.open(filename);
	// for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
	// 	for (int j = 0; j < PML3DVISCOUS_NUM_DOF; j++) {
	// 		myfile << H[i*PML3DVISCOUS_NUM_DOF + j] << " ";
	// 	}
	// 	myfile << "\n";
	// }
	// myfile.close();

}

// =======================================================================
// update
// =======================================================================
int PML3DVISCOUS::update(void)
{
	dt = Domainptr->getDT();
	// opserr << "dt = " << dt << "\n";	
	// get u, v, a from nodes and calculate the ubar vector
	int loc = 0;
	double c1 = dt;
	double c2 = dt * dt * 0.5;
	double c3 = dt*dt*dt*((1.0/6.0)-eta);
	double c4 = dt*dt*dt*eta;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getDisp();
		const Vector& vNode = nodePointers[i]->getVel();
		const Vector& aNode = nodePointers[i]->getAccel();
		const Vector& atpdt = nodePointers[i]->getTrialAccel();
		for (int j = 0; j < 9; j++) {
			ubar(loc) = ubart(loc) + uNode(j)*c1 + vNode(j)*c2 + aNode(j)*c3 + atpdt(j)*c4; 
			loc++;
		}
	}

    double keisi = 1.0/48.0;

    loc = 0;
    c1 = dt;
    c2 = dt * dt * 0.5;
    c3 = dt * dt * dt /6.0;
    c4 = dt * dt * dt * dt * (1.0/24.0 - keisi);
    int c5 = dt * dt * dt * dt * keisi;
    for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
        const Vector& uNode = nodePointers[i]->getDisp();
        const Vector& vNode = nodePointers[i]->getVel();
        const Vector& aNode = nodePointers[i]->getAccel();
        const Vector& atpdt = nodePointers[i]->getTrialAccel();
        for (int j = 0; j < 9; j++) {
            ubarbar(loc) = ubarbart(loc) + ubart(loc)*c1 + uNode(j)*c2 + vNode(j)*c3 + aNode(j)*c4 + atpdt(j)*c5; 
            loc++;
        }
    }

	return 0;
}

// =======================================================================
//	return stiffness matrix 
// =======================================================================
const Matrix& PML3DVISCOUS::getTangentStiff()
{
	// check if the dt is changed to update the tangent stiffness matrix
	double cg = eta*dt/beta;
    double keisi = 1.0/48.0;
    double ch = dt * dt * keisi/beta;
	//keff = k + cg*g( k and g are symmetric matrices)
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF*PML3DVISCOUS_NUM_DOF; i++) {
		Keff[i] = K[i] + cg*G[i] + ch*H[i];
	}
	tangent.setData(Keff, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
	return tangent;
}

// =======================================================================
//	return initial stiffness matrix 
// =======================================================================
const Matrix& PML3DVISCOUS::getInitialStiff()
{
	return this->getTangentStiff();
}

// =======================================================================
//	return mass matrix
// =======================================================================
const Matrix& PML3DVISCOUS::getMass()
{
	mass.setData(M, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
	// mass.Zero();
	return mass;
}

// =======================================================================
//	return damping matrix
// =======================================================================
const Matrix& PML3DVISCOUS::getDamp()
{
	damping.setData(C, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
	// damping.Zero();
	return damping;
}

// =======================================================================
// Ressisting force
// =======================================================================
//get residual
const Vector& PML3DVISCOUS::getResistingForce()
{
	// if (innertag==14) {
	// 	opserr << "getResistingForce function is called\n";
	// }
	int numNodalDOF = 9;
	static Vector theVector(PML3DVISCOUS_NUM_DOF);

	// get K into stiff
	tangent.setData(K, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);

	//
	// perform: R = K * u
	//

	int loc = 0;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < numNodalDOF; j++)
			theVector(loc++) = uNode(j);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);
	return resid;
}


// =======================================================================
//
// =======================================================================
//get residual with inertia terms
const Vector&
PML3DVISCOUS::getResistingForceIncInertia()
{
    // R += K*u
	static Vector theVector(PML3DVISCOUS_NUM_DOF);
	tangent.setData(K, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);

	int loc = 0;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < 9; j++)
			theVector(loc++) = uNode(j);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);



	// R += M*a
	loc = 0;
	Node** theNodes = this->getNodePtrs();
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& acc = theNodes[i]->getTrialAccel();
		for (int j = 0; j < 9; j++) {
			theVector(loc++) = acc(j);
		}
	}
	resid.addMatrixVector(1.0, this->getMass(), theVector, 1.0);

	// R += C*v
	loc = 0;
	for (int i = 0; i < PML3DVISCOUS_NUM_NODES; i++) {
		const Vector& vel = theNodes[i]->getTrialVel();
		for (int j = 0; j < 9; j++) {
			theVector(loc++) = vel[j];
		}
	}
	resid.addMatrixVector(1.0, this->getDamp(), theVector, 1.0);


	// R += G*ubar
	tangent.setData(G, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
	resid.addMatrixVector(1.0, tangent, ubar, 1.0);

    // R += H*ubarbar
    tangent.setData(H, PML3DVISCOUS_NUM_DOF, PML3DVISCOUS_NUM_DOF);
    resid.addMatrixVector(1.0, tangent, ubarbar, 1.0);
	
	
	return resid;
}

// =======================================================================
// get the number of external nodes
// =======================================================================
int  PML3DVISCOUS::getNumExternalNodes() const
{
	return PML3DVISCOUS_NUM_NODES;
}

// =======================================================================
// return connected external nodes
// =======================================================================
const ID& PML3DVISCOUS::getExternalNodes()
{
	return connectedExternalNodes;
}

// =======================================================================
// return node pointers
// =======================================================================
Node** PML3DVISCOUS::getNodePtrs(void)
{
	return nodePointers;
}

// =======================================================================
// return number of dofs
// =======================================================================
int  PML3DVISCOUS::getNumDOF()
{
	return PML3DVISCOUS_NUM_DOF;
}

// =======================================================================
// commit state
// =======================================================================
int  PML3DVISCOUS::commitState()
{
	int success = 0;
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PML3DVISCOUS::commitState () - failed in base class";
	}

	// set ubart to ubar
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
		ubart(i) = ubar(i);
	}

	updateflag = 0;
	return success;
}

// =======================================================================
// revert to last commit 
// =======================================================================
int  PML3DVISCOUS::revertToLastCommit()
{
	int success = 0;

	// set ubar to ubart
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
		ubar(i) = ubart(i);
        ubarbar(i) = ubarbart(i);
	}

	return success;
}

// =======================================================================
// revert to start
// =======================================================================
int  PML3DVISCOUS::revertToStart()
{
	int success = 0;

	// set ubar and ubart to zero
	for (int i = 0; i < PML3DVISCOUS_NUM_DOF; i++) {
		ubar(i) = 0.0;
		ubart(i) = 0.0;
        ubarbar(i) = 0.0;
        ubarbart(i) = 0.0;
	}

	return success;
}

// =======================================================================
// add load
// =======================================================================
int PML3DVISCOUS::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return -1;
}

// =======================================================================
// add zero load
// =======================================================================
void  PML3DVISCOUS::zeroLoad()
{
	return;
}

// =======================================================================
// senself
// =======================================================================
int  PML3DVISCOUS::sendSelf(int commitTag,
	Channel& theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// PML3DVISCOUS packs its data into a Vector and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
	static Vector data(PML3DVISCOUS_NUM_PROPS + 4);
	data(0) = this->getTag();

	for (int ii = 1; ii <= PML3DVISCOUS_NUM_PROPS; ii++) {
		data(ii) = props[ii - 1];
	}
	data(PML3DVISCOUS_NUM_PROPS+1) = eta;
	data(PML3DVISCOUS_NUM_PROPS+2) = beta;
	data(PML3DVISCOUS_NUM_PROPS+3) = gamma;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML3DVISCOUS::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}


	// PML3DVISCOUS then sends the tags of its four nodes
	res += theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML3DVISCOUS::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	return res;
}

// =======================================================================
// recvself
// =======================================================================
int  PML3DVISCOUS::recvSelf(int commitTag,
	Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// PML3DVISCOUS creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(PML3DVISCOUS_NUM_PROPS + 4);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML3DVISCOUS::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));


	for (int ii = 1; ii <= PML3DVISCOUS_NUM_PROPS; ii++) {
		props[ii - 1] = data(ii);
	}

	eta   = data(PML3DVISCOUS_NUM_PROPS+1);
	beta  = data(PML3DVISCOUS_NUM_PROPS+2);
	gamma = data(PML3DVISCOUS_NUM_PROPS+3);

	// PML3DVISCOUS now receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML3DVISCOUS::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	return res;
}


// =======================================================================
// display
// =======================================================================
int PML3DVISCOUS::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
	// Get the end point display coords
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	static Vector v5(3);
	static Vector v6(3);
	static Vector v7(3);
	static Vector v8(3);
	nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
	nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
	nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
	nodePointers[3]->getDisplayCrds(v4, fact, displayMode);
	nodePointers[4]->getDisplayCrds(v5, fact, displayMode);
	nodePointers[5]->getDisplayCrds(v6, fact, displayMode);
	nodePointers[6]->getDisplayCrds(v7, fact, displayMode);
	nodePointers[7]->getDisplayCrds(v8, fact, displayMode);

	// place values in coords matrix
	static Matrix coords(8, 3);
	for (int i = 0; i < 3; i++) {
		coords(0, i) = v1(i);
		coords(1, i) = v2(i);
		coords(2, i) = v3(i);
		coords(3, i) = v4(i);
		coords(4, i) = v5(i);
		coords(5, i) = v6(i);
		coords(6, i) = v7(i);
		coords(7, i) = v8(i);
	}

	// fill RGB vector
	static Vector values(8);
	for (int i = 0; i < 8; i++)
		values(i) = 1.0;

	// draw the cube
	return theViewer.drawCube(coords, values, this->getTag());
}

// =======================================================================
// setresponse
// =======================================================================
Response* PML3DVISCOUS::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;

	// char outputData[32];

	// output.tag("ElementOutput");
	// output.attr("eleType", "PML3DVISCOUS");
	// output.attr("eleTag", this->getTag());
	// for (int i = 1; i <= 8; i++) {
	// 	sprintf(outputData, "node%d", i);
	// 	output.attr(outputData, nodePointers[i - 1]->getTag());
	// }

	// if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0) {

	// 	for (int i = 1; i <= 8; i++) {
	// 		sprintf(outputData, "P1_%d", i);
	// 		output.tag("ResponseType", outputData);
	// 		sprintf(outputData, "P2_%d", i);
	// 		output.tag("ResponseType", outputData);
	// 		sprintf(outputData, "P3_%d", i);
	// 		output.tag("ResponseType", outputData);
	// 	}

	// 	theResponse = new ElementResponse(this, 1, resid);
	// }
	// output.endTag(); // ElementOutput
	return theResponse;
}

// =======================================================================
// getresponse
// =======================================================================
int PML3DVISCOUS::getResponse(int responseID, Information& eleInfo)
{
	// static Vector stresses(48);

	// if (responseID == 1)
	// 	return eleInfo.setVector(this->getResistingForce());

	return -1;
}

// =======================================================================
// set parameter
// =======================================================================
int PML3DVISCOUS::setParameter(const char** argv, int argc, Parameter& param)
{
	int res = -1;
	return res;
}

// =======================================================================
// update parameter
// =======================================================================
int PML3DVISCOUS::updateParameter(int parameterID, Information& info)
{
	int res = -1;
	return res;
}

// =======================================================================
// print
// =======================================================================
void  PML3DVISCOUS::Print(OPS_Stream &s, int flag) {
	
  	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "Element: " << this->getTag() << endln;
		s << "type: PML3DVISCOUS \n";
		s << "Nodes: " << connectedExternalNodes;
		s << "eta: " << eta << " beta: " << beta << " gamma: " << gamma << endln;
		s << endln;
		s << "Resisting Force (no inertia): " << this->getResistingForce();
  	}
    
 	 if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"PML3DVISCOUS\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
		for (int i = 1; i < 7; i++)
		s << connectedExternalNodes(i) << ", ";
		s << connectedExternalNodes(7) << "], ";
  	}
	return;
}
