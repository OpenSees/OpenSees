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
// Eight node PML2DVISCOUS element .. a c++ wrapper to fortran routine 
// provided by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2023


#include "PML2DVISCOUS.h"

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
// PML2DVISCOUS element tcl command
// =======================================================================
void* OPS_PML2DVISCOUS()
{
	// check if the total number of arguments passed is correct
	if (OPS_GetNumRemainingInputArgs() < (5 + PML2DVISCOUS_NUM_PROPS + 3)) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element PML2DVISCOUS eleTag? [4 integer nodeTags] [PML2DVISCOUS_NUM_PARAMS material properties]\n";
		return 0;
	}

	// reading element tag and node numbers 
	int idata[5];
	int num = 5;
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
	double dData[PML2DVISCOUS_NUM_PROPS]; num = PML2DVISCOUS_NUM_PROPS;
	if (OPS_GetDoubleInput(&num, dData) < 0) {
		opserr << "WARNING: invalid double data\n";
		return 0;
	}
	// create a new PML2DVISCOUS element and add it to the Domain
	return new PML2DVISCOUS(idata[0], &idata[1], Newmark, dData);
}

// =======================================================================
// static data
// =======================================================================
Matrix  PML2DVISCOUS::tangent(PML2DVISCOUS_NUM_DOF, PML2DVISCOUS_NUM_DOF);
Matrix  PML2DVISCOUS::mass(PML2DVISCOUS_NUM_DOF, PML2DVISCOUS_NUM_DOF);
Matrix  PML2DVISCOUS::damping(PML2DVISCOUS_NUM_DOF, PML2DVISCOUS_NUM_DOF);
Vector  PML2DVISCOUS::resid(PML2DVISCOUS_NUM_DOF);
double  PML2DVISCOUS::eta = 0.;
double  PML2DVISCOUS::beta = 0.;
double  PML2DVISCOUS::gamma = 0.;
double  PML2DVISCOUS::dt = 0.;
int     PML2DVISCOUS::eleCount = 0;
// int     PML2DVISCOUS::numberOfElements = 0;


// =======================================================================
// null constructor
// =======================================================================
PML2DVISCOUS::PML2DVISCOUS()
	:Element(0, ELE_TAG_PML2DVISCOUS),
	connectedExternalNodes(PML2DVISCOUS_NUM_NODES),
	ubar(PML2DVISCOUS_NUM_DOF),
	ubart(PML2DVISCOUS_NUM_DOF)
{
	for (int i = 0; i < PML2DVISCOUS_NUM_NODES; i++) {
		nodePointers[i] = 0;
	}
	dt = 0;
	ubar.Zero();
	ubart.Zero();
	updateflag = 0;
	update_dt = 0;
	eta = 0;
	beta = 0;
	gamma = 0;
	
}

// =======================================================================
// Full constructor
// =======================================================================
PML2DVISCOUS::PML2DVISCOUS(int tag, int* nodeTags, double* nemwarks, double* eleData)
	:Element(tag, ELE_TAG_PML2DVISCOUS),
	connectedExternalNodes(PML2DVISCOUS_NUM_NODES),
	ubar(PML2DVISCOUS_NUM_DOF),
	ubart(PML2DVISCOUS_NUM_DOF)
{
	eleCount++;
	if (eleCount == 1) {
		opserr << "Perfectly Matched Layer 2D (PML) element with Viscous damping -  Written: W. Zhang, E. Taciroglu , A. Pakzad, P. Arduino, UCLA, UCLA, U.Washington, U.Washington 12/2020\n ";
	}
	// initialize node pointers
	for (int i = 0; i < PML2DVISCOUS_NUM_NODES; i++) {
		connectedExternalNodes(i) = nodeTags[i];
		nodePointers[i] = 0;
	}

	// initialize Newmark parameters
	eta   = nemwarks[0];
	beta  = nemwarks[1];
	gamma = nemwarks[2];

	// initialize material properties
	for (int i = 0; i < PML2DVISCOUS_NUM_PROPS; i++)
		props[i] = eleData[i];

	// initialize the ubar and ubart vectors to zero
	ubart.Zero();
	ubar.Zero();
	updateflag = 0;
	update_dt = 0;
}

// =======================================================================
//  destructor
// ======================================================================= 
PML2DVISCOUS::~PML2DVISCOUS()
{

}

// =======================================================================
// Set Domain
// =======================================================================
void  PML2DVISCOUS::setDomain(Domain* theDomain)
{
	int i ;
	Domainptr = theDomain;

	//node pointers
	for ( i=0; i<4; i++ ) 
		nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;

	this->DomainComponent::setDomain(theDomain);


	// 
	// set constant matrices by invoking fortran routine
	//

	//  double coords[PML2D_NUM_NODES][2];
	double coords[8];

	for (int i = 0; i < 4; i++ )  {
		const Vector &loc = nodePointers[i]->getCrds();
		coords[i*2] = loc(0);      
		coords[i*2+1] = loc(1);      
	} 

	int NDOFEL = 20;
	int NPROPS = 11;
	int MCRD = 2; 
	int NNODE = 4;

	pml2d_(K, 
		C, 
		M,  
		G, 
		&NDOFEL, 
		props, 
		&NPROPS, 
		coords, 
		&MCRD,
		&NNODE);
	if (this->getTag() == 51) {
		// save k matrix in file k.txt 
		std::ofstream kfile;
		kfile.open("k.txt");
		for (int i = 0; i < PML2DVISCOUS_NUM_DOF; i++) {
			for (int j = 0; j < PML2DVISCOUS_NUM_DOF; j++) {
				if (j== PML2DVISCOUS_NUM_DOF-1) {
					kfile << G[i*PML2DVISCOUS_NUM_DOF + j] << "\n";
				} else {
					kfile << G[i*PML2DVISCOUS_NUM_DOF + j] << ",";
				}
			}

		}
	}
	
}

// =======================================================================
// update
// =======================================================================
int PML2DVISCOUS::update(void)
{
	// check if the dt has changed
	if (fabs(Domainptr->getDT() - dt) > 1e-10) {
		update_dt = 1;
		dt = Domainptr->getDT();
	} else {
		update_dt = 0;
	}

	if (updateflag == 1) {
		// get u, v, a from nodes and calculate the ubar vector
		int loc = 0;
		double c1 = dt;
		double c2 = dt * dt * 0.5;
		double c3 = dt*dt*dt*((1.0/6.0)-eta);
		double c4 = dt*dt*dt*eta;
		for (int i = 0; i < PML2DVISCOUS_NUM_NODES; i++) {
			const Vector& uNode = nodePointers[i]->getDisp();
			const Vector& vNode = nodePointers[i]->getVel();
			const Vector& aNode = nodePointers[i]->getAccel();
			const Vector& atpdt = nodePointers[i]->getTrialAccel();
			for (int j = 0; j < 5; j++) {
				ubar(loc) = ubart(loc) + uNode(j)*c1 + vNode(j)*c2 + aNode(j)*c3 + atpdt(j)*c4; 
				loc++;
			}
		}
	}
	updateflag = 1;
	return 0;
}

// =======================================================================
//	return stiffness matrix 
// =======================================================================
const Matrix& PML2DVISCOUS::getTangentStiff()
{
	// check if the dt is changed to update the tangent stiffness matrix
	if (update_dt == 1) {
		double cg = eta*dt/beta;
		//keff = k + cg*g( k and g are symmetric matrices)
		for (int i = 0; i < PML2DVISCOUS_NUM_DOF*PML2DVISCOUS_NUM_DOF; i++) {
				Keff[i]= K[i] + cg * G[i];
		}
		
	}
	tangent.setData(K, PML2DVISCOUS_NUM_DOF, PML2DVISCOUS_NUM_DOF);
	return tangent;
}

// =======================================================================
//	return initial stiffness matrix 
// =======================================================================
const Matrix& PML2DVISCOUS::getInitialStiff()
{
	return this->getTangentStiff();
}

// =======================================================================
//	return mass matrix
// =======================================================================
const Matrix& PML2DVISCOUS::getMass()
{
	mass.setData(M, PML2DVISCOUS_NUM_DOF, PML2DVISCOUS_NUM_DOF);
	// mass.Zero();
	return mass;
}

// =======================================================================
//	return damping matrix
// =======================================================================
const Matrix& PML2DVISCOUS::getDamp()
{
	damping.setData(C, PML2DVISCOUS_NUM_DOF, PML2DVISCOUS_NUM_DOF);
	// damping.Zero();
	return damping;
}

// =======================================================================
// Ressisting force
// =======================================================================
//get residual
const Vector& PML2DVISCOUS::getResistingForce()
{
	// if (innertag==14) {
	// 	opserr << "getResistingForce function is called\n";
	// }
	int numNodalDOF = 5;
	static Vector theVector(PML2DVISCOUS_NUM_DOF);

	// get K into stiff
	tangent.setData(K, PML2DVISCOUS_NUM_DOF, PML2DVISCOUS_NUM_DOF);

	//
	// perform: R = K * u
	//

	int loc = 0;
	for (int i = 0; i < PML2DVISCOUS_NUM_NODES; i++) {
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
PML2DVISCOUS::getResistingForceIncInertia()
{
    // R += K*u
	static Vector theVector(PML2DVISCOUS_NUM_DOF);
	tangent.setData(K, PML2DVISCOUS_NUM_DOF, PML2DVISCOUS_NUM_DOF);

	int loc = 0;
	for (int i = 0; i < PML2DVISCOUS_NUM_NODES; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < 5; j++)
			theVector(loc++) = uNode(j);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);



	// R += M*a
	loc = 0;
	Node** theNodes = this->getNodePtrs();
	for (int i = 0; i < PML2DVISCOUS_NUM_NODES; i++) {
		const Vector& acc = theNodes[i]->getTrialAccel();
		for (int j = 0; j < 5; j++) {
			theVector(loc++) = acc(j);
		}
	}
	resid.addMatrixVector(1.0, this->getMass(), theVector, 1.0);

	// R += C*v
	loc = 0;
	for (int i = 0; i < PML2DVISCOUS_NUM_NODES; i++) {
		const Vector& vel = theNodes[i]->getTrialVel();
		for (int j = 0; j < 5; j++) {
			theVector(loc++) = vel[j];
		}
	}


	resid.addMatrixVector(1.0, this->getDamp(), theVector, 1.0);


	// R += G*ubar
	mass.setData(G, PML2DVISCOUS_NUM_DOF, PML2DVISCOUS_NUM_DOF);
	resid.addMatrixVector(1.0, mass, ubar, 1.0);
	
	
	return resid;
}

// =======================================================================
// get the number of external nodes
// =======================================================================
int  PML2DVISCOUS::getNumExternalNodes() const
{
	return PML2DVISCOUS_NUM_NODES;
}

// =======================================================================
// return connected external nodes
// =======================================================================
const ID& PML2DVISCOUS::getExternalNodes()
{
	return connectedExternalNodes;
}

// =======================================================================
// return node pointers
// =======================================================================
Node** PML2DVISCOUS::getNodePtrs(void)
{
	return nodePointers;
}

// =======================================================================
// return number of dofs
// =======================================================================
int  PML2DVISCOUS::getNumDOF()
{
	return PML2DVISCOUS_NUM_DOF;
}

// =======================================================================
// commit state
// =======================================================================
int  PML2DVISCOUS::commitState()
{
	int success = 0;
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PML2DVISCOUS::commitState () - failed in base class";
	}

	// set ubart to ubar
	for (int i = 0; i < PML2DVISCOUS_NUM_DOF; i++) {
		ubart(i) = ubar(i);
	}

	updateflag = 0;
	return success;
}

// =======================================================================
// revert to last commit 
// =======================================================================
int  PML2DVISCOUS::revertToLastCommit()
{
	int success = 0;

	// set ubar to ubart
	for (int i = 0; i < PML2DVISCOUS_NUM_DOF; i++) {
		ubar(i) = ubart(i);
	}

	return success;
}

// =======================================================================
// revert to start
// =======================================================================
int  PML2DVISCOUS::revertToStart()
{
	int success = 0;

	// set ubar and ubart to zero
	// for (int i = 0; i < PML2DVISCOUS_NUM_DOF; i++) {
	// 	ubar[i] = 0.0;
	// 	ubart[i] = 0.0;
	// }
	ubar.Zero();
	ubart.Zero();

	return success;
}

// =======================================================================
// add load
// =======================================================================
int PML2DVISCOUS::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return -1;
}

// =======================================================================
// add zero load
// =======================================================================
void  PML2DVISCOUS::zeroLoad()
{
	return;
}

// =======================================================================
// senself
// =======================================================================
int  PML2DVISCOUS::sendSelf(int commitTag,
	Channel& theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// PML2DVISCOUS packs its data into a Vector and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
	static Vector data(PML2DVISCOUS_NUM_PROPS + 4);
	data(0) = this->getTag();

	for (int ii = 1; ii <= PML2DVISCOUS_NUM_PROPS; ii++) {
		data(ii) = props[ii - 1];
	}
	data(PML2DVISCOUS_NUM_PROPS+1) = eta;
	data(PML2DVISCOUS_NUM_PROPS+2) = beta;
	data(PML2DVISCOUS_NUM_PROPS+3) = gamma;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML2DVISCOUS::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
	}


	// PML2DVISCOUS then sends the tags of its four nodes
	res += theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML2DVISCOUS::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	return res;
}

// =======================================================================
// recvself
// =======================================================================
int  PML2DVISCOUS::recvSelf(int commitTag,
	Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// PML2DVISCOUS creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(PML2DVISCOUS_NUM_PROPS + 4);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML2DVISCOUS::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));


	for (int ii = 1; ii <= PML2DVISCOUS_NUM_PROPS; ii++) {
		props[ii - 1] = data(ii);
	}

	eta   = data(PML2DVISCOUS_NUM_PROPS+1);
	beta  = data(PML2DVISCOUS_NUM_PROPS+2);
	gamma = data(PML2DVISCOUS_NUM_PROPS+3);

	// PML2DVISCOUS now receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML2DVISCOUS::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	return res;
}


// =======================================================================
// display
// =======================================================================
int PML2DVISCOUS::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
	// Get the end point display coords
	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);
	nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
	nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
	nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
	nodePointers[3]->getDisplayCrds(v4, fact, displayMode);


	// place values in coords matrix
	static Matrix coords(4, 2);
	for (int i = 0; i < 3; i++) {
		coords(0, i) = v1(i);
		coords(1, i) = v2(i);
		coords(2, i) = v3(i);
		coords(3, i) = v4(i);
	}

	// fill RGB vector
	static Vector values(4);
	for (int i = 0; i < 8; i++)
		values(i) = 1.0;

	// draw the cube
	return theViewer.drawCube(coords, values, this->getTag());
}

// =======================================================================
// setresponse
// =======================================================================
Response* PML2DVISCOUS::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;

	// char outputData[32];

	// output.tag("ElementOutput");
	// output.attr("eleType", "PML2DVISCOUS");
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
int PML2DVISCOUS::getResponse(int responseID, Information& eleInfo)
{
	// static Vector stresses(48);

	// if (responseID == 1)
	// 	return eleInfo.setVector(this->getResistingForce());

	return -1;
}

// =======================================================================
// set parameter
// =======================================================================
int PML2DVISCOUS::setParameter(const char** argv, int argc, Parameter& param)
{
	int res = -1;
	return res;
}

// =======================================================================
// update parameter
// =======================================================================
int PML2DVISCOUS::updateParameter(int parameterID, Information& info)
{
	int res = -1;
	return res;
}

// =======================================================================
// print
// =======================================================================
void  PML2DVISCOUS::Print(OPS_Stream &s, int flag) {
	
  	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "Element: " << this->getTag() << endln;
		s << "type: PML2DVISCOUS \n";
		s << "Nodes: " << connectedExternalNodes;
		s << "eta: " << eta << " beta: " << beta << " gamma: " << gamma << endln;
		s << endln;
		s << "Resisting Force (no inertia): " << this->getResistingForce();
  	}
    
 	 if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"PML2DVISCOUS\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
		for (int i = 1; i < 7; i++)
		s << connectedExternalNodes(i) << ", ";
		s << connectedExternalNodes(7) << "], ";
  	}
	return;
}