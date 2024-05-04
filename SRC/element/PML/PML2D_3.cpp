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

// Written by: Amin Pakzad, Pedro Arduino (parduino@uw.edu), and Adriano Trono
//
// Four node PML2D_3 element based on the Adriano Torino PhD thesis


#include "PML2D_3.h"

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
// PML2D_3 element tcl command
// =======================================================================
void* OPS_PML2D_3()
{
	// check if the total number of arguments passed is correct
	if (OPS_GetNumRemainingInputArgs() < (13)) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: element PML2D_3 eleTag? [5 integer nodeTags] [7 PML properties]\n";
		return 0;
	}

	// reading element tag and node numbers 
	int idata[6];
	int num = 6;
	if (OPS_GetIntInput(&num, idata) < 0) {
		opserr << "WARNING: invalid integer data : could be the tag or the node numbers \n";
		return 0;
	}

	// reading material properties
	double dData[7]; num = 7;
	if (OPS_GetDoubleInput(&num, dData) < 0) {
		opserr << "WARNING: invalid double data\n";
		return 0;
	}
	double E                = dData[0];
	double nu				= dData[1];
	double rho		     	= dData[2]; 
	double pmlthicknessx    = dData[3];
	double pmlthicknessy    = dData[4];
	double Halfwidth	    = dData[5];
	double Depth		    = dData[6];


	double G = E / (2 * (1 + nu));
	double Vs = sqrt(G / rho);
	double data[3] = {1.25,1.e-3,Vs};
    num = OPS_GetNumRemainingInputArgs();

	if (num > 3) {num = 3;}
	if (num > 0 ) {
		if (OPS_GetDoubleInput(&num,data) < 0) {
	    	opserr<<"WARNING: invalid double data\n";
	    	return 0;
		}
	}
	double r0 = data[0];
	double R  = data[1];
	double Vc = data[2];
	// create a new PML2D_3 element and add it to the Domain

	return new PML2D_3(idata[0], &idata[1], E, nu, rho, pmlthicknessx, pmlthicknessy, Halfwidth, Depth, r0, R, Vc);
}

// =======================================================================
// static data
// =======================================================================
Matrix  PML2D_3::tangent(PML2D_3_NUM_DOF, PML2D_3_NUM_DOF);
Vector  PML2D_3::resid(PML2D_3_NUM_DOF);
int     PML2D_3::eleCount = 0;

// =======================================================================
// null constructor
// =======================================================================
PML2D_3::PML2D_3()
	:Element(0, ELE_TAG_PML2D_3),
	connectedExternalNodes(PML2D_3_NUM_NODES)
{
	for (int i = 0; i < PML2D_3_NUM_NODES; i++) {
		nodePointers[i] = 0;
	}
}

// =======================================================================
// Full constructor
// =======================================================================
PML2D_3::PML2D_3(int tag, int* nodeTags, double E, double nu,
	double rho, double pmlthicknessx, double pmlthicknessy, double Halfwidth, 
	double Depth, double r0, double R, double Vc)
	:Element(tag, ELE_TAG_PML2D_3),
	connectedExternalNodes(PML2D_3_NUM_NODES), E(E), nu(nu), rho(rho), 
	pmlthicknessx(pmlthicknessx), pmlthicknessy(pmlthicknessy), Halfwidth(Halfwidth),
	Depth(Depth), r0(r0), R(R), Vc(Vc)
{
	eleCount++;
	if (eleCount == 1) {
	opserr << "Perfectly Matched Layer 2D_3 (PML) element -  Written: A. Trono, A. Pakzad, P. Arduino, National University of Córdoba, U.Washington, U.Washington\n ";
	}
	// initialize node pointers
	for (int i = 0; i < PML2D_3_NUM_NODES; i++) {
		connectedExternalNodes(i) = nodeTags[i];
		nodePointers[i] = 0;
	}
}

// =======================================================================
//  destructor
// ======================================================================= 
PML2D_3::~PML2D_3()
{

}
// =======================================================================
// Set Domain
// =======================================================================
void  PML2D_3::setDomain(Domain* theDomain)
{
	// find xi and yi which are the smallest x and y coordinate among nodes
	double xi = 1.0e+10;
	double yj = 1.0e+10;
	for (int i = 0; i < PML2D_3_NUM_NODES; i++) {
		const Vector &loc = theDomain->getNode(connectedExternalNodes(i))->getCrds();
		if (loc(0) < xi) { xi = loc(0); }
		if (loc(1) < yj) { yj = loc(1); }	
	}

	// rearanging the nodetags
	int nodetags[PML2D_3_NUM_NODES];
	double eps = 1.0e-4;
	for (int i = 0; i < PML2D_3_NUM_NODES; i++) {
		const Vector &loc = theDomain->getNode(connectedExternalNodes(i))->getCrds();
		if ((loc(0) < (xi +eps)) && (loc(0) > (xi-eps))) {
			if ((loc(1) < (yj +eps)) && (loc(1) > (yj-eps))) {
				nodetags[1] = connectedExternalNodes(i);
			} else {
				nodetags[0] = connectedExternalNodes(i);
			}
		} else {
			if ((loc(1) < (yj +eps)) && (loc(1) > (yj-eps))) {
			nodetags[2] = connectedExternalNodes(i);
			} else {
				if (theDomain->getNode(connectedExternalNodes(i))->getNumberDOF() == 2) {
					nodetags[3] = connectedExternalNodes(i);
				} else {
					nodetags[4] = connectedExternalNodes(i);
				}
			}
		}
 	}
	

	// set connected external nodes
	for (int i = 0; i < PML2D_3_NUM_NODES; i++) {
		connectedExternalNodes(i) = nodetags[i];
	}

	//node pointers
	for ( int i=0; i<PML2D_3_NUM_NODES; i++ ) {
		nodePointers[i] = theDomain->getNode( connectedExternalNodes(i)) ;
	}
	this->DomainComponent::setDomain(theDomain);


	// obtaining the zone 
	int    WidthSign  = 0;
	int    DepthSign  = 0;
	double x = nodePointers[4]->getCrds()(0);
	double y = nodePointers[4]->getCrds()(1);
	double LXPML = 0.0;
	double LYPML = 0.0;
	if ( x > (Halfwidth - eps))   {WidthSign = 1;}
	if ( x < (-Halfwidth + eps))  {WidthSign = -1;}
	if ( y < (-Depth))            {DepthSign = -1;}
	if ( y >= 0.0) {
		opserr << "PML2D_3 element is not defined for y>=0.0\n";
		exit(-1);
	}
	int zone; 
	if (DepthSign==0  && WidthSign==1 ) {LXPML = pmlthicknessx; LYPML = Depth;        zone = 2;} // zone 2
	if (DepthSign==0  && WidthSign==-1) {LXPML =-pmlthicknessx; LYPML = Depth;        zone = 1;} // zone 1
	if (DepthSign==-1 && WidthSign==0 ) {LXPML = 2*Halfwidth  ; LYPML =-pmlthicknessy;zone = 3;} // zone 3
	if (DepthSign==-1 && WidthSign==1 ) {LXPML = pmlthicknessx; LYPML =-pmlthicknessy;zone = 5;} // zone 5         
	if (DepthSign==-1 && WidthSign==-1) {LXPML =-pmlthicknessx; LYPML =-pmlthicknessy;zone = 4;} // zone 4

	// calculate alpha0 and beta0
	int m = 2;
	int n = 2;
	double alpha0_y = 0;
	double beta0_y  = 0;
	double alpha0_x = 0;
	double beta0_x  = 0;

	if (zone >= 3) {
		alpha0_y = (((n+1)*r0)/(2*pmlthicknessy))*log(1./R);
		beta0_y  = (m+1)*Vc/(2*pmlthicknessy)*log(1./R);
	}
	if (zone !=3) {
		alpha0_x = (((n+1)*r0)/(2*pmlthicknessx))*log(1./R);
		beta0_x  = (m+1)*Vc/(2*pmlthicknessx)*log(1./R);
	}

	//  double coords[PML2D_NUM_NODES][2];
	double coords[8];

	for (int i = 0; i < 4; i++ )  {
		const Vector &loc = nodePointers[i]->getCrds();
		coords[i*2] = loc(0);      
		coords[i*2+1] = loc(1);      
	} 
	this->ComputeK(K,coords, beta0_x, beta0_y, LXPML, LYPML, xi, yj, rho, E, nu);
	this->ComputeM(M,coords, alpha0_x, alpha0_y, LXPML, LYPML, xi, yj, rho, E, nu);
	this->ComputeC(C,coords, alpha0_x, alpha0_y, beta0_x, beta0_y, LXPML, LYPML, xi, yj, rho, E, nu);
	// opserr << "\n\nElement Tag = "<<this->getTag() << endln;
	// opserr << "E = " << E << endln;
	// opserr << "nu = " << nu << endln;
	// opserr << "rho = " << rho << endln;
	// opserr << "pmlthicknessx = " << pmlthicknessx << endln;
	// opserr << "pmlthicknessy = " << pmlthicknessy << endln;
	// opserr << "Halfwidth = " << Halfwidth << endln;
	// opserr << "Depth = " << Depth << endln;
	// opserr << "r0 = " << r0 << endln;
	// opserr << "R = " << R << endln;
	// opserr << "Vc = " << Vc << endln;
	// opserr << "LXPML = " << LXPML << endln;
	// opserr << "LYPML = " << LYPML << endln;
	// opserr << "xi = " << xi << endln;
	// opserr << "yj = " << yj << endln;
	// opserr << "alpha0_x = " << alpha0_x << endln;
	// opserr << "alpha0_y = " << alpha0_y << endln;
	// opserr << "beta0_x = " << beta0_x << endln;
	// opserr << "beta0_y = " << beta0_y << endln;
	// opserr << "zone = " << zone << endln;
}

// =======================================================================
// update
// =======================================================================
int PML2D_3::update(void)
{
	return 0;
}

// =======================================================================
//	return stiffness matrix 
// =======================================================================
const Matrix& PML2D_3::getTangentStiff()
{
	tangent.setData(K, PML2D_3_NUM_DOF, PML2D_3_NUM_DOF);
	return tangent;
}

// =======================================================================
//	return initial stiffness matrix 
// =======================================================================
const Matrix& PML2D_3::getInitialStiff()
{
	return this->getTangentStiff();
}

// =======================================================================
//	return mass matrix
// =======================================================================
const Matrix& PML2D_3::getMass()
{
	tangent.setData(M, PML2D_3_NUM_DOF, PML2D_3_NUM_DOF);
	return tangent;
}

// =======================================================================
//	return damping matrix
// =======================================================================
const Matrix& PML2D_3::getDamp()
{
	tangent.setData(C, PML2D_3_NUM_DOF, PML2D_3_NUM_DOF);
	return tangent;
}

// =======================================================================
// Ressisting force
// =======================================================================
//get residual
const Vector& PML2D_3::getResistingForce()
{
	// if (innertag==14) {
	// 	opserr << "getResistingForce function is called\n";
	// }
	int numNodalDOF = 2;
	static Vector theVector(PML2D_3_NUM_DOF);

	// get K into stiff
	tangent.setData(K, PML2D_3_NUM_DOF, PML2D_3_NUM_DOF);

	//
	// perform: R = K * u
	//

	int loc = 0;
	for (int i = 0; i < 4; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < numNodalDOF; j++)
			theVector(loc++) = uNode(j);
	}
	{
		const Vector& uNode = nodePointers[4]->getTrialDisp();
		theVector(loc++) = uNode(0);
		theVector(loc++) = uNode(1);
		theVector(loc++) = uNode(2);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);
	return resid;
}


// =======================================================================
//
// =======================================================================
//get residual with inertia terms
const Vector&
PML2D_3::getResistingForceIncInertia()
{
    // R += K*u
	static Vector theVector(PML2D_3_NUM_DOF);
	tangent.setData(K, PML2D_3_NUM_DOF, PML2D_3_NUM_DOF);

	int loc = 0;
	for (int i = 0; i < 4; i++) {
		const Vector& uNode = nodePointers[i]->getTrialDisp();
		for (int j = 0; j < 2; j++)
			theVector(loc++) = uNode(j);
	}
	{
		const Vector& uNode = nodePointers[4]->getTrialDisp();
		theVector(loc++) = uNode(0);
		theVector(loc++) = uNode(1);
		theVector(loc++) = uNode(2);
	}
	resid.addMatrixVector(0.0, tangent, theVector, 1.0);


	// R += M*a
	loc = 0;
	for (int i = 0; i < 4; i++) {
		const Vector& uNode = nodePointers[i]->getTrialAccel();
		for (int j = 0; j < 2; j++)
			theVector(loc++) = uNode(j);
	}
	{
		const Vector& uNode = nodePointers[4]->getTrialAccel();
		theVector(loc++) = uNode(0);
		theVector(loc++) = uNode(1);
		theVector(loc++) = uNode(2);
	}
	resid.addMatrixVector(1.0, this->getMass(), theVector, 1.0);


	// R += C*v
	loc = 0;
	for (int i = 0; i < 4; i++) {
		const Vector& uNode = nodePointers[i]->getTrialVel();
		for (int j = 0; j < 2; j++)
			theVector(loc++) = uNode(j);
	}
	{
		const Vector& uNode = nodePointers[4]->getTrialVel();
		theVector(loc++) = uNode(0);
		theVector(loc++) = uNode(1);
		theVector(loc++) = uNode(2);
	}

	resid.addMatrixVector(1.0, this->getDamp(), theVector, 1.0);

	return resid;
}

// =======================================================================
// get the number of external nodes
// =======================================================================
int  PML2D_3::getNumExternalNodes() const
{
	// opserr << "PML2D_3::getNumExternalNodes is called\n";
	return PML2D_3_NUM_NODES;
}

// =======================================================================
// return connected external nodes
// =======================================================================
const ID& PML2D_3::getExternalNodes()
{
	return connectedExternalNodes;
}

// =======================================================================
// return node pointers
// =======================================================================
Node** PML2D_3::getNodePtrs(void)
{
	return nodePointers;
}

// =======================================================================
// return number of dofs
// =======================================================================
int  PML2D_3::getNumDOF()
{
	return PML2D_3_NUM_DOF;
}

// =======================================================================
// commit state
// =======================================================================
int  PML2D_3::commitState()
{
	int success = 0;
	if ((success = this->Element::commitState()) != 0) {
		opserr << "PML2D_3::commitState () - failed in base class";
	}
	return success;
}

// =======================================================================
// revert to last commit 
// =======================================================================
int  PML2D_3::revertToLastCommit()
{
	int success = 0;

	return success;
}

// =======================================================================
// revert to start
// =======================================================================
int  PML2D_3::revertToStart()
{
	int success = 0;

	return success;
}

// =======================================================================
// add load
// =======================================================================
int PML2D_3::addLoad(ElementalLoad* theLoad, double loadFactor)
{
	return -1;
}

// =======================================================================
// add zero load
// =======================================================================
void  PML2D_3::zeroLoad()
{
	return;
}

// =======================================================================
// senself
// =======================================================================
int  PML2D_3::sendSelf(int commitTag,
	Channel& theChannel)
{
	int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();
	static Vector data(11);
    data(0) = this->getTag();
	data(1) = E;
	data(2) = nu;
	data(3) = rho;
	data(4) = pmlthicknessx;
	data(5) = pmlthicknessy;
	data(6) = Halfwidth;
	data(7) = Depth;
	data(8) = r0;
	data(9) = R;
	data(10) = Vc;
	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML2D_3::sendSelf() - " << this->getTag() << " failed to send Vector\n";
		return res;
  	}
  	// PML2D_5 then sends the tags of its four nodes
	res += theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML2D_3::sendSelf() - " << this->getTag() << " failed to send ID\n";
		return res;
	}

	return res;
}

// =======================================================================
// recvself
// =======================================================================
int  PML2D_3::recvSelf(int commitTag,
	Channel& theChannel,
	FEM_ObjectBroker& theBroker)
{
    int res = 0;
	int dataTag = this->getDbTag();
	static Vector data(11);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		opserr << "WARNING PML2D_3::recvSelf() - failed to receive Vector\n";
		return res;
	}

	this->setTag((int)data(0));
	E = data(1);
	nu = data(2);
	rho = data(3);
	pmlthicknessx = data(4);
	pmlthicknessy = data(5);
	Halfwidth = data(6);
	Depth = data(7);
	r0 = data(8);
	R = data(9);
	Vc = data(10);

	res += theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		opserr << "WARNING PML2D_3::recvSelf() - " << this->getTag() << " failed to receive ID\n";
		return res;
	}

	return res;
}


// =======================================================================
// display
// =======================================================================
int PML2D_3::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
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
Response* PML2D_3::setResponse(const char** argv, int argc, OPS_Stream& output)
{
	Response* theResponse = 0;

	// char outputData[32];

	// output.tag("ElementOutput");
	// output.attr("eleType", "PML2D_3");
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
int PML2D_3::getResponse(int responseID, Information& eleInfo)
{
	// static Vector stresses(48);

	// if (responseID == 1)
	// 	return eleInfo.setVector(this->getResistingForce());

	return -1;
}

// =======================================================================
// set parameter
// =======================================================================
int PML2D_3::setParameter(const char** argv, int argc, Parameter& param)
{
	int res = -1;
	return res;
}

// =======================================================================
// update parameter
// =======================================================================
int PML2D_3::updateParameter(int parameterID, Information& info)
{
	int res = -1;
	return res;
}

// =======================================================================
// print
// =======================================================================
void  PML2D_3::Print(OPS_Stream &s, int flag) {
	
  	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "Element: " << this->getTag() << endln;
		s << "type: PML2D_3 \n";
		s << "Nodes: " << connectedExternalNodes;
		s << endln;
		s << "Resisting Force (no inertia): " << this->getResistingForce();
  	}
    
 	 if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{";
		s << "\"name\": " << this->getTag() << ", ";
		s << "\"type\": \"PML2D_3\", ";
		s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
		for (int i = 1; i < 7; i++)
		s << connectedExternalNodes(i) << ", ";
		s << connectedExternalNodes(7) << "], ";
  	}
	return;
}
// =======================================================================
// Compute K matrix 
// =======================================================================
void PML2D_3::ComputeK(double* K,double* XYelement, double beta_0_x, double beta_0_y, double L_PML_x,
                                double L_PML_y, double xi, double yj, double rho, double E, double nu) {
    // Define variables
    double a, b;  // Variables a and b

    // Calculate variables a and b
    a = fabs(XYelement[6] - XYelement[0]);
    b = fabs(XYelement[1] - XYelement[3]);

    // Continue with the rest of the function...
    // Add your C++ code here
	K[0] = (a * b * beta_0_x * beta_0_y * rho * (a * a + 5 * a * xi + 10 * xi * xi) * (6 * b * b + 15 * b * yj + 10 * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[1] = 0;
	K[2] = (a * b * beta_0_x * beta_0_y * rho * (a * a + 5 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[3] = 0;
	K[4] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[5] = 0;
	K[6] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (6 * b * b + 15 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[7] = 0;
	K[8] =  - (b * beta_0_y * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[9] = 0;
	K[10] = (a * beta_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[11] = 0;
	K[12] = (a * b * beta_0_x * beta_0_y * rho * (a * a + 5 * a * xi + 10 * xi * xi) * (6 * b * b + 15 * b * yj + 10 * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[13] = 0;
	K[14] = (a * b * beta_0_x * beta_0_y * rho * (a * a + 5 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[15] = 0;
	K[16] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[17] = 0;
	K[18] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (6 * b * b + 15 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[19] = 0;
	K[20] = (a * beta_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[21] =  - (b * beta_0_y * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[22] = (a * b * beta_0_x * beta_0_y * rho * (a * a + 5 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[23] = 0;
	K[24] = (a * b * beta_0_x * beta_0_y * rho * (a * a + 5 * a * xi + 10 * xi * xi) * (b * b + 5 * b * yj + 10 * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[25] = 0;
	K[26] = (a * b * beta_0_x * beta_0_y * rho * (b * b + 5 * b * yj + 10 * yj * yj) * (3 * a * a + 10 * a * xi + 10 * xi * xi)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[27] = 0;
	K[28] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[29] = 0;
	K[30] =  - (b * beta_0_y * (b * b + 4 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[31] = 0;
	K[32] =  - (a * beta_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[33] = 0;
	K[34] = (a * b * beta_0_x * beta_0_y * rho * (a * a + 5 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[35] = 0;
	K[36] = (a * b * beta_0_x * beta_0_y * rho * (a * a + 5 * a * xi + 10 * xi * xi) * (b * b + 5 * b * yj + 10 * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[37] = 0;
	K[38] = (a * b * beta_0_x * beta_0_y * rho * (b * b + 5 * b * yj + 10 * yj * yj) * (3 * a * a + 10 * a * xi + 10 * xi * xi)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[39] = 0;
	K[40] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[41] = 0;
	K[42] =  - (a * beta_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[43] =  - (b * beta_0_y * (b * b + 4 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[44] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[45] = 0;
	K[46] = (a * b * beta_0_x * beta_0_y * rho * (b * b + 5 * b * yj + 10 * yj * yj) * (3 * a * a + 10 * a * xi + 10 * xi * xi)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[47] = 0;
	K[48] = (a * b * beta_0_x * beta_0_y * rho * (b * b + 5 * b * yj + 10 * yj * yj) * (6 * a * a + 15 * a * xi + 10 * xi * xi)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[49] = 0;
	K[50] = (a * b * beta_0_x * beta_0_y * rho * (6 * a * a + 15 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[51] = 0;
	K[52] = (b * beta_0_y * (b * b + 4 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[53] = 0;
	K[54] =  - (a * beta_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[55] = 0;
	K[56] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[57] = 0;
	K[58] = (a * b * beta_0_x * beta_0_y * rho * (b * b + 5 * b * yj + 10 * yj * yj) * (3 * a * a + 10 * a * xi + 10 * xi * xi)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[59] = 0;
	K[60] = (a * b * beta_0_x * beta_0_y * rho * (b * b + 5 * b * yj + 10 * yj * yj) * (6 * a * a + 15 * a * xi + 10 * xi * xi)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[61] = 0;
	K[62] = (a * b * beta_0_x * beta_0_y * rho * (6 * a * a + 15 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[63] = 0;
	K[64] =  - (a * beta_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[65] = (b * beta_0_y * (b * b + 4 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[66] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (6 * b * b + 15 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[67] = 0;
	K[68] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[69] = 0;
	K[70] = (a * b * beta_0_x * beta_0_y * rho * (6 * a * a + 15 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[71] = 0;
	K[72] = (a * b * beta_0_x * beta_0_y * rho * (6 * a * a + 15 * a * xi + 10 * xi * xi) * (6 * b * b + 15 * b * yj + 10 * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[73] = 0;
	K[74] = (b * beta_0_y * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[75] = 0;
	K[76] = (a * beta_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[77] = 0;
	K[78] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (6 * b * b + 15 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[79] = 0;
	K[80] = (a * b * beta_0_x * beta_0_y * rho * (3 * a * a + 10 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[81] = 0;
	K[82] = (a * b * beta_0_x * beta_0_y * rho * (6 * a * a + 15 * a * xi + 10 * xi * xi) * (3 * b * b + 10 * b * yj + 10 * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[83] = 0;
	K[84] = (a * b * beta_0_x * beta_0_y * rho * (6 * a * a + 15 * a * xi + 10 * xi * xi) * (6 * b * b + 15 * b * yj + 10 * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[85] = 0;
	K[86] = (a * beta_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[87] = (b * beta_0_y * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[88] =  - (b * beta_0_y * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[89] = 0;
	K[90] =  - (b * beta_0_y * (b * b + 4 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[91] = 0;
	K[92] = (b * beta_0_y * (b * b + 4 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[93] = 0;
	K[94] = (b * beta_0_y * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[95] = 0;
	K[96] = (a * b * beta_0_x * beta_0_y * (nu * nu - 1) * (a * a + 3 * a * xi + 3 * xi * xi) * (b * b + 3 * b * yj + 3 * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[97] = (a * b * beta_0_x * beta_0_y * nu * (nu + 1) * (a * a + 3 * a * xi + 3 * xi * xi) * (b * b + 3 * b * yj + 3 * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[98] = 0;
	K[99] = 0;
	K[100] = (a * beta_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[101] = 0;
	K[102] =  - (a * beta_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[103] = 0;
	K[104] =  - (a * beta_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[105] = 0;
	K[106] = (a * beta_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[107] = (a * b * beta_0_x * beta_0_y * nu * (nu + 1) * (a * a + 3 * a * xi + 3 * xi * xi) * (b * b + 3 * b * yj + 3 * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[108] = (a * b * beta_0_x * beta_0_y * (nu * nu - 1) * (a * a + 3 * a * xi + 3 * xi * xi) * (b * b + 3 * b * yj + 3 * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	K[109] = 0;
	K[110] = (a * beta_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[111] =  - (b * beta_0_y * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[112] =  - (a * beta_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[113] =  - (b * beta_0_y * (b * b + 4 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[114] =  - (a * beta_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[115] = (b * beta_0_y * (b * b + 4 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[116] = (a * beta_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	K[117] = (b * beta_0_y * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	K[118] = 0;
	K[119] = 0;
	K[120] =  - (2 * a * b * beta_0_x * beta_0_y * (nu + 1) * (a * a + 3 * a * xi + 3 * xi * xi) * (b * b + 3 * b * yj + 3 * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);

}

// =======================================================================
// Compute M matrix
// =======================================================================

void PML2D_3::ComputeM(double* M,double* XYelement, double alpha_0_x, double alpha_0_y, double L_PML_x,
                                double L_PML_y, double xi, double yj, double rho, double E, double nu) {
	// Define variables
    double a, b;  // Variables a and b

    // Calculate variables a and b
    a = fabs(XYelement[6] - XYelement[0]);
    b = fabs(XYelement[1] - XYelement[3]);


	// Calculate the mass matrix
	M[0] = (a * b * rho * (10 * L_PML_x * L_PML_x + alpha_0_x * a * a + 5 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 6 * alpha_0_y * b * b + 15 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[1] = 0;
	M[2] = (a * b * rho * (10 * L_PML_x * L_PML_x + alpha_0_x * a * a + 5 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[3] = 0;
	M[4] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[5] = 0;
	M[6] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 6 * alpha_0_y * b * b + 15 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[7] = 0;
	M[8] = 0;
	M[9] = 0;
	M[10] = 0;
	M[11] = 0;
	M[12] = (a * b * rho * (10 * L_PML_x * L_PML_x + alpha_0_x * a * a + 5 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 6 * alpha_0_y * b * b + 15 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[13] = 0;
	M[14] = (a * b * rho * (10 * L_PML_x * L_PML_x + alpha_0_x * a * a + 5 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[15] = 0;
	M[16] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[17] = 0;
	M[18] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 6 * alpha_0_y * b * b + 15 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[19] = 0;
	M[20] = 0;
	M[21] = 0;
	M[22] = (a * b * rho * (10 * L_PML_x * L_PML_x + alpha_0_x * a * a + 5 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[23] = 0;
	M[24] = (a * b * rho * (10 * L_PML_x * L_PML_x + alpha_0_x * a * a + 5 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + alpha_0_y * b * b + 5 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[25] = 0;
	M[26] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + alpha_0_y * b * b + 5 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[27] = 0;
	M[28] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[29] = 0;
	M[30] = 0;
	M[31] = 0;
	M[32] = 0;
	M[33] = 0;
	M[34] = (a * b * rho * (10 * L_PML_x * L_PML_x + alpha_0_x * a * a + 5 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[35] = 0;
	M[36] = (a * b * rho * (10 * L_PML_x * L_PML_x + alpha_0_x * a * a + 5 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + alpha_0_y * b * b + 5 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[37] = 0;
	M[38] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + alpha_0_y * b * b + 5 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[39] = 0;
	M[40] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[41] = 0;
	M[42] = 0;
	M[43] = 0;
	M[44] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[45] = 0;
	M[46] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + alpha_0_y * b * b + 5 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[47] = 0;
	M[48] = (a * b * rho * (10 * L_PML_x * L_PML_x + 6 * alpha_0_x * a * a + 15 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + alpha_0_y * b * b + 5 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[49] = 0;
	M[50] = (a * b * rho * (10 * L_PML_x * L_PML_x + 6 * alpha_0_x * a * a + 15 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[51] = 0;
	M[52] = 0;
	M[53] = 0;
	M[54] = 0;
	M[55] = 0;
	M[56] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[57] = 0;
	M[58] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + alpha_0_y * b * b + 5 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[59] = 0;
	M[60] = (a * b * rho * (10 * L_PML_x * L_PML_x + 6 * alpha_0_x * a * a + 15 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + alpha_0_y * b * b + 5 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[61] = 0;
	M[62] = (a * b * rho * (10 * L_PML_x * L_PML_x + 6 * alpha_0_x * a * a + 15 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[63] = 0;
	M[64] = 0;
	M[65] = 0;
	M[66] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 6 * alpha_0_y * b * b + 15 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[67] = 0;
	M[68] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[69] = 0;
	M[70] = (a * b * rho * (10 * L_PML_x * L_PML_x + 6 * alpha_0_x * a * a + 15 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[71] = 0;
	M[72] = (a * b * rho * (10 * L_PML_x * L_PML_x + 6 * alpha_0_x * a * a + 15 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 6 * alpha_0_y * b * b + 15 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[73] = 0;
	M[74] = 0;
	M[75] = 0;
	M[76] = 0;
	M[77] = 0;
	M[78] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 6 * alpha_0_y * b * b + 15 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[79] = 0;
	M[80] = (a * b * rho * (10 * L_PML_x * L_PML_x + 3 * alpha_0_x * a * a + 10 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (3600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[81] = 0;
	M[82] = (a * b * rho * (10 * L_PML_x * L_PML_x + 6 * alpha_0_x * a * a + 15 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 3 * alpha_0_y * b * b + 10 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (1800 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[83] = 0;
	M[84] = (a * b * rho * (10 * L_PML_x * L_PML_x + 6 * alpha_0_x * a * a + 15 * alpha_0_x * a * xi + 10 * alpha_0_x * xi * xi) * (10 * L_PML_y * L_PML_y + 6 * alpha_0_y * b * b + 15 * alpha_0_y * b * yj + 10 * alpha_0_y * yj * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[85] = 0;
	M[86] = 0;
	M[87] = 0;
	M[88] = 0;
	M[89] = 0;
	M[90] = 0;
	M[91] = 0;
	M[92] = 0;
	M[93] = 0;
	M[94] = 0;
	M[95] = 0;
	M[96] = (a * b * (nu * nu - 1) * (3 * L_PML_x * L_PML_x + alpha_0_x * a * a + 3 * alpha_0_x * a * xi + 3 * alpha_0_x * xi * xi) * (3 * L_PML_y * L_PML_y + alpha_0_y * b * b + 3 * alpha_0_y * b * yj + 3 * alpha_0_y * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[97] = (a * b * nu * (nu + 1) * (3 * L_PML_x * L_PML_x + alpha_0_x * a * a + 3 * alpha_0_x * a * xi + 3 * alpha_0_x * xi * xi) * (3 * L_PML_y * L_PML_y + alpha_0_y * b * b + 3 * alpha_0_y * b * yj + 3 * alpha_0_y * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[98] = 0;
	M[99] = 0;
	M[100] = 0;
	M[101] = 0;
	M[102] = 0;
	M[103] = 0;
	M[104] = 0;
	M[105] = 0;
	M[106] = 0;
	M[107] = (a * b * nu * (nu + 1) * (3 * L_PML_x * L_PML_x + alpha_0_x * a * a + 3 * alpha_0_x * a * xi + 3 * alpha_0_x * xi * xi) * (3 * L_PML_y * L_PML_y + alpha_0_y * b * b + 3 * alpha_0_y * b * yj + 3 * alpha_0_y * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[108] = (a * b * (nu * nu - 1) * (3 * L_PML_x * L_PML_x + alpha_0_x * a * a + 3 * alpha_0_x * a * xi + 3 * alpha_0_x * xi * xi) * (3 * L_PML_y * L_PML_y + alpha_0_y * b * b + 3 * alpha_0_y * b * yj + 3 * alpha_0_y * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	M[109] = 0;
	M[110] = 0;
	M[111] = 0;
	M[112] = 0;
	M[113] = 0;
	M[114] = 0;
	M[115] = 0;
	M[116] = 0;
	M[117] = 0;
	M[118] = 0;
	M[119] = 0;
	M[120] =  - (2 * a * b * (nu + 1) * (3 * L_PML_x * L_PML_x + alpha_0_x * a * a + 3 * alpha_0_x * a * xi + 3 * alpha_0_x * xi * xi) * (3 * L_PML_y * L_PML_y + alpha_0_y * b * b + 3 * alpha_0_y * b * yj + 3 * alpha_0_y * yj * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
}


// =================================================================================================
// Compute the matrix C
// =================================================================================================
void PML2D_3::ComputeC(double* C,double* XYelement, double alpha_0_x, double alpha_0_y,
							  double beta_0_x, double beta_0_y, double L_PML_x,double L_PML_y, 
								double xi, double yj, double rho, double E, double nu) { 
    // Define variables
    double a, b;  // Variables a and b

    // Calculate variables a and b
    a = fabs(XYelement[6] - XYelement[0]);
    b = fabs(XYelement[1] - XYelement[3]);

	C[0] = (a * b * rho * (L_PML_y * L_PML_y * a * a * beta_0_x + 10 * L_PML_y * L_PML_y * beta_0_x * xi * xi + 10 * L_PML_x * L_PML_x * beta_0_y * yj * yj + 5 * L_PML_y * L_PML_y * a * beta_0_x * xi + a * a * alpha_0_x * beta_0_y * yj * yj + a * a * alpha_0_y * beta_0_x * yj * yj + 10 * alpha_0_x * beta_0_y * xi * xi * yj * yj + 10 * alpha_0_y * beta_0_x * xi * xi * yj * yj + 5 * a * alpha_0_x * beta_0_y * xi * yj * yj + 5 * a * alpha_0_y * beta_0_x * xi * yj * yj)) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * b * b * rho * (20 * L_PML_x * L_PML_x * beta_0_y * yj + 2 * a * a * alpha_0_x * beta_0_y * yj + 2 * a * a * alpha_0_y * beta_0_x * yj + 20 * alpha_0_x * beta_0_y * xi * xi * yj + 20 * alpha_0_y * beta_0_x * xi * xi * yj + 10 * a * alpha_0_x * beta_0_y * xi * yj + 10 * a * alpha_0_y * beta_0_x * xi * yj)) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * b * b * b * rho * (10 * L_PML_x * L_PML_x * beta_0_y + a * a * alpha_0_x * beta_0_y + a * a * alpha_0_y * beta_0_x + 10 * alpha_0_x * beta_0_y * xi * xi + 10 * alpha_0_y * beta_0_x * xi * xi + 5 * a * alpha_0_x * beta_0_y * xi + 5 * a * alpha_0_y * beta_0_x * xi)) / (150 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[1] = 0;
	C[2] = (a * a * a * b * beta_0_x * rho) / (180 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (60 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[3] = 0;
	C[4] = (a * a * a * b * beta_0_x * rho) / (120 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (120 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (36 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (36 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[5] = 0;
	C[6] = (a * a * a * b * beta_0_x * rho) / (60 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (30 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (12 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[7] = 0;
	C[8] =  - b / 2 - (alpha_0_y * b * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[9] = 0;
	C[10] = a / 2 + (a * alpha_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[11] = 0;
	C[12] = (a * b * rho * (L_PML_y * L_PML_y * a * a * beta_0_x + 10 * L_PML_y * L_PML_y * beta_0_x * xi * xi + 10 * L_PML_x * L_PML_x * beta_0_y * yj * yj + 5 * L_PML_y * L_PML_y * a * beta_0_x * xi + a * a * alpha_0_x * beta_0_y * yj * yj + a * a * alpha_0_y * beta_0_x * yj * yj + 10 * alpha_0_x * beta_0_y * xi * xi * yj * yj + 10 * alpha_0_y * beta_0_x * xi * xi * yj * yj + 5 * a * alpha_0_x * beta_0_y * xi * yj * yj + 5 * a * alpha_0_y * beta_0_x * xi * yj * yj)) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * b * b * rho * (20 * L_PML_x * L_PML_x * beta_0_y * yj + 2 * a * a * alpha_0_x * beta_0_y * yj + 2 * a * a * alpha_0_y * beta_0_x * yj + 20 * alpha_0_x * beta_0_y * xi * xi * yj + 20 * alpha_0_y * beta_0_x * xi * xi * yj + 10 * a * alpha_0_x * beta_0_y * xi * yj + 10 * a * alpha_0_y * beta_0_x * xi * yj)) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * b * b * b * rho * (10 * L_PML_x * L_PML_x * beta_0_y + a * a * alpha_0_x * beta_0_y + a * a * alpha_0_y * beta_0_x + 10 * alpha_0_x * beta_0_y * xi * xi + 10 * alpha_0_y * beta_0_x * xi * xi + 5 * a * alpha_0_x * beta_0_y * xi + 5 * a * alpha_0_y * beta_0_x * xi)) / (150 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[13] = 0;
	C[14] = (a * a * a * b * beta_0_x * rho) / (180 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (60 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[15] = 0;
	C[16] = (a * a * a * b * beta_0_x * rho) / (120 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (120 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (36 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (36 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[17] = 0;
	C[18] = (a * a * a * b * beta_0_x * rho) / (60 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (30 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (12 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[19] = 0;
	C[20] = a / 2 + (a * alpha_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[21] =  - b / 2 - (alpha_0_y * b * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[22] = (a * a * a * b * beta_0_x * rho) / (180 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (60 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[23] = 0;
	C[24] = (a * a * a * b * beta_0_x * rho) / (90 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (90 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (9 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (9 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (9 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (9 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[25] = 0;
	C[26] = (a * a * a * b * beta_0_x * rho) / (60 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (180 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[27] = 0;
	C[28] = (a * a * a * b * beta_0_x * rho) / (120 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (120 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (36 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (36 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[29] = 0;
	C[30] =  - b / 2 - (b * (alpha_0_y * b * b + 4 * alpha_0_y * b * yj + 6 * alpha_0_y * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[31] = 0;
	C[32] =  - a / 2 - (a * alpha_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[33] = 0;
	C[34] = (a * a * a * b * beta_0_x * rho) / (180 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (60 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[35] = 0;
	C[36] = (a * a * a * b * beta_0_x * rho) / (90 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (90 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (9 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (9 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (9 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (9 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[37] = 0;
	C[38] = (a * a * a * b * beta_0_x * rho) / (60 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (180 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[39] = 0;
	C[40] = (a * a * a * b * beta_0_x * rho) / (120 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (120 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (36 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (36 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[41] = 0;
	C[42] =  - a / 2 - (a * alpha_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[43] =  - b / 2 - (b * (alpha_0_y * b * b + 4 * alpha_0_y * b * yj + 6 * alpha_0_y * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[44] = (a * a * a * b * beta_0_x * rho) / (120 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (120 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (36 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (36 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[45] = 0;
	C[46] = (a * a * a * b * beta_0_x * rho) / (60 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (180 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[47] = 0;
	C[48] = (a * a * a * b * beta_0_x * rho) / (15 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (90 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (9 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (6 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (9 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (150 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (150 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (15 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (15 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (9 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (9 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (6 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (6 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[49] = 0;
	C[50] = (a * a * a * b * beta_0_x * rho) / (30 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (60 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (12 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[51] = 0;
	C[52] = b / 2 + (b * (alpha_0_y * b * b + 4 * alpha_0_y * b * yj + 6 * alpha_0_y * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[53] = 0;
	C[54] =  - a / 2 - (a * alpha_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[55] = 0;
	C[56] = (a * a * a * b * beta_0_x * rho) / (120 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (120 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (36 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (36 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[57] = 0;
	C[58] = (a * a * a * b * beta_0_x * rho) / (60 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (180 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (600 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (180 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[59] = 0;
	C[60] = (a * a * a * b * beta_0_x * rho) / (15 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (90 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (9 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (6 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (9 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (150 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (150 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (90 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (15 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (15 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (9 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (9 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (6 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (6 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[61] = 0;
	C[62] = (a * a * a * b * beta_0_x * rho) / (30 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (60 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (12 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[63] = 0;
	C[64] =  - a / 2 - (a * alpha_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[65] = b / 2 + (b * (alpha_0_y * b * b + 4 * alpha_0_y * b * yj + 6 * alpha_0_y * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[66] = (a * a * a * b * beta_0_x * rho) / (60 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (30 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (12 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[67] = 0;
	C[68] = (a * a * a * b * beta_0_x * rho) / (120 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (120 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (36 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (36 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[69] = 0;
	C[70] = (a * a * a * b * beta_0_x * rho) / (30 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (60 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (12 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[71] = 0;
	C[72] = (a * b * rho * (60 * L_PML_y * L_PML_y * a * a * beta_0_x + 60 * L_PML_x * L_PML_x * b * b * beta_0_y + 100 * L_PML_y * L_PML_y * beta_0_x * xi * xi + 100 * L_PML_x * L_PML_x * beta_0_y * yj * yj + 150 * L_PML_y * L_PML_y * a * beta_0_x * xi + 150 * L_PML_x * L_PML_x * b * beta_0_y * yj + 36 * a * a * alpha_0_x * b * b * beta_0_y + 36 * a * a * alpha_0_y * b * b * beta_0_x + 60 * alpha_0_x * b * b * beta_0_y * xi * xi + 60 * alpha_0_y * b * b * beta_0_x * xi * xi + 60 * a * a * alpha_0_x * beta_0_y * yj * yj + 60 * a * a * alpha_0_y * beta_0_x * yj * yj + 100 * alpha_0_x * beta_0_y * xi * xi * yj * yj + 100 * alpha_0_y * beta_0_x * xi * xi * yj * yj + 90 * a * alpha_0_x * b * b * beta_0_y * xi + 90 * a * alpha_0_y * b * b * beta_0_x * xi + 90 * a * a * alpha_0_x * b * beta_0_y * yj + 90 * a * a * alpha_0_y * b * beta_0_x * yj + 150 * a * alpha_0_x * beta_0_y * xi * yj * yj + 150 * a * alpha_0_y * beta_0_x * xi * yj * yj + 150 * alpha_0_x * b * beta_0_y * xi * xi * yj + 150 * alpha_0_y * b * beta_0_x * xi * xi * yj + 225 * a * alpha_0_x * b * beta_0_y * xi * yj + 225 * a * alpha_0_y * b * beta_0_x * xi * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[73] = 0;
	C[74] = b / 2 + (alpha_0_y * b * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[75] = 0;
	C[76] = a / 2 + (a * alpha_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[77] = 0;
	C[78] = (a * a * a * b * beta_0_x * rho) / (60 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (30 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (18 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (12 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[79] = 0;
	C[80] = (a * a * a * b * beta_0_x * rho) / (120 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (120 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (36 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (36 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (36 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (36 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (400 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (120 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (36 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[81] = 0;
	C[82] = (a * a * a * b * beta_0_x * rho) / (30 * L_PML_x * L_PML_x) + (a * b * b * b * beta_0_y * rho) / (60 * L_PML_y * L_PML_y) + (a * b * beta_0_x * rho * xi * xi) / (18 * L_PML_x * L_PML_x) + (a * a * b * beta_0_x * rho * xi) / (12 * L_PML_x * L_PML_x) + (a * b * beta_0_y * rho * yj * yj) / (18 * L_PML_y * L_PML_y) + (a * b * b * beta_0_y * rho * yj) / (18 * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * b * beta_0_y * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * b * beta_0_x * rho) / (100 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * b * beta_0_y * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * b * beta_0_x * rho * xi * xi) / (60 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * b * beta_0_y * rho * xi) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * b * beta_0_x * rho * xi) / (40 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * beta_0_y * rho * yj * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_x * b * b * beta_0_y * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * beta_0_x * rho * yj * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * a * alpha_0_y * b * b * beta_0_x * rho * yj) / (30 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * beta_0_y * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_x * b * b * beta_0_y * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * beta_0_x * rho * xi * xi * yj * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * alpha_0_y * b * b * beta_0_x * rho * xi * xi * yj) / (18 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * beta_0_y * rho * xi * yj * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_x * b * b * beta_0_y * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * beta_0_x * rho * xi * yj * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y) + (a * a * alpha_0_y * b * b * beta_0_x * rho * xi * yj) / (12 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[83] = 0;
	C[84] = (a * b * rho * (60 * L_PML_y * L_PML_y * a * a * beta_0_x + 60 * L_PML_x * L_PML_x * b * b * beta_0_y + 100 * L_PML_y * L_PML_y * beta_0_x * xi * xi + 100 * L_PML_x * L_PML_x * beta_0_y * yj * yj + 150 * L_PML_y * L_PML_y * a * beta_0_x * xi + 150 * L_PML_x * L_PML_x * b * beta_0_y * yj + 36 * a * a * alpha_0_x * b * b * beta_0_y + 36 * a * a * alpha_0_y * b * b * beta_0_x + 60 * alpha_0_x * b * b * beta_0_y * xi * xi + 60 * alpha_0_y * b * b * beta_0_x * xi * xi + 60 * a * a * alpha_0_x * beta_0_y * yj * yj + 60 * a * a * alpha_0_y * beta_0_x * yj * yj + 100 * alpha_0_x * beta_0_y * xi * xi * yj * yj + 100 * alpha_0_y * beta_0_x * xi * xi * yj * yj + 90 * a * alpha_0_x * b * b * beta_0_y * xi + 90 * a * alpha_0_y * b * b * beta_0_x * xi + 90 * a * a * alpha_0_x * b * beta_0_y * yj + 90 * a * a * alpha_0_y * b * beta_0_x * yj + 150 * a * alpha_0_x * beta_0_y * xi * yj * yj + 150 * a * alpha_0_y * beta_0_x * xi * yj * yj + 150 * alpha_0_x * b * beta_0_y * xi * xi * yj + 150 * alpha_0_y * b * beta_0_x * xi * xi * yj + 225 * a * alpha_0_x * b * beta_0_y * xi * yj + 225 * a * alpha_0_y * b * beta_0_x * xi * yj)) / (900 * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[85] = 0;
	C[86] = a / 2 + (a * alpha_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[87] = b / 2 + (alpha_0_y * b * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[88] =  - b / 2 - (alpha_0_y * b * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[89] = 0;
	C[90] =  - b / 2 - (b * (alpha_0_y * b * b + 4 * alpha_0_y * b * yj + 6 * alpha_0_y * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[91] = 0;
	C[92] = b / 2 + (b * (alpha_0_y * b * b + 4 * alpha_0_y * b * yj + 6 * alpha_0_y * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[93] = 0;
	C[94] = b / 2 + (alpha_0_y * b * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[95] = 0;
	C[96] = (a * b * (nu * nu - 1) * (3 * L_PML_y * L_PML_y * a * a * beta_0_x + 3 * L_PML_x * L_PML_x * b * b * beta_0_y + 9 * L_PML_y * L_PML_y * beta_0_x * xi * xi + 9 * L_PML_x * L_PML_x * beta_0_y * yj * yj + 9 * L_PML_y * L_PML_y * a * beta_0_x * xi + 9 * L_PML_x * L_PML_x * b * beta_0_y * yj + a * a * alpha_0_x * b * b * beta_0_y + a * a * alpha_0_y * b * b * beta_0_x + 3 * alpha_0_x * b * b * beta_0_y * xi * xi + 3 * alpha_0_y * b * b * beta_0_x * xi * xi + 3 * a * a * alpha_0_x * beta_0_y * yj * yj + 3 * a * a * alpha_0_y * beta_0_x * yj * yj + 9 * alpha_0_x * beta_0_y * xi * xi * yj * yj + 9 * alpha_0_y * beta_0_x * xi * xi * yj * yj + 3 * a * alpha_0_x * b * b * beta_0_y * xi + 3 * a * alpha_0_y * b * b * beta_0_x * xi + 3 * a * a * alpha_0_x * b * beta_0_y * yj + 3 * a * a * alpha_0_y * b * beta_0_x * yj + 9 * a * alpha_0_x * beta_0_y * xi * yj * yj + 9 * a * alpha_0_y * beta_0_x * xi * yj * yj + 9 * alpha_0_x * b * beta_0_y * xi * xi * yj + 9 * alpha_0_y * b * beta_0_x * xi * xi * yj + 9 * a * alpha_0_x * b * beta_0_y * xi * yj + 9 * a * alpha_0_y * b * beta_0_x * xi * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[97] = (a * b * nu * (nu + 1) * (3 * L_PML_y * L_PML_y * a * a * beta_0_x + 3 * L_PML_x * L_PML_x * b * b * beta_0_y + 9 * L_PML_y * L_PML_y * beta_0_x * xi * xi + 9 * L_PML_x * L_PML_x * beta_0_y * yj * yj + 9 * L_PML_y * L_PML_y * a * beta_0_x * xi + 9 * L_PML_x * L_PML_x * b * beta_0_y * yj + a * a * alpha_0_x * b * b * beta_0_y + a * a * alpha_0_y * b * b * beta_0_x + 3 * alpha_0_x * b * b * beta_0_y * xi * xi + 3 * alpha_0_y * b * b * beta_0_x * xi * xi + 3 * a * a * alpha_0_x * beta_0_y * yj * yj + 3 * a * a * alpha_0_y * beta_0_x * yj * yj + 9 * alpha_0_x * beta_0_y * xi * xi * yj * yj + 9 * alpha_0_y * beta_0_x * xi * xi * yj * yj + 3 * a * alpha_0_x * b * b * beta_0_y * xi + 3 * a * alpha_0_y * b * b * beta_0_x * xi + 3 * a * a * alpha_0_x * b * beta_0_y * yj + 3 * a * a * alpha_0_y * b * beta_0_x * yj + 9 * a * alpha_0_x * beta_0_y * xi * yj * yj + 9 * a * alpha_0_y * beta_0_x * xi * yj * yj + 9 * alpha_0_x * b * beta_0_y * xi * xi * yj + 9 * alpha_0_y * b * beta_0_x * xi * xi * yj + 9 * a * alpha_0_x * b * beta_0_y * xi * yj + 9 * a * alpha_0_y * b * beta_0_x * xi * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[98] = 0;
	C[99] = 0;
	C[100] = a / 2 + (a * alpha_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[101] = 0;
	C[102] =  - a / 2 - (a * alpha_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[103] = 0;
	C[104] =  - a / 2 - (a * alpha_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[105] = 0;
	C[106] = a / 2 + (a * alpha_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[107] = (a * b * nu * (nu + 1) * (3 * L_PML_y * L_PML_y * a * a * beta_0_x + 3 * L_PML_x * L_PML_x * b * b * beta_0_y + 9 * L_PML_y * L_PML_y * beta_0_x * xi * xi + 9 * L_PML_x * L_PML_x * beta_0_y * yj * yj + 9 * L_PML_y * L_PML_y * a * beta_0_x * xi + 9 * L_PML_x * L_PML_x * b * beta_0_y * yj + a * a * alpha_0_x * b * b * beta_0_y + a * a * alpha_0_y * b * b * beta_0_x + 3 * alpha_0_x * b * b * beta_0_y * xi * xi + 3 * alpha_0_y * b * b * beta_0_x * xi * xi + 3 * a * a * alpha_0_x * beta_0_y * yj * yj + 3 * a * a * alpha_0_y * beta_0_x * yj * yj + 9 * alpha_0_x * beta_0_y * xi * xi * yj * yj + 9 * alpha_0_y * beta_0_x * xi * xi * yj * yj + 3 * a * alpha_0_x * b * b * beta_0_y * xi + 3 * a * alpha_0_y * b * b * beta_0_x * xi + 3 * a * a * alpha_0_x * b * beta_0_y * yj + 3 * a * a * alpha_0_y * b * beta_0_x * yj + 9 * a * alpha_0_x * beta_0_y * xi * yj * yj + 9 * a * alpha_0_y * beta_0_x * xi * yj * yj + 9 * alpha_0_x * b * beta_0_y * xi * xi * yj + 9 * alpha_0_y * b * beta_0_x * xi * xi * yj + 9 * a * alpha_0_x * b * beta_0_y * xi * yj + 9 * a * alpha_0_y * b * beta_0_x * xi * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[108] = (a * b * (nu * nu - 1) * (3 * L_PML_y * L_PML_y * a * a * beta_0_x + 3 * L_PML_x * L_PML_x * b * b * beta_0_y + 9 * L_PML_y * L_PML_y * beta_0_x * xi * xi + 9 * L_PML_x * L_PML_x * beta_0_y * yj * yj + 9 * L_PML_y * L_PML_y * a * beta_0_x * xi + 9 * L_PML_x * L_PML_x * b * beta_0_y * yj + a * a * alpha_0_x * b * b * beta_0_y + a * a * alpha_0_y * b * b * beta_0_x + 3 * alpha_0_x * b * b * beta_0_y * xi * xi + 3 * alpha_0_y * b * b * beta_0_x * xi * xi + 3 * a * a * alpha_0_x * beta_0_y * yj * yj + 3 * a * a * alpha_0_y * beta_0_x * yj * yj + 9 * alpha_0_x * beta_0_y * xi * xi * yj * yj + 9 * alpha_0_y * beta_0_x * xi * xi * yj * yj + 3 * a * alpha_0_x * b * b * beta_0_y * xi + 3 * a * alpha_0_y * b * b * beta_0_x * xi + 3 * a * a * alpha_0_x * b * beta_0_y * yj + 3 * a * a * alpha_0_y * b * beta_0_x * yj + 9 * a * alpha_0_x * beta_0_y * xi * yj * yj + 9 * a * alpha_0_y * beta_0_x * xi * yj * yj + 9 * alpha_0_x * b * beta_0_y * xi * xi * yj + 9 * alpha_0_y * b * beta_0_x * xi * xi * yj + 9 * a * alpha_0_x * b * beta_0_y * xi * yj + 9 * a * alpha_0_y * b * beta_0_x * xi * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);
	C[109] = 0;
	C[110] = a / 2 + (a * alpha_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[111] =  - b / 2 - (alpha_0_y * b * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[112] =  - a / 2 - (a * alpha_0_x * (a * a + 4 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[113] =  - b / 2 - (b * (alpha_0_y * b * b + 4 * alpha_0_y * b * yj + 6 * alpha_0_y * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[114] =  - a / 2 - (a * alpha_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[115] = b / 2 + (b * (alpha_0_y * b * b + 4 * alpha_0_y * b * yj + 6 * alpha_0_y * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[116] = a / 2 + (a * alpha_0_x * (3 * a * a + 8 * a * xi + 6 * xi * xi)) / (12 * L_PML_x * L_PML_x);
	C[117] = b / 2 + (alpha_0_y * b * (3 * b * b + 8 * b * yj + 6 * yj * yj)) / (12 * L_PML_y * L_PML_y);
	C[118] = 0;
	C[119] = 0;
	C[120] =  - (2 * a * b * (nu + 1) * (3 * L_PML_y * L_PML_y * a * a * beta_0_x + 3 * L_PML_x * L_PML_x * b * b * beta_0_y + 9 * L_PML_y * L_PML_y * beta_0_x * xi * xi + 9 * L_PML_x * L_PML_x * beta_0_y * yj * yj + 9 * L_PML_y * L_PML_y * a * beta_0_x * xi + 9 * L_PML_x * L_PML_x * b * beta_0_y * yj + a * a * alpha_0_x * b * b * beta_0_y + a * a * alpha_0_y * b * b * beta_0_x + 3 * alpha_0_x * b * b * beta_0_y * xi * xi + 3 * alpha_0_y * b * b * beta_0_x * xi * xi + 3 * a * a * alpha_0_x * beta_0_y * yj * yj + 3 * a * a * alpha_0_y * beta_0_x * yj * yj + 9 * alpha_0_x * beta_0_y * xi * xi * yj * yj + 9 * alpha_0_y * beta_0_x * xi * xi * yj * yj + 3 * a * alpha_0_x * b * b * beta_0_y * xi + 3 * a * alpha_0_y * b * b * beta_0_x * xi + 3 * a * a * alpha_0_x * b * beta_0_y * yj + 3 * a * a * alpha_0_y * b * beta_0_x * yj + 9 * a * alpha_0_x * beta_0_y * xi * yj * yj + 9 * a * alpha_0_y * beta_0_x * xi * yj * yj + 9 * alpha_0_x * b * beta_0_y * xi * xi * yj + 9 * alpha_0_y * b * beta_0_x * xi * xi * yj + 9 * a * alpha_0_x * b * beta_0_y * xi * yj + 9 * a * alpha_0_y * b * beta_0_x * xi * yj)) / (9 * E * L_PML_x * L_PML_x * L_PML_y * L_PML_y);

}
