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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/FourNodeQuad.cpp,v $
                                                                        
                                                                        
// File: ~/element/FourNodeQuad.C
//
// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class definition for FourNodeQuad.
//
// What: "@(#) FourNodeQuad.C, revA"



#include <FourNodeQuad.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <GaussQuadRule1d.h>

#include <G3Globals.h>

Matrix FourNodeQuad::N (2,8);

FourNodeQuad::FourNodeQuad (int tag, int nd1, int nd2, int nd3, int nd4,
	NDMaterial &m, const char *type, double t,
	double p, double r, double b1, double b2)
:Element (tag, ELE_TAG_FourNodeQuad), thickness(t), rho(r),
 M(8,8), K(8,8), C(8,8), P(8), Q(8), b(2), pressure(p), pressureLoad(8),
 B(3,8), J(2,2), L(2,2), connectedExternalNodes(4)
{
	b(0) = b1;
	b(1) = b2;

    // Change this later so object comes through constructor
    theQuadRule = new GaussQuadRule1d (2);
    // theQuadRule = rule.getCopy(); ???
	order = theQuadRule->getOrder();
    
    if (theQuadRule == 0) {
	g3ErrorHandler->fatal("FATAL ERROR FourNodeQuad - failed to get a copy of quadrature rule");
    }
    
    // Get the material model
    NDMaterial *theModel = m.getCopy (type);
    
    if (theModel == 0){
	g3ErrorHandler->fatal("FATAL ERROR FourNodeQuad - failed acquire material model prototype");
    }

    int i, j;
    
    // Allocate arrays of pointers to NDMaterials
    theMaterial = new NDMaterial * *[order];
    
    for (i = 0; i < order; i++){
	// Allocate NDMaterial pointers for each array
	theMaterial[i] = new NDMaterial *[order];

	// Check allocation
	if (theMaterial[i] == 0){
	    g3ErrorHandler->fatal("FATAL ERROR FourNodeQuad - failed allocate material model pointers");
	}

	// Get copies of the material model for each integration point
	for (j = 0; j < order; j++){
	    theMaterial[i][j] = theModel->getCopy();
	    
	    // Check allocation
	    if (theMaterial[i][j] == 0) {
		g3ErrorHandler->fatal("FATAL ERROR FourNodeQuad - failed to get a copy of material model");
	    }
	}
    }

    // Set connected external node IDs
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;
    connectedExternalNodes(2) = nd3;
    connectedExternalNodes(3) = nd4;
}

FourNodeQuad::FourNodeQuad ():Element (0,ELE_TAG_FourNodeQuad),
thickness(0.0), rho(0.0), order(0),
M(8,8), K(8,8), C(8,8), P(8), Q(8), b(2), pressure(0.0),
B(3,8), J(2,2), L(2,2), connectedExternalNodes(4), theMaterial(0), 
pressureLoad(8), theQuadRule(0)
{

}

FourNodeQuad::~FourNodeQuad ()
{    
    for (int i = 0; i < order; i++) {
	for (int j = 0; j < order; j++)
	    // Delete the NDMaterials at each integration point
	    if (theMaterial[i][j])
		delete theMaterial[i][j];
	
	// Delete each array of NDMaterial pointers
	if (theMaterial[i])
	    delete [] theMaterial[i];
    }	

    // Delete the array of pointers to NDMaterial pointer arrays
    if (theMaterial)
	delete [] theMaterial;

    // Delete the quadrature rule
    if (theQuadRule)
	delete theQuadRule;
}

int
FourNodeQuad::getNumExternalNodes () const
{
    return 4;
}

const ID&
FourNodeQuad::getExternalNodes ()
{
    return connectedExternalNodes;
}

int
FourNodeQuad::getNumDOF ()
{
    return 8;
}

void
FourNodeQuad::setDomain (Domain *theDomain)
{
	// Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
	nd1Ptr = 0;
	nd2Ptr = 0;
	nd3Ptr = 0;
	nd4Ptr = 0;
	return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    int Nd3 = connectedExternalNodes(2);
    int Nd4 = connectedExternalNodes(3);

    nd1Ptr = theDomain->getNode(Nd1);
    nd2Ptr = theDomain->getNode(Nd2);
    nd3Ptr = theDomain->getNode(Nd3);
    nd4Ptr = theDomain->getNode(Nd4);

    if (nd1Ptr == 0 || nd2Ptr == 0 || nd3Ptr == 0 || nd4Ptr == 0) {
	//g3ErrorHandler->fatal("FATAL ERROR FourNodeQuad (tag: %d), node not found in domain",
	//	this->getTag());
	
	return;
    }

    int dofNd1 = nd1Ptr->getNumberDOF();
    int dofNd2 = nd2Ptr->getNumberDOF();
    int dofNd3 = nd3Ptr->getNumberDOF();
    int dofNd4 = nd4Ptr->getNumberDOF();
    
    if (dofNd1 != 2 || dofNd2 != 2 || dofNd3 != 2 || dofNd4 != 2) {
	//g3ErrorHandler->fatal("FATAL ERROR FourNodeQuad (tag: %d), has differing number of DOFs at its nodes",
	//	this->getTag());
	
	return;
    }
    this->DomainComponent::setDomain(theDomain);

	// Set mass matrix ... it won't change
	this->getMass();

	// Compute consistent nodal loads due to pressure
	this->setPressureLoadAtNodes();
}

int
FourNodeQuad::commitState ()
{
    int i, j;
    int retVal = 0;

    // Loop over the integration points and commit the material states
    for (i = 0; i < order; i++)
	for (j = 0; j < order; j++)
	    retVal += (theMaterial[i][j])->commitState();

    return retVal;
}

int
FourNodeQuad::revertToLastCommit ()
{
    int i, j;
    int retVal = 0;

    // Loop over the integration points and revert to last committed material states
    for (i = 0; i < order; i++)
	for (j = 0; j < order; j++)
	    retVal += (theMaterial[i][j])->revertToLastCommit();

    return retVal;
}

int
FourNodeQuad::revertToStart ()
{
    int i, j;
    int retVal = 0;

    // Loop over the integration points and revert to initial material states
    for (i = 0; i < order; i++)
	for (j = 0; j < order; j++)
	    retVal += (theMaterial[i][j])->revertToStart();

    return retVal;
}

const Matrix&
FourNodeQuad::getTangentStiff ()
{
	const Matrix &intPt = theQuadRule->getIntegrPointCoords();
	const Vector &intWt = theQuadRule->getIntegrPointWeights();

	const Vector &disp1 = nd1Ptr->getTrialDisp();
	const Vector &disp2 = nd2Ptr->getTrialDisp();
	const Vector &disp3 = nd3Ptr->getTrialDisp();
	const Vector &disp4 = nd4Ptr->getTrialDisp();
	
	static Vector u(8);

	u(0) = disp1(0);
	u(1) = disp1(1);
	u(2) = disp2(0);
	u(3) = disp2(1);
	u(4) = disp3(0);
	u(5) = disp3(1);
	u(6) = disp4(0);
	u(7) = disp4(1);

	static Vector eps (3);

	K.Zero();

	// Loop over the integration points
	for (int i = 0; i < order; i++)
	{
		for (int j = 0; j < order; j++)
		{
			// Determine Jacobian for this integration point
			this->setJacobian (intPt(i,0), intPt(j,0));

			// Interpolate strains
			this->formBMatrix (intPt(i,0), intPt(j,0));
			//eps = B*u;
			eps.addMatrixVector(0.0, B, u, 1.0);

			// Set the material strain
			(theMaterial[i][j])->setTrialStrain (eps);

			// Get the material tangent
			const Matrix &D = (theMaterial[i][j])->getTangent();

			// Form the Jacobian of the coordinate transformation
			double detJ = this->formDetJ (intPt(i,0), intPt(j,0));

			// Perform numerical integration
			//K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
			K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
		}
	}

	K = K * thickness;

	return K;
}

const Matrix&
FourNodeQuad::getSecantStiff ()
{
	return K;
}

const Matrix&
FourNodeQuad::getDamp ()
{
	return C;
}

const Matrix&
FourNodeQuad::getMass ()
{
	M.Zero();

	if (rho == 0.0)
		return M;

	const Matrix &intPt = theQuadRule->getIntegrPointCoords();
	const Vector &intWt = theQuadRule->getIntegrPointWeights();

	int i, j;

	// Loop over the integration points
	for (i = 0; i < order; i++)
	{
		for (j = 0; j < order; j++)
		{
			// Determine Jacobian for this integration point
			this->setJacobian (intPt(i,0), intPt(j,0));

			// Interpolate strains
			this->formNMatrix (intPt(i,0), intPt(j,0));

			// Form the Jacobian of the coordinate transformation
			double detJ = this->formDetJ (intPt(i,0), intPt(j,0));

			// Perform numerical integration
			//M = M + (N^ N) * intWt(i)*intWt(j) * detJ;
			M.addMatrix(1.0, N^ N, intWt(i)*intWt(j)*detJ);
		}
	}

	M = M * (thickness * rho);

	// Lumped mass ... can be optional
	for (j = 0; j < 8; j++) {
		double sum = 0.0;
		// Lump each column onto its diagonal
		for (i = 0; i < 8; i++) {
			sum += M(i,j);
			M(i,j) = 0.0;
		}
		M(j,j) = sum;
	}

	return M;
}

void
FourNodeQuad::zeroLoad(void)
{
	Q.Zero();

	return;
}

int 
FourNodeQuad::addLoad(const Vector &addLoad)
{
	if (addLoad.Size() != 8) {
		g3ErrorHandler->warning("FourNodeQuad::addLoad %s\n",
				"Vector not of correct size");
		return -1;
	}

	// Add to the external nodal loads
	//Q += addLoad;
	Q.addVector(1.0, addLoad, 1.0);

	return 0;
}

int 
FourNodeQuad::addInertiaLoadToUnbalance(const Vector &accel)
{
	// Check for a quick return
	if (rho == 0.0) 
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = nd1Ptr->getRV(accel);
	const Vector &Raccel2 = nd2Ptr->getRV(accel);
	const Vector &Raccel3 = nd3Ptr->getRV(accel);
	const Vector &Raccel4 = nd4Ptr->getRV(accel);

    if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size() ||
		2 != Raccel4.Size()) {
		g3ErrorHandler->warning("FourNodeQuad::addInertiaLoadToUnbalance %s\n",
				"matrix and vector sizes are incompatable");
		return -1;
    }

	static Vector ra(8);

	ra(0) = Raccel1(0);
	ra(1) = Raccel1(1);
	ra(2) = Raccel2(0);
	ra(3) = Raccel2(1);
	ra(4) = Raccel3(0);
	ra(5) = Raccel3(1);
	ra(6) = Raccel4(0);
	ra(7) = Raccel4(1);
    
    // Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	// Mass matrix is computed in setDomain()
    for (int i = 0; i < 8; i++)
		Q(i) += -M(i,i)*ra(i);

    return 0;
}

const Vector&
FourNodeQuad::getResistingForce ()
{
	const Matrix &intPt = theQuadRule->getIntegrPointCoords();
	const Vector &intWt = theQuadRule->getIntegrPointWeights();
	
	const Vector &disp1 = nd1Ptr->getTrialDisp();
	const Vector &disp2 = nd2Ptr->getTrialDisp();
	const Vector &disp3 = nd3Ptr->getTrialDisp();
	const Vector &disp4 = nd4Ptr->getTrialDisp();

	static Vector u(8);

	u(0) = disp1(0);
	u(1) = disp1(1);
	u(2) = disp2(0);
	u(3) = disp2(1);
	u(4) = disp3(0);
	u(5) = disp3(1);
	u(6) = disp4(0);
	u(7) = disp4(1);

	static Vector eps (3);

	P.Zero();

	// Loop over the integration points
	for (int i = 0; i < order; i++)
	{
		for (int j = 0; j < order; j++)
		{
			// Determine Jacobian for this integration point
			this->setJacobian (intPt(i,0), intPt(j,0));

			// Interpolate strains
			this->formBMatrix (intPt(i,0), intPt(j,0));
			//eps = B*u;
			eps.addMatrixVector(0.0, B, u, 1.0);

			// Set the material strain
			(theMaterial[i][j])->setTrialStrain (eps);

			// Get material stress response
			const Vector &sigma = (theMaterial[i][j])->getStress();

			// Form the Jacobian of the coordinate transformation
			double detJ = this->formDetJ (intPt(i,0), intPt(j,0));

			// Perform numerical integration on internal force
			//P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
			P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);

			// Form displacement interpolation for equiv. body forces
			this->formNMatrix(intPt(i,0), intPt(j,0));

			// Subtract equiv. body forces from the nodes
			//P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
			P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
		}
	}

	P = P * thickness;

	// Subtract pressure loading from resisting force
	if (pressure != 0.0) {
		//P = P - pressureLoad;
		P.addVector(1.0, pressureLoad, -1.0);
	}
	
	// Subtract other external nodal loads ... P_res = P_int - P_ext
	//P = P - Q;
	P.addVector(1.0, Q, -1.0);

	return P;
}

const Vector&
FourNodeQuad::getResistingForceIncInertia ()
{
	// Check for a quick return
	if (rho == 0.0)
		return this->getResistingForce();

	const Vector &accel1 = nd1Ptr->getTrialAccel();
	const Vector &accel2 = nd2Ptr->getTrialAccel();
	const Vector &accel3 = nd3Ptr->getTrialAccel();
	const Vector &accel4 = nd4Ptr->getTrialAccel();
	
	static Vector a(8);

	a(0) = accel1(0);
	a(1) = accel1(1);
	a(2) = accel2(0);
	a(3) = accel2(1);
	a(4) = accel3(0);
	a(5) = accel3(1);
	a(6) = accel4(0);
	a(7) = accel4(1);

	// Compute the current resisting force
	this->getResistingForce();

	// Take advantage of lumped mass matrix
	// Mass matrix is computed in setDomain()
	for (int i = 0; i < 8; i++)
		P(i) += M(i,i)*a(i);

	return P;
}

int
FourNodeQuad::sendSelf (int commitTag, Channel &theChannel)
{
	int res = 0;

	// note: we don't check for dataTag == 0 for Element
	// objects as that is taken care of in a commit by the Domain
	// object - don't want to have to do the check if sending data
	int dataTag = this->getDbTag();

	// Quad packs its data into a Vector and sends this to theChannel
	// along with its dbTag and the commitTag passed in the arguments
	static Vector data(7);
	data(0) = this->getTag();
	data(1) = thickness;
	data(2) = rho;
	data(3) = b(0);
	data(4) = b(1);
	data(5) = pressure;
	data(6) = order;

	res += theChannel.sendVector(dataTag, commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("WARNING FourNodeQuad::sendSelf() - %d failed to send Vector\n",this->getTag());
		return res;
	}	      

	// Quad then sends the tags of its four end nodes
	res += theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		g3ErrorHandler->warning("WARNING FourNodeQuad::sendSelf() - %d failed to send ID\n",this->getTag());
		return res;
	}

	// Now quad sends the ids of its materials
	int matDbTag;
	int numMats = order*order;
	ID classTags(2*numMats);

	int i,j;
	for (i = 0; i < order; i++) {
		for (j = 0; j < order; j++) {
			int k = i*order+j;
			classTags(k) = theMaterial[i][j]->getClassTag();
			matDbTag = theMaterial[i][j]->getDbTag();
			// NOTE: we do have to ensure that the material has a database
			// tag if we are sending to a database channel.
			if (matDbTag == 0) {
				matDbTag = theChannel.getDbTag();
				if (matDbTag != 0)
					theMaterial[i][j]->setDbTag(matDbTag);
			}
			classTags(k+numMats) = matDbTag;
		}
	}

	res += theChannel.sendID(dataTag, commitTag, classTags);
	if (res < 0) {
		g3ErrorHandler->warning("WARNING FourNodeQuad::sendSelf() - %d failed to send ID\n",
			this->getTag());
		return res;
	}

	// Finally, quad asks its material objects to send themselves
	for (i = 0; i < order; i++) {
		for (j = 0; j < order; j++) {
			res += theMaterial[i][j]->sendSelf(commitTag, theChannel);
			if (res < 0) {
				g3ErrorHandler->warning("WARNING FourNodeQuad::sendSelf() - %d failed to send its Material\n",this->getTag());
				return res;
			}
		}
	}

	return res;
}

int
FourNodeQuad::recvSelf (int commitTag, Channel &theChannel,
						FEM_ObjectBroker &theBroker)
{
	int res = 0;

	int dataTag = this->getDbTag();

	// Quad creates a Vector, receives the Vector and then sets the 
	// internal data with the data in the Vector
	static Vector data(7);
	res += theChannel.recvVector(dataTag, commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("WARNING FourNodeQuad::recvSelf() - failed to receive Vector\n");
		return res;
	}

	this->setTag((int)data(0));
	thickness = data(1);
	rho = data(2);
	b(0) = data(3);
	b(1) = data(4);
	pressure = data(5);
  
	// Quad now receives the tags of its four external nodes
	res += theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
	if (res < 0) {
		g3ErrorHandler->warning("WARNING FourNodeQuad::recvSelf() - %d failed to receive ID\n", this->getTag());
		return res;
	}

	// Quad now receives the ids of its materials
	int newOrder = (int)data(6);
	int numMats = newOrder*newOrder;
	ID classTags(2*numMats);

	res += theChannel.recvID(dataTag, commitTag, classTags);
	if (res < 0)  {
		g3ErrorHandler->warning("FourNodeQuad::recvSelf() - %s\n",
			    "failed to recv ID data");
		return res;
	}    

	int i,j,k;

	// If the number of materials (quadrature order) is not the same,
	// delete the old materials, allocate new ones and then receive
	if (order != newOrder) {
		// Delete the materials
		for (i = 0; i < order; i++) {
			for (j = 0; j < order; j++) {
				if (theMaterial[i][j])
					delete theMaterial[i][j];
			}
			if (theMaterial[i])
				delete [] theMaterial[i];
		}
		if (theMaterial)
			delete [] theMaterial;

		// Allocate new materials
		order = newOrder;
		theMaterial = new NDMaterial * *[order];
		if (theMaterial == 0) {
			g3ErrorHandler->warning("FourNodeQuad::recvSelf() - %s\n",
				"Could not allocate NDMaterial** pointer");
			return -1;
		}
		for (i = 0; i < order; i++) {
			theMaterial[i] = new NDMaterial *[order];
			if (theMaterial[i] == 0) {
				g3ErrorHandler->warning("FourNodeQuad::recvSelf() - %s\n",
				"Could not allocate NDMaterial* pointer");
				return -1;
			}
			for (j = 0; j < order; j++) {
				k = i*order + j;
				int matClassTag = classTags(k);
				int matDbTag = classTags(k+numMats);
				// Allocate new material with the sent class tag
				theMaterial[i][j] = theBroker.getNewNDMaterial(matClassTag);
				if (theMaterial[i][j] == 0) {
					g3ErrorHandler->warning("FourNodeQuad::recvSelf() - %s %d\n",
						"Broker could not create NDMaterial of class type",matClassTag);
					return -1;
				}
				// Now receive materials into the newly allocated space
				theMaterial[i][j]->setDbTag(matDbTag);
				res += theMaterial[i][j]->recvSelf(commitTag, theChannel, theBroker);
				if (res < 0) {
					g3ErrorHandler->warning("NLBeamColumn3d::recvSelf() - material %d,%d %s\n",
						i,j,"failed to recv itself");
					return res;
				}
			}
		}
	}
	// Number of materials is the same, receive materials into current space
	else {
		order = (int)data(6);
		for (i = 0; i < order; i++) {
			for (j = 0; j < order; j++) {
				k = i*order + j;
				int matClassTag = classTags(k);
				int matDbTag = classTags(k+numMats);
				// Check that material is of the right type; if not,
				// delete it and create a new one of the right type
				if (theMaterial[i][j]->getClassTag() != matClassTag) {
					delete theMaterial[i][j];
					theMaterial[i][j] = theBroker.getNewNDMaterial(matClassTag);
					if (theMaterial[i][j] == 0) {
						g3ErrorHandler->fatal("FourNodeQuad::recvSelf() - %s %d\n",
							"Broker could not create NDMaterial of class type",matClassTag);
						return -1;
					}
				}
				// Receive the material
				theMaterial[i][j]->setDbTag(matDbTag);
				res += theMaterial[i][j]->recvSelf(commitTag, theChannel, theBroker);
				if (res < 0) {
					g3ErrorHandler->warning("FourNodeQuad::recvSelf() - material %d,%d %s\n",
						i,j,"failed to recv itself");
					return res;
				}
			}
		}
	}

	return res;
}

void
FourNodeQuad::Print (ostream &s, int flag)
{
	s << "\nFourNodeQuad, element id:  " << this->getTag() << endl;
	s << "\tConnected external nodes:  " << connectedExternalNodes;
	s << "\tthickness:  " << thickness << endl;
	s << "\tmass density:  " << rho << endl;
	s << "\tsurface pressure:  " << pressure << endl;
	s << "\tMaterial: " << theMaterial[0][0]->getType() << endl;
}

int
FourNodeQuad::displaySelf (Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = nd1Ptr->getCrds();
    const Vector &end2Crd = nd2Ptr->getCrds();	
	const Vector &end3Crd = nd3Ptr->getCrds();	
	const Vector &end4Crd = nd4Ptr->getCrds();	

    const Vector &end1Disp = nd1Ptr->getDisp();
    const Vector &end2Disp = nd2Ptr->getDisp();
    const Vector &end3Disp = nd3Ptr->getDisp();
	const Vector &end4Disp = nd4Ptr->getDisp();

	static Vector v1(3);
	static Vector v2(3);
	static Vector v3(3);
	static Vector v4(3);

	for (int i = 0; i < 2; i++)
	{
		v1(i) = end1Crd(i) + end1Disp(i)*fact;
		v2(i) = end2Crd(i) + end2Disp(i)*fact;    
		v3(i) = end3Crd(i) + end3Disp(i)*fact;    
		v4(i) = end4Crd(i) + end4Disp(i)*fact;    
	}
	
	int error = 0;

	error += theViewer.drawLine (v1, v2, 1.0, 1.0);
	error += theViewer.drawLine (v2, v3, 1.0, 1.0);
	error += theViewer.drawLine (v3, v4, 1.0, 1.0);
	error += theViewer.drawLine (v4, v1, 1.0, 1.0);

	return error;
}

int 
FourNodeQuad::setResponse (char **argv, int argc, Information &eleInformation)
{
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
		Vector *newVector = new Vector(8);
		if (newVector == 0) {
		g3ErrorHandler->warning("WARNING FourNodeQuad::setResponse() - %d out of memory creating vector\n",
				  this->getTag());
		return -1;
		}
		eleInformation.theVector = newVector;
		eleInformation.theType = VectorType;
		return 1;
    } 
    
    else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0) {
		Matrix *newMatrix = new Matrix(8,8);
		if (newMatrix == 0) {
		g3ErrorHandler->warning("WARNING FourNodeQuad::setResponse() - %d out of memory creating matrix\n",
				  this->getTag());
		return -1;
		}
		eleInformation.theMatrix = newMatrix;
		eleInformation.theType = MatrixType;
		return 2;
    } 

	else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {
		int pointNum = atoi(argv[1]);
		int order = theQuadRule->getOrder();
		if (pointNum > 0 && pointNum <= order*order) {
			int i,j;
			this->getMaterialIndices(pointNum,i,j);
			int ok = theMaterial[i][j]->setResponse(&argv[2], argc-2, eleInformation);
			if (ok < 0)
				return -1;
		    else if (ok >= 0 && ok < 100)
				return pointNum*100 + ok;
		}
	    else 
			return -1;
	}
 
    // otherwise response quantity is unknown for the quad class
    else
		return -1;
}

int 
FourNodeQuad::getResponse (int responseID, Information &eleInformation)
{
	switch (responseID) {
		case -1:
			return -1;
      
		case 1:
			if (eleInformation.theVector != 0)
				*(eleInformation.theVector) = P;    
			return 0;
      
		case 2:
			if (eleInformation.theMatrix != 0)
				*(eleInformation.theMatrix) = K;
			return 0;      
		default: 
			if (responseID >= 100) { // material quantity
				int pointNum = responseID/100;
				int order = theQuadRule->getOrder();
				if (pointNum > 0 && pointNum <= order*order) {
					int i,j;
					this->getMaterialIndices(pointNum,i,j);
					return theMaterial[i][j]->getResponse(responseID-100*pointNum, eleInformation);
				}
				else
					return -1;
			} else // unknown
				return -1;
	}
}

int
FourNodeQuad::setParameter(char **argv, int argc, Information &info)
{
	// quad mass density per unit volume
	if (strcmp(argv[0],"rho") == 0) {
		info.theType = DoubleType;
		info.theDouble = rho;
		return 1;
	}
	// quad pressure loading
	if (strcmp(argv[0],"pressure") == 0) {
		info.theType = DoubleType;
		info.theDouble = pressure;
		return 2;
	}
    // a material parameter
    else if (strcmp(argv[0],"material") == 0) {
		int pointNum = atoi(argv[1]);
		int order = theQuadRule->getOrder();
		if (pointNum > 0 && pointNum <= order*order) {
			int i,j;
			this->getMaterialIndices(pointNum,i,j);
			int ok = theMaterial[i][j]->setParameter(&argv[2], argc-2, info);
			if (ok < 0)
				return -1;
		    else if (ok >= 0 && ok < 100)
				return pointNum*100 + ok;
		}
	    else 
			return -1;
	}
    
    // otherwise parameter is unknown for the Truss class
    else
		return -1;

}
    
int
FourNodeQuad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
    case -1:
      return -1;
      
	case 1:
		rho = info.theDouble;
		this->getMass();	// update mass matrix
		return 0;
	case 2:
		pressure = info.theDouble;
		this->setPressureLoadAtNodes();	// update consistent nodal loads
		return 0;
	default: 
		if (parameterID >= 100) { // material parameter
			int pointNum = parameterID/100;
			int order = theQuadRule->getOrder();
			if (pointNum > 0 && pointNum <= order*order) {
				int i,j;
				this->getMaterialIndices(pointNum,i,j);
				return theMaterial[i][j]->updateParameter(parameterID-100*pointNum, info);
			}
			else
				return -1;
		} else // unknown
			return -1;
  }
}

void
FourNodeQuad::getMaterialIndices(int pointNum, int &i, int &j)
{
	int order = theQuadRule->getOrder();

	if (pointNum < 1 || pointNum > order*order) {
		i = 0; j = 0;
		return;
	}

	// Set integration point indices corresponding to nodes in isoparametric domain
	if (order == 1) {
		i = 0; j = 0; 
		return;
	} else if (order == 2) {
		switch (pointNum) {
			case 1:
				i = 0; j = 0; return;
			case 2:
				i = 1; j = 0; return;
			case 3:
				i = 1; j = 1; return;
			case 4:
				i = 0; j = 1; return;
			default:
				i = 0; j = 0; return;
		}
	} else if (order == 3) {
		switch (pointNum) {
			case 1:
				i = 0; j = 0; return;
			case 2:
				i = 2; j = 0; return;
			case 3:
				i = 2; j = 2; return;
			case 4:
				i = 0; j = 2; return;
			case 5:
				i = 1; j = 0; return;
			case 6:
				i = 2; j = 1; return;
			case 7:
				i = 1; j = 2; return;
			case 8:
				i = 0; j = 1; return;
			case 9:
				i = 1; j = 1; return;
			default:
				i = 0; j = 0; return;
		}
	
	} else {	// Add more cases for order 4 or higher,
		i = 0; j = 0;
		return;
	}
}

void
FourNodeQuad::formNMatrix (double xi, double eta)
{
	N.Zero();

	N(0,0) = N(1,1) = 0.25*(1.0-xi)*(1.0-eta);		// N_1
	N(0,2) = N(1,3) = 0.25*(1.0+xi)*(1.0-eta);		// N_2
	N(0,4) = N(1,5) = 0.25*(1.0+xi)*(1.0+eta);		// N_3
	N(0,6) = N(1,7) = 0.25*(1.0-xi)*(1.0+eta);		// N_4
}

void
FourNodeQuad::setJacobian (double xi, double eta)
{
	const Vector &nd1Crds = nd1Ptr->getCrds();
	const Vector &nd2Crds = nd2Ptr->getCrds();
	const Vector &nd3Crds = nd3Ptr->getCrds();
	const Vector &nd4Crds = nd4Ptr->getCrds();

	J(0,0) = -nd1Crds(0)*(1.0-eta) + nd2Crds(0)*(1.0-eta) +
				nd3Crds(0)*(1.0+eta) - nd4Crds(0)*(1.0+eta);

	J(0,1) = -nd1Crds(0)*(1.0-xi) - nd2Crds(0)*(1.0+xi) +
				nd3Crds(0)*(1.0+xi) + nd4Crds(0)*(1.0-xi);

	J(1,0) = -nd1Crds(1)*(1.0-eta) + nd2Crds(1)*(1.0-eta) +
				nd3Crds(1)*(1.0+eta) - nd4Crds(1)*(1.0+eta);

	J(1,1) = -nd1Crds(1)*(1.0-xi) - nd2Crds(1)*(1.0+xi) +
				nd3Crds(1)*(1.0+xi) + nd4Crds(1)*(1.0-xi);

	J = J * 0.25;

	// L = inv(J)
	L(0,0) = J(1,1);
	L(1,0) = -J(0,1);
	L(0,1) = -J(1,0);
	L(1,1) = J(0,0);

	L = L / formDetJ (xi, eta);
}

void
FourNodeQuad::formBMatrix (double xi, double eta)
{
    B.Zero();

    double L00 = L(0,0);
    double L10 = L(1,0);
    double L01 = L(0,1);
    double L11 = L(1,1);

    // See Cook, Malkus, Plesha p. 169 for the derivation of these terms
    B(0,0) = L00*-0.25*(1.0-eta) + L01*-0.25*(1.0-xi);		// N_1,1
    B(0,2) = L00*0.25*(1.0-eta) + L01*-0.25*(1.0+xi);		// N_2,1
    B(0,4) = L00*0.25*(1.0+eta) + L01*0.25*(1.0+xi);		// N_3,1
    B(0,6) = L00*-0.25*(1.0+eta) + L01*0.25*(1.0-xi);		// N_4,1
	
    B(1,1) = L10*-0.25*(1.0-eta) + L11*-0.25*(1.0-xi);	// N_1,2
    B(1,3) = L10*0.25*(1.0-eta) + L11*-0.25*(1.0+xi);		// N_2,2
    B(1,5) = L10*0.25*(1.0+eta) + L11*0.25*(1.0+xi);		// N_3,2
    B(1,7) = L10*-0.25*(1.0+eta) + L11*0.25*(1.0-xi);		// N_4,2

    B(2,0) = B(1,1);
    B(2,1) = B(0,0);
    B(2,2) = B(1,3);
    B(2,3) = B(0,2);
    B(2,4) = B(1,5);
    B(2,5) = B(0,4);
    B(2,6) = B(1,7);
    B(2,7) = B(0,6);
}

double
FourNodeQuad::formDetJ (double xi, double eta)
{
    return J(0,0)*J(1,1) - J(0,1)*J(1,0);
}

void 
FourNodeQuad::setPressureLoadAtNodes(void)
{
	pressureLoad.Zero();

	if (pressure == 0.0)
		return;

	const Vector &node1 = nd1Ptr->getCrds();
	const Vector &node2 = nd2Ptr->getCrds();
	const Vector &node3 = nd3Ptr->getCrds();
	const Vector &node4 = nd4Ptr->getCrds();

	double x1 = node1(0);
	double y1 = node1(1);
	double x2 = node2(0);
	double y2 = node2(1);
	double x3 = node3(0);
	double y3 = node3(1);
	double x4 = node4(0);
	double y4 = node4(1);

	double dx12 = x2-x1;
	double dy12 = y2-y1;
	double dx23 = x3-x2;
	double dy23 = y3-y2;
	double dx34 = x4-x3;
	double dy34 = y4-y3;
	double dx41 = x1-x4;
	double dy41 = y1-y4;

	double pressureOver2 = pressure/2;

	// Contribution from side 12
	pressureLoad(0) += pressureOver2*dy12;
	pressureLoad(2) += pressureOver2*dy12;
	pressureLoad(1) += pressureOver2*-dx12;
	pressureLoad(3) += pressureOver2*-dx12;

	// Contribution from side 23
	pressureLoad(2) += pressureOver2*dy23;
	pressureLoad(4) += pressureOver2*dy23;
	pressureLoad(3) += pressureOver2*-dx23;
	pressureLoad(5) += pressureOver2*-dx23;

	// Contribution from side 34
	pressureLoad(4) += pressureOver2*dy34;
	pressureLoad(6) += pressureOver2*dy34;
	pressureLoad(5) += pressureOver2*-dx34;
	pressureLoad(7) += pressureOver2*-dx34;

	// Contribution from side 41
	pressureLoad(6) += pressureOver2*dy41;
	pressureLoad(0) += pressureOver2*dy41;
	pressureLoad(7) += pressureOver2*-dx41;
	pressureLoad(1) += pressureOver2*-dx41;

	//pressureLoad = pressureLoad*thickness;
}

