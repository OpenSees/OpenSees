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
                                                                        
// $Revision: 1.2 $
// $Date: 2000-12-18 10:38:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beam2d/ElasticBeam2d.cpp,v $
                                                                        
                                                                        
// File: ~/model/ElasticBeam2d.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ElasticBeam2d.
// ElasticBeam2d is a 3d beam element. As such it can only
// connect to a node with 6-dof. 

#include <ElasticBeam2d.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf2d.h>
#include <Information.h>
#include <ElementResponse.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>

ElasticBeam2d::ElasticBeam2d()
:Element(0,ELE_TAG_ElasticBeam2d), 
 A(0), E(0), I(0), L(0.0), rho(0.0),
 connectedExternalNodes(2), theCoordTransf(0),
 m(6,6), d(6,6), Pinert(6), Q(6), kb(3,3), q(3)
{
    // does nothing
}

ElasticBeam2d::ElasticBeam2d(int tag, double a, double e, double i, 
		     int Nd1, int Nd2, 
		     CrdTransf2d &coordTransf, double r)
:Element(tag,ELE_TAG_ElasticBeam2d), 
 A(a), E(e), I(i), L(0.0), rho(r),
 connectedExternalNodes(2), theCoordTransf(0),
 m(6,6), d(6,6), Pinert(6), Q(6), kb(3,3), q(3)
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    theCoordTransf = coordTransf.getCopy();
    
    if (!theCoordTransf)
	g3ErrorHandler->fatal("ElasticBeam2d::ElasticBeam2d -- failed to get copy of coordinate transformation");
}

ElasticBeam2d::~ElasticBeam2d()
{
    if (theCoordTransf)
	delete theCoordTransf;
}

int
ElasticBeam2d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ElasticBeam2d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

int
ElasticBeam2d::getNumDOF(void)
{
    return 6;
}

void
ElasticBeam2d::setDomain(Domain *theDomain)
{
    if (theDomain == 0)
	g3ErrorHandler->fatal("ElasticBeam2d::setDomain -- Domain is null");
    
    node1Ptr = theDomain->getNode(connectedExternalNodes(0));
    node2Ptr = theDomain->getNode(connectedExternalNodes(1));    
    
    if (node1Ptr == 0)
	g3ErrorHandler->fatal("ElasticBeam2d::setDomain -- Node 1: %i does not exist",
			      connectedExternalNodes(0));
    
    if (node2Ptr == 0)
	g3ErrorHandler->fatal("ElasticBeam2d::setDomain -- Node 2: %i does not exist",
			      connectedExternalNodes(1));
 
    int dofNd1 = node1Ptr->getNumberDOF();
    int dofNd2 = node2Ptr->getNumberDOF();    
    
    if (dofNd1 != 3)
	g3ErrorHandler->fatal("ElasticBeam2d::setDomain -- Node 1: %i has incorrect number of DOF",
			      connectedExternalNodes(0));
    
    if (dofNd2 != 3)
	g3ErrorHandler->fatal("ElasticBeam2d::setDomain -- Node 2: %i has incorrect number of DOF",
			      connectedExternalNodes(1));	
	
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(node1Ptr, node2Ptr) != 0)
	g3ErrorHandler->fatal("ElasticBeam2d::setDomain -- Error initializing coordinate transformation");
    
    L = theCoordTransf->getInitialLength();

    if (L == 0.0)
	g3ErrorHandler->fatal("ElasticBeam2d::setDomain -- Element has zero length");
    
    kb(0,0) = E*A/L;
    kb(1,1) = kb(2,2) = 4*E*I/L;
    kb(1,2) = kb(2,1) = 2*E*I/L;
    
    m(0,0) = m(1,1) = m(3,3) = m(4,4) = rho*L/2;
}

int
ElasticBeam2d::commitState()
{
    return theCoordTransf->commitState();
}

int
ElasticBeam2d::revertToLastCommit()
{
    return theCoordTransf->revertToLastCommit();
}

int
ElasticBeam2d::revertToStart()
{
    return theCoordTransf->revertToStart();
}

const Matrix &
ElasticBeam2d::getTangentStiff(void)
{
    theCoordTransf->update();
    
    const Vector &v = theCoordTransf->getBasicTrialDisp();
    
    q(0) = kb(0,0)*v(0);
    q(1) = kb(1,1)*v(1) + kb(1,2)*v(2);
    q(2) = kb(2,1)*v(1) + kb(2,2)*v(2);
    
    //cerr << theCoordTransf->getGlobalStiffMatrix(kb,q);
	return theCoordTransf->getGlobalStiffMatrix(kb,q);
}

const Matrix &
ElasticBeam2d::getSecantStiff(void)
{
    return this->getTangentStiff();
}

const Matrix &
ElasticBeam2d::getDamp(void)
{
    return d;
}

const Matrix &
ElasticBeam2d::getMass(void)
{ 
    return m;
}

void 
ElasticBeam2d::zeroLoad(void)
{
	Q.Zero();

    return;
}

int
ElasticBeam2d::addLoad(const Vector &moreLoad)
{
    if (moreLoad.Size() != 6) {
	g3ErrorHandler->warning("ElasticBeam2d::addLoad: vector not of correct size");
	return -1;
    }

	//Q += moreLoad;
	Q.addVector(1.0, moreLoad, 1.0);

    return 0;
}

int
ElasticBeam2d::addInertiaLoadToUnbalance(const Vector &accel)
{
	if (rho == 0.0)
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = node1Ptr->getRV(accel);
	const Vector &Raccel2 = node2Ptr->getRV(accel);
	
    if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
		g3ErrorHandler->warning("ElasticBeam2d::addInertiaLoadToUnbalance %s\n",
				"matrix and vector sizes are incompatable");
		return -1;
    }
    
	// Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	// Mass matrix is computed in setDomain()
    Q(0) += -m(0,0) * Raccel1(0);
    Q(1) += -m(1,1) * Raccel1(1);
    
    Q(3) += -m(3,3) * Raccel2(0);    
    Q(4) += -m(4,4) * Raccel2(1);    

    return 0;
}

const Vector &
ElasticBeam2d::getResistingForceIncInertia()
{	
    Pinert = this->getResistingForce();
    
    const Vector &accel1 = node1Ptr->getTrialAccel();
    const Vector &accel2 = node2Ptr->getTrialAccel();    
    
    Pinert(0) += m(0,0) * accel1(0);
    Pinert(1) += m(1,1) * accel1(1);
    
    Pinert(3) += m(3,3) * accel2(0);    
    Pinert(4) += m(4,4) * accel2(1);

	return Pinert;
}

const Vector &
ElasticBeam2d::getResistingForce()
{
    theCoordTransf->update();

    const Vector &v = theCoordTransf->getBasicTrialDisp();
    
    q(0) = kb(0,0)*v(0);
    q(1) = kb(1,1)*v(1) + kb(1,2)*v(2);
    q(2) = kb(2,1)*v(1) + kb(2,2)*v(2);
    
    static Vector dummy(2);
    
	Pinert = theCoordTransf->getGlobalResistingForce(q, dummy);

	//Pinert = Pinert - Q;
	Pinert.addVector(1.0, Q, -1.0);

	//cerr << "getResistingForce -- " << Pinert << endl;

    return Pinert;
}

int
ElasticBeam2d::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;

    static Vector data(9);
    
	data(0) = A;
    data(1) = E; 
    data(2) = I; 
    data(3) = rho;
    data(4) = this->getTag();
    data(5) = connectedExternalNodes(0);
    data(6) = connectedExternalNodes(1);
	data(7) = theCoordTransf->getClassTag();
    	
	int dbTag = theCoordTransf->getDbTag();

	if (dbTag == 0) {
		dbTag = theChannel.getDbTag();
		if (dbTag != 0)
			theCoordTransf->setDbTag(dbTag);
	}

	data(8) = dbTag;

	// Send the data vector
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send data Vector",
			"ElasticBeam2d::sendSelf");
		return res;
    }

    // Ask the CoordTransf to send itself
	res += theCoordTransf->sendSelf(cTag, theChannel);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send CoordTransf",
			"ElasticBeam2d::sendSelf");
		return res;
	}

    return res;
}

int
ElasticBeam2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
	
	static Vector data(9);

    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive data Vector",
			"ElasticBeam2d::recvSelf");
		return res;
    }

    A = data(0);
    E = data(1); 
    I = data(2); 
    rho = data(3);
    this->setTag((int)data(4));
    connectedExternalNodes(0) = (int)data(5);
    connectedExternalNodes(1) = (int)data(6);

	// Check if the CoordTransf is null; if so, get a new one
	int crdTag = (int)data(7);
	if (theCoordTransf == 0) {
		theCoordTransf = theBroker.getNewCrdTransf2d(crdTag);
		if (theCoordTransf == 0) {
			g3ErrorHandler->warning("%s -- could not get a CrdTransf2d",
				"ElasticBeam2d::recvSelf");
			return -1;
		}
	}

	// Check that the CoordTransf is of the right type; if not, delete
	// the current one and get a new one of the right type
	if (theCoordTransf->getClassTag() != crdTag) {
		delete theCoordTransf;
		theCoordTransf = theBroker.getNewCrdTransf2d(crdTag);
		if (theCoordTransf == 0) {
			g3ErrorHandler->warning("%s -- could not get a CrdTransf2d",
				"ElasticBeam2d::recvSelf");
			return -1;
		}
	}

	// Now, receive the CoordTransf
	theCoordTransf->setDbTag((int)data(8));
	res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive CoordTransf",
			"ElasticBeam2d::recvSelf");
		return res;
	}

	// Revert the crdtrasf to its last committed state
	theCoordTransf->revertToLastCommit();

    return res;
}

void
ElasticBeam2d::Print(ostream &s, int flag)
{
    s << "\nElasticBeam2d: " << this->getTag() << endl;
    s << "\tConnected Nodes: " << connectedExternalNodes ;
	s << "\tCoordTransf: " << theCoordTransf->getTag() << endl;
}

int
ElasticBeam2d::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = node1Ptr->getCrds();
    const Vector &end2Crd = node2Ptr->getCrds();	

    const Vector &end1Disp = node1Ptr->getDisp();
    const Vector &end2Disp = node2Ptr->getDisp();

	static Vector v1(3);
	static Vector v2(3);

	for (int i = 0; i < 2; i++) {
		v1(i) = end1Crd(i) + end1Disp(i)*fact;
		v2(i) = end2Crd(i) + end2Disp(i)*fact;    
	}
	
	return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response*
ElasticBeam2d::setResponse(char **argv, int argc, Information &info)
{
    // stiffness
    if (strcmp(argv[0],"stiffness") == 0)
		return new ElementResponse(this, 1, Matrix(6,6));

    // forces
    else if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
		return new ElementResponse(this, 2, Vector(6));

    else
		return 0;
}

int
ElasticBeam2d::getResponse (int responseID, Information &eleInfo)
{
    switch (responseID) {
		case 1: // stiffness
			return eleInfo.setMatrix(this->getTangentStiff());

		case 2: // forces
			return eleInfo.setVector(this->getResistingForce());

		default:
			return -1;
    }
}

