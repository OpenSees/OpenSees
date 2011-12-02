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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beam3d/ElasticBeam3d.cpp,v $
                                                                        
                                                                        
// File: ~/model/ElasticBeam3d.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for ElasticBeam3d.
// ElasticBeam3d is a 3d beam element. As such it can only
// connect to a node with 6-dof. 

#include <ElasticBeam3d.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf3d.h>
#include <Information.h>

#include <math.h>
#include <stdlib.h>

ElasticBeam3d::ElasticBeam3d()
:Element(0,ELE_TAG_ElasticBeam3d), 
 A(0), E(0), G(0), Jx(0), Iy(0), Iz(0), L(0.0), rho(0.0),
 connectedExternalNodes(2), theCoordTransf(0),
 m(12,12), d(12,12), Pinert(12), Q(12), kb(6,6), q(6)
{
    // does nothing
}

ElasticBeam3d::ElasticBeam3d(int tag, double a, double e, double g, 
		   double jx, double iy, double iz, int Nd1, int Nd2, 
		   CrdTransf3d &coordTransf, double r)
:Element(tag,ELE_TAG_ElasticBeam3d), 
 A(a), E(e), G(g), Jx(jx), Iy(iy), Iz(iz), L(0.0), rho(r),
 connectedExternalNodes(2), theCoordTransf(0),
 m(12,12), d(12,12), Pinert(12), Q(12), kb(6,6), q(6)
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;
    
    theCoordTransf = coordTransf.getCopy();
    
    if (!theCoordTransf)
	g3ErrorHandler->fatal("ElasticBeam3d::ElasticBeam3d -- failed to get copy of coordinate transformation");
}

ElasticBeam3d::~ElasticBeam3d()
{
    if (theCoordTransf)
	delete theCoordTransf;
}

int
ElasticBeam3d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ElasticBeam3d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

int
ElasticBeam3d::getNumDOF(void)
{
    return 12;
}

void
ElasticBeam3d::setDomain(Domain *theDomain)
{
    if (theDomain == 0)
	g3ErrorHandler->fatal("ElasticBeam3d::setDomain -- Domain is null");
    
    node1Ptr = theDomain->getNode(connectedExternalNodes(0));
    node2Ptr = theDomain->getNode(connectedExternalNodes(1));    
    
    if (node1Ptr == 0)
	g3ErrorHandler->fatal("ElasticBeam3d::setDomain -- Node 1: %i does not exist",
			      connectedExternalNodes(0));
    
    if (node2Ptr == 0)
	g3ErrorHandler->fatal("ElasticBeam3d::setDomain -- Node 2: %i does not exist",
			      connectedExternalNodes(1));
 
    int dofNd1 = node1Ptr->getNumberDOF();
    int dofNd2 = node2Ptr->getNumberDOF();    
    
    if (dofNd1 != 6)
	g3ErrorHandler->fatal("ElasticBeam3d::setDomain -- Node 1: %i has incorrect number of DOF",
			      connectedExternalNodes(0));
    
    if (dofNd2 != 6)
	g3ErrorHandler->fatal("ElasticBeam3d::setDomain -- Node 2: %i has incorrect number of DOF",
			      connectedExternalNodes(1));	
	
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(node1Ptr, node2Ptr) != 0)
	g3ErrorHandler->fatal("ElasticBeam3d::setDomain -- Error initializing coordinate transformation");
    
    L = theCoordTransf->getInitialLength();

    if (L == 0.0)
	g3ErrorHandler->fatal("ElasticBeam3d::setDomain -- Element has zero length");
    
    kb(0,0) = E*A/L;
    kb(1,1) = kb(2,2) = 4*E*Iz/L;
    kb(1,2) = kb(2,1) = 2*E*Iz/L;
    kb(3,3) = kb(4,4) = 4*E*Iy/L;
    kb(3,4) = kb(4,3) = 2*E*Iy/L;
    kb(5,5) = G*Jx/L;
    
    m(0,0) = m(1,1) = m(2,2) = m(6,6) = m(7,7) = m(8,8) = rho*L/2;
}

int
ElasticBeam3d::commitState()
{
    return theCoordTransf->commitState();
}

int
ElasticBeam3d::revertToLastCommit()
{
    return theCoordTransf->revertToLastCommit();
}

int
ElasticBeam3d::revertToStart()
{
    return theCoordTransf->revertToStart();
}

const Matrix &
ElasticBeam3d::getTangentStiff(void)
{
    theCoordTransf->update();
    
    const Vector &v = theCoordTransf->getBasicTrialDisp();
    
    q(0) = kb(0,0)*v(0);
    q(1) = kb(1,1)*v(1) + kb(1,2)*v(2);
    q(2) = kb(2,1)*v(1) + kb(2,2)*v(2);
    q(3) = kb(3,3)*v(3) + kb(3,4)*v(4);
    q(4) = kb(4,3)*v(3) + kb(4,4)*v(4);    
    q(5) = kb(5,5)*v(5);    
    
    return theCoordTransf->getGlobalStiffMatrix(kb,q);
}

const Matrix &
ElasticBeam3d::getSecantStiff(void)
{
    return this->getTangentStiff();
}

const Matrix &
ElasticBeam3d::getDamp(void)
{
    return d;
}

const Matrix &
ElasticBeam3d::getMass(void)
{ 
    return m;
}

void 
ElasticBeam3d::zeroLoad(void)
{
	Q.Zero();

    return;
}

int
ElasticBeam3d::addLoad(const Vector &moreLoad)
{
    if (moreLoad.Size() != 12) {
	g3ErrorHandler->warning("ElasticBeam3d::addLoad: vector not of correct size");
	return -1;
    }

	//Q += moreLoad;
	Q.addVector(1.0, moreLoad, 1.0);

    return 0;
}

int
ElasticBeam3d::addInertiaLoadToUnbalance(const Vector &accel)
{
	if (rho == 0.0)
		return 0;

	// Get R * accel from the nodes
	const Vector &Raccel1 = node1Ptr->getRV(accel);
	const Vector &Raccel2 = node2Ptr->getRV(accel);
	
    if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
		g3ErrorHandler->warning("ElasticBeam3d::addInertiaLoadToUnbalance %s\n",
				"matrix and vector sizes are incompatable");
		return -1;
    }

	// Want to add ( - fact * M R * accel ) to unbalance
	// Take advantage of lumped mass matrix
	// Mass matrix is computed in setDomain()
    Q(0) += -m(0,0) * Raccel1(0);
    Q(1) += -m(1,1) * Raccel1(1);
    Q(2) += -m(2,2) * Raccel1(2);
    
    Q(6) += -m(6,6) * Raccel2(0);    
    Q(7) += -m(7,7) * Raccel2(1);
    Q(8) += -m(8,8) * Raccel2(2);    

    return 0;
}

const Vector &
ElasticBeam3d::getResistingForceIncInertia()
{	
    Pinert = this->getResistingForce();
    
    const Vector &accel1 = node1Ptr->getTrialAccel();
    const Vector &accel2 = node2Ptr->getTrialAccel();    
    
    Pinert(0) += m(0,0) * accel1(0);
    Pinert(1) += m(1,1) * accel1(1);
    Pinert(2) += m(2,2) * accel1(2);
    
    Pinert(6) += m(6,6) * accel2(0);    
    Pinert(7) += m(7,7) * accel2(1);
    Pinert(8) += m(8,8) * accel2(2);    

	return Pinert;
}

const Vector &
ElasticBeam3d::getResistingForce()
{
    theCoordTransf->update();

    const Vector &v = theCoordTransf->getBasicTrialDisp();
    
    q(0) = kb(0,0)*v(0);
    q(1) = kb(1,1)*v(1) + kb(1,2)*v(2);
    q(2) = kb(2,1)*v(1) + kb(2,2)*v(2);
    q(3) = kb(3,3)*v(3) + kb(3,4)*v(4);
    q(4) = kb(4,3)*v(3) + kb(4,4)*v(4);    
    q(5) = kb(5,5)*v(5);
    
    static Vector dummy(3);
    
	Pinert = theCoordTransf->getGlobalResistingForce(q, dummy);

	//Pinert = Pinert - Q;
	Pinert.addVector(1.0, Q, -1.0);

    return Pinert;

}

int
ElasticBeam3d::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;

    static Vector data(12);
    
    data(0) = A;
    data(1) = E; 
    data(2) = G; 
    data(3) = Jx; 
    data(4) = Iy; 
    data(5) = Iz;     
    data(6) = rho;
    data(7) = this->getTag();
    data(8) = connectedExternalNodes(0);
    data(9) = connectedExternalNodes(1);
	data(10) = theCoordTransf->getClassTag();    	
	
	int dbTag = theCoordTransf->getDbTag();

	if (dbTag == 0) {
		dbTag = theChannel.getDbTag();
		if (dbTag != 0)
			theCoordTransf->setDbTag(dbTag);
	}

	data(11) = dbTag;

	// Send the data vector
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send data Vector",
			"ElasticBeam3d::sendSelf");
		return res;
    }

    // Ask the CoordTransf to send itself
	res += theCoordTransf->sendSelf(cTag, theChannel);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send CoordTransf",
			"ElasticBeam3d::sendSelf");
		return res;
	}

    return res;
}

int
ElasticBeam3d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;
	
	static Vector data(12);

    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive data Vector",
			"ElasticBeam3d::recvSelf");
		return res;
    }

    A = data(0);
    E = data(1); 
    G = data(2); 
    Jx = data(3); 
    Iy = data(4); 
    Iz = data(5);     
    rho = data(6);
    this->setTag((int)data(7));
    connectedExternalNodes(0) = (int)data(8);
    connectedExternalNodes(1) = (int)data(9);

	// Check if the CoordTransf is null; if so, get a new one
	int crdTag = (int)data(10);
	if (theCoordTransf == 0) {
		theCoordTransf = theBroker.getNewCrdTransf3d(crdTag);
		if (theCoordTransf == 0) {
			g3ErrorHandler->warning("%s -- could not get a CrdTransf3d",
				"ElasticBeam3d::recvSelf");
			return -1;
		}
	}

	// Check that the CoordTransf is of the right type; if not, delete
	// the current one and get a new one of the right type
	if (theCoordTransf->getClassTag() != crdTag) {
		delete theCoordTransf;
		theCoordTransf = theBroker.getNewCrdTransf3d(crdTag);
		if (theCoordTransf == 0) {
			g3ErrorHandler->warning("%s -- could not get a CrdTransf2d",
				"ElasticBeam3d::recvSelf");
			return -1;
		}
	}

	// Now, receive the CoordTransf
	theCoordTransf->setDbTag((int)data(11));
	res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive CoordTransf",
			"ElasticBeam3d::recvSelf");
		return res;
	}

	// Revert the crdtrasf to its last committed state
	theCoordTransf->revertToLastCommit();

    return res;
}

void
ElasticBeam3d::Print(ostream &s, int flag)
{
    s << "\nElasticBeam3d: " << this->getTag() << endl;
    s << "\tConnected Nodes: " << connectedExternalNodes ;
	s << "\tCoordTransf: " << theCoordTransf->getTag() << endl;
}

int
ElasticBeam3d::setResponse (char **argv, int argc, Information &info)
{
    // Basic stiffness
    if (strcmp(argv[0],"stiffness") == 0) {
        Matrix *newMatrix = new Matrix(6,6);
        if (newMatrix == 0) {
            g3ErrorHandler->warning("WARNING ElasticBeam3d::setResponse() - %d out of memory creating matrix\n",
                                    this->getTag());
            return -1;
        }
        info.theMatrix = newMatrix;
        info.theType = MatrixType;
        return 1;
    }

    // Basic forces
    else if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
	Vector *newVector = new Vector(6);
	if (newVector == 0) {
	    g3ErrorHandler->warning("WARNING ElasticBeam3d::setResponse() - %d out of memory creating vector\n",
                                    this->getTag());
            return -1;
	}
	info.theVector = newVector;
	info.theType = VectorType;
	return 2;
    }
    else
	return -1;
}

int
ElasticBeam3d::getResponse (int responseID, Information &info)
{
    switch (responseID) {
      case -1:
	return -1;
      case 1: // Basic stiffness
	if (info.theMatrix != 0)
	    *(info.theMatrix) = kb;
	return 0;
      case 2: // Basic forces
	if (info.theVector != 0)
	    *(info.theVector) = q;
	return 0;
      default:
	return -1;
    }
}

