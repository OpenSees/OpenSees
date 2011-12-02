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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-14 23:01:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beam3d/beam3d01.cpp,v $
                                                                        
                                                                        
// File: ~/model/beam3d01.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for beam3d01.
// beam3d01 is a 3d beam element. As such it can only
// connect to a node with 6-dof. 

#include <beam3d01.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <math.h>
#include <stdlib.h>
#include <Renderer.h>

Matrix beam3d01::m(12,12);  // these beam members have no mass or damping matrice
Matrix beam3d01::d(12,12);


beam3d01::beam3d01()
:Element(0,ELE_TAG_beam3d01), 
 A(0), E(0), G(0), Jx(0), Iy(0), Iz(0),theta(0),
 k(12,12), rForce(12), load(12), 
 connectedExternalNodes(2), isStiffFormed(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;
}

beam3d01::beam3d01(int tag, double a, double e, double g, 
		   double jx, double iy, double iz, int Nd1, int Nd2, 
		   double Theta)

:Element(tag,ELE_TAG_beam3d01), 
 A(a), E(e), G(g), Jx(jx), Iy(iy), Iz(iz), theta(Theta), L(0),
 k(12,12), rForce(12), load(12), 
 connectedExternalNodes(2), isStiffFormed(0)
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;    

  theNodes[0] = 0;
  theNodes[1] = 0;
}


// ~beam3d01():
// 	destructor

beam3d01::~beam3d01()
{

}



int
beam3d01::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
beam3d01::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
beam3d01::getNodePtrs(void) 
{
  return theNodes;
}

int
beam3d01::getNumDOF(void) {
    int i =12;
    return i;
}


int
beam3d01::revertToLastCommit()
{
    return 0;
    // linear element - nothing to commit
}

const Matrix &
beam3d01::getTangentStiff(void)
{
    return this->getStiff();
}

const Matrix &
beam3d01::getInitialStiff(void)
{
    return this->getStiff();
}

// const Matrix &getStiff():
//	Method to return the stiffness matrix.

const Matrix &
beam3d01::getStiff(void)
{
    // compute the stiffness
    if (isStiffFormed == 0) { 
    
	// first check element has correct number of DOF at attached nodes
	int Nd1, Nd2;
	Nd1 = connectedExternalNodes(0);
	Nd2 = connectedExternalNodes(1);

	Domain *theDomain = this->getDomain();
	Node *end1Ptr = theDomain->getNode(Nd1);
	Node *end2Ptr = theDomain->getNode(Nd2);	

	if (end1Ptr == 0) {
	    opserr << "beam3d01::getStiff: Nd1: ";
	    opserr << Nd1 << "does not exist in model\n";
	    exit(0);
	}
	if (end2Ptr == 0) {
	    opserr << "beam3d01::getStiff: Nd2: ";
	    opserr << Nd2 << "does not exist in model\n";
	    exit(0);
	}
	
	theNodes[0] = end1Ptr;
	theNodes[1] = end2Ptr;

	double dx,dy,dz;
	const Vector &end1Crd = end1Ptr->getCrds();
	const Vector &end2Crd = end2Ptr->getCrds();	
    
	dx = end2Crd(0)-end1Crd(0);
	dy = end2Crd(1)-end1Crd(1);	
	dz = end2Crd(2)-end1Crd(2);		
    
	L = sqrt(dx*dx + dy*dy + dz*dz);
	double L2 = L*L;
	double L3 = L*L*L;
	if (L == 0.0) {
	    opserr << "Element: " << this->getTag();
	    opserr << " beam3d01::getStiff: 0 length\n";
	    return k;  
	}
	
	double EA = E*A/L;
	double twoE = 2*E/L;
	double fourE = 4*E/L;
	double twelveE = 12*E/L3;
	double sixE = 6*E/L2;
	
	if (dy == 0.0 && dz == 0.0 && dx > 0.0 && theta == 90) {
	    // local y in y and local z in z
	    k(0,0) = EA;  
	    k(6,0) = -EA;

	    k(1,1) = twelveE*Iz;
	    k(5,1) = sixE*Iz;
	    k(7,1) = -twelveE*Iz;
	    k(11,1) = sixE*Iz;
	    
	    k(2,2) = twelveE*Iy;
	    k(4,2) = -sixE*Iy;
	    k(8,2) = -twelveE*Iy;
	    k(10,2) = -sixE*Iy;	    
	    
	    k(3,3) = G*Jx/L;
	    k(9,3) = -G*Jx/L;
	    
	    k(2,4) = -sixE*Iy;
	    k(4,4) = fourE*Iy;
	    k(8,4) = sixE*Iy;
	    k(10,4) = twoE*Iy;
	    
	    k(1,5) = sixE*Iz;
	    k(5,5) = fourE*Iz;
	    k(7,5) = -sixE*Iz;
	    k(11,5) = twoE*Iz;	    
	    
	    k(0,6) = -EA;  
	    k(6,6) = EA;

	    k(1,7) = -twelveE*Iz;
	    k(5,7) = -sixE*Iz;
	    k(7,7) = twelveE*Iz;
	    k(11,7) = -sixE*Iz;
	    
	    k(2,8) = -twelveE*Iy;
	    k(4,8) = sixE*Iy;
	    k(8,8) = twelveE*Iy;
	    k(10,8) = sixE*Iy;	    
	    
	    k(3,9) = -G*Jx/L;
	    k(9,9) = G*Jx/L;
	    
	    k(2,10) = -sixE*Iy;
	    k(4,10) = twoE*Iy;
	    k(8,10) = sixE*Iy;
	    k(10,10) = fourE*Iy;
	    
	    k(1,11) = sixE*Iz;
	    k(5,11) = twoE*Iz;
	    k(7,11) = -sixE*Iz;
	    k(11,11) = fourE*Iz;	    	    
	}
	
	else if (dx == 0.0 && dz == 0.0 && dy > 0.0 && theta == 90) {

	    k(0,0) = twelveE*Iz;
	    k(5,0) = -sixE*Iz;
	    k(6,0) = -twelveE*Iz;
	    k(11,0) = -sixE*Iz;
	    
	    k(1,1) = EA;  
	    k(7,1) = -EA;

	    k(2,2) = twelveE*Iy;
	    k(3,2) = sixE*Iy;
	    k(8,2) = -twelveE*Iy;
	    k(9,2) = sixE*Iy;	    

	    k(2,3) = sixE*Iy;
	    k(3,3) = fourE*Iy;
	    k(8,3) = -sixE*Iy;
	    k(9,3) = twoE*Iy;
	    
	    k(4,4) = G*Jx/L;
	    k(10,4) = -G*Jx/L;
	    
	    k(0,5) = -sixE*Iz;
	    k(5,5) = fourE*Iz;
	    k(6,5) = sixE*Iz;
	    k(11,5) = twoE*Iz;	    
	    
	    k(0,6) = -twelveE*Iz;
	    k(5,6) = sixE*Iz;
	    k(6,6) = twelveE*Iz;
	    k(11,6) = sixE*Iz;
	    
	    k(1,7) = -EA;  
	    k(7,7) = EA;

	    k(2,8) = -twelveE*Iy;
	    k(3,8) = -sixE*Iy;
	    k(8,8) = twelveE*Iy;
	    k(9,8) = -sixE*Iy;	    

	    k(2,9) = sixE*Iy;
	    k(3,9) = twoE*Iy;
	    k(8,9) = -sixE*Iy;
	    k(9,9) = fourE*Iy;
	    
	    k(4,10) = -G*Jx/L;
	    k(10,10) = G*Jx/L;
	    
	    k(0,11) = -sixE*Iz;
	    k(5,11) = twoE*Iz;
	    k(6,11) = sixE*Iz;
	    k(11,11) = fourE*Iz;	    	    
	}	    
	
	else if (dx == 0.0 && dy == 0.0 && dz > 0.0 && theta == 90) {
	    // local y of columns in x dirn, local z in y dirn
	    k(0,0) = twelveE*Iz;
	    k(4,0) = sixE*Iz;
	    k(6,0) = -twelveE*Iz;
	    k(10,0) = sixE*Iz;	    

	    k(1,1) = twelveE*Iy;
	    k(3,1) = -sixE*Iy;
	    k(7,1) = -twelveE*Iy;
	    k(9,1) = -sixE*Iy;

	    k(2,2) = EA;  
	    k(8,2) = -EA;

	    k(1,3) = -sixE*Iy;
	    k(3,3) = fourE*Iy;
	    k(7,3) = sixE*Iy;
	    k(9,3) = twoE*Iy;	    	    

	    k(0,4) = sixE*Iz;
	    k(4,4) = fourE*Iz;
	    k(6,4) = -sixE*Iz;
	    k(10,4) = twoE*Iz;

	    k(5,5) = G*Jx/L;
	    k(11,5) = -G*Jx/L;	    
	    
	    k(0,6) = -twelveE*Iz;
	    k(4,6) = -sixE*Iz;
	    k(6,6) = twelveE*Iz;
	    k(10,6) = -sixE*Iz;	    

	    k(1,7) = -twelveE*Iy;
	    k(3,7) = sixE*Iy;
	    k(7,7) = twelveE*Iy;
	    k(9,7) = sixE*Iy;

	    k(2,8) = -EA;  
	    k(8,8) = EA;

	    k(1,9) = -sixE*Iy;
	    k(3,9) = twoE*Iy;
	    k(7,9) = sixE*Iy;
	    k(9,9) = fourE*Iy;	    	    

	    k(0,10) = sixE*Iz;
	    k(4,10) = twoE*Iz;
	    k(6,10) = -sixE*Iz;
	    k(10,10) = fourE*Iz;

	    k(5,11) = -G*Jx/L;
	    k(11,11) = G*Jx/L;	    
	}	    
	else {
	    opserr << "beam3d01::getStiff - NOT FINISHED";
	    opserr << " members not located along global axis directions\n";
	    exit(0);
	    
	}
	
	isStiffFormed = 1;
    }

    return k;
}
    
void 
beam3d01::zeroLoad(void)
{
    load.Zero();
}

int 
beam3d01::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "beam3d01::addLoad() - beam " << this->getTag() << "load type unknown\n";
  return -1;
}

int
beam3d01::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}


const Vector &
beam3d01::getResistingForceIncInertia()
{	
    this->getResistingForce();

    // add rayleigh damping force if factors present    
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	rForce += this->getRayleighDampingForces();

    return rForce;
}

const Vector &
beam3d01::getResistingForce()
{	
    if (isStiffFormed == 0)
	this->getStiff();
    
    // compute the residual Res = k*uTrial
    int Nd1, Nd2;
    Nd1 = connectedExternalNodes(0);
    Nd2 = connectedExternalNodes(1);
    Domain *theDomain = this->getDomain();
    Node *end1Ptr = theDomain->getNode(Nd1);
    Node *end2Ptr = theDomain->getNode(Nd2);	
    
    const Vector &end1Disp = end1Ptr->getTrialDisp();
    const Vector &end2Disp = end2Ptr->getTrialDisp();    
    rForce(0) = end1Disp(0);
    rForce(1) = end1Disp(1);
    rForce(2) = end1Disp(2);    
    rForce(3) = end1Disp(3);
    rForce(4) = end1Disp(4);
    rForce(5) = end1Disp(5);    
    rForce(6) = end2Disp(0);
    rForce(7) = end2Disp(1);
    rForce(8) = end2Disp(2);    
    rForce(9) = end2Disp(3);
    rForce(10) = end2Disp(4);
    rForce(11) = end2Disp(5);        

    rForce = k * rForce;
    
    // add any applied load
    rForce -= load;
    
    return rForce;
}


int
beam3d01::displaySelf(Renderer &theViewer, int displayMode, float fact)
{

    Domain *theDomain = this->getDomain();
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    Node *end1Ptr = theDomain->getNode(Nd1);
    Node *end2Ptr = theDomain->getNode(Nd2);	
    
    const Vector &end1Crd = end1Ptr->getCrds();
    const Vector &end2Crd = end2Ptr->getCrds();	
    const Vector &end1Disp = end1Ptr->getDisp();
    const Vector &end2Disp = end2Ptr->getDisp();    

    Vector v1(3);
    Vector v2(3);
    for (int i=0; i<3; i++) {
	v1(i) = end1Crd(i)+end1Disp(i)*fact;
	v2(i) = end2Crd(i)+end2Disp(i)*fact;    
    }
    return theViewer.drawLine(v1,v2,1.0,1.0);
}


int
beam3d01::sendSelf(int commitTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    Vector data(10);
    data(0) = A; data(1) = E; data(2) = G; 
    data(3) = Jx; data(4) = Iy; data(5) = Iz;     
    data(6) = this->getTag();
    data(7) = connectedExternalNodes(0);
    data(8) = connectedExternalNodes(1);
    data(9) = theta;    
    int result = 0;
    result = theChannel.sendVector(dataTag, commitTag, data);
    if (result < 0) {
	opserr << "beam3d01::sendSelf - failed to send data\n";
	return -1;
    }
    
    return 0;
}

int
beam3d01::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(10);
    int dataTag = this->getDbTag();
    int result = 0;

    result = theChannel.recvVector(dataTag, commitTag, data);
    if (result < 0) {
	opserr << "beam3d01::recvSelf - failed to recv data\n";
	return -1;
    }

    A = data(0); E = data(1); G=data(2); 
    Jx = data(3); Iy = data(4); Iz=data(5);     
    theta = data(9);
    int tag = data(6);
    this->setTag(tag);
    int nd1 = data(7);
    int nd2 = data(8);
    connectedExternalNodes(0) = nd1;
    connectedExternalNodes(1) = nd2;    

    return 0;
}


void
beam3d01::Print(OPS_Stream &s, int flag)
{
    s << "Element: " << this->getTag(); 
    s << " type: beam3d01  iNode: " << connectedExternalNodes(0);
    s << " jNode: " << connectedExternalNodes(1) << " Length: " << L << endln;
//    s << "\tStiffness Matrix:\n" << k;
    s << "\tResisting Force: " << rForce;
//    s << "\tElemt End Force: " << eForce;    
}






