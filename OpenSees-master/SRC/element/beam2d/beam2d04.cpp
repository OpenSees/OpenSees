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
// $Source: /usr/local/cvs/OpenSees/SRC/element/beam2d/beam2d04.cpp,v $
                                                                        
                                                                        
// File: ~/element/beam2d04.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for beam2d04.
// beam2d04 is a 2d plane frame bending member. As such it can only
// connect to a node with 3-dof. It only uses the x and y coordinates
// at the nodes, a z coordinate is not used. 
//
//                                      5
//        2                             |<
//        |                   .=========+-)-4
//       3|      .=============         | 6
// 1 ---(-+======
//       >|
//
// The interface:
//


#include "beam2d04.h"
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <math.h>
#include <stdlib.h>

Matrix beam2d04::k(6,6);
Matrix beam2d04::trans(6,6);

// beam2d04(int tag, double A, double E, double I, int Nd1, int Nd2);
//	constructor which takes the unique element tag, the elements A,E and
//	I and the node ID's of it's nodal end points. 

beam2d04::beam2d04()
   :Element(0,ELE_TAG_beam2d04), A(0), E(0), I(0), L(0),
    connectedExternalNodes(2), 
    rForce(6), load(6), isStiffFormed(0)
{
  theNodes[0] = 0;
  theNodes[1] = 0;
}

beam2d04::beam2d04(int tag, double a, double e, double i, int Nd1, int Nd2)
   :Element(tag,ELE_TAG_beam2d04), A(a), E(e), I(i), L(0),
    connectedExternalNodes(2), 
    rForce(6), load(6), isStiffFormed(0)
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;    
    
    theNodes[0] = 0;
    theNodes[1] = 0;
}


// ~beam2d04():
// 	destructor

beam2d04::~beam2d04()
{
}



int
beam2d04::getNumExternalNodes(void) const
{
    return connectedExternalNodes.Size();
}

const ID &
beam2d04::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
beam2d04::getNodePtrs(void) 
{
  return theNodes;
}

int
beam2d04::getNumDOF(void) {
    int i =6;
    return i;
}


void
beam2d04::formVar(void)
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
	theNodes[0] = end1Ptr;
	theNodes[1] = end2Ptr;

	if (end1Ptr == 0) {
	    opserr << "beam2d04::formVar: Nd1: ";
	    opserr << Nd1 << "does not exist in model\n";
	    exit(0);
	}
	if (end2Ptr == 0) {
	    opserr << "beam2d04::formVar: Nd2: ";
	    opserr << Nd2 << "does not exist in model\n";
	    exit(0);
	}

	double dx,dy;
	const Vector &end1Crd = end1Ptr->getCrds();
	const Vector &end2Crd = end2Ptr->getCrds();	
    
	dx = end2Crd(0)-end1Crd(0);
	dy = end2Crd(1)-end1Crd(1);	
    
	L = sqrt(dx*dx + dy*dy);
	double L2 = L*L;
	double L3 = L*L*L;
	if (L == 0.0) {
	    opserr << "Element: " << this->getTag();
	    opserr << " beam2d04::formVar: 0 length\n";
	    exit(-1);
	}
	
	cs = dx/L;
	sn = dy/L;

	oneEA = E*A/L;
	twoEI = 2*E*I/L;
	fourEI = 4*E*I/L;
	twelveEI = 12*E*I/L3;
	sixEI = 6*E*I/L2;
    }
    isStiffFormed = 1;
}    


int
beam2d04::revertToLastCommit()
{
    return 0;
    // linear element - nothing to commit
}

const Matrix &
beam2d04::getTangentStiff(void)
{
    return this->getStiff();
}

const Matrix &
beam2d04::getInitialStiff(void)
{
    return this->getStiff();
}




// const Matrix &getStiff():
//	Method to return the stiffness matrix.

const Matrix &
beam2d04::getStiff(void)
{
    if (isStiffFormed == 0)
        this->formVar();

    if (cs == 1.0) {
	k(0,0) = oneEA;
	k(1,0) = 0;
	k(2,0) = 0;
	k(3,0) = -oneEA;
	k(4,0) = 0;
	k(5,0) = 0;

	k(0,1) = 0;
	k(1,1) = twelveEI;
	k(2,1)= sixEI;
	k(3,1) = 0;
	k(4,1) = -twelveEI;
	k(5,1) = sixEI;

	k(0,2) = 0;
	k(1,2) = sixEI;
	k(2,2) = fourEI;
	k(3,2) = 0;
	k(4,2) = -sixEI;
	k(5,2) = twoEI;
	
	k(0,3) = -oneEA;
	k(1,3) = 0;
	k(2,3) = 0;
	k(3,3) = oneEA;
	k(4,3) = 0;
	k(5,3) = 0;

	k(0,4) = 0;
	k(1,4) = -twelveEI;
	k(2,4) = -sixEI;
	k(3,4) = 0;
	k(4,4) = twelveEI;
	k(5,4)  = -sixEI;
	
	k(0,5) = 0;
	k(1,5) = sixEI;
	k(2,5) = twoEI;
	k(3,5) = 0;
	k(4,5) = -sixEI;
	k(5,5) = fourEI;
    }
    else if (sn == 1.0) {
	k(0,0) = twelveEI;
	k(1,0) = 0;
	k(2,0) = -sixEI;
	k(3,0) = -twelveEI;;
	k(4,0) = 0;
	k(5,0) = -sixEI;;

	k(0,1) = 0;
	k(1,1) = oneEA;
	k(2,1) = 0;
	k(3,1) = 0;
	k(4,1) = -oneEA;
	k(5,1) = 0;
	
	k(0,2) = -sixEI;
	k(1,2) = 0;
	k(2,2) = fourEI;
	k(3,2) = sixEI;
	k(4,2) = 0;
	k(5,2) = twoEI;
	
	k(0,3) = -twelveEI;
	k(1,3) = 0;
	k(2,3) = sixEI;
	k(3,3) = twelveEI;
	k(4,3) = 0;
	k(5,3) = sixEI;

	k(0,4) = 0;
	k(1,4) = -oneEA;
	k(2,4) = 0;
	k(3,4) = 0;
	k(4,4) = oneEA;
	k(5,4) = 0;
	
	k(0,5) = -sixEI;
	k(1,5) = 0;
	k(2,5) = twoEI;
	k(3,5) = sixEI;
	k(4,5) = 0;
	k(5,5) = fourEI;
    }
    else {
	opserr << "beam2d04::getStiff - more WORK \n";
	exit(0);
    }

    return k;
}
    


void 
beam2d04::zeroLoad(void)
{
    load.Zero();
}

int 
beam2d04::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "beam2d04::addLoad() - beam " << this->getTag() << ",load type unknown\n"; 
  return -1;
}

int
beam2d04::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
beam2d04::getResistingForceIncInertia()
{	
    this->getResistingForce();

    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      rForce += this->getRayleighDampingForces();
    
    return rForce;
}




const Vector &
beam2d04::getResistingForce()
{	
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
    rForce(3) = end2Disp(0);
    rForce(4) = end2Disp(1);
    rForce(5) = end2Disp(2);    

    rForce = k * rForce;
    
    // add any applied load
    rForce -= load;
    
    return rForce;
}

int
beam2d04::sendSelf(int commitTag, Channel &theChannel)
{
  int dataTag = this->getDbTag();

  Vector data(4);
  data(0) = A; data(1) = E; data(2) = I; data(3) = this->getTag();
    
  int result = 0;
  result = theChannel.sendVector(dataTag, commitTag, data);
  if (result < 0) {
    opserr << "beam2d04::sendSelf - failed to send data\n";
    return -1;
  }
  
  result = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (result < 0) {
    opserr << "beam2d04::sendSelf - failed to send data\n";
    return -1;
  }
    
  return 0;
}

int
beam2d04::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    int result = 0;
    int dataTag = this->getDbTag();

    result = theChannel.recvVector(dataTag, commitTag, data);
    if (result < 0) {
	opserr << "beam2d04::recvSelf - failed to recv data\n";
	return -1;
    }

    A = data(0); E = data(1); I=data(2); 
    this->setTag((int)data(3));

    result = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
    if (result < 0) {
	opserr << "beam2d04::recvSelf - failed to recv data\n";
	return -1;
    }
    
    return 0;
}

void
beam2d04::Print(OPS_Stream &s, int flag)
{
    s << "\nElement: " << this->getTag() << " Type: beam2d04 ";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
//    s << "\tStiffness Matrix:\n" << k;
    s << "\tResisting Force: " << rForce;
//    s << "\tElemt End Force: " << load;        
//    s << "\tElemt End Force: " << eForce;    
}





