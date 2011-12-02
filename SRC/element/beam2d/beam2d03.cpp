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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-03-07 17:14:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beam2d/beam2d03.cpp,v $
                                                                        
                                                                        
// File: ~/element/beam2d03.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the class definition for beam2d03.
// beam2d03 is a 2d plane frame bending member. As such it can only
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


#include "beam2d03.h"
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <math.h>
#include <stdlib.h>

// beam2d03(int tag, double A, double E, double I, int Nd1, int Nd2);
//	constructor which takes the unique element tag, the elements A,E and
//	I and the node ID's of it's nodal end points. 

beam2d03::beam2d03()
   :Element(0,ELE_TAG_beam2d03), A(0), E(0), I(0), L(0),
    connectedExternalNodes(2),
    k(6,6), rForce(6), load(6), 
    trans(6,6)
{
  theNodes[0] = 0;
  theNodes[1] = 0;
}

beam2d03::beam2d03(int tag, double a, double e, double i, int Nd1, int Nd2)
   :Element(tag,ELE_TAG_beam2d03), A(a),  E(e), I(i), 
    L(0), connectedExternalNodes(2),
    k(6,6), rForce(6), 
    load(6), trans(6,6)
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;    

    theNodes[0] = 0;
    theNodes[1] = 0;
}


// ~beam2d03():
// 	destructor

beam2d03::~beam2d03()
{
}



int
beam2d03::getNumExternalNodes(void) const
{
    return connectedExternalNodes.Size();
}

const ID &
beam2d03::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
beam2d03::getNodePtrs(void) 
{
  return theNodes;
}

int
beam2d03::getNumDOF(void) {
    int i =6;
    return i;
}


int
beam2d03::revertToLastCommit()
{
    return 0;
    // linear element - nothing to commit
}



void
beam2d03::setDomain(Domain *theDomain)
{
    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	
    if (theNodes[0] == 0) {
	opserr << "WARNING beam2d02::setDomain(): Nd1: ";
	opserr << Nd1 << "does not exist in model for beam \n" << *this;
	return;
    }
    if (theNodes[1] == 0) {
	opserr << "WARNING beam2d02::setDomain(): Nd2: ";
	opserr << Nd2 << "does not exist in model for beam\n" << *this;
	return;
    }	
    
    // now verify the number of dof at node ends
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	
    if (dofNd1 != 3 && dofNd2 != 3) {
	opserr << "WARNING beam2d02::setDomain(): node " << Nd1;
	opserr << " and/or node " << Nd2 << " have/has incorrect number ";
	opserr << "of dof's at end for beam\n " << *this;
	return;
    }	

    // call the base class method
    this->DomainComponent::setDomain(theDomain);

    double dx,dy;
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    
    dx = end2Crd(0)-end1Crd(0);
    dy = end2Crd(1)-end1Crd(1);	
    
    L = sqrt(dx*dx + dy*dy);
    double L2 = L*L;
    double L3 = L*L*L;
    if (L == 0.0) {
      opserr << "Element: " << this->getTag();
      opserr << " beam2d03::getStiff: 0 length\n";
    }
    
    cs = dx/L;
    sn = dy/L;
    
    trans(0,0) = trans(3,3) = cs;
    trans(0,1) = trans(3,4) = sn;
    trans(1,0) = trans(4,3) = -sn;
    trans(1,1) = trans(4,4) = cs;
    trans(2,2) = trans(5,5) = 1.0;
    
    double oneEA = E*A/L;
    double twoEI = 2*E*I/L;
    double fourEI = 4*E*I/L;
    double twelveEI = 12*E*I/L3;
    double sixEI = 6*E*I/L2;
    
    if (sn == 1.0) {
      k(0,0) = twelveEI;
      k(2,0) = -sixEI;
      k(3,0) = -twelveEI;;
      k(5,0) = -sixEI;;
      
      k(1,1) = oneEA;
      k(4,1) = -oneEA;
      
      k(0,2) = -sixEI;
      k(2,2) = fourEI;
      k(3,2) = sixEI;
      k(5,2) = twoEI;
      
      k(0,3) = -twelveEI;
      k(2,3) = sixEI;
      k(3,3) = twelveEI;
      k(5,3) = sixEI;
      
      k(1,4) = -oneEA;
      k(4,4) = oneEA;
      
      k(0,5) = -sixEI;
      k(2,5) = twoEI;
      k(3,5) = sixEI;
      k(5,5) = fourEI;
    }
    else {
      k(0,0) = oneEA;
      k(3,0) = -oneEA;
      
      k(1,1) = twelveEI;
      k(2,1)= sixEI;
      k(4,1) = -twelveEI;
      k(5,1) = sixEI;
      
      k(1,2) = sixEI;
      k(2,2) = fourEI;
      k(4,2) = -sixEI;
      k(5,2) = twoEI;
      
      k(0,3) = -oneEA;
      k(3,3) = oneEA;
      
      k(1,4) = -twelveEI;
      k(2,4) = -sixEI;
      k(4,4) = twelveEI;
      k(5,4)  = -sixEI;
      
      k(1,5) = sixEI;
      k(2,5) = twoEI;
      k(4,5) = -sixEI;
      k(5,5) = fourEI;
      
      if (cs != 1.0)
	k = trans^ k * trans;	    
    }
}



const Matrix &
beam2d03::getTangentStiff(void)
{
    return k;
}


const Matrix &
beam2d03::getInitialStiff(void)
{
    return k;
}
    
void 
beam2d03::zeroLoad(void)
{
    load.Zero();
}

int 
beam2d03::addLoad(ElementalLoad *theLoad, double loadFactor)
{

  opserr << "beam2d03::addLoad() - beam " << this->getTag() << "load type unknown\n";
  return -1;
    
  return 0;
}

int
beam2d03::addInertiaLoadToUnbalance(const Vector &accel)
{
  return 0;
}

const Vector &
beam2d03::getResistingForce()
{	
    // compute the residual Res = k*uTrial
    const Vector &end1Disp = theNodes[0]->getTrialDisp();
    const Vector &end2Disp = theNodes[1]->getTrialDisp();    
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

const Vector &
beam2d03::getResistingForceIncInertia()
{	
    this->getResistingForce();  

    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	rForce += this->getRayleighDampingForces();

    return rForce;
}


int
beam2d03::sendSelf(int commitTag, Channel &theChannel)
{
  int dataTag = this->getDbTag();

  Vector data(4);
  data(0) = A; data(1) = E; data(2) = I; data(3) = this->getTag();
    
  int result = 0;
  result = theChannel.sendVector(dataTag, commitTag, data);
  if (result < 0) {
    opserr << "beam2d03::sendSelf - failed to send data\n";
    return -1;
  }
  
  result = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (result < 0) {
    opserr << "beam2d03::sendSelf - failed to send data\n";
    return -1;
  }
    
  return 0;
}

int
beam2d03::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    int result = 0;
    int dataTag = this->getDbTag();

    result = theChannel.recvVector(dataTag, commitTag, data);
    if (result < 0) {
	opserr << "beam2d03::recvSelf - failed to recv data\n";
	return -1;
    }

    A = data(0); E = data(1); I=data(2); 
    this->setTag((int)data(3));

    result = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
    if (result < 0) {
	opserr << "beam2d03::recvSelf - failed to recv data\n";
	return -1;
    }
    
    return 0;
}



void
beam2d03::Print(OPS_Stream &s, int flag)
{
  //    s << "\nElement: " << this->getTag() << " Type: beam2d03 ";
  //    s << "\tConnected Nodes: " << connectedExternalNodes ;
//    s << "\tStiffness Matrix:\n" << k;
  //    s << "\tResisting Force: " << rForce;
//    s << "\tElemt End Force: " << eForce;    
}




