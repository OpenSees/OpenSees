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
// $Source: /usr/local/cvs/OpenSees/SRC/element/beam2d/beam2d02.cpp,v $
                                                                        
                                                                        
// File: ~/element/beam2d02.C
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the implementation for the beam2d02 class.

#include "beam2d02.h"
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <CrdTransf2d.h>

beam2d02::beam2d02()
   :Element(0,ELE_TAG_beam2d02), 
    A(0.0), E(0.0), I(0.0), M(0.0), L(0),
    connectedExternalNodes(2),
    Kd(3,3), m(6,6),
    q(3), rForce(6), load(6), 
    theCoordTrans(0)
{

}

// beam2d02(int tag, double A, double E, double I, int Nd1, int Nd2);
//	constructor which takes the unique element tag, the elements A,E and
//	I and the node ID's of it's nodal end points. 
beam2d02::beam2d02(int tag, double a, double e, double i, int Nd1, int Nd2, 
		   CrdTransf2d &theTrans, double rho)
   :Element(tag,ELE_TAG_beam2d02), 
    A(a), E(e), I(i), M(rho), L(0),
    connectedExternalNodes(2),
    Kd(3,3), m(6,6), 
    q(3), rForce(6), load(6), 
    theCoordTrans(0)    
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;    

    theNodes[0] = 0;
    theNodes[1] = 0;

    theCoordTrans = theTrans.getCopy();
}



// ~beam2d02():
// 	destructor

beam2d02::~beam2d02()
{
    if (theCoordTrans != 0)
	delete theCoordTrans;
}


int
beam2d02::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
beam2d02::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
beam2d02::getNodePtrs(void) 
{
  return theNodes;
}

int
beam2d02::getNumDOF(void) {
    return 6;
}


void
beam2d02::setDomain(Domain *theDomain)
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

    // determine length and direction cosines
    double dx,dy;
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    
    dx = end2Crd(0)-end1Crd(0);
    dy = end2Crd(1)-end1Crd(1);	
    
    L = sqrt(dx*dx + dy*dy);
    if (L == 0.0) {
	opserr << "WARNING beam2d02::setDomain(): beam " << this->getTag();
	opserr << " has zero length for beam\n" << *this;
	return;  
    }

    Kd(0,0) = E*A/L;
    Kd(0,1) = 0.0;
    Kd(0,2) = 0.0;    

    Kd(1,0) = 0.0;        
    Kd(1,1) = 4.0*E*I/L;
    Kd(1,2) = 2.0*E*I/L;

    Kd(2,0) = 0.0;    
    Kd(2,1) = 2.0*E*I/L;
    Kd(2,2) = 4.0*E*I/L;

    cs = dx/L;
    sn = dy/L;
    
    // set the mass variable equal to 1/2 the mass of the beam = 0.5 * rho*A*L
    M = 0.5*M*A*L;    
}


int
beam2d02::commitState()
{
  int retVal = 0;
  
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "beam2d02::commitState () - failed in base class";
  }    

  retVal = theCoordTrans->commitState();
  return retVal;
}

int
beam2d02::revertToLastCommit()
{
    return theCoordTrans->revertToLastCommit();    
}

int
beam2d02::revertToStart()
{
    return theCoordTrans->revertToStart();    
}

const Matrix &
beam2d02::getTangentStiff(void)
{
    return this->getStiff();
}

const Matrix &
beam2d02::getInitialStiff(void)
{
    return this->getStiff();
}


// const Matrix &getStiff():
//	Method to return the stiffness matrix.

const Matrix &
beam2d02::getStiff(void)
{
    const Vector &v = theCoordTrans->getBasicTrialDisp();
    q.addMatrixVector(0.0,Kd,v,1.0);
    
    return theCoordTrans->getGlobalStiffMatrix(Kd, q);
    //    return k;
}
    
const Matrix &
beam2d02::getMass(void)
{ 
    // lumped mass matrix
    m(0,0) = M;
    m(1,1) = M;
    m(3,3) = M;
    m(4,4) = M;
    return m;
}



void 
beam2d02::zeroLoad(void)
{
    load.Zero();
}


int 
beam2d02::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  opserr << "ElasticBeam2d::addLoad() - beam " << this->getTag() << ", does not handle ele loads\n";
  return -1;
}


int
beam2d02::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (M == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
	
  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  load(0) -= M * Raccel1(0);
  load(1) -= M * Raccel1(1);
    
  load(3) -= M * Raccel2(0);    
  load(4) -= M * Raccel2(1);    

  return 0;
}


const Vector &
beam2d02::getResistingForce()
{	
    // compute the residual Res = k*uTrial
    const Vector &v = theCoordTrans->getBasicTrialDisp();
    q.addMatrixVector(0.0,Kd,v,1.0);

    static Vector uniLoad(2);

    rForce = theCoordTrans->getGlobalResistingForce(q, uniLoad);
    
    // add any applied load
    rForce -= load;
    return rForce;
}


const Vector &
beam2d02::getResistingForceIncInertia()
{	
    this->getResistingForce();     
    
    if (M != 0.0) {
      const Vector &end1Accel = theNodes[0]->getTrialAccel();
      const Vector &end2Accel = theNodes[1]->getTrialAccel();    
      Vector inertia(6);
      rForce(0) += M*end1Accel(0);
      rForce(1) += M*end1Accel(1);
      rForce(3) += M*end2Accel(0);
      rForce(4) += M*end2Accel(1);

      // add rayleigh damping forces
      if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	rForce += this->getRayleighDampingForces();

    } else {

      // add rayleigh damping forces
      if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
	rForce += this->getRayleighDampingForces();

    }

    return rForce;    
}



int
beam2d02::sendSelf(int commitTag, Channel &theChannel)
{
  int dataTag = this->getDbTag();

  Vector data(4);
  data(0) = A; data(1) = E; data(2) = I; data(3) = this->getTag();
    
  int result = 0;
  result = theChannel.sendVector(dataTag, commitTag, data);
  if (result < 0) {
    opserr << "beam2d02::sendSelf - failed to send data\n";
    return -1;
  }
  
  result = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  if (result < 0) {
    opserr << "beam2d02::sendSelf - failed to send data\n";
    return -1;
  }
    
  return 0;
}

int
beam2d02::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(4);
    int result = 0;
    int dataTag = this->getDbTag();

    result = theChannel.recvVector(dataTag, commitTag, data);
    if (result < 0) {
	opserr << "beam2d02::recvSelf - failed to recv data\n";
	return -1;
    }

    A = data(0); E = data(1); I=data(2); 
    this->setTag((int)data(3));

    result = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
    if (result < 0) {
	opserr << "beam2d02::recvSelf - failed to recv data\n";
	return -1;
    }
    
    return 0;
}

int
beam2d02::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
    // first determine the two end points of the truss based on
    // the display factor (a measure of the distorted image)
    // store this information in 2 3d vectors v1 and v2
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();    

    if (displayMode == 1 || displayMode == 2) {
	Vector v1(3);
	Vector v2(3);
	for (int i=0; i<2; i++) {
	    v1(i) = end1Crd(i)+end1Disp(i)*fact;
	    v2(i) = end2Crd(i)+end2Disp(i)*fact;    
	}
	
	return theViewer.drawLine(v1, v2, 1.0,1.0);	
    } else 
	return 0;
}

void
beam2d02::Print(OPS_Stream &s, int flag)
{
    // compute current state
    this->getResistingForce();
    
    s << "Element: " << this->getTag(); 
    s << " type: beam2d02  iNode: " << connectedExternalNodes(0);    
    s << " jNode: " << connectedExternalNodes(1);
    s << " Area: " << A << " E: " << E << " I: " << I << endln;
    s << " resisting Force: " << rForce;
}







