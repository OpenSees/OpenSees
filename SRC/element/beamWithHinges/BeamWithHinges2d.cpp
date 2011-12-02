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
// $Source: /usr/local/cvs/OpenSees/SRC/element/beamWithHinges/BeamWithHinges2d.cpp,v $
                                                                        
                                                                        
///////////////////////////////////////////////////////////
// File:  ~/Src/element/beamWithHinges/BeamWithHinges2d.cpp
//
// Written by Matthew Peavy
//
// Written:	 Feb 13, 2000
// Debugged: Feb 15, 2000
// Revised:        , 200x
//
//
// Purpose:  This cpp file contains the definitions
// for BeamWithHinges.  The element uses hinges (passed as
// object pointers) at the element ends.  The middle section
// is analyzed elastically, while the hinge ends return
// a section flexibility at their middle points.  The element
// calculates a frame (6x6) stiffness matrix from this info.

#include <BeamWithHinges2d.h>
#include <Element.h>
#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Matrix.h>
#include <Vector.h>						
#include <Node.h>
#include <MatrixUtil.h>						//For 3x3 matrix inversion
#include <math.h>							//For sqrt function
#include <stdlib.h>						
#include <iostream.h>						//Not needed
#include <string.h>
#include <stdio.h>

#include <SectionForceDeformation.h>
#include <CrdTransf2d.h>

#include <Information.h>

BeamWithHinges2d::BeamWithHinges2d(void)
	:Element(0, ELE_TAG_BeamWithHinges2d),
	E(0), I(0), A(0), L(0),
	G(0), alpha(0),
	sectionI(0), sectionJ(0),
	connectedExternalNodes(2),
	node1Ptr(0), node2Ptr(0),
	hingeIlen(0), hingeJlen(0),
	K(6,6),	m(6,6),	d(6,6), fElastic(3,3), vElastic(3,2),
	UePrev(6), P(6), Pinert(6), kb(3,3), q(3),
	load(6), prevDistrLoad(2), 
	distrLoadCommit(2),	UeCommit(6), 
	kbCommit(3,3), qCommit(3), massDens(0),
	shearLength(0), shearIkey(-1), shearJkey(-1),
	theCoordTransf(0), maxIter(1), tolerance(1.0e-10), initialFlag(false)
{

}

BeamWithHinges2d::BeamWithHinges2d(int tag, int nodeI, int nodeJ,
						double E_in, double I_in, double A_in,
						double G_in, double alpha_in,
						SectionForceDeformation &sectionRefI, double I_length,
						SectionForceDeformation &sectionRefJ, double J_length,
						CrdTransf2d &coordTransf, double shearL,
						double massDensityPerUnitLength, int max, double tol)
	:Element(tag, ELE_TAG_BeamWithHinges2d),
	E(E_in), I(I_in), A(A_in),
	G(G_in), alpha(alpha_in),
	sectionI(0), sectionJ(0),
	connectedExternalNodes(2),
	node1Ptr(0), node2Ptr(0),
	hingeIlen(I_length), hingeJlen(J_length),
	K(6,6), m(6,6), d(6,6), fElastic(3,3), vElastic(3,2),
	UePrev(6), P(6), Pinert(6), kb(3,3), q(3),
	load(6), prevDistrLoad(2), 
	distrLoadCommit(2),	UeCommit(6), 
	kbCommit(3,3), qCommit(3),
	massDens(massDensityPerUnitLength),
	shearLength(shearL), shearIkey(-1), shearJkey(-1),
	theCoordTransf(0), maxIter(max), tolerance(tol), initialFlag(false)
{
	if (E <= 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Paramater E is zero or negative.";
		//exit(-1);
	}
	
	if (I <= 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Parameter I is zero or negative.";
		//exit(-1);
	}
	
	if (A <= 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Parameter A is zero or negative.";
		//exit(-1); 
	}

	if (G <= 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Parameter G is zero or negative.";
		//exit(-1); 
	}

	if (alpha < 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Parameter alpha is negative.";
		//exit(-1); 
	}

	// Get copies of sections
	sectionI = sectionRefI.getCopy();
	
	if (sectionI == 0)
		g3ErrorHandler->fatal("%s -- failed to get copy of section I",
			"BeamWithHinges2d::BeamWithHinges2d");

	sectionJ = sectionRefJ.getCopy();

	if (sectionJ == 0)
		g3ErrorHandler->fatal("%s -- failed to get copy of section J",
			"BeamWithHinges2d::BeamWithHinges2d");

	theCoordTransf = coordTransf.getCopy();
	
	if (theCoordTransf == 0)
	    g3ErrorHandler->fatal("%s -- failed to get copy of coordinate transformation",
			"BeamWithHinges2d::BeamWithHinges2d");
	
	connectedExternalNodes(0) = nodeI;
	connectedExternalNodes(1) = nodeJ;
}

BeamWithHinges2d::BeamWithHinges2d(int tag, int nodeI, int nodeJ,
				   double E_in, double I_in, double A_in,
				   double G_in, double alpha_in,
				   SectionForceDeformation &sectionRefI, double I_length,
				   SectionForceDeformation &sectionRefJ, double J_length,
				   CrdTransf2d &coordTransf, const Vector &distrLoad,
				   double shearL, double massDensityPerUnitLength,
				   int max, double tol)
	:Element(tag, ELE_TAG_BeamWithHinges2d),
	E(E_in), I(I_in), A(A_in),
	G(G_in), alpha(alpha_in),
	sectionI(0), sectionJ(0),
	connectedExternalNodes(2),
	node1Ptr(0), node2Ptr(0),
	hingeIlen(I_length), hingeJlen(J_length),
	K(6,6), m(6,6), d(6,6), fElastic(3,3), vElastic(3,2),
	UePrev(6), P(6), Pinert(6), kb(3,3), q(3),
	load(6), prevDistrLoad(2), 
	distrLoadCommit(2),	UeCommit(6), 
	kbCommit(3,3), qCommit(3),
	massDens(massDensityPerUnitLength),
	shearLength(shearL), shearIkey(-1), shearJkey(-1),
	theCoordTransf(0), maxIter(max), tolerance(tol), initialFlag(false)
{
	if (E <= 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Paramater E is zero or negative.";
		//exit(-1);
	}
	
	if (I <= 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Parameter I is zero or negative.";
		//exit(-1);
	}
	
	if (A <= 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Parameter A is zero or negative.";
		//exit(-1); 
	}

	if (G <= 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Parameter G is zero or negative.";
		//exit(-1); 
	}

	if (alpha < 0.0)  {
		cerr << "FATAL - BeamWithHinges2d::BeamWithHinges2d() - Parameter alpha is negative.";
		//exit(-1); 
	}

	// Get copies of sections
	sectionI = sectionRefI.getCopy();
	
	if (sectionI == 0)
		g3ErrorHandler->fatal("%s -- failed to get copy of section I",
			"BeamWithHinges2d::BeamWithHinges2d");

	sectionJ = sectionRefJ.getCopy();

	if (sectionJ == 0)
		g3ErrorHandler->fatal("%s -- failed to get copy of section J",
			"BeamWithHinges2d::BeamWithHinges2d");

	theCoordTransf = coordTransf.getCopy();
	
	if (theCoordTransf == 0)
	    g3ErrorHandler->fatal("%s -- failed to get copy of coordinate transformation",
			"BeamWithHinges2d::BeamWithHinges2d");
	
	connectedExternalNodes(0) = nodeI;					//node tags
	connectedExternalNodes(1) = nodeJ;
}

BeamWithHinges2d::~BeamWithHinges2d(void)
{
	//Only hinges must be destroyed
	if (sectionI)
		delete sectionI;

	if (sectionJ)
		delete sectionJ;
	
	if (theCoordTransf)
	    delete theCoordTransf;
}

int 
BeamWithHinges2d::getNumExternalNodes(void) const
{
	return 2;
}

const ID &
BeamWithHinges2d::getExternalNodes(void)
{
	return connectedExternalNodes;
}

int 
BeamWithHinges2d::getNumDOF(void)
{
	return 6;
}

void
BeamWithHinges2d::setDomain (Domain *theDomain)
{
	//This function may be called after a beam is constructed, so
	//geometry may change.  Therefore calculate all beam geometry here.

	if(theDomain == 0) {
		node1Ptr = 0;
		node2Ptr = 0;
		return;
	}
	
	// set the node pointers and verify them
	this->setNodePtrs(theDomain);

	// call the DomainComponent version of the function
	this->DomainComponent::setDomain(theDomain);

	if (theCoordTransf->initialize(node1Ptr, node2Ptr) != 0) {
	    cerr << "BeamWithHinges2d::setDomain(): Error initializing coordinate transformation";  
	    exit(-1);
	}
	
	// get element length
	L = theCoordTransf->getInitialLength();
	
	if (L == 0.0) {
	    cerr << "BeamWithHinges2d::setDomain(): Zero element length:" << this->getTag();  
	    exit(-1);
	}
	
	UePrev.Zero();					//clear the previous nodal displacements
	initialFlag = false;			//used in getTangentStiff to force initialization

	// Set up section interpolation and hinge lengths
	this->setHinges();

	// Calculate the elastic contribution to the element flexibility
	this->setElasticFlex();

	// Set the element lumped mass matrix
	this->setMass();
}

int
BeamWithHinges2d::commitState(void)
{
	int err = 0;
	
	err += sectionI->commitState();
	err += sectionJ->commitState();

	err += theCoordTransf->commitState();
	
	distrLoadCommit = prevDistrLoad;
	UeCommit = UePrev;
	kbCommit = kb;
	qCommit = q;

	initialFlag = false;
	
	return err;
}

int
BeamWithHinges2d::revertToLastCommit(void)
{
	int err = 0;

	// Revert the sections and then get their last commited
	// deformations, stress resultants, and flexibilities
	err += sectionI->revertToLastCommit();
	e1 = sectionI->getSectionDeformation();
	sr1 = sectionI->getStressResultant();
	fs1 = sectionI->getSectionFlexibility();

	err += sectionJ->revertToLastCommit();
	e3 = sectionJ->getSectionDeformation();
	sr3 = sectionJ->getStressResultant();
	fs3 = sectionJ->getSectionFlexibility();

	// Commit the coordinate transformation
	err += theCoordTransf->revertToLastCommit();

	prevDistrLoad = distrLoadCommit;
	UePrev = UeCommit;
	kb = kbCommit;
	q = qCommit;

	initialFlag = false;

	return err;
}

int
BeamWithHinges2d::revertToStart(void)
{
    int err = 0;
    
    err += sectionI->revertToStart();
    err += sectionJ->revertToStart();    
    
    err += theCoordTransf->revertToStart();
    
    fs1.Zero();
    fs3.Zero();    
    
    e1.Zero();
    e3.Zero();	

	sr1.Zero();
	sr3.Zero();
    
    prevDistrLoad.Zero();
    UePrev.Zero();
    kb.Zero();
    q.Zero();

	return err;
}

const Matrix &
BeamWithHinges2d::getTangentStiff(void)
{
	// Update element state
	this->setStiffMatrix();
	//return K;

	// Return transformed basic stiffness
	return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &
BeamWithHinges2d::getSecantStiff(void)
{
	return this->getTangentStiff();
}

const Matrix &
BeamWithHinges2d::getDamp(void)
{
	return d;			//zero matrix
}

const Matrix &
BeamWithHinges2d::getMass(void)
{
	// Computed in setDomain()
	return m;
}

void 
BeamWithHinges2d::zeroLoad(void)
{
	load.Zero();
}

int
BeamWithHinges2d::addLoad(const Vector &moreLoad)
{
	if(moreLoad.Size() != 6) {
		cerr << "BeamWithHinges2d::addLoad - Vector not of correct size.\n";
		return -1;
	}

	//load += moreLoad;
	load.addVector(1.0, moreLoad, 1.0);
	
	return 0;
}

int
BeamWithHinges2d::addInertiaLoadToUnbalance(const Vector &accel)
{
	if (massDens == 0.0)
		return 0;

	static Vector ag(6);

	const Vector &Raccel1 = node1Ptr->getRV(accel);
	const Vector &Raccel2 = node2Ptr->getRV(accel);

	int i,j;
	for (i = 0, j = 3; i < 2; i++, j++) {
		load(i) += m(i,i) * Raccel1(i);
		load(j) += m(j,j) * Raccel2(i);	// Yes, this should be 'i'
	}
	
	return 0;
}

const Vector &
BeamWithHinges2d::getResistingForce(void)
{
	this->setStiffMatrix();
	//return P;
	static Vector dum(2);
	P = theCoordTransf->getGlobalResistingForce(q, dum);

	P.addVector(1.0, load, -1.0);

	return P;
}

const Vector &
BeamWithHinges2d::getResistingForceIncInertia(void)
{
	if (massDens == 0.0)
		return this->getResistingForce();

	static Vector ag(6);

	this->getGlobalAccels(ag);

	Pinert = this->getResistingForce();

	int i,j;	
	for (i = 0, j = 3; i < 2; i++, j++) {
		Pinert(i) += m(i,i) * ag(i);
		Pinert(j) += m(j,j) * ag(j);
	}

	return Pinert;
}

int
BeamWithHinges2d::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	int elemDbTag = this->getDbTag();

	static ID idData(13);
	
	int orderI = sectionI->getOrder();
	int orderJ = sectionJ->getOrder();

	idData(0) = this->getTag();
	idData(1) = sectionI->getClassTag();
	idData(2) = sectionJ->getClassTag();
	idData(5) = theCoordTransf->getClassTag();
	idData(7) = orderI;
	idData(8) = orderJ;
	idData(9) = maxIter;
	idData(10) = (initialFlag) ? 1 : 0;
	idData(11) = connectedExternalNodes(0);
	idData(12) = connectedExternalNodes(1);

	int sectDbTag;

	sectDbTag = sectionI->getDbTag();
    if (sectDbTag == 0) {
		sectDbTag = theChannel.getDbTag();
		sectionI->setDbTag(sectDbTag);
    }
    idData(3) = sectDbTag;

	sectDbTag = sectionJ->getDbTag();
    if (sectDbTag == 0) {
		sectDbTag = theChannel.getDbTag();
		sectionJ->setDbTag(sectDbTag);
    }
    idData(4) = sectDbTag;

	sectDbTag = theCoordTransf->getDbTag();
    if (sectDbTag == 0) {
		sectDbTag = theChannel.getDbTag();
		theCoordTransf->setDbTag(sectDbTag);
    }
    idData(6) = sectDbTag;

	// Send the data ID
	res += theChannel.sendID(elemDbTag, commitTag, idData);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- failed to send data ID",
			"BeamWithHinges2d::sendSelf");
		return res;
	}

	// Ask the sectionI to send itself
	res += sectionI->sendSelf(commitTag, theChannel);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- failed to send section I",
			"BeamWithHinges2d::sendSelf");
		return res;
	}

	// Ask the sectionJ to send itself
	res += sectionJ->sendSelf(commitTag, theChannel);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- failed to send section J",
			"BeamWithHinges2d::sendSelf");
		return res;
	}

	// data, kbcommit, qcommit, Uecommit, distrLoadcommit
	static Vector data(11 + 3*3 + 3 + 6 + 2);
	
	data(0) = L;
	data(1) = E;  
	data(2) = A;
	data(3) = I;  
	data(4) = G;
	data(5) = alpha;
	data(6) = massDens;
	data(7) = hingeIlen/L;		// convert back to ratios
	data(8) = hingeJlen/L;
	data(9) = tolerance;
	data(10) = shearLength/L;

	// Jam committed history variables into the data vector
	int i,j;
	int loc = 11;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			data(loc++) = kbCommit(i,j);

	for (i = 0; i < 3; i++)
		data(loc++) = qCommit(i);

	for (i = 0; i < 6; i++)
		data(loc++) = UeCommit(i);

	for (i = 0; i < 2; i++)
		data(loc++) = distrLoadCommit(i);

	// Send the data vector
	res += theChannel.sendVector(elemDbTag, commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- failed to send tags ID",
			"BeamWithHinges2d::sendSelf");
		return res;
	}

	// Ask the coordinate transformation to send itself
	res += theCoordTransf->sendSelf(commitTag, theChannel);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- failed to send coordinate transformation",
			"BeamWithHinges2d::sendSelf");
		return res;
	}

	return res;
}

int 
BeamWithHinges2d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	int elemDbTag = this->getDbTag();

	static ID idData(13);

	// Receive the ID data
	res += theChannel.recvID(elemDbTag, commitTag, idData);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- failed to received data ID",
			"BeamWithHinges2d::recvSelf");
		return res;
	}

	this->setTag(idData(0));
	maxIter = idData(9);
	initialFlag = (idData(10) == 1) ? true : false;
	connectedExternalNodes(0) = idData(11);
	connectedExternalNodes(1) = idData(12);

	int secTag;

	// Receive section I
	secTag = idData(1);

	if (sectionI == 0)
		sectionI = theBroker.getNewSection(secTag);

	// Check that the section is of the right type; if not, delete
	// the current one and get a new one of the right type
	else if (sectionI->getClassTag() != secTag) {
		delete sectionI;
		sectionI = theBroker.getNewSection(secTag);
	}

	if (sectionI == 0) {
		g3ErrorHandler->warning("%s -- could not get a Section",
			"BeamWithHinges2d::recvSelf");
		return -1;
	}

	// Now, receive the section
	sectionI->setDbTag(idData(3));
	res += sectionI->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Section I",
			"BeamWithHinges2d::recvSelf");
		return res;
	}

	// Receive section J
	secTag = idData(2);

	if (sectionJ == 0)
		sectionJ = theBroker.getNewSection(secTag);

	// Check that the section is of the right type; if not, delete
	// the current one and get a new one of the right type
	else if (sectionJ->getClassTag() != secTag) {
		delete sectionJ;
		sectionJ = theBroker.getNewSection(secTag);
	}

	if (sectionJ == 0) {
		g3ErrorHandler->warning("%s -- could not get a Section",
			"BeamWithHinges2d::recvSelf");
		return -1;
	}

	// Now, receive the section
	sectionJ->setDbTag(idData(4));
	res += sectionJ->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Section J",
			"BeamWithHinges2d::recvSelf");
		return res;
	}

	int orderI = idData(7);
	int orderJ = idData(8);

	static Vector data(11 + 3*3 + 3 + 6 + 2);

	res += theChannel.recvVector(elemDbTag, commitTag, data);
	if(res < 0) {
		g3ErrorHandler->warning("%s -- failed to receive Vector data",
			"BeamWithHinges2d::recvSelf");
		return res;
	}

	L = data(0);
	E = data(1);
	A = data(2);
	I = data(3);
	G = data(4);
	alpha = data(5);
	massDens = data(6);
	hingeIlen = data(7);
	hingeJlen = data(8);
	tolerance = data(9);
	shearLength = data(10);

	// Set up the section force interpolation
	this->setHinges();

	// Compute the elastic flexibility
	this->setElasticFlex();

	// Get committed history variables from the data vector
	int i,j;
	int loc = 11;
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			kbCommit(i,j) = data(loc++);

	for (i = 0; i < 3; i++)
		qCommit(i) = data(loc++);

	for (i = 0; i < 6; i++)
		UeCommit(i) = data(loc++);

	for (i = 0; i < 2; i++)
		distrLoadCommit(i) = data(loc++);

	// Receive the coordinate transformation
	int crdTag = idData(5);

	// Allocate new CoordTransf if null
	if (theCoordTransf == 0)
		theCoordTransf = theBroker.getNewCrdTransf2d(crdTag);

	// Check that the CoordTransf is of the right type; if not, delete
	// the current one and get a new one of the right type
	else if (theCoordTransf->getClassTag() != crdTag) {
		delete theCoordTransf;
		theCoordTransf = theBroker.getNewCrdTransf2d(crdTag);
	}

	// Check if either allocation failed
	if (theCoordTransf == 0) {
		g3ErrorHandler->warning("%s -- could not get a CrdTransf3d",
			"BeamWithHinges2d::recvSelf");
		return -1;
	}

	// Now, receive the CoordTransf
	theCoordTransf->setDbTag(idData(6));
	res += theCoordTransf->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive CoordTransf",
			"BeamWithHinges2d::recvSelf");
		return res;
	}
	
	// Revert the element to the last committed state that was just received
	// Can't call this->revertToLastCommit() because this->setDomain() may
	// not have been called yet
	sectionI->revertToLastCommit();
	e1 = sectionI->getSectionDeformation();
	sectionI->setTrialSectionDeformation(e1);
	sr1 = sectionI->getStressResultant();
	fs1 = sectionI->getSectionFlexibility();

	sectionJ->revertToLastCommit();
	e3 = sectionJ->getSectionDeformation();
	sectionJ->setTrialSectionDeformation(e3);
	sr3 = sectionJ->getStressResultant();
	fs3 = sectionJ->getSectionFlexibility();

	prevDistrLoad = distrLoadCommit;
	UePrev = UeCommit;
	kb = kbCommit;
	q = qCommit;

	initialFlag = false;

	return res;
}

void 
BeamWithHinges2d::Print(ostream &s, int flag)
{
	s << "\nBeamWithHinges2d, tag: " << this->getTag() << endl;
	s << "\tConnected Nodes: " << connectedExternalNodes;
	s << "\tE: " << E << endl;
	s << "\tA: " << A << endl;
	s << "\tI: " << I << endl;

	if (sectionI) {
		s << "\tHinge I, section tag: " << sectionI->getTag() << 
			", length: " << hingeIlen << endl;
		sectionI->Print(s,flag);
	}

	if (sectionJ)
		s << "\tHinge J, section tag: " << sectionJ->getTag() << 
			", length: " << hingeJlen << endl;
}

//////////////////////////////
//Private Function Definitions

void 
BeamWithHinges2d::setNodePtrs(Domain *theDomain)
{
	node1Ptr = theDomain->getNode(connectedExternalNodes(0));
	node2Ptr = theDomain->getNode(connectedExternalNodes(1));

	if(node1Ptr == 0) {
		cerr << "ERROR: BeamWithHinges2d::setNodePtrs() - Node 1: ";
		cerr << connectedExternalNodes(0) << " does not exist.\n";
		exit(-1);
	}
	
	if(node2Ptr == 0) {
		cerr << "ERROR: BeamWithHinges2d::setNodePtrs() - Node 2: ";
		cerr << connectedExternalNodes(1) << " does not exist.\n";
		exit(-1);
	}
	
	// check for correct # of DOF's
	int dofNd1 = node1Ptr->getNumberDOF();
	int dofNd2 = node2Ptr->getNumberDOF();
	if ((dofNd1 != 3) || (dofNd2 != 3))  {
		cerr << "ERROR: BeamWithHinges2d::setNodePtrs() - Number of ";
		cerr << "DOF's in nodes is incorrect.\n";
		exit(-1);
	}
}

void
BeamWithHinges2d::getGlobalDispls(Vector &dg)    //get nodal displacements
{
    const Vector &disp1 = node1Ptr->getTrialDisp();
    const Vector &disp2 = node2Ptr->getTrialDisp();

    dg(0) = disp1(0);		//left side
    dg(1) = disp1(1);
    dg(2) = disp1(2);
    dg(3) = disp2(0);		//right side
    dg(4) = disp2(1);
    dg(5) = disp2(2);
}

void
BeamWithHinges2d::getGlobalAccels(Vector &ag)	//for dynamic analysis
{
    const Vector &accel1 = node1Ptr->getTrialAccel();
    const Vector &accel2 = node2Ptr->getTrialAccel();

    ag(0) = accel1(0);
    ag(1) = accel1(1);
    ag(2) = accel1(2);
    ag(3) = accel2(0);
    ag(4) = accel2(1);
    ag(5) = accel2(2);
}

void
BeamWithHinges2d::setElasticFlex ()
{
	//calculates the middle beam flexibility matrix
	//this is obtained by integrating (from lp1 to L-lp2) the matrix
	//b ^ fs ^ b, where fs is a constant linear section flexibility matrix.
	
	fElastic.Zero();
	vElastic.Zero();

	double a = hingeIlen;
	double b = L-hingeJlen;
	
	// Integrate elastic portion numerically 
	// Create elastic section
	static Matrix fElas(2,2);
	fElas(0,0) = 1.0/(E*A);
	fElas(1,1) = 1.0/(E*I);

	// Section code passed to force interpolation
	static ID code(2);
	code(0) = SECTION_RESPONSE_P;
	code(1) = SECTION_RESPONSE_MZ;

	// Section force interpolation matrix
	static Matrix bElas(2,3);
	
	// Mapping from [-1,1] to [a,b]
	double alpha = 0.5*(b-a);	// Jacobian
	double beta = 0.5*(a+b);

	// Integration points ... the weights are 1.0*(0.5*(b-a)) ==> alpha
	const double xsi = 1.0/sqrt(3.0);
	double x1 = alpha*(-xsi) + beta;
	double x2 = alpha*xsi + beta;

	int dum;	
	// Get interpolations at sample points and integrate
	this->getForceInterpMatrix(bElas, x1, code, dum);
	fElastic.addMatrixTripleProduct(1.0, bElas, fElas, alpha);
	vElastic.addMatrix(1.0, bElas^ fElas, alpha);

	this->getForceInterpMatrix(bElas, x2, code, dum);
	fElastic.addMatrixTripleProduct(1.0, bElas, fElas, alpha);
	vElastic.addMatrix(1.0, bElas^ fElas, alpha);

	/*
	double shearFlex = 0.0;

	if (alpha > 0.0)
		shearFlex = 1/((G*A*alpha)*L*L) * (b-a);

	fElastic(0,0) = 1/(E*A) * ( b-a );
	fElastic(1,1) = 1/(E*I) * ( 1/(3*L*L)*(b*b*b-a*a*a) - 1/L*(b*b-a*a) + b-a ) + shearFlex;
	fElastic(2,2) = 1/(E*I) * ( 1/(3*L*L)*(b*b*b-a*a*a) ) + shearFlex;
	fElastic(1,2) = 1/(E*I) * ( 1/(3*L*L)*(b*b*b-a*a*a) - 1/(2*L)*(b*b-a*a) ) + shearFlex;
	fElastic(2,1) = fElastic(1,2);

	vElastic(0,0) = 1/(E*A) * ( b-a - 1/(2*L)*(b*b-a*a) );
	vElastic(1,1) = 1/(E*I) * ( 1/(4*L*L*L)*(b*b*b*b-a*a*a*a) - 1/(3*L*L)*(b*b*b-a*a*a) );

	if (alpha > 0.0)
		vElastic(2,1) = 1/(G*A*alpha) * ( 1/(2*L*L)*(b*b-a*a) - 1/(2*L)*(b-a) );
	*/
}

void
BeamWithHinges2d::setStiffMatrix(void)
{
    // Get element global end displacements
    static Vector Ue(6);
    this->getGlobalDispls(Ue);
    
    // Compute global end deformation increments
    static Vector dUe(6);
    //dUe = Ue - UePrev;
	dUe = Ue;
	dUe.addVector(1.0, UePrev, -1.0);
    
    // Check if global end displacements have been changed
    if (dUe.Norm() != 0.0 || initialFlag == false) {
	// Compute distributed loads and increments
	static Vector currDistrLoad(2);
	static Vector distrLoadIncr(2); 
	
	currDistrLoad.Zero();							
	//distrLoadIncr = currDistrLoad - prevDistrLoad;
	distrLoadIncr = currDistrLoad;
	distrLoadIncr.addVector(1.0, prevDistrLoad, -1.0);
	prevDistrLoad = currDistrLoad;
	
	// Update the end deformations
	UePrev = Ue;
	
	theCoordTransf->update();
	
	// Convert to basic system from local coord's (eleminate rb-modes)
	static Vector dv(3);
	static Vector v(3);				// basic system deformations
	
	v = theCoordTransf->getBasicTrialDisp();
	dv = theCoordTransf->getBasicIncrDeltaDisp();
	
	// calculate nodal force increments and update nodal forces
	static Vector dq(3);
	//dq = kb * dv;					// using previous stiff matrix k,i
	dq.addMatrixVector(0.0, kb, dv, 1.0);

	// Element properties
	static Matrix f(3,3);	// Element flexibility
	static Vector vr(3);	// Residual element deformations
	
	Vector s1(sectionI->getStressResultant());
	Vector ds1(s1);
	Vector de1(s1);
	
	Vector s3(sectionJ->getStressResultant());
	Vector ds3(s3);	
	Vector de3(s3);
	
	for (int j = 0; j < maxIter; j++) {
	    
	    //q += dq;
		q.addVector(1.0, dq, 1.0);
	    
	    // Section forces
	    //s1 = b1*q + bp1*currDistrLoad;
		s1.addMatrixVector(0.0, b1, q, 1.0);
		s1.addMatrixVector(1.0, bp1, currDistrLoad, 1.0);

	    //s3 = b3*q + bp3*currDistrLoad;
		s3.addMatrixVector(0.0, b3, q, 1.0);
		s3.addMatrixVector(1.0, bp3, currDistrLoad, 1.0);
	    
		// Increment in section forces
		// ds1 = s1 - sr1
		ds1 = s1;
		ds1.addVector(1.0, sr1, -1.0);

		// ds3 = s3 - sr3
		ds3 = s3;
		ds3.addVector(1.0, sr3, -1.0);

	    // compute section deformation increments and update current deformations
	    //e1 += fs1 * ds1;
		//e3 += fs3 * ds3;
		de1.addMatrixVector(0.0, fs1, ds1, 1.0);
		if (initialFlag != false)
			e1.addVector(1.0, de1, 1.0);
		
	    de3.addMatrixVector(0.0, fs3, ds3, 1.0);
		if (initialFlag != false)
			e3.addVector(1.0, de3, 1.0);
	    
	    /*** Dependent on section type ***/
	    // set section deformations
	    sectionI->setTrialSectionDeformation(e1);
	    sectionJ->setTrialSectionDeformation(e3);
	    
	    /*** Dependent on section type ***/
	    // get section resisting forces
	    sr1 = sectionI->getStressResultant();
	    sr3 = sectionJ->getStressResultant();
	    
	    /*** Dependent on section type ***/
	    // get section flexibility matrix
	    fs1 = sectionI->getSectionFlexibility();
	    fs3 = sectionJ->getSectionFlexibility();
	    
		// ds1 = s1 - sr1;
		// ds3 = s3 - sr3;
		ds1 = s1;
		ds1.addVector(1.0, sr1, -1.0);
		ds3 = s3;
		ds3.addVector(1.0, sr3, -1.0);

		de1.addMatrixVector(0.0, fs1, ds1, 1.0);
		de3.addMatrixVector(0.0, fs3, ds3, 1.0);
	    
	    // Use special integration on the diagonal shear flexibility term and
	    // the residual shear deformation
	    double shearWeight;

		if (hingeIlen > 0.0) {
			int orderI = sectionI->getOrder();
			const ID &codeI = sectionI->getType();
			shearWeight = shearLength/hingeIlen;
			for (int i = 0; i < orderI; i++) {
				if (codeI(i) == SECTION_RESPONSE_VY) {
					fs1(i,i) *= shearWeight;
					e1(i) *= shearWeight;
				}
			}
		}

		if (hingeJlen > 0.0) {
			int orderJ = sectionJ->getOrder();
			const ID &codeJ = sectionJ->getType();
			shearWeight = shearLength/hingeIlen;
			for (int i = 0; i < orderJ; i++) {
				if (codeJ(i) == SECTION_RESPONSE_VY) {
					fs3(i,i) *= shearWeight;
					e3(i) *= shearWeight;
				}
			}
		}

		f = fElastic;

	    // integrate section flexibility matrix
	    // Matrix f1 = (b1^ fs1 * b1) * hingeIlen;
		f.addMatrixTripleProduct(1.0, b1, fs1, hingeIlen);

	    // Matrix f3 = (b3^ fs3 * b3) * hingeJlen;
		f.addMatrixTripleProduct(1.0, b3, fs3, hingeJlen);

		// vr = fElastic * q;
		vr.addMatrixVector(0.0, fElastic, q, 1.0);

		// vr = vr + vElastic*currDistrLoad;
		vr.addMatrixVector(1.0, vElastic, currDistrLoad, 1.0);

	    // vr += vr1, where vr1 = (b1^ (e1+de1)) * hingeIlen;
		vr.addMatrixTransposeVector(1.0, b1, e1 + de1, hingeIlen);

	    // vr += vr3, where vr3 = (b3^ (e3+de3)) * hingeJlen;
	    vr.addMatrixTransposeVector(1.0, b3, e3 + de3, hingeJlen);
	    
	    // Undo the temporary change for integrating shear terms
		if (hingeIlen > 0.0) {
			int orderI = sectionI->getOrder();
			const ID &codeI = sectionI->getType();
			shearWeight = shearLength/hingeIlen;
			for (int i = 0; i < orderI; i++) {
				if (codeI(i) == SECTION_RESPONSE_VY) {
					fs1(i,i) /= shearWeight;
					e1(i) /= shearWeight;
				}
			}
		}

		if (hingeJlen > 0.0) {
			int orderJ = sectionJ->getOrder();
			const ID &codeJ = sectionJ->getType();
			shearWeight = shearLength/hingeIlen;
			for (int i = 0; i < orderJ; i++) {
				if (codeJ(i) == SECTION_RESPONSE_VY) {
					fs3(i,i) /= shearWeight;
					e3(i) /= shearWeight;
				}
			}
		}

		// calculate element stiffness matrix
	    invert3by3Matrix(f, kb);

	    //dv = v - vr;
		dv = v;
		dv.addVector(1.0, vr, -1.0);
	
	    // determine resisting forces
	    // dq = kb * dv;
		dq.addMatrixVector(0.0, kb, dv, 1.0);
	
	    double dW = dv^ dq;

	    if (dW < tolerance)
			break;
	}	
    
	// q += dq;
	q.addVector(1.0, dq, 1.0);
    
	//P = theCoordTransf->getGlobalResistingForce (q, currDistrLoad);
	//K = theCoordTransf->getGlobalStiffMatrix (kb, q);

	initialFlag = true;
    }

    return;
}

void
BeamWithHinges2d::setMass(void)
{
    //Note:  Lumped mass matrix only!
    //Note:  It is assumed that section's massDens is the same as beam's.
    m.Zero();
    m(0,0) = m(1,1) = m(3,3) = m(4,4) = massDens * L / 2.0;
}

void
BeamWithHinges2d::setHinges (void)
{
	// Get the number of section response quantities
	int orderI = sectionI->getOrder();
	int orderJ = sectionJ->getOrder();

	// Set interpolation matrices
	b1 = Matrix(orderI,3);
	b3 = Matrix(orderJ,3);

	// Set distributed load interpolation matrices
	bp1 = Matrix(orderI,2);
	bp3 = Matrix(orderJ,2);

	fs1 = Matrix(orderI,orderI);
	fs3 = Matrix(orderJ,orderJ);

	e1 = Vector(orderI);
	e3 = Vector(orderJ);
	
	sr1 = Vector(orderI);
	sr3 = Vector(orderJ);

	// Turn the hinge length ratios into actual lengths since L is not
	// known until setDomain() is called
	hingeIlen *= L;
	hingeJlen *= L;
	shearLength *= L;

	// Interpolation points
	double x1 = hingeIlen/2.0;
	double x3 = L - hingeJlen/2.0;

	// Get codes which indicate the ordering of section response quantities
	const ID &Icode = sectionI->getType();
	const ID &Jcode = sectionJ->getType();

	// Get the force interpolation matrix for each section
	this->getForceInterpMatrix (b1, x1, Icode, shearIkey);
	this->getForceInterpMatrix (b3, x3, Jcode, shearJkey);

	// Get the distributed load interpolation matrix for each section
	this->getDistrLoadInterpMatrix (bp1, x1, Icode);
	this->getDistrLoadInterpMatrix (bp3, x3, Jcode);
}

void
BeamWithHinges2d::getForceInterpMatrix (Matrix &b, double x, const ID &code, int &shearKey)
{			
    b.Zero();

	shearKey = -1;

    double xsi = x/L;

    for (int i = 0; i < code.Size(); i++) {
	switch (code(i)) {
	  case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
	    b(i,1) = xsi - 1.0;
	    b(i,2) = xsi;
	    break;
	  case SECTION_RESPONSE_P:		// Axial, P, interpolation
	    b(i,0) = 1.0;
	    break;
	  case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
	    b(i,1) = 1.0/L;
	    b(i,2) = 1.0/L;
	    shearKey = i;
	    break;
	  default:
	    break;
	}
    }
}

void
BeamWithHinges2d::getDistrLoadInterpMatrix (Matrix &bp, double x, const ID & code)
{
    bp.Zero();

    double xsi = x/L;

    for (int i = 0; i < code.Size(); i++) {
	switch (code(i)) {
	  case SECTION_RESPONSE_MZ:		// Moment, Mz, interpolation
	    bp(i,1) = 0.5*xsi*(xsi-1);
	    break;
	  case SECTION_RESPONSE_P:		// Axial, P, interpolation
	    bp(i,0) = 1 - xsi;
	    break;
	  case SECTION_RESPONSE_VY:		// Shear, Vy, interpolation
	    bp(i,1) = xsi - 0.5;
	    break;
	  default:
	    break;
	}
    }
}

int
BeamWithHinges2d::setResponse (char **argv, int argc, Information &info)
{
    //
    // we compare argv[0] for known response types for BeamWithHinges
    //

    // hinge rotations
    if (strcmp(argv[0],"rotation") == 0) {
	Vector *newVector = new Vector(2);
	if (newVector == 0) {
	    g3ErrorHandler->warning("WARNING BeamWithHinges2d::setResponse() - %d out of memory creating vector\n",
				    this->getTag());
	    return -1;
	}
	info.theVector = newVector;
	info.theType = VectorType;
	return 1;
    }

    // basic forces
    else if (strcmp(argv[0],"force") == 0) {
	Vector *newVector = new Vector(6);
	if (newVector == 0) {
	    g3ErrorHandler->warning("WARNING BeamWithHinges2d::setResponse() - %d out of memory creating vector\n",
				    this->getTag());
	    return -1;
	}
	info.theVector = newVector;
	info.theType = VectorType;
	return 2;
    }
    
    // basic stiffness
    else if (strcmp(argv[0],"stiffness") == 0) {
	Matrix *newMatrix = new Matrix(3,3);
	if (newMatrix == 0) {
	    g3ErrorHandler->warning("WARNING BeamWithHinges2d::setResponse() - %d out of memory creating matrix\n",
				    this->getTag());
	    return -1;
	}
	info.theMatrix = newMatrix;
	info.theType = MatrixType;
	return 3;
    }

    // section response
    else if (strcmp(argv[0],"section") == 0) {
	int sectionNum = atoi(argv[1]);
	
	int ok = 0;
	
	if (sectionNum == 1)
	    ok += sectionI->setResponse(&argv[2], argc-2, info);
	else if (sectionNum == 2)
	    ok += sectionJ->setResponse(&argv[2], argc-2, info);
	else
	    return -1;
	
	if (ok < 0)
	    return -1;
	else if (ok >= 0 && ok < MAX_SECTION_RESPONSE_ID)
	    return sectionNum*MAX_SECTION_RESPONSE_ID + ok;
    }

    else
	return -1;
}

int
BeamWithHinges2d::getResponse (int responseID, Information &info)
{
    switch (responseID) {
      case -1:
	return -1;
      case 1: { // hinge rotations
	  int i;
	  Vector &temp = *(info.theVector);
	  temp.Zero();
	  const Vector &defI = sectionI->getSectionDeformation();
	  int orderI = sectionI->getOrder();
	  const ID &codeI = sectionI->getType();
	  for (i = 0; i < orderI; i++)
		  if (codeI(i) == SECTION_RESPONSE_MZ)
			  temp(0) += defI(i);
	  const Vector &defJ = sectionJ->getSectionDeformation();
	  int orderJ = sectionJ->getOrder();
	  const ID &codeJ = sectionJ->getType();
	  for (i = 0; i < orderJ; i++)
		  if (codeJ(i) == SECTION_RESPONSE_MZ)
			  temp(1) += defJ(i);
	  return 0;
      }
      case 2: { // basic forces
	  if (info.theVector != 0)
	      *(info.theVector) = P;
	  return 0;
      }
      case 3: { // basic stiffness
	  if (info.theMatrix != 0)
	      *(info.theMatrix) = kbCommit;
	  return 0;
      }
      default: {
	  if (responseID >= MAX_SECTION_RESPONSE_ID) { // section quantity
	      
	      int sectionNum = responseID/MAX_SECTION_RESPONSE_ID; 
	      
	      if (sectionNum > 0 && sectionNum <= 2) {
		  if (sectionNum == 1)
		      return sectionI->getResponse (responseID-MAX_SECTION_RESPONSE_ID, info);
		  else if (sectionNum == 2)
		      return sectionJ->getResponse (responseID-MAX_SECTION_RESPONSE_ID*2, info);
	      }
	      else // unknown section
		  return -1;
	  }
	  else // unknown response quantity
	      return -1;
      }
    }
}

int
BeamWithHinges2d::setParameter (char **argv, int argc, Information &info)
{
    // E of the beam interior
    if (strcmp(argv[0],"E") == 0) {
	info.theType = DoubleType;
	return 1;
    }
    
    // G of the beam interior
    else if (strcmp(argv[0],"G") == 0) {
	info.theType = DoubleType;
	return 2;
    }
    
    // A of the beam interior
    else if (strcmp(argv[0],"A") == 0) {
	info.theType = DoubleType;
	return 3;
    }
    
    // I of the beam interior
    else if (strcmp(argv[0],"I") == 0) {
	info.theType = DoubleType;
	return 4;
    }
    
    // alpha of the beam interior
    else if (strcmp(argv[0],"a") == 0 || strcmp(argv[0],"alpha") == 0) {
	info.theType = DoubleType;
	return 5;
    }
    
    // Section parameter
    else if (strcmp(argv[0],"section") ==0) {
	if (argc <= 2)
	    return -1;
	
	int sectionNum = atoi(argv[1]);
	
	int ok = -1;
	
	if (sectionNum == 1)
	    ok = sectionI->setParameter (&argv[2], argc-2, info);
	if (sectionNum == 2)
	    ok = sectionJ->setParameter (&argv[2], argc-2, info);
	
	if (ok < 0)
	    return -1;
	else if (ok < 100)
	    return sectionNum*100 + ok;
	else 
	    return -1;
    }
    
    // Unknown parameter
    else
	return -1;
}	

int
BeamWithHinges2d::updateParameter (int parameterID, Information &info)
{
    switch (parameterID) {
      case 1:
	this->E = info.theDouble;
	return 0;
      case 2:
	this->G = info.theDouble;
	return 0;
      case 3:
	this->A = info.theDouble;
	return 0;
      case 4:
	this->I = info.theDouble;
	return 0;
      case 5:
	this->alpha = info.theDouble;
	return 0;
      default:
	if (parameterID >= 100) { // section quantity
	    int sectionNum = parameterID/100; 
	    if (sectionNum == 1)
		return sectionI->updateParameter (parameterID-100, info);
	    else if (sectionNum == 2)
		return sectionJ->updateParameter (parameterID-2*100, info);
	    else
		return -1;
	}
	else // unknown
	    return -1;
	}	
}	



