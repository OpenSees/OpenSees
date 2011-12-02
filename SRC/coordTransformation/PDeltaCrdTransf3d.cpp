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
                                                                        
// $Revision: 1.11 $
// $Date: 2004-07-21 00:06:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/PDeltaCrdTransf3d.cpp,v $
                                                                        
                                                                        
// File: ~/crdTransf/PDeltaCrdTransf3d.C
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
// 
// Purpose: This file contains the implementation for the 
// PDeltaCrdTransf3d class. PDeltaCrdTransf3d is a linear
// transformation for a planar frame between the global 
// and basic coordinate systems


#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>

#include <PDeltaCrdTransf3d.h>


// constructor:
PDeltaCrdTransf3d::PDeltaCrdTransf3d(int tag, const Vector &vecInLocXZPlane):
  CrdTransf3d(tag, CRDTR_TAG_PDeltaCrdTransf3d),
  nodeIPtr(0), nodeJPtr(0),
  nodeIOffset(0), nodeJOffset(0),
  L(0), ul17(0), ul28(0)
{
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 3; j++)
			R[i][j] = 0.0;

	R[2][0] = vecInLocXZPlane(0);
	R[2][1] = vecInLocXZPlane(1);
	R[2][2] = vecInLocXZPlane(2);

	// Does nothing
}

// constructor:
PDeltaCrdTransf3d::PDeltaCrdTransf3d(int tag, const Vector &vecInLocXZPlane,
										   const Vector &rigJntOffset1,
										   const Vector &rigJntOffset2):
  CrdTransf3d(tag, CRDTR_TAG_PDeltaCrdTransf3d),
  nodeIPtr(0), nodeJPtr(0),
  nodeIOffset(0), nodeJOffset(0),
  L(0), ul17(0), ul28(0)
{
  	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 3; j++)
			R[i][j] = 0.0;

	R[2][0] = vecInLocXZPlane(0);
	R[2][1] = vecInLocXZPlane(1);
	R[2][2] = vecInLocXZPlane(2);

	// check rigid joint offset for node I
	if (&rigJntOffset1 == 0 || rigJntOffset1.Size() != 3 ) {
		opserr << "PDeltaCrdTransf3d::PDeltaCrdTransf3d:  Invalid rigid joint offset vector for node I\n";
		opserr << "Size must be 3\n";      
	}
	else if (rigJntOffset1.Norm() > 0.0) {
		nodeIOffset = new double[3];
		nodeIOffset[0] = rigJntOffset1(0);
		nodeIOffset[1] = rigJntOffset1(1);
		nodeIOffset[2] = rigJntOffset1(2);
	}
   
   // check rigid joint offset for node J
	if (&rigJntOffset2 == 0 || rigJntOffset2.Size() != 3 ) {
		opserr << "PDeltaCrdTransf3d::PDeltaCrdTransf3d:  Invalid rigid joint offset vector for node J\n";
		opserr << "Size must be 3\n";      
	}
	else if (rigJntOffset2.Norm() > 0.0) {
		nodeJOffset = new double[3];
		nodeJOffset[0] = rigJntOffset2(0);
		nodeJOffset[1] = rigJntOffset2(1);
		nodeJOffset[2] = rigJntOffset2(2);
	}
}



 
// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
PDeltaCrdTransf3d::PDeltaCrdTransf3d():
  CrdTransf3d(0, CRDTR_TAG_PDeltaCrdTransf3d),
  nodeIPtr(0), nodeJPtr(0),
  nodeIOffset(0), nodeJOffset(0),
  L(0), ul17(0), ul28(0)
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			R[i][j] = 0.0;
}



// destructor:
PDeltaCrdTransf3d::~PDeltaCrdTransf3d() 
{
	if (nodeIOffset)
		delete [] nodeIOffset;
	if (nodeJOffset)
		delete [] nodeJOffset;
}


int
PDeltaCrdTransf3d::commitState(void)
{
   return 0;
}


int
PDeltaCrdTransf3d::revertToLastCommit(void)
{
   return 0;
}


int
PDeltaCrdTransf3d::revertToStart(void)
{
   return 0;
}


int 
PDeltaCrdTransf3d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
   int error;

   nodeIPtr = nodeIPointer;
   nodeJPtr = nodeJPointer;

   if ((!nodeIPtr) || (!nodeJPtr))
   {
      opserr << "\nPDeltaCrdTransf3d::initialize";
      opserr << "\ninvalid pointers to the element nodes\n";
      return -1;
   }
       
   // get element length and orientation
   if ((error = this->computeElemtLengthAndOrient()))
      return error;

   static Vector XAxis(3);
   static Vector YAxis(3);
   static Vector ZAxis(3);
      
   // get 3by3 rotation matrix
   if ((error = this->getLocalAxes(XAxis, YAxis, ZAxis)))      
      return error;

   return 0;
}


int
PDeltaCrdTransf3d::update(void)
{
	const Vector &disp1 = nodeIPtr->getTrialDisp();
	const Vector &disp2 = nodeJPtr->getTrialDisp();

	static double ug[12];
	for (int i = 0; i < 6; i++) {
		ug[i]   = disp1(i);
		ug[i+6] = disp2(i);
	}

	double ul1, ul7, ul2, ul8;

	ul1 = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
	ul2 = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];

	ul7 = R[1][0]*ug[6] + R[1][1]*ug[7] + R[1][2]*ug[8];
	ul8 = R[2][0]*ug[6] + R[2][1]*ug[7] + R[2][2]*ug[8];

	static double Wu[3];

	if (nodeIOffset) {
		Wu[0] =  nodeIOffset[2]*ug[4] - nodeIOffset[1]*ug[5];
		Wu[1] = -nodeIOffset[2]*ug[3] + nodeIOffset[0]*ug[5];
		Wu[2] =  nodeIOffset[1]*ug[3] - nodeIOffset[0]*ug[4];

		ul1 += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul2 += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}

	if (nodeJOffset) {
		Wu[0] =  nodeJOffset[2]*ug[10] - nodeJOffset[1]*ug[11];
		Wu[1] = -nodeJOffset[2]*ug[9]  + nodeJOffset[0]*ug[11];
		Wu[2] =  nodeJOffset[1]*ug[9]  - nodeJOffset[0]*ug[10];

		ul7 += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul8 += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}

	ul17 = ul1-ul7;
	ul28 = ul2-ul8;

	return 0;
}


int 
PDeltaCrdTransf3d::computeElemtLengthAndOrient()
{
   // element projection
   static Vector dx(3);

   const Vector &ndICoords = nodeIPtr->getCrds();
   const Vector &ndJCoords = nodeJPtr->getCrds();

   dx(0) = ndJCoords(0) - ndICoords(0);
   dx(1) = ndJCoords(1) - ndICoords(1);
   dx(2) = ndJCoords(2) - ndICoords(2);
   
   if (nodeJOffset != 0) {
     dx(0) += nodeJOffset[0];
     dx(1) += nodeJOffset[1];
     dx(2) += nodeJOffset[2];
   }

   if (nodeIOffset != 0) {
     dx(0) -= nodeIOffset[0];
     dx(1) -= nodeIOffset[1];
     dx(2) -= nodeIOffset[2];
   }
   
   // calculate the element length
   L = dx.Norm();

   if (L == 0.0) {
      opserr << "\nPDeltaCrdTransf3d::computeElemtLengthAndOrien: 0 length\n";
      return -2;  
   }

   // calculate the element local x axis components (direction cossines)
   // wrt to the global coordinates
   R[0][0] = dx(0)/L;
   R[0][1] = dx(1)/L;
   R[0][2] = dx(2)/L;
   
   return 0;
}


int
PDeltaCrdTransf3d::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
	// Compute y = v cross x
	// Note: v(i) is stored in R[2][i]
	static Vector vAxis(3);
	vAxis(0) = R[2][0];	vAxis(1) = R[2][1];	vAxis(2) = R[2][2];

	static Vector xAxis(3);
	xAxis(0) = R[0][0];	xAxis(1) = R[0][1];	xAxis(2) = R[0][2];
	XAxis(0) = xAxis(0);    XAxis(1) = xAxis(1);    XAxis(2) = xAxis(2);

	static Vector yAxis(3);

	yAxis(0) = vAxis(1)*xAxis(2) - vAxis(2)*xAxis(1);
	yAxis(1) = vAxis(2)*xAxis(0) - vAxis(0)*xAxis(2);
	yAxis(2) = vAxis(0)*xAxis(1) - vAxis(1)*xAxis(0);

	double ynorm = yAxis.Norm();

	if (ynorm == 0) {
		opserr << "\nPDeltaCrdTransf3d::getLocalAxes";
		opserr << "\nvector v that defines plane xz is parallel to x axis\n";
		return -3;
	}

	yAxis /= ynorm;

	YAxis(0) = yAxis(0);    YAxis(1) = yAxis(1);    YAxis(2) = yAxis(2);

	// Compute z = x cross y
	static Vector zAxis(3);

	zAxis(0) = xAxis(1)*yAxis(2) - xAxis(2)*yAxis(1);
	zAxis(1) = xAxis(2)*yAxis(0) - xAxis(0)*yAxis(2);
	zAxis(2) = xAxis(0)*yAxis(1) - xAxis(1)*yAxis(0);
	ZAxis(0) = zAxis(0);    ZAxis(1) = zAxis(1);    ZAxis(2) = zAxis(2);

	// Fill in transformation matrix
	R[1][0] = yAxis(0);
	R[1][1] = yAxis(1);
	R[1][2] = yAxis(2);

	R[2][0] = zAxis(0);
	R[2][1] = zAxis(1);
	R[2][2] = zAxis(2);

	return 0;
}



double 
PDeltaCrdTransf3d::getInitialLength(void)
{
   return L;
}


double 
PDeltaCrdTransf3d::getDeformedLength(void)
{
   return L;
}


const Vector &
PDeltaCrdTransf3d::getBasicTrialDisp (void)
{
	// determine global displacements
	const Vector &disp1 = nodeIPtr->getTrialDisp();
	const Vector &disp2 = nodeJPtr->getTrialDisp();

	static double ug[12];
	for (int i = 0; i < 6; i++) {
		ug[i]   = disp1(i);
		ug[i+6] = disp2(i);
	}

	double oneOverL = 1.0/L;

	static Vector ub(6);

	static double ul[12];

	ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
	ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
	ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];

	ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
	ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
	ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];

	ul[6]  = R[0][0]*ug[6] + R[0][1]*ug[7] + R[0][2]*ug[8];
	ul[7]  = R[1][0]*ug[6] + R[1][1]*ug[7] + R[1][2]*ug[8];
	ul[8]  = R[2][0]*ug[6] + R[2][1]*ug[7] + R[2][2]*ug[8];

	ul[9]  = R[0][0]*ug[9] + R[0][1]*ug[10] + R[0][2]*ug[11];
	ul[10] = R[1][0]*ug[9] + R[1][1]*ug[10] + R[1][2]*ug[11];
	ul[11] = R[2][0]*ug[9] + R[2][1]*ug[10] + R[2][2]*ug[11];

	static double Wu[3];
	if (nodeIOffset) {
		Wu[0] =  nodeIOffset[2]*ug[4] - nodeIOffset[1]*ug[5];
		Wu[1] = -nodeIOffset[2]*ug[3] + nodeIOffset[0]*ug[5];
		Wu[2] =  nodeIOffset[1]*ug[3] - nodeIOffset[0]*ug[4];

		ul[0] += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
		ul[1] += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul[2] += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}

	if (nodeJOffset) {
		Wu[0] =  nodeJOffset[2]*ug[10] - nodeJOffset[1]*ug[11];
		Wu[1] = -nodeJOffset[2]*ug[9]  + nodeJOffset[0]*ug[11];
		Wu[2] =  nodeJOffset[1]*ug[9]  - nodeJOffset[0]*ug[10];

		ul[6] += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
		ul[7] += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul[8] += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}

	ub(0) = ul[6] - ul[0];
	double tmp;
	tmp = oneOverL*(ul[1]-ul[7]);
	ub(1) = ul[5] + tmp;
	ub(2) = ul[11] + tmp;
	tmp = oneOverL*(ul[8]-ul[2]);
	ub(3) = ul[4] + tmp;
	ub(4) = ul[10] + tmp;
	ub(5) = ul[9] - ul[3];

	return ub;
}


const Vector &
PDeltaCrdTransf3d::getBasicIncrDisp (void)
{
	// determine global displacements
	const Vector &disp1 = nodeIPtr->getIncrDisp();
	const Vector &disp2 = nodeJPtr->getIncrDisp();

	static double ug[12];
	for (int i = 0; i < 6; i++) {
		ug[i]   = disp1(i);
		ug[i+6] = disp2(i);
	}

	double oneOverL = 1.0/L;

	static Vector ub(6);

	static double ul[12];

	ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
	ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
	ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];

	ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
	ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
	ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];

	ul[6]  = R[0][0]*ug[6] + R[0][1]*ug[7] + R[0][2]*ug[8];
	ul[7]  = R[1][0]*ug[6] + R[1][1]*ug[7] + R[1][2]*ug[8];
	ul[8]  = R[2][0]*ug[6] + R[2][1]*ug[7] + R[2][2]*ug[8];

	ul[9]  = R[0][0]*ug[9] + R[0][1]*ug[10] + R[0][2]*ug[11];
	ul[10] = R[1][0]*ug[9] + R[1][1]*ug[10] + R[1][2]*ug[11];
	ul[11] = R[2][0]*ug[9] + R[2][1]*ug[10] + R[2][2]*ug[11];

	static double Wu[3];
	if (nodeIOffset) {
		Wu[0] =  nodeIOffset[2]*ug[4] - nodeIOffset[1]*ug[5];
		Wu[1] = -nodeIOffset[2]*ug[3] + nodeIOffset[0]*ug[5];
		Wu[2] =  nodeIOffset[1]*ug[3] - nodeIOffset[0]*ug[4];

		ul[0] += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
		ul[1] += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul[2] += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}

	if (nodeJOffset) {
		Wu[0] =  nodeJOffset[2]*ug[10] - nodeJOffset[1]*ug[11];
		Wu[1] = -nodeJOffset[2]*ug[9]  + nodeJOffset[0]*ug[11];
		Wu[2] =  nodeJOffset[1]*ug[9]  - nodeJOffset[0]*ug[10];

		ul[6] += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
		ul[7] += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul[8] += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}

	ub(0) = ul[6] - ul[0];
	double tmp;
	tmp = oneOverL*(ul[1]-ul[7]);
	ub(1) = ul[5] + tmp;
	ub(2) = ul[11] + tmp;
	tmp = oneOverL*(ul[8]-ul[2]);
	ub(3) = ul[4] + tmp;
	ub(4) = ul[10] + tmp;
	ub(5) = ul[9] - ul[3];

	return ub;
}


const Vector &
PDeltaCrdTransf3d::getBasicIncrDeltaDisp(void)
{
	// determine global displacements
	const Vector &disp1 = nodeIPtr->getIncrDeltaDisp();
	const Vector &disp2 = nodeJPtr->getIncrDeltaDisp();

	static double ug[12];
	for (int i = 0; i < 6; i++) {
		ug[i]   = disp1(i);
		ug[i+6] = disp2(i);
	}

	double oneOverL = 1.0/L;

	static Vector ub(6);

	static double ul[12];

	ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
	ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
	ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];

	ul[3]  = R[0][0]*ug[3] + R[0][1]*ug[4] + R[0][2]*ug[5];
	ul[4]  = R[1][0]*ug[3] + R[1][1]*ug[4] + R[1][2]*ug[5];
	ul[5]  = R[2][0]*ug[3] + R[2][1]*ug[4] + R[2][2]*ug[5];

	ul[6]  = R[0][0]*ug[6] + R[0][1]*ug[7] + R[0][2]*ug[8];
	ul[7]  = R[1][0]*ug[6] + R[1][1]*ug[7] + R[1][2]*ug[8];
	ul[8]  = R[2][0]*ug[6] + R[2][1]*ug[7] + R[2][2]*ug[8];

	ul[9]  = R[0][0]*ug[9] + R[0][1]*ug[10] + R[0][2]*ug[11];
	ul[10] = R[1][0]*ug[9] + R[1][1]*ug[10] + R[1][2]*ug[11];
	ul[11] = R[2][0]*ug[9] + R[2][1]*ug[10] + R[2][2]*ug[11];

	static double Wu[3];
	if (nodeIOffset) {
		Wu[0] =  nodeIOffset[2]*ug[4] - nodeIOffset[1]*ug[5];
		Wu[1] = -nodeIOffset[2]*ug[3] + nodeIOffset[0]*ug[5];
		Wu[2] =  nodeIOffset[1]*ug[3] - nodeIOffset[0]*ug[4];

		ul[0] += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
		ul[1] += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul[2] += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}

	if (nodeJOffset) {
		Wu[0] =  nodeJOffset[2]*ug[10] - nodeJOffset[1]*ug[11];
		Wu[1] = -nodeJOffset[2]*ug[9]  + nodeJOffset[0]*ug[11];
		Wu[2] =  nodeJOffset[1]*ug[9]  - nodeJOffset[0]*ug[10];

		ul[6] += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
		ul[7] += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul[8] += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}

	ub(0) = ul[6] - ul[0];
	double tmp;
	tmp = oneOverL*(ul[1]-ul[7]);
	ub(1) = ul[5] + tmp;
	ub(2) = ul[11] + tmp;
	tmp = oneOverL*(ul[8]-ul[2]);
	ub(3) = ul[4] + tmp;
	ub(4) = ul[10] + tmp;
	ub(5) = ul[9] - ul[3];

	return ub;
}


const Vector &
PDeltaCrdTransf3d::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
	// transform resisting forces from the basic system to local coordinates
	static double pl[12];

	double q0 = pb(0);
	double q1 = pb(1);
	double q2 = pb(2);
	double q3 = pb(3);
	double q4 = pb(4);
	double q5 = pb(5);

	double oneOverL = 1.0/L;

	pl[0]  = -q0;
	pl[1]  =  oneOverL*(q1+q2);
	pl[2]  = -oneOverL*(q3+q4);
	pl[3]  = -q5;
	pl[4]  =  q3;
	pl[5]  =  q1;
	pl[6]  =  q0;
	pl[7]  = -pl[1];
	pl[8]  = -pl[2];
	pl[9]  =  q5;
	pl[10] =  q4;
	pl[11] =  q2;

	pl[0] += p0(0);
	pl[1] += p0(1);
	pl[7] += p0(2);
	pl[2] += p0(3);
	pl[8] += p0(4);

	// Include leaning column effects (P-Delta)
	double NoverL;
	NoverL = ul17*q0*oneOverL;             
	pl[1] += NoverL;
	pl[7] -= NoverL;
	NoverL = ul28*q0*oneOverL;
	pl[2] += NoverL;
	pl[8] -= NoverL;

	// transform resisting forces  from local to global coordinates
	static Vector pg(12);

	pg(0)  = R[0][0]*pl[0] + R[1][0]*pl[1] + R[2][0]*pl[2];
	pg(1)  = R[0][1]*pl[0] + R[1][1]*pl[1] + R[2][1]*pl[2];
	pg(2)  = R[0][2]*pl[0] + R[1][2]*pl[1] + R[2][2]*pl[2];

	pg(3)  = R[0][0]*pl[3] + R[1][0]*pl[4] + R[2][0]*pl[5];
	pg(4)  = R[0][1]*pl[3] + R[1][1]*pl[4] + R[2][1]*pl[5];
	pg(5)  = R[0][2]*pl[3] + R[1][2]*pl[4] + R[2][2]*pl[5];

	pg(6)  = R[0][0]*pl[6] + R[1][0]*pl[7] + R[2][0]*pl[8];
	pg(7)  = R[0][1]*pl[6] + R[1][1]*pl[7] + R[2][1]*pl[8];
	pg(8)  = R[0][2]*pl[6] + R[1][2]*pl[7] + R[2][2]*pl[8];

	pg(9)  = R[0][0]*pl[9] + R[1][0]*pl[10] + R[2][0]*pl[11];
	pg(10) = R[0][1]*pl[9] + R[1][1]*pl[10] + R[2][1]*pl[11];
	pg(11) = R[0][2]*pl[9] + R[1][2]*pl[10] + R[2][2]*pl[11];

	if (nodeIOffset) {
		pg(3) += -nodeIOffset[2]*pg(1) + nodeIOffset[1]*pg(2);
		pg(4) +=  nodeIOffset[2]*pg(0) - nodeIOffset[0]*pg(2);
		pg(5) += -nodeIOffset[1]*pg(0) + nodeIOffset[0]*pg(1);
	}

	if (nodeJOffset) {
		pg(9)  += -nodeJOffset[2]*pg(7) + nodeJOffset[1]*pg(8);
		pg(10) +=  nodeJOffset[2]*pg(6) - nodeJOffset[0]*pg(8);
		pg(11) += -nodeJOffset[1]*pg(6) + nodeJOffset[0]*pg(7);
	}
	
	return pg;
}


const Matrix &
PDeltaCrdTransf3d::getGlobalStiffMatrix (const Matrix &KB, const Vector &pb)
{
	static Matrix kg(12,12);	// Global stiffness for return
	static double kb[6][6];		// Basic stiffness
	static double kl[12][12];	// Local stiffness
	static double tmp[12][12];	// Temporary storage

	double oneOverL = 1.0/L;

	int i,j;
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
			kb[i][j] = KB(i,j);

	// Transform basic stiffness to local system
	// First compute kb*T_{bl}
	for (i = 0; i < 6; i++) {
		tmp[i][0]  = -kb[i][0];
		tmp[i][1]  =  oneOverL*(kb[i][1]+kb[i][2]);
		tmp[i][2]  = -oneOverL*(kb[i][3]+kb[i][4]);
		tmp[i][3]  = -kb[i][5];
		tmp[i][4]  =  kb[i][3];
		tmp[i][5]  =  kb[i][1];
		tmp[i][6]  =  kb[i][0];
		tmp[i][7]  = -tmp[i][1];
		tmp[i][8]  = -tmp[i][2];
		tmp[i][9]  =  kb[i][5];
		tmp[i][10] =  kb[i][4];
		tmp[i][11] =  kb[i][2];
	}
	
	// Now compute T'_{bl}*(kb*T_{bl})
	for (i = 0; i < 12; i++) {
		kl[0][i]  = -tmp[0][i];
		kl[1][i]  =  oneOverL*(tmp[1][i]+tmp[2][i]);
		kl[2][i]  = -oneOverL*(tmp[3][i]+tmp[4][i]);
		kl[3][i]  = -tmp[5][i];
		kl[4][i]  =  tmp[3][i];
		kl[5][i]  =  tmp[1][i];
		kl[6][i]  =  tmp[0][i];
		kl[7][i]  = -kl[1][i];
		kl[8][i]  = -kl[2][i];
		kl[9][i]  =  tmp[5][i];
		kl[10][i] =  tmp[4][i];
		kl[11][i] =  tmp[2][i];
	}

	// Include geometric stiffness effects in local system
	double NoverL = pb(0)*oneOverL;
	kl[1][1] += NoverL;
	kl[2][2] += NoverL;
	kl[7][7] += NoverL;
	kl[8][8] += NoverL;
	kl[1][7] -= NoverL;
	kl[7][1] -= NoverL;
	kl[2][8] -= NoverL;
	kl[8][2] -= NoverL;
	
	static double RWI[3][3];

	if (nodeIOffset) {
		// Compute RWI
		RWI[0][0] = -R[0][1]*nodeIOffset[2] + R[0][2]*nodeIOffset[1];
		RWI[1][0] = -R[1][1]*nodeIOffset[2] + R[1][2]*nodeIOffset[1];
		RWI[2][0] = -R[2][1]*nodeIOffset[2] + R[2][2]*nodeIOffset[1];

		RWI[0][1] =  R[0][0]*nodeIOffset[2] - R[0][2]*nodeIOffset[0];
		RWI[1][1] =  R[1][0]*nodeIOffset[2] - R[1][2]*nodeIOffset[0];
		RWI[2][1] =  R[2][0]*nodeIOffset[2] - R[2][2]*nodeIOffset[0];

		RWI[0][2] = -R[0][0]*nodeIOffset[1] + R[0][1]*nodeIOffset[0];
		RWI[1][2] = -R[1][0]*nodeIOffset[1] + R[1][1]*nodeIOffset[0];
		RWI[2][2] = -R[2][0]*nodeIOffset[1] + R[2][1]*nodeIOffset[0];
	}

	static double RWJ[3][3];

	if (nodeJOffset) {
		// Compute RWJ
		RWJ[0][0] = -R[0][1]*nodeJOffset[2] + R[0][2]*nodeJOffset[1];
		RWJ[1][0] = -R[1][1]*nodeJOffset[2] + R[1][2]*nodeJOffset[1];
		RWJ[2][0] = -R[2][1]*nodeJOffset[2] + R[2][2]*nodeJOffset[1];

		RWJ[0][1] =  R[0][0]*nodeJOffset[2] - R[0][2]*nodeJOffset[0];
		RWJ[1][1] =  R[1][0]*nodeJOffset[2] - R[1][2]*nodeJOffset[0];
		RWJ[2][1] =  R[2][0]*nodeJOffset[2] - R[2][2]*nodeJOffset[0];

		RWJ[0][2] = -R[0][0]*nodeJOffset[1] + R[0][1]*nodeJOffset[0];
		RWJ[1][2] = -R[1][0]*nodeJOffset[1] + R[1][1]*nodeJOffset[0];
		RWJ[2][2] = -R[2][0]*nodeJOffset[1] + R[2][1]*nodeJOffset[0];
	}

	// Transform local stiffness to global system
	// First compute kl*T_{lg}
	int m;
	for (m = 0; m < 12; m++) {
		tmp[m][0] = kl[m][0]*R[0][0] + kl[m][1]*R[1][0]  + kl[m][2]*R[2][0];
		tmp[m][1] = kl[m][0]*R[0][1] + kl[m][1]*R[1][1]  + kl[m][2]*R[2][1];
		tmp[m][2] = kl[m][0]*R[0][2] + kl[m][1]*R[1][2]  + kl[m][2]*R[2][2];

		tmp[m][3] = kl[m][3]*R[0][0] + kl[m][4]*R[1][0]  + kl[m][5]*R[2][0];
		tmp[m][4] = kl[m][3]*R[0][1] + kl[m][4]*R[1][1]  + kl[m][5]*R[2][1];
		tmp[m][5] = kl[m][3]*R[0][2] + kl[m][4]*R[1][2]  + kl[m][5]*R[2][2];

		if (nodeIOffset) {
			tmp[m][3]  += kl[m][0]*RWI[0][0]  + kl[m][1]*RWI[1][0]  + kl[m][2]*RWI[2][0];
			tmp[m][4]  += kl[m][0]*RWI[0][1]  + kl[m][1]*RWI[1][1]  + kl[m][2]*RWI[2][1];
			tmp[m][5]  += kl[m][0]*RWI[0][2]  + kl[m][1]*RWI[1][2]  + kl[m][2]*RWI[2][2];
		}

		tmp[m][6] = kl[m][6]*R[0][0] + kl[m][7]*R[1][0]  + kl[m][8]*R[2][0];
		tmp[m][7] = kl[m][6]*R[0][1] + kl[m][7]*R[1][1]  + kl[m][8]*R[2][1];
		tmp[m][8] = kl[m][6]*R[0][2] + kl[m][7]*R[1][2]  + kl[m][8]*R[2][2];

		tmp[m][9]  = kl[m][9]*R[0][0] + kl[m][10]*R[1][0] + kl[m][11]*R[2][0];
		tmp[m][10] = kl[m][9]*R[0][1] + kl[m][10]*R[1][1] + kl[m][11]*R[2][1];
		tmp[m][11] = kl[m][9]*R[0][2] + kl[m][10]*R[1][2] + kl[m][11]*R[2][2];

		if (nodeJOffset) {
			tmp[m][9]   += kl[m][6]*RWJ[0][0]  + kl[m][7]*RWJ[1][0]  + kl[m][8]*RWJ[2][0];
			tmp[m][10]  += kl[m][6]*RWJ[0][1]  + kl[m][7]*RWJ[1][1]  + kl[m][8]*RWJ[2][1];
			tmp[m][11]  += kl[m][6]*RWJ[0][2]  + kl[m][7]*RWJ[1][2]  + kl[m][8]*RWJ[2][2];
		}

	}

	// Now compute T'_{lg}*(kl*T_{lg})
	for (m = 0; m < 12; m++) {
		kg(0,m) = R[0][0]*tmp[0][m] + R[1][0]*tmp[1][m]  + R[2][0]*tmp[2][m];
		kg(1,m) = R[0][1]*tmp[0][m] + R[1][1]*tmp[1][m]  + R[2][1]*tmp[2][m];
		kg(2,m) = R[0][2]*tmp[0][m] + R[1][2]*tmp[1][m]  + R[2][2]*tmp[2][m];

		kg(3,m) = R[0][0]*tmp[3][m] + R[1][0]*tmp[4][m]  + R[2][0]*tmp[5][m];
		kg(4,m) = R[0][1]*tmp[3][m] + R[1][1]*tmp[4][m]  + R[2][1]*tmp[5][m];
		kg(5,m) = R[0][2]*tmp[3][m] + R[1][2]*tmp[4][m]  + R[2][2]*tmp[5][m];

		if (nodeIOffset) {
			kg(3,m) += RWI[0][0]*tmp[0][m]  + RWI[1][0]*tmp[1][m] + RWI[2][0]*tmp[2][m];
			kg(4,m) += RWI[0][1]*tmp[0][m]  + RWI[1][1]*tmp[1][m] + RWI[2][1]*tmp[2][m];
			kg(5,m) += RWI[0][2]*tmp[0][m]  + RWI[1][2]*tmp[1][m] + RWI[2][2]*tmp[2][m];
		}

		kg(6,m) = R[0][0]*tmp[6][m] + R[1][0]*tmp[7][m]  + R[2][0]*tmp[8][m];
		kg(7,m) = R[0][1]*tmp[6][m] + R[1][1]*tmp[7][m]  + R[2][1]*tmp[8][m];
		kg(8,m) = R[0][2]*tmp[6][m] + R[1][2]*tmp[7][m]  + R[2][2]*tmp[8][m];

		kg(9,m)  = R[0][0]*tmp[9][m] + R[1][0]*tmp[10][m] + R[2][0]*tmp[11][m];
		kg(10,m) = R[0][1]*tmp[9][m] + R[1][1]*tmp[10][m] + R[2][1]*tmp[11][m];
		kg(11,m) = R[0][2]*tmp[9][m] + R[1][2]*tmp[10][m] + R[2][2]*tmp[11][m];

		if (nodeJOffset) {
			kg(9,m)  += RWJ[0][0]*tmp[6][m]  + RWJ[1][0]*tmp[7][m] + RWJ[2][0]*tmp[8][m];
			kg(10,m) += RWJ[0][1]*tmp[6][m]  + RWJ[1][1]*tmp[7][m] + RWJ[2][1]*tmp[8][m];
			kg(11,m) += RWJ[0][2]*tmp[6][m]  + RWJ[1][2]*tmp[7][m] + RWJ[2][2]*tmp[8][m];
		}
	}

	return kg;
}


const Matrix &
PDeltaCrdTransf3d::getInitialGlobalStiffMatrix (const Matrix &KB)
{
	static Matrix kg(12,12);	// Global stiffness for return
	static double kb[6][6];		// Basic stiffness
	static double kl[12][12];	// Local stiffness
	static double tmp[12][12];	// Temporary storage

	double oneOverL = 1.0/L;

	int i,j;
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
			kb[i][j] = KB(i,j);

	// Transform basic stiffness to local system
	// First compute kb*T_{bl}
	for (i = 0; i < 6; i++) {
		tmp[i][0]  = -kb[i][0];
		tmp[i][1]  =  oneOverL*(kb[i][1]+kb[i][2]);
		tmp[i][2]  = -oneOverL*(kb[i][3]+kb[i][4]);
		tmp[i][3]  = -kb[i][5];
		tmp[i][4]  =  kb[i][3];
		tmp[i][5]  =  kb[i][1];
		tmp[i][6]  =  kb[i][0];
		tmp[i][7]  = -tmp[i][1];
		tmp[i][8]  = -tmp[i][2];
		tmp[i][9]  =  kb[i][5];
		tmp[i][10] =  kb[i][4];
		tmp[i][11] =  kb[i][2];
	}
	
	// Now compute T'_{bl}*(kb*T_{bl})
	for (i = 0; i < 12; i++) {
		kl[0][i]  = -tmp[0][i];
		kl[1][i]  =  oneOverL*(tmp[1][i]+tmp[2][i]);
		kl[2][i]  = -oneOverL*(tmp[3][i]+tmp[4][i]);
		kl[3][i]  = -tmp[5][i];
		kl[4][i]  =  tmp[3][i];
		kl[5][i]  =  tmp[1][i];
		kl[6][i]  =  tmp[0][i];
		kl[7][i]  = -kl[1][i];
		kl[8][i]  = -kl[2][i];
		kl[9][i]  =  tmp[5][i];
		kl[10][i] =  tmp[4][i];
		kl[11][i] =  tmp[2][i];
	}

	// Include geometric stiffness effects in local system
	//	double NoverL = pb(0)*oneOverL;
	//kl[1][1] += NoverL;
	//kl[2][2] += NoverL;
	//kl[7][7] += NoverL;
	//kl[8][8] += NoverL;
	//kl[1][7] -= NoverL;
	//kl[7][1] -= NoverL;
	//kl[2][8] -= NoverL;
	//kl[8][2] -= NoverL;

	
	static double RWI[3][3];

	if (nodeIOffset) {
		// Compute RWI
		RWI[0][0] = -R[0][1]*nodeIOffset[2] + R[0][2]*nodeIOffset[1];
		RWI[1][0] = -R[1][1]*nodeIOffset[2] + R[1][2]*nodeIOffset[1];
		RWI[2][0] = -R[2][1]*nodeIOffset[2] + R[2][2]*nodeIOffset[1];

		RWI[0][1] =  R[0][0]*nodeIOffset[2] - R[0][2]*nodeIOffset[0];
		RWI[1][1] =  R[1][0]*nodeIOffset[2] - R[1][2]*nodeIOffset[0];
		RWI[2][1] =  R[2][0]*nodeIOffset[2] - R[2][2]*nodeIOffset[0];

		RWI[0][2] = -R[0][0]*nodeIOffset[1] + R[0][1]*nodeIOffset[0];
		RWI[1][2] = -R[1][0]*nodeIOffset[1] + R[1][1]*nodeIOffset[0];
		RWI[2][2] = -R[2][0]*nodeIOffset[1] + R[2][1]*nodeIOffset[0];
	}

	static double RWJ[3][3];

	if (nodeJOffset) {
		// Compute RWJ
		RWJ[0][0] = -R[0][1]*nodeJOffset[2] + R[0][2]*nodeJOffset[1];
		RWJ[1][0] = -R[1][1]*nodeJOffset[2] + R[1][2]*nodeJOffset[1];
		RWJ[2][0] = -R[2][1]*nodeJOffset[2] + R[2][2]*nodeJOffset[1];

		RWJ[0][1] =  R[0][0]*nodeJOffset[2] - R[0][2]*nodeJOffset[0];
		RWJ[1][1] =  R[1][0]*nodeJOffset[2] - R[1][2]*nodeJOffset[0];
		RWJ[2][1] =  R[2][0]*nodeJOffset[2] - R[2][2]*nodeJOffset[0];

		RWJ[0][2] = -R[0][0]*nodeJOffset[1] + R[0][1]*nodeJOffset[0];
		RWJ[1][2] = -R[1][0]*nodeJOffset[1] + R[1][1]*nodeJOffset[0];
		RWJ[2][2] = -R[2][0]*nodeJOffset[1] + R[2][1]*nodeJOffset[0];
	}

	// Transform local stiffness to global system
	// First compute kl*T_{lg}

	int m;
	for (m = 0; m < 12; m++) {
		tmp[m][0] = kl[m][0]*R[0][0] + kl[m][1]*R[1][0]  + kl[m][2]*R[2][0];
		tmp[m][1] = kl[m][0]*R[0][1] + kl[m][1]*R[1][1]  + kl[m][2]*R[2][1];
		tmp[m][2] = kl[m][0]*R[0][2] + kl[m][1]*R[1][2]  + kl[m][2]*R[2][2];

		tmp[m][3] = kl[m][3]*R[0][0] + kl[m][4]*R[1][0]  + kl[m][5]*R[2][0];
		tmp[m][4] = kl[m][3]*R[0][1] + kl[m][4]*R[1][1]  + kl[m][5]*R[2][1];
		tmp[m][5] = kl[m][3]*R[0][2] + kl[m][4]*R[1][2]  + kl[m][5]*R[2][2];

		if (nodeIOffset) {
			tmp[m][3]  += kl[m][0]*RWI[0][0]  + kl[m][1]*RWI[1][0]  + kl[m][2]*RWI[2][0];
			tmp[m][4]  += kl[m][0]*RWI[0][1]  + kl[m][1]*RWI[1][1]  + kl[m][2]*RWI[2][1];
			tmp[m][5]  += kl[m][0]*RWI[0][2]  + kl[m][1]*RWI[1][2]  + kl[m][2]*RWI[2][2];
		}

		tmp[m][6] = kl[m][6]*R[0][0] + kl[m][7]*R[1][0]  + kl[m][8]*R[2][0];
		tmp[m][7] = kl[m][6]*R[0][1] + kl[m][7]*R[1][1]  + kl[m][8]*R[2][1];
		tmp[m][8] = kl[m][6]*R[0][2] + kl[m][7]*R[1][2]  + kl[m][8]*R[2][2];

		tmp[m][9]  = kl[m][9]*R[0][0] + kl[m][10]*R[1][0] + kl[m][11]*R[2][0];
		tmp[m][10] = kl[m][9]*R[0][1] + kl[m][10]*R[1][1] + kl[m][11]*R[2][1];
		tmp[m][11] = kl[m][9]*R[0][2] + kl[m][10]*R[1][2] + kl[m][11]*R[2][2];

		if (nodeJOffset) {
			tmp[m][9]   += kl[m][6]*RWJ[0][0]  + kl[m][7]*RWJ[1][0]  + kl[m][8]*RWJ[2][0];
			tmp[m][10]  += kl[m][6]*RWJ[0][1]  + kl[m][7]*RWJ[1][1]  + kl[m][8]*RWJ[2][1];
			tmp[m][11]  += kl[m][6]*RWJ[0][2]  + kl[m][7]*RWJ[1][2]  + kl[m][8]*RWJ[2][2];
		}

	}

	// Now compute T'_{lg}*(kl*T_{lg})
	for (m = 0; m < 12; m++) {
		kg(0,m) = R[0][0]*tmp[0][m] + R[1][0]*tmp[1][m]  + R[2][0]*tmp[2][m];
		kg(1,m) = R[0][1]*tmp[0][m] + R[1][1]*tmp[1][m]  + R[2][1]*tmp[2][m];
		kg(2,m) = R[0][2]*tmp[0][m] + R[1][2]*tmp[1][m]  + R[2][2]*tmp[2][m];

		kg(3,m) = R[0][0]*tmp[3][m] + R[1][0]*tmp[4][m]  + R[2][0]*tmp[5][m];
		kg(4,m) = R[0][1]*tmp[3][m] + R[1][1]*tmp[4][m]  + R[2][1]*tmp[5][m];
		kg(5,m) = R[0][2]*tmp[3][m] + R[1][2]*tmp[4][m]  + R[2][2]*tmp[5][m];

		if (nodeIOffset) {
			kg(3,m) += RWI[0][0]*tmp[0][m]  + RWI[1][0]*tmp[1][m] + RWI[2][0]*tmp[2][m];
			kg(4,m) += RWI[0][1]*tmp[0][m]  + RWI[1][1]*tmp[1][m] + RWI[2][1]*tmp[2][m];
			kg(5,m) += RWI[0][2]*tmp[0][m]  + RWI[1][2]*tmp[1][m] + RWI[2][2]*tmp[2][m];
		}

		kg(6,m) = R[0][0]*tmp[6][m] + R[1][0]*tmp[7][m]  + R[2][0]*tmp[8][m];
		kg(7,m) = R[0][1]*tmp[6][m] + R[1][1]*tmp[7][m]  + R[2][1]*tmp[8][m];
		kg(8,m) = R[0][2]*tmp[6][m] + R[1][2]*tmp[7][m]  + R[2][2]*tmp[8][m];

		kg(9,m)  = R[0][0]*tmp[9][m] + R[1][0]*tmp[10][m] + R[2][0]*tmp[11][m];
		kg(10,m) = R[0][1]*tmp[9][m] + R[1][1]*tmp[10][m] + R[2][1]*tmp[11][m];
		kg(11,m) = R[0][2]*tmp[9][m] + R[1][2]*tmp[10][m] + R[2][2]*tmp[11][m];

		if (nodeJOffset) {
			kg(9,m)  += RWJ[0][0]*tmp[6][m]  + RWJ[1][0]*tmp[7][m] + RWJ[2][0]*tmp[8][m];
			kg(10,m) += RWJ[0][1]*tmp[6][m]  + RWJ[1][1]*tmp[7][m] + RWJ[2][1]*tmp[8][m];
			kg(11,m) += RWJ[0][2]*tmp[6][m]  + RWJ[1][2]*tmp[7][m] + RWJ[2][2]*tmp[8][m];
		}
	}

	return kg;
}
  

CrdTransf3d *
PDeltaCrdTransf3d::getCopy(void)
{
  // create a new instance of PDeltaCrdTransf3d 

  PDeltaCrdTransf3d *theCopy;

  static Vector xz(3);
  xz(0) = R[2][0];
  xz(1) = R[2][1];
  xz(2) = R[2][2];

  Vector offsetI(3);
  Vector offsetJ(3);

  if (nodeIOffset) {
	  offsetI(0) = nodeIOffset[0];
	  offsetI(1) = nodeIOffset[1];
	  offsetI(2) = nodeIOffset[2];
  }

  if (nodeJOffset) {
	  offsetJ(0) = nodeJOffset[0];
	  offsetJ(1) = nodeJOffset[1];
	  offsetJ(2) = nodeJOffset[2];
  }
  
  theCopy = new PDeltaCrdTransf3d(this->getTag(), xz, offsetI, offsetJ);

  theCopy->nodeIPtr = nodeIPtr;
  theCopy->nodeJPtr = nodeJPtr;
  theCopy->L = L;
  theCopy->ul17 = ul17;
  theCopy->ul28 = ul28;
  for (int i = 0; i < 3; i++)
	  for (int j = 0; j < 3; j++)
		  theCopy->R[i][j] = R[i][j];

  
  return theCopy;
}


int 
PDeltaCrdTransf3d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(13);
  data(0) = this->getTag();
  data(1) = L;
  if (nodeIOffset != 0) {
    data(2) = 1.0;
    data(3) = nodeIOffset[0];
    data(4) = nodeIOffset[1];
    data(5) = nodeIOffset[2];
  } else
    data(2) = 0.0;
  
  if (nodeJOffset != 0) {
    data(6) = 1.0;
    data(7) = nodeJOffset[0];
    data(8) = nodeJOffset[1];
    data(9) = nodeJOffset[2];
  } else
    data(6) = 0.0;

  data(10) = R[2][0];
  data(11) = R[2][1];
  data(12) = R[2][2];
  
  res += theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "PDeltaCrdTransf3d::recvSelf - failed to send Vector\n";
    return res;
  }
  
  return res;
}


int 
PDeltaCrdTransf3d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(13);
  
  res += theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "PDeltaCrdTransf3d::recvSelf - failed to receive Vector\n";
    return res;
  }
  
  this->setTag((int)data(0));
  L = data(1);
  data(0) = this->getTag();
  data(1) = L;
  if (data(2) == 1.0) {
    if (nodeIOffset == 0)
      nodeIOffset = new double[3];
    nodeIOffset[0] = data(3);
    nodeIOffset[1] = data(4);
    nodeIOffset[2] = data(5);
  } 
  
  if (data(6) == 1.0) {
    if (nodeJOffset == 0)
      nodeJOffset = new double[3];
    nodeJOffset[0] = data(7);
    nodeJOffset[1] = data(8);
    nodeJOffset[2] = data(9);
  } 

  R[2][0] = data(10);
  R[2][1] = data(11);
  R[2][2] = data(12);

  return res;
}
 	
const Vector &
PDeltaCrdTransf3d::getPointGlobalCoordFromLocal(const Vector &xl)
{
   static Vector xg(3);

   //xg = nodeIPtr->getCrds() + nodeIOffset;
   xg = nodeIPtr->getCrds();
   
   if (nodeIOffset) {
	   xg(0) += nodeIOffset[0];
	   xg(1) += nodeIOffset[1];
	   xg(2) += nodeIOffset[2];
   }

   // xg = xg + Rlj'*xl
   //xg.addMatrixTransposeVector(1.0, Rlj, xl, 1.0);
   xg(0) += R[0][0]*xl(0) + R[1][0]*xl(1) + R[2][0]*xl(2);
   xg(1) += R[0][1]*xl(0) + R[1][1]*xl(1) + R[2][1]*xl(2);
   xg(2) += R[0][2]*xl(0) + R[1][2]*xl(1) + R[2][2]*xl(2);
     
   return xg;  
}

   
const Vector &
PDeltaCrdTransf3d::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
   // determine global displacements
   const Vector &disp1 = nodeIPtr->getTrialDisp();
   const Vector &disp2 = nodeJPtr->getTrialDisp();

   static double ug[12];
   for (int i = 0; i < 6; i++)
   {
      ug[i]   = disp1(i);
      ug[i+6] = disp2(i);
   }

   // transform global end displacements to local coordinates
   //ul.addMatrixVector(0.0, Tlg,  ug, 1.0);       //  ul = Tlg *  ug;
	static double ul[12];

	ul[0]  = R[0][0]*ug[0] + R[0][1]*ug[1] + R[0][2]*ug[2];
	ul[1]  = R[1][0]*ug[0] + R[1][1]*ug[1] + R[1][2]*ug[2];
	ul[2]  = R[2][0]*ug[0] + R[2][1]*ug[1] + R[2][2]*ug[2];

	ul[7]  = R[1][0]*ug[6] + R[1][1]*ug[7] + R[1][2]*ug[8];
	ul[8]  = R[2][0]*ug[6] + R[2][1]*ug[7] + R[2][2]*ug[8];

	static double Wu[3];
	if (nodeIOffset) {
		Wu[0] =  nodeIOffset[2]*ug[4] - nodeIOffset[1]*ug[5];
		Wu[1] = -nodeIOffset[2]*ug[3] + nodeIOffset[0]*ug[5];
		Wu[2] =  nodeIOffset[1]*ug[3] - nodeIOffset[0]*ug[4];

		ul[0] += R[0][0]*Wu[0] + R[0][1]*Wu[1] + R[0][2]*Wu[2];
		ul[1] += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul[2] += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}

	if (nodeJOffset) {
		Wu[0] =  nodeJOffset[2]*ug[10] - nodeJOffset[1]*ug[11];
		Wu[1] = -nodeJOffset[2]*ug[9]  + nodeJOffset[0]*ug[11];
		Wu[2] =  nodeJOffset[1]*ug[9]  - nodeJOffset[0]*ug[10];

		ul[7] += R[1][0]*Wu[0] + R[1][1]*Wu[1] + R[1][2]*Wu[2];
		ul[8] += R[2][0]*Wu[0] + R[2][1]*Wu[1] + R[2][2]*Wu[2];
	}
	
   // compute displacements at point xi, in local coordinates
   static double uxl[3];
   static Vector uxg(3);

   uxl[0] = uxb(0) +        ul[0];
   uxl[1] = uxb(1) + (1-xi)*ul[1] + xi*ul[7];
   uxl[2] = uxb(2) + (1-xi)*ul[2] + xi*ul[8];
  
   // rotate displacements to global coordinates
   // uxg = Rlj'*uxl
   //uxg.addMatrixTransposeVector(0.0, Rlj, uxl, 1.0);
   uxg(0) = R[0][0]*uxl[0] + R[1][0]*uxl[1] + R[2][0]*uxl[2];
   uxg(1) = R[0][1]*uxl[0] + R[1][1]*uxl[1] + R[2][1]*uxl[2];
   uxg(2) = R[0][2]*uxl[0] + R[1][2]*uxl[1] + R[2][2]*uxl[2];
     
   return uxg;  
}



void
PDeltaCrdTransf3d::Print(OPS_Stream &s, int flag)
{
   s << "\nCrdTransf: " << this->getTag() << " Type: PDeltaCrdTransf3d" << endln;
   if (nodeIOffset)
	   s << "\tNode I offset: " << nodeIOffset[0] << " " << nodeIOffset[1] << " "<< nodeIOffset[2] << endln;
   if (nodeJOffset)
	   s << "\tNode J offset: " << nodeJOffset[0] << " " << nodeJOffset[1] << " "<< nodeJOffset[2] << endln;
}

