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
// $Date: 2007-05-15 22:21:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/PDeltaCrdTransf2d.cpp,v $


// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Modified: 04/2005 Andreas Schellenberg (getBasicTrialVel, getBasicTrialAccel)
// 
// Purpose: This file contains the implementation for the 
// PDeltaCrdTransf2d class. PDeltaCrdTransf2d is a linear
// transformation for a planar frame between the global 
// and basic coordinate systems


#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>
#include <elementAPI.h>
#include <string>
#include <PDeltaCrdTransf2d.h>

// initialize static variables
Matrix PDeltaCrdTransf2d::Tlg(6,6);
Matrix PDeltaCrdTransf2d::kg(6,6);

void* OPS_PDeltaCrdTransf2d()
{
    if(OPS_GetNumRemainingInputArgs() < 1) {
	opserr<<"insufficient arguments for PDeltaCrdTransf2d\n";
	return 0;
    }

    // get tag
    int tag;
    int numData = 1;
    if(OPS_GetIntInput(&numData,&tag) < 0) return 0;

    // get option
    Vector jntOffsetI(2), jntOffsetJ(2);
    double *iptr=&jntOffsetI(0), *jptr=&jntOffsetJ(0);
    while(OPS_GetNumRemainingInputArgs() > 4) {
	std::string type = OPS_GetString();
	if(type == "-jntOffset") {
	    numData = 2;
	    if(OPS_GetDoubleInput(&numData,iptr) < 0) return 0;
	    if(OPS_GetDoubleInput(&numData,jptr) < 0) return 0;
	}
    }

    return new PDeltaCrdTransf2d(tag,jntOffsetI,jntOffsetJ);
}


// constructor:
PDeltaCrdTransf2d::PDeltaCrdTransf2d(int tag)
:CrdTransf(tag, CRDTR_TAG_PDeltaCrdTransf2d),
nodeIPtr(0), nodeJPtr(0),
nodeIOffset(0), nodeJOffset(0),
cosTheta(0), sinTheta(0), L(0), ul14(0),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    // Does nothing
}


// constructor:
PDeltaCrdTransf2d::PDeltaCrdTransf2d(int tag,
                                     const Vector &rigJntOffset1,
                                     const Vector &rigJntOffset2)
  :CrdTransf(tag, CRDTR_TAG_PDeltaCrdTransf2d),
   nodeIPtr(0), nodeJPtr(0),
   nodeIOffset(0), nodeJOffset(0),
   cosTheta(0), sinTheta(0), L(0), ul14(0),
   nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    // check rigid joint offset for node I
    if (&rigJntOffset1 == 0 || rigJntOffset1.Size() != 2 ) {
        opserr << "PDeltaCrdTransf2d::PDeltaCrdTransf2d:  Invalid rigid joint offset vector for node I\n";
        opserr << "Size must be 2\n";      
    }
    else if (rigJntOffset1.Norm() > 0.0) {
        nodeIOffset = new double[2];
        nodeIOffset[0] = rigJntOffset1(0);
        nodeIOffset[1] = rigJntOffset1(1);
    }
    
    // check rigid joint offset for node J
    if (&rigJntOffset2 == 0 || rigJntOffset2.Size() != 2 ) {
        opserr << "PDeltaCrdTransf2d::PDeltaCrdTransf2d:  Invalid rigid joint offset vector for node J\n";
        opserr << "Size must be 2\n";      
    }
    else if (rigJntOffset2.Norm() > 0.0) {
        nodeJOffset = new double[2];
        nodeJOffset[0] = rigJntOffset2(0);
        nodeJOffset[1] = rigJntOffset2(1);
    }
}


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
PDeltaCrdTransf2d::PDeltaCrdTransf2d()
:CrdTransf(0, CRDTR_TAG_PDeltaCrdTransf2d),
nodeIPtr(0), nodeJPtr(0),
nodeIOffset(0), nodeJOffset(0),
cosTheta(0), sinTheta(0), L(0), ul14(0),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    
}


// destructor:
PDeltaCrdTransf2d::~PDeltaCrdTransf2d() 
{
    if (nodeIOffset)
        delete [] nodeIOffset;
    if (nodeJOffset)
        delete [] nodeJOffset;
    if (nodeIInitialDisp != 0)
        delete [] nodeIInitialDisp;
    if (nodeJInitialDisp != 0)
        delete [] nodeJInitialDisp;
}


int
PDeltaCrdTransf2d::commitState(void)
{
    return 0;
}


int
PDeltaCrdTransf2d::revertToLastCommit(void)
{
    return 0;
}


int
PDeltaCrdTransf2d::revertToStart(void)
{
    return 0;
}


int 
PDeltaCrdTransf2d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
    int error;
    
    nodeIPtr = nodeIPointer;
    nodeJPtr = nodeJPointer;
    
    if ((!nodeIPtr) || (!nodeJPtr))
    {
        opserr << "\nPDeltaCrdTransf2d::initialize";
        opserr << "\ninvalid pointers to the element nodes\n";
        return -1;
    }
    
    // see if there is some initial displacements at nodes
    if (initialDispChecked == false) {
        
        const Vector &nodeIDisp = nodeIPtr->getDisp();
        const Vector &nodeJDisp = nodeJPtr->getDisp();
        for (int i=0; i<3; i++)
            if (nodeIDisp(i) != 0.0) {
                nodeIInitialDisp = new double [3];
                for (int j=0; j<3; j++)
                    nodeIInitialDisp[j] = nodeIDisp(j);
                i = 3;
            }
            for (int j=0; j<3; j++)
                if (nodeJDisp(j) != 0.0) {
                    nodeJInitialDisp = new double [3];
                    for (int i=0; i<3; i++)
                        nodeJInitialDisp[i] = nodeJDisp(i);
                    j = 3;
                }
                
                initialDispChecked = true;
    }
    
    // get element length and orientation
    if ((error = this->computeElemtLengthAndOrient()))
        return error;
    
    
    return 0;
}


int
PDeltaCrdTransf2d::update(void)
{
    static Vector nodeIDisp(3);
    static Vector nodeJDisp(3);
    nodeIDisp = nodeIPtr->getTrialDisp();
    nodeJDisp = nodeJPtr->getTrialDisp();
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<3; j++)
            nodeIDisp(j) -= nodeIInitialDisp[j];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<3; j++)
            nodeJDisp(j) -= nodeJInitialDisp[j];
    }
    
    double ul1;
    double ul4;
    
    ul1 = -sinTheta*nodeIDisp(0) + cosTheta*nodeIDisp(1);
    ul4 = -sinTheta*nodeJDisp(0) + cosTheta*nodeJDisp(1);
    
    if (nodeIOffset != 0) {
        double t12 = sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
        ul1 += t12*nodeIDisp(2);
    }
    
    if (nodeJOffset != 0) {
        double t45 = sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
        ul4 += t45*nodeJDisp(2);
    }
    
    ul14 = ul1-ul4;
    
    return 0;
}


int 
PDeltaCrdTransf2d::computeElemtLengthAndOrient()
{
    // element projection
    static Vector dx(2);
    
    const Vector &ndICoords = nodeIPtr->getCrds();
    const Vector &ndJCoords = nodeJPtr->getCrds();
    
    dx(0) = ndJCoords(0) - ndICoords(0);
    dx(1) = ndJCoords(1) - ndICoords(1);
    
    if (nodeIInitialDisp != 0) {
        dx(0) -= nodeIInitialDisp[0];
        dx(1) -= nodeIInitialDisp[1];
    }
    
    if (nodeJInitialDisp != 0) {
        dx(0) += nodeJInitialDisp[0];
        dx(1) += nodeJInitialDisp[1];
    }
    
    if (nodeJOffset != 0) {
        dx(0) += nodeJOffset[0];
        dx(1) += nodeJOffset[1];
    }
    
    if (nodeIOffset != 0) {
        dx(0) -= nodeIOffset[0];
        dx(1) -= nodeIOffset[1];
    }
    
    // calculate the element length
    L = dx.Norm();
    
    if (L == 0.0) 
    {
        opserr << "\nPDeltaCrdTransf2d::computeElemtLengthAndOrien: 0 length\n";
        return -2;  
    }
    
    // calculate the element local x axis components (direction cosines)
    // wrt to the global coordinates 
    cosTheta = dx(0)/L;
    sinTheta = dx(1)/L;
    
    return 0;
}


void
PDeltaCrdTransf2d::compTransfMatrixLocalGlobal(Matrix &Tlg) 
{
    // setup transformation matrix from global to local coordinates
    Tlg.Zero();
    
    Tlg(0,0) = Tlg(3,3) =  cosTheta;
    Tlg(0,1) = Tlg(3,4) =  sinTheta;
    Tlg(1,0) = Tlg(4,3) = -sinTheta;
    Tlg(1,1) = Tlg(4,4) =  cosTheta;
    Tlg(2,2) = Tlg(5,5) =  1.0;
}


double 
PDeltaCrdTransf2d::getInitialLength(void)
{
    return L;
}


double 
PDeltaCrdTransf2d::getDeformedLength(void)
{
    return L;
}


const Vector &
PDeltaCrdTransf2d::getBasicTrialDisp(void)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getTrialDisp();
    const Vector &disp2 = nodeJPtr->getTrialDisp();
    
    static double ug[6];
    for (int i = 0; i < 3; i++) {
        ug[i]   = disp1(i);
        ug[i+3] = disp2(i);
    }
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<3; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<3; j++)
            ug[j+3] -= nodeJInitialDisp[j];
    }
    
    static Vector ub(3);
    
    double oneOverL = 1.0/L;
    double sl = sinTheta*oneOverL;
    double cl = cosTheta*oneOverL;
    
    ub(0) = -cosTheta*ug[0] - sinTheta*ug[1] +
        cosTheta*ug[3] + sinTheta*ug[4];
    
    ub(1) = -sl*ug[0] + cl*ug[1] + ug[2] +
        sl*ug[3] - cl*ug[4];
    
    if (nodeIOffset != 0) {
        double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
        double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
        ub(0) -= t02*ug[2];
        ub(1) += oneOverL*t12*ug[2];
    }
    
    if (nodeJOffset != 0) {
        double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
        double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
        ub(0) += t35*ug[5];
        ub(1) -= oneOverL*t45*ug[5];
    }
    
    ub(2) = ub(1) + ug[5] - ug[2];
    
    return ub;
}


const Vector &
PDeltaCrdTransf2d::getBasicIncrDisp(void)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getIncrDisp();
    const Vector &disp2 = nodeJPtr->getIncrDisp();
    
    static double dug[6];
    for (int i = 0; i < 3; i++) {
        dug[i]   = disp1(i);
        dug[i+3] = disp2(i);
    }
    
    static Vector dub(3);
    
    double oneOverL = 1.0/L;
    double sl = sinTheta*oneOverL;
    double cl = cosTheta*oneOverL;
    
    dub(0) = -cosTheta*dug[0] - sinTheta*dug[1] +
        cosTheta*dug[3] + sinTheta*dug[4];
    
    dub(1) = -sl*dug[0] + cl*dug[1] + dug[2] +
        sl*dug[3] - cl*dug[4];
    
    if (nodeIOffset != 0) {
        double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
        double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
        dub(0) -= t02*dug[2];
        dub(1) += oneOverL*t12*dug[2];
    }
    
    if (nodeJOffset != 0) {
        double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
        double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
        dub(0) += t35*dug[5];
        dub(1) -= oneOverL*t45*dug[5];
    }
    
    dub(2) = dub(1) + dug[5] - dug[2];
    
    return dub;
}


const Vector &
PDeltaCrdTransf2d::getBasicIncrDeltaDisp(void)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getIncrDeltaDisp();
    const Vector &disp2 = nodeJPtr->getIncrDeltaDisp();
    
    static double Dug[6];
    for (int i = 0; i < 3; i++) {
        Dug[i]   = disp1(i);
        Dug[i+3] = disp2(i);
    }
    
    static Vector Dub(3);
    
    double oneOverL = 1.0/L;
    double sl = sinTheta*oneOverL;
    double cl = cosTheta*oneOverL;
    
    Dub(0) = -cosTheta*Dug[0] - sinTheta*Dug[1] +
        cosTheta*Dug[3] + sinTheta*Dug[4];
    
    Dub(1) = -sl*Dug[0] + cl*Dug[1] + Dug[2] +
        sl*Dug[3] - cl*Dug[4];
    
    if (nodeIOffset != 0) {
        double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
        double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
        Dub(0) -= t02*Dug[2];
        Dub(1) += oneOverL*t12*Dug[2];
    }
    
    if (nodeJOffset != 0) {
        double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
        double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
        Dub(0) += t35*Dug[5];
        Dub(1) -= oneOverL*t45*Dug[5];
    }
    
    Dub(2) = Dub(1) + Dug[5] - Dug[2];
    
    return Dub;
}


const Vector &
PDeltaCrdTransf2d::getBasicTrialVel(void)
{
	// determine global velocities
	const Vector &vel1 = nodeIPtr->getTrialVel();
	const Vector &vel2 = nodeJPtr->getTrialVel();
	
	static double vg[6];
	for (int i = 0; i < 3; i++) {
		vg[i]   = vel1(i);
		vg[i+3] = vel2(i);
	}
	
	static Vector vb(3);
	
	double oneOverL = 1.0/L;
	double sl = sinTheta*oneOverL;
	double cl = cosTheta*oneOverL;
	
	vb(0) = -cosTheta*vg[0] - sinTheta*vg[1] +
		cosTheta*vg[3] + sinTheta*vg[4];
	
	vb(1) = -sl*vg[0] + cl*vg[1] + vg[2] +
		sl*vg[3] - cl*vg[4];
	
	if (nodeIOffset != 0) {
		double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
		double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
		vb(0) -= t02*vg[2];
		vb(1) += oneOverL*t12*vg[2];
	}
	
	if (nodeJOffset != 0) {
		double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
		double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
		vb(0) += t35*vg[5];
		vb(1) -= oneOverL*t45*vg[5];
	}
	
	vb(2) = vb(1) + vg[5] - vg[2];
	
	return vb;
}


const Vector &
PDeltaCrdTransf2d::getBasicTrialAccel(void)
{
	// determine global accelerations
	const Vector &accel1 = nodeIPtr->getTrialAccel();
	const Vector &accel2 = nodeJPtr->getTrialAccel();
	
	static double ag[6];
	for (int i = 0; i < 3; i++) {
		ag[i]   = accel1(i);
		ag[i+3] = accel2(i);
	}
	
	static Vector ab(3);
	
	double oneOverL = 1.0/L;
	double sl = sinTheta*oneOverL;
	double cl = cosTheta*oneOverL;
	
	ab(0) = -cosTheta*ag[0] - sinTheta*ag[1] +
		cosTheta*ag[3] + sinTheta*ag[4];
	
	ab(1) = -sl*ag[0] + cl*ag[1] + ag[2] +
		sl*ag[3] - cl*ag[4];
	
	if (nodeIOffset != 0) {
		double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
		double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
		ab(0) -= t02*ag[2];
		ab(1) += oneOverL*t12*ag[2];
	}
	
	if (nodeJOffset != 0) {
		double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
		double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
		ab(0) += t35*ag[5];
		ab(1) -= oneOverL*t45*ag[5];
	}
	
	ab(2) = ab(1) + ag[5] - ag[2];
	
	return ab;
}


const Vector &
PDeltaCrdTransf2d::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
    // transform resisting forces from the basic system to local coordinates
    static double pl[6];
    
    double q0 = pb(0);
    double q1 = pb(1);
    double q2 = pb(2);
    
    double oneOverL = 1.0/L;
    
    double V = oneOverL*(q1+q2);
    pl[0] = -q0;
    pl[1] =  V;
    pl[2] =  q1;
    pl[3] =  q0;
    pl[4] = -V;
    pl[5] =  q2;
    
    // add end forces due to element p0 loads
    pl[0] += p0(0);
    pl[1] += p0(1);
    pl[4] += p0(2);
    
    // Include leaning column effects (P-Delta)
    double NoverL = ul14*q0*oneOverL;             
    pl[1] += NoverL;
    pl[4] -= NoverL;
    
    // transform resisting forces  from local to global coordinates
    static Vector pg(6);
    
    pg(0) = cosTheta*pl[0] - sinTheta*pl[1];
    pg(1) = sinTheta*pl[0] + cosTheta*pl[1];
    
    pg(3) = cosTheta*pl[3] - sinTheta*pl[4];
    pg(4) = sinTheta*pl[3] + cosTheta*pl[4];
    
    pg(2) = pl[2];
    pg(5) = pl[5];	
    
    if (nodeIOffset != 0) {
        double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
        double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
        
        pg(2) += t02*pl[0] + t12*pl[1];
    }
    
    if (nodeJOffset != 0) {
        double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
        double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
        
        pg(5) += t35*pl[3] + t45*pl[4];
    }
    
    return pg;
}


const Matrix &
PDeltaCrdTransf2d::getGlobalStiffMatrix(const Matrix &kb, const Vector &pb)
{
    static double kl[6][6];
    static double tmp[6][6];
    double oneOverL = 1.0/L;
    
    // Basic stiffness
    double kb00, kb01, kb02, kb10, kb11, kb12, kb20, kb21, kb22;
    kb00 = kb(0,0);		kb01 = kb(0,1);		kb02 = kb(0,2);
    kb10 = kb(1,0);		kb11 = kb(1,1);		kb12 = kb(1,2);
    kb20 = kb(2,0);		kb21 = kb(2,1);		kb22 = kb(2,2);
    
    // Transform basic stiffness to local system
    kl[0][0] =  kb00;
    kl[1][0] = -oneOverL*(kb10+kb20);
    kl[2][0] = -kb10;
    kl[3][0] = -kb00;
    kl[4][0] = -kl[1][0];
    kl[5][0] = -kb20;
    
    kl[0][1] = -oneOverL*(kb01+kb02);
    kl[1][1] =  oneOverL*oneOverL*(kb11+kb12+kb21+kb22);
    kl[2][1] =  oneOverL*(kb11+kb12);
    kl[3][1] = -kl[0][1];
    kl[4][1] = -kl[1][1];
    kl[5][1] =  oneOverL*(kb21+kb22);
    
    kl[0][2] = -kb01;
    kl[1][2] =  oneOverL*(kb11+kb21);
    kl[2][2] =  kb11;
    kl[3][2] =  kb01;
    kl[4][2] = -kl[1][2];
    kl[5][2] =  kb21;
    
    kl[0][3] = -kl[0][0];
    kl[1][3] = -kl[1][0];
    kl[2][3] = -kl[2][0];
    kl[3][3] = -kl[3][0];
    kl[4][3] = -kl[4][0];
    kl[5][3] = -kl[5][0];
    
    kl[0][4] = -kl[0][1];
    kl[1][4] = -kl[1][1];
    kl[2][4] = -kl[2][1];
    kl[3][4] = -kl[3][1];
    kl[4][4] = -kl[4][1];
    kl[5][4] = -kl[5][1];
    
    kl[0][5] = -kb02;
    kl[1][5] =  oneOverL*(kb12+kb22);
    kl[2][5] =  kb12;
    kl[3][5] =  kb02;
    kl[4][5] = -kl[1][5];
    kl[5][5] =  kb22;
    
    // Include geometric stiffness effects in local system
    double NoverL = pb(0)*oneOverL;
    kl[1][1] += NoverL;
    kl[4][4] += NoverL;
    kl[1][4] -= NoverL;
    kl[4][1] -= NoverL;
    
    double t02 = 0.0;
    double t12 = 0.0;
    double t35 = 0.0;
    double t45 = 0.0;
    
    if (nodeIOffset != 0) {
        t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
        t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
    }
    
    if (nodeJOffset != 0) {
        t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
        t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
    }
    
    // Now transform from local to global ... compute kl*T
    tmp[0][0] = kl[0][0]*cosTheta - kl[0][1]*sinTheta;
    tmp[1][0] = kl[1][0]*cosTheta - kl[1][1]*sinTheta;
    tmp[2][0] = kl[2][0]*cosTheta - kl[2][1]*sinTheta;
    tmp[3][0] = kl[3][0]*cosTheta - kl[3][1]*sinTheta;
    tmp[4][0] = kl[4][0]*cosTheta - kl[4][1]*sinTheta;
    tmp[5][0] = kl[5][0]*cosTheta - kl[5][1]*sinTheta;
    
    tmp[0][1] = kl[0][0]*sinTheta + kl[0][1]*cosTheta;
    tmp[1][1] = kl[1][0]*sinTheta + kl[1][1]*cosTheta;
    tmp[2][1] = kl[2][0]*sinTheta + kl[2][1]*cosTheta;
    tmp[3][1] = kl[3][0]*sinTheta + kl[3][1]*cosTheta;
    tmp[4][1] = kl[4][0]*sinTheta + kl[4][1]*cosTheta;
    tmp[5][1] = kl[5][0]*sinTheta + kl[5][1]*cosTheta;
    
    if (nodeIOffset) {
        tmp[0][2] = kl[0][0]*t02 + kl[0][1]*t12 + kl[0][2];
        tmp[1][2] = kl[1][0]*t02 + kl[1][1]*t12 + kl[1][2];
        tmp[2][2] = kl[2][0]*t02 + kl[2][1]*t12 + kl[2][2];
        tmp[3][2] = kl[3][0]*t02 + kl[3][1]*t12 + kl[3][2];
        tmp[4][2] = kl[4][0]*t02 + kl[4][1]*t12 + kl[4][2];
        tmp[5][2] = kl[5][0]*t02 + kl[5][1]*t12 + kl[5][2];
    }
    else {
        tmp[0][2] = kl[0][2];
        tmp[1][2] = kl[1][2];
        tmp[2][2] = kl[2][2];
        tmp[3][2] = kl[3][2];
        tmp[4][2] = kl[4][2];
        tmp[5][2] = kl[5][2];
    }
    
    tmp[0][3] = kl[0][3]*cosTheta - kl[0][4]*sinTheta;
    tmp[1][3] = kl[1][3]*cosTheta - kl[1][4]*sinTheta;
    tmp[2][3] = kl[2][3]*cosTheta - kl[2][4]*sinTheta;
    tmp[3][3] = kl[3][3]*cosTheta - kl[3][4]*sinTheta;
    tmp[4][3] = kl[4][3]*cosTheta - kl[4][4]*sinTheta;
    tmp[5][3] = kl[5][3]*cosTheta - kl[5][4]*sinTheta;
    
    tmp[0][4] = kl[0][3]*sinTheta + kl[0][4]*cosTheta;
    tmp[1][4] = kl[1][3]*sinTheta + kl[1][4]*cosTheta;
    tmp[2][4] = kl[2][3]*sinTheta + kl[2][4]*cosTheta;
    tmp[3][4] = kl[3][3]*sinTheta + kl[3][4]*cosTheta;
    tmp[4][4] = kl[4][3]*sinTheta + kl[4][4]*cosTheta;
    tmp[5][4] = kl[5][3]*sinTheta + kl[5][4]*cosTheta;
    
    if (nodeJOffset) {
        tmp[0][5] = kl[0][3]*t35 + kl[0][4]*t45 + kl[0][5];
        tmp[1][5] = kl[1][3]*t35 + kl[1][4]*t45 + kl[1][5];
        tmp[2][5] = kl[2][3]*t35 + kl[2][4]*t45 + kl[2][5];
        tmp[3][5] = kl[3][3]*t35 + kl[3][4]*t45 + kl[3][5];
        tmp[4][5] = kl[4][3]*t35 + kl[4][4]*t45 + kl[4][5];
        tmp[5][5] = kl[5][3]*t35 + kl[5][4]*t45 + kl[5][5];
    }
    else {
        tmp[0][5] = kl[0][5];
        tmp[1][5] = kl[1][5];
        tmp[2][5] = kl[2][5];
        tmp[3][5] = kl[3][5];
        tmp[4][5] = kl[4][5];
        tmp[5][5] = kl[5][5];
    }
    
    // Now compute T'*(kl*T)
    kg(0,0) = cosTheta*tmp[0][0] - sinTheta*tmp[1][0];
    kg(0,1) = cosTheta*tmp[0][1] - sinTheta*tmp[1][1];
    kg(0,2) = cosTheta*tmp[0][2] - sinTheta*tmp[1][2];
    kg(0,3) = cosTheta*tmp[0][3] - sinTheta*tmp[1][3];
    kg(0,4) = cosTheta*tmp[0][4] - sinTheta*tmp[1][4];
    kg(0,5) = cosTheta*tmp[0][5] - sinTheta*tmp[1][5];
    
    kg(1,0) = sinTheta*tmp[0][0] + cosTheta*tmp[1][0];
    kg(1,1) = sinTheta*tmp[0][1] + cosTheta*tmp[1][1];
    kg(1,2) = sinTheta*tmp[0][2] + cosTheta*tmp[1][2];
    kg(1,3) = sinTheta*tmp[0][3] + cosTheta*tmp[1][3];
    kg(1,4) = sinTheta*tmp[0][4] + cosTheta*tmp[1][4];
    kg(1,5) = sinTheta*tmp[0][5] + cosTheta*tmp[1][5];
    
    if (nodeIOffset) {
        kg(2,0) = t02*tmp[0][0] + t12*tmp[1][0] + tmp[2][0];
        kg(2,1) = t02*tmp[0][1] + t12*tmp[1][1] + tmp[2][1];
        kg(2,2) = t02*tmp[0][2] + t12*tmp[1][2] + tmp[2][2];
        kg(2,3) = t02*tmp[0][3] + t12*tmp[1][3] + tmp[2][3];
        kg(2,4) = t02*tmp[0][4] + t12*tmp[1][4] + tmp[2][4];
        kg(2,5) = t02*tmp[0][5] + t12*tmp[1][5] + tmp[2][5];
    }
    else {
        kg(2,0) = tmp[2][0];
        kg(2,1) = tmp[2][1];
        kg(2,2) = tmp[2][2];
        kg(2,3) = tmp[2][3];
        kg(2,4) = tmp[2][4];
        kg(2,5) = tmp[2][5];
    }
    
    kg(3,0) = cosTheta*tmp[3][0] - sinTheta*tmp[4][0];
    kg(3,1) = cosTheta*tmp[3][1] - sinTheta*tmp[4][1];
    kg(3,2) = cosTheta*tmp[3][2] - sinTheta*tmp[4][2];
    kg(3,3) = cosTheta*tmp[3][3] - sinTheta*tmp[4][3];
    kg(3,4) = cosTheta*tmp[3][4] - sinTheta*tmp[4][4];
    kg(3,5) = cosTheta*tmp[3][5] - sinTheta*tmp[4][5];
    
    kg(4,0) = sinTheta*tmp[3][0] + cosTheta*tmp[4][0];
    kg(4,1) = sinTheta*tmp[3][1] + cosTheta*tmp[4][1];
    kg(4,2) = sinTheta*tmp[3][2] + cosTheta*tmp[4][2];
    kg(4,3) = sinTheta*tmp[3][3] + cosTheta*tmp[4][3];
    kg(4,4) = sinTheta*tmp[3][4] + cosTheta*tmp[4][4];
    kg(4,5) = sinTheta*tmp[3][5] + cosTheta*tmp[4][5];
    
    if (nodeJOffset) {
        kg(5,0) = t35*tmp[3][0] + t45*tmp[4][0] + tmp[5][0];
        kg(5,1) = t35*tmp[3][1] + t45*tmp[4][1] + tmp[5][1];
        kg(5,2) = t35*tmp[3][2] + t45*tmp[4][2] + tmp[5][2];
        kg(5,3) = t35*tmp[3][3] + t45*tmp[4][3] + tmp[5][3];
        kg(5,4) = t35*tmp[3][4] + t45*tmp[4][4] + tmp[5][4];
        kg(5,5) = t35*tmp[3][5] + t45*tmp[4][5] + tmp[5][5];
    }
    else {
        kg(5,0) = tmp[5][0];
        kg(5,1) = tmp[5][1];
        kg(5,2) = tmp[5][2];
        kg(5,3) = tmp[5][3];
        kg(5,4) = tmp[5][4];
        kg(5,5) = tmp[5][5];
    }
    
    return kg;
}


const Matrix &
PDeltaCrdTransf2d::getInitialGlobalStiffMatrix(const Matrix &kb)
{
    static double tmp [6][6];
    double oneOverL = 1.0/L;
    double kb00, kb01, kb02, kb10, kb11, kb12, kb20, kb21, kb22;
    
    kb00 = kb(0,0);		kb01 = kb(0,1);		kb02 = kb(0,2);
    kb10 = kb(1,0);		kb11 = kb(1,1);		kb12 = kb(1,2);
    kb20 = kb(2,0);		kb21 = kb(2,1);		kb22 = kb(2,2);
    
    double t02 = 0.0;
    double t12 = 1.0;
    double t22 = 0.0;
    
    if (nodeIOffset != 0) {
        t02 =  cosTheta*nodeIOffset[1] - sinTheta*nodeIOffset[0];
        t12 =  oneOverL*(sinTheta*nodeIOffset[1]+cosTheta*nodeIOffset[0]) + 1.0;
        t22 =  oneOverL*(sinTheta*nodeIOffset[1]+cosTheta*nodeIOffset[0]);
    }
    
    double t05 = 0.0;
    double t15 = 0.0;
    double t25 = 1.0;
    
    if (nodeJOffset != 0) {
        t05 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
        t15 = -oneOverL*(sinTheta*nodeJOffset[1]+cosTheta*nodeJOffset[0]);
        t25 = -oneOverL*(sinTheta*nodeJOffset[1]+cosTheta*nodeJOffset[0]) + 1.0;
    }
    
    double sl = sinTheta*oneOverL;
    double cl = cosTheta*oneOverL;
    
    tmp[0][0] = -cosTheta*kb00 - sl*(kb01+kb02);
    tmp[0][1] = -sinTheta*kb00 + cl*(kb01+kb02);
    tmp[0][2] = (nodeIOffset) ? t02*kb00 + t12*kb01 + t22*kb02 : kb01;
    tmp[0][3] = -tmp[0][0];
    tmp[0][4] = -tmp[0][1];
    tmp[0][5] = (nodeJOffset) ? t05*kb00 + t15*kb01 + t25*kb02 : kb02;
    
    tmp[1][0] = -cosTheta*kb10 - sl*(kb11+kb12);
    tmp[1][1] = -sinTheta*kb10 + cl*(kb11+kb12);
    tmp[1][2] = (nodeIOffset) ? t02*kb10 + t12*kb11 + t22*kb12 : kb11;
    tmp[1][3] = -tmp[1][0];
    tmp[1][4] = -tmp[1][1];
    tmp[1][5] = (nodeJOffset) ? t05*kb10 + t15*kb11 + t25*kb12 : kb12;
    
    tmp[2][0] = -cosTheta*kb20 - sl*(kb21+kb22);
    tmp[2][1] = -sinTheta*kb20 + cl*(kb21+kb22);
    tmp[2][2] = (nodeIOffset) ? t02*kb20 + t12*kb21 + t22*kb22 : kb21;
    tmp[2][3] = -tmp[2][0];
    tmp[2][4] = -tmp[2][1];
    tmp[2][5] = (nodeJOffset) ? t05*kb20 + t15*kb21 + t25*kb22 : kb22;
    
    kg(0,0) = -cosTheta*tmp[0][0] - sl*(tmp[1][0]+tmp[2][0]);
    kg(0,1) = -cosTheta*tmp[0][1] - sl*(tmp[1][1]+tmp[2][1]);
    kg(0,2) = -cosTheta*tmp[0][2] - sl*(tmp[1][2]+tmp[2][2]);
    kg(0,3) = -cosTheta*tmp[0][3] - sl*(tmp[1][3]+tmp[2][3]);
    kg(0,4) = -cosTheta*tmp[0][4] - sl*(tmp[1][4]+tmp[2][4]);
    kg(0,5) = -cosTheta*tmp[0][5] - sl*(tmp[1][5]+tmp[2][5]);
    
    kg(1,0) = -sinTheta*tmp[0][0] + cl*(tmp[1][0]+tmp[2][0]);
    kg(1,1) = -sinTheta*tmp[0][1] + cl*(tmp[1][1]+tmp[2][1]);
    kg(1,2) = -sinTheta*tmp[0][2] + cl*(tmp[1][2]+tmp[2][2]);
    kg(1,3) = -sinTheta*tmp[0][3] + cl*(tmp[1][3]+tmp[2][3]);
    kg(1,4) = -sinTheta*tmp[0][4] + cl*(tmp[1][4]+tmp[2][4]);
    kg(1,5) = -sinTheta*tmp[0][5] + cl*(tmp[1][5]+tmp[2][5]);
    
    if (nodeIOffset) {
        kg(2,0) =  t02*tmp[0][0] + t12*tmp[1][0] + t22*tmp[2][0];
        kg(2,1) =  t02*tmp[0][1] + t12*tmp[1][1] + t22*tmp[2][1];
        kg(2,2) =  t02*tmp[0][2] + t12*tmp[1][2] + t22*tmp[2][2];
        kg(2,3) =  t02*tmp[0][3] + t12*tmp[1][3] + t22*tmp[2][3];
        kg(2,4) =  t02*tmp[0][4] + t12*tmp[1][4] + t22*tmp[2][4];
        kg(2,5) =  t02*tmp[0][5] + t12*tmp[1][5] + t22*tmp[2][5];
    }
    else {
        kg(2,0) = tmp[1][0];
        kg(2,1) = tmp[1][1];
        kg(2,2) = tmp[1][2];
        kg(2,3) = tmp[1][3];
        kg(2,4) = tmp[1][4];
        kg(2,5) = tmp[1][5];
    }
    
    kg(3,0) = -kg(0,0);
    kg(3,1) = -kg(0,1);
    kg(3,2) = -kg(0,2);
    kg(3,3) = -kg(0,3);
    kg(3,4) = -kg(0,4);
    kg(3,5) = -kg(0,5);
    
    kg(4,0) = -kg(1,0);
    kg(4,1) = -kg(1,1);
    kg(4,2) = -kg(1,2);
    kg(4,3) = -kg(1,3);
    kg(4,4) = -kg(1,4);
    kg(4,5) = -kg(1,5);
    
    if (nodeJOffset) {
        kg(5,0) =  t05*tmp[0][0] + t15*tmp[1][0] + t25*tmp[2][0];
        kg(5,1) =  t05*tmp[0][1] + t15*tmp[1][1] + t25*tmp[2][1];
        kg(5,2) =  t05*tmp[0][2] + t15*tmp[1][2] + t25*tmp[2][2];
        kg(5,3) =  t05*tmp[0][3] + t15*tmp[1][3] + t25*tmp[2][3];
        kg(5,4) =  t05*tmp[0][4] + t15*tmp[1][4] + t25*tmp[2][4];
        kg(5,5) =  t05*tmp[0][5] + t15*tmp[1][5] + t25*tmp[2][5];
    }
    else {
        kg(5,0) =  tmp[2][0];
        kg(5,1) =  tmp[2][1];
        kg(5,2) =  tmp[2][2];
        kg(5,3) =  tmp[2][3];
        kg(5,4) =  tmp[2][4];
        kg(5,5) =  tmp[2][5];
    }
    
    return kg;
}


CrdTransf *
PDeltaCrdTransf2d::getCopy2d(void)
{
    // create a new instance of PDeltaCrdTransf2d 
    
    PDeltaCrdTransf2d *theCopy;
    
    Vector offsetI(2);
    Vector offsetJ(2);
    
    if (nodeIOffset != 0) {
        offsetI(0) = nodeIOffset[0];
        offsetI(1) = nodeIOffset[1];
    }
    
    if (nodeJOffset != 0) {
        offsetJ(0) = nodeJOffset[0];
        offsetJ(1) = nodeJOffset[1];
    }
    
    theCopy = new PDeltaCrdTransf2d(this->getTag(), offsetI, offsetJ);
    
    theCopy->nodeIPtr = nodeIPtr;
    theCopy->nodeJPtr = nodeJPtr;
    theCopy->cosTheta = cosTheta;
    theCopy->sinTheta = sinTheta;
    theCopy->L = L;
    theCopy->ul14 = ul14;
    
    return theCopy;
}


int 
PDeltaCrdTransf2d::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    
    static Vector data(12);
    data(0) = this->getTag();
    data(1) = L;
    if (nodeIOffset != 0) {
        data(2) = nodeIOffset[0];
        data(3) = nodeIOffset[1];
    } else {
        data(2) = 0.0;
        data(3) = 0.0;
    }
    
    if (nodeJOffset != 0) {
        data(4) = nodeJOffset[0];
        data(5) = nodeJOffset[1];
    } else {
        data(4) = 0.0;
        data(5) = 0.0;
    }
    
    if (nodeIInitialDisp != 0) {
        data(6) = nodeIInitialDisp[0];
        data(7) = nodeIInitialDisp[1];
        data(8) = nodeIInitialDisp[2];
    } else {
        data(6) = 0.0;
        data(7) = 0.0;
        data(8) = 0.0;
    }
    
    if (nodeJInitialDisp != 0) {
        data(9) = nodeJInitialDisp[0];
        data(10) = nodeJInitialDisp[1];
        data(11) = nodeJInitialDisp[2];
    } else {
        data(9) = 0.0;
        data(10) = 0.0;
        data(11) = 0.0;
    }
    
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
        opserr << "PDeltaCrdTransf2d2d::sendSelf - failed to send Vector\n";
        return res;
    }
    
    return res;
}


int 
PDeltaCrdTransf2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    
    static Vector data(12);
    
    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
        opserr << "PDeltaCrdTransf2d2d::recvSelf - failed to receive Vector\n";
        
        return res;
    }
    
    
    this->setTag((int)data(0));
    L = data(1);
    data(0) = this->getTag();
    data(1) = L;
    
    int flag;
    int i,j;
    
    flag = 0;
    for (i=2; i<=3; i++)
        if (data(i) != 0.0)
            flag = 1;
        if (flag == 1) {
            if (nodeIOffset == 0)
                nodeIOffset = new double[2];
            for (i=2, j=0; i<=3; i++, j++)
                nodeIOffset[j] = data(i);
        }
        
        flag = 0;
        for (i=4; i<=5; i++)
            if (data(i) != 0.0)
                flag = 1;
            if (flag == 1) {
                if (nodeJOffset == 0)
                    nodeJOffset = new double[2];
                for (i=4, j=0; i<=5; i++, j++)
                    nodeJOffset[j] = data(i);
            }
            
            flag = 0;
            for (i=6; i<=8; i++)
                if (data(i) != 0.0)
                    flag = 1;
                if (flag == 1) {
                    if (nodeIInitialDisp == 0)
                        nodeIInitialDisp = new double[3];
                    for (i=6, j=0; i<=7; i++, j++)
                        nodeIInitialDisp[j] = data(i);
                }
                
                flag = 0;
                for (i=9; i<=11; i++)
                    if (data(i) != 0.0)
                        flag = 1;
                    if (flag == 1) {
                        if (nodeJInitialDisp == 0)
                            nodeJInitialDisp = new double [3];
                        for (i=9, j=0; i<=11; i++, j++)
                            nodeJInitialDisp[j] = data(i);
                    }
                    
                    initialDispChecked = true;
                    
                    return res;
}


const Matrix &
PDeltaCrdTransf2d::getGlobalMatrixFromLocal(const Matrix &ml)
{
    this->compTransfMatrixLocalGlobal(Tlg);  // OPTIMIZE LATER
    kg.addMatrixTripleProduct(0.0, Tlg, ml, 1.0);  // OPTIMIZE LATER

    return kg;
}


const Vector &
PDeltaCrdTransf2d::getPointGlobalCoordFromLocal(const Vector &xl)
{
    static Vector xg(2);
    
    const Vector &nodeICoords = nodeIPtr->getCrds();
    xg(0) = nodeICoords(0);
    xg(1) = nodeICoords(1);
    
    if (nodeIOffset) {
        xg(0) += nodeIOffset[0];
        xg(1) += nodeIOffset[1];
    }
    
    // xg = xg + Rlj'*xl
    xg(0) += cosTheta*xl(0) - sinTheta*xl(1);
    xg(1) += sinTheta*xl(0) + cosTheta*xl(1);
    
    return xg;  
}


const Vector &
PDeltaCrdTransf2d::getPointGlobalDisplFromBasic(double xi, const Vector &uxb)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getTrialDisp();
    const Vector &disp2 = nodeJPtr->getTrialDisp();
    
    static Vector ug(6);
    for (int i = 0; i < 3; i++)
    {
        ug(i)   = disp1(i);
        ug(i+3) = disp2(i);
    }
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<3; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<3; j++)
            ug[j+3] -= nodeJInitialDisp[j];
    }
    
    // transform global end displacements to local coordinates
    static Vector ul(6);      // total displacements
    
    ul(0) =  cosTheta*ug(0) + sinTheta*ug(1);
    ul(1) = -sinTheta*ug(0) + cosTheta*ug(1);
    ul(2) =  ug(2);
    ul(3) =  cosTheta*ug(3) + sinTheta*ug(4);
    ul(4) = -sinTheta*ug(3) + cosTheta*ug(4);
    ul(5) =  ug(5);
    
    if (nodeIOffset != 0) {
        double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
        double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
        
        ul(0) += t02*ug(2);
        ul(1) += t12*ug(2);
    }
    
    if (nodeJOffset != 0) {
        double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
        double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
        
        ul(3) += t35*ug(5);
        ul(4) += t45*ug(5);
    }
    
    // compute displacements at point xi, in local coordinates
    static Vector uxl(2),  uxg(2);
    
    uxl(0) = uxb(0) +        ul(0);
    uxl(1) = uxb(1) + (1-xi)*ul(1) + xi*ul(4);
    
    // rotate displacements to global coordinates
    // uxg = RljT*uxl
    uxg(0) = cosTheta*uxl(0) - sinTheta*uxl(1);
    uxg(1) = sinTheta*uxl(0) + cosTheta*uxl(1);
    
    return uxg;  
}


const Vector &
PDeltaCrdTransf2d::getPointLocalDisplFromBasic(double xi, const Vector &uxb)
{
    // determine global displacements
    const Vector &disp1 = nodeIPtr->getTrialDisp();
    const Vector &disp2 = nodeJPtr->getTrialDisp();
    
    static Vector ug(6);
    for (int i = 0; i < 3; i++)
    {
        ug(i)   = disp1(i);
        ug(i+3) = disp2(i);
    }
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<3; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<3; j++)
            ug[j+3] -= nodeJInitialDisp[j];
    }
    
    // transform global end displacements to local coordinates
    static Vector ul(6);      // total displacements
    
    ul(0) =  cosTheta*ug(0) + sinTheta*ug(1);
    ul(1) = -sinTheta*ug(0) + cosTheta*ug(1);
    ul(2) =  ug(2);
    ul(3) =  cosTheta*ug(3) + sinTheta*ug(4);
    ul(4) = -sinTheta*ug(3) + cosTheta*ug(4);
    ul(5) =  ug(5);
    
    if (nodeIOffset != 0) {
        double t02 = -cosTheta*nodeIOffset[1] + sinTheta*nodeIOffset[0];
        double t12 =  sinTheta*nodeIOffset[1] + cosTheta*nodeIOffset[0];
        
        ul(0) += t02*ug(2);
        ul(1) += t12*ug(2);
    }
    
    if (nodeJOffset != 0) {
        double t35 = -cosTheta*nodeJOffset[1] + sinTheta*nodeJOffset[0];
        double t45 =  sinTheta*nodeJOffset[1] + cosTheta*nodeJOffset[0];
        
        ul(3) += t35*ug(5);
        ul(4) += t45*ug(5);
    }
    
    // compute displacements at point xi, in local coordinates
    static Vector uxl(2);
    
    uxl(0) = uxb(0) +        ul(0);
    uxl(1) = uxb(1) + (1-xi)*ul(1) + xi*ul(4);
    
    return uxl;  
}


void
PDeltaCrdTransf2d::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "\nCrdTransf: " << this->getTag() << " Type: PDeltaCrdTransf2d";
		if (nodeIOffset != 0)
			s << "\tnodeI Offset: " << nodeIOffset[0] << ' ' << nodeIOffset[1] << endln;
		if (nodeJOffset != 0)
			s << "\tnodeJ Offset: " << nodeJOffset[0] << ' ' << nodeJOffset[1] << endln;
	}
	
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"PDeltaCrdTransf2d\"";
		if (nodeIOffset != 0)
			s << ", \"iOffset\": [" << nodeIOffset[0] << ", " << nodeIOffset[1] << "]";
		if (nodeJOffset != 0)
			s << ", \"jOffset\": [" << nodeJOffset[0] << ", " << nodeJOffset[1] << "]";
		s << "}";
	}
}

int
PDeltaCrdTransf2d::getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis)
{
  xAxis(0) = cosTheta;
  xAxis(1) = sinTheta;
  xAxis(2) = 0;

  yAxis(0) = -sinTheta;
  yAxis(1) =  cosTheta;
  yAxis(2) =  0;    
  
  zAxis(0) = 0.0;
  zAxis(1) = 0.0;
  zAxis(2) = 1.0;

  return 0;
}
