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

// $Revision: 1.4 $
// $Date: 2008-12-03 23:40:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CorotCrdTransfWarping2d.cpp,v $

// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 05/2000
// Revision: rms 06/2000 (using Assemble, and AssembleTranspose)
//
// Modified: 04/2005 Andreas Schellenberg (getBasicTrialVel, getBasicTrialAccel)
// 
// Purpose: This file contains the implementation for the 
// CorotCrdTransfWarping2d class. CorotCrdTransfWarping2d is a Corot
// transformation for a planar frame between the global 
// and basic coordinate systems.

/*
 * References
 *

General Formulation and Analytical Response Sensitivity
---
Scott, M. H. and F. C. Filippou (2007).
"Response Gradients for Nonlinear Beam-Column Elements under Large Displacements."
Journal of Structural Engineering, 133(2):155-165.

 *
 */

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>

#include <CorotCrdTransfWarping2d.h>

// initialize static variables
Matrix CorotCrdTransfWarping2d::Tlg(8,8);
Matrix CorotCrdTransfWarping2d::Tbl(5,8);
Vector CorotCrdTransfWarping2d::uxg(5); 
Vector CorotCrdTransfWarping2d::pg(8); 
Vector CorotCrdTransfWarping2d::dub(5); 
Vector CorotCrdTransfWarping2d::Dub(5); 
Matrix CorotCrdTransfWarping2d::kg(8,8);


// constructor:
CorotCrdTransfWarping2d::CorotCrdTransfWarping2d(int tag, 
                                   const Vector &rigJntOffsetI,
                                   const Vector &rigJntOffsetJ):
CrdTransf(tag, CRDTR_TAG_CorotCrdTransfWarping2d),
nodeIOffset(2), nodeJOffset(2), cosTheta(0), sinTheta(0),
cosAlpha(0), sinAlpha(0),
nodeIPtr(0), nodeJPtr(0), L(0), Ln(0), ub(5), ubcommit(5), ubpr(5),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{ 
    // check rigid joint offset for node I
    if (&rigJntOffsetI == 0 || rigJntOffsetI.Size() != 2 )
    {
        opserr << "CorotCrdTransfWarping2d::CorotCrdTransfWarping2d:  Invalid rigid joint offset vector for node I\n";
        opserr << "Size must be 2\n";      
        nodeIOffset.Zero();      
	}
    else
      nodeIOffset = rigJntOffsetI;
    
    // check rigid joint offset for node J
    if (&rigJntOffsetJ == 0 || rigJntOffsetJ.Size() != 2 )
    {
        opserr << "CorotCrdTransfWarping2d::CorotCrdTransfWarping2d:  Invalid rigid joint offset vector for node J\n";
        opserr << "Size must be 2\n";      
        nodeJOffset.Zero(); 
    }
    else
        nodeJOffset = rigJntOffsetJ;
    
    // temporary
    if (nodeIOffset.Norm() != 0 || nodeJOffset.Norm() != 0) {
      nodeOffsets = true;
    } else
      nodeOffsets = false;
}  


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
CorotCrdTransfWarping2d::CorotCrdTransfWarping2d():
CrdTransf(0, CRDTR_TAG_CorotCrdTransfWarping2d),
nodeIOffset(2), nodeJOffset(2), cosTheta(0), sinTheta(0),
cosAlpha(0), sinAlpha(0),
nodeIPtr(0), nodeJPtr(0), L(0), Ln(0), ub(5), ubcommit(5), ubpr(5),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    
}


// destructor:
CorotCrdTransfWarping2d::~CorotCrdTransfWarping2d() 
{
    if (nodeIInitialDisp != 0)
        delete [] nodeIInitialDisp;
    if (nodeJInitialDisp != 0)
        delete [] nodeJInitialDisp;
}


int
CorotCrdTransfWarping2d::commitState(void)
{
    ubcommit = ub;
    return 0;
}


int
CorotCrdTransfWarping2d::revertToLastCommit(void)
{
    ub = ubcommit;
    
    this->update();
    return 0;
}


int
CorotCrdTransfWarping2d::revertToStart(void)
{
    ub.Zero();
    this->update();
    return 0;
}


int 
CorotCrdTransfWarping2d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
    int error;
    
    nodeIPtr = nodeIPointer;
    nodeJPtr = nodeJPointer;
    
    if ((!nodeIPtr) || (!nodeJPtr))
    {
        opserr << "\nCorotCrdTransfWarping2d::initialize";
        opserr << "\ninvalid pointers to the element nodes\n";
        return -1;
    }
    
    // see if there is some initial displacements at nodes
    if (initialDispChecked == false) {
        const Vector &nodeIDisp = nodeIPtr->getDisp();
        const Vector &nodeJDisp = nodeJPtr->getDisp();
        for (int i=0; i<4; i++)
            if (nodeIDisp(i) != 0.0) {
                nodeIInitialDisp = new double [4];
                for (int j=0; j<4; j++)
                    nodeIInitialDisp[j] = nodeIDisp(j);
                i = 4;
            }
            
            for (int j=0; j<4; j++)
                if (nodeJDisp(j) != 0.0) {
                    nodeJInitialDisp = new double [4];
                    for (int i=0; i<4; i++)
                        nodeJInitialDisp[i] = nodeJDisp(i);
                    j = 8;
                }
                
                initialDispChecked = true;
    }
    
    // get element length and orientation
    if ((error = this->compElemtLengthAndOrient()))
        return error;
    
    return 0;
}


int  
CorotCrdTransfWarping2d::update(void)
{       
    // get global displacements 
    const Vector &dispI = nodeIPtr->getTrialDisp();
    const Vector &dispJ = nodeJPtr->getTrialDisp();
    
    static Vector ug(8);    
    for (int i = 0; i < 4; i++) {
        ug(i  ) = dispI(i);
        ug(i+4) = dispJ(i);
    }
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<4; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<4; j++)
            ug[j+4] -= nodeJInitialDisp[j];
    }

    // account for rigid offsets
    if (nodeOffsets == true) {
      ug(0) = ug(0) - ug(2) * nodeIOffset(1);
      ug(1) = ug(1) + ug(2) * nodeIOffset(0);
      
      ug(4) = ug(4) - ug(6) * nodeJOffset(1);
      ug(5) = ug(5) + ug(6) * nodeJOffset(0);
    }
    
    // transform global end displacements to local coordinates
    static Vector ul(8);
    
    ul(0) = cosTheta*ug(0) + sinTheta*ug(1);
    ul(1) = cosTheta*ug(1) - sinTheta*ug(0);
    ul(2) = ug(2);
	ul(3) = ug(3);
    ul(4) = cosTheta*ug(4) + sinTheta*ug(5);
    ul(5) = cosTheta*ug(5) - sinTheta*ug(4);
    ul(6) = ug(6);
    ul(7) = ug(7);

    // get deformed element length and orientation with respect to the local system
    this->compElemtLengthAndOrientWRTLocalSystem(ul);
    
    // determine displacements in the basic system eliminating rigid body modes 
    ubpr = ub;
    this->transfLocalDisplsToBasic(ul);
    
    // compute the transformation matrix from local to the basic system
    this->getTransfMatrixBasicLocal(Tbl);
    
    return 0;
}


int 
CorotCrdTransfWarping2d::compElemtLengthAndOrient(void)
{
    // element projection
    static Vector dx(2);
    
    if (nodeOffsets == true) 
      dx = (nodeJPtr->getCrds() + nodeJOffset) - (nodeIPtr->getCrds() + nodeIOffset);  
    else
      dx = nodeJPtr->getCrds() - nodeIPtr->getCrds();  
    
    if (nodeIInitialDisp != 0) {
        dx(0) -= nodeIInitialDisp[0];
        dx(1) -= nodeIInitialDisp[1];
    }
    
    if (nodeJInitialDisp != 0) {
        dx(0) += nodeJInitialDisp[0];
        dx(1) += nodeJInitialDisp[1];
    }


    // calculate the element length
    L = dx.Norm();
    
    if (L == 0.0) 
    {
        opserr << "\nCorotCrdTransfWarping2d::compElemtLengthAndOrien: 0 length\n";
        return -2;  
    }
    
    // calculate the element local x axis components (direction cosines)
    // wrt to the global coordinates 
    cosTheta = dx(0)/L;
    sinTheta = dx(1)/L;
   return 0;
}


int
CorotCrdTransfWarping2d::compElemtLengthAndOrientWRTLocalSystem(const Vector &ul)
{
    // get length and chord orientation of the deformed element with respect 
    // to the local coordinate system 
    // (deformed chord corresponding to the basic system)   
    
    double dulx, duly;
    
    dulx = ul(4) - ul(0);           // horizontal relative displacement  
    duly = ul(5) - ul(1);           // vertical relative displacement
    
    Lx = L + dulx;                  // horizontal component of the deformed member
    Ly = duly;                      // vertical component of the deformed member
    
    Ln = sqrt (Lx*Lx + Ly*Ly);      // chord length of the deformed configuration
    
    if (Ln == 0.0) {
        opserr << "\nCorotCrdTransfWarping2d::compElemtLengthAndOrientWRTLocalSystem: 0 length\n";
        return -2;  
    }
    
    cosAlpha = Lx / Ln;             // horizontal projection of the displaced chord 
    
    sinAlpha = Ly / Ln;             // vertical projection of the displaced chord
    
    return 0;
} 


void CorotCrdTransfWarping2d::transfLocalDisplsToBasic(const Vector &ul)
{
    // Eliminate rigid body modes, determining displacements wrt the basic system
    double alpha;
    
    alpha = atan2 (sinAlpha, cosAlpha);
    
    ub(0) = Ln - L;
    ub(1) = ul(2) - alpha;
	ub(2) = ul(3);
    ub(3) = ul(6) - alpha;
	ub(4) = ul(7); 	 
}


void
CorotCrdTransfWarping2d::getTransfMatrixLocalGlobal (Matrix &Tlg) 
{
    // setup transformation matrix
    Tlg.Zero();
    
	Tlg(0,0) = cosTheta;   
	Tlg(0,1) = sinTheta;  
	Tlg(1,0) =-sinTheta;    
	Tlg(1,1) = cosTheta;  
	Tlg(2,2) = 1.0; 
	Tlg(3,3) = 1.0;
	Tlg(4,4) = cosTheta;   
	Tlg(5,4) =-sinTheta;	
    Tlg(4,5) = sinTheta;	
	Tlg(5,5) = cosTheta;	
    Tlg(6,6) = 1.0; 
    Tlg(7,7) = 1.0; 	
}

void
CorotCrdTransfWarping2d::getTransfMatrixBasicLocal(Matrix &Tbl)
{
    // set up exact force transformation matrix from basic to local coordinates
	
    Tbl(0,0) = -cosAlpha;      
    Tbl(1,0) = -sinAlpha/Ln; 
    Tbl(2,0) = 0.0;
	Tbl(3,0) = -sinAlpha/Ln;
	Tbl(4,0) = 0.0;
    
    Tbl(0,1) = -sinAlpha;
    Tbl(1,1) =  cosAlpha/Ln;
    Tbl(2,1) =  0.0;
	Tbl(3,1) =  cosAlpha/Ln;
	Tbl(4,1) =  0.0;
    
    Tbl(0,2) =  0;
    Tbl(1,2) =  1;  
    Tbl(2,2) =  0;  
	Tbl(3,2) =  0;
	Tbl(4,2) =  0;

	Tbl(0,3) =  0;
    Tbl(1,3) =  0;  
    Tbl(2,3) =  1;  
	Tbl(3,3) =  0;
	Tbl(4,3) =  0;
    
    Tbl(0,4) =  cosAlpha;
    Tbl(1,4) =  sinAlpha/Ln;
    Tbl(2,4) =  0.0;
	Tbl(3,4) =  sinAlpha/Ln;
    Tbl(4,4) =  0.0;
    
    Tbl(0,5) =  sinAlpha;
    Tbl(1,5) = -cosAlpha/Ln;
    Tbl(2,5) = 0.0;
	Tbl(3,5) = -cosAlpha/Ln;
    Tbl(4,5) = 0.0;
    
    Tbl(0,6) =  0;
    Tbl(1,6) =  0; 
    Tbl(2,6) =  0; 
	Tbl(3,6) =  1;
    Tbl(4,6) =  0; 

    Tbl(0,7) =  0;
    Tbl(1,7) =  0; 
    Tbl(2,7) =  0; 
	Tbl(3,7) =  0;
    Tbl(4,7) =  1; 
}


const Vector &
CorotCrdTransfWarping2d::getBasicTrialDisp (void)
{
    return ub;    
}


const Vector &
CorotCrdTransfWarping2d::getBasicIncrDeltaDisp (void)
{
    // dub = ub - ubpr;
    dub = ub;
    dub.addVector (1.0, ubpr, -1.0);
    
    return dub;        
}


const Vector &
CorotCrdTransfWarping2d::getBasicIncrDisp(void)
{
    // Dub = ub - ubcommit;
    Dub = ub;
    Dub.addVector(1.0, ubcommit, -1.0);
    
    return Dub;        
}


const Vector &
CorotCrdTransfWarping2d::getBasicTrialVel(void)
{
	// determine global velocities
	const Vector &vel1 = nodeIPtr->getTrialVel();
	const Vector &vel2 = nodeJPtr->getTrialVel();
	
	static double vg[8];
	for (int i = 0; i < 4; i++) {
		vg[i]   = vel1(i);
		vg[i+4] = vel2(i);
	}
	
    // transform global end velocities to local coordinates
    static Vector vl(8);

    vl(0) = cosTheta*vg[0] + sinTheta*vg[1];
    vl(1) = cosTheta*vg[1] - sinTheta*vg[0];
    vl(2) = vg[2];
	vl(3) = vg[3];
    vl(4) = cosTheta*vg[4] + sinTheta*vg[5];
    vl(5) = cosTheta*vg[5] - sinTheta*vg[4];
    vl(6) = vg[6];
	vl(7) = vg[7];

    Lxdot = vl(4) - vl(0);
    Lydot = vl(5) - vl(1);

    // transform local velocities to basic coordinates
    static Vector vb(5);
	
    vb(0) = (Lx*Lxdot + Ly*Lydot)/Ln;
    vb(1) = vl(2) - (Lx*Lydot - Ly*Lxdot)/Ln/Ln;
    vb(2) = vl(3);
	vb(3) = vl(6) - (Lx*Lydot - Ly*Lxdot)/Ln/Ln;
    vb(4) = vl(7);
	
	return vb;
}


const Vector &
CorotCrdTransfWarping2d::getBasicTrialAccel(void)
{
	// determine global velocities
	const Vector &vel1 = nodeIPtr->getTrialVel();
	const Vector &vel2 = nodeJPtr->getTrialVel();
	
	static double vg[8];
	int i;
	for (i = 0; i < 4; i++) {
		vg[i]   = vel1(i);
		vg[i+4] = vel2(i);
	}
	
    // transform global end velocities to local coordinates
    static Vector vl(8);
    vl(0) = cosTheta*vg[0] + sinTheta*vg[1];
    vl(1) = cosTheta*vg[1] - sinTheta*vg[0];
    vl(2) = vg[2];
	vl(3) = vg[3];
    vl(4) = cosTheta*vg[4] + sinTheta*vg[5];
    vl(5) = cosTheta*vg[5] - sinTheta*vg[4];
    vl(6) = vg[6];
	vl(7) = vg[7];

    Lxdot = vl(4) - vl(0);
    Lydot = vl(5) - vl(1);
    
    // determine global accelerations
	const Vector &accel1 = nodeIPtr->getTrialAccel();
	const Vector &accel2 = nodeJPtr->getTrialAccel();
	
	static double ag[8];
	for (i = 0; i < 4; i++) {
		ag[i]   = accel1(i);
		ag[i+4] = accel2(i);
	}
	
    // transform global end accelerations to local coordinates
    static Vector al(8);

	al(0) = cosTheta*ag[0] + sinTheta*ag[1];
    al(1) = cosTheta*ag[1] - sinTheta*ag[0];
    al(2) = ag[2];
	al(3) = ag[3];
    al(4) = cosTheta*ag[4] + sinTheta*ag[5];
    al(5) = cosTheta*ag[5] - sinTheta*ag[4];
    al(6) = ag[6];
	al(7) = ag[7];

    Lxdotdot = al(4) - al(0);
    Lydotdot = al(5) - al(1);

    // transform local accelerations to basic coordinates
    static Vector ab(5);
	
    ab(0) = (Lxdot*Lxdot + Lx*Lxdotdot + Ly*Lydotdot + Lydot*Lydot)/Ln
          - pow(Lx*Lxdot + Ly*Lydot,2)/pow(Ln,3);
    ab(1) = al(2) - (Lx*Lydotdot - Ly*Lxdotdot)/pow(Ln,2)
          + 2*(Lx*Lydot - Ly*Lxdot)*(Lx*Lxdot + Ly*Lydot)/pow(Ln,4);
	ab(2) = al(3);
    ab(3) = al(6) + ab(1) - al(2);
	ab(4) = al(7);
	
	return ab;
}


const Vector &
CorotCrdTransfWarping2d::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
    
    // transform resisting forces from the basic system to local coordinates
    this->getTransfMatrixBasicLocal(Tbl);
    static Vector pl(8);
    pl.addMatrixTransposeVector(0.0, Tbl, pb, 1.0);    // pl = Tbl ^ pb;
    
    // add end forces due to element p0 loads
    // This assumes member loads are in local system
    pl(0) += p0(0);
    pl(1) += p0(1);
    pl(5) += p0(2);

    /*     // This assumes member loads are in basic system
    pl(0) += p0(0)*cosAlpha - p0(1)*sinAlpha;
    pl(1) += p0(0)*sinAlpha + p0(1)*cosAlpha;
    pl(3) -= p0(2)*sinAlpha;
    pl(4) += p0(2)*cosAlpha;
    */

    // transform resisting forces  from local to global coordinates
    //this->getTransfMatrixLocalGlobal(Tlg);     // OPTIMIZE LATER
    //pg.addMatrixTransposeVector(0.0, Tlg, pl, 1.0);   // pg = Tlg ^ pl; residual

    pg(0) = cosTheta*pl[0] - sinTheta*pl[1];
    pg(1) = sinTheta*pl[0] + cosTheta*pl[1];
    
    pg(4) = cosTheta*pl[4] - sinTheta*pl[5];
    pg(5) = sinTheta*pl[4] + cosTheta*pl[5];
    
    pg(2) = pl[2];
    pg(6) = pl[6];

	pg(3) = pl[3];
    pg(7) = pl[7];

    // account for rigid offsets
    if (nodeOffsets == true) {
      pg(2) += -pg(0) * nodeIOffset(1) + pg(1) * nodeIOffset(0);
      pg(6) += -pg(4) * nodeJOffset(1) + pg(5) * nodeJOffset(0);
    }
    

	return pg;
}


const Matrix &
CorotCrdTransfWarping2d::getGlobalStiffMatrix (const Matrix &kb, const Vector &pb)
{
    // transform tangent stiffness matrix from the basic system to local coordinates
    static Matrix kl(8,8);
    this->getTransfMatrixBasicLocal(Tbl);
    kl.addMatrixTripleProduct(0.0, Tbl, kb, 1.0);      // kl = Tbl ^ kb * Tbl;
    
    // add geometric stiffness matrix
    kl.addMatrix(1.0, this->getGeomStiffMatrix(pb), 1.0);
    
    // transform tangent  stiffness matrix from local to global coordinates
   this->getTransfMatrixLocalGlobal(Tlg);
   kg.addMatrixTripleProduct(0.0, Tlg, kl, 1.0); 
    return kg;
}


const Matrix &
CorotCrdTransfWarping2d::getInitialGlobalStiffMatrix (const Matrix &kb)
{
    // transform tangent stiffness matrix from the basic system to local coordinates
    static Matrix kl(8,8);
    static Matrix T(5,8);
    
	int nn = 3;

    T(0,0) = -1.0;
    T(1,0) = 0;
    T(2,0) = 0;
	T(3,0) = 0;
    T(4,0) = 0;
    
    T(0,1) =  0;
    T(1,1) =  1/L;
    T(2,1) =  0.0;
	T(3,1) =  1/L;
    T(4,1) =  0.0;
    
    T(0,2) =  0;
    T(1,2) =  1;  
    T(2,2) =  0;  
	T(3,2) =  0;  
    T(4,2) =  0; 

    T(0,3) =  0;
    T(1,3) =  0;  
    T(2,3) =  1;  
	T(3,3) =  0;  
    T(4,3) =  0; 
    
    T(0,4) =  1;
    T(1,4) =  0;
    T(2,4) =  0;
	T(3,4) =  0;
    T(4,4) =  0;
    
    T(0,5) =  0;
    T(1,5) = -1/L;
    T(2,5) = 0.0;
	T(3,5) = -1/L;
    T(4,5) = 0.0;
    
    T(0,6) =  0;
    T(1,6) =  0; 
    T(2,6) =  0;  
	T(3,6) =  1; 
    T(4,6) =  0;  

	T(0,7) =  0;
    T(1,7) =  0; 
    T(2,7) =  0;  
	T(3,7) =  0; 
    T(4,7) =  1;
    
    kl.addMatrixTripleProduct(0.0, T, kb, 1.0);    // kl = Tbl ^ kb * Tbl;
    
    // add geometric stiffness matrix
   //kl.addMatrix(1.0, this->getGeomStiffMatrix(pb), 1.0);
    
    // transform tangent  stiffness matrix from local to global coordinates
    this->getTransfMatrixLocalGlobal(Tlg);
    kg.addMatrixTripleProduct(0.0, Tlg, kl, 1.0); 
    return kg;
    }

    



const Matrix &
CorotCrdTransfWarping2d::getGeomStiffMatrix(const Vector &pb) const
{
    // get  geometric stiffness matrix present in the transformation 
    // from basic to local system
    double s2, c2, cs;  
    
    s2 = sinAlpha*sinAlpha;
    c2 = cosAlpha*cosAlpha;
    cs = sinAlpha*cosAlpha;
    
    static Matrix kg0(8,8), kg12(8,8);
    kg0.Zero();
    
    kg12.Zero();
    
    kg0(0,0) = kg0(4,4) =  s2;
    kg0(0,1) = kg0(4,5) = -cs;
    kg0(1,0) = kg0(5,4) = -cs;
    kg0(1,1) = kg0(5,5) =  c2;
    
    kg0(0,4) = kg0(4,0) = -s2;
    kg0(0,5) = kg0(4,1) =  cs;
    kg0(1,4) = kg0(5,0) =  cs;
    kg0(1,5) = kg0(5,1) = -c2;
    
    kg0 *= pb(0)/Ln;
    
    kg12(0,0) = kg12(4,4) = -2*cs;
    kg12(0,1) = kg12(4,5) =  c2-s2;
    kg12(1,0) = kg12(5,4) =  c2-s2;
    kg12(1,1) = kg12(5,5) =  2*cs;
    
    kg12(0,4) = kg12(4,0) =  2*cs;
    kg12(0,5) = kg12(4,1) = -c2+s2;
    kg12(1,4) = kg12(5,0) = -c2+s2;
    kg12(1,5) = kg12(5,1) = -2*cs;
    
    kg12 *= (pb(1)+pb(3))/(Ln*Ln);
    
    static Matrix kg(8,8);
    // kg = kg0 + kg12;
    kg = kg0;
    kg.addMatrix(1.0, kg12, 1.0); 
    
    return kg;
}


double 
CorotCrdTransfWarping2d::getInitialLength(void)
{
    return L;
}


double 
CorotCrdTransfWarping2d::getDeformedLength(void)
{
    return Ln;
}


CrdTransf *
CorotCrdTransfWarping2d::getCopy2d(void)
{
    // create a new instance of CorotCrdTransfWarping2d 
    CorotCrdTransfWarping2d *theCopy = new CorotCrdTransfWarping2d (this->getTag(), nodeIOffset, nodeJOffset);
    
    if (!theCopy)
    {
        opserr << "CorotCrdTransfWarping2d::getCopy() - out of memory creating copy\n";
	return 0;
    }    
    
    theCopy->nodeIPtr = nodeIPtr;
    theCopy->nodeJPtr = nodeJPtr;
    theCopy->cosTheta = cosTheta;
    theCopy->sinTheta = sinTheta;
    theCopy->cosAlpha = cosAlpha;
    theCopy->sinAlpha = sinAlpha;
    theCopy->L = L;
    theCopy->Ln = Ln;
    theCopy->ub = ub;
    theCopy->ubcommit = ubcommit;
    
    return theCopy;
}


int 
CorotCrdTransfWarping2d::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(15);
    data(0) = ubcommit(0);
    data(1) = ubcommit(1);
    data(2) = ubcommit(2);
	data(3) = ubcommit(3);
	data(4) = ubcommit(4);
    data(5) = nodeIOffset(0);
    data(6) = nodeIOffset(1);
    data(7) = nodeJOffset(0);
    data(8) = nodeJOffset(1);
    
    if (nodeIInitialDisp != 0) {
        data(9) = nodeIInitialDisp[0];
        data(10) = nodeIInitialDisp[1];
        data(11) = nodeIInitialDisp[2];
    } else {
        data(9) = 0.0;
        data(10) = 0.0;
        data(11) = 0.0;
    }
    
    if (nodeJInitialDisp != 0) {
        data(12) = nodeJInitialDisp[0];
        data(13) = nodeJInitialDisp[1];
        data(14) = nodeJInitialDisp[2];
    } else {
        data(12) = 0.0;
        data(13) = 0.0;
        data(14) = 0.0;
    }
    
    if (theChannel.sendVector(this->getTag(), cTag, data) < 0) {
        opserr << " CorotCrdTransfWarping2d::sendSelf() - data could not be sent\n" ;
        return -1;
    }
    return 0;
}


int 
CorotCrdTransfWarping2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(15);
    if (theChannel.recvVector(this->getTag(), cTag, data) < 0) {
        opserr << " CorotCrdTransfWarping2d::recvSelf() - data could not be received\n" ;
        return -1;
    }
    
    ubcommit(0) = data(0);
    ubcommit(1) = data(1);
    ubcommit(2) = data(2);
	ubcommit(3) = data(3);
	ubcommit(4) = data(4);

    nodeIOffset(0) =data(5);
    nodeIOffset(1) =data(6);
    nodeJOffset(0) =data(7);
    nodeJOffset(1) =data(8);
    
    int flag, i, j;
    flag = 0;
    for (i=9; i<=11; i++)
        if (data(i) != 0.0)
            flag = 1;
        if (flag == 1) {
            if (nodeIInitialDisp == 0)
                nodeIInitialDisp = new double[3];
            for (i=9, j=0; i<=11; i++, j++)
                nodeIInitialDisp[j] = data(i);
        }
        
        flag = 0;
        for (i=12; i<=14; i++)
            if (data(i) != 0.0)
                flag = 1;
            if (flag == 1) {
                if (nodeJInitialDisp == 0)
                    nodeJInitialDisp = new double [3];
                for (i=10, j=0; i<=13; i++, j++)
                    nodeJInitialDisp[j] = data(i);
            }
            
            ub = ubcommit;
            initialDispChecked = true;
            return 0;
}

const Matrix &
CorotCrdTransfWarping2d::getGlobalMatrixFromLocal(const Matrix &ml)
{
  return ml; // FIX LATER
}


const Vector &
CorotCrdTransfWarping2d::getPointGlobalCoordFromLocal(const Vector &xl)
{
    static Vector xg(5);
    opserr << " CorotCrdTransfWarping2d::getPointGlobalCoordFromLocal: not implemented yet" ;
    
    return xg;  
}


const Vector &
CorotCrdTransfWarping2d::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
    opserr << " CorotCrdTransfWarping2d::getPointGlobalDisplFromBasic: not implemented yet" ;
    
    return uxg;  
}


const Vector &
CorotCrdTransfWarping2d::getPointLocalDisplFromBasic (double xi, const Vector &uxb)
{
    opserr << " CorotCrdTransfWarping2d::getPointLocalDisplFromBasic: not implemented yet" ;
    
    return uxg;  
}


void
CorotCrdTransfWarping2d::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "\nCrdTransf: " << this->getTag() << " Type: CorotCrdTransfWarping2d";
		s << "\tnodeI Offset: " << nodeIOffset;
		s << "\tnodeJ Offset: " << nodeJOffset;
	}
	
	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"CorotCrdTransfWarping2d\"";
		if (nodeIOffset != 0)
			s << ", \"iOffset\": [" << nodeIOffset[0] << ", " << nodeIOffset[1] << "]";
		if (nodeJOffset != 0)
			s << ", \"jOffset\": [" << nodeJOffset[0] << ", " << nodeJOffset[1] << "]";
		s << "}";
	}
}

const Vector &
CorotCrdTransfWarping2d::getGlobalResistingForceShapeSensitivity(const Vector &q,
							  const Vector &p0,
							  int gradNumber)
{
  static Vector dpgdh(8);
  dpgdh.Zero();

  int nodeIid = nodeIPtr->getCrdsSensitivity();
  int nodeJid = nodeJPtr->getCrdsSensitivity();
  
  if (nodeIid == 0 && nodeJid == 0)
    return dpgdh;
    
  this->update();

  if (nodeIOffset.Norm() != 0.0 || nodeJOffset.Norm() != 0.0) {
    opserr << "ERROR: Currently a node offset cannot be used in " << endln
	   << " conjunction with random nodal coordinates." << endln;
  }
  
  double dcosThetadh = 0.0;
  double dsinThetadh = 0.0;

  double dx = cosTheta*L;
  double dy = sinTheta*L;	
  
  double dLdh = this->getdLdh();

  if (nodeIid == 1) { // here x1 is random
    //dcosThetadh = (-L+dx*dx/L)/(L*L);
    //dsinThetadh = dx*dy/(L*L*L);
    dcosThetadh = -1/L-cosTheta/L*dLdh;
    dsinThetadh = -sinTheta/L*dLdh;
  }
  if (nodeIid == 2) { // here y1 is random
    //dsinThetadh = (-L+dy*dy/L)/(L*L);
    //dcosThetadh = dx*dy/(L*L*L);
    dcosThetadh = -cosTheta/L*dLdh;
    dsinThetadh = -1/L-sinTheta/L*dLdh;
  }
  
  if (nodeJid == 1) { // here x2 is random
    //dcosThetadh = (L-dx*dx/L)/(L*L);
    //dsinThetadh = -dx*dy/(L*L*L);
    dcosThetadh = 1/L-cosTheta/L*dLdh;
    dsinThetadh = -sinTheta/L*dLdh;
  }
  if (nodeJid == 2) { // here y2 is random
    //dsinThetadh = (L-dy*dy/L)/(L*L);
    //dcosThetadh = -dx*dy/(L*L*L);
    dcosThetadh = -cosTheta/L*dLdh;
    dsinThetadh = 1/L-sinTheta/L*dLdh;
  }
  
  const Vector &disp1 = nodeIPtr->getTrialDisp();
  const Vector &disp2 = nodeJPtr->getTrialDisp();

  static Vector U(6);
  for (int i = 0; i < 4; i++) {
    U(i)   = disp1(i);
    U(i+4) = disp2(i);
  }
  
  static Vector u(8);

  double dux =  cosTheta*(U(4)-U(0)) + sinTheta*(U(5)-U(1));
  double duy = -sinTheta*(U(4)-U(0)) + cosTheta*(U(5)-U(1));

  //double dLdh = this->getdLdh();

  double dcosAlphadh =  sinAlpha*sinAlpha/Ln;
  double dsinAlphadh = -cosAlpha*sinAlpha/Ln;

  double dcosAlphaOverLndh = (2*sinAlpha*sinAlpha-1.0)/(Ln*Ln) ;
  double dsinAlphaOverLndh = -2*cosAlpha*sinAlpha/(Ln*Ln);

  double q0 = q(0);
  double q1 = q(1);
  double q2 = q(2);
  double q3 = q(3);
  double q4 = q(4);

  static Vector dpldh(8);
  dpldh.Zero();

  dpldh(0) = (-dcosAlphadh*q0 - dsinAlphaOverLndh*(q1+q2+q3+q4) )*dLdh;
  dpldh(1) = (-dsinAlphadh*q0 + dcosAlphaOverLndh*(q1+q2+q3+q4) )*dLdh;
  dpldh(2) = 0.0;
  dpldh(3) = 0.0;
  dpldh(4) = ( dcosAlphadh*q0 + dsinAlphaOverLndh*(q1+q2+q3+q4) )*dLdh;
  dpldh(5) = ( dsinAlphadh*q0 - dcosAlphaOverLndh*(q1+q2+q3+q4) )*dLdh;
  dpldh(6) = 0.0;
  dpldh(7) = 0.0;

  this->getTransfMatrixLocalGlobal(Tlg);     // OPTIMIZE LATER
  dpgdh.addMatrixTransposeVector(0.0, Tlg, dpldh, 1.0);   // pg = Tlg ^ pl; residual

  static Vector pl(8);
  pl.Zero();

  static Matrix Abl(5,8);
  this->getTransfMatrixBasicLocal(Abl);

  pl.addMatrixTransposeVector(0.0, Abl, q, 1.0); // OPTIMIZE LATER

  dpgdh(0) += dcosThetadh*pl(0)-dsinThetadh*pl(1);
  dpgdh(1) += dsinThetadh*pl(0)+dcosThetadh*pl(1);
  dpgdh(2) += 0.0;
  dpgdh(3) += 0.0;
  dpgdh(4) += dcosThetadh*pl(4)-dsinThetadh*pl(5);
  dpgdh(5) += dsinThetadh*pl(4)+dcosThetadh*pl(5);
  dpgdh(6) += 0.0;
  dpgdh(7) += 0.0;

  return dpgdh;
}

const Vector&
CorotCrdTransfWarping2d::getBasicDisplSensitivity(int gradNumber)
{
  static Vector dvdh(5);
  dvdh.Zero();

  int nodeIid = nodeIPtr->getCrdsSensitivity();
  int nodeJid = nodeJPtr->getCrdsSensitivity();
  
  this->update();

  double dcosThetadh = 0.0;
  double dsinThetadh = 0.0;

  double dx = cosTheta*L;
  double dy = sinTheta*L;	
  
  double dLdh = this->getdLdh();

  if (nodeIid == 1) { // here x1 is random
    //dcosThetadh = (-L+dx*dx/L)/(L*L);
    //dsinThetadh = dx*dy/(L*L*L);
    dcosThetadh = -1/L-cosTheta/L*dLdh;
    dsinThetadh = -sinTheta/L*dLdh;
  }
  if (nodeIid == 2) { // here y1 is random
    //dsinThetadh = (-L+dy*dy/L)/(L*L);
    //dcosThetadh = dx*dy/(L*L*L);
    dcosThetadh = -cosTheta/L*dLdh;
    dsinThetadh = -1/L-sinTheta/L*dLdh;
  }
  
  if (nodeJid == 1) { // here x2 is random
    //dcosThetadh = (L-dx*dx/L)/(L*L);
    //dsinThetadh = -dx*dy/(L*L*L);
    dcosThetadh = 1/L-cosTheta/L*dLdh;
    dsinThetadh = -sinTheta/L*dLdh;
  }
  if (nodeJid == 2) { // here y2 is random
    //dsinThetadh = (L-dy*dy/L)/(L*L);
    //dcosThetadh = -dx*dy/(L*L*L);
    dcosThetadh = -cosTheta/L*dLdh;
    dsinThetadh = 1/L-sinTheta/L*dLdh;
  }
  
  static Vector U(8);
  static Vector dUdh(8);

  const Vector &disp1 = nodeIPtr->getTrialDisp();
  const Vector &disp2 = nodeJPtr->getTrialDisp();
  for (int i = 0; i < 4; i++) {
    U(i)   = disp1(i);
    U(i+4) = disp2(i);
    dUdh(i)   = nodeIPtr->getDispSensitivity((i+1),gradNumber);
    dUdh(i+4) = nodeJPtr->getDispSensitivity((i+1),gradNumber);
  }

  static Vector dudh(8);

  dudh(0) =  cosTheta*dUdh(0) + sinTheta*dUdh(1);
  dudh(1) = -sinTheta*dUdh(0) + cosTheta*dUdh(1);
  dudh(2) =  dUdh(2);
  dudh(3) =  dUdh(3);
  dudh(4) =  cosTheta*dUdh(4) + sinTheta*dUdh(5);
  dudh(5) = -sinTheta*dUdh(4) + cosTheta*dUdh(5);
  dudh(6) =  dUdh(6);
  dudh(7) =  dUdh(7);

  if (nodeIid != 0 || nodeJid != 0) {
    dudh(0) +=  dcosThetadh*U(0) + dsinThetadh*U(1);
    dudh(1) += -dsinThetadh*U(0) + dcosThetadh*U(1);
    dudh(3) +=  dcosThetadh*U(4) + dsinThetadh*U(5);
    dudh(4) += -dsinThetadh*U(4) + dcosThetadh*U(5);
  }

  double duxdh = dudh(4)-dudh(0);
  double duydh = dudh(5)-dudh(1);

  //double dLdh  = this->getdLdh();
  double dLndh = cosAlpha*(dLdh+duxdh) + sinAlpha*duydh;

  double dalphadh = (cosAlpha*duydh - sinAlpha*(dLdh+duxdh))/Ln;

  // direct differentiation of v(u) wrt theta
  dvdh(0) = dLndh-dLdh;
  dvdh(1) = dudh(2)-dalphadh;
  dvdh(2) = dudh(5)-dalphadh;

  return dvdh;
}

const Vector&
CorotCrdTransfWarping2d::getBasicTrialDispShapeSensitivity(void)
{
  static Vector dvdh(5);
  dvdh.Zero();

  int nodeIid = nodeIPtr->getCrdsSensitivity();
  int nodeJid = nodeJPtr->getCrdsSensitivity();
  
  if (nodeIid == 0 && nodeJid == 0)
    return dvdh;

  static Matrix Abl(5,8);

  this->update();
  this->getTransfMatrixBasicLocal(Abl);

  double dcosThetadh = 0.0;
  double dsinThetadh = 0.0;

  double dx = cosTheta*L;
  double dy = sinTheta*L;	
  
  double dLdh = this->getdLdh();

  if (nodeIid == 1) { // here x1 is random
    //dcosThetadh = (-L+dx*dx/L)/(L*L);
    //dsinThetadh = dx*dy/(L*L*L);
    dcosThetadh = -1/L-cosTheta/L*dLdh;
    dsinThetadh = -sinTheta/L*dLdh;
  }
  if (nodeIid == 2) { // here y1 is random
    //dsinThetadh = (-L+dy*dy/L)/(L*L);
    //dcosThetadh = dx*dy/(L*L*L);
    dcosThetadh = -cosTheta/L*dLdh;
    dsinThetadh = -1/L-sinTheta/L*dLdh;
  }
  
  if (nodeJid == 1) { // here x2 is random
    //dcosThetadh = (L-dx*dx/L)/(L*L);
    //dsinThetadh = -dx*dy/(L*L*L);
    dcosThetadh = 1/L-cosTheta/L*dLdh;
    dsinThetadh = -sinTheta/L*dLdh;
  }
  if (nodeJid == 2) { // here y2 is random
    //dsinThetadh = (L-dy*dy/L)/(L*L);
    //dcosThetadh = -dx*dy/(L*L*L);
    dcosThetadh = -cosTheta/L*dLdh;
    dsinThetadh = 1/L-sinTheta/L*dLdh;
  }
  
  const Vector &disp1 = nodeIPtr->getTrialDisp();
  const Vector &disp2 = nodeJPtr->getTrialDisp();

  static Vector U(8);
  for (int i = 0; i < 4; i++) {
    U(i)   = disp1(i);
    U(i+4) = disp2(i);
  }

  dvdh(0) = (cosAlpha-1.0)*dLdh;
  dvdh(1) =  (sinAlpha/Ln)*dLdh;
  dvdh(2) =  (sinAlpha/Ln)*dLdh;

  static Vector dAdh_U(8);
  // dAdh * U
  dAdh_U(0) =  dcosThetadh*U(0) + dsinThetadh*U(1);
  dAdh_U(1) = -dsinThetadh*U(0) + dcosThetadh*U(1);
  dAdh_U(2) = 0.0;
  dAdh_U(3) = 0.0;
  dAdh_U(4) =  dcosThetadh*U(4) + dsinThetadh*U(5);
  dAdh_U(5) = -dsinThetadh*U(4) + dcosThetadh*U(5);
  dAdh_U(6) = 0.0;
  dAdh_U(7) = 0.0;

  dvdh += Abl*dAdh_U;

  return dvdh;
}

bool
CorotCrdTransfWarping2d::isShapeSensitivity(void)
{
  int nodeIid = nodeIPtr->getCrdsSensitivity();
  int nodeJid = nodeJPtr->getCrdsSensitivity();
  
  return (nodeIid != 0 || nodeJid != 0);
}

double
CorotCrdTransfWarping2d::getdLdh(void)
{
  int nodeIid = nodeIPtr->getCrdsSensitivity();
  int nodeJid = nodeJPtr->getCrdsSensitivity();
  
  if (nodeIid == 0 && nodeJid == 0) 
    return 0.0;
    
  if (nodeIOffset.Norm() != 0.0 || nodeJOffset.Norm() != 0.0) {
    opserr << "ERROR: Currently a node offset cannot be used in " << endln
	   << " conjunction with random nodal coordinates." << endln;
  }
  
  if (nodeIid == 1) // here x1 is random
    return -cosTheta;
  if (nodeIid == 2) // here y1 is random
    return -sinTheta;
  
  if (nodeJid == 1) // here x2 is random
    return cosTheta;
  if (nodeJid == 2) // here y2 is random
    return sinTheta;

  return 0.0;
}

double
CorotCrdTransfWarping2d::getd1overLdh(void)
{
  int nodeIid = nodeIPtr->getCrdsSensitivity();
  int nodeJid = nodeJPtr->getCrdsSensitivity();
  
  if (nodeIid == 0 && nodeJid == 0)
    return 0.0;
    
  if (nodeIOffset.Norm() != 0.0 || nodeJOffset.Norm() != 0.0) {
    opserr << "ERROR: Currently a node offset cannot be used in " << endln
	   << " conjunction with random nodal coordinates." << endln;
  }
  
  if (nodeIid == 1) // here x1 is random
    return cosTheta/(L*L);
  if (nodeIid == 2) // here y1 is random
    return sinTheta/(L*L);
  
  if (nodeJid == 1) // here x2 is random
    return -cosTheta/(L*L);
  if (nodeJid == 2) // here y2 is random
    return -sinTheta/(L*L);

  return 0.0;
}

// AddingSensitivity:END /////////////////////////////////////
