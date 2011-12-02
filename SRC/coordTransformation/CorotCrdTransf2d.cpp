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
// $Date: 2007-08-07 16:52:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CorotCrdTransf2d.cpp,v $

// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 05/2000
// Revision: rms 06/2000 (using Assemble, and AssembleTranspose)
//
// Modified: 04/2005 Andreas Schellenberg (getBasicTrialVel, getBasicTrialAccel)
// 
// Purpose: This file contains the implementation for the 
// CorotCrdTransf2d class. CorotCrdTransf2d is a Corot
// transformation for a planar frame between the global 
// and basic coordinate systems.


#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>

#include <CorotCrdTransf2d.h>

// initialize static variables
Matrix CorotCrdTransf2d::Tlg(6,6);
Matrix CorotCrdTransf2d::Tbl(3,6);
Vector CorotCrdTransf2d::uxg(3); 
Vector CorotCrdTransf2d::pg(6); 
Vector CorotCrdTransf2d::dub(3); 
Vector CorotCrdTransf2d::Dub(3); 
Matrix CorotCrdTransf2d::kg(6,6);


// constructor:
CorotCrdTransf2d::CorotCrdTransf2d(int tag, 
                                   const Vector &rigJntOffsetI,
                                   const Vector &rigJntOffsetJ):
CrdTransf2d(tag, CRDTR_TAG_CorotCrdTransf2d),
nodeIOffset(2), nodeJOffset(2), cosTheta(0), sinTheta(0),
cosAlpha(0), sinAlpha(0),
nodeIPtr(0), nodeJPtr(0), L(0), Ln(0), ub(3), ubcommit(3), ubpr(3),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    // check rigid joint offset for node I
    if (&rigJntOffsetI == 0 || rigJntOffsetI.Size() != 2 )
    {
        opserr << "CorotCrdTransf2d::CorotCrdTransf2d:  Invalid rigid joint offset vector for node I\n";
        opserr << "Size must be 2\n";      
        nodeIOffset.Zero();      
    }
    else
      nodeIOffset = rigJntOffsetI;
    
    // check rigid joint offset for node J
    if (&rigJntOffsetJ == 0 || rigJntOffsetJ.Size() != 2 )
    {
        opserr << "CorotCrdTransf2d::CorotCrdTransf2d:  Invalid rigid joint offset vector for node J\n";
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
CorotCrdTransf2d::CorotCrdTransf2d():
CrdTransf2d(0, CRDTR_TAG_CorotCrdTransf2d),
nodeIOffset(2), nodeJOffset(2), cosTheta(0), sinTheta(0),
cosAlpha(0), sinAlpha(0),
nodeIPtr(0), nodeJPtr(0), L(0), Ln(0), ub(3), ubcommit(3), ubpr(3),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    
}


// destructor:
CorotCrdTransf2d::~CorotCrdTransf2d() 
{
    if (nodeIInitialDisp != 0)
        delete [] nodeIInitialDisp;
    if (nodeJInitialDisp != 0)
        delete [] nodeJInitialDisp;
}


int
CorotCrdTransf2d::commitState(void)
{
    ubcommit = ub;
    return 0;
}


int
CorotCrdTransf2d::revertToLastCommit(void)
{
    ub = ubcommit;
    
    this->update();
    return 0;
}


int
CorotCrdTransf2d::revertToStart(void)
{
    ub.Zero();
    this->update();
    return 0;
}


int 
CorotCrdTransf2d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
    int error;
    
    nodeIPtr = nodeIPointer;
    nodeJPtr = nodeJPointer;
    
    if ((!nodeIPtr) || (!nodeJPtr))
    {
        opserr << "\nCorotCrdTransf2d::initialize";
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
                    j = 6;
                }
                
                initialDispChecked = true;
    }
    
    // get element length and orientation
    if ((error = this->compElemtLengthAndOrient()))
        return error;
    
    return 0;
}


int  
CorotCrdTransf2d::update(void)
{       
    // get global displacements 
    const Vector &dispI = nodeIPtr->getTrialDisp();
    const Vector &dispJ = nodeJPtr->getTrialDisp();
    
    static Vector ug(6);    
    for (int i = 0; i < 3; i++) {
        ug(i  ) = dispI(i);
        ug(i+3) = dispJ(i);
    }
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<3; j++)
            ug[j] -= nodeIInitialDisp[j];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<3; j++)
            ug[j+3] -= nodeJInitialDisp[j];
    }

    // account for rigid offsets
    if (nodeOffsets == true) {
      ug(0) = ug(0) - ug(2) * nodeIOffset(1);
      ug(1) = ug(1) + ug(2) * nodeIOffset(0);
      
      ug(3) = ug(3) - ug(5) * nodeJOffset(1);
      ug(4) = ug(4) + ug(5) * nodeJOffset(0);
    }
    
    // transform global end displacements to local coordinates
    static Vector ul(6);
    
    ul(0) = cosTheta*ug(0) + sinTheta*ug(1);
    ul(1) = cosTheta*ug(1) - sinTheta*ug(0);
    ul(2) = ug(2);
    ul(3) = cosTheta*ug(3) + sinTheta*ug(4);
    ul(4) = cosTheta*ug(4) - sinTheta*ug(3);
    ul(5) = ug(5);
    
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
CorotCrdTransf2d::compElemtLengthAndOrient(void)
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
        opserr << "\nCorotCrdTransf2d::compElemtLengthAndOrien: 0 length\n";
        return -2;  
    }
    
    // calculate the element local x axis components (direction cosines)
    // wrt to the global coordinates 
    cosTheta = dx(0)/L;
    sinTheta = dx(1)/L;
    
    return 0;
}


int
CorotCrdTransf2d::compElemtLengthAndOrientWRTLocalSystem(const Vector &ul)
{
    // get length and chord orientation of the deformed element with respect 
    // to the local coordinate system 
    // (deformed chord corresponding to the basic system)   
    
    double dulx, duly;
    
    dulx = ul(3) - ul(0);           // horizontal relative displacement  
    duly = ul(4) - ul(1);           // vertical relative displacement
    
    Lx = L + dulx;                  // horizontal component of the deformed member
    Ly = duly;                      // vertical component of the deformed member
    
    Ln = sqrt (Lx*Lx + Ly*Ly);      // chord length of the deformed configuration
    
    if (Ln == 0.0) {
        opserr << "\nCorotCrdTransf2d::compElemtLengthAndOrientWRTLocalSystem: 0 length\n";
        return -2;  
    }
    
    cosAlpha = Lx / Ln;             // horizontal projection of the displaced chord 
    
    sinAlpha = Ly / Ln;             // vertical projection of the displaced chord
    
    return 0;
} 


void CorotCrdTransf2d::transfLocalDisplsToBasic(const Vector &ul)
{
    // Eliminate rigid body modes, determining displacements wrt the basic system
    double alpha;
    
    alpha = atan2 (sinAlpha, cosAlpha);
    
    ub(0) = Ln - L;
    ub(1) = ul(2) - alpha;
    ub(2) = ul(5) - alpha;
}


void
CorotCrdTransf2d::getTransfMatrixLocalGlobal (Matrix &Tlg) 
{
    // setup transformation matrix
    Tlg.Zero();
    
    Tlg(0,0) = Tlg(3,3) =  cosTheta;           
    Tlg(0,1) = Tlg(3,4) =  sinTheta;
    Tlg(1,0) = Tlg(4,3) = -sinTheta;
    Tlg(1,1) = Tlg(4,4) =  cosTheta;
    Tlg(2,2) = Tlg(5,5) =  1.0;
}


void
CorotCrdTransf2d::getTransfMatrixBasicLocal(Matrix &Tbl)
{
    // set up exact force transformation matrix from basic to local coordinates
    Tbl(0,0) = -cosAlpha;      
    Tbl(1,0) = -sinAlpha/Ln;
    Tbl(2,0) = -sinAlpha/Ln;
    
    Tbl(0,1) = -sinAlpha;
    Tbl(1,1) =  cosAlpha/Ln;
    Tbl(2,1) =  cosAlpha/Ln;
    
    Tbl(0,2) =  0;
    Tbl(1,2) =  1;  
    Tbl(2,2) =  0;  
    
    Tbl(0,3) =  cosAlpha;
    Tbl(1,3) =  sinAlpha/Ln;
    Tbl(2,3) =  sinAlpha/Ln;
    
    Tbl(0,4) =  sinAlpha;
    Tbl(1,4) = -cosAlpha/Ln;
    Tbl(2,4) = -cosAlpha/Ln;
    
    Tbl(0,5) =  0;
    Tbl(1,5) =  0; 
    Tbl(2,5) =  1;  
}


const Vector &
CorotCrdTransf2d::getBasicTrialDisp (void)
{
    return ub;    
}


const Vector &
CorotCrdTransf2d::getBasicIncrDeltaDisp (void)
{
    // dub = ub - ubpr;
    dub = ub;
    dub.addVector (1.0, ubpr, -1.0);
    
    return dub;        
}


const Vector &
CorotCrdTransf2d::getBasicIncrDisp(void)
{
    // Dub = ub - ubcommit;
    Dub = ub;
    Dub.addVector(1.0, ubcommit, -1.0);
    
    return Dub;        
}


const Vector &
CorotCrdTransf2d::getBasicTrialVel(void)
{
	// determine global velocities
	const Vector &vel1 = nodeIPtr->getTrialVel();
	const Vector &vel2 = nodeJPtr->getTrialVel();
	
	static double vg[6];
	for (int i = 0; i < 3; i++) {
		vg[i]   = vel1(i);
		vg[i+3] = vel2(i);
	}
	
    // transform global end velocities to local coordinates
    static Vector vl(6);

    vl(0) = cosTheta*vg[0] + sinTheta*vg[1];
    vl(1) = cosTheta*vg[1] - sinTheta*vg[0];
    vl(2) = vg[2];
    vl(3) = cosTheta*vg[3] + sinTheta*vg[4];
    vl(4) = cosTheta*vg[4] - sinTheta*vg[3];
    vl(5) = vg[5];

    Lxdot = vl(3) - vl(0);
    Lydot = vl(4) - vl(1);

    // transform local velocities to basic coordinates
    static Vector vb(3);
	
    vb(0) = (Lx*Lxdot + Ly*Lydot)/Ln;
    vb(1) = vl(2) - (Lx*Lydot - Ly*Lxdot)/pow(Ln,2);
    vb(2) = vl(5) + vb(1) - vl(2);
	
	return vb;
}


const Vector &
CorotCrdTransf2d::getBasicTrialAccel(void)
{
	// determine global velocities
	const Vector &vel1 = nodeIPtr->getTrialVel();
	const Vector &vel2 = nodeJPtr->getTrialVel();
	
	static double vg[6];
	int i;
	for (i = 0; i < 3; i++) {
		vg[i]   = vel1(i);
		vg[i+3] = vel2(i);
	}
	
    // transform global end velocities to local coordinates
    static Vector vl(6);

    vl(0) = cosTheta*vg[0] + sinTheta*vg[1];
    vl(1) = cosTheta*vg[1] - sinTheta*vg[0];
    vl(2) = vg[2];
    vl(3) = cosTheta*vg[3] + sinTheta*vg[4];
    vl(4) = cosTheta*vg[4] - sinTheta*vg[3];
    vl(5) = vg[5];

    Lxdot = vl(3) - vl(0);
    Lydot = vl(4) - vl(1);
    
    // determine global accelerations
	const Vector &accel1 = nodeIPtr->getTrialAccel();
	const Vector &accel2 = nodeJPtr->getTrialAccel();
	
	static double ag[6];
	for (i = 0; i < 3; i++) {
		ag[i]   = accel1(i);
		ag[i+3] = accel2(i);
	}
	
    // transform global end accelerations to local coordinates
    static Vector al(6);

    al(0) = cosTheta*ag[0] + sinTheta*ag[1];
    al(1) = cosTheta*ag[1] - sinTheta*ag[0];
    al(2) = ag[2];
    al(3) = cosTheta*ag[3] + sinTheta*ag[4];
    al(4) = cosTheta*ag[4] - sinTheta*ag[3];
    al(5) = ag[5];

    Lxdotdot = al(3) - al(0);
    Lydotdot = al(4) - al(1);

    // transform local accelerations to basic coordinates
    static Vector ab(3);
	
    ab(0) = (Lxdot*Lxdot + Lx*Lxdotdot + Ly*Lydotdot + Lydot*Lydot)/Ln
          - pow(Lx*Lxdot + Ly*Lydot,2)/pow(Ln,3);
    ab(1) = al(2) - (Lx*Lydotdot - Ly*Lxdotdot)/pow(Ln,2)
          + 2*(Lx*Lydot - Ly*Lxdot)*(Lx*Lxdot + Ly*Lydot)/pow(Ln,4);
    ab(2) = al(5) + ab(1) - al(2);
	
	return ab;
}


const Vector &
CorotCrdTransf2d::getGlobalResistingForce(const Vector &pb, const Vector &unifLoad)
{
    
    // transform resisting forces from the basic system to local coordinates
    this->getTransfMatrixBasicLocal(Tbl);
    static Vector pl(6);
    pl.addMatrixTransposeVector(0.0, Tbl, pb, 1.0);    // pl = Tbl ^ pb;
    
    ///opserr << "pl: " << pl;
    // check distributed load is zero (not implemented yet)
    
    if (unifLoad.Norm() != 0.0)
    {
        opserr << "CorotCrdTransf2d::getGlobalResistingForce: affect of Po not implemented yet.";
        opserr << "using zero value";
    }  
    
    // transform resisting forces  from local to global coordinates
    this->getTransfMatrixLocalGlobal(Tlg);     // OPTIMIZE LATER
    pg.addMatrixTransposeVector(0.0, Tlg, pl, 1.0);   // pg = Tlg ^ pl; residual


    // account for rigid offsets
    if (nodeOffsets == true) {
      pg(2) += -pg(0) * nodeIOffset(1) + pg(1) * nodeIOffset(0);
      pg(5) += -pg(3) * nodeJOffset(1) + pg(4) * nodeJOffset(0);
    }
    


    return pg;
}


const Matrix &
CorotCrdTransf2d::getGlobalStiffMatrix (const Matrix &kb, const Vector &pb)
{
    // transform tangent stiffness matrix from the basic system to local coordinates
    static Matrix kl(6,6);
    this->getTransfMatrixBasicLocal(Tbl);
    kl.addMatrixTripleProduct(0.0, Tbl, kb, 1.0);      // kl = Tbl ^ kb * Tbl;
    
    // add geometric stiffness matrix
    kl.addMatrix(1.0, this->getGeomStiffMatrix(pb), 1.0);
    
    // transform tangent  stiffness matrix from local to global coordinates
    
    // kg.addMatrixTripleProduct(0.0, Tlg, kl, 1.0);
    double s2, c2, cs;  
    
    s2 = sinTheta*sinTheta;
    c2 = cosTheta*cosTheta;
    cs = sinTheta*cosTheta;
    
    double k11, k12, k13, k21, k22, k23, k31, k32, k33;
    
    k11 = kl(0,0);    k12 = kl(0,1);    k13 = kl(0,2);
    k21 = kl(1,0);    k22 = kl(1,1);    k23 = kl(1,2);
    k31 = kl(2,0);    k32 = kl(2,1);    k33 = kl(2,2);
    
    kg(0,0) = c2*k11+s2*k22-cs*(k21+k12); 
    kg(1,0) = c2*k21-s2*k12+cs*(k11-k22);
    kg(2,0) = cosTheta*k31-sinTheta*k32;
    
    kg(0,1) = c2*k12-s2*k21+cs*(k11-k22);
    kg(1,1) = c2*k22+s2*k11+cs*(k21+k12);
    kg(2,1) = sinTheta*k31+cosTheta*k32;
    
    kg(0,2) = cosTheta*k13 - sinTheta*k23;
    kg(1,2) = sinTheta*k13 + cosTheta*k23;
    kg(2,2) = k33;
    
    k11 = kl(0,3);    k12 = kl(0,4);    k13 = kl(0,5);
    k21 = kl(1,3);    k22 = kl(1,4);    k23 = kl(1,5);
    k31 = kl(2,3);    k32 = kl(2,4);    k33 = kl(2,5);
    
    kg(0,3) = c2*k11+s2*k22-cs*(k21+k12); 
    kg(1,3) = c2*k21-s2*k12+cs*(k11-k22);
    kg(2,3) = cosTheta*k31-sinTheta*k32;
    
    kg(0,4) = c2*k12-s2*k21+cs*(k11-k22);
    kg(1,4) = c2*k22+s2*k11+cs*(k21+k12);
    kg(2,4) = sinTheta*k31+cosTheta*k32;
    
    kg(0,5) = cosTheta*k13 - sinTheta*k23;
    kg(1,5) = sinTheta*k13 + cosTheta*k23;
    kg(2,5) = k33;
    
    k11 = kl(3,0);    k12 = kl(3,1);    k13 = kl(3,2);
    k21 = kl(4,0);    k22 = kl(4,1);    k23 = kl(4,2);
    k31 = kl(5,0);    k32 = kl(5,1);    k33 = kl(5,2);
    
    kg(3,0) = c2*k11+s2*k22-cs*(k21+k12); 
    kg(4,0) = c2*k21-s2*k12+cs*(k11-k22);
    kg(5,0) = cosTheta*k31-sinTheta*k32;
    
    kg(3,1) = c2*k12-s2*k21+cs*(k11-k22);
    kg(4,1) = c2*k22+s2*k11+cs*(k21+k12);
    kg(5,1) = sinTheta*k31+cosTheta*k32;
    
    kg(3,2) = cosTheta*k13 - sinTheta*k23;
    kg(4,2) = sinTheta*k13 + cosTheta*k23;
    kg(5,2) = k33;
    
    k11 = kl(3,3);    k12 = kl(3,4);    k13 = kl(3,5);
    k21 = kl(4,3);    k22 = kl(4,4);    k23 = kl(4,5);
    k31 = kl(5,3);    k32 = kl(5,4);    k33 = kl(5,5);
    
    kg(3,3) = c2*k11+s2*k22-cs*(k21+k12); 
    kg(4,3) = c2*k21-s2*k12+cs*(k11-k22);
    kg(5,3) = cosTheta*k31-sinTheta*k32;
    
    kg(3,4) = c2*k12-s2*k21+cs*(k11-k22);
    kg(4,4) = c2*k22+s2*k11+cs*(k21+k12);
    kg(5,4) = sinTheta*k31+cosTheta*k32;
    
    kg(3,5) = cosTheta*k13 - sinTheta*k23;
    kg(4,5) = sinTheta*k13 + cosTheta*k23;
    kg(5,5) = k33;


    if (nodeOffsets == true) {
      double X1 = nodeIOffset(0);
      double Y1 = nodeIOffset(1);
      double X2 = nodeJOffset(0);
      double Y2 = nodeJOffset(1);

      double k11 = kg(0,0);

      double k12 = kg(0,1);
      double k22 = kg(1,1);

      double k13 = kg(0,2);
      double k23 = kg(1,2);
      double k33 = kg(2,2);

      double k14 = kg(0,3);
      double k24 = kg(1,3);
      double k34 = kg(2,3);
      double k44 = kg(3,3);

      double k15 = kg(0,4);
      double k25 = kg(1,4);
      double k35 = kg(2,4);
      double k45 = kg(3,4);
      double k55 = kg(3,4);

      double k16 = kg(0,5);
      double k26 = kg(1,5);
      double k36 = kg(2,5);
      double k46 = kg(3,5);
      double k56 = kg(4,5);
      double k66 = kg(5,5);


      double K13 = -k11*Y1 + k12*X1 + k13;
      double K23 = -k12*Y1 + k22*X1 + k23;
      kg(0,2) = K13;
      kg(2,0) = K13;
      kg(1,2) = K23;
      kg(2,1) = K23;
      kg(2,2) = -Y1*K13 + X1*K23 -k13*Y1 + k23*X1 + k33;
      

      double K16 = -k14*Y2 + k15*X2 + k16;
      double K26 = -k24*Y2 + k25*X2 + k26;
      kg(0,5) = K16;
      kg(5,0) = K16;
      kg(1,5) = K26;
      kg(5,1) = K26;
      kg(2,5) = -Y2*K16 + X2*K26 -k16*Y1 + k26*X1 + k36;      
      kg(5,2) = kg(2,5);

      double K46 = -k44*Y2 + k45*X2 + k46;
      double K56 = -k45*Y2 + k55*X2 + k56;
      kg(3,5) = K46;
      kg(5,3) = K46;
      kg(4,5) = K56;
      kg(5,4) = K56;
      kg(5,5) = -Y2*K46 + X2*K56 -k46*Y2 + k56*X2 + k66;      

      double K34 = -k14*Y1 + k24*X1 + k34;
      double K35 = -k15*Y1 + k25*X1 + k35;
      kg(2,3) = K34;
      kg(3,2) = K34;
      kg(2,4) = K35;
      kg(4,2) = K35;

    }
    
    return kg;
}


const Matrix &
CorotCrdTransf2d::getInitialGlobalStiffMatrix (const Matrix &kb)
{
    // transform tangent stiffness matrix from the basic system to local coordinates
    static Matrix kl(6,6);
    static Matrix T(3,6);
    
    T(0,0) = -1.0;
    T(1,0) = 0;
    T(2,0) = 0;
    
    T(0,1) = 0;
    T(1,1) = 1/L;
    T(2,1) = 1/L;
    
    T(0,2) =  0;
    T(1,2) =  1;  
    T(2,2) =  0;  
    
    T(0,3) =  1;
    T(1,3) =  0;
    T(2,3) =  0;
    
    T(0,4) =  0;
    T(1,4) = -1/L;
    T(2,4) = -1/L;
    
    T(0,5) =  0;
    T(1,5) =  0; 
    T(2,5) =  1;  
    
    kl.addMatrixTripleProduct(0.0, T, kb, 1.0);      // kl = Tbl ^ kb * Tbl;
    
    // add geometric stiffness matrix
    // kl.addMatrix(1.0, this->getGeomStiffMatrix(pb), 1.0);
    
    // transform tangent  stiffness matrix from local to global coordinates
    
    // kg.addMatrixTripleProduct(0.0, Tlg, kl, 1.0);
    double s2, c2, cs;  
    
    s2 = sinTheta*sinTheta;
    c2 = cosTheta*cosTheta;
    cs = sinTheta*cosTheta;
    
    double k11, k12, k13, k21, k22, k23, k31, k32, k33;
    
    k11 = kl(0,0);    k12 = kl(0,1);    k13 = kl(0,2);
    k21 = kl(1,0);    k22 = kl(1,1);    k23 = kl(1,2);
    k31 = kl(2,0);    k32 = kl(2,1);    k33 = kl(2,2);
    
    kg(0,0) = c2*k11+s2*k22-cs*(k21+k12); 
    kg(1,0) = c2*k21-s2*k12+cs*(k11-k22);
    kg(2,0) = cosTheta*k31-sinTheta*k32;
    
    kg(0,1) = c2*k12-s2*k21+cs*(k11-k22);
    kg(1,1) = c2*k22+s2*k11+cs*(k21+k12);
    kg(2,1) = sinTheta*k31+cosTheta*k32;
    
    kg(0,2) = cosTheta*k13 - sinTheta*k23;
    kg(1,2) = sinTheta*k13 + cosTheta*k23;
    kg(2,2) = k33;
    
    k11 = kl(0,3);    k12 = kl(0,4);    k13 = kl(0,5);
    k21 = kl(1,3);    k22 = kl(1,4);    k23 = kl(1,5);
    k31 = kl(2,3);    k32 = kl(2,4);    k33 = kl(2,5);
    
    kg(0,3) = c2*k11+s2*k22-cs*(k21+k12); 
    kg(1,3) = c2*k21-s2*k12+cs*(k11-k22);
    kg(2,3) = cosTheta*k31-sinTheta*k32;
    
    kg(0,4) = c2*k12-s2*k21+cs*(k11-k22);
    kg(1,4) = c2*k22+s2*k11+cs*(k21+k12);
    kg(2,4) = sinTheta*k31+cosTheta*k32;
    
    kg(0,5) = cosTheta*k13 - sinTheta*k23;
    kg(1,5) = sinTheta*k13 + cosTheta*k23;
    kg(2,5) = k33;
    
    k11 = kl(3,0);    k12 = kl(3,1);    k13 = kl(3,2);
    k21 = kl(4,0);    k22 = kl(4,1);    k23 = kl(4,2);
    k31 = kl(5,0);    k32 = kl(5,1);    k33 = kl(5,2);
    
    kg(3,0) = c2*k11+s2*k22-cs*(k21+k12); 
    kg(4,0) = c2*k21-s2*k12+cs*(k11-k22);
    kg(5,0) = cosTheta*k31-sinTheta*k32;
    
    kg(3,1) = c2*k12-s2*k21+cs*(k11-k22);
    kg(4,1) = c2*k22+s2*k11+cs*(k21+k12);
    kg(5,1) = sinTheta*k31+cosTheta*k32;
    
    kg(3,2) = cosTheta*k13 - sinTheta*k23;
    kg(4,2) = sinTheta*k13 + cosTheta*k23;
    kg(5,2) = k33;
    
    k11 = kl(3,3);    k12 = kl(3,4);    k13 = kl(3,5);
    k21 = kl(4,3);    k22 = kl(4,4);    k23 = kl(4,5);
    k31 = kl(5,3);    k32 = kl(5,4);    k33 = kl(5,5);
    
    kg(3,3) = c2*k11+s2*k22-cs*(k21+k12); 
    kg(4,3) = c2*k21-s2*k12+cs*(k11-k22);
    kg(5,3) = cosTheta*k31-sinTheta*k32;
    
    kg(3,4) = c2*k12-s2*k21+cs*(k11-k22);
    kg(4,4) = c2*k22+s2*k11+cs*(k21+k12);
    kg(5,4) = sinTheta*k31+cosTheta*k32;
    
    kg(3,5) = cosTheta*k13 - sinTheta*k23;
    kg(4,5) = sinTheta*k13 + cosTheta*k23;
    kg(5,5) = k33;

    if (nodeOffsets == true) {
      double X1 = nodeIOffset(0);
      double Y1 = nodeIOffset(1);
      double X2 = nodeJOffset(0);
      double Y2 = nodeJOffset(1);

      double k11 = kg(0,0);

      double k12 = kg(0,1);
      double k22 = kg(1,1);

      double k13 = kg(0,2);
      double k23 = kg(1,2);
      double k33 = kg(2,2);

      double k14 = kg(0,3);
      double k24 = kg(1,3);
      double k34 = kg(2,3);
      double k44 = kg(3,3);

      double k15 = kg(0,4);
      double k25 = kg(1,4);
      double k35 = kg(2,4);
      double k45 = kg(3,4);
      double k55 = kg(3,4);

      double k16 = kg(0,5);
      double k26 = kg(1,5);
      double k36 = kg(2,5);
      double k46 = kg(3,5);
      double k56 = kg(4,5);
      double k66 = kg(5,5);


      double K13 = -k11*Y1 + k12*X1 + k13;
      double K23 = -k12*Y1 + k22*X1 + k23;
      kg(0,2) = K13;
      kg(2,0) = K13;
      kg(1,2) = K23;
      kg(2,1) = K23;
      kg(2,2) = -Y1*K13 + X1*K23 -k13*Y1 + k23*X1 + k33;
      

      double K16 = -k14*Y2 + k15*X2 + k16;
      double K26 = -k24*Y2 + k25*X2 + k26;
      kg(0,5) = K16;
      kg(5,0) = K16;
      kg(1,5) = K26;
      kg(5,1) = K26;
      kg(2,5) = -Y2*K16 + X2*K26 -k16*Y1 + k26*X1 + k36;      
      kg(5,2) = kg(2,5);

      double K46 = -k44*Y2 + k45*X2 + k46;
      double K56 = -k45*Y2 + k55*X2 + k56;
      kg(3,5) = K46;
      kg(5,3) = K46;
      kg(4,5) = K56;
      kg(5,4) = K56;
      kg(5,5) = -Y2*K46 + X2*K56 -k46*Y2 + k56*X2 + k66;      

      double K34 = -k14*Y1 + k24*X1 + k34;
      double K35 = -k15*Y1 + k25*X1 + k35;
      kg(2,3) = K34;
      kg(3,2) = K34;
      kg(2,4) = K35;
      kg(4,2) = K35;

    }

    
    return kg;
}


const Matrix &
CorotCrdTransf2d::getGeomStiffMatrix(const Vector &pb) const
{
    // get  geometric stiffness matrix present in the transformation 
    // from basic to local system
    double s2, c2, cs;  
    
    s2 = sinAlpha*sinAlpha;
    c2 = cosAlpha*cosAlpha;
    cs = sinAlpha*cosAlpha;
    
    static Matrix kg0(6,6), kg12(6,6);
    kg0.Zero();
    
    kg12.Zero();
    
    kg0(0,0) = kg0(3,3) =  s2;
    kg0(0,1) = kg0(3,4) = -cs;
    kg0(1,0) = kg0(4,3) = -cs;
    kg0(1,1) = kg0(4,4) =  c2;
    
    kg0(0,3) = kg0(3,0) = -s2;
    kg0(0,4) = kg0(3,1) =  cs;
    kg0(1,3) = kg0(4,0) =  cs;
    kg0(1,4) = kg0(4,1) = -c2;
    
    kg0 *= pb(0)/Ln;
    
    kg12(0,0) = kg12(3,3) = -2*cs;
    kg12(0,1) = kg12(3,4) =  c2-s2;
    kg12(1,0) = kg12(4,3) =  c2-s2;
    kg12(1,1) = kg12(4,4) =  2*cs;
    
    kg12(0,3) = kg12(3,0) =  2*cs;
    kg12(0,4) = kg12(3,1) = -c2+s2;
    kg12(1,3) = kg12(4,0) = -c2+s2;
    kg12(1,4) = kg12(4,1) = -2*cs;
    
    kg12 *= (pb(1)+pb(2))/(Ln*Ln);
    
    static Matrix kg(6,6);
    // kg = kg0 + kg12;
    kg = kg0;
    kg.addMatrix(1.0, kg12, 1.0);
    
    return kg;
}


double 
CorotCrdTransf2d::getInitialLength(void)
{
    return L;
}


double 
CorotCrdTransf2d::getDeformedLength(void)
{
    return Ln;
}


CrdTransf2d *
CorotCrdTransf2d::getCopy(void)
{
    // create a new instance of CorotCrdTransf2d 
    CorotCrdTransf2d *theCopy = new CorotCrdTransf2d (this->getTag(), nodeIOffset, nodeJOffset);
    
    if (!theCopy)
    {
        opserr << "CorotCrdTransf2d::getCopy() - out of memory creating copy\n";
        exit(-1);
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
CorotCrdTransf2d::sendSelf(int cTag, Channel &theChannel)
{
    Vector data(13);
    data(0) = ubcommit(0);
    data(1) = ubcommit(1);
    data(2) = ubcommit(2);
    data(3) = nodeIOffset(0);
    data(4) = nodeIOffset(1);
    data(5) = nodeJOffset(0);
    data(6) = nodeJOffset(1);
    
    if (nodeIInitialDisp != 0) {
        data(7) = nodeIInitialDisp[0];
        data(8) = nodeIInitialDisp[1];
        data(9) = nodeIInitialDisp[2];
    } else {
        data(7) = 0.0;
        data(8) = 0.0;
        data(9) = 0.0;
    }
    
    if (nodeJInitialDisp != 0) {
        data(10) = nodeJInitialDisp[0];
        data(11) = nodeJInitialDisp[1];
        data(12) = nodeJInitialDisp[2];
    } else {
        data(10) = 0.0;
        data(11) = 0.0;
        data(12) = 0.0;
    }
    
    if (theChannel.sendVector(this->getTag(), cTag, data) < 0) {
        opserr << " CorotCrdTransf2d::sendSelf() - data could not be sent\n" ;
        return -1;
    }
    return 0;
}


int 
CorotCrdTransf2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    Vector data(13);
    if (theChannel.recvVector(this->getTag(), cTag, data) < 0) {
        opserr << " CorotCrdTransf2d::recvSelf() - data could not be received\n" ;
        return -1;
    }
    
    ubcommit(0) = data(0);
    ubcommit(1) = data(1);
    ubcommit(2) = data(2);
    nodeIOffset(0) =data(4);
    nodeIOffset(1) =data(4);
    nodeJOffset(0) =data(5);
    nodeJOffset(1) =data(6);
    
    int flag, i, j;
    flag = 0;
    for (i=7; i<=9; i++)
        if (data(i) != 0.0)
            flag = 1;
        if (flag == 1) {
            if (nodeIInitialDisp == 0)
                nodeIInitialDisp = new double[3];
            for (i=7, j=0; i<=9; i++, j++)
                nodeIInitialDisp[j] = data(i);
        }
        
        flag = 0;
        for (i=10; i<=13; i++)
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


const Vector &
CorotCrdTransf2d::getPointGlobalCoordFromLocal(const Vector &xl)
{
    static Vector xg(3);
    opserr << " CorotCrdTransf2d::getPointGlobalCoordFromLocal: not implemented yet" ;
    
    return xg;  
}


const Vector &
CorotCrdTransf2d::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
    opserr << " CorotCrdTransf2d::getPointGlobalDisplFromBasic: not implemented yet" ;
    
    return uxg;  
}


void
CorotCrdTransf2d::Print(OPS_Stream &s, int flag)
{
    s << "\nCrdTransf: " << this->getTag() << " Type: LinearCrdTransf2d";
    s << "\tnodeI Offset: " << nodeIOffset;
    s << "\tnodeJ Offset: " << nodeJOffset;
}
