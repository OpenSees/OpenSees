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
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CorotCrdTransf3d.cpp,v $

// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 05/2000
// Revision: rms 06/2000 (using Assemble, and AssembleTranspose)
//
// Modified: 04/2005 Andreas Schellenberg (getBasicTrialVel, getBasicTrialAccel)
// 
// Purpose: This file contains the implementation for the 
// CorotCrdTransf3d class. CorotCrdTransf3d is a Corot
// transformation for a spatial frame between the global 
// and basic coordinate systems.

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>
#include <elementAPI.h>
#include <string>
#include <CorotCrdTransf3d.h>

// initialize static variables
Matrix CorotCrdTransf3d::Tp(6,7); 
Matrix CorotCrdTransf3d::Tlg(12,12);
Matrix CorotCrdTransf3d::TlgInv(12, 12);
Matrix CorotCrdTransf3d::Tbl(6,12);
Matrix CorotCrdTransf3d::kg(12,12);

void* OPS_CorotCrdTransf3d()
{
    if(OPS_GetNumRemainingInputArgs() < 4) {
	opserr<<"insufficient arguments for CorotCrdTransf3d\n";
	return 0;
    }

    // get tag
    int tag;
    int numData = 1;
    if(OPS_GetIntInput(&numData,&tag) < 0) return 0;

    // get vector
    Vector vec(3);
    double* vptr = &vec(0);
    numData = 3;
    if(OPS_GetDoubleInput(&numData,vptr) < 0) return 0;

    // get option
    Vector jntOffsetI(3), jntOffsetJ(3);
    double *iptr=&jntOffsetI(0), *jptr=&jntOffsetJ(0);
    while(OPS_GetNumRemainingInputArgs() > 6) {
	std::string type = OPS_GetString();
	if(type == "-jntOffset") {
	    if(OPS_GetDoubleInput(&numData,iptr) < 0) return 0;
	    if(OPS_GetDoubleInput(&numData,jptr) < 0) return 0;
	}
    }

    return new CorotCrdTransf3d(tag,vec,jntOffsetI,jntOffsetJ);
}


// constructor:
CorotCrdTransf3d::CorotCrdTransf3d(int tag, const Vector &vecInLocXZPlane,
                                   const Vector &rigJntOffsetI,
                                   const Vector &rigJntOffsetJ):
CrdTransf(tag, CRDTR_TAG_CorotCrdTransf3d),
vAxis(3), nodeIOffset(3), nodeJOffset(3), xAxis(3),
nodeIPtr(0), nodeJPtr(0), R0(3,3), L(0), Ln(0), 
alphaIq(4), alphaJq(4), 
alphaIqcommit(4), alphaJqcommit(4), alphaI(3), alphaJ(3),
ulcommit(7), ul(7),  ulpr(7),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false),
RI(3,3), RJ(3,3), Rbar(3,3), e(3,3), T(7,12), T0(7,12), Lr2(12,3), Lr3(12,3), A(3,3)
{
    // check vector that defines local xz plane
    if (vecInLocXZPlane.Size() != 3 )
    {
        opserr << "CorotCrdTransf3d::CorotCrdTransf3d:  Vector that defines local xz plane is invalid\n";
        opserr << "Size must be 3\n. Using (0,0,1)";      
        vAxis(0) = 0;       vAxis(1) = 0;      vAxis(2) = 1;
    }
    else
        vAxis = vecInLocXZPlane;
    
    // check rigid joint offset for node I
    if (rigJntOffsetI.Size() != 3 )
    {
        opserr << "CorotCrdTransf3d::CorotCrdTransf3d:  Invalid rigid joint offset vector for node I\n";
        opserr << "Size must be 3\n";      
        nodeIOffset.Zero();      
    }
    else
        nodeIOffset = rigJntOffsetI;
    
    // check rigid joint offset for node J
    if (rigJntOffsetJ.Size() != 3 )
    {
        opserr << "CorotCrdTransf3d::CorotCrdTransf3d:  Invalid rigid joint offset vector for node J\n";
        opserr << "Size must be 3\n";      
        nodeJOffset.Zero(); 
    }
    else
        nodeJOffset = rigJntOffsetJ;
    
    //opserr << "nodeIOffset:" << nodeIOffset;
    //opserr << "nodeJOffset:" << nodeJOffset;
    
    // temporary
    if (nodeIOffset.Norm() != 0 || nodeJOffset.Norm() != 0)
    {
        opserr << "CorotCrdTransf3d::CorotCrdTransf3d: rigid joint zones not implemented yet\n";
        opserr << "Using zero values\n"; 
        nodeIOffset.Zero();
        nodeJOffset.Zero();
    }
    
    // Permutation matrix (to renumber basic dof's)
    
    //       0 1  2 3 4  5 6  
    //
    // Tp=  [0 0  0 0 0  0 1;  0
    //       0 1  0 0 0  0 0;  1            
    //       0 0  0 0 1  0 0;  2
    //       0 0 -1 0 0  0 0;  3
    //       0 0  0 0 0 -1 0;  4
    //      -1 0  0 1 0  0 0]; 5 
    
    // using static matrix (one constant matrix for all objects)
    
    if (Tp(0, 6) == 0)  // initialize only once
    {
        Tp(0, 6) =  1;
        Tp(1, 1) =  1;
        Tp(2, 4) =  1;
        Tp(3, 2) = -1;
        Tp(4, 5) = -1;
        Tp(5, 0) = -1;
        Tp(5, 3) =  1;
    }
    
    //opserr << "Tp: " << Tp;
    
}  


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
CorotCrdTransf3d::CorotCrdTransf3d():
CrdTransf(0, CRDTR_TAG_CorotCrdTransf3d),
vAxis(3), nodeIOffset(3), nodeJOffset(3), xAxis(3),
nodeIPtr(0), nodeJPtr(0), R0(3,3), L(0), Ln(0), 
alphaIq(4), alphaJq(4), 
alphaIqcommit(4), alphaJqcommit(4), alphaI(3), alphaJ(3),
ulcommit(7), ul(7),  ulpr(7),
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false),
RI(3,3), RJ(3,3), Rbar(3,3), e(3,3), T(7,12), T0(7,12), Lr2(12,3), Lr3(12,3), A(3,3)
{
    // Permutation matrix (to renumber basic dof's)
    
    //       0 1  2 3 4  5 6  
    //
    // Tp=  [0 0  0 0 0  0 1;  0
    //       0 1  0 0 0  0 0;  1            
    //       0 0  0 0 1  0 0;  2
    //       0 0 -1 0 0  0 0;  3
    //       0 0  0 0 0 -1 0;  4
    //      -1 0  0 1 0  0 0]; 5 
    
    // using static matrix (one constant matrix for all objects)
    
    if (Tp(0, 6) == 0)  // initialize only once
    {
        Tp(0, 6) =  1;
        Tp(1, 1) =  1;
        Tp(2, 4) =  1;
        Tp(3, 2) = -1;
        Tp(4, 5) = -1;
        Tp(5, 0) = -1;
	Tp(5, 3) =  1;
    }
    
    //opserr << "Tp: " << Tp;
}


// destructor:
CorotCrdTransf3d::~CorotCrdTransf3d() 
{
    if (nodeIInitialDisp != 0)
        delete [] nodeIInitialDisp;
    if (nodeJInitialDisp != 0)
        delete [] nodeJInitialDisp;
}


int
CorotCrdTransf3d::commitState(void)
{
    ulcommit      = ul;
    alphaIqcommit = alphaIq;
    alphaJqcommit = alphaJq;
    
    return 0;
}


int
CorotCrdTransf3d::revertToLastCommit(void)
{
    // determine global displacement increments from last iteration
    const Vector &dispI = nodeIPtr->getTrialDisp();
    const Vector &dispJ = nodeJPtr->getTrialDisp();
    
    for (int k = 0; k < 3; k++) {
        alphaI(k) =  dispI(k+3);
        alphaJ(k) =  dispJ(k+3);
    }   
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<3; j++)
            alphaI(j) -= nodeIInitialDisp[j+3];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<3; j++)
            alphaJ(j) -= nodeJInitialDisp[j+3];
    }   
    
    ul      = ulcommit;
    alphaIq = alphaIqcommit;
    alphaJq = alphaJqcommit;
    
    this->update();
    
    return 0;
}


int
CorotCrdTransf3d::revertToStart(void)
{
    ul.Zero();
    alphaIq = this->getQuaternionFromRotMatrix(R0);    // pseudo-vector for node 1
    alphaJq = this->getQuaternionFromRotMatrix(R0);    // pseudo-vector for node J
    
    alphaI.Zero();
    alphaJ.Zero();

    // also reset the committed state, otherwise getBasicIncrDisp() and a
    // later revertToLastCommit() would use the old committed values
    alphaIqcommit = alphaIq;
    alphaJqcommit = alphaJq;
    ulcommit.Zero();
    ulpr.Zero();

    this->update();
    return 0;
}


int 
CorotCrdTransf3d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
    int error;
    
    nodeIPtr = nodeIPointer;
    nodeJPtr = nodeJPointer;
    
    if ((!nodeIPtr) || (!nodeJPtr))
    {
        opserr << "\nCorotCrdTransf3d::initialize";
        opserr << "\ninvalid pointers to the element nodes\n";
        return -1;
    }
    
    // see if there is some initial displacements at nodes
    if (initialDispChecked == false) {
        const Vector &nodeIDisp = nodeIPtr->getDisp();
        const Vector &nodeJDisp = nodeJPtr->getDisp();
        for (int i=0; i<6; i++)
	  if (nodeIDisp(i) != 0.0) {
	    nodeIInitialDisp = new double [6];
	    for (int j=0; j<6; j++)
	      nodeIInitialDisp[j] = nodeIDisp(j);
	    i = 6;
	  }
            
	for (int j=0; j<6; j++)
	  if (nodeJDisp(j) != 0.0) {
	    nodeJInitialDisp = new double [6];
	    for (int i=0; i<6; i++)
	      nodeJInitialDisp[i] = nodeJDisp(i);
	    j = 6;
	  }
	initialDispChecked = true;
    }
    
    static Vector XAxis(3);
    static Vector YAxis(3);
    static Vector ZAxis(3);
    
    // get 3by3 rotation matrix
    if ((error = this->getLocalAxes(XAxis, YAxis, ZAxis)))
      return error;
    
    // transformation matrix in the undeformed configuration (used by the
    // initial stiffness): all triads equal R0, Ln = L, ul = 0
    {
      static Vector ulSave(7);
      ulSave = ul;
      ul.Zero();
      RI = R0;  RJ = R0;  Rbar = R0;  e = R0;
      Ln = L;
      this->compTransfMatrixBasicGlobal();
      T0 = T;
      ul = ulSave;
    }
    
    // compute initial pseudo-vectors for nodal triads
    //opserr << setiosflags(ios::scientific);
    //opserr << setiosflags(ios::showpos);
    //opserr << setprecision(16);
    
    //opserr << "L: " << L;
    //opserr << "R0: " << R0;
    
    // If the committed nodal quaternions are already set (state restored by
    // recvSelf(), or initialize() called again), keep that state instead of
    // resetting the nodal triads to R0. A valid quaternion is never zero, so a
    // zero alphaIqcommit means this object has not been initialized yet.
    if (alphaIqcommit.Norm() != 0.0 && alphaJqcommit.Norm() != 0.0) {
        alphaIq = alphaIqcommit;
        alphaJq = alphaJqcommit;
        ul      = ulcommit;
        return 0;
    }

    alphaIq = this->getQuaternionFromRotMatrix(R0);    // pseudo-vector for node I
    alphaJq = this->getQuaternionFromRotMatrix(R0);    // pseudo-vector for node J

    //opserr << "alphaIq: " << alphaIq;
    //opserr << "alphaJq: " << alphaJq;

    this->commitState();

    return 0;
}


int
CorotCrdTransf3d::update(void)
{
    int i, j, k;
    
    /********* OLD REMO - REPLACED BELOW TO FIX BUG ***************
    // determine global displacement increments from last iteration
    const Vector &dispIncrI = nodeIPtr->getIncrDeltaDisp();
    const Vector &dispIncrJ = nodeJPtr->getIncrDeltaDisp();
    
     // get the iterative spins dAlphaI and dAlphaJ 
     // (rotational displacement increments at both nodes)
     
      static Vector dAlphaI(3);
      static Vector dAlphaJ(3);
      
       
        for (k = 0; k < 3; k++)
        {
        dAlphaI(k) = dispIncrI(k+3);
        dAlphaJ(k) = dispIncrJ(k+3);
        }   
    **************************************************************/
    
    // determine global displacement increments from last iteration
    static Vector dispI(6);
    static Vector dispJ(6);
    dispI = nodeIPtr->getTrialDisp();
    dispJ = nodeJPtr->getTrialDisp();
    
    if (nodeIInitialDisp != 0) {
        for (int j=0; j<6; j++)
            dispI(j) -= nodeIInitialDisp[j];
    }
    
    if (nodeJInitialDisp != 0) {
        for (int j=0; j<6; j++)
            dispJ(j) -= nodeJInitialDisp[j];
    }
    
    // get the iterative spins dAlphaI and dAlphaJ 
    // (rotational displacement increments at both nodes)
    
    static Vector dAlphaI(3);
    static Vector dAlphaJ(3);
    
    for (k = 0; k < 3; k++) {
        dAlphaI(k) = dispI(k+3) - alphaI(k);
        dAlphaJ(k) = dispJ(k+3) - alphaJ(k);
        alphaI(k) =  dispI(k+3);
        alphaJ(k) =  dispJ(k+3);
    }   
    
    /************** END OF REPLACEMENT **************************/
    
    // update the nodal triads TI and RJ using quaternions
    static Vector dAlphaIq(4);
    static Vector dAlphaJq(4);

    dAlphaIq = this->getQuaternionFromPseudoRotVector (dAlphaI);
    dAlphaJq = this->getQuaternionFromPseudoRotVector (dAlphaJ);

    alphaIq = this->quaternionProduct (alphaIq, dAlphaIq);
    alphaJq = this->quaternionProduct (alphaJq, dAlphaJq);

    RI = this->getRotationMatrixFromQuaternion (alphaIq);
    RJ = this->getRotationMatrixFromQuaternion (alphaJq);

    // compute the mean nodal triad
    static Matrix dRgamma(3,3); 
    static Vector gammaq(4);
    static Vector gammaw(3);
    
    dRgamma.Zero();
    
    //dRgamma = RJ * RIt;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            for (k = 0; k < 3; k++) 
                dRgamma(i,j) += RJ(i,k) * RI(j,k);
            
            gammaq = this->getQuaternionFromRotMatrix(dRgamma);    // pseudo-vector for node J
            gammaw = this->getTangScaledPseudoVectorFromQuaternion(gammaq);
            
            // opserr << "dRgamma: " << dRgamma;
            // opserr << "gammaq: " << gammaq << endln;
            // opserr << "gammaw: " << gammaw << endln;
            
            dRgamma = this->getRotMatrixFromTangScaledPseudoVector(gammaw/2);
            
            Rbar.addMatrixProduct(0.0, dRgamma, RI, 1.0);
            
            // compute the base vectors e1, e2, e3
            static Vector e1(3);
            static Vector e2(3);
            static Vector e3(3);
            
            // relative translation displacements
            static Vector dJI(3);    
            for (int kk = 0; kk < 3; kk++)
                dJI(kk) = dispJ(kk) - dispI(kk);
            
            // element projection
            static Vector xJI(3);
            xJI = nodeJPtr->getCrds() - nodeIPtr->getCrds();
            
            if (nodeIInitialDisp != 0) {
                xJI(0) -= nodeIInitialDisp[0];
                xJI(1) -= nodeIInitialDisp[1];
                xJI(2) -= nodeIInitialDisp[2];
            }
            
            if (nodeJInitialDisp != 0) {
                xJI(0) += nodeJInitialDisp[0];
                xJI(1) += nodeJInitialDisp[1];
                xJI(2) += nodeJInitialDisp[2];
            }
            
            static Vector dx(3);
            // dx = xJI + dJI;  
            dx = xJI;
            dx.addVector (1.0, dJI, 1.0);
            
            // calculate the deformed element length
            Ln = dx.Norm();
            
            
            if (Ln == 0.0) 
            {
                opserr << "\nCorotCrdTransf3d::update: 0 deformed length\n";
                return -2;  
            }
            
            // compute the base vector e1
            e1 = dx/Ln;
            
            // 'rotate' the mean rotation matrix Rbar on to e1 to 
            // obtain e2 and e3 (using the 'mid-point' procedure)
            static Vector r1(3);
            static Vector r2(3);
            static Vector r3(3);
            
            for (k = 0; k < 3; k ++)
            {
                r1(k) = Rbar(k,0);
                r2(k) = Rbar(k,1);
                r3(k) = Rbar(k,2);
            }
            
            //    e2 = r2 - (e1 + r1)*((r2^ e1)*0.5);
            // e3 = r3 - (e1 + r1)*((r3^ e1)*0.5);
            
            static Vector tmp(3);
            tmp = e1;
            tmp += r1;
            
            e2=tmp;
            e3=tmp;
            
            e2 *= (r2^e1)*0.5;
            e2.addVector(-1.0,  r2, 1.0);
            
            // e2 = r2 - (e1 + r1)*((r2^ e1)*0.5);
            
            e3 *= (r3^e1)*0.5;
            e3.addVector(-1.0,  r3, 1.0);
            
            // compute the basic rotations
            static Vector rI1(3), rI2(3), rI3(3);
            static Vector rJ1(3), rJ2(3), rJ3(3);
            
            for (k = 0; k < 3; k ++)
            {
                e(k,0) = e1(k);
                e(k,1) = e2(k);
                e(k,2) = e3(k);
                
                rI1(k) = RI(k,0);
                rI2(k) = RI(k,1);
                rI3(k) = RI(k,2);
                
                rJ1(k) = RJ(k,0);
                rJ2(k) = RJ(k,1);
                rJ3(k) = RJ(k,2);
            }
            
            // compute the basic displacements
            ulpr = ul;
            ul(0) = asin (((rI2^ e3) - (rI3^ e2))*0.5);
            ul(1) = asin (((rI1^ e2) - (rI2^ e1))*0.5);
            ul(2) = asin (((rI1^ e3) - (rI3^ e1))*0.5);
            ul(3) = asin (((rJ2^ e3) - (rJ3^ e2))*0.5);
            ul(4) = asin (((rJ1^ e2) - (rJ2^ e1))*0.5);
            ul(5) = asin (((rJ1^ e3) - (rJ3^ e1))*0.5);		   
            
            // ul = Ln - L;
            // ul(6) = 2 * ((xJI + dJI/2)^ dJI) / (Ln + L);  //mid-point formula   
            xJI.addVector(1.0, dJI, 0.5);
            ul(6) = 2 * (xJI ^ dJI) / (Ln + L);  //mid-point formula   
            
            // compute the transformation matrix
            this->compTransfMatrixBasicGlobal();
            
            return 0;
}


void
CorotCrdTransf3d::compTransfMatrixBasicGlobal(void)
{
    // Computes A, Lr2, Lr3 and the transformation matrix T (basic <- global)
    // from RI, RJ, Rbar, e, Ln and ul (fixed-size arrays; same expressions
    // as the original Matrix/Vector implementation).
    int i, j, k;

    double r1[3], r2[3], r3[3], e1[3], e2[3], e3[3];
    double rI1[3], rI2[3], rI3[3], rJ1[3], rJ2[3], rJ3[3];
    for (k = 0; k < 3; k++) {
        r1[k] = Rbar(k,0);  r2[k] = Rbar(k,1);  r3[k] = Rbar(k,2);
        e1[k] = e(k,0);     e2[k] = e(k,1);     e3[k] = e(k,2);
        rI1[k] = RI(k,0);   rI2[k] = RI(k,1);   rI3[k] = RI(k,2);
        rJ1[k] = RJ(k,0);   rJ2[k] = RJ(k,1);   rJ3[k] = RJ(k,2);
    }

    //   A = (1/Ln)*(I - e1*e1');
    double a[3][3];
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            a[i][j] = ((i == j ? 1.0 : 0.0) - e1[i]*e1[j])/Ln;
            A(i,j) = a[i][j];
        }

    // L(ri) = [L1; L2; -L1; L2] with
    //   L1 = ri'*e1 * A/2 + A*ri*(e1 + r1)'/2
    //   L2 = Sri/2 - ri'*e1*S(r1)/4 - Sri*e1*(e1 + r1)'/4
    double lr2[12][3], lr3[12][3];
    {
        const double *ris[2] = {r2, r3};
        double (*Ls[2])[3] = {lr2, lr3};
        double e1r1[3];
        for (k = 0; k < 3; k++) e1r1[k] = e1[k] + r1[k];
        for (int q = 0; q < 2; q++) {
            const double *ri = ris[q];
            double (*L)[3] = Ls[q];
            const double rie1 = ri[0]*e1[0] + ri[1]*e1[1] + ri[2]*e1[2];
            double Ari[3], rixe1[3];
            for (i = 0; i < 3; i++)
                Ari[i] = a[i][0]*ri[0] + a[i][1]*ri[1] + a[i][2]*ri[2];
            rixe1[0] = ri[1]*e1[2] - ri[2]*e1[1];
            rixe1[1] = ri[2]*e1[0] - ri[0]*e1[2];
            rixe1[2] = ri[0]*e1[1] - ri[1]*e1[0];
            // skew matrices S(ri), S(r1)
            const double Sri[3][3] = {{0.0, -ri[2], ri[1]}, {ri[2], 0.0, -ri[0]}, {-ri[1], ri[0], 0.0}};
            const double Sr1[3][3] = {{0.0, -r1[2], r1[1]}, {r1[2], 0.0, -r1[0]}, {-r1[1], r1[0], 0.0}};
            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++) {
                    const double L1 = a[i][j]*rie1*0.5 + Ari[i]*e1r1[j]*0.5;
                    const double L2 = Sri[i][j]*0.5 - Sr1[i][j]*rie1/4.0 - rixe1[i]*e1r1[j]*0.25;
                    L[i][j]   =  L1;
                    L[3+i][j] =  L2;
                    L[6+i][j] = -L1;
                    L[9+i][j] =  L2;
                }
        }
    }
    for (i = 0; i < 12; i++)
        for (j = 0; j < 3; j++) {
            Lr2(i,j) = lr2[i][j];
            Lr3(i,j) = lr3[i][j];
        }

    double t[7][12];
    for (i = 0; i < 7; i++)
        for (j = 0; j < 12; j++)
            t[i][j] = 0.0;

    auto cross = [](const double *x, const double *y, double *z) {
        z[0] = x[1]*y[2] - x[2]*y[1];
        z[1] = x[2]*y[0] - x[0]*y[2];
        z[2] = x[0]*y[1] - x[1]*y[0];
    };
    double c1[3], c2[3], At[3];

    //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';
    cross(rI3, e2, c1);  cross(rI2, e3, c2);
    for (i = 0; i < 3; i++) t[0][i+3] = -c1[i] + c2[i];

    //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';
    cross(rI2, e1, c1);  cross(rI1, e2, c2);
    for (i = 0; i < 3; i++) At[i] = a[i][0]*rI2[0] + a[i][1]*rI2[1] + a[i][2]*rI2[2];
    for (i = 0; i < 3; i++) { t[1][i] = At[i]; t[1][i+3] = -c1[i] + c2[i]; t[1][i+6] = -At[i]; }

    //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';
    cross(rI3, e1, c1);  cross(rI1, e3, c2);
    for (i = 0; i < 3; i++) At[i] = a[i][0]*rI3[0] + a[i][1]*rI3[1] + a[i][2]*rI3[2];
    for (i = 0; i < 3; i++) { t[2][i] = At[i]; t[2][i+3] = -c1[i] + c2[i]; t[2][i+6] = -At[i]; }

    //   T4 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
    cross(rJ3, e2, c1);  cross(rJ2, e3, c2);
    for (i = 0; i < 3; i++) t[3][i+9] = -c1[i] + c2[i];

    //   T5 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
    cross(rJ2, e1, c1);  cross(rJ1, e2, c2);
    for (i = 0; i < 3; i++) At[i] = a[i][0]*rJ2[0] + a[i][1]*rJ2[1] + a[i][2]*rJ2[2];
    for (i = 0; i < 3; i++) { t[4][i] = At[i]; t[4][i+6] = -At[i]; t[4][i+9] = -c1[i] + c2[i]; }

    //   T6 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']';
    cross(rJ3, e1, c1);  cross(rJ1, e3, c2);
    for (i = 0; i < 3; i++) At[i] = a[i][0]*rJ3[0] + a[i][1]*rJ3[1] + a[i][2]*rJ3[2];
    for (i = 0; i < 3; i++) { t[5][i] = At[i]; t[5][i+6] = -At[i]; t[5][i+9] = -c1[i] + c2[i]; }

    // T(:,1) += Lr3*rI2 - Lr2*rI3;  T(:,2) += Lr2*rI1;  T(:,3) += Lr3*rI1;
    // T(:,4) += Lr3*rJ2 - Lr2*rJ3;  T(:,5) += Lr2*rJ1;  T(:,6) += Lr3*rJ1;
    for (i = 0; i < 12; i++) {
        const double *l2 = lr2[i], *l3 = lr3[i];
        t[0][i] += (l3[0]*rI2[0] + l3[1]*rI2[1] + l3[2]*rI2[2]) - (l2[0]*rI3[0] + l2[1]*rI3[1] + l2[2]*rI3[2]);
        t[1][i] +=  l2[0]*rI1[0] + l2[1]*rI1[1] + l2[2]*rI1[2];
        t[2][i] +=  l3[0]*rI1[0] + l3[1]*rI1[1] + l3[2]*rI1[2];
        t[3][i] += (l3[0]*rJ2[0] + l3[1]*rJ2[1] + l3[2]*rJ2[2]) - (l2[0]*rJ3[0] + l2[1]*rJ3[1] + l2[2]*rJ3[2]);
        t[4][i] +=  l2[0]*rJ1[0] + l2[1]*rJ1[1] + l2[2]*rJ1[2];
        t[5][i] +=  l3[0]*rJ1[0] + l3[1]*rJ1[1] + l3[2]*rJ1[2];
    }

    for (j = 0; j < 6; j++) {
        const double c = 2 * cos(ul(j));
        for (i = 0; i < 12; i++)
            t[j][i] /= c;
    }

    // T(:,7) = [-e1' O' e1' O']';
    for (i = 0; i < 3; i++) {
        t[6][i  ] = -e1[i];
        t[6][i+6] =  e1[i];
    }

    for (i = 0; i < 7; i++)
        for (j = 0; j < 12; j++)
            T(i,j) = t[i][j];

    (void)r2; (void)r3;
}


// NOTE: not called anywhere; experimental/incomplete alternative to
// compTransfMatrixBasicGlobal(). Do not use without verification.
void
CorotCrdTransf3d::compTransfMatrixBasicGlobalNew(void)
{
    // extract columns of rotation matrices
    int i, j, k;
    
    //opserr << "comprTransfMatrixBasicGlobal: *****************************\n";
    static Vector r1(3), r2(3), r3(3);
    static Vector e1(3), e2(3), e3(3);
    static Vector rI1(3), rI2(3), rI3(3);
    static Vector rJ1(3), rJ2(3), rJ3(3);
    
    for (k = 0; k < 3; k ++)
    {
        r1(k) = Rbar(k,0);
        r2(k) = Rbar(k,1);
        r3(k) = Rbar(k,2);
        
        e1(k) = e(k,0);
        e2(k) = e(k,1);
        e3(k) = e(k,2);
        
        rI1(k) = RI(k,0);
        rI2(k) = RI(k,1);
        rI3(k) = RI(k,2);
        
        rJ1(k) = RJ(k,0);
        rJ2(k) = RJ(k,1);
        rJ3(k) = RJ(k,2);
    }
    
    // compute the transformation matrix from the basic to the
    // global system
    static Matrix I(3,3);
    
    //   A = (1/Ln)*(I - e1*e1');
    for (i = 0; i < 3; i++)
        I(i,i) = 1;
    
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            A(i,j) = (I(i,j) - e1(i)*e1(j))/Ln;
        
        
            /*
            opserr << "r1: " << r1;
            opserr << "r2: " << r2;
            opserr << "r3: " << r3;
            opserr << "e1: " << e1;
            opserr << "e2: " << e2;
            opserr << "e3: " << e3;
            opserr << "rI1: " << rI1;
            opserr << "rI2: " << rI2;
            opserr << "rI3: " << rI3;
            opserr << "rJ1: " << rJ1;
            opserr << "rJ2: " << rJ2;
            opserr << "rJ3: " << rJ3;
            opserr << "A: " << A;
        */
        
        Lr2 = this->getLMatrix (r2);
        Lr3 = this->getLMatrix (r3);
        
        // opserr << "Lr2: " << Lr2;
        // opserr << "Lr3: " << Lr3;
        
        static Matrix Sr1(3,3), Sr2(3,3), Sr3(3,3);
        static Vector Se(3), At(3);
        
        
        // O = zeros(3,1);
        // hI1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';
        // hI2 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';
        // hI3 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';
        // hJ1 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
        // hJ2 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']';
        // hJ3 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
        
        static Vector hI1(12);
        static Vector hI2(12);
        static Vector hI3(12);
        static Vector hJ1(12);
        static Vector hJ2(12);
        static Vector hJ3(12);
        
        Sr1 = this->getSkewSymMatrix(rI1);
        Sr2 = this->getSkewSymMatrix(rI2);
        Sr3 = this->getSkewSymMatrix(rI3);
        
        // hI1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';
        Se.addMatrixVector(0.0, Sr3, e2, -1.0);     // (-S(rI3)*e2 + S(rI2)*e3)
        Se.addMatrixVector(1.0, Sr2, e3,  1.0);    
        for (i=0; i<3; i++)
            hI1(i+3) = Se(i);
        
        // hI2 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';
        At.addMatrixVector(0.0, A, rI3, 1.0);   
        Se.addMatrixVector(0.0, Sr3, e1, -1.0);     // (-S(rI3)*e1 + S(rI1)*e3)
        Se.addMatrixVector(1.0, Sr1, e3,  1.0);
        
        for (i = 0; i < 3; i++) {
            hI2(i  ) =  At(i);  
            hI2(i+3) =  Se(i);
            hI2(i+6) = -At(i);
        }
        
        // hI3 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';
        At.addMatrixVector(0.0, A, rI2, 1.0);   
        Se.addMatrixVector(0.0, Sr2, e1, -1.0);     // (-S(rI2)*e1 + S(rI1)*e2)'
        Se.addMatrixVector(1.0, Sr1, e2,  1.0);
        
        for (i = 0; i < 3; i++) {
            hI3(i  ) =  At(i);  
            hI3(i+3) =  Se(i);
            hI3(i+6) = -At(i);
        }
        
        Sr1 = this->getSkewSymMatrix(rJ1);
        Sr2 = this->getSkewSymMatrix(rJ2);
        Sr3 = this->getSkewSymMatrix(rJ3);
        
        // hJ1 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';  
        Se.addMatrixVector(0.0, Sr3, e2, -1.0);    // -S(rJ3)*e2 + S(rJ2)*e3
        Se.addMatrixVector(1.0, Sr2, e3,  1.0);    
        for (i = 0; i < 3; i++)
            hJ1(i+9) =  Se(i);
        
        
        // hJ2 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']';
        At.addMatrixVector(0.0, A, rJ3, 1.0);   
        Se.addMatrixVector(0.0, Sr3, e1, -1.0);     // (-S(rJ3)*e1 + S(rJ1)*e3)
        Se.addMatrixVector(1.0, Sr1, e3,  1.0);
        
        for (i = 0; i < 3; i++) {
            hJ2(i  ) =  At(i);  
            hJ2(i+6) = -At(i);
            hJ2(i+9) =  Se(i);
        }
        
        // hJ3 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';   
        At.addMatrixVector(0.0, A, rJ2, 1.0);   
        Se.addMatrixVector(0.0, Sr2, e1, -1.0);     // (-S(rJ2)*e1 + S(rJ1)*e2)
        Se.addMatrixVector(1.0, Sr1, e2,  1.0);
        
        for (i = 0; i < 3; i++) {
            hJ3(i  ) =  At(i);  
            hJ3(i+6) = -At(i);
            hJ3(i+9) =  Se(i);
        }
        
        /*
        opserr << "hI1: " << hI1;
        opserr << "hI2: " << hI2;
        opserr << "hI3: " << hI3;
        opserr << "hj1: " << hJ1;
        opserr << "hj2: " << hJ2;
        opserr << "hj3: " << hJ3;
        */
        
        // f1 =  [-e1' O' e1' O'];
        // f2  = ( Lr2*rI1 + hI3)'./(2*(cos(thetalI(3))));
        // f3  = ( Lr2*rJ1 + hJ3)'./(2*(cos(thetalJ(3))));
        // f4  = (-Lr3*rI1 - hI2)'./(2*(cos(thetalI(2))));
        // f5  = (-Lr3*rJ1 - hJ2)'./(2*(cos(thetalJ(2))));
        // f6I = ( Lr3*rI2 - Lr2*rI3 + hI1)'./(2*(cos(thetalI(1))));
        // f6J = ( Lr3*rJ2 - Lr2*rJ3 + hJ1)'./(2*(cos(thetalJ(1))));
        
        // F = [f1;
        //     f2;
        //     f3;
        //     f4;
        //     f5;
        //     f6J-f6I];      
        
        // T = F'
        T.Zero();
        static Vector Lr(12);
        
        // f1 =  [-e1' O' e1' O'];
        for (i=0; i<3; i++) {
            T(i  ,0) = -e1(i);
            T(i+3,0) = e1(i);
        }
        
        static Vector thetaI(3);
        static Vector thetaJ(3);
        
        
        thetaI(0) = ul(0);
        thetaI(1) = -ul(2);
        thetaI(2) = ul(1);
        
        thetaJ(0) = ul(3);
        thetaJ(1) = -ul(5);
        thetaJ(2) = ul(4);
        
        
        
        double c;
        
        // f2  = ( Lr2*rI1 + hI3)'./(2*(cos(thetalI(3))));
        Lr.addMatrixVector(0.0, Lr2, rI1,  1.0);  
        Lr += hI3;
        c = 1.0/(2.0*cos(thetaI(2)));
        for (i=0; i<12; i++)
            T(1,i) = Lr(i)*c;
        
        // f3  = ( Lr2*rJ1 + hJ3)'./(2*(cos(thetalJ(3))));
        Lr.addMatrixVector(0.0, Lr2, rJ1,  1.0);  
        Lr += hJ3;
        c = 1.0/(2.0*cos(thetaJ(2)));   
        for (i=0; i<12; i++)
            T(2,i) = Lr(i)*c;
        
        // f4  = (-Lr3*rI1 - hI2)'./(2*(cos(thetalI(2))));
        Lr.addMatrixVector(0.0, Lr3, rI1,  -1.0);  
        Lr -= hI2;
        c = 1.0/(2.0*cos(thetaI(1)));   
        for (i=0; i<12; i++)
            T(3,i) = Lr(i)*c;
        
        // f5  = (-Lr3*rJ1 - hJ2)'./(2*(cos(thetalJ(2))));
        Lr.addMatrixVector(0.0, Lr3, rJ1,  -1.0);  
        Lr -= hJ2;
        c = 1.0/(2.0*cos(thetaJ(1)));   
        for (i=0; i<12; i++)
            T(4,i) = Lr(i)*c;
        
        // f6I = ( Lr3*rI2 - Lr2*rI3 + hI1)'./(2*(cos(thetalI(1))));
        Lr.addMatrixVector(0.0, Lr3, rI2,  1.0);
        Lr.addMatrixVector(1.0, Lr2, rI3, -1.0);
        Lr += hI1;
        c = 1.0/(2.0*cos(thetaI(0)));   
        for (i=0; i<12; i++)
            T(5,i) = Lr(i)*c;
        
        // f6J = ( Lr3*rJ2 - Lr2*rJ3 + hJ1)'./(2*(cos(thetalJ(1))));
        Lr.addMatrixVector(0.0, Lr3, rJ2,  1.0);
        Lr.addMatrixVector(1.0, Lr2, rJ3, -1.0);
        Lr += hJ1;
        c = 1.0/(2.0*cos(thetaJ(0)));
        for (i=0; i<12; i++)
            T(6,i) -= Lr(i)*c;
}


void
CorotCrdTransf3d::compTransfMatrixLocalGlobal(Matrix &Tlg) 
{
    // setup transformation matrix from local to global
    Tlg.Zero();
    
    Tlg(0,0) = Tlg(3,3) = Tlg(6,6) = Tlg(9,9)   = R0(0,0);
    Tlg(0,1) = Tlg(3,4) = Tlg(6,7) = Tlg(9,10)  = R0(1,0);
    Tlg(0,2) = Tlg(3,5) = Tlg(6,8) = Tlg(9,11)  = R0(2,0);
    Tlg(1,0) = Tlg(4,3) = Tlg(7,6) = Tlg(10,9)  = R0(0,1);
    Tlg(1,1) = Tlg(4,4) = Tlg(7,7) = Tlg(10,10) = R0(1,1);
    Tlg(1,2) = Tlg(4,5) = Tlg(7,8) = Tlg(10,11) = R0(2,1);
    Tlg(2,0) = Tlg(5,3) = Tlg(8,6) = Tlg(11,9)  = R0(0,2);
    Tlg(2,1) = Tlg(5,4) = Tlg(8,7) = Tlg(11,10) = R0(1,2);
    Tlg(2,2) = Tlg(5,5) = Tlg(8,8) = Tlg(11,11) = R0(2,2);
}


void
CorotCrdTransf3d::compTransfMatrixBasicLocal(Matrix &Tbl)
{
    // setup transformation matrix from basic to local
    Tbl.Zero();

    // first get transformation matrix from basic to global 
    static Matrix Tbg(6, 12);
    Tbg.addMatrixProduct(0.0, Tp, T, 1.0);

    // get inverse of transformation matrix from local to global
    this->compTransfMatrixLocalGlobal(Tlg);
    // Tlg.Invert(TlgInv);
    TlgInv.addMatrixTranspose(0.0, Tlg, 1.0);  // for square rot-matrix: Tlg^-1 = Tlg'

    // finally get transformation matrix from basic to local
    Tbl.addMatrixProduct(0.0, Tbg, TlgInv, 1.0);
}


const Vector &
CorotCrdTransf3d::getBasicTrialDisp(void)
{
    static Vector ub(6);
    
    // use transformation matrix to renumber the degrees of freedom
    //ub.addMatrixVector(0.0, Tp, ul, 1.0);
    ub(0) = ul(6);
	ub(1) = ul(1);
	ub(2) = ul(4);
	ub(3) = -ul(2);
	ub(4) = -ul(5);
	ub(5) = ul(3)-ul(0);
	
    return ub;    
}


const Vector &
CorotCrdTransf3d::getBasicIncrDeltaDisp(void)
{
    static Vector dub(6);
    static Vector dul(7);
    
    // dul = ul - ulpr;
    dul = ul;
    dul.addVector (1.0, ulpr, -1.0);
    
    // use transformation matrix to renumber the degrees of freedom
    //dub.addMatrixVector(0.0, Tp, dul, 1.0);
    dub(0) = dul(6);
	dub(1) = dul(1);
	dub(2) = dul(4);
	dub(3) = -dul(2);
	dub(4) = -dul(5);
	dub(5) = dul(3)-dul(0);
	
    return dub;        
}


const Vector &
CorotCrdTransf3d::getBasicIncrDisp(void)
{
    static Vector Dub(6);
    static Vector Dul(7);
    
    // Dul = ul - ulcommit;
    Dul = ul;
    Dul.addVector (1.0, ulcommit, -1.0);
    
    // use transformation matrix to renumber the degrees of freedom
    //Dub.addMatrixVector(0.0, Tp, Dul, 1.0);
    Dub(0) = Dul(6);
	Dub(1) = Dul(1);
	Dub(2) = Dul(4);
	Dub(3) = -Dul(2);
	Dub(4) = -Dul(5);
	Dub(5) = Dul(3)-Dul(0);
	
    return Dub;        
}


const Vector &
CorotCrdTransf3d::getBasicTrialVel(void)
{
    opserr << "WARNING CorotCrdTransf3d::getBasicTrialVel()"
        << " - has not been implemented yet. Returning zeros." << endln;
    
    static Vector dummy(6);
    return dummy;
}


const Vector &
CorotCrdTransf3d::getBasicTrialAccel(void)
{
    opserr << "WARNING CorotCrdTransf3d::getBasicTrialAccel()"
        << " - has not been implemented yet. Returning zeros." << endln;
    
    static Vector dummy(6);
    return dummy;
}


const Vector &
CorotCrdTransf3d::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
    // T, RI, RJ, Rbar, e are per-object and current (set by update())

    static Vector pg(12);
    pg.Zero();
    
    // if there are no element loads present
    if (p0.Norm() == 0.0) {
        // transform resisting forces from the basic system to local coordinates
        static Vector pl(7);
        //pl.addMatrixTransposeVector(0.0, Tp, pb, 1.0);    // pl = Tp ^ pb;
		pl(0) = -pb(5);
		pl(1) = pb(1);
		pl(2) = -pb(3);
		pl(3) = pb(5);
		pl(4) = pb(2);
		pl(5) = -pb(4);
		pl(6) = pb(0);

        // transform resisting forces from local to global coordinates
        pg.addMatrixTransposeVector(0.0, T, pl, 1.0);   // pg = T ^ pl; residual
    }

    // if there are element loads present
    else {
        // SLOWER!!!! THIS IS SAME APPROACH AS 2D CASE
        // ===========================================
        /* transform resisting forces from the basic system to local coordinates
        this->compTransfMatrixBasicLocal(Tbl);
        static Vector pl(12);
        pl.addMatrixTransposeVector(0.0, Tbl, pb, 1.0);    // pl = Tbl ^ pb;

        // add end forces due to element p0 loads
        // assuming member loads are in local system
        pl(0) += p0(0);
        pl(1) += p0(1);
        pl(7) += p0(2);
        pl(2) += p0(3);
        pl(8) += p0(4);

        // transform resisting forces from local to global coordinates
        pg.addMatrixTransposeVector(0.0, Tlg, pl, 1.0);   // pg = Tlg ^ pl;*/

        // FASTER!!!! TRANSFORM REACTIONS AND ADD AT END
        // =============================================
        // transform resisting forces from the basic system to local coordinates
        static Vector pl(7);
        //pl.addMatrixTransposeVector(0.0, Tp, pb, 1.0);    // pl = Tp ^ pb;
		pl(0) = -pb(5);
		pl(1) = pb(1);
		pl(2) = -pb(3);
		pl(3) = pb(5);
		pl(4) = pb(2);
		pl(5) = -pb(4);
		pl(6) = pb(0);
		
        // transform resisting forces from local to global coordinates
        pg.addMatrixTransposeVector(0.0, T, pl, 1.0);   // pg = T ^ pl; residual

        // add end forces due to element p0 loads
        // assuming member loads are in local system
        static Vector pl0(12), pg0(12);
        pl0.Zero();
        pl0(0) = p0(0);
        pl0(1) = p0(1);
        pl0(7) = p0(2);
        pl0(2) = p0(3);
        pl0(8) = p0(4);
        this->compTransfMatrixLocalGlobal(Tlg);
        pg0.addMatrixTransposeVector(0.0, Tlg, pl0, 1.0);
        pg.addVector(1.0, pg0, 1.0);
    }
    
    return pg;
}


const Matrix &
CorotCrdTransf3d::getGlobalStiffMatrix(const Matrix &kb, const Vector &pb)
{
    // T, RI, RJ, Rbar, e, A, Lr2, Lr3 are per-object and current (set by update())
    //
    //   ks = T'*kl*T + ks1 + T * diag (m .* tan(thetal))*T' + ...
    //        m(4)*(ks2r2t3_u3 + ks2r3u2_t2) + ...
    //        m(2)*ks2r2t1 + m(3)*ks2r3t1 + ...
    //        m(5)*ks2r2u1 + m(6)*ks2r3u1 + ...
    //        ks3 + ks3' + ks4 + ks5;
    //
    // (same terms as the original implementation, evaluated with fixed-size
    // arrays; skew-matrix products use x'*S(r) = (x cross r)' and
    // S(a)*S(b) = b*a' - (a.b)*I)

    int i, j, k;

    // basic -> local: Tp has one nonzero (+-1) per column
    //   local dof a = (basic dof tpRow[a]) * tpSgn[a]
    static const int    tpRow[7] = { 5,   1,  3,   5,   2,  4,   0 };
    static const double tpSgn[7] = {-1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0};

    double kl[7][7], pl[7];
    for (i = 0; i < 7; i++) {
        pl[i] = tpSgn[i]*pb(tpRow[i]);
        for (j = 0; j < 7; j++)
            kl[i][j] = tpSgn[i]*tpSgn[j]*kb(tpRow[i], tpRow[j]);
    }

    double t[7][12];
    for (i = 0; i < 7; i++)
        for (j = 0; j < 12; j++)
            t[i][j] = T(i,j);

    // K = T' * kl * T
    double K[12][12];
    {
        double W[7][12];
        for (i = 0; i < 7; i++)
            for (j = 0; j < 12; j++) {
                double s = 0.0;
                for (k = 0; k < 7; k++)
                    s += kl[i][k]*t[k][j];
                W[i][j] = s;
            }
        for (i = 0; i < 12; i++)
            for (j = 0; j < 12; j++) {
                double s = 0.0;
                for (k = 0; k < 7; k++)
                    s += t[k][i]*W[k][j];
                K[i][j] = s;
            }
    }

    double m[6];
    for (i = 0; i < 6; i++)
        m[i] = pl[i]/(2*cos(ul(i)));

    double e1[3], e2[3], e3[3], rI1[3], rI2[3], rI3[3], rJ1[3], rJ2[3], rJ3[3];
    double a[3][3], lr2[12][3], lr3[12][3];
    for (k = 0; k < 3; k++) {
        e1[k] = e(k,0);   e2[k] = e(k,1);   e3[k] = e(k,2);
        rI1[k] = RI(k,0); rI2[k] = RI(k,1); rI3[k] = RI(k,2);
        rJ1[k] = RJ(k,0); rJ2[k] = RJ(k,1); rJ3[k] = RJ(k,2);
        for (j = 0; j < 3; j++)
            a[k][j] = A(k,j);
    }
    for (i = 0; i < 12; i++)
        for (j = 0; j < 3; j++) {
            lr2[i][j] = Lr2(i,j);
            lr3[i][j] = Lr3(i,j);
        }

    // ksigma1:  ks1_11 = A*pl(6);  ks1 = [ks1_11 o -ks1_11 o; ...]
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            const double v = a[i][j]*pl[6];
            K[i][j]     += v;   K[i][j+6]   -= v;
            K[i+6][j]   -= v;   K[i+6][j+6] += v;
        }

    // ksigma3:
    //  kbar2 = -Lr2*(m(3)*S(rI3) + m(1)*S(rI1)) + Lr3*(m(3)*S(rI2) - m(2)*S(rI1))
    //  kbar4 =  Lr2*(m(3)*S(rJ3) - m(4)*S(rJ1)) - Lr3*(m(3)*S(rJ2) + m(5)*S(rJ1))
    //  ks3 = [o kbar2 o kbar4] (+ transpose);  (Lr*S(w)) row i = Lr_i x w
    {
        double wa[3], wb[3], wc[3], wd[3];
        for (k = 0; k < 3; k++) {
            wa[k] = m[3]*rI3[k] + m[1]*rI1[k];
            wb[k] = m[3]*rI2[k] - m[2]*rI1[k];
            wc[k] = m[3]*rJ3[k] - m[4]*rJ1[k];
            wd[k] = m[3]*rJ2[k] + m[5]*rJ1[k];
        }
        for (i = 0; i < 12; i++) {
            const double *l2 = lr2[i], *l3 = lr3[i];
            double kb2[3], kb4[3];
            kb2[0] = -(l2[1]*wa[2] - l2[2]*wa[1]) + (l3[1]*wb[2] - l3[2]*wb[1]);
            kb2[1] = -(l2[2]*wa[0] - l2[0]*wa[2]) + (l3[2]*wb[0] - l3[0]*wb[2]);
            kb2[2] = -(l2[0]*wa[1] - l2[1]*wa[0]) + (l3[0]*wb[1] - l3[1]*wb[0]);
            kb4[0] =  (l2[1]*wc[2] - l2[2]*wc[1]) - (l3[1]*wd[2] - l3[2]*wd[1]);
            kb4[1] =  (l2[2]*wc[0] - l2[0]*wc[2]) - (l3[2]*wd[0] - l3[0]*wd[2]);
            kb4[2] =  (l2[0]*wc[1] - l2[1]*wc[0]) - (l3[0]*wd[1] - l3[1]*wd[0]);
            for (j = 0; j < 3; j++) {
                K[i][3+j] += kb2[j];   K[3+j][i] += kb2[j];
                K[i][9+j] += kb4[j];   K[9+j][i] += kb4[j];
            }
        }
    }

    // ksigma4:
    //  Ks4_22 =  m(3)*( S(e2)*S(rI3) - S(e3)*S(rI2)) + m(1)*(-S(e1)*S(rI2) + S(e2)*S(rI1))
    //          + m(2)*(-S(e1)*S(rI3) + S(e3)*S(rI1));
    //  Ks4_44 = -m(3)*( S(e2)*S(rJ3) - S(e3)*S(rJ2)) + m(4)*(-S(e1)*S(rJ2) + S(e2)*S(rJ1))
    //          + m(5)*(-S(e1)*S(rJ3) + S(e3)*S(rJ1));
    //  with S(x)*S(y) = y x' - (x.y) I
    {
        const double e2rI3 = e2[0]*rI3[0]+e2[1]*rI3[1]+e2[2]*rI3[2];
        const double e3rI2 = e3[0]*rI2[0]+e3[1]*rI2[1]+e3[2]*rI2[2];
        const double e2rI1 = e2[0]*rI1[0]+e2[1]*rI1[1]+e2[2]*rI1[2];
        const double e1rI2 = e1[0]*rI2[0]+e1[1]*rI2[1]+e1[2]*rI2[2];
        const double e3rI1 = e3[0]*rI1[0]+e3[1]*rI1[1]+e3[2]*rI1[2];
        const double e1rI3 = e1[0]*rI3[0]+e1[1]*rI3[1]+e1[2]*rI3[2];
        const double e2rJ3 = e2[0]*rJ3[0]+e2[1]*rJ3[1]+e2[2]*rJ3[2];
        const double e3rJ2 = e3[0]*rJ2[0]+e3[1]*rJ2[1]+e3[2]*rJ2[2];
        const double e2rJ1 = e2[0]*rJ1[0]+e2[1]*rJ1[1]+e2[2]*rJ1[2];
        const double e1rJ2 = e1[0]*rJ2[0]+e1[1]*rJ2[1]+e1[2]*rJ2[2];
        const double e3rJ1 = e3[0]*rJ1[0]+e3[1]*rJ1[1]+e3[2]*rJ1[2];
        const double e1rJ3 = e1[0]*rJ3[0]+e1[1]*rJ3[1]+e1[2]*rJ3[2];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                const double d = (i == j) ? 1.0 : 0.0;
                const double kI =
                     m[3]*((rI3[i]*e2[j] - e2rI3*d) - (rI2[i]*e3[j] - e3rI2*d))
                   + m[1]*((rI1[i]*e2[j] - e2rI1*d) - (rI2[i]*e1[j] - e1rI2*d))
                   + m[2]*((rI1[i]*e3[j] - e3rI1*d) - (rI3[i]*e1[j] - e1rI3*d));
                const double kJ =
                    -m[3]*((rJ3[i]*e2[j] - e2rJ3*d) - (rJ2[i]*e3[j] - e3rJ2*d))
                   + m[4]*((rJ1[i]*e2[j] - e2rJ1*d) - (rJ2[i]*e1[j] - e1rJ2*d))
                   + m[5]*((rJ1[i]*e3[j] - e3rJ1*d) - (rJ3[i]*e1[j] - e1rJ3*d));
                K[3+i][3+j] += kI;
                K[9+i][9+j] += kJ;
            }
    }

    // ksigma5:
    //  v = (1/Ln)*(m(2)*rI2 + m(3)*rI3 + m(5)*rJ2 + m(6)*rJ3)
    //  Ks5_11 = A*v*e1' + e1*v'*A + (e1'*v)*A
    //  Ks5_12 = -(m(2)*A*S(rI2) + m(3)*A*S(rI3)),  Ks5_14 = -(m(5)*A*S(rJ2) + m(6)*A*S(rJ3))
    {
        double v[3], Av[3];
        for (k = 0; k < 3; k++)
            v[k] = (m[1]*rI2[k] + m[2]*rI3[k] + m[4]*rJ2[k] + m[5]*rJ3[k])/Ln;
        const double e1tv = e1[0]*v[0] + e1[1]*v[1] + e1[2]*v[2];
        for (i = 0; i < 3; i++)
            Av[i] = a[i][0]*v[0] + a[i][1]*v[1] + a[i][2]*v[2];
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++) {
                // A is symmetric: e1*v'*A = e1*(A*v)'
                const double ks = e1tv*a[i][j] + Av[i]*e1[j] + e1[i]*Av[j];
                K[i][j]     += ks;   K[i][j+6]   -= ks;
                K[i+6][j]   -= ks;   K[i+6][j+6] += ks;
            }
        double wI[3], wJ[3];
        for (k = 0; k < 3; k++) {
            wI[k] = m[1]*rI2[k] + m[2]*rI3[k];
            wJ[k] = m[4]*rJ2[k] + m[5]*rJ3[k];
        }
        for (i = 0; i < 3; i++) {
            // (A*S(w)) row i = a_i x w
            double k12[3], k14[3];
            k12[0] = -(a[i][1]*wI[2] - a[i][2]*wI[1]);
            k12[1] = -(a[i][2]*wI[0] - a[i][0]*wI[2]);
            k12[2] = -(a[i][0]*wI[1] - a[i][1]*wI[0]);
            k14[0] = -(a[i][1]*wJ[2] - a[i][2]*wJ[1]);
            k14[1] = -(a[i][2]*wJ[0] - a[i][0]*wJ[2]);
            k14[2] = -(a[i][0]*wJ[1] - a[i][1]*wJ[0]);
            for (j = 0; j < 3; j++) {
                K[i][3+j]   += k12[j];  K[i+6][3+j] -= k12[j];
                K[3+j][i]   += k12[j];  K[3+j][i+6] -= k12[j];
                K[i][9+j]   += k14[j];  K[i+6][9+j] -= k14[j];
                K[9+j][i]   += k14[j];  K[9+j][i+6] -= k14[j];
            }
        }
    }

    // ksigma2
    {
        static Vector r2(3), r3(3), rIv1(3), rJv1(3), rm(3);
        for (k = 0; k < 3; k++) {
            r2(k)   = Rbar(k,1);
            r3(k)   = Rbar(k,2);
            rIv1(k) = rI1[k];
            rJv1(k) = rJ1[k];
        }
        for (k = 0; k < 3; k++) rm(k) = rI3[k] - rJ3[k];
        this->addKs2Matrix(r2, rm, m[3], K);
        for (k = 0; k < 3; k++) rm(k) = rJ2[k] - rI2[k];
        this->addKs2Matrix(r3, rm, m[3], K);
        this->addKs2Matrix(r2, rIv1, m[1], K);
        this->addKs2Matrix(r3, rIv1, m[2], K);
        this->addKs2Matrix(r2, rJv1, m[4], K);
        this->addKs2Matrix(r3, rJv1, m[5], K);
    }

    //  T * diag (M .* tan(thetal))*T'
    for (k = 0; k < 6; k++) {
        const double factor = pl[k] * tan(ul(k));
        if (factor == 0.0) continue;
        for (i = 0; i < 12; i++) {
            const double fi = t[k][i]*factor;
            if (fi == 0.0) continue;
            for (j = 0; j < 12; j++)
                K[i][j] += fi*t[k][j];
        }
    }

    for (i = 0; i < 12; i++)
        for (j = 0; j < 12; j++)
            kg(i,j) = K[i][j];

    return kg;
}


const Matrix &
CorotCrdTransf3d::getInitialGlobalStiffMatrix(const Matrix &kb)
{
    // transform tangent stiffness matrix from the basic system to local coordinates
    static Matrix kl(7,7);
    kl.addMatrixTripleProduct(0.0, Tp, kb, 1.0);      // kl = Tp ^ kb * Tp;
    
    // transform tangent  stiffness matrix from local to global coordinates
    //static Matrix kg(12,12);
    
    // compute the tangent stiffness matrix in global coordinates, using the
    // transformation of the undeformed configuration
    kg.addMatrixTripleProduct(0.0, T0, kl, 1.0);
    
    return kg;
}


int 
CorotCrdTransf3d::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
    // element projection
    
    static Vector dx(3);
    
    dx = (nodeJPtr->getCrds() + nodeJOffset) - (nodeIPtr->getCrds() + nodeIOffset);  
    if (nodeIInitialDisp != 0) {
        dx(0) -= nodeIInitialDisp[0];
        dx(1) -= nodeIInitialDisp[1];
        dx(2) -= nodeIInitialDisp[2];
    }
    
    if (nodeJInitialDisp != 0) {
        dx(0) += nodeJInitialDisp[0];
        dx(1) += nodeJInitialDisp[1];
        dx(2) += nodeJInitialDisp[2];
    }
    
    // calculate the element length
    
    L = dx.Norm();
    
    if (L == 0.0) 
    {
      opserr << "\nCorotCrdTransf3d::getLocalAxes transfTag = " << this->getTag();
      opserr << "\nelement has zero length" << endln;
      return -2;  
    }
    
    // calculate the element local x axis components (direction cossines)
    // wrt to the global coordinates 
    xAxis = dx/L;
    
    XAxis(0) = xAxis(0);    XAxis(1) = xAxis(1);    XAxis(2) = xAxis(2);
    
    // calculate the cross-product y = v * x   
    static Vector yAxis(3), zAxis(3);
    
    yAxis(0) = vAxis(1)*xAxis(2) - vAxis(2)*xAxis(1);
    yAxis(1) = vAxis(2)*xAxis(0) - vAxis(0)*xAxis(2);
    yAxis(2) = vAxis(0)*xAxis(1) - vAxis(1)*xAxis(0);
    
    double ynorm = yAxis.Norm();
    
    if (ynorm == 0)
    {
      opserr << "\nCorotCrdTransf3d::getLocalAxes transfTag = " << this->getTag();
      opserr << "\nvector v that defines plane xz is parallel to x axis" << endln;
      return -3;
    }
    
    yAxis /= ynorm;
    YAxis(0) = yAxis(0);    YAxis(1) = yAxis(1);    YAxis(2) = yAxis(2);
    
    // calculate the cross-product z = x * y 
    
    zAxis(0) = xAxis(1)*yAxis(2) - xAxis(2)*yAxis(1);
    zAxis(1) = xAxis(2)*yAxis(0) - xAxis(0)*yAxis(2);
    zAxis(2) = xAxis(0)*yAxis(1) - xAxis(1)*yAxis(0);
    ZAxis(0) = zAxis(0);    ZAxis(1) = zAxis(1);    ZAxis(2) = zAxis(2);
    
    for (int i=0; i < 3; i++) {
        R0(i,0) = xAxis(i);
        R0(i,1) = yAxis(i);
        R0(i,2) = zAxis(i);
    }
    
    return 0;
}

int
CorotCrdTransf3d::getRigidOffsets(Vector &offsets)
{
  offsets(0) = nodeIOffset(0);
  offsets(1) = nodeIOffset(1);
  offsets(2) = nodeIOffset(2);
  offsets(3) = nodeJOffset(0);
  offsets(4) = nodeJOffset(1);
  offsets(5) = nodeJOffset(2);

  return 0;
}

double 
CorotCrdTransf3d::getInitialLength(void)
{
    return L;
}


double 
CorotCrdTransf3d::getDeformedLength(void)
{
    return Ln;
}


const Vector &
CorotCrdTransf3d::getQuaternionFromRotMatrix(const Matrix &R) const
{
    // obtains the normalised quaternion from the rotation matrix
    int i, j, k;
    double trR;              // trace of R
    double a    ;
    static Vector q(4);      // normalized quaternion
    
    trR = R(0,0) + R(1,1) + R(2,2);    
    
    // a = max ([trR R(0,0) R(1,1) R(2,2)]);
    a = trR;
    for (i = 0; i < 3; i++)
        if (R(i,i) > a)
            a = R(i,i);
        
        if (a == trR)
        {
            q(3) = sqrt(1+a)*0.5;
            
            for (i = 0; i < 3; i++)
            {
                j = (i+1)%3;
                k = (i+2)%3;
                q(i) = (R(k,j) - R(j,k))/(4*q(3));
            }
        }  
        else
        {
            for (i = 0; i < 3; i++)
                if (a == R(i,i))
                {
                    j = (i+1)%3;
                    k = (i+2)%3;
                    
                    q(i) = sqrt(a*0.5 + (1 - trR)/4.0);
                    q(3)   = (R(k,j) - R(j,k))/(4*q(i));
                    q(j) = (R(j,i) + R(i,j))/(4*q(i));
                    q(k) = (R(k,i) + R(i,k))/(4*q(i));
                }
        }
        
        return q;
}


const Vector &
CorotCrdTransf3d::getQuaternionFromPseudoRotVector(const Vector  &theta) const
{
    double t;                // norm of the pseudo rotation vector
    double factor;
    static Vector q(4);      // normalized quaternion
    
    t = theta.Norm();
    
    if (t == 0)
        q.Zero();
    else
    {
        factor = sin(t*0.5)/ t;
        for (int i = 0; i < 3; i++)
            q(i) = theta(i) * factor;
    }
    
    q(3) = cos(t*0.5);
    
    return q;
}


const Vector &
CorotCrdTransf3d::quaternionProduct(const Vector &q1, const Vector &q2) const
{
    
    static Vector q12(4);
    int i;
    double q1Tq2= 0;  // dot product
    static Vector q1xq2(3);     // cross product
    
    // calculate the dot product q1.q2
    for (i = 0; i < 3; i++)       // NOTE i <3, not i<4
        q1Tq2 += q1(i)* q2(i);
    
    // calculate the cross-product q1 x q2
    q1xq2(0) = q1(1)*q2(2) - q1(2)*q2(1);
    q1xq2(1) = q1(2)*q2(0) - q1(0)*q2(2);
    q1xq2(2) = q1(0)*q2(1) - q1(1)*q2(0);
    
    // calculate the quaternion product
    for (i = 0; i < 3; i++)
        q12(i) = q1(3)*q2(i) + q2(3)*q1(i) - q1xq2(i);  
    
    q12(3) = q1(3)*q2(3) - q1Tq2;
    
    return q12;
}


const Matrix &
CorotCrdTransf3d::getRotationMatrixFromQuaternion(const Vector &q) const
{ 
    int i, j;
    double factor;
    static Matrix I(3,3); // identity matrix
    static Matrix qqT(3,3); 
    static Matrix S(3,3);
    static Matrix R(3,3);
    
    // R = (q0^2 - q' * q) * I + 2 * q * q' + 2*q0*S(q);
    
    factor = q(3)*q(3) - (q(0)*q(0) + q(1)*q(1) + q(2)*q(2));
    
    
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            qqT(i,j) = q(i) * q(j);
        
        // get skew symmetric matrix	     
        S = this->getSkewSymMatrix (q);
        
        R.Zero();
        
        for (i = 0; i < 3; i++)
            R(i,i) = factor;
        
        R.addMatrix (1.0, qqT, 2.0);
        R.addMatrix (1.0, S, 2.0*q(3));
        
        return R;
}


const Vector &
CorotCrdTransf3d::getTangScaledPseudoVectorFromQuaternion(const Vector &q) const
{ 
    static Vector w(3);
    
    for (int i = 0; i < 3; i++)
        w(i) = 2.0 * q(i)/q(3);
    
    return w;
}


const Matrix &
CorotCrdTransf3d::getRotMatrixFromTangScaledPseudoVector(const Vector &w) const
{ 
    // Rotation matrix in terms of the tangent-scaled pseudo-vector
    static Matrix S(3,3);
    static Matrix S2(3,3);
    static Matrix R(3,3);
    double normw2;
    
    S = this->getSkewSymMatrix(w);
    
    // R = I + (S + S*S/2)/(1 + w' * w / 4);
    
    R.Zero();
    
    for (int i = 0; i < 3; i++)
        R(i,i) = 1;
    
    S2 = S;
    S2.addMatrixProduct(1.0, S, S, 0.5);
    
    normw2 = w^ w;
    
    R.addMatrix(1.0, S2, 1.0/(1 + normw2/4.0));
    
    return R;
}


const Matrix &
CorotCrdTransf3d::getSkewSymMatrix(const Vector &theta) const
{
    static Matrix S(3,3);
    
    //  St = [   0       -theta(2)  theta(1);
    //         theta(2)     0      -theta(0);
    //        -theta(1)   theta(0)      0   ];
    
    S(0,0) =  0.0;
    S(0,1) = -theta(2);
    S(0,2) =  theta(1);
    
    S(1,0) =  theta(2);
    S(1,1) =  0.0;
    S(1,2) = -theta(0);
    
    S(2,0) = -theta(1);
    S(2,1) =  theta(0);
    S(2,2) =  0.0;
    
    //opserr << "getSkew  theta: " << theta;
    //opserr << "getSkew  S: " << S;
    
    return S;
}   


const Matrix &
CorotCrdTransf3d::getLMatrix(const Vector &ri) const
{
    static Matrix L1(3,3), L2(3,3);
    static Vector r1(3), e1(3);
    double rie1, e1r1k;
    static Matrix rie1r1(3,3);
    static Matrix e1e1r1(3,3);
    static Matrix Sri(3,3);
    static Matrix Sr1(3,3);
    static Matrix L(12,3);
    
    int j, k;
    
    for (j = 0; j < 3; j++)
    {
        e1(j) =    e(j,0);
        r1(j) = Rbar(j,0);
    }  
    
    rie1 = ri^ e1;
    
    //opserr << "rie1: " << rie1 << endln;
    
    for (k = 0; k < 3; k++)
    {
        e1r1k = (e1(k) + r1(k));
        for (j = 0; j < 3; j++)
        {
            rie1r1(j,k) = ri(j)*e1r1k;
            e1e1r1(j,k) = e1(j)*e1r1k;
        }
    }
    
    //opserr << "rie1r1: " << rie1r1;
    
    //opserr << "e1e1r1: " << e1e1r1;
    
    //L1  = ri'*e1 * A/2 + A*ri*(e1 + r1)'/2;
    
    L1.addMatrix (0.0, A, rie1*0.5);
    L1.addMatrixProduct (1.0, A, rie1r1, 0.5);
    
    //opserr << "L1: " << L1;
    
    Sri = this->getSkewSymMatrix(ri);       // Change ??????????????????????
    Sr1 = this->getSkewSymMatrix(r1);       // Change ??????????????????????
    
    // L2  = Sri/2 - ri'*e1*S(r1)/4 - Sri*e1*(e1 + r1)'/4;
    L2.addMatrix (0.0, Sri, 0.5);
    L2.addMatrix (1.0, Sr1, -rie1/4.0);
    L2.addMatrixProduct (1.0, Sri, e1e1r1, -0.25);
    
    //opserr << "L2: " << L2;
    
    // L = [L1
    //      L2 
    //     -L1 
    //      L2 ;
    
    L.Zero();
    L.Assemble(L1, 0, 0,  1.0);
    L.Assemble(L2, 3, 0,  1.0);
    L.Assemble(L1, 6, 0, -1.0);
    L.Assemble(L2, 9, 0,  1.0);
    
    //opserr << "L: " << L;
    
    return L;
}


// Adds factor * Ksigma2(ri, z) directly into kg (Crisfield's Ks2 term).
//
//  Ksigma2 = [ K11   K12 -K11   K12;
//              K12t  K22 -K12t  K22;
//             -K11  -K12  K11  -K12;
//              K12t  K22 -K12t  K22];
//
//  U   = (-1/2)*A*z*ri'*A + ri'*e1*A*z*e1'/(2*Ln) + z'*(e1+r1)*A*ri*e1'/(2*Ln)
//  K11 = U + U' + ri'*e1*(2*(e1'*z)+z'*r1)*A/(2*Ln)
//  K12 = (1/4)*(-A*z*e1'*Sri - A*ri*z'*Sr1 - z'*(e1+r1)*A*Sri)
//  K22 = (1/8)*((-ri'*e1)*Sz*Sr1 + Sr1*z*e1'*Sri + Sri*e1*z'*Sr1
//               - (e1+r1)'*z*S(e1)*Sri + 2*Sz*Sri)
//
// evaluated with 3-vectors using  a'*S(b) = (a x b)',  S(a)*b = a x b  and
// S(a)*S(b) = b*a' - (a.b)*I  (same expressions as the original
// getKs2Matrix, without 12x12 temporaries).
void
CorotCrdTransf3d::addKs2Matrix(const Vector &riV, const Vector &zV, double factor,
                               double K[12][12]) const
{
    if (factor == 0.0)
        return;

    double ri[3], z[3], e1[3], r1[3];
    for (int i = 0; i < 3; i++) {
        ri[i] = riV(i);
        z[i]  = zV(i);
        e1[i] = e(i,0);
        r1[i] = Rbar(i,0);
    }
    double a[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            a[i][j] = A(i,j);

    const double rite1 = ri[0]*e1[0] + ri[1]*e1[1] + ri[2]*e1[2];
    const double zte1  = z[0]*e1[0]  + z[1]*e1[1]  + z[2]*e1[2];
    const double ztr1  = z[0]*r1[0]  + z[1]*r1[1]  + z[2]*r1[2];
    const double s     = zte1 + ztr1;
    const double rz    = ri[0]*z[0]  + ri[1]*z[1]  + ri[2]*z[2];
    const double r1ri  = r1[0]*ri[0] + r1[1]*ri[1] + r1[2]*ri[2];
    const double e1ri  = rite1;
    const double zr1   = ztr1;

    // A*z, A*ri
    double Az[3], Ari[3];
    for (int i = 0; i < 3; i++) {
        Az[i]  = a[i][0]*z[0]  + a[i][1]*z[1]  + a[i][2]*z[2];
        Ari[i] = a[i][0]*ri[0] + a[i][1]*ri[1] + a[i][2]*ri[2];
    }
    // cross products
    double e1xri[3], zxr1[3], r1xz[3], rixe1[3];
    e1xri[0] = e1[1]*ri[2] - e1[2]*ri[1];  e1xri[1] = e1[2]*ri[0] - e1[0]*ri[2];  e1xri[2] = e1[0]*ri[1] - e1[1]*ri[0];
    zxr1[0]  = z[1]*r1[2]  - z[2]*r1[1];   zxr1[1]  = z[2]*r1[0]  - z[0]*r1[2];   zxr1[2]  = z[0]*r1[1]  - z[1]*r1[0];
    for (int i = 0; i < 3; i++) {
        r1xz[i]  = -zxr1[i];
        rixe1[i] = -e1xri[i];
    }

    // K11
    double U[3][3], K11[3][3];
    const double c1 = rite1/(2.0*Ln), c2 = s/(2.0*Ln);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            U[i][j] = -0.5*Az[i]*Ari[j] + c1*Az[i]*e1[j] + c2*Ari[i]*e1[j];
    const double c3 = rite1*(2.0*zte1 + ztr1)/(2.0*Ln);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            K11[i][j] = U[i][j] + U[j][i] + c3*a[i][j];

    // K12 = 1/4 (-(Az)(e1 x ri)' - (A ri)(z x r1)' - s A S(ri))
    double K12[3][3];
    for (int i = 0; i < 3; i++) {
        // row i of A*S(ri) = (A row i) x ... : (A S(ri))_ij = sum_k a_ik S(ri)_kj
        double ASri[3];
        ASri[0] =  a[i][1]*ri[2] - a[i][2]*ri[1];
        ASri[1] = -a[i][0]*ri[2] + a[i][2]*ri[0];
        ASri[2] =  a[i][0]*ri[1] - a[i][1]*ri[0];
        for (int j = 0; j < 3; j++)
            K12[i][j] = 0.25*(-Az[i]*e1xri[j] - Ari[i]*zxr1[j] - s*ASri[j]);
    }

    // K22 = 1/8 ( -rite1 S(z)S(r1) + (r1 x z)(e1 x ri)' + (ri x e1)(z x r1)'
    //             - s S(e1)S(ri) + 2 S(z)S(ri) ),  S(a)S(b) = b a' - (a.b) I
    double K22[3][3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            const double d = (i == j) ? 1.0 : 0.0;
            const double SzSr1 = r1[i]*z[j]  - zr1*d;
            const double Se1Sri = ri[i]*e1[j] - e1ri*d;
            const double SzSri = ri[i]*z[j]  - rz*d;
            K22[i][j] = 0.125*(-rite1*SzSr1 + r1xz[i]*e1xri[j] + rixe1[i]*zxr1[j]
                               - s*Se1Sri + 2.0*SzSri);
        }
    (void)r1ri;

    // assemble factor * Ksigma2 into K
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            const double k11 = factor*K11[i][j];
            const double k12 = factor*K12[i][j];
            const double k21 = factor*K12[j][i];
            const double k22 = factor*K22[i][j];
            K[i][j] += k11;  K[i][j+6] -= k11;
            K[i+6][j] -= k11;  K[i+6][j+6] += k11;
            K[i][j+3] += k12;  K[i][j+9] += k12;
            K[i+6][j+3] -= k12;  K[i+6][j+9] -= k12;
            K[i+3][j] += k21;  K[i+9][j] += k21;
            K[i+3][j+6] -= k21;  K[i+9][j+6] -= k21;
            K[i+3][j+3] += k22;  K[i+3][j+9] += k22;
            K[i+9][j+3] += k22;  K[i+9][j+9] += k22;
        }
}


CrdTransf *
CorotCrdTransf3d::getCopy3d(void)
{
    // create a new instance of CorotCrdTransf3d 
    
    CorotCrdTransf3d *theCopy = new CorotCrdTransf3d (this->getTag(), vAxis, nodeIOffset, nodeJOffset);
    
    if (!theCopy)
    {
        opserr << "CorotCrdTransf3d::getCopy() - out of memory creating copy\n";
        return 0;
    }    
    
    theCopy->nodeIPtr = nodeIPtr;
    theCopy->nodeJPtr = nodeJPtr;
    theCopy->xAxis = xAxis;
    theCopy->L = L;
    theCopy->Ln = Ln;
    theCopy->R0 = R0;
    theCopy->alphaIq = alphaIq;
    theCopy->alphaJq = alphaJq;
    theCopy->alphaIqcommit = alphaIqcommit;
    theCopy->alphaJqcommit = alphaJqcommit;
    theCopy->ul = ul;
    theCopy->ulcommit = ulcommit;
    theCopy->ulpr = ulpr;

    // last trial rotations: without these the copy would treat the full
    // current nodal rotation as a new increment on its next update()
    theCopy->alphaI = alphaI;
    theCopy->alphaJ = alphaJ;

    // initial displacements
    if (nodeIInitialDisp != 0) {
        theCopy->nodeIInitialDisp = new double [6];
        for (int i=0; i<6; i++)
            theCopy->nodeIInitialDisp[i] = nodeIInitialDisp[i];
    }
    if (nodeJInitialDisp != 0) {
        theCopy->nodeJInitialDisp = new double [6];
        for (int i=0; i<6; i++)
            theCopy->nodeJInitialDisp[i] = nodeJInitialDisp[i];
    }
    theCopy->initialDispChecked = initialDispChecked;

    return theCopy;
}


int 
CorotCrdTransf3d::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(48);
  for (int i=0; i<7; i++) 
    data(i) = ulcommit(i);
  for (int j=0; j<4; j++) {
    data(7+j) = alphaIqcommit(j);
    data(11+j) = alphaJqcommit(j);
  }
  
  for (int k=0; k<3; k++) {
    data(15+k) = xAxis(k);
    data(18+k) = vAxis(k);
    data(21+k) = nodeIOffset(k);
    data(24+k) = nodeJOffset(k);
    data(27+k) = alphaI(k);
    data(30+k) = alphaJ(k);
  }
  
  if (nodeIInitialDisp != 0) {
    data(34) = nodeIInitialDisp[0];
    data(35) = nodeIInitialDisp[1];
    data(36) = nodeIInitialDisp[2];
    data(37) = nodeIInitialDisp[3];
    data(38) = nodeIInitialDisp[4];
    data(39) = nodeIInitialDisp[5];
  } else {
    data(34)  = 0.0;
    data(35)  = 0.0;
    data(36) = 0.0;
    data(37) = 0.0;
    data(38) = 0.0;
    data(39) = 0.0;
  }
  
  if (nodeJInitialDisp != 0) {
    data(40) = nodeJInitialDisp[0];
    data(41) = nodeJInitialDisp[1];
    data(42) = nodeJInitialDisp[2];
    data(43) = nodeJInitialDisp[3];
    data(44) = nodeJInitialDisp[4];
    data(45) = nodeJInitialDisp[5];
  } else {
    data(40) = 0.0;
    data(41) = 0.0;
    data(42) = 0.0;
    data(43) = 0.0;
    data(44) = 0.0;
    data(45) = 0.0;
  }
  data(46) = L;
  data(47) = Ln;
  
  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << " CorotCrdTransf3d::sendSelf() - data could not be sent\n" ;
    return -1;
  }


  return 0;
}


int 
CorotCrdTransf3d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(48);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << " CorotCrdTransf3d::recvSelf() - data could not be received\n" ;
    return -1;
    }
  int i,j;
  for (i=0; i<7; i++) 
    ulcommit(i) = data(i);
  for (j=0; j<4; j++) {
    alphaIqcommit(j) = data(7+j);
    alphaJqcommit(j) = data(11+j);
  }
  
  for (int k=0; k<3; k++) {
    xAxis(k) = data(15+k);
    vAxis(k) = data(18+k); 
    nodeIOffset(k) = data(21+k);
    nodeJOffset(k) = data(24+k);
    alphaI(k) = data(27+k);
    alphaJ(k) = data(30+k);
  }
  
  int flag;
  flag = 0;
  for (i=34; i<=39; i++)
    if (data(i) != 0.0)
      flag = 1;
  if (flag == 1) {
    if (nodeIInitialDisp == 0)
      nodeIInitialDisp = new double[6];
    for (i=34, j=0; i<=39; i++, j++)
      nodeIInitialDisp[j] = data(i);
  }
  
  flag = 0;
  for (i=40; i<=45; i++)
    if (data(i) != 0.0)
      flag = 1;
  if (flag == 1) {
    if (nodeJInitialDisp == 0)
      nodeJInitialDisp = new double[6];
    for (i=40, j=0; i<=45; i++, j++)
      nodeJInitialDisp[j] = data(i);
  }
  L = data(46);
  Ln = data(47);
  
  ul = ulcommit;
  alphaIq = alphaIqcommit;
  alphaJq = alphaJqcommit;
  
  initialDispChecked = true;
  return 0;  
}


const Matrix &
CorotCrdTransf3d::getGlobalMatrixFromLocal(const Matrix &local)
{
    this->compTransfMatrixLocalGlobal(Tlg);  // OPTIMIZE LATER
    kg.addMatrixTripleProduct(0.0, Tlg, local, 1.0);  // OPTIMIZE LATER

    return kg;
}


const Vector &
CorotCrdTransf3d::getPointGlobalCoordFromLocal(const Vector &xl)
{
    static Vector xg(3);
    opserr << " CorotCrdTransf3d::getPointGlobalCoordFromLocal: not implemented yet" ;
    
    return xg;  
}


const Vector &
CorotCrdTransf3d::getPointGlobalDisplFromBasic(double xi, const Vector &uxb)
{
    static Vector uxg(3);
    opserr << " CorotCrdTransf3d::getPointGlobalDisplFromBasic: not implemented yet" ;
    
    
    return uxg;  
}


const Vector &
CorotCrdTransf3d::getPointLocalDisplFromBasic(double xi, const Vector &uxb)
{
    static Vector uxg(3);
    opserr << " CorotCrdTransf3d::getPointLocalDisplFromBasic: not implemented yet" ;
    
    
    return uxg;  
}


void
CorotCrdTransf3d::Print(OPS_Stream &s, int flag)
{
	if (flag == OPS_PRINT_CURRENTSTATE) {
		s << "\nCrdTransf: " << this->getTag() << " Type: CorotCrdTransf3d";
		s << "\tvAxis: " << vAxis;
		s << "\tnodeI Offset: " << nodeIOffset;
		s << "\tnodeJ Offset: " << nodeJOffset;
	}

	if (flag == OPS_PRINT_PRINTMODEL_JSON) {
		s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"CorotCrdTransf3d\"";
        s << ", \"vecInLocXZPlane\": [" << vAxis(0) << ", " << vAxis(1) << ", " << vAxis(2) << "]";
        if (nodeIOffset != 0)
            s << ", \"iOffset\": [" << nodeIOffset[0] << ", " << nodeIOffset[1] << ", " << nodeIOffset[2] << "]";
        if (nodeJOffset != 0)
            s << ", \"jOffset\": [" << nodeJOffset[0] << ", " << nodeJOffset[1] << ", " << nodeJOffset[2] << "]";
        s << "}";
    }
}
