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
Matrix CorotCrdTransf3d::RI(3,3); 
Matrix CorotCrdTransf3d::RJ(3,3); 
Matrix CorotCrdTransf3d::Rbar(3,3); 
Matrix CorotCrdTransf3d::e(3,3); 
Matrix CorotCrdTransf3d::Tp(6,7); 
Matrix CorotCrdTransf3d::T(7,12);
Matrix CorotCrdTransf3d::Tlg(12,12);
Matrix CorotCrdTransf3d::TlgInv(12, 12);
Matrix CorotCrdTransf3d::Tbl(6,12);
Matrix CorotCrdTransf3d::kg(12,12);
Matrix CorotCrdTransf3d::Lr2(12,3);
Matrix CorotCrdTransf3d::Lr3(12,3);
Matrix CorotCrdTransf3d::A(3,3);

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
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    // check vector that defines local xz plane
    if (&vecInLocXZPlane == 0 || vecInLocXZPlane.Size() != 3 )
    {
        opserr << "CorotCrdTransf3d::CorotCrdTransf3d:  Vector that defines local xz plane is invalid\n";
        opserr << "Size must be 3\n. Using (0,0,1)";      
        vAxis(0) = 0;       vAxis(1) = 0;      vAxis(2) = 1;
    }
    else
        vAxis = vecInLocXZPlane;
    
    // check rigid joint offset for node I
    if (&rigJntOffsetI == 0 || rigJntOffsetI.Size() != 3 )
    {
        opserr << "CorotCrdTransf3d::CorotCrdTransf3d:  Invalid rigid joint offset vector for node I\n";
        opserr << "Size must be 3\n";      
        nodeIOffset.Zero();      
    }
    else
        nodeIOffset = rigJntOffsetI;
    
    // check rigid joint offset for node J
    if (&rigJntOffsetJ == 0 || rigJntOffsetJ.Size() != 3 )
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
nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
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
    
    // compute initial pseudo-vectors for nodal triads
    //opserr << setiosflags(ios::scientific);
    //opserr << setiosflags(ios::showpos);
    //opserr << setprecision(16);
    
    //opserr << "L: " << L;
    //opserr << "R0: " << R0;
    
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
        
        Lr2 = this->getLMatrix (r2);
        Lr3 = this->getLMatrix (r3);
        
        static Matrix Sr1(3,3), Sr2(3,3), Sr3(3,3);
        static Vector Se(3), At(3);
        
        //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';
        //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';
        //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';
        //  
        //   T4 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
        //   T5 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
        //   T6 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']';
        
        T.Zero();
        
        Sr1 = this->getSkewSymMatrix(rI1);
        Sr2 = this->getSkewSymMatrix(rI2);
        Sr3 = this->getSkewSymMatrix(rI3);
        
        //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';
        
        Se.addMatrixVector(0.0, Sr3, e2, -1.0);     // (-S(rI3)*e2 + S(rI2)*e3)
        Se.addMatrixVector(1.0, Sr2, e3,  1.0);    
        
        for (i = 0; i < 3; i++)
            T(0,i+3) =  Se(i);
        
        //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';
        
        At.addMatrixVector(0.0, A, rI2, 1.0);   
        
        Se.addMatrixVector(0.0, Sr2, e1, -1.0);     // (-S(rI2)*e1 + S(rI1)*e2)'
        Se.addMatrixVector(1.0, Sr1, e2,  1.0);
        
        for (i = 0; i < 3; i++)
        {
            T(1,i  ) =  At(i);  
            T(1,i+3) =  Se(i);
            T(1,i+6) = -At(i);
        }
        
        //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';
        
        At.addMatrixVector(0.0, A, rI3, 1.0);   
        
        Se.addMatrixVector(0.0, Sr3, e1, -1.0);     // (-S(rI3)*e1 + S(rI1)*e3)
        Se.addMatrixVector(1.0, Sr1, e3,  1.0);
        
        for (i = 0; i < 3; i++)
        {
            T(2,i  ) =  At(i);  
            T(2,i+3) =  Se(i);
            T(2,i+6) = -At(i);
        }
        
        Sr1 = this->getSkewSymMatrix(rJ1);
        Sr2 = this->getSkewSymMatrix(rJ2);
        Sr3 = this->getSkewSymMatrix(rJ3);
        
        //   T4 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
        
        Se.addMatrixVector(0.0, Sr3, e2, -1.0);    // -S(rJ3)*e2 + S(rJ2)*e3
        Se.addMatrixVector(1.0, Sr2, e3,  1.0);    
        
        for (i = 0; i < 3; i++)
            T(3,i+9) =  Se(i);
        
        //   T5 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
        
        At.addMatrixVector(0.0, A, rJ2, 1.0);   
        
        Se.addMatrixVector(0.0, Sr2, e1, -1.0);     // (-S(rJ2)*e1 + S(rJ1)*e2)
        Se.addMatrixVector(1.0, Sr1, e2,  1.0);
        
        for (i = 0; i < 3; i++)
        {
            T(4,i  ) =  At(i);  
            T(4,i+6) = -At(i);
            T(4,i+9) =  Se(i);
        }
        
        //   T6 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']';
        
        At.addMatrixVector(0.0, A, rJ3, 1.0);   
        
        Se.addMatrixVector(0.0, Sr3, e1, -1.0);     // (-S(rJ3)*e1 + S(rJ1)*e3)
        Se.addMatrixVector(1.0, Sr1, e3,  1.0);
        
        for (i = 0; i < 3; i++)
        {
            T(5,i  ) =  At(i);  
            T(5,i+6) = -At(i);
            T(5,i+9) =  Se(i);
        }
        
        // setup tranformation matrix
        static Vector Lr(12);
        
        // T(:,1) += Lr3*rI2 - Lr2*rI3;
        // T(:,2) +=           Lr2*rI1;
        // T(:,3) += Lr3*rI1          ;
        
        // T(:,4) += Lr3*rJ2 - Lr2*rJ3;
        // T(:,5) += Lr2*rJ1          ;      // ?????? check sign
        // T(:,6) += Lr3*rJ1          ;      // ?????? check sign
        
        Lr.addMatrixVector(0.0, Lr3, rI2,  1.0);  //  T(:,1) += Lr3*rI2 - Lr2*rI3 
        Lr.addMatrixVector(1.0, Lr2, rI3, -1.0);
        
        for (i = 0; i < 12; i++)
            T(0,i) += Lr(i);
        
        Lr.addMatrixVector(0.0, Lr2, rI1,  1.0);  //  T(:,2) +=           Lr2*rI1
        
        for (i = 0; i < 12; i++)
            T(1,i) += Lr(i);
        
        Lr.addMatrixVector(0.0, Lr3, rI1,  1.0);  //  T(:,3) += Lr3*rI1 
        
        for (i = 0; i < 12; i++)
            T(2,i) += Lr(i);
        
        Lr.addMatrixVector(0.0, Lr3, rJ2,  1.0);  //  T(:,4) += Lr3*rJ2 - Lr2*rJ3;
        Lr.addMatrixVector(1.0, Lr2, rJ3, -1.0);
        
        for (i = 0; i < 12; i++)
            T(3,i) += Lr(i);
        
        Lr.addMatrixVector(0.0, Lr2, rJ1,  1.0);  //  T(:,5) += Lr2*rJ1   
        
        for (i = 0; i < 12; i++)
            T(4,i) += Lr(i);
        
        Lr.addMatrixVector(0.0, Lr3, rJ1,  1.0);  //   T(:,6) += Lr3*rJ1
        
        for (i = 0; i < 12; i++)
            T(5,i) += Lr(i);
        
        double c;
        for (j = 0; j < 6; j++)
        {
            c = 2 * cos(ul(j));
            
            for (i = 0; i < 12; i++) 
                T(j,i) /= c;
        }
        
        // T(:,7) = [-e1' O' e1' O']';
        for (i = 0; i < 3; i++)
        {
            T(6,i  ) = -e1(i);
            T(6,i+6) =  e1(i);
        }      
}


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
        
        
        opserr << "thetaI: " << thetaI;
        opserr << "thetaJ: " << thetaJ;
        
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
        c = 1.0/(2.0*cos(thetaI(0)));   
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
    ub.addMatrixVector(0.0, Tp, ul, 1.0);
    
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
    dub.addMatrixVector(0.0, Tp, dul, 1.0);
    
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
    Dub.addMatrixVector(0.0, Tp, Dul, 1.0);
    
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
    this->update();
    
    static Vector pg(12);
    pg.Zero();
    
    // if there are no element loads present
    if (p0 == 0.0) {
        // transform resisting forces from the basic system to local coordinates
        static Vector pl(7);
        pl.addMatrixTransposeVector(0.0, Tp, pb, 1.0);    // pl = Tp ^ pb;

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
        pl.addMatrixTransposeVector(0.0, Tp, pb, 1.0);    // pl = Tp ^ pb;

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
    this->update();
    
    int i, j, k;   
    // transform tangent stiffness matrix from the basic system to local coordinates
    static Matrix kl(7,7);
    kl.addMatrixTripleProduct(0.0, Tp, kb, 1.0);      // kl = Tp ^ kb * Tp;

    //    opserr << "kb: " << kb;
    //    opserr << "Tp: " << Tp;
    
    // transform resisting forces from the basic system to local coordinates
    static Vector pl(7);
    pl.addMatrixTransposeVector(0.0, Tp, pb, 1.0);    // pl = Tp ^ pb;
    
    // transform tangent  stiffness matrix from local to global coordinates
    //static Matrix kg(12,12);
    
    // compute the tangent stiffness matrix in global coordinates
    kg.addMatrixTripleProduct(0.0, T, kl, 1.0);
    
    static Vector m(6);
    for (i = 0; i < 6; i++)
        m(i) = pl(i)/(2*cos(ul(i)));
    
    // compute the basic rotations
    
    static Vector e1(3), e2(3), e3(3);
    static Vector r1(3), r2(3), r3(3);
    static Vector rI1(3), rI2(3), rI3(3);
    static Vector rJ1(3), rJ2(3), rJ3(3);
    
    for (k = 0; k < 3; k ++)
    {
        e1(k) = e(k,0);
        e2(k) = e(k,1);
        e3(k) = e(k,2);
        
        r1(k) = Rbar(k,0);
        r2(k) = Rbar(k,1);
        r3(k) = Rbar(k,2);
        
        rI1(k) = RI(k,0);
        rI2(k) = RI(k,1);
        rI3(k) = RI(k,2);
        
        rJ1(k) = RJ(k,0);
        rJ2(k) = RJ(k,1);
        rJ3(k) = RJ(k,2);
    }
    
    //   ks = t'*kl*t + ks1 + t * diag (m .* tan(thetal))*t' + ...
    //        m(4)*(ks2r2t3_u3 + ks2r3u2_t2) + ...
    //        m(2)*ks2r2t1 + m(3)*ks2r3t1 + ...
    //        m(5)*ks2r2u1 + m(6)*ks2r3u1 + ...
    //        ks3 + ks3' + ks4 + ks5;
    
    static Matrix Se1(3,3), Se2(3,3), Se3(3,3);
    static Matrix SrI1(3,3), SrI2(3,3), SrI3(3,3);
    static Matrix SrJ1(3,3), SrJ2(3,3), SrJ3(3,3);
    
    Se1 = this->getSkewSymMatrix(e1);
    Se2 = this->getSkewSymMatrix(e2);
    Se3 = this->getSkewSymMatrix(e3);
    
    SrJ1 = this->getSkewSymMatrix(rJ1);
    SrJ2 = this->getSkewSymMatrix(rJ2);
    SrJ3 = this->getSkewSymMatrix(rJ3);
    
    SrI1 = this->getSkewSymMatrix(rI1);
    SrI2 = this->getSkewSymMatrix(rI2);
    SrI3 = this->getSkewSymMatrix(rI3);
    
    
    // ksigma1 -------------------------------
    //   ks1_11 =  a*pl(6);
    //   ks1 = [ ks1_11  o  -ks1_11  o;
    //             o     o      o    o;
    //          -ks1_11  o   ks1_11  o;
    //             o     o      o    o];
    
    double factor;
    
    kg.Assemble(A, 0, 0,  pl(6));
    kg.Assemble(A, 0, 6, -pl(6));
    kg.Assemble(A, 6, 0, -pl(6));
    kg.Assemble(A, 6, 6,  pl(6));
    
    //opserr << "kg += ksigma1: " << kg;
    
    // ksigma3 -------------------------------
    //  kbar2 = -Lr2*(m(3)*S(rI3) + m(1)*S(rI1)) + ...
    //           Lr3*(m(3)*S(rI2) - m(2)*S(rI1)) ;
    
    //  kbar4 =  Lr2*(m(3)*S(rJ3) - m(4)*S(rJ1)) - ...
    //           Lr3*(m(3)*S(rJ2) + m(5)*S(rJ1));
    
    //     ks3 = [o kbar2 o kbar4];
    
    static Matrix Sm(3,3);
    static Matrix kbar(12,3);
    
    Sm.addMatrix(0.0, SrI3,  m(3));
    Sm.addMatrix(1.0, SrI1,  m(1));
    
    kbar.addMatrixProduct(0.0, Lr2, Sm, -1.0);
    
    Sm.addMatrix(0.0, SrI2,  m(3));
    Sm.addMatrix(1.0, SrI1, -m(2));
    
    kbar.addMatrixProduct(1.0, Lr3, Sm,  1.0);
    
    kg.Assemble(kbar, 0, 3, 1.0);
    kg.AssembleTranspose(kbar, 3, 0, 1.0);
    
    Sm.addMatrix(0.0, SrJ3,  m(3));
    Sm.addMatrix(1.0, SrJ1, -m(4));
    
    kbar.addMatrixProduct(0.0, Lr2, Sm, 1.0);
    
    Sm.addMatrix(0.0, SrJ2, m(3));
    Sm.addMatrix(1.0, SrJ1, m(5));
    
    kbar.addMatrixProduct(1.0, Lr3, Sm,  -1.0);
    
    kg.Assemble(kbar, 0, 9, 1.0);
    kg.AssembleTranspose(kbar, 9, 0, 1.0);
    
    //opserr << "kg += ksigma3: " << kg;
    
    // Ksigma4 -------------------------------
    // Ks4_22 =  m(3)*( S(e2)*S(rI3) - S(e3)*S(rI2)) + ...
    //           m(1)*(-S(e1)*S(rI2) + S(e2)*S(rI1)) + ...
    //           m(2)*(-S(e1)*S(rI3) + S(e3)*S(rI1));
    
    // Ks4_44 = -m(3)*( S(e2)*S(rJ3) - S(e3)*S(rJ2)) + ...
    //           m(4)*(-S(e1)*S(rJ2) + S(e2)*S(rJ1)) + ...
    //           m(5)*(-S(e1)*S(rJ3) + S(e3)*S(rJ1));
    
    // Ks4 = [   O    O     O    O;
    //           O  Ks4_22  O    O;
    //           O    O     O    O;
    //           O    O     O  Ks4_44];
    
    static Matrix ks33(3,3);
    
    ks33.addMatrixProduct(0.0, Se2, SrI3,  m(3));
    ks33.addMatrixProduct(1.0, Se3, SrI2, -m(3));
    
    ks33.addMatrixProduct(1.0, Se2, SrI1,  m(1));
    ks33.addMatrixProduct(1.0, Se1, SrI2, -m(1));
    
    ks33.addMatrixProduct(1.0, Se3, SrI1,  m(2));
    ks33.addMatrixProduct(1.0, Se1, SrI3, -m(2));
    
    kg.Assemble(ks33, 3, 3, 1.0);
    
    ks33.addMatrixProduct(0.0, Se2, SrJ3, -m(3));
    ks33.addMatrixProduct(1.0, Se3, SrJ2,  m(3));
    
    ks33.addMatrixProduct(1.0, Se2, SrJ1,  m(4));
    ks33.addMatrixProduct(1.0, Se1, SrJ2, -m(4));
    
    ks33.addMatrixProduct(1.0, Se3, SrJ1,  m(5));
    ks33.addMatrixProduct(1.0, Se1, SrJ3, -m(5));
    
    kg.Assemble(ks33, 9, 9, 1.0);
    
    //opserr << "kg += ksigma4: " << kg;
    
    // Ksigma5 -------------------------------
    //
    //  Ks5 = [ Ks5_11   Ks5_12 -Ks5_11   Ks5_14;
    //          Ks5_12t    O    -Ks5_12t   O;
    //         -Ks5_11  -Ks5_12  Ks5_11  -Ks5_14;
    //          Ks5_14t     O   -Ks5_14t   O];
    
    // v = (1/Ln)*(m(2)*rI2 + m(3)*rI3 + m(5)*rJ2 + m(6)*rJ3);
    static Vector v(3);
    v.addVector (0.0, rI2, m(1));
    v.addVector (1.0, rI3, m(2));
    v.addVector (1.0, rJ2, m(4));
    v.addVector (1.0, rJ3, m(5));
    
    v /= Ln;
    
    //Ks5_11 = A*v*e1' + e1*v'*A + (e1'*v)*A;
    static Matrix m33(3,3);
    double  e1tv = 0;   // dot product e1. v
    
    for (i = 0; i < 3; i++)
        e1tv += e1(i) * v(i);
    
    ks33.addMatrix (0.0, A, e1tv);
    
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            m33(i,j) = v(i)*e1(j);
        
        ks33.addMatrixProduct (1.0, A, m33, 1.0);
        
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                m33(i,j) = e1(i)*v(j);
            
            ks33.addMatrixProduct (1.0, m33, A, 1.0);
            
            kg.Assemble(ks33, 0, 0,  1.0);
            kg.Assemble(ks33, 0, 6, -1.0);
            kg.Assemble(ks33, 6, 0, -1.0);
            kg.Assemble(ks33, 6, 6,  1.0);
            
            //Ks5_12 = -(m(2)*A*S(rI2) + m(3)*A*S(rI3));
            
            ks33.addMatrixProduct(0.0, A, SrI2, -m(1));
            ks33.addMatrixProduct(1.0, A, SrI3, -m(2));
            
            kg.Assemble(ks33, 0, 3,  1.0);
            kg.Assemble(ks33, 6, 3, -1.0);
            
            kg.AssembleTranspose(ks33, 3, 0,  1.0);
            kg.AssembleTranspose(ks33, 3, 6, -1.0);
            
            //  Ks5_14 = -(m(5)*A*S(rJ2) + m(6)*A*S(rJ3));
            
            ks33.addMatrixProduct(0.0, A, SrJ2, -m(4));
            ks33.addMatrixProduct(1.0, A, SrJ3, -m(5));
            
            kg.Assemble(ks33, 0, 9,  1.0);
            kg.Assemble(ks33, 6, 9, -1.0);
            
            kg.AssembleTranspose(ks33, 9, 0,  1.0);
            kg.AssembleTranspose(ks33, 9, 6, -1.0);
            
            //opserr << "kg += ksigma5: " << kg;
            
            // Ksigma -------------------------------
            static Vector rm(3);
            
            rm = rI3;
            rm.addVector (1.0, rJ3, -1.0); 
            //opserr << "ks2(r2,rI3-rJ3):\n "; 
            kg.addMatrix (1.0, this->getKs2Matrix(r2, rm), m(3));
            
            rm = rJ2;
            rm.addVector (1.0, rI2, -1.0); 
            //opserr << "ks2(r3,rJ2-rI2):\n "; 
            kg.addMatrix (1.0, this->getKs2Matrix(r3, rm), m(3));
            //opserr << "ks2(r2,rI1):\n "; 
            kg.addMatrix (1.0, this->getKs2Matrix(r2, rI1), m(1));
            //opserr << "ks2(r3,rI1):\n "; 
            kg.addMatrix (1.0, this->getKs2Matrix(r3, rI1), m(2));
            //opserr << "ks2(r2,rJ1):\n "; 
            kg.addMatrix (1.0, this->getKs2Matrix(r2, rJ1), m(4));
            //opserr << "ks2(r3,rJ1):\n "; 
            kg.addMatrix (1.0, this->getKs2Matrix(r3, rJ1), m(5));
            
            //opserr << "kg += ksigma2: " << kg;
            
            //  T * diag (M .* tan(thetal))*T' 
            
            for (k = 0; k < 6; k++)
            {
                factor = pl(k) * tan(ul(k));
                for (i = 0; i < 12; i++)
                    for (j = 0; j < 12; j++)
                        kg(i,j) += T(k,i) * factor * T(k,j);
            }
            
	    //            opserr << "COROATIONAL 3d: kg final: " << kg;
            
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
    
    // compute the tangent stiffness matrix in global coordinates
    kg.addMatrixTripleProduct(0.0, T, kl, 1.0);
    
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
        opserr << "\nCorotCrdTransf3d::computeElemtLengthAndOrien: 0 length\n";
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
        opserr << "\nCorotCrdTransf3d::getElementLengthAndOrientation";
        opserr << "\nvector v that defines plane xz is parallel to x axis\n";
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


const Matrix &
CorotCrdTransf3d::getKs2Matrix(const Vector &ri, const Vector &z) const
{
    static Matrix ks2(12,12);
    static Vector e1(3), r1(3);
    
    //opserr << "\ngetKs2Matrix:\n";
    //opserr << "ri: " << ri;
    
    //opserr << "z: " << z;
    
    //  Ksigma2 = [ K11   K12 -K11   K12;
    //              K12t  K22 -K12t  K22;
    //             -K11  -K12  K11  -K12;
    //              K12t  K22 -K12t  K22];
    
    // U = (-1/2)*A*z*ri'*A + ri'*e1*A*z*e1'/(2*Ln)+...
    //      z'*(e1+r1)*A*ri*e1'/(2*Ln);
    
    int i, j;
    double rite1 = 0;   // dot product ri . e1
    double zte1  = 0;   // dot product z  . e1
    double ztr1  = 0;   // dot product z  . r1
    
    for (i = 0; i < 3; i++)
    {
        e1(i) = e(i,0);
        r1(i) = Rbar(i,0);
        rite1 += ri(i)*e1(i);
        zte1  += z(i)*e1(i);
        ztr1  += z(i)*r1(i);
    }
    
    static Matrix zrit(3,3), ze1t(3,3);
    static Matrix rizt(3,3), r1e1t(3,3), rie1t(3,3);
    static Matrix e1zt(3,3);
    
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
        {
            zrit(i,j) = z(i)*ri(j);
            rizt(i,j) = ri(i)*z(j);
            ze1t(i,j) = z(i)*e1(j);
            e1zt(i,j) = e1(i)*z(j);
            r1e1t(i,j) = r1(i)*e1(j);
            rie1t(i,j) = ri(i)*e1(j);
        }
        
        static Matrix U(3,3);
        //opserr << " rite1: "<< rite1;
        //opserr << " zte1: "<< zte1;
        //opserr << " ztr1: "<< ztr1;
        
        //opserr << " A: "<< A;
        //opserr << " zrit: "<< zrit;
        //opserr << " ze1t: "<< ze1t;
        //opserr << " rie1t: "<< rie1t;
        
        //opserr << " Ln: "<< Ln;
        
        U.addMatrixTripleProduct(0.0, A, zrit, -0.5);
        
        //opserr << "U: A*zrit*zrit " << U;
        U.addMatrixProduct (1.0, A, ze1t,   rite1/(2*Ln));
        //opserr << "U: A*ze1t*rite1/(2*Ln) " << U;
        U.addMatrixProduct (1.0, A, rie1t, (zte1 + ztr1)/(2*Ln));
        
        //opserr << "U: " << U;
        static Matrix ks(3,3);
        
        //K11 = U + U' + ri'*e1*(2*(e1'*z)+z'*r1)*A/(2*Ln);
        
        ks.addMatrix (0.0, U, 1.0);
        
        // add matrix U transpose
        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                ks(i,j) += U(j,i);
            
            ks.addMatrix(1.0, A, rite1*(2*zte1 + ztr1)/(2*Ln));
            
            //opserr << "Ks211: " << ks;
            
            ks2.Zero();
            
            ks2.Assemble(ks, 0, 0,  1.0);
            ks2.Assemble(ks, 0, 6, -1.0);
            ks2.Assemble(ks, 6, 0, -1.0);
            ks2.Assemble(ks, 6, 6,  1.0);
            
            static Matrix Sri(3,3), Sr1(3,3), Sz(3,3), Se1(3,3);
            
            Sri = this->getSkewSymMatrix(ri);  
            Sr1 = this->getSkewSymMatrix(r1);
            Sz  = this->getSkewSymMatrix(z); 
            Se1 = this->getSkewSymMatrix(e1); 
            
            //K12 = (1/4)*(-A*z*e1'*Sri - A*ri*z'*Sr1 - z'*(e1+r1)*A*Sri);
            
            static Matrix m1(3,3);
            
            m1.addMatrixProduct(0.0, A, ze1t, -1.0);
            ks.addMatrixProduct(0.0, m1, Sri, 0.25);
            
            m1.addMatrixProduct(0.0, A, rizt, -1.0);
            ks.addMatrixProduct(1.0, m1, Sr1, 0.25);
            
            ks.addMatrixProduct(1.0, A, Sri, -0.25*(zte1+ztr1));
            
            //opserr << "Ks2_12: " << ks;
            
            ks2.Assemble(ks, 0, 3,  1.0);
            ks2.Assemble(ks, 0, 9,  1.0);
            ks2.Assemble(ks, 6, 3, -1.0);
            ks2.Assemble(ks, 6, 9, -1.0);
            
            ks2.AssembleTranspose(ks, 3, 0,  1.0);
            ks2.AssembleTranspose(ks, 3, 6, -1.0);
            ks2.AssembleTranspose(ks, 9, 0,  1.0);
            ks2.AssembleTranspose(ks, 9, 6, -1.0);
            
            //K22 = (1/8)*((-ri'*e1)*Sz*Sr1 + Sr1*z*e1'*Sri + ...
            //      Sri*e1*z'*Sr1 - (e1+r1)'*z*S(e1)*Sri + 2*Sz*Sri);
            
            ks.addMatrixProduct (0.0, Sz, Sr1, -0.125*(rite1));
            
            m1.addMatrixProduct (0.0, Sr1, ze1t, 1.0);
            ks.addMatrixProduct (1.0, m1, Sri, 0.125);
            
            m1.addMatrixProduct (0.0, Sri, e1zt, 1.0);
            ks.addMatrixProduct (1.0, m1, Sr1, 0.125);
            
            ks.addMatrixProduct (1.0, Se1, Sri, -0.125*(zte1 + ztr1));
            ks.addMatrixProduct (1.0, Sz, Sri, 0.25);
            
            //opserr << "Ks2_22: " << ks;
            
            //  Ksigma2 = [ K11   K12 -K11   K12;
            //              K12t  K22 -K12t  K22;
            //             -K11  -K12  K11  -K12;
            //              K12t  K22 -K12t  K22];
            
            ks2.Assemble(ks, 3, 3, 1.0);
            ks2.Assemble(ks, 3, 9, 1.0);
            ks2.Assemble(ks, 9, 3, 1.0);
            ks2.Assemble(ks, 9, 9, 1.0);
            
            return ks2;
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
