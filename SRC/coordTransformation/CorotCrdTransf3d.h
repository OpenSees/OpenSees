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
                                                                        
// $Revision: 1.9 $                                                              
// $Date: 2006-01-13 01:07:48 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CorotCrdTransf3d.h,v $      


// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the class definition for
// CorotCrdTransf3d.h. CorotCrdTransf3d provides the
// abstraction of a corotation transformation for a spatial frame element

// What: "@(#) CorotCrdTransf3d.h, revA"

#ifndef CorotCrdTransf3d_h
#define CorotCrdTransf3d_h

#include <CrdTransf.h>
#include <Vector.h>
#include <Matrix.h>

class CorotCrdTransf3d: public CrdTransf
{
public:
    CorotCrdTransf3d(int tag, const Vector &vecInLocXZPlane,
        const Vector &rigJntOffsetI, const Vector &rigJntOffsetJ);
    
    CorotCrdTransf3d();
    ~CorotCrdTransf3d();
    
    const char *getClassType() const {return "CorotCrdTransf3d";};
    
    int initialize(Node *nodeIPointer, Node *nodeJPointer);
    int update(void);
    double getInitialLength(void);
    double getDeformedLength(void);
    
    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);
    
    const Vector &getBasicTrialDisp(void);
    const Vector &getBasicIncrDisp(void);
    const Vector &getBasicIncrDeltaDisp(void);
	const Vector &getBasicTrialVel(void);
	const Vector &getBasicTrialAccel(void);
    
    const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &p0);
    const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce);
    const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff);
    
    CrdTransf *getCopy3d(void);
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);
    
    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps);    
    
    int  getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);
  int getRigidOffsets(Vector &offsets);
  
private:
    void compTransfMatrixBasicGlobal(void);
    void compTransfMatrixBasicGlobalNew(void);
    void compTransfMatrixLocalGlobal(Matrix &Tlg);
    void compTransfMatrixBasicLocal(Matrix &Tbl);
    const Vector &getQuaternionFromRotMatrix(const Matrix &RotMatrix) const;
    const Vector &getQuaternionFromPseudoRotVector(const Vector &theta) const;
    const Vector &getTangScaledPseudoVectorFromQuaternion(const Vector &theta) const;
    const Vector &quaternionProduct(const Vector &q1, const Vector &q2) const;
    const Matrix &getRotationMatrixFromQuaternion(const Vector &q) const;
    const Matrix &getRotMatrixFromTangScaledPseudoVector(const Vector &w) const;
    const Matrix &getSkewSymMatrix(const Vector &theta) const;
    const Matrix &getLMatrix(const Vector &ri) const;
    const Matrix &getKs2Matrix(const Vector &ri, const Vector &z) const;
    
    // internal data
    Node *nodeIPtr, *nodeJPtr;  // pointers to the element two endnodes
    
    Vector vAxis;               // Vector that lies in local plane xz
    Vector nodeIOffset, nodeJOffset;    // rigid joint offsets
    
    Vector xAxis;               // local x axis
    double L;                   // undeformed element length
    double Ln;                  // deformed element length
    
    Matrix R0;                  // rotation matrix from local to global coordinates
                                // (the columns of which are the element local axes)
    Vector alphaIq;             // quaternion for node I
    Vector alphaJq;             // quaternion for node I
    
    Vector alphaIqcommit;       // committed quaternion for node I
    Vector alphaJqcommit;       // committed quaternion for node J
    Vector alphaI;              // last trial rotations end i
    Vector alphaJ;              // last trial rotatations end j
    
    Vector ul;                  // local displacements
    Vector ulcommit;            // committed local displacements
    Vector ulpr;                // previous local displacements
    
    static Matrix RI;           // nodal triad for node 1
    static Matrix RJ;           // nodal triad for node 2
    static Matrix Rbar;         // mean nodal triad 
    static Matrix e;            // base vectors
    static Matrix Tp;           // transformation matrix to renumber dofs
    static Matrix T;            // transformation matrix from basic to global system
    static Matrix Tlg;          // transformation matrix from global to local system
    static Matrix TlgInv;       // inverse of transformation matrix from global to local system
    static Matrix Tbl;          // transformation matrix from local to basic system
    static Matrix kg;           // global stiffness matrix
    static Matrix Lr2, Lr3, A;  // auxiliary matrices
    
    double *nodeIInitialDisp, *nodeJInitialDisp;
    bool initialDispChecked;
};
#endif
