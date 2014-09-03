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

// $Revision: 1.7 $
// $Date: 2005-12-15 00:30:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/PDeltaCrdTransf3d.h,v $

// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the class definition for
// PDeltaCrdTransf3d.h. PDeltaCrdTransf3d provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems

// What: "@(#) PDeltaCrdTransf3d.h, revA"

#ifndef PDeltaCrdTransf3d_h
#define PDeltaCrdTransf3d_h

#include <CrdTransf.h>
#include <Vector.h>
#include <Matrix.h>

class PDeltaCrdTransf3d: public CrdTransf
{
public:
    PDeltaCrdTransf3d(int tag, const Vector &vecInLocXZPlane);
    PDeltaCrdTransf3d(int tag, const Vector &vecInLocXZPlane,
        const Vector &rigJntOffsetI,
        const Vector &rigJntOffsetJ);
    
    PDeltaCrdTransf3d();
    ~PDeltaCrdTransf3d();
    
    const char *getClassType() const {return "PDeltaCrdTransf3d";};
    
    int initialize(Node *node1Pointer, Node *node2Pointer);
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
    
    int getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);
    
private:
    int computeElemtLengthAndOrient(void);
    void compTransfMatrixLocalGlobal(Matrix &Tlg);
    
    // internal data
    Node *nodeIPtr, *nodeJPtr;  // pointers to the element two endnodes
    
    double *nodeIOffset, *nodeJOffset;	// rigid joint offsets
    
    double R[3][3];	// rotation matrix
    double L;       // undeformed element length
    double ul17;	// Transverse local displacement offsets of P-Delta
    double ul28;

    static Matrix Tlg;  // matrix that transforms from global to local coordinates
    static Matrix kg;   // global stiffness matrix

    double *nodeIInitialDisp, *nodeJInitialDisp;
    bool initialDispChecked;
};

#endif
