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
// $Date: 2008-12-03 23:40:07 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CorotCrdTransfWarping2d.h,v $ 
                                        
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the class definition for
// CorotCrdTransfWarping2d.h. CorotCrdTransfWarping2d provides the
// abstraction of a corotation transformation for a planar frame element

// What: "@(#) CorotCrdTransfWarping2d.h, revA"

#ifndef CorotCrdTransfWarping2d_h
#define CorotCrdTransfWarping2d_h

#include <CrdTransf.h>
#include <Vector.h>
#include <Matrix.h>

class CorotCrdTransfWarping2d: public CrdTransf
{
public:
    CorotCrdTransfWarping2d(int tag, const Vector &rigJntOffsetI, const Vector &rigJntOffsetJ);
    
    CorotCrdTransfWarping2d();
    ~CorotCrdTransfWarping2d();
    
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
    
// AddingSensitivity:BEGIN //////////////////////////////////
    const Vector &getBasicDisplSensitivity(int gradNumber);
    const Vector &getGlobalResistingForceShapeSensitivity(const Vector &q,
							  const Vector &p0,
							  int gradNumber);
    const Vector &getBasicTrialDispShapeSensitivity(void);
    bool isShapeSensitivity(void);
    double getdLdh(void);
    double getd1overLdh(void);
// AddingSensitivity:END //////////////////////////////////

    CrdTransf *getCopy2d(void);
    
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
    const Matrix &getGlobalMatrixFromLocal(const Matrix &ml);
    // functions used in post-processing only    
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps);    
    
private:
    int compElemtLengthAndOrient(void);
    int compElemtLengthAndOrientWRTLocalSystem(const Vector &ul);
    void transfLocalDisplsToBasic(const Vector &ul);
    void getTransfMatrixLocalGlobal(Matrix &Tlg);
    void getTransfMatrixBasicLocal(Matrix &Tbl);
    const Matrix &getGeomStiffMatrix(const Vector &pb) const;
    
    // internal data
    Node *nodeIPtr, *nodeJPtr; // pointers to the element two endnodes
    
    Vector nodeIOffset, nodeJOffset;  // rigid joint offsets
    
    double cosTheta, sinTheta; // direction cossines of undeformed element wrt to global system 
    double cosAlpha, sinAlpha; // direction cossines of deformed element wrt to local system
    double L;                  // undeformed element length
    double Ln;                 // deformed element length
    double Lx, Ly;             // components of the deformed member
    double Lxdot, Lydot;       // time derivatives of components of the deformed member
    double Lxdotdot, Lydotdot; // double time derivatives of components of the deformed member
    
    Vector ub;                 // basic displacements
    Vector ubcommit;           // committed basic displacements
    Vector ubpr;               // previous basic displacements
    
    static Matrix Tlg;         // matrix that transforms from global to local coordinates
    static Matrix Tbl;         // matrix that transforms from local  to basic coordinates
    static Matrix kg;     
    static Vector uxg;     
    static Vector pg;     
    static Vector dub;     
    static Vector Dub;     
    
    double *nodeIInitialDisp, *nodeJInitialDisp;
    bool initialDispChecked;
    bool nodeOffsets;
};
#endif
