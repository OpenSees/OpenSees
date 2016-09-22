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

// $Revision: 1.6 $
// $Date: 2009-08-19 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CrdTransf.h,v $

// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the class definition for 
// CrdTransf.h. CrdTransf provides the abstraction of a frame 
// coordinate transformation. It is an abstract base class and 
// thus no objects of  it's type can be instatiated. It has pure 
// virtual functions which  must be implemented in it's derived classes.
//
// What: "@(#) CrdTransf.h, revA"

#ifndef CrdTransf_h
#define CrdTransf_h

#include <MovableObject.h>
#include <TaggedObject.h>

class Vector;
class Matrix;
class Node;

// class definition

class CrdTransf: public TaggedObject, public MovableObject
{
public:
    CrdTransf(int tag, int classTag);
    CrdTransf();
    virtual ~CrdTransf();

    virtual CrdTransf *getCopy2d(void) {return 0;};
    virtual CrdTransf *getCopy3d(void) {return 0;};
    virtual int getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis) {return -1;};
    
    virtual int    initialize(Node *node1Pointer, Node *node2Pointer) = 0;
    virtual int    update(void) = 0;
    virtual double getInitialLength(void) = 0;
    virtual double getDeformedLength(void) = 0;
    
    virtual int commitState(void) = 0;
    virtual int revertToLastCommit(void) = 0;
    virtual int revertToStart(void) = 0;
    
    virtual const Vector &getBasicTrialDisp(void) = 0;
    virtual const Vector &getBasicIncrDisp(void) = 0;
    virtual const Vector &getBasicIncrDeltaDisp(void) = 0;
    virtual const Vector &getBasicTrialVel(void) = 0;
    virtual const Vector &getBasicTrialAccel(void) = 0;
    
    // AddingSensitivity:BEGIN //////////////////////////////////
    virtual const Vector &getBasicDisplSensitivity(int gradNumber);
    virtual const Vector &getBasicDisplSensitivity(int gradNumber,int); // used by Quan 
    //virtual const Vector &getGlobalResistingForceShapeSensitivity(const Vector &basicForce, const Vector &uniformLoad);
    virtual const Vector &getGlobalResistingForceShapeSensitivity(const Vector &pb, const Vector &p0, int gradNumber);
    virtual const Vector &getGlobalResistingForceShapeSensitivity(const Vector &pb, const Vector &p0);
    virtual const Vector &getBasicTrialDispShapeSensitivity(void);
    virtual bool isShapeSensitivity(void) {return false;}
    virtual double getdLdh(void) {return 0.0;}
    virtual double getd1overLdh(void) {return 0.0;}
    // AddingSensitivity:END //////////////////////////////////
    
    virtual const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &uniformLoad) = 0;
    virtual const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce) = 0;
    virtual const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff) = 0;
    
    // method used to rotate consistent mass matrix
    virtual const Matrix &getGlobalMatrixFromLocal(const Matrix &local) = 0;
    
    // methods used in post-processing only
    virtual const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords) = 0;
    virtual const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps) = 0;
    
protected:
    
private:
};

// some additional methods related to prototypes created for copy constructors
extern bool       OPS_addCrdTransf(CrdTransf *newComponent);
extern CrdTransf *OPS_getCrdTransf(int tag);
extern void       OPS_clearAllCrdTransf(void);

#endif
