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
// $Date: 2008-04-14 21:22:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/dispBeamColumnInt/LinearCrdTransf2dInt.h,v $
                                                                        
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the class definition for
// LinearCrdTransf2d.h. LinearCrdTransf2d provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems

// What: "@(#) LinearCrdTransf2d.h, revA"

#ifndef LinearCrdTransf2dInt_h
#define LinearCrdTransf2dInt_h

#include <CrdTransf.h>
#include <Vector.h>
#include <Matrix.h>


class LinearCrdTransf2dInt: public CrdTransf
{
 public:
  LinearCrdTransf2dInt (int tag);
  LinearCrdTransf2dInt (int tag,
			const Vector &rigJntOffsetI,
			const Vector &rigJntOffsetJ);
  
  LinearCrdTransf2dInt();
  ~LinearCrdTransf2dInt();

  CrdTransf *getCopy2d(void);
  
  int    initialize(Node *node1Pointer, Node *node2Pointer);
  int    update(void);
  double getInitialLength(void);
  double getDeformedLength(void);
  
  int commitState(void);
  int revertToLastCommit(void);        
  int revertToStart(void);
  
  // AddingSensitivity:BEGIN //////////////////////////////////
  const Vector &getBasicDisplSensitivity (int gradNumber);
  const Vector &getGlobalResistingForceShapeSensitivity (const Vector &basicForce, const Vector &p0);
  const Vector &getBasicTrialDispShapeSensitivity     (void);
  // AddingSensitivity:END //////////////////////////////////

  const Vector &getBasicTrialDisp(void);
  const Vector &getBasicTrialDispInt(void);
  const Vector &getBasicTrialVel(void);
  const Vector &getBasicTrialAccel(void);
  
  const Vector &getBasicIncrDisp      (void);
  const Vector &getBasicIncrDeltaDisp (void);
  
  const Vector &getGlobalResistingForce (const Vector &basicForce, const Vector &p0);
  const Vector &getGlobalResistingForceInt (const Vector &basicForce, const Vector &p0);
  
  const Matrix &getGlobalStiffMatrix    (const Matrix &basicStiff, const Vector &basicForce);
  const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff);
  
  const Matrix &getGlobalStiffMatrixInt    (const Matrix &basicStiff, const Vector &basicForce);
  const Matrix &getInitialGlobalStiffMatrixInt(const Matrix &basicStiff);
  
  
  int sendSelf(int cTag, Channel &theChannel);
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  
  void Print(OPS_Stream &s, int flag =0);
  
  // method used to rotate consistent mass matrix
  const Matrix &getGlobalMatrixFromLocal(const Matrix &local);
  
  // methods used in post-processing only
  const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
  const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
  
 private:
  int  computeElemtLengthAndOrient(void);
  void compTransfMatrixLocalGlobal(Matrix &Tlg);
  
  // internal data
  Node *nodeIPtr;
  Node *nodeJPtr;          // pointers to the element two endnodes
  
  double *nodeIOffset, *nodeJOffset;	// rigid joint offsets
  
  double cosTheta, sinTheta;  // direction cosines of undeformed element wrt to global system 
  double L;  // undeformed element length
  
  static Matrix Tlg;  // matrix that transforms from global to local coordinates
  static Matrix kg;   // global stiffness matrix
};

#endif


