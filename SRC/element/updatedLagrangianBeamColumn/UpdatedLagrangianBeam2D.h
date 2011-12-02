/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** See file 'COPYRIGHT'  in main directory for information on usage   **
** and redistribution of OpenSees, and for a DISCLAIMER OF ALL        **
** WARRANTIES.                                                        **
**                                                                    **
** UpdatedLagrangianBeam2D.h: interface for UpdatedLagrangianBeam2D   **
** Developed by:                                                      **
**    Rohit Kaul       (rkaul@stanford.edu)                           **
**    Greg Deierlein   (ggd@stanford.edu)                             **
**                                                                    **
**           John A. Blume Earthquake Engineering Center              **
**                    Stanford University                             **
** ****************************************************************** **/


// UpdatedLagrangianBeam2D.h: interface for the UpdatedLagrangianBeam2D class
// Written: rkaul
//
// Description: This file contains the class definition for UpdatedLagrangianBeam2D.

// UpdatedLagrangianBeam2D is an abstract class providing most methods required by
// the base class "Element", for 2D beam-column elements. Geometric
// nonlinearity is incorporated at this level using updated lagrangian
// formulation

#ifndef UpdatedLagrangianElement2D
#define UpdatedLagrangianElement2D


#include <bool.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Element.h>
#include <ID.h>

#define min_(a,b)  (((a) < (b)) ? (a) : (b))
#define max_(a,b)  (((a) > (b)) ? (a) : (b))
#define sign(a)    (((a) < (0)) ?(-1) : (1))

class Response;

class UpdatedLagrangianBeam2D : public Element
{
public:

  UpdatedLagrangianBeam2D(int classTag);
  UpdatedLagrangianBeam2D(int tag, int classTag, int nd1, int nd2, bool islinear = false);

  virtual ~UpdatedLagrangianBeam2D();

  //////////////////////////////////////////////////////////////////////
  // Overridden public methods, defined in Element class
  // (see ~/OpenSees/SRC/element/Element.h)
  //////////////////////////////////////////////////////////////////////
 public:
  virtual int update(void);
  // void	setEndRelease(ID &g);
  
  // methods dealing with nodes and number of external dof
  int		getNumExternalNodes(void) const;
  const	ID &getExternalNodes(void);
  Node **getNodePtrs(void);

  int		getNumDOF(void);
  virtual void setDomain(Domain *theDomain);
  
  virtual int commitState(void);
  virtual int revertToLastCommit(void);
  
  // methods to return the current linearized stiffness,
  // damping and mass matrices
  virtual const	Matrix &getTangentStiff(void);
  virtual const Matrix &getInitialStiff(void);
  virtual const Matrix &getMass(void);
  
  // methods for returning and applying loads
  virtual Vector &getUVLoadVector(double q1, double q2);
  void	zeroLoad(void);
  int		addLoad(const Vector &load);
  int     addLoad(ElementalLoad *theLoad, double loadFactor)
    { return -1;}
  int     addInertiaLoadToUnbalance(const Vector &accel)
    { return -1;}
  
  virtual const Vector &getResistingForce(void);
  const	Vector &getResistingForceIncInertia(void);
  
  // method for graphics
  virtual int displaySelf(Renderer &theViewer,
			  int displayMode, float fact);
  
  virtual Response *setResponse(const char **argv, int argc,
				Information &eleInformation);
  virtual int getResponse(int responseID, Information &eleInformation);
  
 protected:
  
  virtual void  getLocalStiff(Matrix &K)=0;
  virtual void  getLocalMass(Matrix &M)=0;
  
  
  void	getIncrLocalDisp(Vector &localDisp);
  void	getTrialNaturalDisp(Vector &localDisp);
  void    getIncrNaturalDisp(Vector &nDisp);
  void    getConvLocalDisp(Vector &lDisp);
  void    getTrialLocalDisp(Vector &lDisp);
  void    getTrialLocalForce(Vector &force);
  
  virtual	void	updateState(void);

  void  addInternalGeomStiff(Matrix &K);
  void  addExternalGeomStiff(Matrix &K);
  
  void  transformToGlobal(Matrix &K);
  
 protected:
  
  bool    isLinear;
  int     numDof;
  double  L, sn, cs, massDof;
  
  ID      connectedExternalNodes;
  Vector  load;
  
  Node    *end1Ptr, *end2Ptr;
  double  L_hist, cs_hist, sn_hist;
  
  Vector  eleForce, eleForce_hist;

  int     nodeRecord, dofRecord;
  int     m_Iter;

  Matrix *Ki;

  static Matrix K, Kg, Kt; // stiffness matrices
  static Matrix M; // mass matrix
  static Matrix D; // damping matrix
  static Matrix T; // transformation matrix
  
  static Vector disp;
  static Vector force;
  
  // ZeroVector and ZeroMatrix should always have zero value,
  // used for quick return
  static Vector ZeroVector;
  static Matrix ZeroMatrix;
  
  // used for temporarily storing nodal displacements
  static Vector end1IncrDisp;
  static Vector end2IncrDisp;

  static Node *theNodes[2];
};


#endif // !defined(UpdatedLagrangianElement2D)
