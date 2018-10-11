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
// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthInterface2D.h,v $
// $Revision: 1.1 $
// $Date: July 02 2010

// Written: Roozbeh Geraili Mikola  (roozbehg@berkeley.edu)
//          Prof. Nicholas Sitar (nsitar@ce.berkeley.edu)

#ifndef ZeroLengthInterface2D_h
#define ZeroLengthInterface2D_h

/*----+----+----+----+----+----+----+----+----+----+----+----+----+----+------*
 |                                                                            |
 |                    ZeroLengthContactNTS2D element                        |
 +                                                                            +
 |----------------------------------------------------------------------------|
 |                                                                            |
 +        Authors: Roozbeh Geraili Mikola  AND  Professor Nicholas Sitar      +
 |                                                                            |
 |		  Department of Civil and Environmental Engineering                   |
 +		  University of California, Berkeley, CA 94720, USA                   +
 |                                                                            |
 |        Email: roozbehg@berkeley.edu                                        |
 +                                                                            +
 |  Disclaimers:                                                              |
 |  (1) Frame of this code is based on zeroLength element                     |
 +  (2) Questions regarding this code should be directed to Roozbeh G. Mikola +
 |  (3) Documentation could be found at                                       |
 |             www2.decf.berkeley.edu/~roozbehg/OpenSees.html                 |
 +                                                                            +
 |  Development History:                                                      |
 |  Created       07/02/2010                                                  |
 |                                                                            |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+------*/

/*
 element zeroLengthInterface2D eleTag? -sNdNum sNdNum? -mNdNum mNdNum? �dof sdof? mdof? -Nodes Nodes? Kn? Kt? phi? 

 Description: This file contains the class definition for ZeroLengthContact2D.

[1] The contact element is node-to-segment (NTS) contact. The relation follows 
    Mohr-coulomb frictional law: T = N * tan($phi), where T is tangential force 
	and N is normal force across the interface. $phi is friction angle.
[2] For 2D contact, slave nodes and master nodes must be 2 DOF and notice note 
    that the slave and master nodes must be entered in counterclockwise order.
[3] The resulting tangent from the contact element is non-symmetric. Switch to 
    the non-symmetric matrix solver if convergence problem is experienced. 
[4] As opposed to node-to-node contact, predefined normal vector for node-to-segment (NTS) 
    element is not required because contact normal will be calculated automatically 
	at each step. And also this element can handle contact between different DOFs such 
	as beam-beam, beam-solid and solid-solid.
[5] The contact element is implemented to handle large deformations.

  References:

[1] P. Wriggers, V.T. Vu and E. Stein, Finite-element formulation of large deformation
    impact�contact problems with friction, Comput. Struct. 37 (1990), pp. 319�331.
[2] Peter Wriggers. Computational Contact Mechanics. John Wiley & Sons Ltd. Chichester, 2002.


  www2.decf.berkeley.edu/~roozbehg/
*/

#include <Element.h>
#include <Matrix.h>

// Tolerance for zero length of element
#define	LENTOL 1.0e-6

class Node;
class Channel;
//class UniaxialMaterial;
class Response;

class ZeroLengthInterface2D: public Element
{
  public:
  // Constructor
  ZeroLengthInterface2D(int tag, int sNdNum, int mNdNum, int sDof, int mDof, const ID& Nodes,
		                   double Kn, double Kt, double fRatio);
  // Null constructor
  ZeroLengthInterface2D();

  ~ZeroLengthInterface2D();

  // public methods to obtain information about dof & connectivity
  int getNumExternalNodes(void) const;
  const ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  int getNumDOF(void);
  void setDomain(Domain *theDomain);

  // public methods to set the state of the element
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  int update(void);

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);
  const Matrix &getDamp(void);
  const Matrix &getMass(void);
  void zeroLoad(void);
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);
  const Vector &getResistingForce(void);
  const Vector &getResistingForceIncInertia(void);

  // public methods for element output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
  void Print(OPS_Stream &s, int flag =0);
  Response *setResponse(const char **argv, int argc, Information &eleInformation);
  int getResponse(int responseID, Information &eleInformation);

  //void updateDir (const Vector& x, const Vector& y);

 protected:

 private:
  //int    directionID;
  ID     connectedExternalNodes;         // contains the tags of the end nodes
  int  numberNodes  ;
  Node **nodePointers;   // node pointer
  // contact forces
  Vector pressure;    // contact pressure at n+1
  // parameters
  Vector normal_gap;  // normal gap of time n+1 step
  Vector shear_gap;   // shear gap of time n+1 step
  double Kn;			// normal penalty
  double Kt;			// tangential penalty
  double fc;			// friction ratio
  // contact point and contact plane stuffs
  Vector stored_shear_gap;  // (keci)
  double xi;       // trial stick point in local coord
  // Normal and Tangental Vectors for Elemental Nodes, (4*1)
  Vector N;
  Vector T;
  Vector ContactNormal;  // out normal of master element
  int ContactFlag;                    // 0: not contact; 1: stick; 2: slide
  int numDOF;	                        // number of dof for ZeroLength
  // detect the contact and set flag
  int contactDetect(int s, int m1, int m2, int stage);
  //form residual and tangent
  void formLocalResidAndTangent( int tang_flag , int slave, int master1, int master2, int stage);
  void formGlobalResidAndTangent(int tang_flag );
  void GlobalResidAndTangentOrder(int slave, int master1, int master2);
  void GlobalResidAndTangentOrder2(int slave, int master1, int master2);
  Matrix *Ki; 	    	// pointer to objects matrix (a class Matrix)
  Vector *load;         	// pointer to objects vector (a class Vector)
  // variables for 2D contact
  Matrix stiff;   // for stiff matrix
  Vector resid;   // for force residual vector
  Matrix zeroMatrix;
  int SlaveNodeNum;
  int MasterNodeNum;
  int SlaveDof;
  int MasterDof;
  double *restore_shear_gap;
  int loctoglob[6];
};

#endif











