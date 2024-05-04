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


// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthContact2D.h,v $
// $Revision: 1.2 $
// $Date: 2007-11-28 00:08:58 $


#ifndef ZeroLengthContact2D_h
#define ZeroLengthContact2D_h


/*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                          |
 |                    ZeroLengthContact2D element                           |
 +                                                                          +
 |--------------------------------------------------------------------------|
 |                                                                          |
 +             Authors: Gang Wang  AND  Professor Nicholas Sitar            +
 |                                                                          |
 |             Department of Civil and Environmental Engineering            |
 +             University of California, Berkeley, CA 94720, USA            +
 |                                                                          |
 |             Email: wang@ce.berkeley.edu (G.W.)                           |
 +                                                                          +
 |  Disclaimers:                                                            |
 |  (1) Frame of this code is based on zeroLength element                   |
 +  (2) Questions regarding this code should be directed to Gang Wang       +
 |  (3) Documentation could be found at                                     |
 |             www.ce.berkeley.edu/~wang/OpenSees.html                      |
 +                                                                          +
 |  Development History:                                                    |
 |  Created       -- Jan/23/2004                                            |
 +  Final Release -- June 2004                                              +
 |  Add Augmented Lagrange  --Oct 21,2004                                   |
 |                                                                          |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/


/*
 element zeroLengthContact2D  $eleID $sNdID $mNdID $Kn $Kt $fs -normal $Nx $Ny

 Description: This file contains the class definition for ZeroLengthContact2D.
 (1) A ZeroLengthContact2D element is defined by two nodes with the same coordinate
     in R^2
 (2) Normal to the contact plane through these contact nodes is preset by the user
     to be either radial or any of x,y,z axis.
 (3) Penalty (Kn,Kt) is used to enforce the constraints, i.e.
           No (in fact, very small) penetration in the normal direction, and
	       Coulomb frictional law in the tangential direction.
 (4) Backward Euler return mapping scheme is used for implicit formulation of
     residual and consistent tangent modulus.
 (5) This element is a modified version of continuum-based implicit formulation
     (Laursen and Simo,1993) for small deformation. Contact detection is simplified
	 to be only related to these two nodes, and contributions of geometric
	 nonlinearities to consistent tangent moudulus are neglected.

  References:
  Laursen,T.A. and Simo, J.C. "A Continuum-based Finite Element Formulation for the
  Implicit Solution of Multibody, Large Deformation Frictional Contact Problems".
  Int. J. Numer. Methods Eng. , 36, 3451-3485 (1993)

  Armero, F. and Petocz, E. "A New Dissipative Time-Stepping Algorithm for Frictional
  Contact Problems: Fromulation and Analysis", Comp. Methods Appl. Mech. Eng., 179,
  151-178 (1999)

  www.ce.berkeley.edu/~wang
*/


#include <Element.h>
#include <Matrix.h>


// Tolerance for zero length of element
#define	LENTOL 1.0e-6

class Node;
class Channel;
//class UniaxialMaterial;
class Response;

class ZeroLengthContact2D: public Element
{
  public:

  // Constructor
  ZeroLengthContact2D(int tag, int Nd1, int Nd2,
          double Kn, double Kt, double fRatio,
          const Vector& normal);

  // Null constructor
  ZeroLengthContact2D();
  ~ZeroLengthContact2D();

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

  Response *setResponse(const char **argv, int argc, OPS_Stream &output);
  int getResponse(int responseID, Information &eleInformation);

  //void updateDir (const Vector& x, const Vector& y);

 protected:


 private:
  //int    directionID;
  ID     connectedExternalNodes;         // contains the tags of the end nodes
  static const int  numberNodes  ;
  Node *nodePointers[2];   // node pointer


  // contact forces
  double pressure;    // contact pressure at n+1
  double lambda;      // lagrange multiplier, used for augmented Lagrange
  double t1;          // friction in local dir1
  double t2;          // friction in local dir2

  // parameters
  double gap;         // gap of time n+1 step
  double gap_n;       // gap of time n step
  double Kn;          // normal penalty
  double Kt;          // tangential penalty
  double fs;          // friction ratio

  // contact point and contact plane stuffs
  double stickPt;  // (keci)
  double xi;       // trial stick point in local coord

  // Normal and Tangental Vectors for Elemental Nodes, (4*1)
  Vector N;
  Vector T;
  Vector ContactNormal;  // out normal of primary element

  int ContactFlag;                    // 0: not contact; 1: stick; 2: slide
  int numDOF;	                        // number of dof for ZeroLength

  // detect the contact and set flag
  int contactDetect();

  //form residual and tangent
  void formResidAndTangent(int tang_flag ) ;

  //static variables for 2D contact
  static Matrix stiff4;   // for stiff matrix
  static Vector resid4;   // for force residual vector

  static Matrix stiff6;   // for stiff matrix
  static Vector resid6;   // for force residual vector
};

#endif
