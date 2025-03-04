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

#ifndef ZeroLengthNew_h
#define ZeroLengthNew_h

/*----+----+----+----+----+----+----+----+----+----+----+----+----+----+------*
 |                                                                            |
 |                    ZeroLengthContactNTS2D element                          |
 +                                                                            +
 |----------------------------------------------------------------------------|
 |                                                                            |
 +        Authors: Roozbeh Geraili Mikola  AND  Professor Nicholas Sitar      +
 |                                                                            |
 |        Department of Civil and Environmental Engineering                   |
 +        University of California, Berkeley, CA 94720, USA                   +
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
  element zeroLengthInterface2D eleTag? -sNdNum sNdNum? -pNdNum pNdNum? ?dof sdof? mdof? -Nodes Nodes? Kn? Kt? phi?

  Description: This file contains the class definition for ZeroLengthContact2D.

 [1] The contact element is node-to-segment (NTS) contact. The relation follows
     Mohr-coulomb frictional law: T = N * tan($phi), where T is tangential force
     and N is normal force across the interface. $phi is friction angle.
 [2] For 2D contact, secondary nodes and primary nodes must be 2 DOF and notice note
     that the secondary and primary nodes must be entered in counterclockwise order.
 [3] The resulting tangent from the contact element is non-symmetric. Switch to
     the non-symmetric matrix solver if convergence problem is experienced.
 [4] As opposed to node-to-node contact, predefined normal vector for node-to-segment (NTS)
     element is not required because contact normal will be calculated automatically
     at each step. And also this element can handle contact between different DOFs such
     as beam-beam, beam-solid and solid-solid.
 [5] The contact element is implemented to handle large deformations.

   References:

 [1] P. Wriggers, V.T. Vu and E. Stein, Finite-element formulation of large deformation
     impact?contact problems with friction, Comput. Struct. 37 (1990), pp. 319?331.
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

class ZeroLengthNew : public Element
{
public:
    // Constructor
    ZeroLengthNew(int tag, int sNdNum, int pNdNum, int sDof, int mDof, const ID& Nodes,
        double Kn, double Kt, double fRatio);
    // Null constructor
    ZeroLengthNew();

    ~ZeroLengthNew();

    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes(void) const;
    const ID& getExternalNodes(void);
    Node** getNodePtrs(void);
    int getNumDOF(void);
    void setDomain(Domain* theDomain);

    // public methods to set the state of the element
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix& getTangentStiff(void);
    const Matrix& getInitialStiff(void);
    const Matrix& getDamp(void);
    const Matrix& getMass(void);
    void zeroLoad(void);
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);
    const Vector& getResistingForce(void);
    const Vector& getResistingForceIncInertia(void);

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
    int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);
    void Print(OPS_Stream& s, int flag = 0);
    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInformation);
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);
    //void updateDir (const Vector& x, const Vector& y);

protected:

private:
    //int    directionID;
    ID     connectedExternalNodes;         // contains the tags of the end nodes
    int  numberNodes;
    Node** nodePointers;   // node pointer
    // contact forces
    Vector pressure;    // contact pressure at n+1
    // parameters
    Vector normal_gap;  // normal gap of time n+1 step, this will be updated under the loop between secondary pt and all primary pt segment
    Vector normal_gap_contact;  // normal gap of time n+1 step, only for the contact pair of sec. and pri.
    Vector shear_gap;   // delta shear gap of time n+1 step this will be updated under the loop between secondary pt and all primary pt segment
    Vector shear_gap_contact;   // delta shear gap time n+1 step, only for the contact pair of sec. and pri.
    Vector shear_force_c;   // the history variable, committed (total) shear force of time n step
    Vector shear_force_delta;  // shear force increment of time n+1 step
    Vector shear_gap_c;   // committed (total) shear gap of time n step.

    //Vector uxs_deltax;   // x trial disp. for secoundary pt after the last committed step
    //Vector uxs_deltay;   // x trial disp. for primary pt1 after the last committed step
    int fricChange;     // flag for chaning friction
    int comitn;     // number of commitment
    double Kn;			// normal penalty
    double Kt;			// tangential penalty
    double fc;			// friction ratio
    // contact point and contact plane stuffs
    Vector stored_shear_gap;  // (keci)
    double xi;       // trial stick point in local coord
    // Normal and Tangental Vectors for Elemental Nodes, (4*1)
    Vector N;
    Vector T;
    Vector ContactNormal;  // out normal of primary element
    int ContactFlag;                    // 0: not contact; 1: stick; 2: slide
    int numDOF;	                        // number of dof for ZeroLength
    // detect the contact and set flag
    int contactDetect(int s, int m1, int m2, int stage);
    //form residual and tangent
    void formLocalResidAndTangent(int tang_flag, int secondary, int primary1, int primary2, int stage);
    void formGlobalResidAndTangent(int tang_flag);
    void GlobalResidAndTangentOrder(int secondary, int primary1, int primary2);
    void GlobalResidAndTangentOrder2(int secondary, int primary1, int primary2);
    Matrix* Ki; 	    	// pointer to objects matrix (a class Matrix)
    Vector* load;         	// pointer to objects vector (a class Vector)
    // variables for 2D contact
    Matrix stiff;   // for stiff matrix
    Vector resid;   // for force residual vector
    Matrix zeroMatrix;
    int SecondaryNodeNum;
    int PrimaryNodeNum;
    int SecondaryDof;
    int PrimaryDof;
    double* restore_shear_gap;
    int loctoglob[6];
};

#endif











#pragma once
