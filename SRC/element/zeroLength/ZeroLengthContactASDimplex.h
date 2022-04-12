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

// $Source: /usr/local/cvs/OpenSees/SRC/element/zeroLength/ZeroLengthContactASDimplex.h,v $
// $Revision: 1.0 $
// $Date: 2020-XX-XX XX:XX:XX $

#ifndef ZeroLengthContactASDimplex_h
#define ZeroLengthContactASDimplex_h

// Written: Onur Deniz Akan         (onur.akan@iusspavia.it)
//          Dr. Massimo Petracca    
//          Prof. Guido Camata      
//          Prof. Enrico Spacone
//          Prof. Carlo G. Lai
//
// Created: May 2020
//

 /*----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*
 |                                                                                |
 |                        ZeroLengthContactASDimplex element                      |
 +                                                                                +
 |--------------------------------------------------------------------------------|
 |                                                                                |
 +         Authors: Onur Deniz Akan (IUSS), Dr. Massimo Petracca (ASDEA)          +
 |                  Prof. Guido Camata, Prof. Enrico Spacone (UNICH) and          |
 |                  Prof. Carlo G. Lai (UNIPV)                                    |
 |                                                                                |
 +      Istituto Universitario di Studi Superiori di Pavia          (IUSS)        +
 |      Advanced Structural Design & Analysis Software Technology   (ASDEA)       |
 |		Universita degli Studi 'G. d'Annunzio' Chieti - Pescara	    (UNICH)       |
 |      Università degli Studi di Pavia                             (UNIPV)       |
 +			                                                                      +
 |                                                                                |
 |                    Email: onur.akan@iusspavia.it (O.D.A.)                      |
 +                                                                                +
 |  Development History:                                                          |
 |  Created       -- May 2020                                                     |
 +  Final Release -- XXX XXXX                                                     +
 |                                                                                |
 |                                                                                |
 +----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----*/

 /* command

  element zeroLengthContactASDimplex $eleID $mNdID $sNdID $Kn $Kt $mu <-orient $x1 $x2 $x3> <-intType $type>

  where:
   $eleID   : element ID of this contact element
   $mNdID   : master node ID
   $sNdID   : slave node ID
   $Kn      : penalty in the normal direction
   $Kt      : penalty in the tangential direction
   $mu      : friction coefficient [mu = tan(Phi) where Phi is the angle of internal friction]
   -orient  : vector components in global coordinates defining local x-axis in 2D or 3D (optional)
                default = <1,0> in 2D and <1,0,0> in 3D
                i.e. normal vector directed from the master towards the slave node
   -intType : type of integration scheme (optional)
                rFlag = 0 implicit, backward-euler integration scheme (default)
                rFlag = 1 impl-ex, implicit/explicit integration scheme

  Description: This file contains the class definition for ZeroLengthContactASDimplex.
  (1) A ZeroLengthContactASDimplex element is defined by two nodes in R^2 (x,y) or R^3 (x,y,z).
  (2) Penalty (Kn,Kt) method is used to constrain displacements in the normal direction,
      and in the tangential direction along with a Mohr-Coulomb frictional yield surface.
  (3) In 2D and 3D, 2, 4, 6 or 12 DOFs with the possibility of having  different DOFsets at
      each of the two nodes is supported automatically. (Example: Solid <--> Shell contact in 3D)
  (4) Optionally, the normal vector of the contact plane, orienting the contact axis can be
      specified by the user. Otherwise, the global X vector is chosen by default.
  (5) Optionally, the user can prefer using the IMPL-EX scheme over the standard Backward Euler
      return mapping scheme for the implicit-explicit integration of the contact material residual
      and tangent modulus.
  (6) Currently, only small deformations are considered, however an update on large deformations
      may be implemented in the near future.
  (7) In IMPL-EX scheme, the material tangent and residual at the current step (n+1) are computed
      by lineary extrapolating the material parameters commited at steps (n) and (n-1). Then, after
      the convergence criteria is achieved, extrapolated material parameters are corrected with a
      step of traditional implicit computation from the last commited step (n). Finally, corrected
      parameters are saved as the current step (n+1) commited variables.
  (8) In a Newton-Raphson scheme, IMPL-EX integration translates to a step-wise linear solution of the material
      non-linearity, and a symmetric, semi-positive definite tangent tensor regardless of what the analytical
      tangent might be, hence improving the robustness of the solution. The order of accuracy is the same with
      the implicit (backward) euler integration scheme, however the commited error at a time step is larger.
      An automatic time-stepping algorithm should be employed in order to control (a-priori) this error. (Oliver et al, 2008)

   References:
   Oliver, J., Huespe, A.E. and Cante, J.C. "An Implicit/Explicit Integration Scheme to
        Increase Computability of Non-linear Material and Contact/Friction Problems", Comp.
        Methods Appl. Mech. Eng., 197, 1865-1889 (2008)
 */

#include <Element.h>
#include <Matrix.h>
#include <array>

class Node;
class Channel;
class Response;

class ZeroLengthContactASDimplex : public Element {

public:
    class StateVariables {
    public:
        // material strain and stress [(0 = normal) (1 = tangent_1) (2 = tangent_2)]
        Vector eps = Vector(3);             // material strain 
        Vector eps_commit = Vector(3);      // commited material strain
        Vector shear = Vector(2);           // effective shear stress
        Vector shear_commit = Vector(2);    // commited effective shear stress
        double xs = 0.0;                    // equivalent plastic strain
        double xs_commit = 0.0;             // commited equivalent plastic strain
        double rs = 0.0;                    // equivalent shear stress
        double rs_commit = 0.0;             // commited equivalent shear stress
        double rs_commit_old = 0.0;         // commited equivalent shear stress (n-1)
        double cres = 0.0;                  // residual shear stress
        double cres_commit = 0.0;           // commited residual shear stress
        double cres_commit_old = 0.0;       // commited residual shear stress (n-1)
        double PC = 1.0;                    // normal stress compressive projector
        double PC_commit = 1.0;             // commited normal stress compressive projector
        double dtime_n = 0.0;               // time factor
        double dtime_n_commit = 0.0;        // commited time factor
        bool dtime_is_user_defined = false;
        bool dtime_first_set = false;
        // modulus
        Matrix C = Matrix(3, 3);            // tangent modulus matrix
        Vector sig = Vector(3);
        Vector sig_implex = Vector(3); // only for output
        // constructor
        StateVariables() = default;
        StateVariables(const StateVariables&) = default;
        StateVariables& operator = (const StateVariables&) = default;
    };

public:
    // constructor
    ZeroLengthContactASDimplex(int tag, int Nd1, int Nd2,
        double Kn, double Kt, double fs,
        int ndm, bool itype, double xN, double yN, double zN);

    // null constructor & destructor
    ZeroLengthContactASDimplex();
    ~ZeroLengthContactASDimplex();

    // public methods to obtain information about dof & connectivity
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();
    void setDomain(Domain* theDomain);

    // public methods to set the state of the element
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getDamp();
    const Matrix& getMass();

    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);
    int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);
    void Print(OPS_Stream& s, int flag = 0);

    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInformation);
    int updateParameter(int parameterID, double value);

private:
    // element info
    ID connectedExternalNodes = ID(2); // contains the tags of the end nodes
    double Knormal = 0.0; // normal penalty
    double Kfriction = 0.0; // tangential penalty
    double mu = 0.0; // friction coefficient
    int numDIM = 0; // model dimension
    std::array<int, 2> numDOF = { { 0, 0 } }; // no. of DOF at 1st and 2nd node
    bool use_implex = false; // integration type flag
    Vector Xorient = Vector(3); // contact axis orientation in global coords.

    // element pointers
    std::array<Node*, 2> theNodes = { { nullptr, nullptr } }; // node pointers

    // initial gap in global coordinates
    // it includes the initial gap in geometry and the initial displacement
    // (if any)
    Vector gap0 = Vector(3);
    bool gap0_initialized = false;

    // state variables
    StateVariables sv; // element & material state variables

private:
    // compute rotation matrix
    const Matrix& getRotationMatrix33(); // 3-by-3 for operations in the local dof-set
    const Matrix& getRotationMatrix66(); // 6-by-6 for operations in the global dof-set
    // compute local initial gap
    const Vector& getInitialGap(); // returns initial gap in local coordinate system
    // compute strain
    void computeStrain();
    // compute material response
    void updateInternal(bool do_implex, bool do_tangent);
    // compute siffness
    void formStiffnessMatrix(const Matrix& C, Matrix& K);
    // compute B matrix
    const Matrix& theBMatrix();
};
#endif