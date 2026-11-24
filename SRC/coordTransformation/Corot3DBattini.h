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
**   Gregory lo. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Purpose: Corot3DBattini is a 3d corotational coordinate transformation for
// two-node frame elements in which the corotated element frame is defined as
// in Battini & Pacoste (2002):
//
//   r1 = deformed chord,
//   q  = (RtI y0 + RtJ y0)/2  (mean of the two nodal y-axes),
//   r3 = r1 x q / |r1 x q|,  r2 = r3 x r1.
//
// The frame is symmetric in the two nodes (results do not depend on the
// element node order). Local nodal rotations are the true logarithm of
// R_r^T R_i (valid up to pi), and the global tangent is the exact consistent
// linearization of the internal force vector (material + geometric parts).
//
// Global rotational dof increments are interpreted as spatial (left) spins
// and accumulated per object in nodal quaternions, as in CorotCrdTransf3d.
//
// Basic system (same as all OpenSees 3d frame transformations):
//   ub = [ ln - lo,  thetaZ_I, thetaZ_J,  thetaY_I, thetaY_J,  thetaX_J - thetaX_I ]
//
// References:
//   Battini, J.-M. and Pacoste, C. (2002) "Co-rotational beam elements with
//   warping effects in instability problems", CMAME 191(17-18):1755-1789.
//   Battini, J.-M. (2002) "Co-rotational beam elements in instability
//   problems", Doctoral thesis, KTH Stockholm: local frame Eqs. (4.25)-(4.28),
//   local rotations Eqs. (4.30)-(4.31).
//   Perez, C.M. and Filippou, F.C. (2024) "On nonlinear geometric
//   transformations of finite elements", IJNME 125(17).
//
// Rigid joint offsets (-jntOffset), [B] Section 4.2.4:
//   The offsets voI, voJ (global, undeformed, from node to element end, as in
//   LinearCrdTransf3d) are rigid links that rotate with their node:
//   end position x_end = x_node + v, v = R_g vo, with R_g = Rt Ro^T the nodal
//   rotation since the start. The element (frame, ub, B, tangent) is formed at
//   the element ends and transformed to the nodes with
//     du_end = du_node - S(v) dw_node   (A = [I -S(v); 0 I] per end)
//     p_node = A^T p_end,   K_node = A^T K_end A + K_off,
//     K_off (rotation block of each node) = v f^T - (v.f) I,  f = end force.
//   These are exact for large rotations (not the small-rotation offsets of the
//   2d transformations).
//
// Limitations:
//   - element loads p0 are applied in the initial local frame
//   - getBasicTrialAccel() returns B*a (the dB/dt*v term is neglected)
//
// Tangent symmetry:
//   The consistent tangent is in general NOT symmetric (rotation increments
//   are spatial spins; see [B] after Eq. 4.94). Use a nonsymmetric solver
//   (e.g. BandGeneral, UmfPack, SparseGeneral). With the optional flag
//   -symmetric the tangent is replaced by its symmetric part (K + K^T)/2;
//   following Nour-Omid & Rankin (1991) this keeps quadratic convergence
//   unless concentrated moments are applied, but it does not make the
//   tangent positive definite (near buckling an SPD solver can still fail).

#ifndef Corot3DBattini_h
#define Corot3DBattini_h

#include <CrdTransf.h>
#include <Vector.h>
#include <Matrix.h>

class Corot3DBattini: public CrdTransf
{
public:
    Corot3DBattini(int tag, const Vector &vecInLocXZPlane,
                   const Vector &rigJntOffsetI, const Vector &rigJntOffsetJ,
                   bool symmetricTangent = false);
    Corot3DBattini();
    ~Corot3DBattini();

    const char *getClassType() const {return "Corot3DBattini";};

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

    int getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);
    int getRigidOffsets(Vector &offsets);

private:
    int  computeInitialGeometry(void);
    int  computeState(void);              // frame, local rotations, ub, GTg, B
    void getNodalDisp(Node *node, const double *initDisp, double *d) const;
    void computeGeometricStiffness(const Vector &pb, double K[12][12]) const;

    // nodes
    Node *nodeIPtr, *nodeJPtr;

    // input
    double vecxz[3];                      // vector in local x-z plane

    // undeformed geometry
    double Ro[3][3];                      // columns are the local axes
    double lo;                            // undeformed length

    // nodal rotations: quaternions (w, x, y, z), trial and committed
    double quatRtI[4], quatRtJ[4];
    double quatRtICommit[4], quatRtJCommit[4];
    double alphaI[3], alphaJ[3];          // last trial rotational dofs seen by update()

    // current configuration (computed in update)
    double RtI[3][3], RtJ[3][3];          // nodal triads (R_g Ro in [B])
    double Rr[3][3];                      // corotated frame, columns r1, r2, r3
    double ln;                            // deformed length
    double qI[3], qJ[3], q[3];            // nodal y-axes and their mean
    double thetaBarI[3], thetaBarJ[3];    // local nodal rotations, log(Rr^T Rt)
    double GTg[3][12];                    // frame spin per end dof (global comps)
    double B[6][12];                      // d(ub)/d(element-end dofs)

    // basic displacements
    Vector ub, ubCommit, ubPrev;

    // initial displacements (treated as stress free)
    double *nodeIInitialDisp, *nodeJInitialDisp;
    bool initialDispChecked;

    // set by recvSelf() so that initialize() keeps the restored state
    bool stateRestored;

    // return the symmetric part of the tangent, (K + K^T)/2
    bool symmetricTangent;

    // initial-configuration B (for getInitialGlobalStiffMatrix)
    double B0[6][12];

    // rigid joint offsets [B 4.2.4]: initial (input) and current, global
    double voI[3], voJ[3];
    double vI[3], vJ[3];
    bool   hasOffsets;
    // end forces from element loads p0 (set in getGlobalResistingForce and
    // used for the offset geometric term K_off in getGlobalStiffMatrix)
    double pEndP0[12];

    // transform end quantities to the nodes for rigid offsets
    void endToNodeVector(const double *pEnd, double *pNode) const;
    void endToNodeMatrix(double K[12][12], const double *pEnd) const;
    void nodeToEnd(const Vector &rI, const Vector &rJ, double *rEnd) const;

    // work areas
    static Matrix Kg;                     // 12x12 global stiffness
    static Vector Pg;                     // 12 global force
    static Matrix Tlg;                    // 12x12 local->global rotation
};

#endif
