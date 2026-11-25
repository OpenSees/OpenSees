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

// Written: Amin Pakzad (University of Washington)
// Created: October 2026
//
// ---------------------------------------------------------------------------
// Credits / origin of the implementation
// ---------------------------------------------------------------------------
// This element combines the ideas of two existing ASDEA elements and does
// NOT modify either of them:
//
//  (1) ASDEmbeddedNodeElement
//      Original implementation: Massimo Petracca, Guido Camata (ASDEA Software
//      Technology). Provides the idea of embedding a constrained node inside a
//      host triangle (2D) or tetrahedron (3D) through the host shape functions
//      evaluated at the constrained node position, and of removing the initial
//      displacements so the element can be activated in a pre-stressed state.
//
//  (2) ZeroLengthContactASDimplex
//      Implemented by: Onur Deniz Akan (IUSS), Massimo Petracca (ASDEA),
//      Guido Camata, Enrico Spacone (UNICH), Carlo G. Lai (UNIPV).
//      Akan, O.D., Petracca, M., Camata, G., Spacone, E. & Lai, C.G. (2020).
//      Provides the penalty frictional (Mohr-Coulomb) contact law with
//      compressive-only normal response and the optional IMPL-EX integration
//      (Oliver et al., 2008).
//
// In EmbeddedNodeContact the "master" node of ZeroLengthContactASDimplex is
// replaced by a virtual point interpolated inside the host element, while the
// "slave" node is the embedded (constrained) node. Hence the embedded node can
// stick, slide, or separate from the same material point of the host element
// without the need of an extra node or of a stiff tie element.
// ---------------------------------------------------------------------------
//
// Command:
//
//   2D: element EmbeddedNodeContact $tag $Cnode $R1 $R2 $R3      $Kn $Kt $mu <-orient $x $y>    <-intType $t>
//   3D: element EmbeddedNodeContact $tag $Cnode $R1 $R2 $R3 $R4  $Kn $Kt $mu <-orient $x $y $z> <-intType $t>
//
//   $Cnode     : embedded (constrained / slave) node
//   $R1..$R3/4 : retained nodes of the host triangle (2D) or tetrahedron (3D)
//   $Kn, $Kt   : normal and tangential penalty stiffness (force/length, NOT scaled)
//   $mu        : friction coefficient
//   -orient    : contact normal in global coordinates, pointing from the host
//                (master point) towards the embedded node. Default = global X.
//   -intType   : 0 = implicit backward Euler (default), 1 = IMPL-EX.
//                (-int_type is accepted as an alias)
//
// Notes / limitations:
//   - Only translational DOFs are coupled. Nodes may have any DOF set supported
//     by the model dimension (2D: 2 or 3; 3D: 3, 4 or 6), also mixed.
//   - Small sliding: the host element and the shape functions are computed
//     once, at setDomain, from the reference configuration.
//   - Local contact strain:  eps = T * ( u_C - sum_i N_i u_Ri + gap0 )
//     eps(0) > 0 means opening (no normal force).

#ifndef EmbeddedNodeContact_h
#define EmbeddedNodeContact_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <vector>

class Node;
class Channel;
class Response;

class EmbeddedNodeContact : public Element
{
public:
    // contact material state variables (same layout as ZeroLengthContactASDimplex)
    class StateVariables {
    public:
        // [(0 = normal) (1 = tangent_1) (2 = tangent_2)]
        Vector eps = Vector(3);
        Vector eps_commit = Vector(3);
        Vector shear = Vector(2);
        Vector shear_commit = Vector(2);
        double xs = 0.0;
        double xs_commit = 0.0;
        double rs = 0.0;
        double rs_commit = 0.0;
        double rs_commit_old = 0.0;
        double cres = 0.0;
        double cres_commit = 0.0;
        double cres_commit_old = 0.0;
        double PC = 1.0;
        double PC_commit = 1.0;
        double dtime_n = 0.0;
        double dtime_n_commit = 0.0;
        bool dtime_is_user_defined = false;
        bool dtime_first_set = false;
        Matrix C = Matrix(3, 3);
        Vector sig = Vector(3);
        Vector sig_implex = Vector(3);
        StateVariables() = default;
        StateVariables(const StateVariables&) = default;
        StateVariables& operator = (const StateVariables&) = default;
    };

public:
    // life cycle
    EmbeddedNodeContact();
    EmbeddedNodeContact(int tag, int ndm, const ID& nodes,
        double Kn, double Kt, double mu, bool use_implex,
        double xN, double yN, double zN);
    virtual ~EmbeddedNodeContact();

    // domain
    const char* getClassType(void) const;
    void setDomain(Domain* theDomain);

    // print
    void Print(OPS_Stream& s, int flag = 0);

    // nodes and dofs
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();

    // state
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // matrices
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getMass();
    const Matrix& getDamp();

    // forces
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // parallel / database
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    // output
    int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);
    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInformation);

private:
    void computeRotationMatrix();
    bool computeShapeFunctions();
    void computeBMatrix();
    void computeStrain();
    void updateInternal(bool do_implex, bool do_tangent);
    const Vector& getGlobalDisplacements() const;

private:
    // nodes: [Cnode, R1, R2, R3, (R4)]
    ID m_node_ids;
    std::vector<Node*> m_nodes;
    // model dimension (2 or 3)
    int m_ndm = 0;
    // total number of element dofs and first dof of each node
    int m_num_dofs = 0;
    ID m_dof_offset;
    // contact parameters
    double m_Kn = 0.0;
    double m_Kt = 0.0;
    double m_mu = 0.0;
    bool m_use_implex = false;
    Vector m_orient = Vector(3);
    // host shape functions at the embedded node (computed at setDomain)
    Vector m_N;
    // rotation matrix (rows = normal, tangent 1, tangent 2)
    Matrix m_T = Matrix(3, 3);
    // global-to-local contact strain operator (3 x num_dofs)
    Matrix m_B;
    // initial gap (global) removing initial displacements
    Vector m_gap0 = Vector(3);
    bool m_gap0_initialized = false;
    // contact state
    StateVariables sv;
};

#endif // EmbeddedNodeContact_h
