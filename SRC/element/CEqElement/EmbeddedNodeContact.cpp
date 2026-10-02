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
// EmbeddedNodeContact: frictional contact between an embedded node and the
// material point of a host triangle (2D) / tetrahedron (3D) where it lies.
//
// The implementation takes its ideas from two ASDEA elements, which are NOT
// modified by this file:
//   - ASDEmbeddedNodeElement      (Massimo Petracca, Guido Camata - ASDEA)
//       -> host shape-function interpolation of the embedded node and
//          removal of initial displacements.
//   - ZeroLengthContactASDimplex  (Onur Deniz Akan - IUSS, Massimo Petracca - ASDEA,
//                                  Guido Camata, Enrico Spacone - UNICH,
//                                  Carlo G. Lai - UNIPV)
//       -> penalty Mohr-Coulomb frictional contact law with compressive-only
//          normal response and optional IMPL-EX integration. The contact law
//          below (updateInternal) follows the one of ZeroLengthContactASDimplex.
// See EmbeddedNodeContact.h for the command syntax and conventions.

#include <EmbeddedNodeContact.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <ElementResponse.h>
#include <Renderer.h>
#include <elementAPI.h>

#include <cmath>
#include <cstring>
#include <limits>
#include <algorithm>

namespace
{
    inline void cross3(const Vector& A, const Vector& B, Vector& C) {
        C(0) = A(1) * B(2) - A(2) * B(1);
        C(1) = A(2) * B(0) - A(0) * B(2);
        C(2) = A(0) * B(1) - A(1) * B(0);
    }
}

void*
OPS_EmbeddedNodeContact(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using EmbeddedNodeContact - Amin Pakzad (UW). Based on ASDEmbeddedNodeElement (Petracca, Camata - ASDEA) "
            "and ZeroLengthContactASDimplex (Akan, Petracca, Camata, Spacone, Lai - 2020)\n";
        first_done = true;
    }

    const char* descr =
        "Want: element EmbeddedNodeContact $tag $Cnode $R1 $R2 $R3 <$R4 (3D only)> $Kn $Kt $mu "
        "<-orient $x $y <$z>> <-intType $type>\n";

    int ndm = OPS_GetNDM();
    if (ndm != 2 && ndm != 3) {
        opserr << "EmbeddedNodeContact ERROR: Unsupported NDM (" << ndm << "). It should be 2 or 3\n";
        return 0;
    }

    // 2D: constrained + 3 retained (triangle); 3D: constrained + 4 retained (tetrahedron)
    int num_nodes = (ndm == 2) ? 4 : 5;
    if (OPS_GetNumRemainingInputArgs() < 1 + num_nodes + 3) {
        opserr << "EmbeddedNodeContact ERROR: Few arguments:\n" << descr;
        return 0;
    }

    int tag;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) {
        opserr << "EmbeddedNodeContact ERROR: invalid element tag\n" << descr;
        return 0;
    }

    int inodes[5];
    numData = num_nodes;
    if (OPS_GetIntInput(&numData, inodes) != 0) {
        opserr << "EmbeddedNodeContact ERROR: element " << tag << " invalid node tags (expected "
            << num_nodes << " integers in " << ndm << "D)\n" << descr;
        return 0;
    }
    ID nodes(num_nodes);
    for (int i = 0; i < num_nodes; ++i)
        nodes(i) = inodes[i];

    double ddata[3];
    numData = 3;
    if (OPS_GetDoubleInput(&numData, ddata) != 0) {
        opserr << "EmbeddedNodeContact ERROR: element " << tag << " invalid Kn Kt mu\n" << descr;
        return 0;
    }
    if (ddata[0] <= 0.0 || ddata[1] <= 0.0 || ddata[2] < 0.0) {
        opserr << "EmbeddedNodeContact ERROR: element " << tag << " requires Kn > 0, Kt > 0, mu >= 0\n";
        return 0;
    }

    double orient[3] = { 1.0, 0.0, 0.0 };
    int integration_type = 0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* what = OPS_GetString();
        if (what == nullptr) {
            opserr << "EmbeddedNodeContact ERROR: element " << tag << " unexpected (non string) optional argument\n" << descr;
            return 0;
        }
        if (strcmp(what, "-orient") == 0) {
            if (OPS_GetNumRemainingInputArgs() < ndm) {
                opserr << "EmbeddedNodeContact ERROR: element " << tag << " -orient needs " << ndm << " values\n";
                return 0;
            }
            numData = ndm;
            orient[2] = 0.0;
            if (OPS_GetDoubleInput(&numData, orient) != 0) {
                opserr << "EmbeddedNodeContact ERROR: element " << tag << " invalid -orient values\n";
                return 0;
            }
        }
        else if (strcmp(what, "-intType") == 0 || strcmp(what, "-int_type") == 0) {
            numData = 1;
            if (OPS_GetIntInput(&numData, &integration_type) != 0 || (integration_type != 0 && integration_type != 1)) {
                opserr << "EmbeddedNodeContact ERROR: element " << tag << " -intType must be 0 (implicit) or 1 (IMPL-EX)\n";
                return 0;
            }
        }
        else {
            opserr << "EmbeddedNodeContact ERROR: element " << tag << " unknown option " << what << "\n" << descr;
            return 0;
        }
    }

    double norm = std::sqrt(orient[0] * orient[0] + orient[1] * orient[1] + orient[2] * orient[2]);
    if (norm < 1.0e-12) {
        opserr << "EmbeddedNodeContact ERROR: element " << tag << " -orient vector cannot be zero\n";
        return 0;
    }

    return new EmbeddedNodeContact(tag, ndm, nodes, ddata[0], ddata[1], ddata[2], integration_type == 1,
        orient[0] / norm, orient[1] / norm, orient[2] / norm);
}

EmbeddedNodeContact::EmbeddedNodeContact()
    : Element(0, ELE_TAG_EmbeddedNodeContact)
{
}

EmbeddedNodeContact::EmbeddedNodeContact(int tag, int ndm, const ID& nodes,
    double Kn, double Kt, double mu, bool use_implex,
    double xN, double yN, double zN)
    : Element(tag, ELE_TAG_EmbeddedNodeContact)
    , m_node_ids(nodes)
    , m_nodes(static_cast<std::size_t>(nodes.Size()), nullptr)
    , m_ndm(ndm)
    , m_Kn(Kn)
    , m_Kt(Kt)
    , m_mu(mu)
    , m_use_implex(use_implex)
{
    m_orient(0) = xN;
    m_orient(1) = yN;
    m_orient(2) = zN;
}

EmbeddedNodeContact::~EmbeddedNodeContact()
{
}

const char* EmbeddedNodeContact::getClassType(void) const
{
    return "EmbeddedNodeContact";
}

void EmbeddedNodeContact::setDomain(Domain* theDomain)
{
    if (theDomain == nullptr) {
        for (auto& n : m_nodes)
            n = nullptr;
        return;
    }

    int num_nodes = m_node_ids.Size();
    if ((m_ndm == 2 && num_nodes != 4) || (m_ndm == 3 && num_nodes != 5)) {
        opserr << "EmbeddedNodeContact Error in setDomain: element " << getTag()
            << " wrong number of nodes for ndm = " << m_ndm << "\n";
        exit(-1);
    }

    m_nodes.resize(static_cast<std::size_t>(num_nodes), nullptr);
    m_dof_offset.resize(num_nodes);
    m_num_dofs = 0;
    for (int i = 0; i < num_nodes; ++i) {
        Node* node = theDomain->getNode(m_node_ids(i));
        if (node == nullptr) {
            opserr << "EmbeddedNodeContact Error in setDomain: element " << getTag()
                << " node " << m_node_ids(i) << " does not exist in the domain\n";
            exit(-1);
        }
        if (node->getCrds().Size() != m_ndm) {
            opserr << "EmbeddedNodeContact Error in setDomain: element " << getTag()
                << " node " << m_node_ids(i) << " has " << node->getCrds().Size()
                << " coordinates, expected " << m_ndm << "\n";
            exit(-1);
        }
        int ndf = node->getNumberDOF();
        bool ok = (m_ndm == 2) ? (ndf == 2 || ndf == 3) : (ndf == 3 || ndf == 4 || ndf == 6);
        if (!ok) {
            opserr << "EmbeddedNodeContact Error in setDomain: element " << getTag()
                << " node " << m_node_ids(i) << " has unsupported ndf = " << ndf << "\n";
            exit(-1);
        }
        m_nodes[static_cast<std::size_t>(i)] = node;
        m_dof_offset(i) = m_num_dofs;
        m_num_dofs += ndf;
    }

    // geometry: shape functions at the embedded node, rotation, B operator
    if (!computeShapeFunctions()) {
        opserr << "EmbeddedNodeContact Error in setDomain: element " << getTag()
            << " degenerate host element (zero area/volume)\n";
        exit(-1);
    }
    computeRotationMatrix();
    computeBMatrix();

    // initial gap (global): geometric gap minus initial relative displacement.
    // the geometric gap is zero by construction (linear interpolation reproduces
    // the embedded node position), but it is computed for generality.
    if (!m_gap0_initialized) {
        m_gap0.Zero();
        const Vector& XC = m_nodes[0]->getCrds();
        const Vector& UC = m_nodes[0]->getTrialDisp();
        for (int j = 0; j < m_ndm; ++j)
            m_gap0(j) = XC(j) - UC(j);
        for (int i = 1; i < num_nodes; ++i) {
            const Vector& XR = m_nodes[static_cast<std::size_t>(i)]->getCrds();
            const Vector& UR = m_nodes[static_cast<std::size_t>(i)]->getTrialDisp();
            double Ni = m_N(i - 1);
            for (int j = 0; j < m_ndm; ++j)
                m_gap0(j) -= Ni * (XR(j) - UR(j));
        }
        m_gap0_initialized = true;
    }

    DomainComponent::setDomain(theDomain);
}

bool EmbeddedNodeContact::computeShapeFunctions()
{
    // Linear triangle (2D) or tetrahedron (3D) shape functions at the embedded
    // node, as in ASDEmbeddedNodeElement:
    //   X_C = X_1 + J * xi,  J = [X_2 - X_1, X_3 - X_1, (X_4 - X_1)]
    //   N = [1 - sum(xi), xi_1, xi_2, (xi_3)]
    int nr = m_ndm + 1;
    m_N.resize(nr);
    Matrix J(m_ndm, m_ndm);
    Vector D(m_ndm);
    Vector xi(m_ndm);
    const Vector& X1 = m_nodes[1]->getCrds();
    const Vector& XC = m_nodes[0]->getCrds();
    for (int k = 0; k < m_ndm; ++k) {
        const Vector& Xk = m_nodes[static_cast<std::size_t>(k + 2)]->getCrds();
        for (int j = 0; j < m_ndm; ++j)
            J(j, k) = Xk(j) - X1(j);
    }
    for (int j = 0; j < m_ndm; ++j)
        D(j) = XC(j) - X1(j);

    double detJ = (m_ndm == 2)
        ? J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0)
        : J(0, 0) * (J(1, 1) * J(2, 2) - J(1, 2) * J(2, 1))
        - J(0, 1) * (J(1, 0) * J(2, 2) - J(1, 2) * J(2, 0))
        + J(0, 2) * (J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0));
    // relative tolerance on the jacobian
    double scale = 0.0;
    for (int j = 0; j < m_ndm; ++j)
        for (int k = 0; k < m_ndm; ++k)
            scale = std::max(scale, std::fabs(J(j, k)));
    if (scale <= 0.0 || std::fabs(detJ) <= 1.0e-12 * std::pow(scale, m_ndm))
        return false;

    if (J.Solve(D, xi) < 0)
        return false;

    double sum = 0.0;
    for (int j = 0; j < m_ndm; ++j) {
        m_N(j + 1) = xi(j);
        sum += xi(j);
    }
    m_N(0) = 1.0 - sum;

    // warn if the embedded node is outside the host element
    constexpr double tol = 1.0e-6;
    for (int i = 0; i < nr; ++i) {
        if (m_N(i) < -tol || m_N(i) > 1.0 + tol) {
            opserr << "EmbeddedNodeContact WARNING: element " << getTag() << " node " << m_node_ids(0)
                << " lies outside its host element (N = " << m_N(i) << "); extrapolation is used\n";
            break;
        }
    }
    return true;
}

void EmbeddedNodeContact::computeRotationMatrix()
{
    // same local frame construction as ZeroLengthContactASDimplex:
    // row 0 = normal (orient), rows 1-2 = tangents.
    // in 2D, orient lies in the XY plane, so tangent 1 is in-plane and
    // tangent 2 is the out-of-plane Z axis (unused).
    static Vector gY(3);
    static Vector gZ(3);
    gY.Zero(); gY(1) = 1.0;
    gZ.Zero(); gZ(2) = 1.0;
    Vector rY(3);
    Vector rZ(3);
    if (std::fabs(m_orient ^ gY) < 0.99) {
        cross3(m_orient, gY, rZ);
        rZ.Normalize();
        cross3(rZ, m_orient, rY);
        rY.Normalize();
    }
    else {
        cross3(m_orient, gZ, rY);
        rY.Normalize();
        cross3(rY, m_orient, rZ);
        rZ.Normalize();
    }
    for (int j = 0; j < 3; ++j) {
        m_T(0, j) = m_orient(j);
        m_T(1, j) = rY(j);
        m_T(2, j) = rZ(j);
    }
}

void EmbeddedNodeContact::computeBMatrix()
{
    // local contact strain:
    //   eps = T * (u_C - sum_i N_i * u_Ri) + T * gap0
    // only translational dofs (first ndm dofs of each node) are coupled
    m_B.resize(3, m_num_dofs);
    m_B.Zero();
    int num_nodes = m_node_ids.Size();
    for (int i = 0; i < num_nodes; ++i) {
        double factor = (i == 0) ? 1.0 : -m_N(i - 1);
        int off = m_dof_offset(i);
        for (int r = 0; r < 3; ++r)
            for (int j = 0; j < m_ndm; ++j)
                m_B(r, off + j) = factor * m_T(r, j);
    }
}

void EmbeddedNodeContact::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"EmbeddedNodeContact\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "], ";
        s << "\"Kn\": " << m_Kn << ", \"Kt\": " << m_Kt << ", \"mu\": " << m_mu << ", ";
        s << "\"orient\": [" << m_orient(0) << ", " << m_orient(1) << ", " << m_orient(2) << "], ";
        s << "\"intType\": " << (m_use_implex ? 1 : 0) << "}";
        return;
    }
    s << "Element: " << getTag() << " type: EmbeddedNodeContact  Cnode: " << m_node_ids(0) << " Rnodes:";
    for (int i = 1; i < m_node_ids.Size(); ++i)
        s << " " << m_node_ids(i);
    s << "  Kn: " << m_Kn << " Kt: " << m_Kt << " mu: " << m_mu
        << " orient: " << m_orient(0) << " " << m_orient(1) << " " << m_orient(2)
        << " intType: " << (m_use_implex ? 1 : 0) << endln;
    if (m_N.Size() > 0) {
        s << "  shape functions at Cnode:";
        for (int i = 0; i < m_N.Size(); ++i)
            s << " " << m_N(i);
        s << endln;
    }
}

int EmbeddedNodeContact::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& EmbeddedNodeContact::getExternalNodes()
{
    return m_node_ids;
}

Node** EmbeddedNodeContact::getNodePtrs()
{
    return m_nodes.data();
}

int EmbeddedNodeContact::getNumDOF()
{
    return m_num_dofs;
}

int EmbeddedNodeContact::commitState()
{
    // implicit correction step for IMPL-EX (as in ZeroLengthContactASDimplex)
    if (m_use_implex)
        updateInternal(false, false);

    sv.eps_commit = sv.eps;
    sv.shear_commit = sv.shear;
    sv.xs_commit = sv.xs;
    sv.rs_commit_old = sv.rs_commit;
    sv.rs_commit = sv.rs;
    sv.cres_commit_old = sv.cres_commit;
    sv.cres_commit = sv.cres;
    sv.PC_commit = sv.PC;
    sv.dtime_n_commit = sv.dtime_n;
    return 0;
}

int EmbeddedNodeContact::revertToLastCommit()
{
    sv.eps = sv.eps_commit;
    sv.shear = sv.shear_commit;
    sv.xs = sv.xs_commit;
    sv.rs = sv.rs_commit;
    sv.cres = sv.cres_commit;
    sv.PC = sv.PC_commit;
    sv.dtime_n = sv.dtime_n_commit;
    return 0;
}

int EmbeddedNodeContact::revertToStart()
{
    sv = StateVariables();
    return 0;
}

int EmbeddedNodeContact::update()
{
    if (!sv.dtime_is_user_defined) {
        sv.dtime_n = ops_Dt;
        if (!sv.dtime_first_set) {
            sv.dtime_n_commit = sv.dtime_n;
            sv.dtime_first_set = true;
        }
    }
    computeStrain();
    if (m_use_implex) {
        updateInternal(true, true);
        sv.sig_implex = sv.sig;
    }
    else {
        // implicit: numerical tangent (same choice as ZeroLengthContactASDimplex)
        static Vector strain(3);
        static Matrix Cnum(3, 3);
        constexpr double pert = 1.0e-9;
        strain = sv.eps;
        for (int j = 0; j < 3; ++j) {
            sv.eps(j) = strain(j) + pert;
            updateInternal(true, false);
            for (int i = 0; i < 3; ++i)
                Cnum(i, j) = sv.sig(i);
            sv.eps(j) = strain(j) - pert;
            updateInternal(true, false);
            for (int i = 0; i < 3; ++i)
                Cnum(i, j) = (Cnum(i, j) - sv.sig(i)) / 2.0 / pert;
            sv.eps(j) = strain(j);
        }
        updateInternal(true, false);
        sv.C = Cnum;
    }
    return 0;
}

const Matrix& EmbeddedNodeContact::getTangentStiff()
{
    static Matrix K;
    K.resize(m_num_dofs, m_num_dofs);
    K.addMatrixTripleProduct(0.0, m_B, sv.C, 1.0);
    return K;
}

const Matrix& EmbeddedNodeContact::getInitialStiff()
{
    static Matrix K0;
    K0.resize(m_num_dofs, m_num_dofs);
    static Matrix C0(3, 3);
    C0.Zero();
    // local initial normal gap
    double gn = 0.0;
    for (int j = 0; j < 3; ++j)
        gn += m_T(0, j) * m_gap0(j);
    if (gn <= 1.0e-10) {
        C0(0, 0) = m_Kn;
        C0(1, 1) = C0(2, 2) = m_Kt;
    }
    K0.addMatrixTripleProduct(0.0, m_B, C0, 1.0);
    return K0;
}

const Matrix& EmbeddedNodeContact::getMass()
{
    static Matrix M;
    M.resize(m_num_dofs, m_num_dofs);
    M.Zero();
    return M;
}

const Matrix& EmbeddedNodeContact::getDamp()
{
    static Matrix D;
    D.resize(m_num_dofs, m_num_dofs);
    D.Zero();
    return D;
}

const Vector& EmbeddedNodeContact::getResistingForce()
{
    static Vector R;
    R.resize(m_num_dofs);
    R.addMatrixTransposeVector(0.0, m_B, sv.sig, 1.0);
    return R;
}

const Vector& EmbeddedNodeContact::getResistingForceIncInertia()
{
    return getResistingForce();
}

const Vector& EmbeddedNodeContact::getGlobalDisplacements() const
{
    static Vector U;
    U.resize(m_num_dofs);
    int counter = 0;
    for (Node* node : m_nodes) {
        const Vector& iu = node->getTrialDisp();
        for (int i = 0; i < iu.Size(); ++i)
            U(counter++) = iu(i);
    }
    return U;
}

void EmbeddedNodeContact::computeStrain()
{
    const Vector& U = getGlobalDisplacements();
    sv.eps.addMatrixVector(0.0, m_B, U, 1.0);
    // add local initial gap: T * gap0
    for (int r = 0; r < 3; ++r) {
        double g = 0.0;
        for (int j = 0; j < 3; ++j)
            g += m_T(r, j) * m_gap0(j);
        sv.eps(r) += g;
    }
    // no out-of-plane component in 2D
    if (m_ndm == 2)
        sv.eps(2) = 0.0;
}

void EmbeddedNodeContact::updateInternal(bool do_implex, bool do_tangent)
{
    // contact law taken from ZeroLengthContactASDimplex::updateInternal
    // (Akan, Petracca, Camata, Spacone, Lai 2020).
    // strain layout: [Normal, Tangential1, Tangential2]

    sv.rs = sv.rs_commit;
    sv.xs = sv.xs_commit;
    sv.shear = sv.shear_commit;
    sv.cres = sv.cres_commit;

    // note: as in the parent element, the IMPL-EX extrapolation uses a unit
    // time factor (avoids issues with pseudo-time in continuation methods)
    const double time_factor = 1.0;

    // elastic trial
    double SN = m_Kn * sv.eps(0);
    double T1 = sv.shear(0) + m_Kt * (sv.eps(1) - sv.eps_commit(1));
    double T2 = sv.shear(1) + m_Kt * (sv.eps(2) - sv.eps_commit(2));
    double SS = std::sqrt(T1 * T1 + T2 * T2);

    // residual shear strength
    if (do_implex && m_use_implex) {
        sv.cres = std::max(0.0, sv.cres_commit + time_factor * (sv.cres_commit - sv.cres_commit_old));
    }
    else {
        if (SN < 0.0)
            sv.cres = -m_mu * SN;
        else {
            if (!m_use_implex && sv.eps(0) < 1.0e-6)
                sv.cres = 1.0e-10;
        }
    }

    // shear stress in total stress space
    double SS_bar = SS + sv.xs * m_Kt;

    // equivalent shear stress exceeding cres
    if (do_implex && m_use_implex) {
        sv.rs = sv.rs_commit + time_factor * (sv.rs_commit - sv.rs_commit_old);
    }
    else {
        double rs_trial = SS_bar - sv.cres;
        sv.rs = std::max(sv.rs, rs_trial);
    }

    // plastic multiplier and apparent damage
    double equivalent_total_strain = sv.rs / m_Kt;
    double rs_effective = (equivalent_total_strain - sv.xs) * m_Kt + sv.cres;
    sv.xs = equivalent_total_strain;
    double damage = 0.0;
    if (sv.xs > std::numeric_limits<double>::epsilon()) {
        damage = 1.0;
        if (rs_effective > std::numeric_limits<double>::epsilon())
            damage = 1.0 - sv.cres / rs_effective;
    }

    // compressive projector
    if (do_implex && m_use_implex)
        sv.PC = sv.PC_commit;
    else
        sv.PC = SN <= 0.0 ? 1.0 : 0.0;
    double SC = sv.PC * SN;

    // effective shear
    sv.shear(0) = (1.0 - damage) * T1;
    sv.shear(1) = (1.0 - damage) * T2;

    // nominal stress (force)
    sv.sig(0) = SC;
    sv.sig(1) = sv.shear(0);
    sv.sig(2) = sv.shear(1);

    // tangent (IMPL-EX: symmetric, step-wise linear). The implicit tangent is
    // computed numerically in update().
    if (do_tangent) {
        sv.C.Zero();
        sv.C(0, 0) = sv.PC * m_Kn;
        sv.C(1, 1) = sv.C(2, 2) = (1.0 - damage) * m_Kt;
    }
}

int EmbeddedNodeContact::sendSelf(int commitTag, Channel& theChannel)
{
    int dataTag = this->getDbTag();

    // INT data
    static ID idata(18);
    idata.Zero();
    idata(0) = getTag();
    idata(1) = m_ndm;
    idata(2) = m_node_ids.Size();
    for (int i = 0; i < m_node_ids.Size(); ++i)
        idata(3 + i) = m_node_ids(i);
    idata(8) = m_num_dofs;
    idata(9) = m_use_implex ? 1 : 0;
    idata(10) = sv.dtime_is_user_defined ? 1 : 0;
    idata(11) = sv.dtime_first_set ? 1 : 0;
    idata(12) = m_gap0_initialized ? 1 : 0;
    for (int i = 0; i < m_dof_offset.Size() && i < 5; ++i)
        idata(13 + i) = m_dof_offset(i);
    if (theChannel.sendID(dataTag, commitTag, idata) < 0) {
        opserr << "WARNING EmbeddedNodeContact::sendSelf() - " << getTag() << " failed to send ID\n";
        return -1;
    }

    // DOUBLE data
    static Vector ddata(35);
    ddata.Zero();
    ddata(0) = m_Kn;
    ddata(1) = m_Kt;
    ddata(2) = m_mu;
    for (int j = 0; j < 3; ++j) {
        ddata(3 + j) = m_orient(j);
        ddata(6 + j) = m_gap0(j);
        ddata(9 + j) = sv.eps(j);
        ddata(12 + j) = sv.eps_commit(j);
    }
    ddata(15) = sv.shear(0);
    ddata(16) = sv.shear(1);
    ddata(17) = sv.shear_commit(0);
    ddata(18) = sv.shear_commit(1);
    ddata(19) = sv.xs;
    ddata(20) = sv.xs_commit;
    ddata(21) = sv.rs;
    ddata(22) = sv.rs_commit;
    ddata(23) = sv.rs_commit_old;
    ddata(24) = sv.cres;
    ddata(25) = sv.cres_commit;
    ddata(26) = sv.cres_commit_old;
    ddata(27) = sv.PC;
    ddata(28) = sv.PC_commit;
    ddata(29) = sv.dtime_n;
    ddata(30) = sv.dtime_n_commit;
    for (int i = 0; i < m_N.Size() && i < 4; ++i)
        ddata(31 + i) = m_N(i);
    if (theChannel.sendVector(dataTag, commitTag, ddata) < 0) {
        opserr << "WARNING EmbeddedNodeContact::sendSelf() - " << getTag() << " failed to send Vector\n";
        return -1;
    }
    return 0;
}

int EmbeddedNodeContact::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int dataTag = this->getDbTag();

    static ID idata(18);
    if (theChannel.recvID(dataTag, commitTag, idata) < 0) {
        opserr << "WARNING EmbeddedNodeContact::recvSelf() - failed to receive ID\n";
        return -1;
    }
    setTag(idata(0));
    m_ndm = idata(1);
    int num_nodes = idata(2);
    m_node_ids.resize(num_nodes);
    for (int i = 0; i < num_nodes; ++i)
        m_node_ids(i) = idata(3 + i);
    m_nodes.assign(static_cast<std::size_t>(num_nodes), nullptr);
    m_num_dofs = idata(8);
    m_use_implex = idata(9) == 1;
    sv.dtime_is_user_defined = idata(10) == 1;
    sv.dtime_first_set = idata(11) == 1;
    m_gap0_initialized = idata(12) == 1;
    m_dof_offset.resize(num_nodes);
    for (int i = 0; i < num_nodes; ++i)
        m_dof_offset(i) = idata(13 + i);

    static Vector ddata(35);
    if (theChannel.recvVector(dataTag, commitTag, ddata) < 0) {
        opserr << "WARNING EmbeddedNodeContact::recvSelf() - failed to receive Vector\n";
        return -1;
    }
    m_Kn = ddata(0);
    m_Kt = ddata(1);
    m_mu = ddata(2);
    for (int j = 0; j < 3; ++j) {
        m_orient(j) = ddata(3 + j);
        m_gap0(j) = ddata(6 + j);
        sv.eps(j) = ddata(9 + j);
        sv.eps_commit(j) = ddata(12 + j);
    }
    sv.shear(0) = ddata(15);
    sv.shear(1) = ddata(16);
    sv.shear_commit(0) = ddata(17);
    sv.shear_commit(1) = ddata(18);
    sv.xs = ddata(19);
    sv.xs_commit = ddata(20);
    sv.rs = ddata(21);
    sv.rs_commit = ddata(22);
    sv.rs_commit_old = ddata(23);
    sv.cres = ddata(24);
    sv.cres_commit = ddata(25);
    sv.cres_commit_old = ddata(26);
    sv.PC = ddata(27);
    sv.PC_commit = ddata(28);
    sv.dtime_n = ddata(29);
    sv.dtime_n_commit = ddata(30);
    m_N.resize(m_ndm + 1);
    for (int i = 0; i < m_N.Size(); ++i)
        m_N(i) = ddata(31 + i);
    // T and B are recomputed in setDomain (gap0 is kept, flag is set)
    return 0;
}

int EmbeddedNodeContact::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
    if (m_nodes.empty() || m_nodes[0] == nullptr)
        return 0;
    static Vector v1(3);
    m_nodes[0]->getDisplayCrds(v1, fact);
    return theViewer.drawPoint(v1, 1.0, 10);
}

Response* EmbeddedNodeContact::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    Response* theResponse = nullptr;
    if (argc < 1)
        return theResponse;

    output.tag("ElementOutput");
    output.attr("eleType", "EmbeddedNodeContact");
    output.attr("eleTag", this->getTag());
    for (int i = 0; i < m_node_ids.Size(); ++i) {
        char buf[32];
        snprintf(buf, sizeof(buf), "node%d", i + 1);
        output.attr(buf, m_node_ids(i));
    }

    auto open_gauss = [&output]() {
        output.tag("GaussPoint");
        output.attr("number", 1);
        output.attr("eta", 0.0);
        output.tag("NdMaterialOutput");
        output.attr("classType", 0);
        output.attr("tag", 0);
    };
    auto close_gauss = [&output]() {
        output.endTag(); // NdMaterialOutput
        output.endTag(); // GaussPoint
    };
    const char* gl[3] = { "x", "y", "z" };
    const char* lc[3] = { "N", "Tx", "Ty" };

    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
        strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            int ndf = (i + 1 < m_node_ids.Size() ? m_dof_offset(i + 1) : m_num_dofs) - m_dof_offset(i);
            for (int j = 0; j < ndf; ++j) {
                char buf[32];
                snprintf(buf, sizeof(buf), "P%d_%d", j + 1, i + 1);
                output.tag("ResponseType", buf);
            }
        }
        theResponse = new ElementResponse(this, 1, Vector(m_num_dofs));
    }
    else if (strcmp(argv[0], "displacement") == 0 || strcmp(argv[0], "dispJump") == 0) {
        open_gauss();
        for (int j = 0; j < m_ndm; ++j) {
            char buf[16];
            snprintf(buf, sizeof(buf), "dU%s", gl[j]);
            output.tag("ResponseType", buf);
        }
        close_gauss();
        theResponse = new ElementResponse(this, 2, Vector(m_ndm));
    }
    else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {
        open_gauss();
        for (int j = 0; j < m_ndm; ++j)
            output.tag("ResponseType", lc[j]);
        close_gauss();
        theResponse = new ElementResponse(this, 3, Vector(m_ndm));
    }
    else if (strcmp(argv[0], "localForceImplex") == 0 || strcmp(argv[0], "localForcesImplex") == 0) {
        open_gauss();
        for (int j = 0; j < m_ndm; ++j)
            output.tag("ResponseType", lc[j]);
        close_gauss();
        theResponse = new ElementResponse(this, 33, Vector(m_ndm));
    }
    else if (strcmp(argv[0], "localDisplacement") == 0 || strcmp(argv[0], "localDispJump") == 0) {
        open_gauss();
        output.tag("ResponseType", "dUN");
        output.tag("ResponseType", "dUTx");
        if (m_ndm == 3)
            output.tag("ResponseType", "dUTy");
        close_gauss();
        theResponse = new ElementResponse(this, 4, Vector(m_ndm));
    }
    else if (strcmp(argv[0], "slip") == 0 || strcmp(argv[0], "slipMultiplier") == 0) {
        open_gauss();
        output.tag("ResponseType", "lambda");
        close_gauss();
        theResponse = new ElementResponse(this, 5, Vector(1));
    }
    else if (strcmp(argv[0], "normalContactForce") == 0 || strcmp(argv[0], "NormalContactForce") == 0) {
        open_gauss();
        output.tag("ResponseType", "N");
        close_gauss();
        theResponse = new ElementResponse(this, 6, Vector(1));
    }
    else if (strcmp(argv[0], "tangentialContactForce") == 0 || strcmp(argv[0], "TangentialContactForce") == 0) {
        open_gauss();
        output.tag("ResponseType", "|T|");
        close_gauss();
        theResponse = new ElementResponse(this, 7, Vector(1));
    }
    else if (strcmp(argv[0], "contactState") == 0 || strcmp(argv[0], "state") == 0) {
        // 0 = open, 1 = stick, 2 = slip
        open_gauss();
        output.tag("ResponseType", "state");
        close_gauss();
        theResponse = new ElementResponse(this, 8, Vector(1));
    }
    else if (strcmp(argv[0], "shapeFunctions") == 0 || strcmp(argv[0], "N") == 0) {
        for (int i = 0; i < m_ndm + 1; ++i) {
            char buf[16];
            snprintf(buf, sizeof(buf), "N%d", i + 1);
            output.tag("ResponseType", buf);
        }
        theResponse = new ElementResponse(this, 9, Vector(m_ndm + 1));
    }
    (void)gl;
    output.endTag(); // ElementOutput
    return theResponse;
}

int EmbeddedNodeContact::getResponse(int responseID, Information& eleInfo)
{
    static Vector small;
    static Vector scalar(1);
    small.resize(m_ndm);

    switch (responseID) {
    case 1:
        return eleInfo.setVector(getResistingForce());
    case 2: {
        // global displacement jump: T^T * eps
        for (int i = 0; i < m_ndm; ++i) {
            double v = 0.0;
            for (int r = 0; r < 3; ++r)
                v += m_T(r, i) * sv.eps(r);
            small(i) = v;
        }
        return eleInfo.setVector(small);
    }
    case 3:
        for (int i = 0; i < m_ndm; ++i)
            small(i) = sv.sig(i);
        return eleInfo.setVector(small);
    case 33:
        for (int i = 0; i < m_ndm; ++i)
            small(i) = sv.sig_implex(i);
        return eleInfo.setVector(small);
    case 4:
        for (int i = 0; i < m_ndm; ++i)
            small(i) = sv.eps(i);
        return eleInfo.setVector(small);
    case 5:
        scalar(0) = sv.xs;
        return eleInfo.setVector(scalar);
    case 6:
        scalar(0) = sv.sig(0);
        return eleInfo.setVector(scalar);
    case 7:
        scalar(0) = std::sqrt(sv.sig(1) * sv.sig(1) + sv.sig(2) * sv.sig(2));
        return eleInfo.setVector(scalar);
    case 8: {
        double state = 0.0; // open
        if (sv.PC > 0.5) {
            state = 1.0; // stick
            // slipping if the equivalent plastic shear grew in the last
            // step (rs_commit_old = value at the beginning of that step)
            if (sv.rs > sv.rs_commit_old + 1.0e-12 * std::max(1.0, std::fabs(sv.rs)))
                state = 2.0; // slip
        }
        scalar(0) = state;
        return eleInfo.setVector(scalar);
    }
    case 9:
        return eleInfo.setVector(m_N);
    default:
        return -1;
    }
}
