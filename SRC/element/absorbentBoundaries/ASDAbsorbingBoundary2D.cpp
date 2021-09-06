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

// $Revision: 1.10 $
// $Date: 2021/08/30 22:51:21 $

// Original implementation: Massimo Petracca (ASDEA)

#include "ASDAbsorbingBoundary2D.h"

#include <Domain.h>
#include <Node.h>
#include <ErrorHandler.h>
#include <ElementResponse.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Renderer.h>
#include <Parameter.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

namespace {

    static Vector RVector;

}

void*
OPS_ASDAbsorbingBoundary2D(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDAbsorbingBoundary2D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    const char* descr = "Want: element ASDAbsorbingBoundary2D $tag $n1 $n2 $n3 $n4 $G $v $rho $thickness\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 9) {
        opserr << "ASDAbsorbingBoundary2D ERROR : Few arguments:\n" << descr;
        return 0;
    }

    // int parameters
    int iData[5];
    int numData = 5;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "ASDAbsorbingBoundary2D ERROR: Invalid integer mandatory values: element ASDAbsorbingBoundary2D wants 5 integer parameters\n" << descr;
        return 0;
    }

    // double parameters
    double dData[4];
    numData = 4;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "ASDAbsorbingBoundary2D ERROR: Invalid double mandatory values: element ASDAbsorbingBoundary2D wants 4 double parameters\n" << descr;
        return 0;
    }

    // done
    return new ASDAbsorbingBoundary2D(iData[0], iData[1], iData[2], iData[3], iData[4], dData[0], dData[1], dData[2], dData[3]);
}

ASDAbsorbingBoundary2D::ASDAbsorbingBoundary2D()
	: Element(0, ELE_TAG_ASDAbsorbingBoundary2D)
{
}

ASDAbsorbingBoundary2D::ASDAbsorbingBoundary2D(int tag, int node1, int node2, int node3, int node4, double G, double v, double rho, double thickness)
	: Element(tag, ELE_TAG_ASDAbsorbingBoundary2D)
    , m_G(G)
    , m_v(v)
    , m_rho(rho)
    , m_thickness(thickness)
{
    // save node ids
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
    m_node_ids(2) = node3;
    m_node_ids(3) = node4;

    // initialize node vector
    m_nodes.resize(4, nullptr);
}

ASDAbsorbingBoundary2D::~ASDAbsorbingBoundary2D()
{
}

const char* ASDAbsorbingBoundary2D::getClassType(void) const
{
    return "ASDAbsorbingBoundary2D";
}

void ASDAbsorbingBoundary2D::setDomain(Domain* theDomain)
{
    // Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        for (int i = 0; i < 4; ++i)
            m_nodes[i] = nullptr;
        return;
    }

    // check nodes and mapping
    int local_pos = 0;
    for (std::size_t i = 0; i < m_nodes.size(); ++i) {

        // check node
        int node_id = m_node_ids(static_cast<int>(i));
        Node* node = theDomain->getNode(node_id);
        if (node == nullptr) {
            opserr << "ASDAbsorbingBoundary2D Error in setDomain: node " << node_id << " does not exist in the domain\n";
            exit(-1);
        }

        // store node
        m_nodes[i] = node;

        // check NDM
        int ndm = node->getCrds().Size();
        if (ndm != 2) {
            opserr << "ASDAbsorbingBoundary2D Error in setDomain: Nodes must have 2 dimensions, not " << ndm << "\n";
            exit(-1);
        }

        // check NDF
        int ndf = node->getNumberDOF();
        if (ndf != 2 && ndf != 3) {
            opserr << "ASDAbsorbingBoundary2D Error in setDomain: In 2D only 2 or 3 DOFs are allowed, not " << ndf << "\n";
            exit(-1);
        }

        // set up mapping
        int index = static_cast<int>(i) * 2;
        m_mapping(index) = local_pos; // Ux
        m_mapping(index + 1) = local_pos + 1;// Uy
        local_pos += ndf;
    }
    m_num_dofs = local_pos;

    // check direction 1->2 and length 1-2
    static Vector D12(2);
    D12 = m_nodes[1]->getCrds();
    D12 -= m_nodes[0]->getCrds();
    double L12 = D12.Norm();
    if (L12 < std::numeric_limits<double>::epsilon()) {
        opserr << "ASDAbsorbingBoundary2D Error in setDomain: The distance between nodes 1-2 cannot be zero.\n";
        exit(-1);
    }
    D12 /= L12;
    if (std::abs(D12(0)) > 0.99) {
        m_boundary = Boundary_Vertical;
    }
    else if (std::abs(D12(1)) > 0.99) {
        m_boundary = Boundary_Horizontal;
    }
    else {
        opserr << "ASDAbsorbingBoundary2D Error in setDomain: The direction 1-2 can only be X or Y, not [" << D12(0) << ", " << D12(1) << "].\n";
        exit(-1);
    }

    // check direction 1->4 and length 1->4 (make sure they are orthogonal)
    static Vector D14(2);
    D14 = m_nodes[3]->getCrds();
    D14 -= m_nodes[0]->getCrds();
    double L14 = D14.Norm();
    if (L14 < std::numeric_limits<double>::epsilon()) {
        opserr << "ASDAbsorbingBoundary2D Error in setDomain: The distance between nodes 1-4 cannot be zero.\n";
        exit(-1);
    }
    D14 /= L14;
    if (m_boundary == Boundary_Vertical) {
        if (std::abs(D14(1)) < 0.99) {
            opserr << "ASDAbsorbingBoundary2D Error in setDomain: The direction 1-4 can only be Y, not [" << D14(0) << ", " << D14(1) << "].\n";
            exit(-1);
        }
    }
    else {
        if (std::abs(D14(0)) < 0.99) {
            opserr << "ASDAbsorbingBoundary2D Error in setDomain: The direction 1-4 can only be X, not [" << D14(0) << ", " << D14(1) << "].\n";
            exit(-1);
        }
    }

    // compute variables which define the geometry
    if (m_boundary == Boundary_Vertical) {
        m_lx = L12;
        m_ly = L14;
    }
    else {
        m_lx = L14;
        m_ly = L12;
    }
    m_n = 1;// D12(0) > 0.0 ? 1.0 : -1.0;

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void ASDAbsorbingBoundary2D::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ASDAbsorbingBoundary2D\t" << eleTag << " :";
        for (int i = 0; i < m_node_ids.Size(); ++i)
            s << "\t" << m_node_ids(i);
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDAbsorbingBoundary2D\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "]}";
    }
}

int ASDAbsorbingBoundary2D::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& ASDAbsorbingBoundary2D::getExternalNodes()
{
    return m_node_ids;
}

Node** ASDAbsorbingBoundary2D::getNodePtrs()
{
    return m_nodes.data();
}

int ASDAbsorbingBoundary2D::getNumDOF()
{
    return m_num_dofs;
}

int ASDAbsorbingBoundary2D::revertToLastCommit()
{
    return 0;
}

const Matrix& ASDAbsorbingBoundary2D::getTangentStiff(void)
{
    // initialize matrix
    static Matrix K;
    K.resize(m_num_dofs, m_num_dofs);
    K.Zero();

    // fill K matrix
    if (m_stage == Stage_StaticConstraint) {
        addKPenaltyStage0(K);
    }
    else {
        addKPenaltyStage1(K);
        addKff(K);
        addKffToSoil(K);
    }

    // done
    return K;
}

const Matrix& ASDAbsorbingBoundary2D::getInitialStiff(void)
{
    return getTangentStiff();
}

const Matrix& ASDAbsorbingBoundary2D::getDamp(void)
{
    // initialize matrix
    static Matrix C;
    C.resize(m_num_dofs, m_num_dofs);
    C.Zero();

    // fill the C matrix
    if (m_stage == Stage_Absorbing) {
        addCff(C);
        addClk(C);
    }

    // done
    return C;
}

const Matrix& ASDAbsorbingBoundary2D::getMass(void)
{
    // initialize matrix
    static Matrix M;
    M.resize(m_num_dofs, m_num_dofs);
    M.Zero();

    if (m_stage == Stage_Absorbing) {
        // free-field mass
        addMff(M);
    }

    // done
    return M;
}

int ASDAbsorbingBoundary2D::addInertiaLoadToUnbalance(const Vector& accel)
{
    // we don't need this!
    return 0;
}

const Vector& ASDAbsorbingBoundary2D::getResistingForce()
{
    // initialize vector
    Vector& R = RVector;
    R.resize(m_num_dofs);
    R.Zero();

    // fill R vector
    if (m_stage == Stage_StaticConstraint) {
        addRPenaltyStage0(R);
    }
    else {
        addRPenaltyStage1(R);
        addRff(R);
        addRffToSoil(R);
        addRReactions(R);
    }

    // done
    return R;
}

const Vector& ASDAbsorbingBoundary2D::getResistingForceIncInertia()
{
    // initialize vector with resisting forces
    getResistingForce();
    Vector& R = RVector;

    // done if in stage 0
    if (m_stage == Stage_StaticConstraint)
        return R;

    // add damping term
    addRCff(R);
    addRlk(R);

    // add mass term
    addRMff(R);

    // done
    return R;
}

int ASDAbsorbingBoundary2D::addResistingForceToNodalReaction(int flag)
{
    m_is_computing_reactions = true;
    int result = Element::addResistingForceToNodalReaction(flag);
    m_is_computing_reactions = false;
    return result;
}

int ASDAbsorbingBoundary2D::sendSelf(int commitTag, Channel& theChannel)
{
    // todo
    return -1;
}

int ASDAbsorbingBoundary2D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    // todo
    return -1;
}

int ASDAbsorbingBoundary2D::setParameter(const char** argv, int argc, Parameter& param)
{
    // initial check
    if (argc < 1)
        return -1;

    // update G?
    if (strcmp(argv[0], "G") == 0) {
        return param.addObject(1, this);
    }
    // update stage?
    else if (strcmp(argv[0], "stage") == 0) {
        return param.addObject(2, this);
    }

    // no valid parameter
    return -1;
}

int ASDAbsorbingBoundary2D::updateParameter(int parameterID, Information& info)
{
    switch (parameterID)
    {
    case 1:
        // G
        m_G = info.theDouble;
        return 0;
    case 2:
        // stage
        if (m_stage == Stage_StaticConstraint) {
            int istage = static_cast<int>(info.theDouble);
            if (istage == 1) {
                // update stage
                updateStage();
                // done
                return 0;
            }
            else {
                opserr << "Error in ASDAbsorbingBoundary2D::updateParameter.\n"
                    "Current stage = 0 (Stage_StaticConstraint).\n"
                    "The next stage can only be 1 (Stage_Absorbing), not " << istage << "!\n";
                exit(-1);
            }
        }
        else {
            // We cannot move from the Stage_Absorbing
            opserr << "Error in ASDAbsorbingBoundary2D::updateParameter.\n"
                "Current stage = 1 (Stage_Absorbing).\n"
                "You cannot change the stage at this point!\n";
            exit(-1);
        }
    default:
        return -1;
    }
}

double ASDAbsorbingBoundary2D::penaltyFactor()
{
    // order of magnitude of the (approximate) max stiffness value
    // of an equivalent 2D solid element with this G
    int oom = static_cast<int>(std::round(std::log10(m_G * m_thickness)));
    // ideal power = oom + 8 (8 being half the number of significant digits of a double floating number)
    int p = oom + 8;
    // compute the penalty factor
    return std::pow(10.0, p);
}

void ASDAbsorbingBoundary2D::addKPenaltyStage0(Matrix& K)
{
    // Enforce constraints in stage = 0.
    // In this stage:
    // 1) all nodes are fixed in either X (Vertical boundary) or Y (Horizontal boundary) [SP constraint].
    // 2) nodes not belonging to the soil domain
    //    are equalDOF'ed in either X (Horizontal boundary) or Y (Vertical boundary)
    //    to the corresponding nodes on the soil domain [MP constraint]

    // penalty factor
    double p = penaltyFactor();

    // SP and diagonal part of the MP
    for (int i = 0; i < 4; ++i) {
        int j1 = i * 2;
        int j2 = j1 + 1;
        j1 = m_mapping(j1);
        j2 = m_mapping(j2);
        K(j1, j1) += p;
        K(j2, j2) += p;
    }

    // off-diagonal part of the MP
    int dof = m_boundary == Boundary_Horizontal ? 0 : 1;
    for (int i = 0; i < 2; ++i) {
        int j1 = i * 2;
        int j2 = j1 + 1;
        int q1 = j1 * 2 + dof;
        int q2 = j2 * 2 + dof;
        q1 = m_mapping(q1);
        q2 = m_mapping(q2);
        K(q1, q2) -= p;
        K(q2, q1) -= p;
    }
}

void ASDAbsorbingBoundary2D::addRPenaltyStage0(Vector& R)
{
    // Enforce constraints in stage = 0.
    // In this stage:
    // 1) all nodes are fixed in either X (Vertical boundary) or Y (Horizontal boundary) [SP constraint].
    // 2) nodes not belonging to the soil domain
    //    are equalDOF'ed in either X (Horizontal boundary) or Y (Vertical boundary)
    //    to the corresponding nodes on the soil domain [MP constraint]

    // skip if computing reactions
    if (m_is_computing_reactions)
        return;

    // penalty factor
    double p = penaltyFactor();

    // SP and diagonal part of the MP
    for (int i = 0; i < 4; ++i) {
        const Vector& iU = m_nodes[static_cast<std::size_t>(i)]->getTrialDisp();
        int j1 = i * 2;
        int j2 = j1 + 1;
        j1 = m_mapping(j1);
        j2 = m_mapping(j2);
        R(j1) += p * iU(0);
        R(j2) += p * iU(1);
    }

    // off-diagonal part of the MP
    int dof = m_boundary == Boundary_Horizontal ? 0 : 1;
    for (int i = 0; i < 2; ++i) {
        int j1 = i * 2;
        int j2 = j1 + 1;
        const Vector& U1 = m_nodes[static_cast<std::size_t>(j1)]->getTrialDisp();
        const Vector& U2 = m_nodes[static_cast<std::size_t>(j2)]->getTrialDisp();
        int q1 = j1 * 2 + dof;
        int q2 = j2 * 2 + dof;
        q1 = m_mapping(q1);
        q2 = m_mapping(q2);
        R(q1) -= p * U2(dof);
        R(q2) -= p * U1(dof);
    }
}

void ASDAbsorbingBoundary2D::updateStage()
{
    // update stage
    m_stage = Stage_Absorbing;

    // save current displacement and velocities (in local dofset)
    for (int i = 0; i < 4; ++i) {
        Node* inode = m_nodes[static_cast<std::size_t>(i)];
        const Vector& iU = inode->getTrialDisp();
        const Vector& iV = inode->getTrialVel();
        int j1 = i * 2;
        int j2 = j1 + 1;
        m_U0(j1) = iU(0);
        m_U0(j2) = iU(1);
        m_V0(j1) = iV(0);
        m_V0(j2) = iV(1);
    }

    // save reactions (in local dofset)
    static Vector R;
    R.resize(m_num_dofs);
    R.Zero();
    addRPenaltyStage0(R);
    for (int i = 0; i < 4; ++i) {
        // local dofset index
        int j1 = i * 2;
        int j2 = j1 + 1;
        // global dofset index
        int q1 = m_mapping(j1);
        int q2 = m_mapping(j2);
        m_R0(j1) = -R(q1);
        m_R0(j2) = -R(q2);
    }
}

void ASDAbsorbingBoundary2D::addKPenaltyStage1(Matrix& K)
{
    // Enforce constraints in stage = 1.
    // In this stage the edge opposed to the soil domain:
    // 1) is fixed in both X and Y only on Horizontal boundary [SP constraint].

    // skip vertical boundary
    if (m_boundary == Boundary_Vertical)
        return;

    // penalty factor
    double p = penaltyFactor();

    // SP Ux and Uy for all nodes not in the soil domain
    // todo: only in horizontal boundary (1 and 4)
    for (int i = 0; i < 2; ++i) {
        int j = i * 3; // i=0 -> j=0, i=1 -> j=3
        int j1 = j * 2;
        int j2 = j1 + 1;
        int q1 = m_mapping(j1);
        int q2 = m_mapping(j2);
        K(q1, q1) += p;
        K(q2, q2) += p;
    }
}

void ASDAbsorbingBoundary2D::addRPenaltyStage1(Vector& R)
{
    // Enforce constraints in stage = 1.
    // In this stage the edge opposed to the soil domain:
    // 1) is fixed in both X and Y only on Horizontal boundary [SP constraint].

    // skip vertical boundary
    if (m_boundary == Boundary_Vertical)
        return;

    // skip if computing reactions
    if (m_is_computing_reactions)
        return;

    // penalty factor
    double p = penaltyFactor();

    // SP Ux and Uy for all nodes not in the soil domain
    // todo: only in horizontal boundary (1 and 4)
    for (int i = 0; i < 2; ++i) {
        int j = i * 3; // i=0 -> j=0, i=1 -> j=3
        const Vector& iU = m_nodes[static_cast<std::size_t>(j)]->getTrialDisp();
        int j1 = j * 2;
        int j2 = j1 + 1;
        int q1 = m_mapping(j1);
        int q2 = m_mapping(j2);
        R(q1) += p * (iU(0) - m_U0(j1));
        R(q2) += p * (iU(1) - m_U0(j2));
    }
}

void ASDAbsorbingBoundary2D::addRReactions(Vector& R)
{
    // In stage 1, nodes on the soil domain side were fixed in either X or Y direction.
    // In stage 2 those constraints are removed, so we need to restore
    // those reactions as external forces.

    // add restoring reaction forces as external loads (-)
    for (int i = 0; i < 4; ++i) {
        // local dofset index
        int j1 = i * 2;
        int j2 = j1 + 1;
        // global dofset index
        int q1 = m_mapping(j1);
        int q2 = m_mapping(j2);
        R(q1) -= m_R0(j1);
        R(q2) -= m_R0(j2);
    }
}

void ASDAbsorbingBoundary2D::addMff(Matrix& M, double scale)
{
    // Add the mass terms due to the Free-Field columns, thus only 
    // on vertical boundaries (nodes 1-4).
    // The optional scale parameter can be used while assembling
    // the damping due to the free-field

    // skip horizontal boundary
    if (m_boundary == Boundary_Horizontal)
        return;

    // lumped mass (half mass of the whole element)
    double hm = scale * m_rho * m_thickness * m_lx * m_ly / 2.0;
    for (int i = 0; i < 2; ++i) {
        int j = i * 3; // i=0 -> j=0, i=1 -> j=3
        int j1 = j * 2;
        int j2 = j1 + 1;
        int q1 = m_mapping(j1);
        int q2 = m_mapping(j2);
        M(q1, q1) += hm;
        M(q2, q2) += hm;
    }
}

void ASDAbsorbingBoundary2D::addRMff(Vector& R)
{
    // Add the mass terms due to the Free-Field columns, thus only 
    // on vertical boundaries (nodes 1-4).

    // skip horizontal boundary
    if (m_boundary == Boundary_Horizontal)
        return;

    // lumped mass (half mass of the whole element)
    double hm = m_rho * m_thickness * m_lx * m_ly / 2.0;
    for (int i = 0; i < 2; ++i) {
        int j = i * 3; // i=0 -> j=0, i=1 -> j=3
        const Vector& jA = m_nodes[static_cast<std::size_t>(j)]->getTrialAccel();
        int j1 = j * 2;
        int j2 = j1 + 1;
        int q1 = m_mapping(j1);
        int q2 = m_mapping(j2);
        R(q1) += hm * jA(0);
        R(q2) += hm * jA(1);
    }
}

void ASDAbsorbingBoundary2D::addKff(Matrix& K, double scale)
{
    // Add the stiffness matrix of the free-field column.
    // Only on vertical boundaries on nodes 1-4.
    // The optional scale parameter can be used while assembling
    // the damping due to the free-field

    // skip horizontal boundary
    if (m_boundary == Boundary_Horizontal)
        return;

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // stiffness coefficients
    double b = m_lx;
    double h = m_ly;
    double t = m_thickness;
    double kx = scale * b * mu * t / h;
    double ky = scale * b * t * (lam + 2.0 * mu) / h;

    // fill
    int x1 = m_mapping(0);
    int y1 = m_mapping(1);
    int x4 = m_mapping(6);
    int y4 = m_mapping(7);
    K(x1, x1) += kx;
    K(x4, x4) += kx;
    K(y1, y1) += ky;
    K(y1, y1) += ky;
    K(x1, x4) -= kx;
    K(x4, x1) -= kx;
    K(y1, y4) -= ky;
    K(y4, y1) -= ky;
}

void ASDAbsorbingBoundary2D::addRff(Vector& R)
{
    // Add the stiffness matrix of the free-field column.
    // Only on vertical boundaries on nodes 1-4

    // skip horizontal boundary
    if (m_boundary == Boundary_Horizontal)
        return;

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // stiffness coefficients
    double b = m_lx;
    double h = m_ly;
    double t = m_thickness;
    double kx = b * mu * t / h;
    double ky = b * t * (lam + 2.0 * mu) / h;

    // displacement of nodes 1-4, subtracting initial
    // values at stage = 0
    const Vector& U1 = m_nodes[0]->getTrialDisp();
    const Vector& U4 = m_nodes[3]->getTrialDisp();
    double Ux1 = U1(0) - m_U0(0);
    double Uy1 = U1(1) - m_U0(1);
    double Ux4 = U4(0) - m_U0(6);
    double Uy4 = U4(1) - m_U0(7);

    // forces
    double fx = kx * (Ux1 - Ux4);
    double fy = ky * (Uy1 - Uy4);

    // fill
    int x1 = m_mapping(0);
    int y1 = m_mapping(1);
    int x4 = m_mapping(6);
    int y4 = m_mapping(7);
    R(x1) += fx;
    R(x4) -= fx;
    R(y1) += fy;
    R(y4) -= fy;
}

void ASDAbsorbingBoundary2D::addKffToSoil(Matrix& K)
{
    // Add the stiffness matrix of the forces transfered from the
    // free-field column to the soil domain.
    // Only on vertical boundaries on nodes 1-4

    // skip horizontal boundary
    if (m_boundary == Boundary_Horizontal)
        return;

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // geometry
    double h = -1.0;// m_ly;
    double t = m_thickness;
    double n = m_n;

    // fill
    K(m_mapping(2), m_mapping(1)) += h * lam * n * t / 2.0;
    K(m_mapping(2), m_mapping(7)) += -h * lam * n * t / 2.0;
    K(m_mapping(3), m_mapping(0)) += h * mu * n * t / 2.0;
    K(m_mapping(3), m_mapping(6)) += -h * mu * n * t / 2.0;
    K(m_mapping(4), m_mapping(1)) += h * lam * n * t / 2.0;
    K(m_mapping(4), m_mapping(7)) += -h * lam * n * t / 2.0;
    K(m_mapping(5), m_mapping(0)) += h * mu * n * t / 2.0;
    K(m_mapping(5), m_mapping(6)) += -h * mu * n * t / 2.0;
}

void ASDAbsorbingBoundary2D::addRffToSoil(Vector& R)
{
    // Add the stiffness matrix of the forces transfered from the
    // free-field column to the soil domain.
    // Only on vertical boundaries on nodes 1-4

    // skip horizontal boundary
    if (m_boundary == Boundary_Horizontal)
        return;

    // get velocity vector
    static Vector U;
    U.resize(m_num_dofs);
    U.Zero();
    for (int i = 0; i < 4; ++i) {
        const Vector& iU = m_nodes[static_cast<std::size_t>(i)]->getTrialDisp();
        int j1 = i * 2;
        int j2 = j1 + 1;
        int q1 = m_mapping(j1);
        int q2 = m_mapping(j2);
        U(q1) = iU(0) - m_U0(j1);
        U(q2) = iU(1) - m_U0(j2);
    }

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // geometry
    double h = -1.0;// m_ly;
    double t = m_thickness;
    double n = m_n;

    // fill
    R(m_mapping(2)) += -h * lam * n * t * (-U(m_mapping(1)) + U(m_mapping(7))) / 2.0;
    R(m_mapping(3)) += -h * mu * n * t * (-U(m_mapping(0)) + U(m_mapping(6))) / 2.0;
    R(m_mapping(4)) += -h * lam * n * t * (-U(m_mapping(1)) + U(m_mapping(7))) / 2.0;
    R(m_mapping(5)) += -h * mu * n * t * (-U(m_mapping(0)) + U(m_mapping(6))) / 2.0;
}

void ASDAbsorbingBoundary2D::getDampParam(double& alpha, double& beta)
{
    alpha = alphaM;
    beta = betaK;
    if (beta == 0.0) {
        beta = betaK0;
        if (beta == 0)
            beta = betaKc;
    }
}

void ASDAbsorbingBoundary2D::addCff(Matrix& C)
{
    // skip horizontal boundary
    if (m_boundary == Boundary_Horizontal)
        return;

    // compute alpha and beta parameters
    double alpha, beta;
    getDampParam(alpha, beta);

    // mass proportional term
    if (alpha != 0.0) {
        addMff(C, alpha);
    }

    // stiffness proportional term
    if (beta != 0.0) {
        addKff(C, beta);
    }
}

void ASDAbsorbingBoundary2D::addRCff(Vector& R)
{
    // skip horizontal boundary
    if (m_boundary == Boundary_Horizontal)
        return;

    // compute alpha and beta parameters
    double alpha, beta;
    getDampParam(alpha, beta);
    if (alpha == 0.0 && beta == 0.0)
        return;

    // compute damping matrix of the free-field
    static Matrix C;
    C.resize(m_num_dofs, m_num_dofs);
    C.Zero();

    // mass proportional term
    if (alpha != 0.0) {
        addMff(C, alpha);
    }

    // stiffness proportional term
    if (beta != 0.0) {
        addKff(C, beta);
    }

    // get velocity vector
    static Vector V;
    V.resize(m_num_dofs);
    V.Zero();
    for (int i = 0; i < 4; ++i) {
        const Vector& iV = m_nodes[static_cast<std::size_t>(i)]->getTrialVel();
        int j1 = i * 2;
        int j2 = j1 + 1;
        int q1 = m_mapping(j1);
        int q2 = m_mapping(j2);
        V(q1) = iV(0) - m_V0(j1);
        V(q2) = iV(1) - m_V0(j2);
    }
    R.addMatrixVector(1.0, C, V, 1.0);
}

void ASDAbsorbingBoundary2D::getLKcoeff(double& ap, double& as)
{
    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // wave velocities
    double cp = std::sqrt((lam + 2.0 * mu) / m_rho);
    double cs = std::sqrt(mu / m_rho);

    // coefficients
    double h = m_ly;
    if (m_boundary == Boundary_Horizontal) {
        // if the boundary is horizontal:
        // 1) the element size is lx
        // 2) the cp and cs should be swapped
        h = m_lx;
        double aux = cp;
        cp = cs;
        cs = aux;
    }
    double t = m_thickness;
    ap = -cp * h * m_rho * t / 2.0;
    as = -cs * h * m_rho * t / 2.0;
}

void ASDAbsorbingBoundary2D::addClk(Matrix& C)
{
    // dashpot coefficients
    double ap, as;
    getLKcoeff(ap, as);

    // compute derivatives of dashpot forces between nodes:
    // 1->2
    // 4->3
    C(m_mapping(2), m_mapping(0)) += ap;
    C(m_mapping(2), m_mapping(2)) += -ap;
    C(m_mapping(3), m_mapping(1)) += as;
    C(m_mapping(3), m_mapping(3)) += -as;
    C(m_mapping(4), m_mapping(4)) += -ap;
    C(m_mapping(4), m_mapping(6)) += ap;
    C(m_mapping(5), m_mapping(5)) += -as;
    C(m_mapping(5), m_mapping(7)) += as;
}

void ASDAbsorbingBoundary2D::addRlk(Vector& R)
{
    // get velocity vector
    static Vector V;
    V.resize(m_num_dofs);
    V.Zero();
    for (int i = 0; i < 4; ++i) {
        const Vector& iV = m_nodes[static_cast<std::size_t>(i)]->getTrialVel();
        int j1 = i * 2;
        int j2 = j1 + 1;
        int q1 = m_mapping(j1);
        int q2 = m_mapping(j2);
        V(q1) = iV(0) - m_V0(j1);
        V(q2) = iV(1) - m_V0(j2);
    }

    // compute dashpot forces between nodes:
    // 1->2
    // 4->3
    double ap, as;
    getLKcoeff(ap, as);

    // X1 = 2 = 0 - 2
    // Y1 = 3 = 1 - 3
    R(m_mapping(2)) += ap * (V(m_mapping(0)) - V(m_mapping(2)));
    R(m_mapping(3)) += as * (V(m_mapping(1)) - V(m_mapping(3)));
    // X2 = 4 = 6 - 4
    // Y2 = 5 = 7 - 5
    R(m_mapping(4)) += ap * (V(m_mapping(6)) - V(m_mapping(4)));
    R(m_mapping(5)) += as * (V(m_mapping(7)) - V(m_mapping(5)));
}





