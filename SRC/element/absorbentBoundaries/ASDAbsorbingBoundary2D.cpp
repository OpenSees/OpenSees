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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

void*
OPS_ASDAbsorbingBoundary2D(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDAbsorbingBoundary2D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    const char* descr = "Want: element ASDAbsorbingBoundary2D $tag $n1 $n2 $n3 $n4 $G $rho $thickness\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 8) {
        opserr << "ASDAbsorbingBoundary2D ERROR : Few arguments:\n" << descr;
        return 0;
    }

    // int parameters
    int iData[5];
    int numData = 5;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "ASDAbsorbingBoundary2D ERROR: Invalid integer mandatory values: element ASDEmbeddedNodeElement wants 5 integer parameters\n" << descr;
        return 0;
    }

    // double parameters
    double dData[3];
    numData = 3;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "ASDAbsorbingBoundary2D ERROR: Invalid double mandatory values: element ASDEmbeddedNodeElement wants 3 double parameters\n" << descr;
        return 0;
    }

    // done
    return new ASDAbsorbingBoundary2D(iData[0], iData[1], iData[2], iData[3], iData[4], dData[0], dData[1], dData[2]);
}

ASDAbsorbingBoundary2D::ASDAbsorbingBoundary2D()
	: Element(0, ELE_TAG_ASDAbsorbingBoundary2D)
{
}

ASDAbsorbingBoundary2D::ASDAbsorbingBoundary2D(int tag, int node1, int node2, int node3, int node4, double G, double rho, double thickness)
	: Element(tag, ELE_TAG_ASDAbsorbingBoundary2D)
    , m_G(G)
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
    m_n = D12(0) > 0.0 ? 1.0 : -1.0;

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
        addKPenalty(K);
    }
    else {
        // todo...
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

    // todo

    // done
    return C;
}

const Matrix& ASDAbsorbingBoundary2D::getMass(void)
{
    // initialize matrix
    static Matrix M;
    M.resize(m_num_dofs, m_num_dofs);
    M.Zero();

    // todo

    // done
    return M;
}

int ASDAbsorbingBoundary2D::addInertiaLoadToUnbalance(const Vector& accel)
{
    // todo
    return 0;
}

const Vector& ASDAbsorbingBoundary2D::getResistingForce()
{
    // initialize vector
    static Vector R;
    R.resize(m_num_dofs);
    R.Zero();

    // fill R vector
    if (m_stage == Stage_StaticConstraint) {
        addRPenalty(R);
    }
    else {

        // todo...
    }
    

    // done
    return R;
}

const Vector& ASDAbsorbingBoundary2D::getResistingForceIncInertia()
{
    // todo
    return getResistingForce();
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
    // todo
    return -1;
}

int ASDAbsorbingBoundary2D::updateParameter(int parameterID, Information& info)
{
    // todo
    return -1;
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

void ASDAbsorbingBoundary2D::addKPenalty(Matrix& K)
{
    // penalty factor
    double p = penaltyFactor();

    // SP and diagonal part of the MP
    for (int i = 0; i < 4; ++i) {
        int j = i * 2;
        K(j, j) += p;
        K(j + 1, j + 1) += p;
    }

    // off-diagonal part of the MP
    int dof = m_boundary == Boundary_Horizontal ? 0 : 1;
    for (int i = 0; i < 2; ++i) {
        int j1 = i * 2;
        int j2 = j1 + 1;
        int q1 = j1 * 2 + dof;
        int q2 = j2 * 2 + dof;
        K(q1, q2) -= p;
        K(q2, q1) -= p;
    }
}

void ASDAbsorbingBoundary2D::addRPenalty(Vector& R)
{
    // penalty factor
    double p = penaltyFactor();

    // SP and diagonal part of the MP
    for (int i = 0; i < 4; ++i) {
        const Vector& iU = m_nodes[static_cast<std::size_t>(i)]->getTrialDisp();
        int j = i * 2;
        R(j) += p * iU(0);
        R(j + 1) += p * iU(1);
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
        R(q1) -= p * U2(dof);
        R(q2) -= p * U1(dof);
    }

}



