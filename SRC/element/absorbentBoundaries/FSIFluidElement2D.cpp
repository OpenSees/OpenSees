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
// $Date: 2023/12/18 22:51:21 $

// Implementation: Massimo Petracca (ASDEA)
// Based on feap fortran element elmt03 by Ushnish Basu

#include "FSIFluidElement2D.h"

#include <Domain.h>
#include <Node.h>
#include <ErrorHandler.h>
#include <ElementResponse.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Renderer.h>
#include <Parameter.h>

#include <limits>
#include <cmath>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>

void*
OPS_FSIFluidElement2D(void)
{
    const char* descr = 
        "Want: element FSIFluidElement2D $tag $n1 $n2 $n3 $n4   $c "
        "<-thickness $thickess>\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 6) {
        opserr << "FSIFluidElement2D ERROR : Few arguments:\n" << descr;
        return 0;
    }

    // mandatory int parameters
    int iData[5];
    int numData = 5;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "FSIFluidElement2D ERROR: Invalid integer mandatory values: element FSIFluidElement2D wants 5 integer parameters\n" << descr;
        return 0;
    }

    // mandatory double parameters
    double dData[1];
    numData = 1;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "FSIFluidElement2D ERROR: Invalid double mandatory values: element FSIFluidElement2D wants 1 double parameter\n" << descr;
        return 0;
    }

    // optional
    numData = 1;
    double thickness = 1.0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();
        if (strcmp(type, "-thickness") == 0) {
            if (OPS_GetNumRemainingInputArgs() == 0) {
                opserr << "thickness not found\n";
                return 0;
            }
            if (OPS_GetDoubleInput(&numData, &thickness) < 0)
                return 0;
        }
    }

    // done
    return new FSIFluidElement2D(
        iData[0], iData[1], iData[2], iData[3], iData[4], 
        dData[0],
        thickness);
}

FSIFluidElement2D::FSIFluidElement2D()
    : Element(0, ELE_TAG_FSIFluidElement2D)
{
}

FSIFluidElement2D::FSIFluidElement2D(
    int tag,
    int node1,
    int node2,
    int node3,
    int node4,
    double c,
    double thickness)
    : Element(tag, ELE_TAG_FSIFluidElement2D)
    , m_c(c)
    , m_thickness(thickness)
{
    // save node ids
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
    m_node_ids(2) = node3;
    m_node_ids(3) = node4;
}

FSIFluidElement2D::~FSIFluidElement2D()
{
    // clean up load vectors
    if (m_load)
        delete m_load;
}

const char* FSIFluidElement2D::getClassType(void) const
{
    return "FSIFluidElement2D";
}

void FSIFluidElement2D::setDomain(Domain* theDomain)
{
    // Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        for (std::size_t i = 0; i < m_nodes.size(); ++i)
            m_nodes[i] = nullptr;
        return;
    }

    // check nodes
    for (std::size_t i = 0; i < m_nodes.size(); ++i) {

        // check node
        int node_id = m_node_ids(static_cast<int>(i));
        Node* node = theDomain->getNode(node_id);
        if (node == nullptr) {
            opserr << "FSIFluidElement2D Error in setDomain: node " << node_id << " does not exist in the domain\n";
            exit(-1);
        }

        // store node
        m_nodes[i] = node;

        // check NDM
        int ndm = node->getCrds().Size();
        if (ndm != 2) {
            opserr << "FSIFluidElement2D Error in setDomain: Nodes must have 2 dimensions, not " << ndm << "\n";
            exit(-1);
        }

        // check NDF
        int ndf = node->getNumberDOF();
        if (ndf != 1) {
            opserr << "FSIFluidElement2D Error in setDomain: 1 DOF is required, not " << ndf << "\n";
            exit(-1);
        }

    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void FSIFluidElement2D::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_FSIFluidElement2D\t" << eleTag << " :";
        for (int i = 0; i < m_node_ids.Size(); ++i)
            s << "\t" << m_node_ids(i);
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"FSIFluidElement2D\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "]}";
    }
}

int FSIFluidElement2D::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& FSIFluidElement2D::getExternalNodes()
{
    return m_node_ids;
}

Node** FSIFluidElement2D::getNodePtrs()
{
    return m_nodes.data();
}

int FSIFluidElement2D::getNumDOF()
{
    return 4;
}

int FSIFluidElement2D::revertToLastCommit()
{
    return 0;
}

const Matrix& FSIFluidElement2D::getTangentStiff(void)
{
    // fluid stiffness
    static Matrix K(4, 4);
    K.Zero();
    static Matrix dNdX(2, 4);
    // gauss loop
    const Matrix& g = q4gauss();
    for (int i = 0; i < g.noRows(); ++i) {
        // initegration
        double xi = g(i, 0);
        double eta = g(i, 1);
        double wi = g(i, 2);
        // shape functions, natural derivatives and jacobian
        const Vector& N = q4N(xi, eta);
        const Matrix& dN = q4dN(xi, eta);
        double detJ;
        const Matrix& invJ = q4Jac(dN, detJ);
        double dV = wi * detJ * m_thickness;
        // shape function cartesian derivatives
        dNdX.addMatrixProduct(0.0, invJ, dN, 1.0);
        // integrate stiffness
        K.addMatrixTransposeProduct(1.0, dNdX, dNdX, dV);
    }
    return K;
}

const Matrix& FSIFluidElement2D::getInitialStiff(void)
{
    return getTangentStiff();
}

const Matrix& FSIFluidElement2D::getDamp(void)
{
    // initialize matrix
    static Matrix C(4, 4);
    C.Zero();
    return C;
}

const Matrix& FSIFluidElement2D::getMass(void)
{
    // fluid mass
    static Matrix M(4, 4);
    M.Zero();
    double c2 = m_c * m_c;
    // gauss loop
    const Matrix& g = q4gauss();
    for (int i = 0; i < g.noRows(); ++i) {
        // initegration
        double xi = g(i, 0);
        double eta = g(i, 1);
        double wi = g(i, 2);
        // shape functions, natural derivatives and jacobian
        const Vector& N = q4N(xi, eta);
        const Matrix& dN = q4dN(xi, eta);
        double detJ;
        const Matrix& invJ = q4Jac(dN, detJ);
        double dV = wi * detJ * m_thickness;
        // integrate Mass: M += N'*N*dV/c^2
        dV /= c2;
        for (int j = 0; j < 4; ++j) {
            double Nj = N(j);
            for (int k = 0; k < 4; ++k) {
                double Nk = N(k);
                M(j, k) += Nj * Nk * dV;
            }
        }
    }
    return M;
}

void FSIFluidElement2D::zeroLoad()
{
    if (m_load)
        m_load->Zero();
}

int FSIFluidElement2D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    opserr << "FSIFluidElement2D::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
}

int FSIFluidElement2D::addInertiaLoadToUnbalance(const Vector& accel)
{
    // allocate load vector if necessary
    if (m_load == nullptr)
        m_load = new Vector(4);
    Vector& F = *m_load;

    // get mass matrix
    const Matrix& M = getMass();

    // add -M*R*acc to unbalance
    static Vector RA(4);
    for (int i = 0; i < 4; ++i) 
        RA(i) = m_nodes[static_cast<std::size_t>(i)]->getRV(accel)(0);
    F.addMatrixVector(1.0, M, RA, -1.0);

    // done
    return 0;
}

const Vector& FSIFluidElement2D::getResistingForce()
{
    // initialize vector
    static Vector R(4);

    // add stiffness terms
    const Matrix& K = getTangentStiff();
    const Vector& U = getU();
    R.addMatrixVector(0.0, K, U, 1.0);

    // subtract external loads if any
    if (m_load) 
        R.addVector(1.0, *m_load, -1.0);

    // done
    return R;
}

const Vector& FSIFluidElement2D::getResistingForceIncInertia()
{
    // initialize vector
    static Vector R(4);

    // add stiffness terms
    const Matrix& K = getTangentStiff();
    const Vector& U = getU();
    R.addMatrixVector(0.0, K, U, 1.0);

    // subtract external loads if any
    if (m_load)
        R.addVector(1.0, *m_load, -1.0);

    // add damping terms
    const Matrix& C = getDamp();
    const Vector& V = getV();
    R.addMatrixVector(1.0, C, V, 1.0);

    // add mass terms
    const Matrix& M = getMass();
    const Vector& A = getA();
    R.addMatrixVector(1.0, M, A, 1.0);

    // done
    return R;
}

int FSIFluidElement2D::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // a counter
    int counter;

    // has load flag
    bool has_load = m_load != nullptr;

    // INT data
    // 1 tag +
    // 4 node tags +
    // 1 has_load_flag
    static ID idData(6);
    counter = 0;
    idData(counter++) = getTag();
    for (int i = 0; i < 4; ++i)
        idData(counter++) = m_node_ids(i);
    idData(counter++) = static_cast<int>(has_load);
    res = theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING FSIFluidElement2D::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    // 1 c +
    // 1 thickness +
    // 4 if has_load else 0
    int NLoad = has_load ? 4 : 0;
    Vector vectData(1 + 1 + NLoad);
    counter = 0;
    vectData(counter++) = m_c;
    vectData(counter++) = m_thickness;
    if (has_load) {
        for (int i = 0; i < 4; ++i)
            vectData(counter++) = (*m_load)(i);
    }
    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING FSIFluidElement2D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // done
    return res;
}

int FSIFluidElement2D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // a counter
    int counter;

    // INT data
    // 1 tag +
    // 4 node tags +
    // 1 has_load_flag
    static ID idData(6);
    res = theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING FSIFluidElement2D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    counter = 0;
    setTag(idData(counter++));
    for (int i = 0; i < 4; ++i)
        m_node_ids(i) = idData(counter++);
    bool has_load = static_cast<bool>(idData(counter++));

    // create load
    if (has_load) {
        if (m_load == nullptr)
            m_load = new Vector(4);
    }
    else {
        if (m_load) {
            delete m_load;
            m_load = nullptr;
        }
    }

    // DOUBLE data
    // 1 c +
    // 1 thickness +
    // 4 if has_load else 0
    int NLoad = has_load ? 4 : 0;
    Vector vectData(1 + 1 + NLoad);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING FSIFluidElement2D::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }

    counter = 0;
    m_c = vectData(counter++);
    m_thickness = vectData(counter++);
    if (has_load) {
        for (int i = 0; i < 4; ++i)
            (*m_load)(i) = vectData(counter++);
    }

    // done
    return res;
}

const Vector& FSIFluidElement2D::getU() const
{
    static Vector U(4);
    for (int j = 0; j < 4; ++j) {
        const Vector& iu = m_nodes[static_cast<std::size_t>(j)]->getTrialDisp();
        U(j) = iu(0);
    }
    return U;
}

const Vector& FSIFluidElement2D::getV() const
{
    static Vector U(4);
    for (int j = 0; j < 4; ++j) {
        const Vector& iu = m_nodes[static_cast<std::size_t>(j)]->getTrialVel();
        U(j) = iu(0);
    }
    return U;
}

const Vector& FSIFluidElement2D::getA() const
{
    static Vector U(4);
    for (int j = 0; j < 4; ++j) {
        const Vector& iu = m_nodes[static_cast<std::size_t>(j)]->getTrialAccel();
        U(j) = iu(0);
    }
    return U;
}

const Matrix& FSIFluidElement2D::q4gauss() const
{
    static Matrix q(4, 3);

    constexpr double GLOC = 0.577350269189626;
    q(0, 0) = -GLOC;  q(0, 1) = -GLOC; q(0, 2) = 1.0;
    q(1, 0) =  GLOC;  q(1, 1) = -GLOC; q(1, 2) = 1.0;
    q(2, 0) =  GLOC;  q(2, 1) =  GLOC; q(2, 2) = 1.0;
    q(3, 0) = -GLOC;  q(3, 1) =  GLOC; q(3, 2) = 1.0;

    return q;
}

const Vector& FSIFluidElement2D::q4N(double xi, double eta) const
{
    static Vector N(4);

    N(0) = 0.25 * (1.0 - xi) * (1.0 - eta);
    N(1) = 0.25 * (1.0 + xi) * (1.0 - eta);
    N(2) = 0.25 * (1.0 + xi) * (1.0 + eta);
    N(3) = 0.25 * (1.0 - xi) * (1.0 + eta);

    return N;
}

const Matrix& FSIFluidElement2D::q4dN(double xi, double eta) const
{
    static Matrix dN(2, 4);

    dN(0, 0) = -(1.0 - eta) * 0.25;
    dN(0, 1) =  (1.0 - eta) * 0.25;
    dN(0, 2) =  (1.0 + eta) * 0.25;
    dN(0, 3) = -(1.0 + eta) * 0.25;

    dN(1, 0) = -(1.0 - xi) * 0.25;
    dN(1, 1) = -(1.0 + xi) * 0.25;
    dN(1, 2) =  (1.0 + xi) * 0.25;
    dN(1, 3) =  (1.0 - xi) * 0.25;

    return dN;
}

const Matrix& FSIFluidElement2D::q4Jac(const Matrix& dN, double& detJ) const
{
    static Matrix J(2, 2);
    static Matrix invJ(2, 2);

    // coordinates
    const Vector& X1 = m_nodes[0]->getCrds();
    const Vector& X2 = m_nodes[1]->getCrds();
    const Vector& X3 = m_nodes[2]->getCrds();
    const Vector& X4 = m_nodes[3]->getCrds();

    // jacobian
    J(0, 0) = dN(0, 0) * X1(0) + dN(0, 1) * X2(0) + dN(0, 2) * X3(0) + dN(0, 3) * X4(0);
    J(1, 0) = dN(0, 0) * X1(1) + dN(0, 1) * X2(1) + dN(0, 2) * X3(1) + dN(0, 3) * X4(1);
    J(0, 1) = dN(1, 0) * X1(0) + dN(1, 1) * X2(0) + dN(1, 2) * X3(0) + dN(1, 3) * X4(0);
    J(1, 1) = dN(1, 0) * X1(1) + dN(1, 1) * X2(1) + dN(1, 2) * X3(1) + dN(1, 3) * X4(1);

    // determinant
    detJ = J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1);
    double mult = 1.0 / detJ;

    // inv(jacobian)
    invJ(0, 0) =  J(1, 1) * mult;
    invJ(1, 1) =  J(0, 0) * mult;
    invJ(1, 0) = -J(0, 1) * mult;
    invJ(0, 1) = -J(1, 0) * mult;

    return invJ;
}
