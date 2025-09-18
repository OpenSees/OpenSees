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
// Based on feap fortran element elmt04 by Ushnish Basu

#include "FSIFluidBoundaryElement2D.h"

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
#include <array>

void*
OPS_FSIFluidBoundaryElement2D(void)
{
    const char* descr =
        "Want: element FSIFluidBoundaryElement2D $tag $n1 $n2   $c $alpha $g "
        "<-thickness $thickess>\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 6) {
        opserr << "FSIFluidBoundaryElement2D ERROR : Few arguments:\n" << descr;
        return 0;
    }

    // mandatory int parameters
    int iData[3];
    int numData = 3;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "FSIFluidBoundaryElement2D ERROR: Invalid integer mandatory values: element FSIFluidBoundaryElement2D wants 3 integer parameters\n" << descr;
        return 0;
    }

    // mandatory double parameters
    double dData[3];
    numData = 3;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "FSIFluidBoundaryElement2D ERROR: Invalid double mandatory values: element FSIFluidBoundaryElement2D wants 3 double parameters\n" << descr;
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
    return new FSIFluidBoundaryElement2D(
        iData[0], iData[1], iData[2],
        dData[0], dData[1], dData[2],
        thickness);
}

FSIFluidBoundaryElement2D::FSIFluidBoundaryElement2D()
    : Element(0, ELE_TAG_FSIFluidBoundaryElement2D)
{
}

FSIFluidBoundaryElement2D::FSIFluidBoundaryElement2D(
    int tag,
    int node1,
    int node2,
    double c,
    double alpha,
    double g,
    double thickness)
    : Element(tag, ELE_TAG_FSIFluidBoundaryElement2D)
    , m_c(c)
    , m_alpha(alpha)
    , m_g(g)
    , m_thickness(thickness)
{
    // save node ids
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
}

FSIFluidBoundaryElement2D::~FSIFluidBoundaryElement2D()
{
    // clean up load vectors
    if (m_load)
        delete m_load;
}

const char* FSIFluidBoundaryElement2D::getClassType(void) const
{
    return "FSIFluidBoundaryElement2D";
}

void FSIFluidBoundaryElement2D::setDomain(Domain* theDomain)
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
            opserr << "FSIFluidBoundaryElement2D Error in setDomain: node " << node_id << " does not exist in the domain\n";
            exit(-1);
        }

        // store node
        m_nodes[i] = node;

        // check NDM
        int ndm = node->getCrds().Size();
        if (ndm != 2) {
            opserr << "FSIFluidBoundaryElement2D Error in setDomain: Nodes must have 2 dimensions, not " << ndm << "\n";
            exit(-1);
        }

        // check NDF
        int ndf = node->getNumberDOF();
        if (ndf != 1) {
            opserr << "FSIFluidBoundaryElement2D Error in setDomain: 1 DOF is required, not " << ndf << "\n";
            exit(-1);
        }

    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void FSIFluidBoundaryElement2D::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_FSIFluidBoundaryElement2D\t" << eleTag << " :";
        for (int i = 0; i < m_node_ids.Size(); ++i)
            s << "\t" << m_node_ids(i);
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"FSIFluidBoundaryElement2D\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "]}";
    }
}

int FSIFluidBoundaryElement2D::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& FSIFluidBoundaryElement2D::getExternalNodes()
{
    return m_node_ids;
}

Node** FSIFluidBoundaryElement2D::getNodePtrs()
{
    return m_nodes.data();
}

int FSIFluidBoundaryElement2D::getNumDOF()
{
    return 2;
}

int FSIFluidBoundaryElement2D::revertToLastCommit()
{
    return 0;
}

const Matrix& FSIFluidBoundaryElement2D::getTangentStiff(void)
{
    static Matrix K(2, 2);
    K.Zero();
    return K;
}

const Matrix& FSIFluidBoundaryElement2D::getInitialStiff(void)
{
    return getTangentStiff();
}

const Matrix& FSIFluidBoundaryElement2D::getDamp(void)
{
    static Matrix C(2, 2);
    C.Zero();
    if (m_c != 0.0) {
        double factor = (1.0 - m_alpha) / (1.0 + m_alpha) / m_c;
        const Matrix& NTN = getNTN();
        C.addMatrix(0.0, NTN, factor);
    }
    return C;
}

const Matrix& FSIFluidBoundaryElement2D::getMass(void)
{
    static Matrix M(2, 2);
    M.Zero();
    if (m_g != 0.0) {
        double factor = 1.0 / m_g;
        const Matrix& NTN = getNTN();
        M.addMatrix(0.0, NTN, factor);
    }
    return M;
}

void FSIFluidBoundaryElement2D::zeroLoad()
{
    if (m_load)
        m_load->Zero();
}

int FSIFluidBoundaryElement2D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    opserr << "FSIFluidBoundaryElement2D::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
}

int FSIFluidBoundaryElement2D::addInertiaLoadToUnbalance(const Vector& accel)
{
    // only if mass terms are included
    if (m_g != 0.0) {

        // allocate load vector if necessary
        if (m_load == nullptr)
            m_load = new Vector(2);
        Vector& F = *m_load;

        // get mass matrix
        const Matrix& M = getMass();

        // add -M*R*acc to unbalance
        static Vector RA(2);
        for (int i = 0; i < 2; ++i) {
            const Vector& iRA = m_nodes[static_cast<std::size_t>(i)]->getRV(accel);
            RA(i) = iRA(0);
        }
        F.addMatrixVector(1.0, M, RA, -1.0);

    }

    // done
    return 0;
}

const Vector& FSIFluidBoundaryElement2D::getResistingForce()
{
    // initialize vector
    static Vector R(2);
    R.Zero();

    // no stiffness term

    // subtract external loads if any
    if (m_load)
        R.addVector(1.0, *m_load, -1.0);

    // done
    return R;
}

const Vector& FSIFluidBoundaryElement2D::getResistingForceIncInertia()
{
    // initialize vector
    static Vector R(2);
    R.Zero();

    // no stiffness term

    // subtract external loads if any
    if (m_load)
        R.addVector(1.0, *m_load, -1.0);

    // add damping terms if included
    if (m_c != 0.0) {
        const Matrix& C = getDamp();
        const Vector& V = getV();
        R.addMatrixVector(0.0, C, V, 1.0);
    }

    // add mass terms if included
    if (m_g != 0.0) {
        const Matrix& M = getMass();
        const Vector& A = getA();
        R.addMatrixVector(1.0, M, A, 1.0);
    }

    // done
    return R;
}

int FSIFluidBoundaryElement2D::sendSelf(int commitTag, Channel& theChannel)
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
    // 2 node tags +
    // 1 has_load_flag
    static ID idData(4);
    counter = 0;
    idData(counter++) = getTag();
    for (int i = 0; i < 2; ++i)
        idData(counter++) = m_node_ids(i);
    idData(counter++) = static_cast<int>(has_load);
    res = theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING FSIFluidBoundaryElement2D::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    // 1 c +
    // 1 alpha +
    // 1 g
    // 1 thickness +
    // 2 if has_load else 0
    int NLoad = has_load ? 2 : 0;
    Vector vectData(4 + NLoad);
    counter = 0;
    vectData(counter++) = m_c;
    vectData(counter++) = m_alpha;
    vectData(counter++) = m_g;
    vectData(counter++) = m_thickness;
    if (has_load) {
        for (int i = 0; i < 2; ++i)
            vectData(counter++) = (*m_load)(i);
    }
    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING FSIFluidBoundaryElement2D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // done
    return res;
}

int FSIFluidBoundaryElement2D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // a counter
    int counter;

    // INT data
    // 1 tag +
    // 2 node tags +
    // 1 has_load_flag
    static ID idData(4);
    res = theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING FSIFluidBoundaryElement2D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    counter = 0;
    setTag(idData(counter++));
    for (int i = 0; i < 2; ++i)
        m_node_ids(i) = idData(counter++);
    bool has_load = static_cast<bool>(idData(counter++));

    // create load
    if (has_load) {
        if (m_load == nullptr)
            m_load = new Vector(2);
    }
    else {
        if (m_load) {
            delete m_load;
            m_load = nullptr;
        }
    }

    // DOUBLE data
    // 1 c +
    // 1 alpha +
    // 1 g
    // 1 thickness +
    // 2 if has_load else 0
    int NLoad = has_load ? 2 : 0;
    Vector vectData(4 + NLoad);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING FSIFluidBoundaryElement2D::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }

    counter = 0;
    m_c = vectData(counter++);
    m_alpha = vectData(counter++);
    m_g = vectData(counter++);
    m_thickness = vectData(counter++);
    if (has_load) {
        for (int i = 0; i < 2; ++i)
            (*m_load)(i) = vectData(counter++);
    }

    // done
    return res;
}

const Matrix& FSIFluidBoundaryElement2D::getNTN() const
{
    // domain size
    double x0 = m_nodes[0]->getCrds()(0);
    double y0 = m_nodes[0]->getCrds()(1);
    double x1 = m_nodes[1]->getCrds()(0);
    double y1 = m_nodes[1]->getCrds()(1);
    double Tx = x1 - x0;
    double Ty = y1 - y0;
    double L = std::sqrt(Tx * Tx + Ty * Ty);
    if (L == 0.0) {
        opserr << "ERROR FSIFluidBoundaryElement2D::getS() - Element " << this->getTag() << " - edge is collapsed (zero length)\n";
        exit(-1);
    }

    // integration
    constexpr double GLOC = 0.577350269189626;
    std::array<double, 2> X = { -GLOC, GLOC };
    double dV = L * m_thickness / 2.0;

    // holds int(N'*N)dV
    static Matrix NTN(2, 2);
    NTN.Zero();

    // holds shape functions
    static Vector N(2);

    // gauss loop
    for (int igauss = 0; igauss < 2; ++igauss) {
        double xi = X[igauss];
        // shape functions
        N(0) = 0.5 * (1.0 - xi);
        N(1) = 0.5 * (1.0 + xi);
        // integrate the N'N matrix
        for (int i = 0; i < 2; ++i) {
            double Ni = N(i);
            for (int j = 0; j < 2; ++j) {
                double Nj = N(j);
                NTN(i, j) += Ni * Nj * dV;
            }
        }
    }

    // done
    return NTN;
}

const Vector& FSIFluidBoundaryElement2D::getV() const
{
    static Vector U(2);
    for (int j = 0; j < 2; ++j) {
        const Vector& iu = m_nodes[static_cast<std::size_t>(j)]->getTrialVel();
        U(j) = iu(0);
    }
    return U;
}

const Vector& FSIFluidBoundaryElement2D::getA() const
{
    static Vector U(2);
    for (int j = 0; j < 2; ++j) {
        const Vector& iu = m_nodes[static_cast<std::size_t>(j)]->getTrialAccel();
        U(j) = iu(0);
    }
    return U;
}
