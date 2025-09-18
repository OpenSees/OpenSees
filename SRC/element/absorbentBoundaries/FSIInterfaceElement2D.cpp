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
// Based on feap fortran element elmt41 by Ushnish Basu

#include "FSIInterfaceElement2D.h"

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
OPS_FSIInterfaceElement2D(void)
{
    const char* descr = 
        "Want: element FSIInterfaceElement2D $tag $n1 $n2   $rho "
        "<-thickness $thickess>\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 4) {
        opserr << "FSIInterfaceElement2D ERROR : Few arguments:\n" << descr;
        return 0;
    }

    // mandatory int parameters
    int iData[3];
    int numData = 3;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "FSIInterfaceElement2D ERROR: Invalid integer mandatory values: element FSIInterfaceElement2D wants 3 integer parameters\n" << descr;
        return 0;
    }

    // mandatory double parameters
    double dData[1];
    numData = 1;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "FSIInterfaceElement2D ERROR: Invalid double mandatory values: element FSIInterfaceElement2D wants 1 double parameter\n" << descr;
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
    return new FSIInterfaceElement2D(
        iData[0], iData[1], iData[2], 
        dData[0],
        thickness);
}

FSIInterfaceElement2D::FSIInterfaceElement2D()
    : Element(0, ELE_TAG_FSIInterfaceElement2D)
{
}

FSIInterfaceElement2D::FSIInterfaceElement2D(
    int tag,
    int node1,
    int node2,
    double rho,
    double thickness)
    : Element(tag, ELE_TAG_FSIInterfaceElement2D)
    , m_rho(rho)
    , m_thickness(thickness)
{
    // save node ids
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
}

FSIInterfaceElement2D::~FSIInterfaceElement2D()
{
    // clean up load vectors
    if (m_load)
        delete m_load;
}

const char* FSIInterfaceElement2D::getClassType(void) const
{
    return "FSIInterfaceElement2D";
}

void FSIInterfaceElement2D::setDomain(Domain* theDomain)
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
            opserr << "FSIInterfaceElement2D Error in setDomain: node " << node_id << " does not exist in the domain\n";
            exit(-1);
        }

        // store node
        m_nodes[i] = node;

        // check NDM
        int ndm = node->getCrds().Size();
        if (ndm != 2) {
            opserr << "FSIInterfaceElement2D Error in setDomain: Nodes must have 2 dimensions, not " << ndm << "\n";
            exit(-1);
        }

        // check NDF
        int ndf = node->getNumberDOF();
        if (ndf != 3) {
            opserr << "FSIInterfaceElement2D Error in setDomain: 3 DOFs are required, not " << ndf << "\n";
            exit(-1);
        }

    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void FSIInterfaceElement2D::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_FSIInterfaceElement2D\t" << eleTag << " :";
        for (int i = 0; i < m_node_ids.Size(); ++i)
            s << "\t" << m_node_ids(i);
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"FSIInterfaceElement2D\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "]}";
    }
}

int FSIInterfaceElement2D::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& FSIInterfaceElement2D::getExternalNodes()
{
    return m_node_ids;
}

Node** FSIInterfaceElement2D::getNodePtrs()
{
    return m_nodes.data();
}

int FSIInterfaceElement2D::getNumDOF()
{
    return 6;
}

int FSIInterfaceElement2D::revertToLastCommit()
{
    return 0;
}

const Matrix& FSIInterfaceElement2D::getTangentStiff(void)
{
    static Matrix K(6, 6);
    K.Zero();
    const Matrix& S = getS();

    // put -S^T in K
    for (int n = 0; n < 2; ++n) {
        for (int dof = 0; dof < 2; ++dof) {
            int jloc = n * 2 + dof;
            int j = n * 3 + dof; // mapped to U DOFs
            for (int iloc = 0; iloc < 2; ++iloc) {
                int i = iloc * 3 + 2; // mapped to P DOFs
                K(j, i) = -S(iloc, jloc);
            }
        }
    }

    return K;
}

const Matrix& FSIInterfaceElement2D::getInitialStiff(void)
{
    return getTangentStiff();
}

const Matrix& FSIInterfaceElement2D::getDamp(void)
{
    // initialize matrix
    static Matrix C(6, 6);
    C.Zero();
    return C;
}

const Matrix& FSIInterfaceElement2D::getMass(void)
{
    static Matrix M(6, 6);
    M.Zero();
    const Matrix& S = getS();

    // put rho*S in the U DOFs
    for (int n = 0; n < 2; ++n) {
        for (int dof = 0; dof < 2; ++dof) {
            int jloc = n * 2 + dof;
            int j = n * 3 + dof; // mapped to U DOFs
            for (int iloc = 0; iloc < 2; ++iloc) {
                int i = iloc * 3 + 2; // mapped to P DOFs
                M(i, j) = m_rho * S(iloc, jloc);
            }
        }
    }

    return M;
}

void FSIInterfaceElement2D::zeroLoad()
{
    if (m_load)
        m_load->Zero();
}

int FSIInterfaceElement2D::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    opserr << "FSIInterfaceElement2D::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
}

int FSIInterfaceElement2D::addInertiaLoadToUnbalance(const Vector& accel)
{
    // allocate load vector if necessary
    if (m_load == nullptr)
        m_load = new Vector(6);
    Vector& F = *m_load;

    // get mass matrix
    const Matrix& M = getMass();

    // add -M*R*acc to unbalance
    static Vector RA(6);
    for (int i = 0; i < 2; ++i) {
        const Vector& iRA = m_nodes[static_cast<std::size_t>(i)]->getRV(accel);
        for (int j = 0; j < 3; ++j) {
            RA(i * 3 + j) = iRA(j);
        }
    }
    F.addMatrixVector(1.0, M, RA, -1.0);

    // done
    return 0;
}

const Vector& FSIInterfaceElement2D::getResistingForce()
{
    // initialize vector
    static Vector R(6);

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

const Vector& FSIInterfaceElement2D::getResistingForceIncInertia()
{
    // initialize vector
    static Vector R(6);

    // add stiffness terms
    const Matrix& K = getTangentStiff();
    const Vector& U = getU();
    R.addMatrixVector(0.0, K, U, 1.0);

    // subtract external loads if any
    if (m_load)
        R.addVector(1.0, *m_load, -1.0);

    // no damping term

    // add mass terms
    const Matrix& M = getMass();
    const Vector& A = getA();
    R.addMatrixVector(1.0, M, A, 1.0);

    // done
    return R;
}

int FSIInterfaceElement2D::sendSelf(int commitTag, Channel& theChannel)
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
        opserr << "WARNING FSIInterfaceElement2D::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    // 1 rho +
    // 1 thickness +
    // 6 if has_load else 0
    int NLoad = has_load ? 6 : 0;
    Vector vectData(1 + 1 + NLoad);
    counter = 0;
    vectData(counter++) = m_rho;
    vectData(counter++) = m_thickness;
    if (has_load) {
        for (int i = 0; i < 6; ++i)
            vectData(counter++) = (*m_load)(i);
    }
    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING FSIInterfaceElement2D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // done
    return res;
}

int FSIInterfaceElement2D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
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
        opserr << "WARNING FSIInterfaceElement2D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
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
            m_load = new Vector(6);
    }
    else {
        if (m_load) {
            delete m_load;
            m_load = nullptr;
        }
    }

    // DOUBLE data
    // 1 rho +
    // 1 thickness +
    // 6 if has_load else 0
    int NLoad = has_load ? 6 : 0;
    Vector vectData(1 + 1 + NLoad);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING FSIInterfaceElement2D::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }

    counter = 0;
    m_rho = vectData(counter++);
    m_thickness = vectData(counter++);
    if (has_load) {
        for (int i = 0; i < 6; ++i)
            (*m_load)(i) = vectData(counter++);
    }

    // done
    return res;
}

const Matrix& FSIInterfaceElement2D::getS() const
{
    // fluid outward normal vector = -cross(Z, dx)
    double x0 = m_nodes[0]->getCrds()(0);
    double y0 = m_nodes[0]->getCrds()(1);
    double x1 = m_nodes[1]->getCrds()(0);
    double y1 = m_nodes[1]->getCrds()(1);
    double Tx = x1 - x0;
    double Ty = y1 - y0;
    double L = std::sqrt(Tx * Tx + Ty * Ty);
    if (L == 0.0) {
        opserr << "ERROR FSIInterfaceElement2D::getS() - Element " << this->getTag() << " - edge is collapsed (zero length)\n";
        exit(-1);
    }
    static Matrix nt(1, 2);
    nt(0, 0) = Ty/L;
    nt(0, 1) = -Tx/L;

    // coupling matrix (already scaled by the interface area)
    // S = int( Nf' * nt * Ns ) * dV
    static Matrix S(2, 4);
    S.Zero();
    // aux matrix to hold the nt*Ns product
    static Matrix ntNs(1, 4);
    // fluid shape functions
    static Matrix Nf(1, 2);
    // structure shape functions
    static Matrix Ns(2, 4);

    // integration
    constexpr double GLOC = 0.577350269189626;
    std::array<double, 2> X = { -GLOC, GLOC };
    double dV = L * m_thickness / 2.0;

    // gauss loop
    for (int igauss = 0; igauss < 2; ++igauss) {
        double xi = X[igauss];
        // shape functions
        double N1 = 0.5 * (1.0 - xi);
        double N2 = 0.5 * (1.0 + xi);
        // for fluid (P DOF)
        Nf(0, 0) = N1; Nf(0, 1) = N2;
        // for structure (Ux,Uy DOFs)
        Ns(0, 0) = N1; Ns(0, 1) = 0.; Ns(0, 2) = N2; Ns(0, 3) = 0.;
        Ns(1, 0) = 0.; Ns(1, 1) = N1; Ns(1, 2) = 0.; Ns(1, 3) = N2;
        // integrate the S coupling matrix
        //S.addMatrixTripleProduct(1.0, Nf, nt, Ns, dV); // there's some bug in here
        ntNs.addMatrixProduct(0.0, nt, Ns, dV);
        S.addMatrixTransposeProduct(1.0, Nf, ntNs, 1.0);
    }

    // done
    return S;
}

const Vector& FSIInterfaceElement2D::getU() const
{
    static Vector U(6);
    for (int i = 0; i < 2; ++i) {
        const Vector& iu = m_nodes[static_cast<std::size_t>(i)]->getTrialDisp();
        for (int j = 0; j < 3; ++j)
            U(i * 3 + j) = iu(j);
    }
    return U;
}

const Vector& FSIInterfaceElement2D::getA() const
{
    static Vector U(6);
    for (int i = 0; i < 2; ++i) {
        const Vector& iu = m_nodes[static_cast<std::size_t>(i)]->getTrialAccel();
        for (int j = 0; j < 3; ++j)
            U(i * 3 + j) = iu(j);
    }
    return U;
}
