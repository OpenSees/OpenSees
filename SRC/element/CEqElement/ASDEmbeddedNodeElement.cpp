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
// $Date: 2021/04/28 22:51:21 $

// Original implementation: Massimo Petracca (ASDEA)
//
//

#include <ASDEmbeddedNodeElement.h>

#include <Domain.h>
#include <Node.h>
#include <ErrorHandler.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Renderer.h>
#include <analysis/dof_grp/DOF_Group.h>

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <limits>

// anonymous namespace for utilities
namespace
{

    double det2(const Matrix& J) {
        return J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
    }

    double det3(const Matrix& J) {
        return  J(0, 0) * J(1, 1) * J(2, 2) - J(0, 0) * J(1, 2) * J(2, 1) - J(0, 1) * J(1, 0) * J(2, 2) +
            J(0, 1) * J(1, 2) * J(2, 0) + J(0, 2) * J(1, 0) * J(2, 1) - J(0, 2) * J(1, 1) * J(2, 0);
    }

    void cross(const Vector& a, const Vector& b, Vector& c) {
        c(0) = a(1) * b(2) - a(2) * b(1);
        c(1) = a(2) * b(0) - a(0) * b(2);
        c(2) = a(0) * b(1) - a(1) * b(0);
    }

    namespace tri {

        double shapeFun(double x, double y, int i) {
            if (i == 0)
                return 1.0 - x - y;
            else if (i == 1)
                return x;
            else if (i == 2)
                return y;
            return 0.0;
        }

        void shapeFunDer(Matrix& dN) {
            dN(0, 0) = -1.0; dN(0, 1) = -1.0;
            dN(1, 0) = 1.0; dN(1, 1) = 0.0;
            dN(2, 0) = 0.0; dN(2, 1) = 1.0;
        }

        void globalCoord(const Matrix& X, double lx, double ly, double& gx, double& gy) {
            gx = gy = 0.0;
            for (int i = 0; i < 3; i++) {
                double N = shapeFun(lx, ly, i);
                gx += N * X(0, i);
                gy += N * X(1, i);
            }
        }

        void globalCoord(const Matrix& X, double lx, double ly, double& gx, double& gy, double& gz) {
            gx = gy = gz = 0.0;
            for (int i = 0; i < 3; i++) {
                double N = shapeFun(lx, ly, i);
                gx += N * X(0, i);
                gy += N * X(1, i);
                gz += N * X(2, i);
            }
        }

        void localCoord(const Matrix& X, const Matrix& invJ, double gx, double gy, double& lx, double& ly) {
            lx = ly = 0.0;
            double px, py;
            globalCoord(X, lx, ly, px, py);
            Vector D(2);
            Vector DL(2);
            D(0) = gx - px;
            D(1) = gy - py;
            DL.addMatrixVector(0.0, invJ, D, 1.0);
            lx = DL(0);
            ly = DL(1);
        }

        void localCoord(const Matrix& X, const Matrix& invJ, double gx, double gy, double gz, double& lx, double& ly) {
            lx = ly = 0.0;
            double px, py, pz;
            globalCoord(X, lx, ly, px, py, pz);
            Vector D(3);
            Vector DL(3);
            D(0) = gx - px;
            D(1) = gy - py;
            D(2) = gz - pz;
            DL.addMatrixVector(0.0, invJ, D, 1.0);
            lx = DL(0);
            ly = DL(1);
        }

        void fillVzInJacobian(Matrix& J) {
            double nx = J(1, 0) * J(2, 1) - J(1, 1) * J(2, 0);
            double ny = J(0, 1) * J(2, 0) - J(0, 0) * J(2, 1);
            double nz = J(0, 0) * J(1, 1) - J(0, 1) * J(1, 0);
            double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
            if (norm > std::numeric_limits<double>::epsilon()) {
                J(0, 2) = nx / norm;
                J(1, 2) = ny / norm;
                J(2, 2) = nz / norm;
            }
        }

    }

    namespace tet {

        double shapeFun(double x, double y, double z, int i) {
            if (i == 0)
                return 1.0 - (x + y + z);
            else if (i == 1)
                return x;
            else if (i == 2)
                return y;
            else if (i == 3)
                return z;
            return 0.0;
        }

        void shapeFunDer(Matrix& dN) {
            dN(0, 0) = -1.0; dN(0, 1) = -1.0; dN(0, 2) = -1.0;
            dN(1, 0) = 1.0; dN(1, 1) = 0.0; dN(1, 2) = 0.0;
            dN(2, 0) = 0.0; dN(2, 1) = 1.0; dN(2, 2) = 0.0;
            dN(3, 0) = 0.0; dN(3, 1) = 0.0; dN(3, 2) = 1.0;
        }

        void globalCoord(const Matrix& X, double lx, double ly, double lz, double& gx, double& gy, double& gz) {
            gx = gy = gz = 0.0;
            for (int i = 0; i < 4; i++) {
                double N = shapeFun(lx, ly, lz, i);
                gx += N * X(0, i);
                gy += N * X(1, i);
                gz += N * X(2, i);
            }
        }

        void localCoord(const Matrix& X, const Matrix& invJ, double gx, double gy, double gz, double& lx, double& ly, double& lz) {
            lx = ly = lz = 0.0;
            double px, py, pz;
            globalCoord(X, lx, ly, lz, px, py, pz);
            Vector D(3);
            Vector DL(3);
            D(0) = gx - px;
            D(1) = gy - py;
            D(2) = gz - pz;
            DL.addMatrixVector(0.0, invJ, D, 1.0);
            lx = DL(0);
            ly = DL(1);
            lz = DL(2);
        }

    }

}

void * OPS_ADD_RUNTIME_VPV(OPS_ASDEmbeddedNodeElement)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDEmbeddedNodeElement - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    const char* descr = "Want: element ASDEmbeddedNodeElement $tag $Cnode $Rnode1 $Rnode2 $Rnode3 <$Rnode4> <-rot> <-K $K>\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 5) {
        opserr << "ASDEmbeddedNodeElement ERROR : Few arguments:\n" << descr;
        return 0;
    }
    
    // mandatory parameters
    int iData[5];
    int numData = 5;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "ASDEmbeddedNodeElement ERROR: Invalid integer mandatory values: element ASDEmbeddedNodeElement wants at least 5 integer parameters\n" << descr;
        return 0;
    }

    // parse optional parameters
    bool rot = false;
    int N4 = 0;
    bool has_N4 = false;
    double K = 1.0e18;
    for (int i = 5; i < numArgs; i++) {
        const char* what = OPS_GetString();
        if (strcmp(what, "-rot") == 0) {
            rot = true;
        }
        else if (strcmp(what, "-K") == 0) {
            if (i == numArgs - 1) {
                opserr << "ASDEmbeddedNodeElement ERROR: The -K keyword should be followed by a floating point number.\n" << descr;
                return 0;
            }
            ++i;
            numData = 1;
            if (OPS_GetDouble(&numData, &K) != 0) {
                opserr << "ASDEmbeddedNodeElement ERROR invalid floating point number for -K keyword.\n";
                return 0;
            }
        }
        else {
            // should be an integer for the only if i == 6
            if (i == 5) {
                try {
                    N4 = std::stoi(std::string(what));
                    has_N4 = true;
                }
                catch (...) {
                    N4 = -1;
                    has_N4 = false;
                }
            }
        }
    }

    // done
    if (has_N4)
        return new ASDEmbeddedNodeElement(iData[0], iData[1], iData[2], iData[3], iData[4], N4, rot, K);
    else
        return new ASDEmbeddedNodeElement(iData[0], iData[1], iData[2], iData[3], iData[4], rot, K);
}

ASDEmbeddedNodeElement::ASDEmbeddedNodeElement() 
    : Element(0, ELE_TAG_ASDEmbeddedNodeElement)
{
}

ASDEmbeddedNodeElement::ASDEmbeddedNodeElement(int tag, int cNode, int rNode1, int rNode2, int rNode3, bool rot_flag, double K)
    : Element(tag, ELE_TAG_ASDEmbeddedNodeElement)
    , m_rot_c_flag(rot_flag)
    , m_K(K)
{
    m_node_ids.resize(4);
    m_node_ids(0) = cNode;
    m_node_ids(1) = rNode1;
    m_node_ids(2) = rNode2;
    m_node_ids(3) = rNode3;
    m_nodes.resize(4, nullptr);
}

ASDEmbeddedNodeElement::ASDEmbeddedNodeElement(int tag, int cNode, int rNode1, int rNode2, int rNode3, int rNode4, bool rot_flag, double K)
    : Element(tag, ELE_TAG_ASDEmbeddedNodeElement)
    , m_rot_c_flag(rot_flag)
    , m_K(K)
{
    m_node_ids.resize(5);
    m_node_ids(0) = cNode;
    m_node_ids(1) = rNode1;
    m_node_ids(2) = rNode2;
    m_node_ids(3) = rNode3;
    m_node_ids(4) = rNode4;
    m_nodes.resize(5, nullptr);
}

ASDEmbeddedNodeElement::~ASDEmbeddedNodeElement( )
{
}

const char* ASDEmbeddedNodeElement::getClassType(void) const
{
    return "ASDEmbeddedNodeElement";
}

void ASDEmbeddedNodeElement::setDomain(Domain* theDomain)
{
    // check nodes
    m_num_dofs = 0;
    int local_dof_counter = 0;
    int local_pos = 0;
    std::vector<ID> aux_mapping(m_nodes.size());
    for (std::size_t i = 0; i < m_nodes.size(); ++i) {

        // check node
        int node_id = m_node_ids(static_cast<int>(i));
        Node* node = theDomain->getNode(node_id);
        if (node == nullptr) {
            opserr << "ASDEmbeddedNodeElement Error in setDomain: node " << node_id << " does not exit in the domain\n";
            exit(-1);
        }

        // store node
        m_nodes[i] = node;

        // check NDM
        int ndm = node->getCrds().Size();
        if (ndm != 2 && ndm != 3) {
            opserr << "ASDEmbeddedNodeElement Error in setDomain: Nodes should have either 2 or 3 dimensions, not " << ndm << "\n";
            exit(-1);
        }
        if (i == 0) {
            // save ndm at first node
            m_ndm = ndm;
        }
        else {
            if (m_ndm != ndm) {
                opserr << "ASDEmbeddedNodeElement Error in setDomain: Nodes should have the same dimension (2 or 3)\n";
                exit(-1);
            }
        }

        // check NDF
        int ndf = node->getNumberDOF();
        if (m_ndm == 2) {
            if (ndf != 2 && ndf != 3) {
                opserr << "ASDEmbeddedNodeElement Error in setDomain: In 2D only 2 or 3 DOFs are allowed, not " << ndf << "\n";
                exit(-1);
            }
            if (i == 0) {
                m_rot_c = (m_rot_c_flag && ndf == 3);
            }
        }
        else {
            if (ndf != 3 && ndf != 4 && ndf != 6) {
                opserr << "ASDEmbeddedNodeElement Error in setDomain: In 3D only 3, 4 or 6 DOFs are allowed, not " << ndf << "\n";
                exit(-1);
            }
            if (i == 0) {
                m_rot_c = (m_rot_c_flag && ndf == 6);
            }
        }

        // set up mapping
        ID& imap = aux_mapping[i];
        imap.resize(m_ndm + ((i == 0 && m_rot_c) ? (m_ndm == 2 ? 1 : 3) : 0));
        imap(0) = local_pos; // Ux
        imap(1) = local_pos + 1; // Uy
        if (m_ndm == 3) {
            imap(2) = local_pos + 2; // Uz
            if (i == 0 && m_rot_c) {
                imap(3) = local_pos + 3; // Rx
                imap(4) = local_pos + 4; // Rx
                imap(5) = local_pos + 5; // Rx
            }
        }
        else {
            if (i == 0 && m_rot_c) {
                imap(2) = local_pos + 2; // Rz
            }
        }
        local_pos += ndf;

        // update total dof counter
        m_num_dofs += ndf;
        // update local dof counter
        local_dof_counter += imap.Size();
    }

    // flatten mapping
    m_mapping.resize(local_dof_counter);
    local_pos = 0;
    for (const ID& imap : aux_mapping) {
        for (int i = 0; i < imap.Size(); ++i) {
            m_mapping(local_pos++) = imap(i);
        }
    }

    // compute initial displacement vector
    if (!m_U0_computed) {
        m_U0.resize(m_num_dofs);
        m_U0 = getGlobalDisplacements();
        m_U0_computed = true;
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void ASDEmbeddedNodeElement::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ASDEmbeddedNodeElement\t" << eleTag << " :";
        for (int i = 0; i < m_node_ids.Size(); ++i)
            s << "\t" << m_node_ids(i);
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDEmbeddedNodeElement\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "]}";
    }
}

int ASDEmbeddedNodeElement::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& ASDEmbeddedNodeElement::getExternalNodes()
{
    return m_node_ids;
}

Node**
ASDEmbeddedNodeElement::getNodePtrs(void)
{
    return m_nodes.data();
}

int ASDEmbeddedNodeElement::getNumDOF()
{
    return m_num_dofs;
}

int ASDEmbeddedNodeElement::update()
{
    return 0;
}

int ASDEmbeddedNodeElement::revertToLastCommit()
{
    return 0;
}

const Matrix& ASDEmbeddedNodeElement::getTangentStiff()
{
    // compute stiffness matrix in reduced local dofset
    auto compute_local = [this]() -> const Matrix& {
        if (m_nodes.size() == 4) {
            // support shape is a triangle ...
            if (m_ndm == 2) {
                // ... in 2D
                if (m_rot_c) {
                    // ... with rotational dofs
                    return TRI_2D_UR();
                }
                else {
                    // ... without rotational dofs
                    return TRI_2D_U();
                }
            }
            else {
                // ... in 3D
                if (m_rot_c) {
                    // ... with rotational dofs
                    return TRI_3D_UR();
                }
                else {
                    // ... without rotational dofs
                    return TRI_3D_U();
                }
            }
        }
        else {
            // support shape is a tetrahedron ...
            if (m_rot_c) {
                // ... with rotational dofs
                return TET_3D_UR();
            }
            else {
                // ... without rotational dofs
                return TET_3D_U();
            }
        }
    };
    const Matrix& KL = compute_local();

    // output matrix
    static Matrix K;
    K.resize(m_num_dofs, m_num_dofs);
    K.Zero();

    // copy in global dofset
    for (int i = 0; i < KL.noRows(); ++i) {
        int ig = m_mapping(i);
        for (int j = 0; j < KL.noCols(); ++j) {
            int jg = m_mapping(j);
            K(ig, jg) = KL(i, j);
        }
    }

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::getInitialStiff()
{
    return getTangentStiff();
}

const Matrix& ASDEmbeddedNodeElement::getMass()
{
    static Matrix M;
    M.resize(m_num_dofs, m_num_dofs);
    M.Zero();
    return M;
}

int
ASDEmbeddedNodeElement::addInertiaLoadToUnbalance(const Vector& accel)
{
    return 0;
}

const Vector& ASDEmbeddedNodeElement::getResistingForce()
{
    static Vector F;
    F.resize(m_num_dofs);
    const Matrix& K = getTangentStiff();
    const Vector& U = getGlobalDisplacements();
    F.addMatrixVector(0.0, K, U, 1.0);
    return F;
}

const Vector& ASDEmbeddedNodeElement::getResistingForceIncInertia()
{
    return getResistingForce();
}

int ASDEmbeddedNodeElement::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // INT data
    // tag
    // number of nodes
    // 5 nodal ids (the last one is 0 if number of nodes = 4)
    // ndm
    // num_dofs
    // rot_c_flag
    // rot_c
    // U0_computed
    // number of local dofs (at most 18)
    // mapping of local dofs (if less then 18, the last ones will be set to 0)
    static ID idData(31);
    idData.Zero();
    idData(0) = getTag();
    idData(1) = m_node_ids.Size();
    idData(2) = m_node_ids(0);
    idData(3) = m_node_ids(1);
    idData(4) = m_node_ids(2);
    idData(5) = m_node_ids(3);
    if (m_node_ids.Size() == 5)
        idData(6) = m_node_ids(4);
    idData(7) = m_ndm;
    idData(8) = m_num_dofs;
    idData(9) = m_rot_c_flag ? 1 : 0;
    idData(10) = m_rot_c ? 1 : 0;
    idData(11) = m_U0_computed ? 1 : 0;
    idData(12) = m_mapping.Size();
    for (int i = 0; i < m_mapping.Size(); ++i)
        idData(12 + i) = m_mapping(i);
    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    // K
    // initial displacement (at most we can have 30 dofs)
    static Vector vectData(31);
    vectData.Zero();
    vectData(0) = m_K;
    for (int i = 0; i < m_num_dofs; ++i)
        vectData(1 + i) = m_U0(i);
    res += theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // done
    return res;
}

int ASDEmbeddedNodeElement::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // INT data
    // tag
    // number of nodes
    // 5 nodal ids (the last one is 0 if number of nodes = 4)
    // ndm
    // num_dofs
    // rot_c_flag
    // rot_c
    // number of local dofs (at most 18)
    // mapping of local dofs (if less then 18, the last ones will be set to 0)
    static ID idData(31);
    res += theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    setTag(idData(0));
    int num_nodes = idData(1);
    m_node_ids.resize(num_nodes);
    m_nodes.resize(static_cast<std::size_t>(num_nodes), nullptr);
    m_node_ids(0) = idData(2);
    m_node_ids(1) = idData(3);
    m_node_ids(2) = idData(4);
    m_node_ids(3) = idData(5);
    if (m_node_ids.Size() == 5)
        m_node_ids(4) = idData(6);
    m_ndm = idData(7);
    m_num_dofs = idData(8);
    m_rot_c_flag = idData(9) == 1;
    m_rot_c = idData(10) == 1;
    m_U0_computed = idData(11) == 1;
    int num_local_dofs = idData(12);
    m_mapping.resize(num_local_dofs);
    for (int i = 0; i < m_mapping.Size(); ++i)
        m_mapping(i) = idData(12 + i);

    // DOUBLE data
    // K
    static Vector vectData(31);
    res += theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDEmbeddedNodeElement::sendSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }
    m_K = vectData(0);
    m_U0.resize(m_num_dofs);
    for (int i = 0; i < m_num_dofs; ++i)
        m_U0(i) = vectData(1 + i);

    // done
    return res;
}

const Vector& ASDEmbeddedNodeElement::getGlobalDisplacements() const
{
    static Vector U(m_num_dofs);
    int counter = 0;
    for (Node* node : m_nodes) {
        const Vector& iu = node->getTrialDisp();
        for (int i = 0; i < iu.Size(); ++i) {
            U(counter++) = iu(i);
        }
    }
    if (m_U0_computed) {
        U.addVector(1.0, m_U0, -1.0);
    }
    return U;
}

const Matrix& ASDEmbeddedNodeElement::TRI_2D_U()
{
    // output
    static Matrix K(8, 8);

    // collect triangle coordinates
    static Matrix X(2, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i-1) = node->getCrds()(0);
        X(1, i-1) = node->getCrds()(1);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 2);
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(2, 2);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det2(J);
    double V = detJ / 2.0;
    static Matrix invJ(2, 2);
    J.Invert(invJ);

    // find local coordinates of constrained node
    double lx, ly;
    tri::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for(int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    static Matrix B(2, 8);
    B.Zero();
    for (int i = 0; i < 2; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 3; i++) {
        int j = 2 + i * 2;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
    }

    // Penalty stiffness
    double iK = m_K * V;

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TRI_2D_UR()
{
    // output
    static Matrix K(9, 9);
    
    // collect triangle coordinates
    static Matrix X(2, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 2);
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(2, 2);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det2(J);
    double V = detJ / 2.0;
    static Matrix invJ(2, 2);
    J.Invert(invJ);

    // shape functions cartesian derivatives
    static Matrix dNdX(3, 2);
    dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

    // find local coordinates of constrained node
    double lx, ly;
    tri::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for (int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // RCz = sum(d_URy_dX - d_URx_dY)/2.0 -> sum(d_URy_dX - d_URx_dY)/2.0 - RCz = 0
    static Matrix B(3, 9);
    B.Zero();
    for (int i = 0; i < 3; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 3; i++) {
        int j = 3 + i * 2;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j) = -dNdX(i, 1)/2.0; B(2, j + 1) = dNdX(i, 0)/2.0;
    }

    // Penalty stiffness
    double iK = m_K* V;

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TRI_3D_U()
{
    // output
    static Matrix K(12, 12);

    // collect triangle coordinates
    static Matrix X(3, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 3);
    dN.Zero(); // note: the tri::shapeFunDer fills the first 2 columns, the 3rd should be zero!
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(3, 3);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    tri::fillVzInJacobian(J); // note: the 3rd column will be zero, so fill it with the unit normal vector
    double detJ = det3(J);
    double V = detJ / 2.0;
    static Matrix invJ(3, 3);
    J.Invert(invJ);
    for (int i = 0; i < 3; ++i) invJ(i, 2) = 0.0; // note: the 3rd column in the inverse should be zero

    // find local coordinates of constrained node
    double lx, ly;
    tri::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), m_nodes[0]->getCrds()(2), lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for (int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    static Matrix B(3, 12);
    B.Zero();
    for (int i = 0; i < 3; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 3; i++) {
        int j = 3 + i * 3;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j + 2) = N(i);
    }

    // Penalty stiffness
    double iK = m_K * V;

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TRI_3D_UR()
{
    // output
    static Matrix K(15, 15);

    // collect triangle coordinates
    // in global coordinates
    static Matrix X(3, 3);
    for (int i = 1; i < 4; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // compute orientation
    static Vector dx(3);
    static Vector dy(3);
    static Vector dz(3);
    for (int i = 0; i < 3; ++i) {
        dx(i) = X(i, 1) - X(i, 0);
        dy(i) = X(i, 2) - X(i, 0);
    }
    dx.Normalize();
    dy.Normalize();
    cross(dx, dy, dz);
    dz.Normalize();
    cross(dz, dx, dy);

    // assemble orientation matrix (transposed of rotation)
    static Matrix R(3, 3);
    for (int i = 0; i < 3; ++i) {
        R(0, i) = dx(i);
        R(1, i) = dy(i);
        R(2, i) = dz(i);
    }

    // triangle coordinates in local coordinates
    static Matrix XL(2, 3);
    for (int i = 0; i < 3; ++i) {
        XL(0, i) = X(0, i) * dx(0) + X(1, i) * dx(1) + X(2, i) * dx(2);
        XL(1, i) = X(0, i) * dy(0) + X(1, i) * dy(1) + X(2, i) * dy(2);
    }

    // shape functions natural derivatives
    static Matrix dN(3, 2);
    tri::shapeFunDer(dN);

    // jacobian
    static Matrix J(2, 2);
    J.addMatrixProduct(0.0, XL, dN, 1.0);
    double detJ = det2(J);
    double V = detJ / 2.0;
    static Matrix invJ(2, 2);
    J.Invert(invJ);

    // shape functions cartesian derivatives
    static Matrix dNdX(3, 2);
    dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

    // find local coordinates of constrained node
    const Vector CPos = m_nodes[0]->getCrds();
    double CPosX = CPos(0) * dx(0) + CPos(1) * dx(1) + CPos(2) * dx(2);
    double CPosY = CPos(0) * dy(0) + CPos(1) * dy(1) + CPos(2) * dy(2);
    double lx, ly;
    tri::localCoord(XL, invJ, CPosX, CPosY, lx, ly);

    // compute shape functions at constrained node
    static Vector N(3);
    for (int i = 0; i < 3; ++i)
        N(i) = tri::shapeFun(lx, ly, i);

    // compute B matrix (in local coordinates)
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    // RCx = sum( d_URz_dY) -> sum( d_URz_dY) - RCx = 0 (local)
    // RCy = sum(-d_URz_dX) -> sum(-d_URz_dX) - RCy = 0 (local)
    // RCz = sum(d_URy_dX - d_URx_dY)/2.0 -> sum(d_URy_dX - d_URx_dY)/2.0 - RCz = 0 (local)
    static Matrix B(6, 15);
    B.Zero();
    // fill the -identity 6x6 block (transformed to global coordinates)
    // for constrained node+
    /*for (int i = 0; i < 6; ++i)
        B(i, i) = -1.0;*/
    for (int i = 0; i < 2; ++i) {
        int j = i * 3;
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                B(j + row, j + col) = -R(row, col);
    }
    // fill the 2 rows of 3 3x3 blocks (transformed to global coordinates)
    static Matrix BL(3, 3);
    static Matrix BG(3, 3);
    for (int i = 0; i < 3; ++i) {
        int j = 6 + i * 3;
        // U block
        BL.Zero();
        BL(0, 0) = N(i);
        BL(1, 1) = N(i);
        BL(2, 2) = N(i);
        BG.addMatrixProduct(0.0, BL, R, 1.0);
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                B(row, j + col) = BG(row, col);
        // R block
        BL.Zero();
        BL(0, 2) = dNdX(i, 1);
        BL(1, 2) = -dNdX(i, 0);
        BL(2, 0) = -dNdX(i, 1) / 2.0; 
        BL(2, 1) = dNdX(i, 0) / 2.0;
        BG.addMatrixProduct(0.0, BL, R, 1.0);
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < 3; ++col)
                B(3 + row, j + col) = BG(row, col);
    }

    // Penalty stiffness
    double iK = m_K * V;

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TET_3D_U()
{
    // output
    static Matrix K(15, 15);

    // collect tetrahedron coordinates
    static Matrix X(3, 4);
    for (int i = 1; i < 5; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // shape functions natural derivatives
    static Matrix dN(4, 3);
    tet::shapeFunDer(dN);

    // jacobian
    static Matrix J(3, 3);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det3(J);
    double V = detJ / 6.0;
    static Matrix invJ(3, 3);
    J.Invert(invJ);

    // find local coordinates of constrained node
    double lx, ly, lz;
    tet::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), m_nodes[0]->getCrds()(2), lx, ly, lz);

    // compute shape functions at constrained node
    static Vector N(4);
    for (int i = 0; i < 4; ++i)
        N(i) = tet::shapeFun(lx, ly, lz, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    static Matrix B(3, 15);
    B.Zero();
    for (int i = 0; i < 3; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 4; i++) {
        int j = 3 + i * 3;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j + 2) = N(i);
    }

    // Penalty stiffness
    double iK = m_K * V;

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}

const Matrix& ASDEmbeddedNodeElement::TET_3D_UR()
{
    // output
    static Matrix K(18, 18);

    // collect tetrahedron coordinates
    static Matrix X(3, 4);
    for (int i = 1; i < 5; i++) {
        const Node* node = m_nodes[static_cast<std::size_t>(i)];
        X(0, i - 1) = node->getCrds()(0);
        X(1, i - 1) = node->getCrds()(1);
        X(2, i - 1) = node->getCrds()(2);
    }

    // shape functions natural derivatives
    static Matrix dN(4, 3);
    tet::shapeFunDer(dN);

    // jacobian
    static Matrix J(3, 3);
    J.addMatrixProduct(0.0, X, dN, 1.0);
    double detJ = det3(J);
    double V = detJ / 6.0;
    static Matrix invJ(3, 3);
    J.Invert(invJ);

    // shape functions cartesian derivatives
    static Matrix dNdX(4, 3);
    dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

    // find local coordinates of constrained node
    double lx, ly, lz;
    tet::localCoord(X, invJ, m_nodes[0]->getCrds()(0), m_nodes[0]->getCrds()(1), m_nodes[0]->getCrds()(2), lx, ly, lz);

    // compute shape functions at constrained node
    static Vector N(4);
    for (int i = 0; i < 4; ++i)
        N(i) = tet::shapeFun(lx, ly, lz, i);

    // compute B matrix
    // UCx = sum(N*URx) -> sum(N*URx) - UCx = 0
    // UCy = sum(N*URy) -> sum(N*URy) - UCy = 0
    // UCz = sum(N*URz) -> sum(N*URz) - UCz = 0
    // RCx = sum(d_URz_dY - d_URy_dZ)/2.0 -> sum(d_URz_dY - d_URy_dZ)/2.0 - RCx = 0
    // RCy = sum(d_URx_dZ - d_URz_dX)/2.0 -> sum(d_URx_dZ - d_URz_dX)/2.0 - RCy = 0
    // RCz = sum(d_URy_dX - d_URx_dY)/2.0 -> sum(d_URy_dX - d_URx_dY)/2.0 - RCz = 0
    static Matrix B(6, 18);
    B.Zero();
    for (int i = 0; i < 6; i++)
        B(i, i) = -1.0;
    for (int i = 0; i < 4; i++) {
        int j = 6 + i * 3;
        B(0, j) = N(i);
        B(1, j + 1) = N(i);
        B(2, j + 2) = N(i);
        B(3, j + 1) = -dNdX(i, 2) / 2.0; B(3, j + 2) = dNdX(i, 1) / 2.0;
        B(4, j) = dNdX(i, 2) / 2.0; B(4, j + 2) = -dNdX(i, 0) / 2.0;
        B(5, j) = -dNdX(i, 1) / 2.0; B(5, j + 1) = dNdX(i, 0) / 2.0;
    }

    // Penalty stiffness
    double iK = m_K * V;

    // compute stiffness
    K.addMatrixTransposeProduct(0.0, B, B, iK);

    // done
    return K;
}
