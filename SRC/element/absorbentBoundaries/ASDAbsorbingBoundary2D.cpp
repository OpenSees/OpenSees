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
#include <TimeSeries.h>

#include <limits>
#include <cmath>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>


namespace {

    // flags for boundary types
    static constexpr int BND_NONE = 0;
    static constexpr int BND_BOTTOM = (1 << 1);
    static constexpr int BND_LEFT = (1 << 2);
    static constexpr int BND_RIGHT = (1 << 3);

    // a simple wrapper for sorting nodes
    struct SortedNode {
        // local position of this node (0 to 3)
        std::size_t pos = 0;
        // cartesian coordinates
        double x = 0.0;
        double y = 0.0;
        // number of dofs at this node
        int ndf = 0;
        // tolerance against round-off errors in formatting of floating-point numbers in input files...
        double tolerance = 1.0e-12;
        SortedNode() = default;
        SortedNode(std::size_t _pos, Node* n)
            : pos(_pos)
        {
            const Vector& v = n->getCrds();
            x = v(0);
            y = v(1);
            ndf = n->getNumberDOF();
        }
    };

    // compute tolerance
    inline void computeTolerance(std::vector<SortedNode>& n) {
        double xmin = std::numeric_limits<double>::max();
        double xmax = -xmin;
        double ymin = xmin;
        double ymax = xmax;
        for (const SortedNode& ni : n) {
            xmin = std::min(xmin, ni.x);
            xmax = std::max(xmax, ni.x);
            ymin = std::min(ymin, ni.y);
            ymax = std::max(ymax, ni.y);
        }
        double dx = std::abs(xmax - xmin);
        double dy = std::abs(ymax - ymin);
        double dmax = std::max(dx, dy);
        double tol = std::max(1.0e-10 * dmax, std::numeric_limits<double>::epsilon());
        for (SortedNode& ni : n)
            ni.tolerance = tol;
    }

    // sorts a node left-to-right (first) and bottom-to-top (second)
    struct SorterLeft {
        bool operator()(const SortedNode& a, const SortedNode& b) const {
            if (a.x < b.x - a.tolerance) return true;
            if (a.x > b.x + a.tolerance) return false;
            return a.y < b.y - a.tolerance;
        }
    };

    // sorts a node right-to-left (first) and bottom-to-top (second)
    struct SorterRight {
        bool operator()(const SortedNode& a, const SortedNode& b) const {
            if (a.x > b.x + a.tolerance) return true;
            if (a.x < b.x - a.tolerance) return false;
            return a.y < b.y - a.tolerance;
        }
    };

    // sort nodes with a custom sorter rule
    template<class T>
    inline void sortNodes(const std::vector<SortedNode>& nodes, std::vector<std::size_t>& ids, ID& dofs, int& ndf) {
        std::vector<SortedNode> aux = nodes;
        computeTolerance(aux);
        std::sort(aux.begin(), aux.end(), T());
        ids.resize(aux.size());
        dofs.resize(static_cast<int>(aux.size()) * 2);
        ndf = 0;
        for (std::size_t i = 0; i < aux.size(); ++i) {
            ids[i] = aux[i].pos;
            int j = static_cast<int>(i) * 2;
            dofs[j] = ndf;
            dofs[j + 1] = ndf + 1;
            ndf += nodes[i].ndf;
        }
        ID aux_dofs = dofs;
        for (std::size_t i = 0; i < aux.size(); ++i) {
            int j = static_cast<int>(i) * 2;
            int q = static_cast<int>(ids[i]) * 2;
            dofs[j] = aux_dofs[q];
            dofs[j + 1] = aux_dofs[q + 1];
        }
    }

}

void*
OPS_ASDAbsorbingBoundary2D(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDAbsorbingBoundary2D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    const char* descr = "Want: element ASDAbsorbingBoundary2D $tag $n1 $n2 $n3 $n4 $G $v $rho $thickness $btype <-fx $tsxTag> <-fy $tsyTag>\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 10) {
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

    // string parameter
    const char* btype = OPS_GetString();
    int bflag = BND_NONE;
    if (strstr(btype, "B")) bflag |= BND_BOTTOM;
    if (strstr(btype, "L")) bflag |= BND_LEFT;
    if (strstr(btype, "R")) bflag |= BND_RIGHT;
    if (bflag == BND_NONE) {
        opserr << "ASDAbsorbingBoundary2D ERROR: Invalid string mandatory value: the $btype "
            "argument should contain at least one of the following characters:\n"
            "'B', 'L', 'R'.\n" << descr;
        return 0;
    }

    // optional time series
    TimeSeries* fx = nullptr;
    TimeSeries* fy = nullptr;
    // only on bottom boundaries
    if (bflag & BND_BOTTOM) {
        numData = 1;
        int tsTag = 0;
        // util: get x
        auto get_fx = [&numData, &tsTag, &fx, descr]() -> bool {
            if (OPS_GetInt(&numData, &tsTag) != 0) {
                opserr << "ASDAbsorbingBoundary2D ERROR: Invalid integer for -fx optional time series.\n" << descr;
                return false;
            }
            fx = OPS_getTimeSeries(tsTag);
            if (fx == nullptr) {
                opserr << "ASDAbsorbingBoundary2D ERROR: Cannot find -fx time series with id = " << tsTag << ".\n" << descr;
                return false;
            }
            return true;
        };
        // util: get y
        auto get_fy = [&numData, &tsTag, &fy, descr]() -> bool {
            if (OPS_GetInt(&numData, &tsTag) != 0) {
                opserr << "ASDAbsorbingBoundary2D ERROR: Invalid integer for -fy optional time series.\n" << descr;
                return false;
            }
            fy = OPS_getTimeSeries(tsTag);
            if (fy == nullptr) {
                opserr << "ASDAbsorbingBoundary2D ERROR: Cannot find -fy time series with id = " << tsTag << ".\n" << descr;
                return false;
            }
            return true;
        };
        // parse first
        numArgs = OPS_GetNumRemainingInputArgs();
        if (numArgs > 1) {
            const char* key = OPS_GetString();
            if (strcmp(key, "-fx") == 0) {
                if (!get_fx()) return 0;
            }
            else if (strcmp(key, "-fy") == 0) {
                if (!get_fy()) return 0;
            }
            else {
                opserr << "ASDAbsorbingBoundary2D ERROR: Invalid optional flag \"" << key << "\".\n" << descr;
                return 0;
            }
        }
        // parse second
        numArgs = OPS_GetNumRemainingInputArgs();
        if (numArgs > 1) {
            const char* key = OPS_GetString();
            if (strcmp(key, "-fx") == 0) {
                if (fx) {
                    opserr << "ASDAbsorbingBoundary2D ERROR: -fx flag specified twice!.\n" << descr;
                    return 0;
                }
                if (!get_fx()) return 0;
            }
            else if (strcmp(key, "-fy") == 0) {
                if (fy) {
                    opserr << "ASDAbsorbingBoundary2D ERROR: -fy flag specified twice!.\n" << descr;
                    return 0;
                }
                if (!get_fy()) return 0;
            }
            else {
                opserr << "ASDAbsorbingBoundary2D ERROR: Invalid optional flag \"" << key << "\".\n" << descr;
                return 0;
            }
        }
    }

    // done
    return new ASDAbsorbingBoundary2D(iData[0], iData[1], iData[2], iData[3], iData[4], dData[0], dData[1], dData[2], dData[3], bflag, fx, fy);
}

ASDAbsorbingBoundary2D::ASDAbsorbingBoundary2D()
    : Element(0, ELE_TAG_ASDAbsorbingBoundary2D)
{
}

ASDAbsorbingBoundary2D::ASDAbsorbingBoundary2D(
    int tag,
    int node1,
    int node2,
    int node3,
    int node4,
    double G,
    double v,
    double rho,
    double thickness,
    int btype,
    TimeSeries* actionx,
    TimeSeries* actiony)
    : Element(tag, ELE_TAG_ASDAbsorbingBoundary2D)
    , m_G(G)
    , m_v(v)
    , m_rho(rho)
    , m_thickness(thickness)
    , m_boundary(btype)
{
    // save node ids
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
    m_node_ids(2) = node3;
    m_node_ids(3) = node4;

    // copy time-series
    if (actionx)
        m_tsx = actionx->getCopy();
    if (actiony)
        m_tsy = actiony->getCopy();
}

ASDAbsorbingBoundary2D::~ASDAbsorbingBoundary2D()
{
    // delete time-series
    if (m_tsx)
        delete m_tsx;
    if (m_tsy)
        delete m_tsy;
}

const char* ASDAbsorbingBoundary2D::getClassType(void) const
{
    return "ASDAbsorbingBoundary2D";
}

void ASDAbsorbingBoundary2D::setDomain(Domain* theDomain)
{
    // Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
        for (std::size_t i = 0; i < m_nodes.size(); ++i)
            m_nodes[i] = nullptr;
        return;
    }

    // check nodes and sort them in the in-coming order.
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
        if (ndf < 2) {
            opserr << "ASDAbsorbingBoundary2D Error in setDomain: In 2D at least 2 DOFs are required, not " << ndf << "\n";
            exit(-1);
        }

    }

    // only on setDomain during creation, not after a recvSelf
    if (!m_initialized) {

        // reorder nodes and compute node and dof mapping
        std::vector<SortedNode> sortednodes(m_nodes.size());
        for (std::size_t i = 0; i < m_nodes.size(); ++i)
            sortednodes[i] = SortedNode(i, m_nodes[i]);
        // all boundaries are sorted as LEFT, unless they are on the RIGHT side
        if (m_boundary & BND_RIGHT) {
            sortNodes<SorterRight>(sortednodes, m_node_map, m_dof_map, m_num_dofs);
        }
        else {
            sortNodes<SorterLeft>(sortednodes, m_node_map, m_dof_map, m_num_dofs);
        }

        // allocate initial displacement vector and also the static constraints
        m_U0.resize(m_num_dofs);
        m_U0.Zero();
        addDisplacement(m_U0);
        m_R0.resize(m_num_dofs);
        m_R0.Zero();

        // done with initialization
        m_initialized = true;

    }

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
    static Vector R;
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
        addBaseActions(R);
    }

    // done
    return R;
}

const Vector& ASDAbsorbingBoundary2D::getResistingForceIncInertia()
{
    // initialize vector
    static Vector R;
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
        addBaseActions(R);
        // add damping term
        addRCff(R);
        addRlk(R);
        // add mass term
        addRMff(R);
    }

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
    int res = 0;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // aux
    int pos = 0;

    // INT data
    static ID idData(28);

    // tag
    idData(pos++) = getTag();
    // nodes
    for (int i = 0; i < 4; ++i)
        idData(pos++) = m_node_ids(i);
    // stage
    idData(pos++) = static_cast<int>(m_stage);
    // boundary
    idData(pos++) = m_boundary;
    // num dofs
    idData(pos++) = m_num_dofs;
    // dof map
    for (int i = 0; i < 8; ++i)
        idData(pos++) = m_dof_map(i);
    // node map
    for (std::size_t i = 0; i < 4; ++i)
        idData(pos++) = static_cast<int>(m_node_map[i]);
    // time series data
    if (m_tsx) {
        idData(pos++) = 1;
        int dbtag = m_tsx->getDbTag();
        int classtag = m_tsx->getClassTag();
        if (dbtag == 0) {
            dbtag = theChannel.getDbTag();
            m_tsx->setDbTag(dbtag);
        }
        idData(pos++) = classtag;
        idData(pos++) = dbtag;
    }
    else {
        idData(pos++) = 0;
        idData(pos++) = 0;
        idData(pos++) = 0;
    }
    if (m_tsy) {
        idData(pos++) = 1;
        int dbtag = m_tsy->getDbTag();
        int classtag = m_tsy->getClassTag();
        if (dbtag == 0) {
            dbtag = theChannel.getDbTag();
            m_tsy->setDbTag(dbtag);
        }
        idData(pos++) = classtag;
        idData(pos++) = dbtag;
    }
    else {
        idData(pos++) = 0;
        idData(pos++) = 0;
        idData(pos++) = 0;
    }
    // initialization flag
    idData(pos++) = static_cast<int>(m_initialized);
    // double data size (not known at compile-time)
    int dsize = 4 + 2 * m_num_dofs;
    idData(pos++) = dsize;

    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDAbsorbingBoundary2D::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    static Vector vectData;
    vectData.resize(dsize);

    pos = 0;
    vectData(pos++) = m_G;
    vectData(pos++) = m_v;
    vectData(pos++) = m_rho;
    vectData(pos++) = m_thickness;
    for (int i = 0; i < m_num_dofs; ++i)
        vectData(pos++) = m_U0(i);
    for (int i = 0; i < m_num_dofs; ++i)
        vectData(pos++) = m_R0(i);

    res += theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDAbsorbingBoundary2D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // send time-series
    if (m_tsx) {
        if (m_tsx->sendSelf(commitTag, theChannel) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary2D::sendSelf() - " << this->getTag() << " failed to send TimeSeries (X)\n";
            return -1;
        }
    }
    if (m_tsy) {
        if (m_tsy->sendSelf(commitTag, theChannel) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary2D::sendSelf() - " << this->getTag() << " failed to send TimeSeries (Y)\n";
            return -1;
        }
    }

    // done
    return res;
}

int ASDAbsorbingBoundary2D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // aux
    int pos = 0;

    // INT data
    static ID idData(28);
    res += theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDAbsorbingBoundary2D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    // tag
    setTag(idData(pos++));
    // nodes
    for (int i = 0; i < 4; ++i)
        m_node_ids(i) = idData(pos++);
    // stage
    m_stage = static_cast<StageType>(idData(pos++));
    // boundary
    m_boundary = idData(pos++);
    // num dofs
    m_num_dofs = idData(pos++);
    // dof map
    for (int i = 0; i < 8; ++i)
        m_dof_map(i) = idData(pos++);
    // node map
    for (std::size_t i = 0; i < 4; ++i)
        m_node_map[i] = static_cast<std::size_t>(idData(pos++));
    // time series
    m_tsx = nullptr;
    m_tsy = nullptr;
    bool has_tsx = idData(pos++) == 1;
    int tsx_classtag = 0;
    int tsx_dbtag = 0;
    if (has_tsx) {
        tsx_classtag = idData(pos++);
        tsx_dbtag = idData(pos++);
    }
    else {
        pos += 2;
    }
    bool has_tsy = idData(pos++) == 1;
    int tsy_classtag = 0;
    int tsy_dbtag = 0;
    if (has_tsy) {
        tsy_classtag = idData(pos++);
        tsy_dbtag = idData(pos++);
    }
    else {
        pos += 2;
    }
    // initialization flag
    m_initialized = static_cast<bool>(idData(pos++));
    // double data size (not known at compile-time)
    int dsize = idData(pos++);

    // DOUBLE data
    static Vector vectData;
    vectData.resize(dsize);

    res += theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDAbsorbingBoundary2D::sendSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }

    pos = 0;
    m_G = vectData(pos++);
    m_v = vectData(pos++);
    m_rho = vectData(pos++);
    m_thickness = vectData(pos++);
    m_U0.resize(m_num_dofs);
    m_R0.resize(m_num_dofs);
    for (int i = 0; i < m_num_dofs; ++i)
        m_U0(i) = vectData(pos++);
    for (int i = 0; i < m_num_dofs; ++i)
        m_R0(i) = vectData(pos++);

    // receive time-series
    if (has_tsx) {
        m_tsx = theBroker.getNewTimeSeries(tsx_classtag);
        if (m_tsx == nullptr) {
            opserr << "WARNING ASDAbsorbingBoundary2D::recvSelf() - " << this->getTag() << " failed to create TimeSeries (X)\n";
            return -1;
        }
        m_tsx->setDbTag(tsx_dbtag);
        if (m_tsx->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary2D::recvSelf() - " << this->getTag() << " failed to recv TimeSeries (X)\n";
            return -1;
        }
    }
    if (has_tsy) {
        m_tsy = theBroker.getNewTimeSeries(tsy_classtag);
        if (m_tsy == nullptr) {
            opserr << "WARNING ASDAbsorbingBoundary2D::recvSelf() - " << this->getTag() << " failed to create TimeSeries (Y)\n";
            return -1;
        }
        m_tsy->setDbTag(tsy_dbtag);
        if (m_tsy->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary2D::recvSelf() - " << this->getTag() << " failed to recv TimeSeries (Y)\n";
            return -1;
        }
    }

    // done
    return res;
}

int ASDAbsorbingBoundary2D::setParameter(const char** argv, int argc, Parameter& param)
{
    // initial check
    if (argc < 1)
        return -1;

    // update stage?
    if (strcmp(argv[0], "stage") == 0) {
        return param.addObject(1, this);
    }
    // update G?
    else if (strcmp(argv[0], "G") == 0) {
        return param.addObject(2, this);
    }
    // update v?
    else if (strcmp(argv[0], "v") == 0) {
        return param.addObject(3, this);
    }
    // update rho?
    else if (strcmp(argv[0], "rho") == 0) {
        return param.addObject(4, this);
    }

    // no valid parameter
    return -1;
}

int ASDAbsorbingBoundary2D::updateParameter(int parameterID, Information& info)
{
    switch (parameterID)
    {
    case 1:
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
                opserr << "Error in ASDAbsorbingBoundary2D::updateParameter (element = " << getTag() << ").\n"
                    "Current stage = 0 (Stage_StaticConstraint).\n"
                    "The next stage can only be 1 (Stage_Absorbing), not " << istage << "!\n";
                exit(-1);
            }
        }
        else {
            // We cannot move from the Stage_Absorbing
            opserr << "Error in ASDAbsorbingBoundary2D::updateParameter (element = " << getTag() << ").\n"
                "Current stage = " << static_cast<int>(m_stage) << " (Stage_Absorbing).\n"
                "You cannot change the stage at this point!\n";
            exit(-1);
        }
    case 2:
        // G
        m_G = info.theDouble;
        return 0;
    case 3:
        // v
        m_v = info.theDouble;
        return 0;
    case 4:
        // rho
        m_rho = info.theDouble;
        return 0;
    default:
        return -1;
    }
}

Response* ASDAbsorbingBoundary2D::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    // check and quick return
    if (argc < 1)
        return nullptr;
    // support responses such as "E" or "material 1 E"
    int iarg = 0;
    if (argc == 3) {
        if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "integrPoint") == 0) {
            int pointNum = atoi(argv[1]);
            if (pointNum == 1) {
                iarg = 2;
            }
        }
    }
    // check argument
    int rtype = 0;
    if (strcmp(argv[iarg], "stage") == 0)
        rtype = 1;
    else if (strcmp(argv[iarg], "G") == 0)
        rtype = 2;
    else if (strcmp(argv[iarg], "v") == 0)
        rtype = 3;
    else if (strcmp(argv[iarg], "rho") == 0)
        rtype = 4;
    else if (strcmp(argv[iarg], "E") == 0)
        rtype = 5;
    if (rtype == 0)
        return Element::setResponse(argv, argc, output);

    // prepare response meta-data
    output.tag("ElementOutput");
    output.attr("eleType", this->getClassType());
    output.attr("eleTag", this->getTag());
    int numNodes = this->getNumExternalNodes();
    const ID& nodes = this->getExternalNodes();
    static char nodeData[32];
    for (int i = 0; i < numNodes; i++) {
        sprintf(nodeData, "node%d", i + 1);
        output.attr(nodeData, nodes(i));
    }
    output.tag("GaussPoint");
    output.attr("number", 1);
    output.attr("eta", 0.0);
    output.attr("neta", 0.0);
    output.attr("zeta", 0.0);
    output.tag("NdMaterialOutput");
    switch (rtype)
    {
    case 1: output.tag("ResponseType", "stage"); break;
    case 2: output.tag("ResponseType", "G"); break;
    case 3: output.tag("ResponseType", "v"); break;
    case 4: output.tag("ResponseType", "rho"); break;
    case 5: output.tag("ResponseType", "E"); break;
    default:
        break;
    }
    output.endTag(); // GaussPoint
    output.endTag(); // NdMaterialOutput
    output.endTag(); // ElementOutput

    // done
    return new ElementResponse(this, rtype, Vector(1));
}

int ASDAbsorbingBoundary2D::getResponse(int responseID, Information& eleInfo)
{
    static Vector r1(1);
    switch (responseID) {
    case 1: r1(0) = static_cast<double>(m_stage); return eleInfo.setVector(r1);
    case 2: r1(0) = m_G; return eleInfo.setVector(r1);
    case 3: r1(0) = m_v; return eleInfo.setVector(r1);
    case 4: r1(0) = m_rho; return eleInfo.setVector(r1);
    case 5: r1(0) = 2.0 * m_G * (1.0 + m_v); return eleInfo.setVector(r1);
    default:
        return Element::getResponse(responseID, eleInfo);
    }
}

void ASDAbsorbingBoundary2D::addDisplacement(Vector& U)
{
    int counter = 0;
    for (Node* node : m_nodes) {
        const Vector& iU = node->getTrialDisp();
        for (int i = 0; i < iU.Size(); ++i)
            U(counter++) += iU(i);
    }
}

const Vector& ASDAbsorbingBoundary2D::getDisplacement()
{
    static Vector U;
    U.resize(m_num_dofs);
    U.Zero();
    addDisplacement(U);
    U.addVector(1.0, m_U0, -1.0);
    return U;
}

const Vector& ASDAbsorbingBoundary2D::getVelocity()
{
    static Vector U;
    U.resize(m_num_dofs);
    int counter = 0;
    for (Node* node : m_nodes) {
        const Vector& iU = node->getTrialVel();
        for (int i = 0; i < iU.Size(); ++i)
            U(counter++) = iU(i);
    }
    return U;
}

const Vector& ASDAbsorbingBoundary2D::getAcceleration()
{
    static Vector U;
    U.resize(m_num_dofs);
    int counter = 0;
    for (Node* node : m_nodes) {
        const Vector& iU = node->getTrialAccel();
        for (int i = 0; i < iU.Size(); ++i)
            U(counter++) = iU(i);
    }
    return U;
}

void ASDAbsorbingBoundary2D::getElementSizes(double& lx, double& ly, double& nx)
{
    // get nodes
    Node* N0 = m_nodes[m_node_map[0]];
    Node* N1 = m_nodes[m_node_map[1]];
    Node* N2 = m_nodes[m_node_map[2]];
    // this should be always positive due to sorting...
    ly = std::abs(N1->getCrds()(1) - N0->getCrds()(1));
    // this can be negative on right element
    lx = std::abs(N2->getCrds()(0) - N0->getCrds()(0));
    // nx positive if the outward normal vector from soil domain
    // points in the positive +X direction.
    // note: change the sign here: external forces on soil domain
    if (m_boundary & BND_RIGHT)
        nx = -1.0;
    else
        nx = 1.0;
}

void ASDAbsorbingBoundary2D::updateStage()
{
    // save reactions
    m_R0.Zero();
    addRPenaltyStage0(m_R0);
    m_R0 *= -1.0;

    // accumulate current displacements
    addDisplacement(m_U0);

    // update stage
    m_stage = Stage_Absorbing;
}

void ASDAbsorbingBoundary2D::penaltyFactor(double& sp, double& mp)
{
    // order of magnitude of the (approximate) max stiffness value
    // of an equivalent 2D solid element with this G
    int oom = static_cast<int>(std::round(std::log10(m_G * m_thickness)));
    // ideal power = oom + 8 (8 being half the number of significant digits of a double floating number)
    int psp = oom + 8;
    int pmp = oom + 3;
    // compute the penalty factor
    // compute the penalty factor
    sp = std::pow(10.0, psp);
    mp = std::pow(10.0, pmp);
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
    double sp, mp;
    penaltyFactor(sp, mp);

    if (m_boundary & BND_BOTTOM) {
        // fix Uy at all nodes
        for (int i = 0; i < 4; ++i) {
            int j1 = i * 2 + 1;
            int q1 = m_dof_map(j1);
            K(q1, q1) += sp;
        }
        // edof Ux at nodes 0-1 and 2-3
        for (int i = 0; i < 2; ++i) {
            int n1 = i * 2;
            int n2 = n1 + 1;
            int j1 = n1 * 2;
            int j2 = n2 * 2;
            int q1 = m_dof_map(j1);
            int q2 = m_dof_map(j2);
            K(q1, q1) += mp;
            K(q2, q2) += mp;
            K(q1, q2) -= mp;
            K(q2, q1) -= mp;
        }
    }
    else {
        // fix Ux at all nodes
        for (int i = 0; i < 4; ++i) {
            int j1 = i * 2;
            int q1 = m_dof_map(j1);
            K(q1, q1) += sp;
        }
        // edof Uy at nodes 0-2 and 1-3
        for (int i = 0; i < 2; ++i) {
            int n1 = i;
            int n2 = i + 2;
            int j1 = n1 * 2 + 1;
            int j2 = n2 * 2 + 1;
            int q1 = m_dof_map(j1);
            int q2 = m_dof_map(j2);
            K(q1, q1) += mp;
            K(q2, q2) += mp;
            K(q1, q2) -= mp;
            K(q2, q1) -= mp;
        }
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
    double sp, mp;
    penaltyFactor(sp, mp);

    // get displacement vector
    const Vector& U = getDisplacement();

    if (m_boundary & BND_BOTTOM) {
        // fix Uy at all nodes
        for (int i = 0; i < 4; ++i) {
            int j1 = i * 2 + 1;
            int q1 = m_dof_map(j1);
            R(q1) += sp * U(q1);
        }
        // edof Ux at nodes 0-1 and 2-3
        for (int i = 0; i < 2; ++i) {
            int n1 = i * 2;
            int n2 = n1 + 1;
            int j1 = n1 * 2;
            int j2 = n2 * 2;
            int q1 = m_dof_map(j1);
            int q2 = m_dof_map(j2);
            R(q1) += mp * (U(q1) - U(q2));
            R(q2) += mp * (U(q2) - U(q1));
        }
    }
    else {
        // fix Ux at all nodes
        for (int i = 0; i < 4; ++i) {
            int j1 = i * 2;
            int q1 = m_dof_map(j1);
            R(q1) += sp * U(q1);
        }
        // edof Uy at nodes 0-2 and 1-3
        for (int i = 0; i < 2; ++i) {
            int n1 = i;
            int n2 = i + 2;
            int j1 = n1 * 2 + 1;
            int j2 = n2 * 2 + 1;
            int q1 = m_dof_map(j1);
            int q2 = m_dof_map(j2);
            R(q1) += mp * (U(q1) - U(q2));
            R(q2) += mp * (U(q2) - U(q1));
        }
    }
}

void ASDAbsorbingBoundary2D::addKPenaltyStage1(Matrix& K)
{
    // Enforce constraints in stage = 1.
    // In this stage the edge opposed to the soil domain:
    // 1) is fixed in both X and Y only on Horizontal boundary [SP constraint].

    // skip vertical boundary
    if (!(m_boundary & BND_BOTTOM))
        return;

    // penalty factor
    double sp, mp;
    penaltyFactor(sp, mp);

    // SP Ux and Uy for all nodes not in the soil domain
    // (nodes 0 and 2 in bottom boundary)
    // nodes:  0     2
    //  dofs: 0 1   4 5
    for (int i = 0; i < 2; ++i) {
        int j = i * 2;
        int j1 = j * 2;
        int j2 = j1 + 1;
        int q1 = m_dof_map(j1);
        int q2 = m_dof_map(j2);
        K(q1, q1) += sp;
        K(q2, q2) += sp;
    }
}

void ASDAbsorbingBoundary2D::addRPenaltyStage1(Vector& R)
{
    // Enforce constraints in stage = 1.
    // In this stage the edge opposed to the soil domain:
    // 1) is fixed in both X and Y only on Horizontal boundary [SP constraint].

    // skip vertical boundary
    if (!(m_boundary & BND_BOTTOM))
        return;

    // skip if computing reactions
    if (m_is_computing_reactions)
        return;

    // penalty factor
    double sp, mp;
    penaltyFactor(sp, mp);

    // get displacement vector
    const Vector& U = getDisplacement();

    // SP Ux and Uy for all nodes not in the soil domain
    // (nodes 0 and 2 in bottom boundary)
    // nodes:  0     2
    //  dofs: 0 1   4 5
    for (int i = 0; i < 2; ++i) {
        int j = i * 2;
        int j1 = j * 2;
        int j2 = j1 + 1;
        int q1 = m_dof_map(j1);
        int q2 = m_dof_map(j2);
        R(q1) += sp * U(q1);
        R(q2) += sp * U(q2);
    }
}

void ASDAbsorbingBoundary2D::addRReactions(Vector& R)
{
    // In stage 1, nodes on the soil domain side were fixed in either X or Y direction.
    // In stage 2 those constraints are removed, so we need to restore
    // those reactions as external forces.

    // todo: skip if computing reactions??
    if (m_is_computing_reactions)
        return;

    // add restoring reaction forces as external loads (-)
    R.addVector(1.0, m_R0, -1.0);
}

void ASDAbsorbingBoundary2D::addMff(Matrix& M, double scale)
{
    // Add the mass terms due to the Free-Field columns, thus only 
    // on vertical boundaries (nodes 0-1).
    // The optional scale parameter can be used while assembling
    // the damping due to the free-field

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // lumped mass (half mass of the whole element)
    double lx, ly, nx;
    getElementSizes(lx, ly, nx);
    double hm = scale * m_rho * m_thickness * lx * ly / 2.0;
    for (int i = 0; i < 2; ++i) {
        int j1 = i * 2;
        int j2 = j1 + 1;
        int q1 = m_dof_map(j1);
        int q2 = m_dof_map(j2);
        M(q1, q1) += hm;
        M(q2, q2) += hm;
    }
}

void ASDAbsorbingBoundary2D::addRMff(Vector& R)
{
    // Add the mass terms due to the Free-Field columns, thus only 
    // on vertical boundaries (nodes 0-1).

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // get acceleration
    const Vector& A = getAcceleration();

    // lumped mass (half mass of the whole element)
    double lx, ly, nx;
    getElementSizes(lx, ly, nx);
    double hm = m_rho * m_thickness * lx * ly / 2.0;
    for (int i = 0; i < 2; ++i) {
        int j1 = i * 2;
        int j2 = j1 + 1;
        int q1 = m_dof_map(j1);
        int q2 = m_dof_map(j2);
        R(q1) += hm * A(q1);
        R(q2) += hm * A(q2);
    }
}

void ASDAbsorbingBoundary2D::addKff(Matrix& K, double scale)
{
    // Add the stiffness matrix of the free-field column.
    // Only on vertical boundaries on nodes 0-1.
    // The optional scale parameter can be used while assembling
    // the damping due to the free-field

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // stiffness coefficients
    double lx, ly, nx;
    getElementSizes(lx, ly, nx);
    double t = m_thickness;
    double kx = scale * lx * mu * t / ly;
    double ky = scale * lx * t * (lam + 2.0 * mu) / ly;

    // fill
    K(m_dof_map(0), m_dof_map(0)) += kx;
    K(m_dof_map(0), m_dof_map(2)) += -kx;
    K(m_dof_map(1), m_dof_map(1)) += ky;
    K(m_dof_map(1), m_dof_map(3)) += -ky;
    K(m_dof_map(2), m_dof_map(0)) += -kx;
    K(m_dof_map(2), m_dof_map(2)) += kx;
    K(m_dof_map(3), m_dof_map(1)) += -ky;
    K(m_dof_map(3), m_dof_map(3)) += ky;
}

void ASDAbsorbingBoundary2D::addRff(Vector& R)
{
    // Add the stiffness matrix of the free-field column.
    // Only on vertical boundaries on nodes 0-1

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // stiffness coefficients
    double lx, ly, nx;
    getElementSizes(lx, ly, nx);
    double t = m_thickness;
    double kx = lx * mu * t / ly;
    double ky = lx * t * (lam + 2.0 * mu) / ly;

    // get displacement
    const Vector& U = getDisplacement();

    // fill
    R(m_dof_map(0)) += kx * (U(m_dof_map(0)) - U(m_dof_map(2)));
    R(m_dof_map(1)) += ky * (U(m_dof_map(1)) - U(m_dof_map(3)));
    R(m_dof_map(2)) += kx * (-U(m_dof_map(0)) + U(m_dof_map(2)));
    R(m_dof_map(3)) += ky * (-U(m_dof_map(1)) + U(m_dof_map(3)));
}

void ASDAbsorbingBoundary2D::addKffToSoil(Matrix& K)
{
    // Add the stiffness matrix of the forces transfered from the
    // free-field column to the soil domain.
    // Only on vertical boundaries

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // geometry
    double lx, ly, nx;
    getElementSizes(lx, ly, nx);
    double t = m_thickness;

    // fill
    K(m_dof_map(4), m_dof_map(1)) += -lam * nx * t / 2.0;
    K(m_dof_map(4), m_dof_map(3)) += lam * nx * t / 2.0;
    K(m_dof_map(5), m_dof_map(0)) += -mu * nx * t / 2.0;
    K(m_dof_map(5), m_dof_map(2)) += mu * nx * t / 2.0;
    K(m_dof_map(6), m_dof_map(1)) += -lam * nx * t / 2.0;
    K(m_dof_map(6), m_dof_map(3)) += lam * nx * t / 2.0;
    K(m_dof_map(7), m_dof_map(0)) += -mu * nx * t / 2.0;
    K(m_dof_map(7), m_dof_map(2)) += mu * nx * t / 2.0;
}

void ASDAbsorbingBoundary2D::addRffToSoil(Vector& R)
{
    // Add the forces transfered from the
    // free-field column to the soil domain.
    // Only on vertical boundaries

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // geometry
    double lx, ly, nx;
    getElementSizes(lx, ly, nx);
    double t = m_thickness;

    // get displacement 
    const Vector& U = getDisplacement();

    // fill
    R(m_dof_map(4)) += lam * nx * t * (-U(m_dof_map(1)) + U(m_dof_map(3))) / 2.0;
    R(m_dof_map(5)) += mu * nx * t * (-U(m_dof_map(0)) + U(m_dof_map(2))) / 2.0;
    R(m_dof_map(6)) += lam * nx * t * (-U(m_dof_map(1)) + U(m_dof_map(3))) / 2.0;
    R(m_dof_map(7)) += mu * nx * t * (-U(m_dof_map(0)) + U(m_dof_map(2))) / 2.0;
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
    // Add the rayleigh damping matrix of the free-field

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
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
    // Add the rayleigh damping forces of the free-field

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
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
    const Vector& V = getVelocity();

    // compute damping forces
    R.addMatrixVector(1.0, C, V, 1.0);
}

void ASDAbsorbingBoundary2D::getLKcoeff(double& ap, double& as)
{
    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // wave velocities
    double vp = std::sqrt((lam + 2.0 * mu) / m_rho);
    double vs = std::sqrt(mu / m_rho);

    // sizes
    double lx, ly, nx;
    getElementSizes(lx, ly, nx);
    double t = m_thickness;
    double h = ly;

    // swap coefficients and lumping size on horizontal
    // boundaries
    if (m_boundary & BND_BOTTOM) {
        // if the boundary is horizontal:
        // 1) the element size is lx
        // 2) the vp and vs should be swapped
        h = lx;
        double aux = vp;
        vp = vs;
        vs = aux;
    }

    // coefficients (put the minus sign here: external forces acting on soil domain)
    ap = -vp * h * m_rho * t / 2.0;
    as = -vs * h * m_rho * t / 2.0;
}

void ASDAbsorbingBoundary2D::addClk(Matrix& C)
{
    // dashpot coefficients
    double ap, as;
    getLKcoeff(ap, as);

    // compute derivatives of dashpot forces
    if (m_boundary & BND_BOTTOM) {
        // if it's on the corners, add only to nodes 1-0 with twice the ap and sp coefficients
        if (m_boundary != BND_BOTTOM) {
            // for BOTTOM-LEFT and BOTTOM-RIGHT (corner) boundaries
            C(m_dof_map(2), m_dof_map(0)) += 2.0 * ap;
            C(m_dof_map(2), m_dof_map(2)) += -2.0 * ap;
            C(m_dof_map(3), m_dof_map(1)) += 2.0 * as;
            C(m_dof_map(3), m_dof_map(3)) += -2.0 * as;
        }
        else {
            // for BOTTOM (horizontal) boundaries
            C(m_dof_map(2), m_dof_map(0)) += ap;
            C(m_dof_map(2), m_dof_map(2)) += -ap;
            C(m_dof_map(3), m_dof_map(1)) += as;
            C(m_dof_map(3), m_dof_map(3)) += -as;
            C(m_dof_map(6), m_dof_map(4)) += ap;
            C(m_dof_map(6), m_dof_map(6)) += -ap;
            C(m_dof_map(7), m_dof_map(5)) += as;
            C(m_dof_map(7), m_dof_map(7)) += -as;
        }
    }
    else {
        // for LEFT and RIGHT (vertical) boundaries
        C(m_dof_map(4), m_dof_map(0)) += ap;
        C(m_dof_map(4), m_dof_map(4)) += -ap;
        C(m_dof_map(5), m_dof_map(1)) += as;
        C(m_dof_map(5), m_dof_map(5)) += -as;
        C(m_dof_map(6), m_dof_map(2)) += ap;
        C(m_dof_map(6), m_dof_map(6)) += -ap;
        C(m_dof_map(7), m_dof_map(3)) += as;
        C(m_dof_map(7), m_dof_map(7)) += -as;
    }
}

void ASDAbsorbingBoundary2D::addRlk(Vector& R)
{
    // get velocity
    const Vector& V = getVelocity();

    // dashpot coefficients
    double ap, as;
    getLKcoeff(ap, as);

    // compute dashpot forces
    if (m_boundary & BND_BOTTOM) {
        // if it's on the corners, add only to nodes 1-0 with twice the ap and sp coefficients
        if (m_boundary != BND_BOTTOM) {
            // for BOTTOM-LEFT and BOTTOM-RIGHT (corner) boundaries
            R(m_dof_map(2)) += 2.0 * ap * (V(m_dof_map(0)) - V(m_dof_map(2)));
            R(m_dof_map(3)) += 2.0 * as * (V(m_dof_map(1)) - V(m_dof_map(3)));
        }
        else {
            // for BOTTOM (horizontal) boundaries
            R(m_dof_map(2)) += ap * (V(m_dof_map(0)) - V(m_dof_map(2)));
            R(m_dof_map(3)) += as * (V(m_dof_map(1)) - V(m_dof_map(3)));
            R(m_dof_map(6)) += ap * (V(m_dof_map(4)) - V(m_dof_map(6)));
            R(m_dof_map(7)) += as * (V(m_dof_map(5)) - V(m_dof_map(7)));
        }
    }
    else {
        // for LEFT and RIGHT (vertical) boundaries
        R(m_dof_map(4)) += ap * (V(m_dof_map(0)) - V(m_dof_map(4)));
        R(m_dof_map(5)) += as * (V(m_dof_map(1)) - V(m_dof_map(5)));
        R(m_dof_map(6)) += ap * (V(m_dof_map(2)) - V(m_dof_map(6)));
        R(m_dof_map(7)) += as * (V(m_dof_map(3)) - V(m_dof_map(7)));
    }
}

void ASDAbsorbingBoundary2D::addBaseActions(Vector& R)
{
    // skip vertical boundaries
    if (!(m_boundary & BND_BOTTOM))
        return;

    // re-use the getLKcoeff utilty
    double ap, as;
    getLKcoeff(ap, as);
    // swap back...
    double aux = ap;
    ap = as;
    as = aux;

    // add forces
    if (m_tsx) {
        Domain* domain = getDomain();
        if (domain == nullptr) {
            opserr << "ASDAbsorbingBoundary2D Error: cannot get domain!\n";
            exit(-1);
        }
        double vel = m_tsx->getFactor(domain->getCurrentTime());
        double fx = 2.0 * vel * as;
        if (m_boundary != BND_BOTTOM) {
            // on corners only on node 1 with twice the intensity
            R(m_dof_map(2)) += 2.0 * fx;
        }
        else {
            // on standard bottom boundary
            R(m_dof_map(2)) += fx;
            R(m_dof_map(6)) += fx;
        }
    }
    if (m_tsy) {
        Domain* domain = getDomain();
        if (domain == nullptr) {
            opserr << "ASDAbsorbingBoundary2D Error: cannot get domain!\n";
            exit(-1);
        }
        double vel = m_tsy->getFactor(domain->getCurrentTime());
        double fy = 2.0 * vel * ap;
        if (m_boundary != BND_BOTTOM) {
            // on corners only on node 1 with twice the intensity
            R(m_dof_map(3)) += 2.0 * fy;
        }
        else {
            // on standard bottom boundary
            R(m_dof_map(3)) += fy;
            R(m_dof_map(7)) += fy;
        }
    }
}





