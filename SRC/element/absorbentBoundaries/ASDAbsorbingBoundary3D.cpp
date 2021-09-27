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

#include "ASDAbsorbingBoundary3D.h"

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <algorithm>
#include <string>
#include <array>

namespace {

    // flags for boundary types
    static constexpr int BND_NONE = 0;
    static constexpr int BND_BOTTOM = (1 << 1);
    static constexpr int BND_LEFT = (1 << 2);
    static constexpr int BND_RIGHT = (1 << 3);
    static constexpr int BND_FRONT = (1 << 4);
    static constexpr int BND_BACK = (1 << 5);

    // a simple wrapper for sorting nodes
    struct SortedNode {
        // local position of this node (0 to 7)
        std::size_t pos = 0;
        // cartesian coordinates
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
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
            z = v(2);
            ndf = n->getNumberDOF();
        }
    };

    // compute tolerance
    inline void computeTolerance(std::vector<SortedNode>& n) {
        double xmin = std::numeric_limits<double>::max();
        double xmax = -xmin;
        double ymin = xmin;
        double ymax = xmax;
        double zmin = xmin;
        double zmax = xmax;
        for (const SortedNode& ni : n) {
            xmin = std::min(xmin, ni.x);
            xmax = std::max(xmax, ni.x);
            ymin = std::min(ymin, ni.y);
            ymax = std::max(ymax, ni.y);
            zmin = std::min(zmin, ni.z);
            zmax = std::max(zmax, ni.z);
        }
        double dx = std::abs(xmax - xmin);
        double dy = std::abs(ymax - ymin);
        double dz = std::abs(zmax - zmin);
        double dmax = std::max(dx, std::max(dy, dz));
        double tol = std::max(1.0e-10 * dmax, std::numeric_limits<double>::epsilon());
        for (SortedNode& ni : n)
            ni.tolerance = tol;
    }

    // sorts a node front-to-back (first) left-to-right (second) and bottom-to-top (third)
    struct SorterLeft {
        bool operator()(const SortedNode& a, const SortedNode& b) const {
            if (a.y < b.y - a.tolerance) return true;
            if (a.y > b.y + a.tolerance) return false;
            if (a.x < b.x - a.tolerance) return true;
            if (a.x > b.x + a.tolerance) return false;
            return a.z < b.z - a.tolerance;
        }
    };

    // sorts a node front-to-back (first) right-to-left (second) and bottom-to-top (third)
    struct SorterRight {
        bool operator()(const SortedNode& a, const SortedNode& b) const {
            if (a.y < b.y - a.tolerance) return true;
            if (a.y > b.y + a.tolerance) return false;
            if (a.x > b.x + a.tolerance) return true;
            if (a.x < b.x - a.tolerance) return false;
            return a.z < b.z - a.tolerance;
        }
    };

    // sorts a node left-to-right (first) front-to-back (second) and bottom-to-top (third)
    struct SorterFront {
        bool operator()(const SortedNode& a, const SortedNode& b) const {
            if (a.x < b.x - a.tolerance) return true;
            if (a.x > b.x + a.tolerance) return false;
            if (a.y < b.y - a.tolerance) return true;
            if (a.y > b.y + a.tolerance) return false;
            return a.z < b.z - a.tolerance;
        }
    };

    // sorts a node left-to-right (first) back-to-front (second) and bottom-to-top (third)
    struct SorterBack {
        bool operator()(const SortedNode& a, const SortedNode& b) const {
            if (a.x < b.x - a.tolerance) return true;
            if (a.x > b.x + a.tolerance) return false;
            if (a.y > b.y + a.tolerance) return true;
            if (a.y < b.y - a.tolerance) return false;
            return a.z < b.z - a.tolerance;
        }
    };

    // sorts a node right-to-left (first) back-to-front (second) and bottom-to-top (third)
    struct SorterBackRight {
        bool operator()(const SortedNode& a, const SortedNode& b) const {
            if (a.x > b.x + a.tolerance) return true;
            if (a.x < b.x - a.tolerance) return false;
            if (a.y > b.y + a.tolerance) return true;
            if (a.y < b.y - a.tolerance) return false;
            return a.z < b.z - a.tolerance;
        }
    };

    // sort nodes with a custom sorter rule
    template<class T>
    inline void sortNodes(const std::vector<SortedNode>& nodes, std::vector<std::size_t>& ids, ID& dofs, int& ndf) {
        std::vector<SortedNode> aux = nodes;
        computeTolerance(aux);
        std::sort(aux.begin(), aux.end(), T());
        ids.resize(aux.size());
        dofs.resize(static_cast<int>(aux.size()) * 3);
        ndf = 0;
        for (std::size_t i = 0; i < aux.size(); ++i) {
            ids[i] = aux[i].pos;
            int j = static_cast<int>(i) * 3;
            dofs[j] = ndf;
            dofs[j + 1] = ndf + 1;
            dofs[j + 2] = ndf + 2;
            ndf += nodes[i].ndf;
        }
        ID aux_dofs = dofs;
        for (std::size_t i = 0; i < aux.size(); ++i) {
            int j = static_cast<int>(i) * 3;
            int q = static_cast<int>(ids[i]) * 3;
            dofs[j] = aux_dofs[q];
            dofs[j + 1] = aux_dofs[q + 1];
            dofs[j + 2] = aux_dofs[q + 2];
        }
    }

    // node locations inside an element
    constexpr std::array<int, 4> N_BOTTOM = { {0,2,4,6} };
    constexpr std::array<int, 4> N_TOP = { {1,3,5,7} };
    constexpr std::array<int, 4> N_FF = { {0,1,4,5} };
    constexpr std::array<int, 4> N_SS = { {2,3,6,7} };
    // for vertical edges
    constexpr std::array<int, 2> NVE_IN = { {6,7} };
    constexpr std::array<std::array<int, 3>, 2> NVE_OUT = { {{0,2,4},{1,3,5}} };
    // for horizontal edges
    constexpr std::array<int, 2> NHE_IN = { {3,7} };
    constexpr std::array<std::array<int, 3>, 2> NHE_OUT = { {{0,1,2},{4,5,6}} };
    // for corners
    constexpr int NC_IN = 7;
    constexpr std::array<int, 7> NC_OUT = { {0,1,2,3,4,5,6} };

    // free-field ids for pure sides
    constexpr std::array<int, 4> FF_ID_SIDE = { {0,1,4,5} };
    constexpr std::array<int, 2> FF_ID_EDGE = { {0,1} };

    // dof labels
    constexpr int UX = 0;
    constexpr int UY = 1;
    constexpr int UZ = 2;

    // constraints utilities
    inline void cfix(int _n1, int _dof, Matrix& K, double p, const ID& m_dof_map) {
        int iL = _n1 * 3 + _dof;
        int iG = m_dof_map(iL);
        K(iG, iG) += p;
    };
    inline void cfix(int _n1, int _dof, Vector& R, const Vector& U, double p, const ID& m_dof_map) {
        int iL = _n1 * 3 + _dof;
        int iG = m_dof_map(iL);
        R(iG) += p * U(iG);
    };
    inline void cedof(int _n1, int _n2, int _dof, Matrix& K, double p, const ID& m_dof_map) {
        int iL = _n1 * 3 + _dof;
        int jL = _n2 * 3 + _dof;
        int iG = m_dof_map(iL);
        int jG = m_dof_map(jL);
        K(iG, iG) += p;
        K(jG, jG) += p;
        K(iG, jG) -= p;
        K(jG, iG) -= p;
    };
    inline void cedof(int _n1, int _n2, int _dof, Vector& R, const Vector& U, double p, const ID& m_dof_map, bool print=false) {
        int iL = _n1 * 3 + _dof;
        int jL = _n2 * 3 + _dof;
        int iG = m_dof_map(iL);
        int jG = m_dof_map(jL);
        R(iG) += p * (U(iG) - U(jG));
        R(jG) += p * (U(jG) - U(iG));
    };

    // dashpot utilities
    struct LKnodes {
        int ff = 0; // node on free-field (or fixed) boundary
        int ss = 0; // node on soild boundary
        double w = 1.0; // weight
        LKnodes() = default;
        LKnodes(int _ff, int _ss, double _w) : ff(_ff), ss(_ss), w(_w) {}
    };
    // for vertical boundaries
    // L, R, F, K
    static std::vector<LKnodes> LK_NODES_V_SIDE = { LKnodes(0,2,1.0), LKnodes(1,3,1.0), LKnodes(4,6,1.0), LKnodes(5,7,1.0) };
    // L-F, L-K, R-F, R-K
    static std::vector<LKnodes> LK_NODES_V_EDGE = { LKnodes(0,2,2.0), LKnodes(1,3,2.0), LKnodes(0,4,2.0), LKnodes(1,5,2.0) };
    // for horizontal boundaries
    // B
    static std::vector<LKnodes> LK_NODES_H_SIDE = { LKnodes(0,1,1.0), LKnodes(2,3,1.0), LKnodes(4,5,1.0), LKnodes(6,7,1.0) };
    // B-F, B-K, B-L, B-R
    static std::vector<LKnodes> LK_NODES_H_EDGE = { LKnodes(0,1,2.0), LKnodes(4,5,2.0) };
    // bottom corners
    static std::vector<LKnodes> LK_NODES_H_CORNER = { LKnodes(0,1,4.0) };
    // select appropriate pairs
    inline const std::vector<LKnodes>& LKselectPairs(int boundary) {
        if (boundary & BND_BOTTOM) {
            if ((boundary == (BND_BOTTOM | BND_LEFT | BND_FRONT)) ||
                (boundary == (BND_BOTTOM | BND_LEFT | BND_BACK)) ||
                (boundary == (BND_BOTTOM | BND_RIGHT | BND_FRONT)) ||
                (boundary == (BND_BOTTOM | BND_RIGHT | BND_BACK))) {
                return LK_NODES_H_CORNER;
            }
            else if (boundary == BND_BOTTOM) {
                return LK_NODES_H_SIDE;
            }
            else {
                return LK_NODES_H_EDGE;
            }
        }
        else {
            if ((boundary == BND_LEFT) ||
                (boundary == BND_RIGHT) ||
                (boundary == BND_FRONT) ||
                (boundary == BND_BACK)) {
                return LK_NODES_V_SIDE;
            }
            else {
                return LK_NODES_V_EDGE;
            }
        }
    }

    // utilities for 8-node hexa
    // H8 quadrature data
    constexpr double H8_G = 0.577350269189626;
    constexpr std::array<double, 8> H8_GX = { {-H8_G, H8_G, H8_G, -H8_G,   -H8_G, H8_G, H8_G, -H8_G} };
    constexpr std::array<double, 8> H8_GY = { {-H8_G, -H8_G, H8_G, H8_G,   -H8_G, -H8_G, H8_G, H8_G} };
    constexpr std::array<double, 8> H8_GZ = { {-H8_G, -H8_G, -H8_G, -H8_G,   H8_G, H8_G, H8_G, H8_G} };
    constexpr std::array<double, 8> H8_GW = { {1.0, 1.0, 1.0, 1.0,   1.0, 1.0, 1.0, 1.0} };
    // H8 node coordinate matrix
    inline void H8_nodeMatrix(const std::vector<Node*>& nodes, const std::vector<std::size_t>& node_map, Matrix& P) {
        for (int j = 0; j < 8; ++j) {
            const Node* n = nodes[static_cast<std::size_t>(j)];
            const Vector& x = n->getCrds();
            for (int i = 0; i < 3; ++i) {
                P(i, j) = x(i);
            }
        }
    }
    // H8 shape function local gradients at a gauss point
    inline void H8_dN(double x, double y, double z, Matrix& dN) {
        dN(0, 0) = -0.125 * (1.0 - y) * (1.0 - z);
        dN(0, 1) = -0.125 * (1.0 - x) * (1.0 - z);
        dN(0, 2) = -0.125 * (1.0 - x) * (1.0 - y);
        dN(1, 0) = 0.125 * (1.0 - y) * (1.0 - z);
        dN(1, 1) = -0.125 * (1.0 + x) * (1.0 - z);
        dN(1, 2) = -0.125 * (1.0 + x) * (1.0 - y);
        dN(2, 0) = 0.125 * (1.0 + y) * (1.0 - z);
        dN(2, 1) = 0.125 * (1.0 + x) * (1.0 - z);
        dN(2, 2) = -0.125 * (1.0 + x) * (1.0 + y);
        dN(3, 0) = -0.125 * (1.0 + y) * (1.0 - z);
        dN(3, 1) = 0.125 * (1.0 - x) * (1.0 - z);
        dN(3, 2) = -0.125 * (1.0 - x) * (1.0 + y);
        dN(4, 0) = -0.125 * (1.0 - y) * (1.0 + z);
        dN(4, 1) = -0.125 * (1.0 - x) * (1.0 + z);
        dN(4, 2) = 0.125 * (1.0 - x) * (1.0 - y);
        dN(5, 0) = 0.125 * (1.0 - y) * (1.0 + z);
        dN(5, 1) = -0.125 * (1.0 + x) * (1.0 + z);
        dN(5, 2) = 0.125 * (1.0 + x) * (1.0 - y);
        dN(6, 0) = 0.125 * (1.0 + y) * (1.0 + z);
        dN(6, 1) = 0.125 * (1.0 + x) * (1.0 + z);
        dN(6, 2) = 0.125 * (1.0 + x) * (1.0 + y);
        dN(7, 0) = -0.125 * (1.0 + y) * (1.0 + z);
        dN(7, 1) = 0.125 * (1.0 - x) * (1.0 + z);
        dN(7, 2) = 0.125 * (1.0 - x) * (1.0 + y);
    }
    // H8 elastic material matrix
    inline void H8_C0(double lam, double mu, Matrix& C) {
        C.Zero();
        C(0, 0) = C(1, 1) = C(2, 2) = lam + 2.0 * mu;
        C(0, 1) = C(1, 0) = lam;
        C(0, 2) = C(2, 0) = lam;
        C(1, 2) = C(2, 1) = lam;
        C(3, 3) = C(4, 4) = C(5, 5) = mu;
    }
    // H8 det 3x3
    inline double H8_det3(const Matrix& J) {
        return  J(0, 0) * J(1, 1) * J(2, 2) - J(0, 0) * J(1, 2) * J(2, 1) - J(0, 1) * J(1, 0) * J(2, 2) +
            J(0, 1) * J(1, 2) * J(2, 0) + J(0, 2) * J(1, 0) * J(2, 1) - J(0, 2) * J(1, 1) * J(2, 0);
    }
    // H8 B matrix
    inline void H8_Bmatrix(const Matrix& dN, Matrix& B) {
        B.Zero();
        for (int node = 0; node < 8; ++node) {
            int j = node * 3;
            B(0, j) = dN(node, 0);
            B(1, j + 1) = dN(node, 1);
            B(2, j + 2) = dN(node, 2);
            B(3, j) = dN(node, 1);
            B(3, j + 1) = dN(node, 0);
            B(4, j + 1) = dN(node, 2);
            B(4, j + 2) = dN(node, 1);
            B(5, j) = dN(node, 2);
            B(5, j + 2) = dN(node, 0);
        }
    }
    

}

void*
OPS_ASDAbsorbingBoundary3D(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDAbsorbingBoundary3D - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    const char* descr = "Want: element ASDAbsorbingBoundary3D $tag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 $G $v $rho $btype <-fx $tsxTag> <-fy $tsyTag> <-fz $tszTag>\n";

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 13) {
        opserr << "ASDAbsorbingBoundary3D ERROR : Few arguments:\n" << descr;
        return 0;
    }

    // int parameters
    int iData[9];
    int numData = 9;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "ASDAbsorbingBoundary3D ERROR: Invalid integer mandatory values: element ASDAbsorbingBoundary3D wants 9 integer parameters\n" << descr;
        return 0;
    }

    // double parameters
    double dData[3];
    numData = 3;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "ASDAbsorbingBoundary3D ERROR: Invalid double mandatory values: element ASDAbsorbingBoundary3D wants 3 double parameters\n" << descr;
        return 0;
    }

    // string parameter
    const char* btype = OPS_GetString();
    int bflag = BND_NONE;
    if (strstr(btype, "B")) bflag |= BND_BOTTOM;
    if (strstr(btype, "L")) bflag |= BND_LEFT;
    if (strstr(btype, "R")) bflag |= BND_RIGHT;
    if (strstr(btype, "F")) bflag |= BND_FRONT;
    if (strstr(btype, "K")) bflag |= BND_BACK;
    if(bflag == BND_NONE) {
        opserr << "ASDAbsorbingBoundary3D ERROR: Invalid string mandatory value: the $btype "
            "argument should contain at least one of the following characters:\n"
            "'B', 'L', 'R', 'F', 'K'.\n" << descr;
        return 0;
    }

    // optional time series
    TimeSeries* fx = nullptr;
    TimeSeries* fy = nullptr;
    TimeSeries* fz = nullptr;
    // only on bottom boundaries
    if (bflag & BND_BOTTOM) {
        numData = 1;
        int tsTag = 0;
        // util: get x
        auto get_fx = [&numData, &tsTag, &fx, descr]() -> bool {
            if (OPS_GetInt(&numData, &tsTag) != 0) {
                opserr << "ASDAbsorbingBoundary3D ERROR: Invalid integer for -fx optional time series.\n" << descr;
                return false;
            }
            fx = OPS_getTimeSeries(tsTag);
            if (fx == nullptr) {
                opserr << "ASDAbsorbingBoundary3D ERROR: Cannot find -fx time series with id = " << tsTag << ".\n" << descr;
                return false;
            }
            return true;
        };
        // util: get y
        auto get_fy = [&numData, &tsTag, &fy, descr]() -> bool {
            if (OPS_GetInt(&numData, &tsTag) != 0) {
                opserr << "ASDAbsorbingBoundary3D ERROR: Invalid integer for -fy optional time series.\n" << descr;
                return false;
            }
            fy = OPS_getTimeSeries(tsTag);
            if (fy == nullptr) {
                opserr << "ASDAbsorbingBoundary3D ERROR: Cannot find -fy time series with id = " << tsTag << ".\n" << descr;
                return false;
            }
            return true;
        };
        // util: get z
        auto get_fz = [&numData, &tsTag, &fz, descr]() -> bool {
            if (OPS_GetInt(&numData, &tsTag) != 0) {
                opserr << "ASDAbsorbingBoundary3D ERROR: Invalid integer for -fz optional time series.\n" << descr;
                return false;
            }
            fz = OPS_getTimeSeries(tsTag);
            if (fz == nullptr) {
                opserr << "ASDAbsorbingBoundary3D ERROR: Cannot find -fz time series with id = " << tsTag << ".\n" << descr;
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
            else if (strcmp(key, "-fz") == 0) {
                if (!get_fz()) return 0;
            }
            else {
                opserr << "ASDAbsorbingBoundary3D ERROR: Invalid optional flag \"" << key << "\".\n" << descr;
                return 0;
            }
        }
        // parse second and third
        for (int i = 0; i < 2; ++i) {
            numArgs = OPS_GetNumRemainingInputArgs();
            if (numArgs > 1) {
                const char* key = OPS_GetString();
                if (strcmp(key, "-fx") == 0) {
                    if (fx) {
                        opserr << "ASDAbsorbingBoundary3D ERROR: -fx flag specified twice!.\n" << descr;
                        return 0;
                    }
                    if (!get_fx()) return 0;
                }
                else if (strcmp(key, "-fy") == 0) {
                    if (fy) {
                        opserr << "ASDAbsorbingBoundary3D ERROR: -fy flag specified twice!.\n" << descr;
                        return 0;
                    }
                    if (!get_fy()) return 0;
                }
                else if (strcmp(key, "-fz") == 0) {
                    if (fz) {
                        opserr << "ASDAbsorbingBoundary3D ERROR: -fz flag specified twice!.\n" << descr;
                        return 0;
                    }
                    if (!get_fz()) return 0;
                }
                else {
                    opserr << "ASDAbsorbingBoundary3D ERROR: Invalid optional flag \"" << key << "\".\n" << descr;
                    return 0;
                }
            }
        }
    }

    // done
    return new ASDAbsorbingBoundary3D(
        iData[0], 
        iData[1], iData[2], iData[3], iData[4], iData[5], iData[6], iData[7], iData[8],
        dData[0], dData[1], dData[2], 
        bflag, fx, fy, fz);
}

ASDAbsorbingBoundary3D::ASDAbsorbingBoundary3D()
	: Element(0, ELE_TAG_ASDAbsorbingBoundary3D)
{
}

ASDAbsorbingBoundary3D::ASDAbsorbingBoundary3D(
    int tag,
    int node1,
    int node2,
    int node3,
    int node4,
    int node5,
    int node6,
    int node7,
    int node8,
    double G,
    double v,
    double rho,
    int btype,
    TimeSeries* actionx,
    TimeSeries* actiony,
    TimeSeries* actionz)
	: Element(tag, ELE_TAG_ASDAbsorbingBoundary3D)
    , m_G(G)
    , m_v(v)
    , m_rho(rho)
    , m_boundary(btype)
{
    // save node ids
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
    m_node_ids(2) = node3;
    m_node_ids(3) = node4;
    m_node_ids(4) = node5;
    m_node_ids(5) = node6;
    m_node_ids(6) = node7;
    m_node_ids(7) = node8;

    // copy time-series
    if (actionx)
        m_tsx = actionx->getCopy();
    if (actiony)
        m_tsy = actiony->getCopy();
    if (actionz)
        m_tsz = actionz->getCopy();
}

ASDAbsorbingBoundary3D::~ASDAbsorbingBoundary3D()
{
    // delete time-series
    if (m_tsx)
        delete m_tsx;
    if (m_tsy)
        delete m_tsy;
    if (m_tsz)
        delete m_tsz;
}

const char* ASDAbsorbingBoundary3D::getClassType(void) const
{
    return "ASDAbsorbingBoundary3D";
}

void ASDAbsorbingBoundary3D::setDomain(Domain* theDomain)
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
            opserr << "ASDAbsorbingBoundary3D Error in setDomain: node " << node_id << " does not exist in the domain\n";
            exit(-1);
        }

        // store node
        m_nodes[i] = node;

        // check NDM
        int ndm = node->getCrds().Size();
        if (ndm != 3) {
            opserr << "ASDAbsorbingBoundary3D Error in setDomain: Nodes must have 3 dimensions, not " << ndm << "\n";
            exit(-1);
        }

        // check NDF
        int ndf = node->getNumberDOF();
        if (ndf < 3) {
            opserr << "ASDAbsorbingBoundary3D Error in setDomain: In 3D at least 3 DOFs are required, not " << ndf << "\n";
            exit(-1);
        }

    }

    // only on setDomain during creation, not after a recvSelf
    if (!m_initialized) {

        // reorder nodes and compute node and dof mapping
        std::vector<SortedNode> sortednodes(m_nodes.size());
        for (std::size_t i = 0; i < m_nodes.size(); ++i)
            sortednodes[i] = SortedNode(i, m_nodes[i]);
        // choose appropriate sorting: 17 combinations
        if (m_boundary & BND_LEFT) {
            if (m_boundary & BND_BACK) {
                // LEFT-BACK, LEFT-BACK-BOTTOM
                sortNodes<SorterBack>(sortednodes, m_node_map, m_dof_map, m_num_dofs);
            }
            else {
                // LEFT, LEFT-FRONT, LEFT-BOTTOM, LEFT-FRONT-BOTTOM
                sortNodes<SorterLeft>(sortednodes, m_node_map, m_dof_map, m_num_dofs);
            }
        }
        else if (m_boundary & BND_RIGHT) {
            if (m_boundary & BND_BACK) {
                // RIGHT-BACK, RIGHT-BACK-BOTTOM
                sortNodes<SorterBackRight>(sortednodes, m_node_map, m_dof_map, m_num_dofs);
            }
            else {
                // RIGHT, RIGHT-FRONT, RIGHT-BOTTOM, RIGHT-FRONT-BOTTOM
                sortNodes<SorterRight>(sortednodes, m_node_map, m_dof_map, m_num_dofs);
            }
        }
        else if (m_boundary & BND_FRONT) {
            // FRONT FRONT-BOTTOM
            sortNodes<SorterFront>(sortednodes, m_node_map, m_dof_map, m_num_dofs);
        }
        else if (m_boundary & BND_BACK) {
            // BACK BACK-BOTTOM
            sortNodes<SorterBack>(sortednodes, m_node_map, m_dof_map, m_num_dofs);
        }
        else if(m_boundary & BND_BOTTOM) {
            // BOTTOM
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

void ASDAbsorbingBoundary3D::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ASDAbsorbingBoundary3D\t" << eleTag << " :";
        for (int i = 0; i < m_node_ids.Size(); ++i)
            s << "\t" << m_node_ids(i);
        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDAbsorbingBoundary3D\", ";
        s << "\"nodes\": [";
        for (int i = 0; i < m_node_ids.Size(); ++i) {
            if (i > 0)
                s << ", ";
            s << m_node_ids(i);
        }
        s << "]}";
    }
}

int ASDAbsorbingBoundary3D::getNumExternalNodes() const
{
    return m_node_ids.Size();
}

const ID& ASDAbsorbingBoundary3D::getExternalNodes()
{
    return m_node_ids;
}

Node** ASDAbsorbingBoundary3D::getNodePtrs()
{
    return m_nodes.data();
}

int ASDAbsorbingBoundary3D::getNumDOF()
{
    return m_num_dofs;
}

int ASDAbsorbingBoundary3D::revertToLastCommit()
{
    return 0;
}

const Matrix& ASDAbsorbingBoundary3D::getTangentStiff(void)
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

const Matrix& ASDAbsorbingBoundary3D::getInitialStiff(void)
{
    return getTangentStiff();
}

const Matrix& ASDAbsorbingBoundary3D::getDamp(void)
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

const Matrix& ASDAbsorbingBoundary3D::getMass(void)
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

int ASDAbsorbingBoundary3D::addInertiaLoadToUnbalance(const Vector& accel)
{
    // we don't need this!
    return 0;
}

const Vector& ASDAbsorbingBoundary3D::getResistingForce()
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

const Vector& ASDAbsorbingBoundary3D::getResistingForceIncInertia()
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

int ASDAbsorbingBoundary3D::addResistingForceToNodalReaction(int flag)
{
    m_is_computing_reactions = true;
    int result = Element::addResistingForceToNodalReaction(flag);
    m_is_computing_reactions = false;
    return result;
}

int ASDAbsorbingBoundary3D::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // aux
    int pos = 0;

    // INT data
    static ID idData(55);

    // tag
    idData(pos++) = getTag();
    // nodes
    for(int i = 0; i < 8; ++i)
        idData(pos++) = m_node_ids(i);
    // stage
    idData(pos++) = static_cast<int>(m_stage);
    // boundary
    idData(pos++) = m_boundary;
    // num dofs
    idData(pos++) = m_num_dofs;
    // dof map
    for (int i = 0; i < 24; ++i)
        idData(pos++) = m_dof_map(i);
    // node map
    for (std::size_t i = 0; i < 8; ++i)
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
    if (m_tsz) {
        idData(pos++) = 1;
        int dbtag = m_tsz->getDbTag();
        int classtag = m_tsz->getClassTag();
        if (dbtag == 0) {
            dbtag = theChannel.getDbTag();
            m_tsz->setDbTag(dbtag);
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
    int dsize = 3 + 2 * m_num_dofs;
    idData(pos++) = dsize;

    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDAbsorbingBoundary3D::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    static Vector vectData;
    vectData.resize(dsize);

    pos = 0;
    vectData(pos++) = m_G;
    vectData(pos++) = m_v;
    vectData(pos++) = m_rho;
    for (int i = 0; i < m_num_dofs; ++i)
        vectData(pos++) = m_U0(i);
    for (int i = 0; i < m_num_dofs; ++i)
        vectData(pos++) = m_R0(i);

    res += theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDAbsorbingBoundary3D::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // send time-series
    if (m_tsx) {
        if (m_tsx->sendSelf(commitTag, theChannel) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary3D::sendSelf() - " << this->getTag() << " failed to send TimeSeries (X)\n";
            return -1;
        }
    }
    if (m_tsy) {
        if (m_tsy->sendSelf(commitTag, theChannel) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary3D::sendSelf() - " << this->getTag() << " failed to send TimeSeries (Y)\n";
            return -1;
        }
    }
    if (m_tsz) {
        if (m_tsz->sendSelf(commitTag, theChannel) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary3D::sendSelf() - " << this->getTag() << " failed to send TimeSeries (Z)\n";
            return -1;
        }
    }

    // done
    return res;
}

int ASDAbsorbingBoundary3D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // aux
    int pos = 0;

    // INT data
    static ID idData(55);
    res += theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDAbsorbingBoundary3D::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    // tag
    setTag(idData(pos++));
    // nodes
    for (int i = 0; i < 8; ++i)
        m_node_ids(i) = idData(pos++);
    // stage
    m_stage = static_cast<StageType>(idData(pos++));
    // boundary
    m_boundary = idData(pos++);
    // num dofs
    m_num_dofs = idData(pos++);
    // dof map
    for (int i = 0; i < 24; ++i)
        m_dof_map(i) = idData(pos++);
    // node map
    for (std::size_t i = 0; i < 8; ++i)
        m_node_map[i] = static_cast<std::size_t>(idData(pos++));
    // time series
    m_tsx = nullptr;
    m_tsy = nullptr;
    m_tsz = nullptr;
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
    bool has_tsz = idData(pos++) == 1;
    int tsz_classtag = 0;
    int tsz_dbtag = 0;
    if (has_tsz) {
        tsz_classtag = idData(pos++);
        tsz_dbtag = idData(pos++);
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
        opserr << "WARNING ASDAbsorbingBoundary3D::sendSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }

    pos = 0;
    m_G = vectData(pos++);
    m_v = vectData(pos++);
    m_rho = vectData(pos++);
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
            opserr << "WARNING ASDAbsorbingBoundary3D::recvSelf() - " << this->getTag() << " failed to create TimeSeries (X)\n";
            return -1;
        }
        m_tsx->setDbTag(tsx_dbtag);
        if (m_tsx->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary3D::recvSelf() - " << this->getTag() << " failed to recv TimeSeries (X)\n";
            return -1;
        }
    }
    if (has_tsy) {
        m_tsy = theBroker.getNewTimeSeries(tsy_classtag);
        if (m_tsy == nullptr) {
            opserr << "WARNING ASDAbsorbingBoundary3D::recvSelf() - " << this->getTag() << " failed to create TimeSeries (Y)\n";
            return -1;
        }
        m_tsy->setDbTag(tsy_dbtag);
        if (m_tsy->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary3D::recvSelf() - " << this->getTag() << " failed to recv TimeSeries (Y)\n";
            return -1;
        }
    }
    if (has_tsz) {
        m_tsz = theBroker.getNewTimeSeries(tsy_classtag);
        if (m_tsz == nullptr) {
            opserr << "WARNING ASDAbsorbingBoundary3D::recvSelf() - " << this->getTag() << " failed to create TimeSeries (Z)\n";
            return -1;
        }
        m_tsz->setDbTag(tsy_dbtag);
        if (m_tsz->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "WARNING ASDAbsorbingBoundary3D::recvSelf() - " << this->getTag() << " failed to recv TimeSeries (Z)\n";
            return -1;
        }
    }

    // done
    return res;
}

int ASDAbsorbingBoundary3D::setParameter(const char** argv, int argc, Parameter& param)
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

int ASDAbsorbingBoundary3D::updateParameter(int parameterID, Information& info)
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
                opserr << "Error in ASDAbsorbingBoundary3D::updateParameter (element = " << getTag() << ").\n"
                    "Current stage = 0 (Stage_StaticConstraint).\n"
                    "The next stage can only be 1 (Stage_Absorbing), not " << istage << "!\n";
                exit(-1);
            }
        }
        else {
            // We cannot move from the Stage_Absorbing
            opserr << "Error in ASDAbsorbingBoundary3D::updateParameter (element = " << getTag() << ").\n"
                "Current stage = " << static_cast<int>(m_stage) << " (Stage_Absorbing).\n"
                "You cannot change the stage at this point!\n";
            exit(-1);
        }
    default:
        return -1;
    }
}

void ASDAbsorbingBoundary3D::addDisplacement(Vector& U)
{
    int counter = 0;
    for (Node* node : m_nodes) {
        const Vector& iU = node->getTrialDisp();
        for (int i = 0; i < iU.Size(); ++i)
            U(counter++) += iU(i);
    }
}

const Vector& ASDAbsorbingBoundary3D::getDisplacement()
{
    static Vector U;
    U.resize(m_num_dofs);
    U.Zero();
    addDisplacement(U);
    U.addVector(1.0, m_U0, -1.0);
    return U;
}

const Vector& ASDAbsorbingBoundary3D::getVelocity()
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

const Vector& ASDAbsorbingBoundary3D::getAcceleration()
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

void ASDAbsorbingBoundary3D::getElementSizes(double& lx, double& ly, double& lz, double& nx, double& ny)
{
    // get nodes
    Node* N0 = m_nodes[m_node_map[0]];
    Node* N1 = m_nodes[m_node_map[1]];
    Node* N2 = m_nodes[m_node_map[2]];
    Node* N4 = m_nodes[m_node_map[4]];

    // this should be always positive due to sorting...
    lz = std::abs(N1->getCrds()(2) - N0->getCrds()(2));

    // nx/ny positive if the outward normal vector from soil domain
    // points in the positive +X/+Y direction.
    // note: change the sign here: external forces on soil domain
    nx = ny = 1.0; // default

    // choose lx,ly,nx,ny: 17 combinations
    if (m_boundary & BND_LEFT) {
        if (m_boundary & BND_BACK) {
            // LEFT-BACK, LEFT-BACK-BOTTOM
            lx = std::abs(N4->getCrds()(0) - N0->getCrds()(0));
            ly = std::abs(N2->getCrds()(1) - N0->getCrds()(1));
            ny = -1.0;
        }
        else {
            // LEFT, LEFT-FRONT, LEFT-BOTTOM, LEFT-FRONT-BOTTOM
            lx = std::abs(N2->getCrds()(0) - N0->getCrds()(0));
            ly = std::abs(N4->getCrds()(1) - N0->getCrds()(1));
        }
    }
    else if (m_boundary & BND_RIGHT) {
        if (m_boundary & BND_BACK) {
            // RIGHT-BACK, RIGHT-BACK-BOTTOM
            lx = std::abs(N4->getCrds()(0) - N0->getCrds()(0));
            ly = std::abs(N2->getCrds()(1) - N0->getCrds()(1));
            ny = -1.0;
        }
        else {
            // RIGHT, RIGHT-FRONT, RIGHT-BOTTOM, RIGHT-FRONT-BOTTOM
            lx = std::abs(N2->getCrds()(0) - N0->getCrds()(0));
            ly = std::abs(N4->getCrds()(1) - N0->getCrds()(1));
            nx = -1.0;
        }
    }
    else if (m_boundary & BND_FRONT) {
        // FRONT FRONT-BOTTOM
        lx = std::abs(N4->getCrds()(0) - N0->getCrds()(0));
        ly = std::abs(N2->getCrds()(1) - N0->getCrds()(1));
    }
    else if (m_boundary & BND_BACK) {
        // BACK BACK-BOTTOM
        lx = std::abs(N4->getCrds()(0) - N0->getCrds()(0));
        ly = std::abs(N2->getCrds()(1) - N0->getCrds()(1));
        ny = -1.0;
    }
    else if (m_boundary & BND_BOTTOM) {
        // BOTTOM
        lx = std::abs(N2->getCrds()(0) - N0->getCrds()(0));
        ly = std::abs(N4->getCrds()(1) - N0->getCrds()(1));
    }

}

void ASDAbsorbingBoundary3D::updateStage()
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

void ASDAbsorbingBoundary3D::penaltyFactor(double& sp, double& mp)
{
    // get characteristic length
    double lx, ly, lz, nx, ny;
    getElementSizes(lx, ly, lz, nx, ny);
    double lch = std::cbrt(lx * ly * lz);
    // order of magnitude of the (approximate) max stiffness value
    // of an equivalent 3D solid element with this G
    int oom = static_cast<int>(std::round(std::log10(m_G * lch)));
    // ideal power = oom + 8 (8 being half the number of significant digits of a double floating number)
    int psp = oom + 8;
    int pmp = oom + 3;
    // compute the penalty factor
    sp = std::pow(10.0, psp);
    mp = std::pow(10.0, pmp);
}

void ASDAbsorbingBoundary3D::addKPenaltyStage0(Matrix& K)
{
    // Enforce constraints in stage = 0.

    // penalty factor
    double sp, mp;
    penaltyFactor(sp, mp);

    if (m_boundary == BND_BOTTOM) { // 1
        for (int i = 0; i < 4; ++i) {
            // fix Uz (interior and exterior)
            cfix(N_TOP[i], UZ, K, sp, m_dof_map);
            cfix(N_BOTTOM[i], UZ, K, sp, m_dof_map);
            // edof Ux and Uy (exterior->interior)
            cedof(N_BOTTOM[i], N_TOP[i], UX, K, mp, m_dof_map);
            cedof(N_BOTTOM[i], N_TOP[i], UY, K, mp, m_dof_map);
        }
    }
    else if (
        (m_boundary == BND_LEFT) || 
        (m_boundary == BND_RIGHT)) { // 3
        for (int i = 0; i < 4; ++i) {
            // fix Ux (interior and exterior)
            cfix(N_SS[i], UX, K, sp, m_dof_map);
            cfix(N_FF[i], UX, K, sp, m_dof_map);
            // edof Uy and Uz (exterior->interior)
            cedof(N_FF[i], N_SS[i], UY, K, mp, m_dof_map);
            cedof(N_FF[i], N_SS[i], UZ, K, mp, m_dof_map);
        }
    }
    else if (
        (m_boundary == BND_FRONT) || 
        (m_boundary == BND_BACK)) { // 5
        for (int i = 0; i < 4; ++i) {
            // fix Uy (interior and exterior)
            cfix(N_SS[i], UY, K, sp, m_dof_map);
            cfix(N_FF[i], UY, K, sp, m_dof_map);
            // edof Ux and Uz (exterior->interior)
            cedof(N_FF[i], N_SS[i], UX, K, mp, m_dof_map);
            cedof(N_FF[i], N_SS[i], UZ, K, mp, m_dof_map);
        }
    }
    else if (
        (m_boundary == (BND_LEFT | BND_FRONT)) ||
        (m_boundary == (BND_LEFT | BND_BACK)) ||
        (m_boundary == (BND_RIGHT | BND_FRONT)) ||
        (m_boundary == (BND_RIGHT | BND_BACK))) { // 9
        for (int i = 0; i < 2; ++i) {
            // fix Ux and Uy (interior)
            cfix(NVE_IN[i], UX, K, sp, m_dof_map);
            cfix(NVE_IN[i], UY, K, sp, m_dof_map);
            for (int j = 0; j < 3; ++j) {
                // fix Ux and Uy (exterior)
                cfix(NVE_OUT[i][j], UX, K, sp, m_dof_map);
                cfix(NVE_OUT[i][j], UY, K, sp, m_dof_map);
                // edof Uz (exterior->interior)
                cedof(NVE_OUT[i][j], NVE_IN[i], UZ, K, mp, m_dof_map);
            }
        }
    }
    else if (
        (m_boundary == (BND_BOTTOM | BND_LEFT)) ||
        (m_boundary == (BND_BOTTOM | BND_RIGHT))) { // 11
        for (int i = 0; i < 2; ++i) {
            // fix Ux and Uz (interior)
            cfix(NHE_IN[i], UX, K, sp, m_dof_map);
            cfix(NHE_IN[i], UZ, K, sp, m_dof_map);
            for (int j = 0; j < 3; ++j) {
                // fix Ux and Uz (exterior)
                cfix(NHE_OUT[i][j], UX, K, sp, m_dof_map);
                cfix(NHE_OUT[i][j], UZ, K, sp, m_dof_map);
                // edof Uy (exterior->interior)
                cedof(NHE_OUT[i][j], NHE_IN[i], UY, K, mp, m_dof_map);
            }
        }
    }
    else if (
        (m_boundary == (BND_BOTTOM | BND_FRONT)) ||
        (m_boundary == (BND_BOTTOM | BND_BACK))) { // 13
        for (int i = 0; i < 2; ++i) {
            // fix Uy and Uz (interior)
            cfix(NHE_IN[i], UY, K, sp, m_dof_map);
            cfix(NHE_IN[i], UZ, K, sp, m_dof_map);
            for (int j = 0; j < 3; ++j) {
                // fix Uy and Uz (exterior)
                cfix(NHE_OUT[i][j], UY, K, sp, m_dof_map);
                cfix(NHE_OUT[i][j], UZ, K, sp, m_dof_map);
                // edof Ux (exterior->interior)
                cedof(NHE_OUT[i][j], NHE_IN[i], UX, K, mp, m_dof_map);
            }
        }
    }
    else if (
        (m_boundary == (BND_BOTTOM | BND_LEFT | BND_FRONT)) ||
        (m_boundary == (BND_BOTTOM | BND_LEFT | BND_BACK)) ||
        (m_boundary == (BND_BOTTOM | BND_RIGHT | BND_FRONT)) ||
        (m_boundary == (BND_BOTTOM | BND_RIGHT | BND_BACK))) { // 17
        for (int i = 0; i < 8; ++i) {
            // fix ALL
            cfix(i, UX, K, sp, m_dof_map);
            cfix(i, UY, K, sp, m_dof_map);
            cfix(i, UZ, K, sp, m_dof_map);
        }
    }
}

void ASDAbsorbingBoundary3D::addRPenaltyStage0(Vector& R)
{
    // Enforce constraints in stage = 0.
    
    // skip if computing reactions
    if (m_is_computing_reactions)
        return;

    // penalty factor
    double sp, mp;
    penaltyFactor(sp, mp);

    // get displacement vector
    const Vector& U = getDisplacement();

    if (m_boundary == BND_BOTTOM) { // 1
        for (int i = 0; i < 4; ++i) {
            // fix Uz (interior and exterior)
            cfix(N_TOP[i], UZ, R, U, sp, m_dof_map);
            cfix(N_BOTTOM[i], UZ, R, U, sp, m_dof_map);
            // edof Ux and Uy (exterior->interior)
            cedof(N_BOTTOM[i], N_TOP[i], UX, R, U, mp, m_dof_map);
            cedof(N_BOTTOM[i], N_TOP[i], UY, R, U, mp, m_dof_map);
        }
    }
    else if (
        (m_boundary == BND_LEFT) ||
        (m_boundary == BND_RIGHT)) { // 3
        for (int i = 0; i < 4; ++i) {
            // fix Ux (interior and exterior)
            cfix(N_SS[i], UX, R, U, sp, m_dof_map);
            cfix(N_FF[i], UX, R, U, sp, m_dof_map);
            // edof Uy and Uz (exterior->interior)
            cedof(N_FF[i], N_SS[i], UY, R, U, mp, m_dof_map);
            cedof(N_FF[i], N_SS[i], UZ, R, U, mp, m_dof_map);
        }
    }
    else if (
        (m_boundary == BND_FRONT) ||
        (m_boundary == BND_BACK)) { // 5
        for (int i = 0; i < 4; ++i) {
            // fix Uy (interior and exterior)
            cfix(N_SS[i], UY, R, U, sp, m_dof_map);
            cfix(N_FF[i], UY, R, U, sp, m_dof_map);
            // edof Ux and Uz (exterior->interior)
            cedof(N_FF[i], N_SS[i], UX, R, U, mp, m_dof_map);
            cedof(N_FF[i], N_SS[i], UZ, R, U, mp, m_dof_map);
        }
    }
    else if (
        (m_boundary == (BND_LEFT | BND_FRONT)) ||
        (m_boundary == (BND_LEFT | BND_BACK)) ||
        (m_boundary == (BND_RIGHT | BND_FRONT)) ||
        (m_boundary == (BND_RIGHT | BND_BACK))) { // 9
        for (int i = 0; i < 2; ++i) {
            // fix Ux and Uy (interior)
            cfix(NVE_IN[i], UX, R, U, sp, m_dof_map);
            cfix(NVE_IN[i], UY, R, U, sp, m_dof_map);
            for (int j = 0; j < 3; ++j) {
                // fix Ux and Uy (exterior)
                cfix(NVE_OUT[i][j], UX, R, U, sp, m_dof_map);
                cfix(NVE_OUT[i][j], UY, R, U, sp, m_dof_map);
                // edof Uz (exterior->interior)
                cedof(NVE_OUT[i][j], NVE_IN[i], UZ, R, U, mp, m_dof_map);
            }
        }
    }
    else if (
        (m_boundary == (BND_BOTTOM | BND_LEFT)) ||
        (m_boundary == (BND_BOTTOM | BND_RIGHT))) { // 11
        for (int i = 0; i < 2; ++i) {
            // fix Ux and Uz (interior)
            cfix(NHE_IN[i], UX, R, U, sp, m_dof_map);
            cfix(NHE_IN[i], UZ, R, U, sp, m_dof_map);
            for (int j = 0; j < 3; ++j) {
                // fix Ux and Uz (exterior)
                cfix(NHE_OUT[i][j], UX, R, U, sp, m_dof_map);
                cfix(NHE_OUT[i][j], UZ, R, U, sp, m_dof_map);
                // edof Uy (exterior->interior)
                cedof(NHE_OUT[i][j], NHE_IN[i], UY, R, U, mp, m_dof_map);
            }
        }
    }
    else if (
        (m_boundary == (BND_BOTTOM | BND_FRONT)) ||
        (m_boundary == (BND_BOTTOM | BND_BACK))) { // 13
        for (int i = 0; i < 2; ++i) {
            // fix Uy and Uz (interior)
            cfix(NHE_IN[i], UY, R, U, sp, m_dof_map);
            cfix(NHE_IN[i], UZ, R, U, sp, m_dof_map);
            for (int j = 0; j < 3; ++j) {
                // fix Uy and Uz (exterior)
                cfix(NHE_OUT[i][j], UY, R, U, sp, m_dof_map);
                cfix(NHE_OUT[i][j], UZ, R, U, sp, m_dof_map);
                // edof Ux (exterior->interior)
                cedof(NHE_OUT[i][j], NHE_IN[i], UX, R, U, mp, m_dof_map);
            }
        }
    }
    else if (
        (m_boundary == (BND_BOTTOM | BND_LEFT | BND_FRONT)) ||
        (m_boundary == (BND_BOTTOM | BND_LEFT | BND_BACK)) ||
        (m_boundary == (BND_BOTTOM | BND_RIGHT | BND_FRONT)) ||
        (m_boundary == (BND_BOTTOM | BND_RIGHT | BND_BACK))) { // 17
        for (int i = 0; i < 8; ++i) {
            // fix ALL
            cfix(i, UX, R, U, sp, m_dof_map);
            cfix(i, UY, R, U, sp, m_dof_map);
            cfix(i, UZ, R, U, sp, m_dof_map);
        }
    }
}

void ASDAbsorbingBoundary3D::addKPenaltyStage1(Matrix& K)
{
    // Enforce constraints in stage = 1.
    // In this stage the edge opposed to the soil domain:
    // 1) is fixed in both X and Y and Z only on Horizontal boundary [SP constraint].

    // skip vertical boundary
    if (!(m_boundary & BND_BOTTOM))
        return;

    // penalty factor
    double sp, mp;
    penaltyFactor(sp, mp);

    // Fix Ux, Uy and Uz for all nodes not in the soil domain
    for (int i = 0; i < 4; ++i) {
        cfix(N_BOTTOM[i], UX, K, sp, m_dof_map);
        cfix(N_BOTTOM[i], UY, K, sp, m_dof_map);
        cfix(N_BOTTOM[i], UZ, K, sp, m_dof_map);
    }
}

void ASDAbsorbingBoundary3D::addRPenaltyStage1(Vector& R)
{
    // Enforce constraints in stage = 1.
    // In this stage the edge opposed to the soil domain:
    // 1) is fixed in both X and Y and Z only on Horizontal boundary [SP constraint].

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

    // Fix Ux, Uy and Uz for all nodes not in the soil domain
    for (int i = 0; i < 4; ++i) {
        cfix(N_BOTTOM[i], UX, R, U, sp, m_dof_map);
        cfix(N_BOTTOM[i], UY, R, U, sp, m_dof_map);
        cfix(N_BOTTOM[i], UZ, R, U, sp, m_dof_map);
    }
}

void ASDAbsorbingBoundary3D::addRReactions(Vector& R)
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

void ASDAbsorbingBoundary3D::addMff(Matrix& M, double scale)
{
    // Add the mass terms due to the Free-Field columns, thus only 
    // on vertical boundaries
    // The optional scale parameter can be used while assembling
    // the damping due to the free-field

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // get the whole element free-field mass
    double lx, ly, lz, nx, ny;
    getElementSizes(lx, ly, lz, nx, ny);
    double hm = scale * m_rho * lx * ly * lz;

    if (
        (m_boundary == BND_LEFT) || 
        (m_boundary == BND_RIGHT) ||
        (m_boundary == BND_FRONT) ||
        (m_boundary == BND_BACK)) {
        // on pure side boundaries we have 4 nodes for free-field
        hm /= 4.0;
        for (int node_id : FF_ID_SIDE) {
            int iL = node_id * 3;
            int iG = m_dof_map(iL);
            for (int dof = 0; dof < 3; ++dof) {
                M(iG + dof, iG + dof) += hm;
            }
        }
    }
    else if (
        (m_boundary == (BND_LEFT | BND_FRONT)) ||
        (m_boundary == (BND_LEFT | BND_BACK)) ||
        (m_boundary == (BND_RIGHT | BND_FRONT)) ||
        (m_boundary == (BND_RIGHT | BND_BACK))) {
        // on vertical edges only 2 nodes have the whole free-field
        hm /= 2.0;
        for (int node_id : FF_ID_EDGE) {
            int iL = node_id * 3;
            int iG = m_dof_map(iL);
            for (int dof = 0; dof < 3; ++dof) {
                M(iG + dof, iG + dof) += hm;
            }
        }
    }

}

void ASDAbsorbingBoundary3D::addRMff(Vector& R)
{
    // Add the mass terms due to the Free-Field columns, thus only 
    // on vertical boundaries

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // get acceleration
    const Vector& A = getAcceleration();

    // get the whole element free-field mass
    double lx, ly, lz, nx, ny;
    getElementSizes(lx, ly, lz, nx, ny);
    double hm = m_rho * lx * ly * lz;

    if (
        (m_boundary == BND_LEFT) ||
        (m_boundary == BND_RIGHT) ||
        (m_boundary == BND_FRONT) ||
        (m_boundary == BND_BACK)) {
        // on pure side boundaries we have 4 nodes for free-field
        hm /= 4.0;
        for (int node_id : FF_ID_SIDE) {
            int iL = node_id * 3;
            int iG = m_dof_map(iL);
            for (int dof = 0; dof < 3; ++dof) {
                R(iG + dof) += hm * A(iG + dof);
            }
        }
    }
    else if (
        (m_boundary == (BND_LEFT | BND_FRONT)) ||
        (m_boundary == (BND_LEFT | BND_BACK)) ||
        (m_boundary == (BND_RIGHT | BND_FRONT)) ||
        (m_boundary == (BND_RIGHT | BND_BACK))) {
        // on vertical edges only 2 nodes have the whole free-field
        hm /= 2.0;
        for (int node_id : FF_ID_EDGE) {
            int iL = node_id * 3;
            int iG = m_dof_map(iL);
            for (int dof = 0; dof < 3; ++dof) {
                R(iG + dof) += hm * A(iG + dof);
            }
        }
    }
}

const ID& ASDAbsorbingBoundary3D::ffMapping()
{
    static ID m(24);

    // h8 is formed with the original order, but with a local dofset
    int counter = 0;
    for (int i = 0; i < 8; ++i) {
        Node* node = m_nodes[static_cast<std::size_t>(i)];
        int j = i * 3;
        m(j) = counter;
        m(j + 1) = counter + 1;
        m(j + 2) = counter + 2;
        counter += node->getNumberDOF();
    }

    int n0 = static_cast<int>(m_node_map[0]);
    int n1 = static_cast<int>(m_node_map[1]);
    int n2 = static_cast<int>(m_node_map[2]);
    int n3 = static_cast<int>(m_node_map[3]);
    int n4 = static_cast<int>(m_node_map[4]);
    int n5 = static_cast<int>(m_node_map[5]);
    int n6 = static_cast<int>(m_node_map[6]);
    int n7 = static_cast<int>(m_node_map[7]);

    auto map = [](int na, int nb) {
        // all DOFs of na are set equal to all DOFs of nb
        m(na * 3) = m(nb * 3);
        m(na * 3 + 1) = m(nb * 3 + 1);
        m(na * 3 + 2) = m(nb * 3 + 2);
    };

    if (
        (m_boundary == (BND_LEFT | BND_FRONT)) ||
        (m_boundary == (BND_LEFT | BND_BACK)) ||
        (m_boundary == (BND_RIGHT | BND_FRONT)) ||
        (m_boundary == (BND_RIGHT | BND_BACK))
        ) {
        map(n2, n0);
        map(n6, n0);
        map(n4, n0);
        map(n3, n1);
        map(n5, n1);
        map(n7, n1);
    }
    else {
        map(n2, n0);
        map(n3, n1);
        map(n6, n4);
        map(n7, n5);
    }

    return m;
}

void ASDAbsorbingBoundary3D::addKff(Matrix& K, double scale)
{
    // Add the stiffness matrix of the free-field column.
    // Only on vertical boundaries.
    // The optional scale parameter can be used while assembling
    // the damping due to the free-field

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // transformation vector for static condensation
    const ID& T = ffMapping();

    // H8 node matrix
    static Matrix P(3, 8);
    H8_nodeMatrix(m_nodes, m_node_map, P);

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);
    static Matrix E(6, 6);
    H8_C0(lam, mu, E);

    // other data
    static Matrix dN(8, 3);
    static Matrix J(3, 3);
    static Matrix invJ(3, 3);
    double detJ = 0.0;
    static Matrix dNdX(8, 3);
    static Matrix B(6, 24);
    static Matrix BB;
    BB.resize(6, m_num_dofs);

    // gauss integration
    for (int gauss_id = 0; gauss_id < H8_GW.size(); ++gauss_id) {

        // quadrature data
        double gx = H8_GX[gauss_id];
        double gy = H8_GY[gauss_id];
        double gz = H8_GZ[gauss_id];
        double gw = H8_GW[gauss_id];

        // shape function derivatives
        H8_dN(gx, gy, gz, dN);
        J.addMatrixProduct(0.0, P, dN, 1.0);
        detJ = H8_det3(J);
        J.Invert(invJ);
        dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

        // b-matrix of the h8
        H8_Bmatrix(dNdX, B);

        // condensate to obtain the b-matrix of the 4-node free-field
        BB.Zero();
        for (int j = 0; j < 24; ++j) {
            int jj = T(j);
            for (int i = 0; i < 6; ++i) {
                BB(i, jj) += B(i, j);
            }
        }

        // integrate
        double dV = detJ * gw * scale;
        K.addMatrixTripleProduct(1.0, BB, E, dV);
    }
}

void ASDAbsorbingBoundary3D::addRff(Vector& R)
{
    // Add the residual of the free-field column.

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // transformation vector for static condensation
    const ID& T = ffMapping();

    // get displacement
    const Vector& U = getDisplacement();

    // H8 node matrix
    static Matrix P(3, 8);
    H8_nodeMatrix(m_nodes, m_node_map, P);

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);
    static Matrix E(6, 6);
    H8_C0(lam, mu, E);

    // other data
    static Matrix dN(8, 3);
    static Matrix J(3, 3);
    static Matrix invJ(3, 3);
    double detJ = 0.0;
    static Matrix dNdX(8, 3);
    static Matrix B(6, 24);
    static Matrix BB;
    BB.resize(6, m_num_dofs);
    static Vector strain(6);
    static Vector stress(6);

    // gauss integration
    for (int gauss_id = 0; gauss_id < H8_GW.size(); ++gauss_id) {

        // quadrature data
        double gx = H8_GX[gauss_id];
        double gy = H8_GY[gauss_id];
        double gz = H8_GZ[gauss_id];
        double gw = H8_GW[gauss_id];

        // shape function derivatives
        H8_dN(gx, gy, gz, dN);
        J.addMatrixProduct(0.0, P, dN, 1.0);
        detJ = H8_det3(J);
        J.Invert(invJ);
        dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

        // b-matrix of the h8
        H8_Bmatrix(dNdX, B);

        // condensate to obtain the b-matrix of the 4-node free-field
        BB.Zero();
        for (int j = 0; j < 24; ++j) {
            int jj = T(j);
            for (int i = 0; i < 6; ++i) {
                BB(i, jj) += B(i, j);
            }
        }

        // strain vector
        strain.addMatrixVector(0.0, BB, U, 1.0);

        // stress vector
        stress.addMatrixVector(0.0, E, strain, 1.0);

        // integrate
        double dV = detJ * gw;
        R.addMatrixTransposeVector(1.0, BB, stress, dV);
    }
}

const Matrix& ASDAbsorbingBoundary3D::computeNmatrix()
{
    // N matrix to obtain the traction vector
    // from the free-field stress 
    // applied on the soil nodes
    static Matrix N;
    N.resize(m_num_dofs, 6);
    N.Zero();

    // aux vector for normal from free-field to soil
    // (it should be soil-to-free-field, but we invert the sign
    // here because these are external forces!
    static Vector n(3);

    // choose ff-ss pairs (same as per LK dashpots)
    const std::vector<LKnodes>& pairs = LKselectPairs(m_boundary);

    // obtain element sizes
    double lx, ly, lz, nx, ny;
    getElementSizes(lx, ly, lz, nx, ny);

    // process each pair
    for (const LKnodes& ip : pairs) {

        // get nodes
        const Node* nff = m_nodes[m_node_map[static_cast<std::size_t>(ip.ff)]];
        const Node* nss = m_nodes[m_node_map[static_cast<std::size_t>(ip.ss)]];

        // get direction (nss - nff)
        n.addVector(0.0, nss->getCrds(), 1.0);
        n.addVector(1.0, nff->getCrds(), -1.0);
        if (n.Normalize() != 0) {
            opserr << "ASDAbsordbinBoundary3D Error: distance between nodes " << nff->getTag() << " and " << nss->getTag() << " is ZERO!\n";
            exit(-1);
        }

        // compute lumping coefficients (area normal to n)
        // based on the direction.
        // also include the weight of the LK pair, to account for side elements with only 2 ff nodes
        double dA = 0.0;
        if (std::abs(n(0)) > 0.99) { // X
            dA = ly * lz * ip.w / 8.0 / 4.0;
        }
        else if (std::abs(n(1)) > 0.99) { // Y
            dA = lx * lz * ip.w / 8.0 / 4.0;
        }
        else { // others
            opserr << "ASDAbsordbinBoundary3D Error: normal vector can be only X or Y, not " << n << "\n";
            exit(-1);
        }

        // scale the normal vector with the lumping coefficient
        // so that we don't do it during gauss integration
        n *= dA;

        // fill the N matrix with the components of the normal vector,
        // properly places to account for voigt notation
        int ix = ip.ss * 3;
        int iy = ix + 1;
        int iz = ix + 2;
        ix = m_dof_map(ix);
        iy = m_dof_map(iy);
        iz = m_dof_map(iz);
        N(ix, 0) += n(0);
        N(ix, 3) += n(1);
        N(ix, 5) += n(2);
        N(iy, 1) += n(1);
        N(iy, 3) += n(0);
        N(iy, 4) += n(2);
        N(iz, 2) += n(2);
        N(iz, 4) += n(1);
        N(iz, 5) += n(0);
    }

    // done
    return N;
}

void ASDAbsorbingBoundary3D::addKffToSoil(Matrix& K)
{
    // Add the stiffness matrix of the forces transfered from the
    // free-field column to the soil domain.
    // Only on vertical boundaries

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // transformation vector for static condensation
    const ID& T = ffMapping();

    // H8 node matrix
    static Matrix P(3, 8);
    H8_nodeMatrix(m_nodes, m_node_map, P);

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);
    static Matrix E(6, 6);
    H8_C0(lam, mu, E);

    // N matrix
    const Matrix& N = computeNmatrix();

    // other data
    static Matrix dN(8, 3);
    static Matrix J(3, 3);
    static Matrix invJ(3, 3);
    double detJ = 0.0;
    static Matrix dNdX(8, 3);
    static Matrix B(6, 24);
    static Matrix BB;
    BB.resize(6, m_num_dofs);
    static Matrix NE;
    NE.resize(m_num_dofs, 6);

    // gauss integration
    for (int gauss_id = 0; gauss_id < H8_GW.size(); ++gauss_id) {

        // quadrature data
        double gx = H8_GX[gauss_id];
        double gy = H8_GY[gauss_id];
        double gz = H8_GZ[gauss_id];
        double gw = H8_GW[gauss_id];

        // shape function derivatives
        H8_dN(gx, gy, gz, dN);
        J.addMatrixProduct(0.0, P, dN, 1.0);
        J.Invert(invJ);
        dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

        // b-matrix of the h8
        H8_Bmatrix(dNdX, B);

        // condensate to obtain the b-matrix of the 4-node free-field
        BB.Zero();
        for (int j = 0; j < 24; ++j) {
            int jj = T(j);
            for (int i = 0; i < 6; ++i) {
                BB(i, jj) += B(i, j);
            }
        }

        // K += N*E*BB
        NE.addMatrixProduct(0.0, N, E, 1.0);
        K.addMatrixProduct(1.0, NE, BB, 1.0);
    }
}

void ASDAbsorbingBoundary3D::addRffToSoil(Vector& R)
{
    // Add the forces transfered from the
    // free-field column to the soil domain.
    // Only on vertical boundaries

    // skip horizontal boundary
    if (m_boundary & BND_BOTTOM)
        return;

    // transformation vector for static condensation
    const ID& T = ffMapping();

    // get displacement
    const Vector& U = getDisplacement();

    // H8 node matrix
    static Matrix P(3, 8);
    H8_nodeMatrix(m_nodes, m_node_map, P);

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);
    static Matrix E(6, 6);
    H8_C0(lam, mu, E);

    // N matrix
    const Matrix& N = computeNmatrix();

    // other data
    static Matrix dN(8, 3);
    static Matrix J(3, 3);
    static Matrix invJ(3, 3);
    double detJ = 0.0;
    static Matrix dNdX(8, 3);
    static Matrix B(6, 24);
    static Matrix BB(6, m_num_dofs);
    static Vector strain(6);
    static Vector stress(6);

    // gauss integration
    for (int gauss_id = 0; gauss_id < H8_GW.size(); ++gauss_id) {

        // quadrature data
        double gx = H8_GX[gauss_id];
        double gy = H8_GY[gauss_id];
        double gz = H8_GZ[gauss_id];
        double gw = H8_GW[gauss_id];

        // shape function derivatives
        H8_dN(gx, gy, gz, dN);
        J.addMatrixProduct(0.0, P, dN, 1.0);
        J.Invert(invJ);
        dNdX.addMatrixProduct(0.0, dN, invJ, 1.0);

        // b-matrix of the h8
        H8_Bmatrix(dNdX, B);

        // condensate to obtain the b-matrix of the 4-node free-field
        BB.Zero();
        for (int j = 0; j < 24; ++j) {
            int jj = T(j);
            for (int i = 0; i < 6; ++i) {
                BB(i, jj) += B(i, j);
            }
        }

        // strain vector
        strain.addMatrixVector(0.0, BB, U, 1.0);

        // stress vector
        stress.addMatrixVector(0.0, E, strain, 1.0);

        // integrate
        R.addMatrixVector(1.0, N, stress, 1.0);
    }
}

void ASDAbsorbingBoundary3D::getDampParam(double& alpha, double& beta)
{
    alpha = alphaM;
    beta = betaK;
    if (beta == 0.0) {
        beta = betaK0;
        if (beta == 0)
            beta = betaKc;
    }
}

void ASDAbsorbingBoundary3D::addCff(Matrix& C)
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

void ASDAbsorbingBoundary3D::addRCff(Vector& R)
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

void ASDAbsorbingBoundary3D::addClk(Matrix& C)
{
    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // wave velocities
    double vp = std::sqrt((lam + 2.0 * mu) / m_rho);
    double vs = std::sqrt(mu / m_rho);

    // sizes
    double lx, ly, lz, nx, ny;
    getElementSizes(lx, ly, lz, nx, ny);

    // divide all sizes for initial lumping
    lx *= 0.5;
    ly *= 0.5;
    lz *= 0.5;

    // choose pairs
    const std::vector<LKnodes>& pairs = LKselectPairs(m_boundary);

    // aux
    static Vector direction(3);
    static Vector coeff(3);

    // process each pair
    for (const LKnodes& ip : pairs) {

        // get nodes
        const Node* nff = m_nodes[m_node_map[static_cast<std::size_t>(ip.ff)]];
        const Node* nss = m_nodes[m_node_map[static_cast<std::size_t>(ip.ss)]];

        // get direction (nss - nff)
        direction.addVector(0.0, nss->getCrds(), 1.0);
        direction.addVector(1.0, nff->getCrds(), -1.0);
        if (direction.Normalize() != 0) {
            opserr << "ASDAbsordbinBoundary3D Error: distance between nodes " << nff->getTag() << " and " << nss->getTag() << " is ZERO!\n";
            exit(-1);
        }

        // compute lumped coefficients based on direction
        // (put the minus sign here: external forces acting on soil domain)
        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;
        if (std::abs(direction(0)) > 0.99) { // X
            double area = ly * lz * ip.w;
            cx = -vp * m_rho * area;
            cy = cz = -vs * m_rho * area;
        }
        else if (std::abs(direction(1)) > 0.99) { // Y
            double area = lx * lz * ip.w;
            cy = -vp * m_rho * area;
            cx = cz = -vs * m_rho * area;
        }
        else { // Z
            double area = lx * ly * ip.w;
            cz = -vp * m_rho * area;
            cx = cy = -vs * m_rho * area;
        }
        coeff(0) = cx;
        coeff(1) = cy;
        coeff(2) = cz;

        // fill
        for (int dof = 0; dof < 3; ++dof) {
            int ff_local = ip.ff * 3 + dof;
            int ss_local = ip.ss * 3 + dof;
            int ff_global = m_dof_map(ff_local);
            int ss_global = m_dof_map(ss_local);
            C(ss_global, ff_global) += coeff(dof);
            C(ss_global, ss_global) -= coeff(dof);
        }
    }
}

void ASDAbsorbingBoundary3D::addRlk(Vector& R)
{
    // get velocity
    const Vector& V = getVelocity();

    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // wave velocities
    double vp = std::sqrt((lam + 2.0 * mu) / m_rho);
    double vs = std::sqrt(mu / m_rho);

    // sizes
    double lx, ly, lz, nx, ny;
    getElementSizes(lx, ly, lz, nx, ny);

    // divide all sizes for initial lumping
    lx *= 0.5;
    ly *= 0.5;
    lz *= 0.5;

    // choose pairs
    const std::vector<LKnodes>& pairs = LKselectPairs(m_boundary);

    // aux
    static Vector direction(3);
    static Vector coeff(3);

    // process each pair
    for (const LKnodes& ip : pairs) {

        // get nodes
        const Node* nff = m_nodes[m_node_map[static_cast<std::size_t>(ip.ff)]];
        const Node* nss = m_nodes[m_node_map[static_cast<std::size_t>(ip.ss)]];

        // get direction (nss - nff)
        direction.addVector(0.0, nss->getCrds(), 1.0);
        direction.addVector(1.0, nff->getCrds(), -1.0);
        if (direction.Normalize() != 0) {
            opserr << "ASDAbsordbinBoundary3D Error: distance between nodes " << nff->getTag() << " and " << nss->getTag() << " is ZERO!\n";
            exit(-1);
        }

        // compute lumped coefficients based on direction
        // (put the minus sign here: external forces acting on soil domain)
        double cx = 0.0;
        double cy = 0.0;
        double cz = 0.0;
        if (std::abs(direction(0)) > 0.99) { // X
            double area = ly * lz * ip.w;
            cx = -vp * m_rho * area;
            cy = cz = -vs * m_rho * area;
        }
        else if (std::abs(direction(1)) > 0.99) { // Y
            double area = lx * lz * ip.w;
            cy = -vp * m_rho * area;
            cx = cz = -vs * m_rho * area;
        }
        else { // Z
            double area = lx * ly * ip.w;
            cz = -vp * m_rho * area;
            cx = cy = -vs * m_rho * area;
        }
        coeff(0) = cx;
        coeff(1) = cy;
        coeff(2) = cz;

        // fill
        for (int dof = 0; dof < 3; ++dof) {
            int ff_local = ip.ff * 3 + dof;
            int ss_local = ip.ss * 3 + dof;
            int ff_global = m_dof_map(ff_local);
            int ss_global = m_dof_map(ss_local);
            R(ss_global) += coeff(dof) * (V(ff_global) - V(ss_global));
        }
    }
}

void ASDAbsorbingBoundary3D::addBaseActions(Vector& R)
{
    // skip vertical boundaries
    if (!(m_boundary & BND_BOTTOM))
        return;

    // compute velocities from time series
    auto get_vel = [this](TimeSeries* ts) {
        if (ts == nullptr)
            return 0.0;
        Domain* domain = getDomain();
        if (domain == nullptr) {
            opserr << "ASDAbsorbingBoundary3D Error: cannot get domain!\n";
            exit(-1);
        }
        return ts->getFactor(domain->getCurrentTime());
    };
    double vx = get_vel(m_tsx);
    double vy = get_vel(m_tsy);
    double vz = get_vel(m_tsz);

    // quick return
    if ((vx == 0.0) && (vy == 0.0) && (vz == 0.0))
        return;
    
    // elasticity constants
    double mu = m_G;
    double lam = 2.0 * mu * m_v / (1.0 - 2.0 * m_v);

    // wave velocities
    double vp = std::sqrt((lam + 2.0 * mu) / m_rho);
    double vs = std::sqrt(mu / m_rho);

    // sizes
    double lx, ly, lz, nx, ny;
    getElementSizes(lx, ly, lz, nx, ny);

    // lumped coefficients:
    // - put the minus sign here: external forces acting on soil domain
    double ap = -vp * m_rho * lx * ly / 4.0;
    double as = -vs * m_rho * lx * ly / 4.0;

    // lumped forces (2: accounts for half force absorbed by dashpots)
    std::array<double, 3> f = { {2.0 * as * vx, 2.0 * as * vy, 2.0 * ap * vz} };

    // add forces on top nodes
    static ID loaded_nodes(4);
    if ((m_boundary & BND_LEFT) || (m_boundary & BND_RIGHT)) {
        if ((m_boundary & BND_FRONT) || (m_boundary & BND_BACK)) {
            // corner: 4 contributions on same node (1)
            loaded_nodes(0) = loaded_nodes(1) = loaded_nodes(2) = loaded_nodes(3) = 1;
        }
        else {
            // side: 2 contributions on each of the 2 nodes (1 and 5)
            loaded_nodes(0) = loaded_nodes(1) = 1;
            loaded_nodes(2) = loaded_nodes(3) = 5;
        }
    }
    else {
        if ((m_boundary & BND_FRONT) || (m_boundary & BND_BACK)) {
            // side: 2 contributions on each of the 2 nodes (1 and 5)
            loaded_nodes(0) = loaded_nodes(1) = 1;
            loaded_nodes(2) = loaded_nodes(3) = 5;
        }
        else {
            // bottom only: 1 contribution on each of the 4 nodes
            loaded_nodes(0) = 1;
            loaded_nodes(1) = 3;
            loaded_nodes(2) = 5;
            loaded_nodes(3) = 7;
        }
    }
    for (int i = 0; i < 4; ++i) {
        int node = loaded_nodes[i];
        for (int j = 0; j < 3; ++j) {
            int iL = node * 3 + j;
            int iG = m_dof_map(iL);
            R(iG) += f[j];
        }
    }
}





