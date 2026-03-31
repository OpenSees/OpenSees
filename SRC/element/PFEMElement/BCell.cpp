/* ******************************************************************
***
**    OpenSees - Open System for Earthquake Engineering
*Simulation    **
**          Pacific Earthquake Engineering Research Center
***
** **
** **
** (C) Copyright 1999, The Regents of the University of
*California    **
** All Rights Reserved. **
** **
** Commercial use of this program without express permission
*of the   **
** University of California, Berkeley, is strictly
*prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on
*usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.
***
** **
** Developed by: **
**   Frank McKenna (fmckenna@ce.berkeley.edu) **
**   Gregory L. Fenves (fenves@ce.berkeley.edu) **
**   Filip C. Filippou (filippou@ce.berkeley.edu) **
** **
** ****************************************************************** */

// Written: Minjie Zhu
//
//

#include "BCell.h"

#include <Domain.h>
#include <Node.h>
#include <Pressure_Constraint.h>
#include <elementAPI.h>

#include "BNode.h"

BCell::BCell()
    : pts(),
      type(BACKGROUND_FLUID),
      bnodes(),
      bindices(),
      center_node(0),
      center_pc(0) {}


void BCell::removeCenterNode() {
    auto* domain = OPS_GetDomain();
    if (center_node != 0) {
        if (domain != 0) {
            Node* nd = domain->removeNode(center_node->getTag());
            if (nd != 0) {
                delete nd;
            }
        }
        center_node = 0;
    }
    if (center_pc != 0) {
        if (domain != 0) {
            auto* pc = domain->removePressure_Constraint(
                center_pc->getTag());
            if (pc != 0) {
                delete pc;
            }
        }
        center_pc = 0;
    }
}

void BCell::add(Particle* pt) { pts.push_back(pt); }

void BCell::setType(BackgroundType t) { type = t; }

BackgroundType BCell::getType() const { return type; }

VParticle& BCell::getPts() { return pts; }

std::vector<BNode*>& BCell::getNodes() { return bnodes; }

std::vector<VInt>& BCell::getIndices() { return bindices; }

void BCell::addNode(BNode* bnode, const VInt& index) {
    bnodes.push_back(bnode);
    bindices.push_back(index);
}

Node* BCell::setCenterNode(int new_tag, int new_p_tag) {
    if (bnodes.size() != 4 && bnodes.size() != 8) {
        return 0;
    }
    removeCenterNode();
    Domain* domain = OPS_GetDomain();
    if (domain == 0) {
        return 0;
    }
    int ndm = OPS_GetNDM();

    // compute the center node
    VDouble crds(ndm), vel(ndm), accel(ndm);
    double pressure = 0.0, dp = 0.0;
    int count = 0;
    for (int i = 0; i < (int)bnodes.size(); ++i) {
        if (bnodes[i] != 0 &&
            bnodes[i]->getTags().size() > 0) {
            crds += bnodes[i]->getCrds()[0];
            vel += bnodes[i]->getVel()[0];
            accel += bnodes[i]->getAccel()[0];
            pressure += bnodes[i]->getPressure()[0];
            dp += bnodes[i]->getPdot()[0];
            count++;
        }
    }
    crds /= count;
    vel /= count;
    accel /= count;
    pressure /= count;
    dp /= count;

    // create the center node
    if (ndm == 2) {
        center_node =
            new Node(new_tag, ndm, crds[0], crds[1]);
    } else if (ndm == 3) {
        center_node = new Node(new_tag, ndm, crds[0],
                               crds[1], crds[2]);
    }
    if (center_node == 0) {
        return 0;
    }
    if (domain->addNode(center_node) == false) {
        delete center_node;
        center_node = 0;
        return 0;
    }

    // set the center node
    Vector vec;
    toVector(vel, vec);
    center_node->setTrialVel(vec);
    toVector(accel, vec);
    center_node->setTrialAccel(vec);
    center_node->commitState();

    // set the pressure constraint
    center_pc = domain->getPressure_Constraint(
        center_node->getTag());
    Node* pnode = 0;
    if (center_pc != 0) {
        pnode = center_pc->getPressureNode();
        if (pnode == 0) {
            return 0;
        }
    } else {
        // create pressure node
        if (ndm == 2) {
            pnode =
                new Node(new_p_tag, 1, crds[0], crds[1]);
        } else if (ndm == 3) {
            pnode = new Node(new_p_tag, 1, crds[0], crds[1],
                             crds[2]);
        }
        if (pnode == 0) {
            return 0;
        }
        if (domain->addNode(pnode) == false) {
            delete pnode;
            pnode = 0;
            return 0;
        }

        // create pressure constraint
        center_pc = new Pressure_Constraint(center_node->getTag(),
                                     pnode->getTag());
        if (center_pc == 0) {
            return 0;
        }
        domain->addPressure_Constraint(center_pc);
    }
    center_pc->setPressure(pressure);
    center_pc->setPdot(dp);

    return center_node;
}

Node* BCell::getCenterNode() { return center_node; }

void BCell::clearParticles() { pts.clear(); }