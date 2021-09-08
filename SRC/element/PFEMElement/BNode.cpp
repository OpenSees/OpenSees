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

// Written: Minjie Zhu
//
//

#include "BNode.h"

BNode::BNode()
    : tags(), crdsn(), vn(), dvn(), pn(), dpn(), type(BACKGROUND_FLUID) {}

void BNode::addNode(int tag, const VDouble& crds, const VDouble& v,
                    const VDouble& dv, double p, double dp,
                    BackgroundType tp, int id) {
    tags.push_back(tag);
    crdsn.push_back(crds);
    vn.push_back(v);
    dvn.push_back(dv);
    pn.push_back(p);
    dpn.push_back(dp);
    type = tp;
    sid.push_back(id);
}

void BNode::clear() {
    tags.clear();
    crdsn.clear();
    vn.clear();
    dvn.clear();
    pn.clear();
    dpn.clear();
    type = BACKGROUND_FLUID;
    sid.clear();
}

void BNode::setType(BackgroundType t) { type = t; }

int BNode::size() const { return (int)tags.size(); }

VInt& BNode::getTags() { return tags; }

VVDouble& BNode::getCrds() { return crdsn; }

VVDouble& BNode::getVel() { return vn; }

VVDouble& BNode::getAccel() { return dvn; }

VDouble& BNode::getPressure() { return pn; }

VDouble& BNode::getPdot() { return dpn; }

BackgroundType BNode::getType() { return type; }

VInt& BNode::getSid() { return sid; }