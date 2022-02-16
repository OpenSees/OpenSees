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

#include "BCell.h"

#include "BNode.h"

BCell::BCell() : pts(), type(BACKGROUND_FLUID), bnodes(), bindices() {}

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

void BCell::clearParticles() { pts.clear(); }