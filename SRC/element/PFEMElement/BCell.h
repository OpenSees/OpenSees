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

#ifndef BCell_h
#define BCell_h

#include "BackgroundDef.h"

class BNode;

// BACKGROUND_FLUID - a grid fluid cell
// BACKGROUND_STRUCTURE - a structural cell, which should have no particles
class BCell {
   private:
    VParticle pts;
    BackgroundType type;
    std::vector<BNode*> bnodes;
    std::vector<VInt> bindices;

   public:
    BCell();
    void add(Particle* pt);
    void setType(BackgroundType t);
    BackgroundType getType() const;
    VParticle& getPts();
    std::vector<BNode*>& getNodes();
    std::vector<VInt>& getIndices();

    void addNode(BNode* bnode, const VInt& index);

    void clearParticles();
};

#endif