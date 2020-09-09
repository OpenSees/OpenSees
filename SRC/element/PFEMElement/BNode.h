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

#ifndef BNode_h
#define BNode_h

#include "BackgroundDef.h"

// BACKGROUND_FLUID - a grid fluid node
// BACKGROUND_STRUCTURE - a structural node
// BACKGROUND_FIXED - a fixed grid fluid node
class BNode {
   private:
    VInt tags;
    VVDouble crdsn;
    VVDouble vn;
    VVDouble dvn;
    VDouble pn;
    VDouble dpn;
    BackgroundType type;
    VInt sid;  // structure id, <0:fluid, >0:structure, =0:not in
               // contact

   public:
    BNode();
    void addNode(int tag, const VDouble& crds, const VDouble& v,
                 const VDouble& dv, double p, double dp, BackgroundType tp,
                 int id = -1);
    void clear();
    void setType(BackgroundType t);

    int size() const;
    VInt& getTags();
    VVDouble& getCrds();
    VVDouble& getVel();
    VVDouble& getAccel();
    VDouble& getPressure();
    VDouble& getPdot();
    BackgroundType getType();
    VInt& getSid();
};

#endif