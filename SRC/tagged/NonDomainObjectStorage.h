/* ******************************************************************
***
**    OpenSees - Open System for Earthquake Engineering Simulation **
**          Pacific Earthquake Engineering Research Center **
** **
** **
** (C) Copyright 1999, The Regents of the University of California **
** All Rights Reserved. **
** **
** Commercial use of this program without express permission of the **
** University of California, Berkeley, is strictly prohibited.  See **
** file 'COPYRIGHT'  in main directory for information on usage and **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES. **
** **
** Developed by: **
**   Frank McKenna (fmckenna@ce.berkeley.edu) **
**   Gregory L. Fenves (fenves@ce.berkeley.edu) **
**   Filip C. Filippou (filippou@ce.berkeley.edu) **
** **
** ******************************************************************
*/

// Minjie

// This is to store all objects that are not stored in the domain
// such as, uniaxial material, nDMaterial, sections, timeSeries,
// crdTransf, beam integration, limitCurves, damage models,
// backbone, etc.
// Because they are currently stored in their own .cpp files
// and are hard to keep track of them.
// The Tcl version will still use the old form of storage
// The Python version will move them here.

#ifndef NonDomainObjectStorage_H
#define NonDomainObjectStorage_H

#include <map>

class NonDomainObjectStorage {
   private:
   public:
    NonDomainObjectStorage() {}
    ~NonDomainObjectStorage() { clear(); }
    void clear() {}
};

#endif