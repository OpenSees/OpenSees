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
#ifndef Pipe_h
#define Pipe_h

#include <Domain.h>
#include <ElasticBeam3d.h>
#include <Matrix.h>
#include <PipeMaterial.h>
#include <PipeSection.h>
#include <Vector.h>
#include <elementAPI.h>

#include <vector>

class Node;

// straight pipe element
class Pipe : public ElasticBeam3d {
   private:
    PipeMaterial *theMat;
    PipeSection *theSect;

    double alp;       // thermal expansion coefficient
    double T0;        // stress free temperature
    double pressure;  // internal pressure

   public:
    Pipe();
    Pipe(int tag, int Nd1, int Nd2, CrdTransf &theTransf,
         PipeMaterial &mat, PipeSection &sect, double to = 0.0,
         double pre = 0.0, int cMass = 0, int releasez = 0,
         int releasey = 0);

    ~Pipe();

    const char *getClassType(void) const;

    void setDomain(Domain *theDomain);

    void zeroLoad(void);

   private:
    double aveTemp();
    int updateSectionData();
    int updateMaterialData();
};

#endif
