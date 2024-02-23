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
#include "CurvedPipe.h"

#include <CrdTransf.h>

void *OPS_CurvedPipeElement() {
    // check inputs
    if (OPS_GetNumRemainingInputArgs() < 10) {
        opserr << "Invalid #args,  want: element CurvedPipe "
                  "tag? nd1? nd2? transfTag? pipeMatTag? pipeSecTag?"
                  "R? xC? yC? zC?"
                  "<-T0 T0? -p p? -cMass? -releasey releasey? "
                  "-releasez releasez?>\n";
        return 0;
    }

    // get tag
    int iData[6];
    int numData = 6;
    if (OPS_GetIntInput(&numData, iData) < 0) {
        opserr << "WARNING invalid integer input for curved pipe "
                  "element\n";
        return 0;
    }

    // get center and radius
    double data[4];
    numData = 4;
    if (OPS_GetDoubleInput(&numData, data) < 0) {
        opserr << "WARNING invalid center or radius input for curved "
                  "pipe element\n";
        return 0;
    }

    // get data
    double T0 = 0.0, pressure = 0.0;
    int cMass = 0;
    int releasez = 0;
    int releasey = 0;
    numData = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char *theType = OPS_GetString();
        if (strcmp(theType, "-T0") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &T0) < 0) {
                    opserr << "WARNING: failed to read T0\n";
                    return 0;
                }
            }
        } else if (strcmp(theType, "-p") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numData, &pressure) < 0) {
                    opserr << "WARNING: failed to read internal "
                              "pressure\n";
                    return 0;
                }
            }
        } else if (strcmp(theType, "-cMass") == 0) {
            cMass = 1;
        } else if (strcmp(theType, "-releasez") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numData, &releasez) < 0) {
                    opserr << "WARNING: failed to get releasez";
                    return 0;
                }
            }
        } else if (strcmp(theType, "-releasey") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numData, &releasey) < 0) {
                    opserr << "WARNING: failed to get releasey";
                    return 0;
                }
            }
        }
    }

    auto *theSect = dynamic_cast<PipeSection *>(
        OPS_getSectionForceDeformation(iData[5]));
    if (theSect == 0) {
        opserr << "WARNING: section " << iData[5]
               << " is not found or not a curved pipe section\n";
        return 0;
    }

    auto *theMat = dynamic_cast<PipeMaterial *>(
        OPS_getUniaxialMaterial(iData[4]));
    if (theMat == 0) {
        opserr << "WARNING: uniaxialMaterial " << iData[4]
               << " is not found or not a curved pipe material\n";
        return 0;
    }

    auto *theTrans = OPS_getCrdTransf(iData[3]);
    if (theTrans == 0) {
        opserr << "WARNING: CrdTransf " << iData[3]
               << " is not found\n";
        return 0;
    }

    auto *ele = new CurvedPipe(iData[0], iData[1], iData[2],
                               *theTrans, *theMat, *theSect, T0,
                               pressure, cMass, releasez, releasey);

    return ele;
}

CurvedPipe::CurvedPipe() : Pipe() {}

CurvedPipe::CurvedPipe(int tag, int nd1, int nd2,
                       CrdTransf &theTransf, PipeMaterial &mat,
                       PipeSection &sect, double to, double pre,
                       int cm, int rz, int ry)
    : Pipe(tag, ELE_TAG_CurvedPipe) {
    if (Pipe::createPipe(nd1, nd2, theTransf, mat, sect, cm, rz, ry) <
        0) {
        opserr << "WARNING: failed to create curved pipe element\n";
        exit(-1);
    }
}