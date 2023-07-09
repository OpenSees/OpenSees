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

//
// Written: Minjie Zhu
//
// Description: A class for line mesh generator
//

#include "LineMeshGenerator.h"

#include <cmath>

LineMeshGenerator::LineMeshGenerator()
    : point(), line(), pointout(), lineout() {}

LineMeshGenerator::~LineMeshGenerator() {}

// mesh
int LineMeshGenerator::mesh(double size) {
    if (size <= 0) return -1;

    // copy point to pointout
    clearOutput();
    pointout = point;

    // mesh each line
    for (int i = 0; i < (int)line.size(); ++i) {
        if (meshLine(size, line[i]) < 0) {
            opserr << "WARNING: failed to mesh line\n";
            return -1;
        }
    }

    return 0;
}

// inputs
int LineMeshGenerator::addPoint(const Vector& crds) {
    point.push_back(crds);
    return 0;
}

int LineMeshGenerator::addLine(const ID& pts) {
    line.push_back(pts);
    return 0;
}

// outputs
int LineMeshGenerator::getNumPoints() const {
    return (int)pointout.size();
}

void LineMeshGenerator::getPoint(int i, Vector& crds) {
    if (i >= 0 && i < getNumPoints()) {
        crds = pointout[i];
    }
}

int LineMeshGenerator::getNumLines() const { return (int)lineout.size(); }

void LineMeshGenerator::getLine(int i, ID& pts) {
    if (i >= 0 && i < getNumLines()) {
        pts = lineout[i];
    }
}

int LineMeshGenerator::meshLine(double size, const ID& pts) {
    // check inputs
    if (pts.Size() != 2) return -1;
    int numpts = getNumPoints();
    if (pts(0) < 0 || pts(0) >= numpts) return -1;
    if (pts(1) < 0 || pts(1) >= numpts) return -1;

    // point coordinates
    const Vector& crds1 = point[pts(0)];
    const Vector& crds2 = point[pts(1)];

    if (crds1.Size() != crds2.Size()) {
        opserr << "WARNING: crds of points not compatible\n";
        return -1;
    }

    // increments of lines
    Vector incr = crds2 - crds1;
    // int nele = ceil(incr.Norm()/size);
    int nele = (int)floor(incr.Norm() / size + 0.5);
    if (nele < 1) {
        return 0;
    } else if (nele == 1) {
        lineout.push_back(pts);
        return 0;
    }
    incr /= nele;

    // first line
    Vector crds = crds1 + incr;
    ID lpts(2);
    lpts(0) = pts(0);
    lpts(1) = (int)pointout.size();
    pointout.push_back(crds);
    lineout.push_back(lpts);

    // create inner points and lines
    for (int i = 1; i < nele - 1; ++i) {
        // next line
        lpts(0) = lpts(1);
        lpts(1) = (int)pointout.size();

        crds += incr;
        pointout.push_back(crds);
        lineout.push_back(lpts);
    }

    // last line
    lpts(0) = lpts(1);
    lpts(1) = pts(1);
    lineout.push_back(lpts);

    return 0;
}

void LineMeshGenerator::clear() {
    point.clear();
    line.clear();
    clearOutput();
}

void LineMeshGenerator::clearOutput() {
    pointout.clear();
    lineout.clear();
}
