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


#ifndef LineMeshGenerator_h
#define LineMeshGenerator_h

#include <vector>
#include <Vector.h>
#include <ID.h>

class LineMeshGenerator
{
public:
    LineMeshGenerator();
    ~LineMeshGenerator();

    // mesh
    int mesh(double size);

    // inputs
    int addPoint(const Vector& crds);
    int addLine(const ID& pts);

    // outputs
    int getNumPoints() const;
    void getPoint(int i, Vector& crds);
    int getNumLines() const;
    void getLine(int i, ID& pts);

    // clear()
    void clear();
    void clearOutput();

private:
    int meshLine(double size, const ID& pts);

private:
    std::vector<Vector> point;
    std::vector<ID> line;

    std::vector<Vector> pointout;
    std::vector<ID> lineout;
};


#endif
