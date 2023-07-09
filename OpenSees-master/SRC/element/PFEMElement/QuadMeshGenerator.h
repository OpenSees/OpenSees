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
// Description: A class for quad mesh generator
//


#ifndef QuadMeshGenerator_h
#define QuadMeshGenerator_h

#include <vector>
#include <Vector.h>
#include <ID.h>

class QuadMeshGenerator
{
public:
    QuadMeshGenerator();
    ~QuadMeshGenerator();

    // mesh
    int mesh(double size);

    // inputs
    int addPoint(const Vector& crds);
    int addLine(const ID& pts);

    // outputs
    int getNumPoints() const;
    void getPoint(int i, Vector& crds);
    int getNumQuads() const;
    void getQuad(int i, ID& pts);

    // clear()
    void clear();
    void clearOutput();

private:
    std::vector<Vector> point;
    std::vector<ID> line;

    std::vector<Vector> pointout;
    std::vector<ID> quadout;
};


#endif
