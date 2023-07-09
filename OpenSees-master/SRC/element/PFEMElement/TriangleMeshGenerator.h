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
                                                                        
// $Revision: 1.0 $
// $Date: 2015-2-17 12:01:59 $
                                                                        
// Written: Minjie Zhu
//
// Description: This file defines the class 'TriangleMeshGenerator', which
//              is a c++ wrapper of 'triangle' program.

#ifndef TriangleMeshGenerator_h
#define TriangleMeshGenerator_h

//#define REAL double
//#define VOID void
extern "C" {
#include <triangle.h>
}

#include <vector>

class TriangleMeshGenerator
{
public:
    TriangleMeshGenerator();
    ~TriangleMeshGenerator();

    // mesh
    int mesh(double size, bool pointOnBoundary=true);
    int remesh(double alpha);

    // inputs
    int addPoint(double x, double y);
    int addSegment(int p1, int p2, int mark);

    // outputs
    int getNumPoints() const;
    void getPoint(int i, double& x, double& y, int& mark);
    int getNumTriangles() const;
    void getTriangle(int i, int& p1, int& p2, int& p3);
    void getNeighbor(int i, int& t1, int& t2, int& t3);

    // clear
    void clear();
    
private:
    void reset();
    void initializeTri(triangulateio& in);
    void freeTri(triangulateio& in);
    void freeTriOut(triangulateio& in);
    
    triangulateio in, out, vout;

    std::vector<double> pointlist;
    std::vector<int> pointmarkerlist;
    std::vector<int> segmentlist;
    std::vector<int> segmentmarkerlist;
    std::vector<int> trianglelist;
    std::vector<int> neighborlist;
    int numberofcorners;
};


#endif
