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
                                                                        
// $Revision $
// $Date $
                                                                        
// Written: Minjie Zhu
//
// Description: This file defines the class 'TetMeshGenerator', which
//              is a c++ wrapper of 'TetGen' program.

#ifndef TetMeshGenerator_h
#define TetMeshGenerator_h

#include <tetgen.h>
#include <vector>

class TetMeshGenerator
{
public:

    typedef std::vector<int> Polygon;
    typedef std::vector<Polygon> Facet;
    
public:
    TetMeshGenerator();
    ~TetMeshGenerator();

    // mesh
    int mesh(double vol, bool pointOnBoundary=true);
    int remesh(double alpha);

    // inputs
    int addPoint(double x, double y, double z, int mark);
    int addHole(double x, double y, double z);
    int addFacet(const Facet& facet, int mark);

    // outputs
    int getNumPoints() const;
    void getPoint(int i, double& x, double& y, double& z, int& mark);
    int getNumTets() const;
    void getTet(int i, int& p1, int& p2, int& p3, int&p4);
    void getNeighbor(int i, int& t1, int& t2, int& t3, int& t4);
    int getNumFaces() const;
    void getTriFace(int i, int& p1, int& p2, int& p3, int& mark);
    int getNumEdges() const;
    void getEdge(int i, int& p1, int& p2, int& mark);

    // clear
    void clear();
    
private:
    void reset();
    
    tetgenio in, out;

    std::vector<double> pointlist;
    std::vector<int> pointmarkerlist;
    std::vector<Facet> facetlist;
    std::vector<int> facetmarkerlist;
    std::vector<int> tetrahedronlist;
    std::vector<int> neighborlist;
    std::vector<double> holelist;
    std::vector<int> trifacelist;
    std::vector<int> trifacemarkerlist;
    std::vector<int> edgelist;
    std::vector<int> edgemarkerlist;
    int numberofcorners;
};


#endif
