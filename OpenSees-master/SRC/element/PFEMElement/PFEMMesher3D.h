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
// $Date: 2013-2-25 11:25:05 $
// $Source: SRC/element/PFEMElement/PFEMMesher3D.h,v $
                                                                        
// Written: Minjie Zhu
// Created: Feb. 25
//
// Description: The class PFEMMesher3D encapsulates all necessary
// routines for PFEM analysis for fluid-structure interaction
// The triangulation uses the 
// tetgen program written by Hang Si http://tetgen.berlios.de
//


//
// What: "@(#) PFEMMesher3D.h, revA"

#ifndef PFEMMesher3D_H
#define PFEMMesher3D_H

#include <tetgen.h>
#include <string>
#include <Vector.h>
#include <vector>
#include <map>
#include <set>

class ID;
class Domain;


class PFEMMesher3D 
{

public:
    // defs
    typedef std::vector<int> Polygon;
    typedef std::vector<Polygon> PolygonVec;
    typedef std::vector<int> ivector;
    typedef std::vector<double> dvector;

    class Point {
    private:
        double x1,x2,x3;
    public:
        Point():x1(0),x2(0),x3(0){}
        Point(const double& xx, const double& yy, const double& zz)
            :x1(xx),x2(yy),x3(zz) {}
        void coord(const double& xx, const double& yy, const double& zz){
            x1=xx;x2=yy;x3=zz;
        }
        const double& x()const {return x1;}
        const double& y()const {return x2;}
        const double& z()const {return x3;}
    };
    typedef std::vector<Point> PointVec;

    class Facet {
    private:
        PolygonVec polys;
        PointVec hls;
    public:
        Facet():polys(),hls(){}
        ~Facet(){}
        void addpolygon(const Polygon& polygon){
            polys.push_back(polygon);
        }
        void addhole(const double& x, const double& y, const double& z){
            hls.push_back(Point(x,y,z));
        }
        const PolygonVec& polygons()const {return polys;}
        const PointVec& holes()const {return hls;}
    };
    typedef std::vector<Facet> FacetVec;

public:

    // constructor
    PFEMMesher3D();

    // destructor
    ~PFEMMesher3D();

    // discretize domain
    int discretize(int startnode, const PointVec& points, const FacetVec facets,
                   const PointVec& holes, double maxvol, int ndf, 
                   const ivector& fix, const dvector& mass,
                   Domain* theDomain, int& endnode); // PLC
    int discretize(int startnode, const Point& pt, const dvector& hs, 
                   const ivector& ns, int ndf, const ivector& fix,
                   const dvector& mass, const dvector& vel,
                   Domain* theDomain, int& endnode); // cube

    // triangulation
    int doTriangulation(const ivector& nodes, double alpha, double volthresh,
                        const ivector& addnodes, Domain* theDomain, 
                        ivector& eles);
    int doTriangulation(int startele, const ivector& nodes, double alpha, 
                        double valthresh, const ivector& addnodes, Domain* theDomain,
                        double rho, double mu, double b1, double b2, double b3);

    // save
    int save(const char* name, const ID& snodes, int step, Domain* theDomain);

    // set boundary of fluids
    void setBoundary(double x1, double y1, double z1, double x2, double y2, double z2);
    void removeOutBoundNodes(const ID& nodes, Domain* theDomain);

private:

    Vector bound;
    double avesize;
};




#endif
