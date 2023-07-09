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

#include "TriangleMeshGenerator.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <OPS_Globals.h>

TriangleMeshGenerator::TriangleMeshGenerator()
    :in(),out(),vout(),pointlist(),pointmarkerlist(),
     segmentlist(),segmentmarkerlist(),trianglelist(),
     neighborlist(),numberofcorners(0)
{
    initializeTri(in);
    initializeTri(out);
    initializeTri(vout);
}

TriangleMeshGenerator::~TriangleMeshGenerator()
{
    this->reset();
}

void
TriangleMeshGenerator::reset()
{
    freeTri(in);
    freeTri(vout);
    freeTriOut(out);
}

int
TriangleMeshGenerator::addPoint(double x, double y)
{
    pointlist.push_back(x);
    pointlist.push_back(y);
    return 0;
}

int
TriangleMeshGenerator::addSegment(int p1, int p2, int mark)
{
    segmentlist.push_back(p1);
    segmentlist.push_back(p2);
    segmentmarkerlist.push_back(mark);
    return 0;
}

int
TriangleMeshGenerator::getNumPoints() const
{
    return pointlist.size()/2;
}

int
TriangleMeshGenerator::getNumTriangles() const
{
    return trianglelist.size()/numberofcorners;
}

void
TriangleMeshGenerator::getPoint(int i, double& x, double& y, int& mark)
{
    if(i<0 || 2*i>=(int)pointlist.size()) return;
    x = pointlist[2*i];
    y = pointlist[2*i+1];
    if(i<0 || i>=(int)pointmarkerlist.size()) return;
    mark = pointmarkerlist[i];
}

void
TriangleMeshGenerator::getTriangle(int i, int& p1, int& p2, int& p3)
{
    if(i<0 || numberofcorners*i>=(int)trianglelist.size()) return;
    p1 = trianglelist[numberofcorners*i];
    p2 = trianglelist[numberofcorners*i+1];
    p3 = trianglelist[numberofcorners*i+2];
}

void
TriangleMeshGenerator::getNeighbor(int i, int& t1, int& t2, int& t3)
{
    int num = 3;
    if(i<0 || num*i>=(int)neighborlist.size()) return;
    t1 = neighborlist[num*i];
    t2 = neighborlist[num*i+1];
    t3 = neighborlist[num*i+2];
}

int
TriangleMeshGenerator::mesh(double size, bool pointOnBoundary)
{
    // reset
    this->reset();

    // get size
    in.numberofpoints = pointlist.size()/2;
    in.numberofsegments = segmentlist.size()/2;
        
    // quick return
    if(in.numberofsegments < 1) return 0;
    if(in.numberofpoints < 1) return 0;
    
    // set pointers
    in.pointlist = &pointlist[0];
    in.segmentlist = &segmentlist[0];
    in.segmentmarkerlist = &segmentmarkerlist[0];

    // meshing
    char s[128];
    if (pointOnBoundary) {
	sprintf(s, "DnQzqpa%.20f", size);
    } else {
	sprintf(s, "DnYYQzqpa%.20f", size);
    }
    triangulate(s, &in, &out, &vout);

    // reset pointers
    in.pointlist = NULL;
    in.segmentlist = NULL;
    in.segmentmarkerlist = NULL;

    // clear input data
    pointlist.clear();
    pointmarkerlist.clear();
    segmentlist.clear();
    segmentmarkerlist.clear();
    trianglelist.clear();
    neighborlist.clear();
    
    // copy output data
    numberofcorners = out.numberofcorners;
    pointlist.assign(out.pointlist,
		     out.pointlist+out.numberofpoints*2);
    pointmarkerlist.assign(out.pointmarkerlist,
		     out.pointmarkerlist+out.numberofpoints);
    trianglelist.assign(out.trianglelist,
			out.trianglelist+out.numberoftriangles*numberofcorners);
    neighborlist.assign(out.neighborlist,
			out.neighborlist+out.numberoftriangles*3);


    // clear memory
    this->reset();
    
    return 0;
}

int
TriangleMeshGenerator::remesh(double alpha)
{
    // reset
    this->reset();

    // get size
    in.numberofpoints = pointlist.size()/2;
        
    // quick return
    if(in.numberofpoints < 3) return 0;
    
    // set pointers
    in.pointlist = &pointlist[0];

    // meshing
    char s[] = "Qnzv";
    triangulate(s, &in, &out, &vout);

    // reset pointers
    in.pointlist = NULL;

    // clear input data
    pointmarkerlist.clear();
    segmentlist.clear();
    segmentmarkerlist.clear();
    trianglelist.clear();
    neighborlist.clear();

    neighborlist.assign(out.neighborlist,
			out.neighborlist+out.numberoftriangles*3);

    // radius and average size of triangles
    numberofcorners = out.numberofcorners;
    int numtri = out.numberoftriangles;
    std::vector<double> radius(numtri);
    std::vector<double> beta(numtri);
    double avesize = 0.0;
    for(int i=0; i<numtri; i++) {

	// triangle circumcenter
	double xc = vout.pointlist[2*i];
	double yc = vout.pointlist[2*i+1];

	// triangle points
	int pt[3];
	for(int j=0; j<3; j++) {
	    pt[j] = out.trianglelist[numberofcorners*i+j];
	}

	// nodal coordinates
	double x[3], y[3];
	for(int j=0; j<3; j++) {
	    x[j] = out.pointlist[2*pt[j]];
	    y[j] = out.pointlist[2*pt[j]+1];
	}

	// size of triangle
	double he = -1.0;
	for(int j=0; j<3; j++) {
	    for(int k=j+1; k<3; k++) {
		double h = (x[j]-x[k])*(x[j]-x[k])+(y[j]-y[k])*(y[j]-y[k]);
		if(h<he || he==-1.0) {
		    he = h;
		}
	    }
	}

	// take average
	double she = sqrt(he);
	avesize += she;

	// radius
	radius[i] = sqrt((xc-x[0])*(xc-x[0])+(yc-y[0])*(yc-y[0]));

	// beta, element distortion
	beta[i] = radius[i]/she*sqrt(3.0);
    }
    avesize /= numtri;

    // alpha shape test
    for(int i=0; i<numtri; i++) {
	
	// if pass the alpha test
    	//if(radius[i] / avesize <= alpha && beta[i] <= 50.*alpha) {
	if(radius[i] / avesize <= alpha || alpha < 0) {
	    for(int j=0; j<3; j++) {
		trianglelist.push_back(out.trianglelist[numberofcorners*i+j]);
	    }
    	}
	// if(beta[i] > 100.*alpha && radius[i] / avesize <= alpha) {
	//     opserr<<"removing distorted triangles with beta = "<<beta[i]<<"\n";
	// }
    }

    // clear memory
    this->reset();

    return 0;
}

void
TriangleMeshGenerator::clear()
{
    pointlist.clear();
    pointmarkerlist.clear();
    segmentlist.clear();
    segmentmarkerlist.clear();
    trianglelist.clear();
    neighborlist.clear();
    numberofcorners = 0;
    this->reset();
}
    

void
TriangleMeshGenerator::initializeTri(triangulateio& in)
{
    in.pointlist = NULL;
    in.pointattributelist = NULL;
    in.pointmarkerlist = NULL;
    in.numberofpoints = 0;
    in.numberofpointattributes = 0;

    in.trianglelist = NULL;
    in.triangleattributelist = NULL;
    in.trianglearealist = NULL;
    in.neighborlist = NULL;
    in.numberoftriangles = 0;
    in.numberofcorners = 0;
    in.numberoftriangleattributes = 0;

    in.segmentlist = NULL;
    in.segmentmarkerlist = NULL;
    in.numberofsegments = 0;

    in.holelist = NULL;
    in.numberofholes = 0;

    in.regionlist = NULL;
    in.numberofregions = 0;

    in.edgelist = NULL;
    in.edgemarkerlist = NULL;
    in.normlist = NULL;
    in.numberofedges = 0;
}

void
TriangleMeshGenerator::freeTri(triangulateio& in)
{
    if(in.pointlist != NULL) {
        free(in.pointlist);
    }
    if(in.pointattributelist != NULL) {
        free(in.pointattributelist);
    }
    if(in.pointmarkerlist != NULL) {
        free(in.pointmarkerlist);
    }

    if(in.trianglelist != NULL) {
        free(in.trianglelist);
    }
    if(in.triangleattributelist != NULL) {
        free(in.triangleattributelist);
    }
    if(in.trianglearealist != NULL) {
        free(in.trianglearealist);
    }
    if(in.neighborlist != NULL) {
        free(in.neighborlist);
    }

    if(in.segmentlist != NULL) {
        free(in.segmentlist);
    }
    if(in.segmentmarkerlist != NULL) {
        free(in.segmentmarkerlist);
    }

    if(in.holelist != NULL) {
        free(in.holelist);
    }

    if(in.regionlist != NULL) {
        free(in.regionlist);
    }

    if(in.edgelist != NULL) {
        free(in.edgelist);
    }
    if(in.edgemarkerlist != NULL) {
        free(in.edgemarkerlist);
    }
    if(in.normlist != NULL) {
        free(in.normlist);
    }

    initializeTri(in);
}

void
TriangleMeshGenerator::freeTriOut(triangulateio& in)
{
    if(in.pointlist != NULL) {
        free(in.pointlist);
    }
    if(in.pointattributelist != NULL) {
        free(in.pointattributelist);
    }
    if(in.pointmarkerlist != NULL) {
        free(in.pointmarkerlist);
    }

    if(in.trianglelist != NULL) {
        free(in.trianglelist);
    }
    if(in.triangleattributelist != NULL) {
        free(in.triangleattributelist);
    }
    if(in.trianglearealist != NULL) {
        free(in.trianglearealist);
    }
    if(in.neighborlist != NULL) {
        free(in.neighborlist);
    }

    if(in.segmentlist != NULL) {
        free(in.segmentlist);
    }
    if(in.segmentmarkerlist != NULL) {
        free(in.segmentmarkerlist);
    }

    if(in.edgelist != NULL) {
        free(in.edgelist);
    }
    if(in.edgemarkerlist != NULL) {
        free(in.edgemarkerlist);
    }
    if(in.normlist != NULL) {
        free(in.normlist);
    }

    initializeTri(in);
}


