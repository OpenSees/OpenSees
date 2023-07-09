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
// Description: This file defines the class 'TetMeshGenerator', which
//              is a c++ wrapper of 'TetGen' program.

#include "TetMeshGenerator.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <OPS_Globals.h>

TetMeshGenerator::TetMeshGenerator()
    :in(),out(), pointlist(), pointmarkerlist(),
     facetlist(), facetmarkerlist(),tetrahedronlist(),
     neighborlist(), holelist(), trifacelist(),
     trifacemarkerlist(), edgelist(), edgemarkerlist(),
     numberofcorners(4)
{
}

TetMeshGenerator::~TetMeshGenerator()
{
    clear();
    reset();
}

void
TetMeshGenerator::clear()
{
    pointlist.clear();
    pointmarkerlist.clear();
    facetlist.clear();
    facetmarkerlist.clear();
    tetrahedronlist.clear();
    neighborlist.clear();
    holelist.clear();
    trifacelist.clear();
    trifacemarkerlist.clear();
    edgelist.clear();
    edgemarkerlist.clear();
    numberofcorners = 4;
}

void
TetMeshGenerator::reset()
{
    in.deinitialize();
    in.initialize();
    out.deinitialize();
    out.initialize();
}

int
TetMeshGenerator::addPoint(double x, double y, double z, int mark)
{
    pointlist.push_back(x);
    pointlist.push_back(y);
    pointlist.push_back(z);
    pointmarkerlist.push_back(mark);
    return 0;
}

int
TetMeshGenerator::addHole(double x, double y, double z)
{
    holelist.push_back(x);
    holelist.push_back(y);
    holelist.push_back(z);
    return 0;
}

int
TetMeshGenerator::addFacet(const Facet& facet, int mark)
{
    facetlist.push_back(facet);
    facetmarkerlist.push_back(mark);
    return 0;
}

int
TetMeshGenerator::getNumPoints() const
{
    return (int)pointlist.size()/3;
}

int
TetMeshGenerator::getNumTets() const
{
    return (int)tetrahedronlist.size()/numberofcorners;
}

void
TetMeshGenerator::getPoint(int i, double& x, double& y, double&z, int& mark)
{
    if(i<0 || 3*i>=(int)pointlist.size()) return;
    x = pointlist[3*i];
    y = pointlist[3*i+1];
    z = pointlist[3*i+2];
    if(i<0 || i>=(int)pointmarkerlist.size()) return;
    mark = pointmarkerlist[i];
}

void
TetMeshGenerator::getTet(int i, int& p1, int& p2, int& p3, int&p4)
{
    if(i<0 || numberofcorners*i>=(int)tetrahedronlist.size()) return;
    p1 = tetrahedronlist[numberofcorners*i];
    p2 = tetrahedronlist[numberofcorners*i+1];
    p3 = tetrahedronlist[numberofcorners*i+2];
    p4 = tetrahedronlist[numberofcorners*i+3];
}

void
TetMeshGenerator::getNeighbor(int i, int& t1, int& t2, int& t3, int& t4)
{
    int num = 4;
    if(i<0 || num*i>=(int)neighborlist.size()) return;
    t1 = neighborlist[num*i];
    t2 = neighborlist[num*i+1];
    t3 = neighborlist[num*i+2];
    t4 = neighborlist[num*i+3];
}

int
TetMeshGenerator::getNumFaces() const
{
    return (int)trifacelist.size()/3;
}

void
TetMeshGenerator::getTriFace(int i, int& p1, int& p2, int& p3, int& mark)
{
    if(i<0 || 3*i>=(int)trifacelist.size()) return;
    p1 = trifacelist[3*i];
    p2 = trifacelist[3*i+1];
    p3 = trifacelist[3*i+2];
    if (i<0 || i>(int)trifacemarkerlist.size()) {
	return;
    }
    mark = trifacemarkerlist[i];
}

int
TetMeshGenerator::getNumEdges() const
{
    return (int)trifacelist.size()/3;
}

void
TetMeshGenerator::getEdge(int i, int& p1, int& p2, int& mark)
{
    if(i<0 || 2*i>=(int)edgelist.size()) return;
    p1 = edgelist[2*i];
    p2 = edgelist[2*i+1];
    if (i<0 || i>(int)edgemarkerlist.size()) {
	return;
    }
    mark = edgemarkerlist[i];
}

int
TetMeshGenerator::mesh(double vol, bool pointOnBoundary)
{
    // reset
    reset();

    // get size
    in.numberofpoints = (int)pointlist.size()/3;
    in.numberoffacets = (int)facetlist.size();
        
    // quick return
    if(in.numberoffacets < 1) return 0;
    if(in.numberofpoints < 1) return 0;
    
    // set pointers
    in.pointlist = &pointlist[0];
    in.pointmarkerlist = &pointmarkerlist[0];
    in.facetmarkerlist = &facetmarkerlist[0];
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    for (int i=0; i<in.numberoffacets; ++i) {
	tetgenio::facet& f = in.facetlist[i];
	f.numberofpolygons = (int)facetlist[i].size();
	f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
	f.numberofholes = 0;
	f.holelist = NULL;
	for (int j=0; j<f.numberofpolygons; ++j) {
	    tetgenio::polygon& p = f.polygonlist[j];
	    p.numberofvertices = (int)facetlist[i][j].size();
	    p.vertexlist = new int[p.numberofvertices];
	    for (int k=0; k<p.numberofvertices; ++k) {
		p.vertexlist[k] = facetlist[i][j][k];
	    }
	}
    }
    if (holelist.empty() == false) {
	in.holelist = &holelist[0];
	in.numberofholes = (int)holelist.size();
    }

    // meshing
    char s[128];
    if (pointOnBoundary) {
	sprintf(s, "nQfezqpa%.20f", vol);
    } else {
	sprintf(s, "nYQfezqpa%.20f", vol);
    }
    tetrahedralize(s, &in, &out);

    // reset pointers
    in.pointlist = NULL;
    in.pointmarkerlist = NULL;
    in.facetmarkerlist = NULL;
    in.holelist = NULL;
    in.numberofpoints = 0;
    in.numberofholes = 0;

    // copy output data
    numberofcorners = out.numberofcorners;
    pointlist.assign(out.pointlist,
		     out.pointlist+out.numberofpoints*3);
    pointmarkerlist.assign(out.pointmarkerlist,
		     out.pointmarkerlist+out.numberofpoints);
    tetrahedronlist.assign(out.tetrahedronlist,
			   out.tetrahedronlist+out.numberoftetrahedra*numberofcorners);
    neighborlist.assign(out.neighborlist,
			out.neighborlist+out.numberoftetrahedra*4);
    trifacelist.assign(out.trifacelist,
		       out.trifacelist+out.numberoftrifaces*3);
    trifacemarkerlist.assign(out.trifacemarkerlist,
			     out.trifacemarkerlist+out.numberoftrifaces);
    edgelist.assign(out.edgelist,
		    out.edgelist+out.numberofedges*2);
    edgemarkerlist.assign(out.edgemarkerlist,
			  out.edgemarkerlist+out.numberofedges);


    // clear memory
    reset();
    
    return 0;
}

int
TetMeshGenerator::remesh(double alpha)
{
    // reset
    reset();

    // get size
    in.numberofpoints = (int)pointlist.size()/3;
        
    // quick return
    if(in.numberofpoints < 4) return 0;
    
    // set pointers
    in.pointlist = &pointlist[0];

    // meshing
    char s[] = "Qnzfev";
    tetrahedralize(s, &in, &out);

    // reset pointers
    in.pointlist = NULL;
    in.numberofpoints = 0;

    // copy output data
    numberofcorners = out.numberofcorners;
    neighborlist.assign(out.neighborlist,
			out.neighborlist+out.numberoftetrahedra*4);
    trifacelist.assign(out.trifacelist,
		       out.trifacelist+out.numberoftrifaces*3);
    trifacemarkerlist.assign(out.trifacemarkerlist,
			     out.trifacemarkerlist+out.numberoftrifaces);
    edgelist.assign(out.edgelist,
		    out.edgelist+out.numberofedges*2);
    edgemarkerlist.assign(out.edgemarkerlist,
			  out.edgemarkerlist+out.numberofedges);

    // radius and average size of triangles
    int numtet = out.numberoftetrahedra;
    std::vector<double> radius(numtet);
    double avesize = 0.0;
    for(int i=0; i<numtet; i++) {

	// triangle circumcenter
	double xc = out.vpointlist[3*i];
	double yc = out.vpointlist[3*i+1];
	double zc = out.vpointlist[3*i+2];

	// tetrahedra points
	int pt[4];
	for(int j=0; j<4; j++) {
	    pt[j] = out.tetrahedronlist[numberofcorners*i+j];
	}

	// nodal coordinates
	double x[4], y[4], z[4];
	for(int j=0; j<4; j++) {
	    x[j] = out.pointlist[3*pt[j]];
	    y[j] = out.pointlist[3*pt[j]+1];
	    z[j] = out.pointlist[3*pt[j]+2];
	}

	// size of tetrahedra
	double he = -1.0;
	for(int j=0; j<4; j++) {
	    for(int k=j+1; k<4; k++) {
		double h = (x[j]-x[k])*(x[j]-x[k])+(y[j]-y[k])*(y[j]-y[k])+(z[j]-z[k])*(z[j]-z[k]);
		if(h<he || he==-1.0) {
		    he = h;
		}
	    }
	}

	// take average
	double she = sqrt(he);
	avesize += she;

	// radius
	radius[i] = sqrt((xc-x[0])*(xc-x[0])+(yc-y[0])*(yc-y[0])+(zc-z[0])*(zc-z[0]));
    }
    avesize /= numtet;

    // alpha shape test
    for(int i=0; i<numtet; i++) {
	
	// if pass the alpha test
	if(radius[i] / avesize <= alpha || alpha < 0) {
	    for(int j=0; j<numberofcorners; j++) {
		tetrahedronlist.push_back(out.tetrahedronlist[numberofcorners*i+j]);
	    }
    	}
    }

    // clear memory
    reset();

    return 0;
}


