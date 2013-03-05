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
// $Date: 2012-10-31 9:56:05 $
// $Source: SRC/element/PFEMElement/PFEMMesher2D.cpp,v $
                                                                        
// Written: Minjie Zhu
// Created: Oct. 31
//
// Description: The class PFEMMesher2D encapsulates all necessary
// routines for PFEM analysis for fluid-structure interaction
// The triangulation uses the 
// Triangle program written by Jonathan Richard Shewchuk
//

//
// What: "@(#) PFEMMesher2D.cpp, revA"

#include "PFEMMesher2D.h"
#include <Matrix.h>
#include <ID.h>
#include <Domain.h>
#include <Node.h>
#include <Element.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <Pressure_Constraint.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <Pressure_ConstraintIter.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <PFEMElement2D.h>
#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include <Timer.h>

PFEMMesher2D::PFEMMesher2D()
    :PI(3.1415926535897932384626433), bound(4)
{
}

PFEMMesher2D::~PFEMMesher2D()
{
}

int 
PFEMMesher2D::discretize(int startnodetag, const Vector& points, const ID& segments, 
                         const Vector& holes, double maxarea, int ndf, const ID& fix, 
                         const Vector& mass, Domain* theDomain)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }

    // structs
    triangulateio in, out, vout;
    initializeTri(in);
    initializeTri(out);
    initializeTri(vout);

    // number of segments
    in.numberofsegments = segments.Size()/2;

    // no segments
    if(in.numberofsegments < 1) {
        return 0;
    }

    // number of points
    in.numberofpoints = points.Size()/2;

    // no points
    if(in.numberofpoints < 2) {
        return 0;
    }

    // number of holes
    in.numberofholes = holes.Size()/2;
    
    // get space for segmentlist 
    in.segmentlist = (int*) calloc(segments.Size(), sizeof(int));
    if(in.segmentlist == NULL) {
        opserr<<"WARNING: no enough memory -- PFEMMesher2D::discretize\n";
        return -1;
    }

    // copy segments
    int maxpt=0;
    for(int i=0; i<segments.Size(); i++) {
        if(segments(i) > maxpt) maxpt = segments(i);
        in.segmentlist[i] = segments(i);
    }
    if(maxpt >= in.numberofpoints) {
        opserr<<"WARNING: no point "<<maxpt<<" in the point list -- ";
        opserr<<"PFEMMesher2D::discretize\n";
        freeTri(in);
        return -1;
    }

    // get space for pointlist
    in.pointlist = (double*) calloc(points.Size(), sizeof(double));
    if(in.pointlist == NULL) {
        opserr<<"WARNING: no enough memory -- PFEMMesher2D::discretize\n";
        return -1;
    }

    // copy points
    for(int i=0; i<points.Size(); i++) {
        in.pointlist[i] = points(i);
    }

    // get space for holelist
    if(in.numberofholes > 0) {
        in.holelist = (double*) calloc(holes.Size(), sizeof(double));
        if(in.holelist == NULL) {
            opserr<<"WARNING: no enough memory -- PFEMMesher2D::discretize\n";
            return -1;
        }
        
        // copy points
        for(int i=0; i<holes.Size(); i++) {
            in.holelist[i] = holes(i);
        }
    }

    // conforming Delaunay triangulation
    //Timer timer;
    //timer.start();
    char s[100];
    sprintf(s,"Qzqpa%.60f",maxarea);
    //opserr<<"discretize : ";
    triangulate(s, &in, &out, &vout);
    //timer.pause();
    //opserr<<timer;
    freeTri(in);

    // read outputs
    int newtag = startnodetag-1;
    for(int i=0; i<out.numberofpoints; i++) {
        const double& x = out.pointlist[2*i];
        const double& y = out.pointlist[2*i+1];

        // create nodes
        Node* theNode = new Node(++newtag, ndf, x, y);
        if (theNode == 0) {
            opserr << "WARNING ran out of memory creating node\n";
            opserr << "node: " << newtag << "\n";
            opserr << "PFEMMesher2D::discretize\n";
            return -1;
        }
        if(theDomain->addNode(theNode) == false) {
            opserr << "WARNING failed to add node to the domain\n";
            opserr << "node: " << newtag << "\n";
            opserr << "PFEMMesher2D::discretize\n";
            delete theNode; // otherwise memory leak
            return -1;
        }
        
        // fix
        int numfix = fix.Size();
        if(numfix > ndf) numfix = ndf;
        for(int i=0; i<numfix; i++) {
            if(fix(i) != 0) {
                SP_Constraint* theSP = new SP_Constraint(newtag, i, 0.0, true);
                if(theSP == 0) {
                    opserr << "WARNING ran out of memory creating SP_Constraint\n";
                    opserr << "node: " << newtag << "\n";
                    opserr << "PFEMMesher2D::discretize\n";
                    return -1;
                }
                if(theDomain->addSP_Constraint(theSP) == false) {
                    opserr << "WARNING failed to add SP_Constraint to the domain\n";
                    opserr << "node: " << newtag << "\n";
                    opserr << "PFEMMesher2D::discretize\n";
                    delete theSP; // otherwise memory leak
                    return -1;
                }
            }
        }

        // mass
        int nummass = mass.Size();
        if(nummass > 0) {
            if(nummass > ndf) nummass = ndf;
            Matrix theMass(ndf, ndf);
            for(int i=0; i<nummass; i++) {
                theMass(i,i) = mass(i);
            }
            if(theNode->setMass(theMass) < 0) {
                opserr<<"WARNING: failed to set mass of node "<<newtag;
                opserr<<" -- PFEMMesher2D::discretize\n";
                return -1;
            }
        }
    }
    

    // free out and vout
    freeTriOut(out);
    freeTri(vout);


    return newtag;
}

int 
PFEMMesher2D::discretize(int startnodetag, double x1, double y1, double hx, double hy, double angle,
                         int nx, int ny, int ndf, const ID& fix, const Vector& mass, Domain* theDomain)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }

    if(nx<0 || ny<0) return 0;
    double sy = sin((angle+90)*PI/180.0);
    double cy = cos((angle+90)*PI/180.0);
    
    // foreach line
    double lhx = cy*hy;
    double lhy = sy*hy;

    int newtag = startnodetag-1;
    for(int i=0; i<=ny; i++) {
        double lx = i*lhx+x1;
        double ly = i*lhy+y1;
        newtag = discretize(++newtag, lx, ly, hx, angle, nx, ndf, fix, mass, theDomain);
    }


    return newtag;
}

int 
PFEMMesher2D::discretize(int startnodetag, double x1, double y1, double h, double angle, 
                         int num, int ndf, const ID& fix, const Vector& mass, Domain* theDomain)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }

    if(num < 0) return 0;

    double s = sin(angle*PI/180.0);
    double c = cos(angle*PI/180.0);

    double hx = c*h;
    double hy = s*h;

    int newtag = startnodetag-1;
    for(int i=0; i<=num; i++) {
        double x = i*hx+x1;
        double y = i*hy+y1;

        // define node
        Node* theNode = new Node(++newtag, ndf, x, y);
        if (theNode == 0) {
            opserr << "WARNING ran out of memory creating node\n";
            opserr << "node: " << newtag << "\n";
            opserr << "PFEMMesher2D::discretize\n";
            return -1;
        }
        if(theDomain->addNode(theNode) == false) {
            opserr << "WARNING failed to add node to the domain\n";
            opserr << "node: " << newtag << "\n";
            opserr << "PFEMMesher2D::discretize\n";
            delete theNode; // otherwise memory leak
            return -1;
        }

        // fix
        int numfix = fix.Size();
        if(numfix > ndf) numfix = ndf;
        for(int i=0; i<numfix; i++) {
            if(fix(i) != 0) {
                SP_Constraint* theSP = new SP_Constraint(newtag, i, 0.0, true);
                if(theSP == 0) {
                    opserr << "WARNING ran out of memory creating SP_Constraint\n";
                    opserr << "node: " << newtag << "\n";
                    opserr << "PFEMMesher2D::discretize\n";
                    return -1;
                }
                if(theDomain->addSP_Constraint(theSP) == false) {
                    opserr << "WARNING failed to add SP_Constraint to the domain\n";
                    opserr << "node: " << newtag << "\n";
                    opserr << "PFEMMesher2D::discretize\n";
                    delete theSP; // otherwise memory leak
                    return -1;
                }
            }
        }

        // mass
        int nummass = mass.Size();
        if(nummass > 0) {
            if(nummass > ndf) nummass = ndf;
            Matrix theMass(ndf, ndf);
            for(int i=0; i<nummass; i++) {
                theMass(i,i) = mass(i);
            }
            if(theNode->setMass(theMass) < 0) {
                opserr<<"WARNING: failed to set mass of node "<<newtag;
                opserr<<" -- PFEMMesher2D::discretize\n";
                return -1;
            }
        }
    }

    return newtag;
}

int
PFEMMesher2D::doTriangulation(int starteletag, int startnodetag, int endnodetag, double alpha, 
                              int startanodetag, int endanodetag, Domain* theDomain,
                              double rho, double mu, double b1, double b2)

{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }
    //Timer timer;
    // do triangulation
    ID eles;
    int res = doTriangulation(startnodetag, endnodetag, alpha,
                              startanodetag, endanodetag, 
                              theDomain, eles);
    if(res < 0) {
        return res;
    }

    // add PFEM elements
    //timer.start();
    int numeles = eles.Size()/3;
    if(numeles == 0) return 0;
    int etag = starteletag-1;
    for(int i=0; i<numeles; i++) {
        PFEMElement2D* theEle = new PFEMElement2D(++etag, eles(3*i), eles(3*i+1), eles(3*i+2),
                                                  rho, mu, b1, b2);

        if(theEle == 0) {
            opserr<<"WARNING: no enough memory -- ";
            opserr<<" -- delaunay2D::doTriangulation\n";
            return -1;
        }
        if(theDomain->addElement(theEle) == false) {
            opserr<<"WARNING: failed to add element to domain -- ";
            opserr<<" -- delaunay2D::doTriangulation\n";
            delete theEle;
            return -1;
        }
    }
    //timer.pause();
    //opserr<<"create PFEM elements :"<<timer;
    return etag;
    
}

int
PFEMMesher2D::doTriangulation(int startnodetag, int endnodetag, double alpha, 
                              int startanodetag, int endanodetag , Domain* theDomain, ID& eles)
{
    //Timer theTimer;
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }

    // structs
    triangulateio in, out, vout;
    initializeTri(in);
    initializeTri(out);
    initializeTri(vout);

    // number of nodes
    int numnodes = endnodetag-startnodetag+1;
    int numadd = endanodetag-startanodetag+1;
    if(numnodes < 0) numnodes = 0;
    if(numadd < 0) numadd = 0;
    in.numberofpoints = numnodes+numadd;
    if(in.numberofpoints < 3) {
        return 0;
    }

    // get space for pointlist
    in.pointlist = (double*) calloc(in.numberofpoints*2, sizeof(double));
    if(in.pointlist == NULL) {
        opserr<<"WARNING: no enough memory -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }

    // copy nodes to pointlist, use current coordinates
    int loc = 0, numpoints = 0;
    ID p2nd(0, in.numberofpoints);
    for(int i=0; i<in.numberofpoints; i++) {

        // pointer to node
        Node* node = 0;
        int ndtag = 0;
        if(i<numnodes) {
            ndtag = startnodetag+i;
        } else {
            ndtag = startanodetag+i-numnodes;
        }
        node = theDomain->getNode(ndtag);
        if(node == 0) continue;

        // get coordinates
        const Vector& coord = node->getCrds();
        if(coord.Size() < 2) {
            //opserr<<"node "<<ndtag<<" has ndm = "<<coord.Size()<<"\n";
            opserr<<"WARNING: 2d is required -- PFEMMesher2D::doTriangulation\n";
            return -1;
        }
        const Vector& disp = node->getTrialDisp();
        if(disp.Size() < 2) {
            opserr<<"WARNING: 2 ndf is required -- PFEMMesher2D::doTriangulation\n";
            return -1;
        }

        // set pointlist
        double x = coord(0);        // initial
        double y = coord(1);
        x += disp(0);               // current
        y += disp(1);

        in.pointlist[loc++] = x;
        in.pointlist[loc++] = y;

        p2nd[numpoints++] = ndtag;

    }
    in.numberofpoints = numpoints;

    // Delaunay Triangulation
    char s[] = "QIzv";
    //theTimer.start();
    triangulate(s, &in, &out, &vout);
    //theTimer.pause();
    //opserr<<"triangulation";
    //opserr<<theTimer;
    // free in
    freeTri(in);

    // no outputs
    if(out.numberoftriangles < 1) {
        return 0;
    }

    // do alpha shape test
    if(alpha > 0) {

        // radius and average size of triangles
        double avesize = 0.0;
        Vector radius(out.numberoftriangles);
        for(int i=0; i<out.numberoftriangles; i++) {

            // circumcenter of traingle
            double& xc = vout.pointlist[2*i];
            double& yc = vout.pointlist[2*i+1];

            // triangle points
            int pt[3];
            for(int j=0; j<3; j++) {
                pt[j] = out.trianglelist[out.numberofcorners*i+j];
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
            avesize += sqrt(he);
            // double h = sqrt((x[1]-x[0])*(x[1]-x[0])+(y[1]-y[0])*(y[1]-y[0]));
            // h += sqrt((x[1]-x[2])*(x[1]-x[2])+(y[1]-y[2])*(y[1]-y[2]));
            // h += sqrt((x[0]-x[2])*(x[0]-x[2])+(y[0]-y[2])*(y[0]-y[2]));
            // h /= 3.0;
            // avesize += h;
            // double area = fabs((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0])) / 2.0;
            // avesize += sqrt(area);

            // radius
            radius(i) = sqrt((xc-x[0])*(xc-x[0])+(yc-y[0])*(yc-y[0]));
        }
        avesize /= out.numberoftriangles;

        // alpha test
        int num = 0;
        eles.resize(out.numberoftriangles*3);
        for(int i=0; i<out.numberoftriangles; i++) {
            if(radius(i) / avesize <= alpha) {
                // pass the test

                // check if all nodes are additional nodes
                bool add = true;
                for(int j=0; j<3; j++) {
                    int tag = p2nd(out.trianglelist[out.numberofcorners*i+j]);
                    if(tag>=startnodetag && tag<=endnodetag) {
                        add = false;
                        break;
                    }
                }
                    
                // add ele
                if(!add) {
                    for(int j=0; j<3; j++) {
                        int tag = p2nd(out.trianglelist[out.numberofcorners*i+j]);
                        eles(num++) = tag;
                    }
                }

            }
        }
        if(num == 0) {
            eles = ID();
        } else {
            eles.resize(num);
        }
        
    } else if(alpha < 0) {

        // eles
        eles.resize(out.numberoftriangles*3);

        // copy triangles
        for(int i=0; i<out.numberoftriangles; i++) {
            for(int j=0; j<3; j++) {
                int tag = p2nd(out.trianglelist[out.numberofcorners*i+j]);
                eles(3*i+j) = tag;
            }
        }

    }
    
    // free vout 
    freeTri(vout);

    // free out
    freeTriOut(out);

    return 0;
}

int 
PFEMMesher2D::save(const char* filename, int startsnode, int endsnode, Domain* theDomain)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::save\n";
        return -1;
    }

    // write nodes
    std::ofstream file2(std::string(filename).append(".node").c_str());
    file2<<theDomain->getCurrentTime()<<" "<<0<<" "<<0<<" ";
    file2<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<"\n";

    std::map<int,int> nodeIndex;
    int index = 1;

    // pressure constraints
    Pressure_ConstraintIter& thePCs = theDomain->getPCs();
    Pressure_Constraint* thePC = 0;
    while((thePC = thePCs()) != 0) {
        int ntag = thePC->getTag();
        Node* node = theDomain->getNode(ntag);
        int ptag = thePC->getPressureNode();
        Node* pnode = theDomain->getNode(ptag);
        if(node!=0 || pnode!=0) {
            nodeIndex[ntag] = index++;
            const Vector& coord = node->getCrds();
            const Vector& disp = node->getDisp();
            const Vector& vel = node->getVel();
            const Vector& accel = node->getAccel();
            const Vector& pvel = pnode->getVel();
            file2<<coord(0)+disp(0)<<" "<<coord(1)+disp(1)<<" "<<1<<" ";
            file2<<vel(0)<<" "<<vel(1)<<" "<<accel(0)<<" "<<accel(1);
            file2<<" "<<pvel(2)<<"\n";
        }
    }

    // nodes
    for(int tag=startsnode; tag<=endsnode; tag++) {
        Node* theNode = theDomain->getNode(tag);
        if(theNode != 0) {
            nodeIndex[tag] = index++;
            const Vector& coord = theNode->getCrds();
            const Vector& disp = theNode->getDisp();
            const Vector& vel = theNode->getVel();
            const Vector& accel = theNode->getAccel();
            file2<<coord(0)+disp(0)<<" "<<coord(1)+disp(1)<<" "<<2<<" ";
            file2<<vel(0)<<" "<<vel(1)<<" "<<accel(0)<<" "<<accel(1);
            file2<<" "<<0<<"\n";

        }
    }
    file2.close();

    // write edges
    std::ofstream file1(std::string(filename).append(".ele").c_str());
    ElementIter& theEles = theDomain->getElements();
    Element* theEle = 0;
    while((theEle = theEles()) != 0) {
        const ID& allntags = theEle->getExternalNodes();
        ID ntags(0,allntags.Size());
        for(int i=0; i<allntags.Size(); i++) {
            if(nodeIndex.find(allntags(i)) != nodeIndex.end()) {
                ntags[ntags.Size()] = allntags(i);
            }
        }
        if(ntags.Size()==1) continue;
        for(int i=0; i<ntags.Size(); i++) {
            file1<<nodeIndex[ntags(i)]<<" ";
        }
        file1<<"\n";
    }
    file1.close();

    // std::ofstream file1(std::string(filename).append(".edge").c_str());
    // typedef std::set< std::pair<int,int> > EdgeSet;
    // typedef EdgeSet::iterator EdgeSetIter;
    // typedef std::pair<EdgeSetIter, bool> EdgeSetRet;
    // EdgeSet edges;
    // EdgeSetRet ret;
    // ElementIter& theEles = theDomain->getElements();
    // Element* theEle = 0;
    // while((theEle = theEles()) != 0) {
    //     const ID& allntags = theEle->getExternalNodes();
    //     ID ntags(0,allntags.Size());
    //     for(int i=0; i<allntags.Size(); i++) {
    //         if(nodeIndex.find(allntags(i)) != nodeIndex.end()) {
    //             ntags[ntags.Size()] = allntags(i);
    //         }
    //     }
    //     if(ntags.Size()==1) continue;
    //     if(ntags.Size()>2) {
    //         ntags[ntags.Size()] = ntags(0);
    //     }
    //     for(int i=0; i<ntags.Size()-1; i++) {
    //         ID nd(0,2);
    //         for(int j=0; j<2; j++) {
    //             nd.insert(ntags(i+j));
    //         }
    //         ret = edges.insert(std::make_pair(nd(0),nd(1)));
    //         if(ret.second == true) {
    //             file1<<nodeIndex[nd(0)]<<" "<<nodeIndex[nd(1)]<<"\n";
    //         }
    //     }
    // }
    // file1.close();

    return 0;
}

void
PFEMMesher2D::setBoundary(double x1, double y1, double x2, double y2)
{
    bound(0) = x1;
    bound(1) = y1;
    bound(2) = x2;
    bound(3) = y2;
}

void 
PFEMMesher2D::removeOutBoundNodes(int startnode, int endnode, Domain* theDomain)
{
    // check nodes
    int numnodes = endnode-startnode+1;
    ID removelist(0, numnodes);
    for(int tag=startnode; tag<=endnode; tag++) {
        Node* theNode = theDomain->getNode(tag);
        if(theNode == 0) continue;
        const Vector& coord = theNode->getCrds();
        if(coord.Size() < 2) {
            continue;
        }
        const Vector& disp = theNode->getTrialDisp();
        if(disp.Size() < 2) {
            continue;
        }
        double x = coord(0);        // initial
        double y = coord(1);
        x += disp(0);               // current
        y += disp(1);
        
        // if out of boundary
        if(x<bound(0) || x>bound(2) || y<bound(1) || y>bound(3)) {
            removelist.insert(tag);
        }
    }

    // remove nodes
    for(int i=0; i<removelist.Size(); i++) {
        Node* theNode = theDomain->removeNode(removelist(i));
        if(theNode != 0) {
            delete theNode;
        }
        Pressure_Constraint* thePC = theDomain->removePressure_Constraint(removelist(i));
        if(thePC != 0) {
            delete thePC;
        }
    }
}

void
PFEMMesher2D::initializeTri(triangulateio& in)
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
PFEMMesher2D::freeTri(triangulateio& in)
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
PFEMMesher2D::freeTriOut(triangulateio& in)
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
