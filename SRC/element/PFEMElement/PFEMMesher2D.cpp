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
                                                                        
// $Revision: 1.1 $
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
#include <PFEMElement2DCompressible.h>
#include <PFEMElement2DBubble.h>
#include <NDMaterial.h>
#include <Tri31.h>
#include <set>
#include <fstream>
#include <iostream>
//#include <Timer.h>
#include <algorithm>
#include <elementAPI.h>
#include <MeshRegion.h>

double PFEMMesher2D::PI = 3.1415926535897932384626433;

PFEMMesher2D::PFEMMesher2D()
{
}

PFEMMesher2D::~PFEMMesher2D()
{
}

// PSPG
int 
PFEMMesher2D::discretize(int startnodetag, const Vector& points, const Vector& segments, 
                         const Vector& holes, double maxarea, int ndf, const Vector& fix, 
                         const Vector& vel, const Vector& mass, Domain* theDomain, int& endnodetag)
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
        if(segments(i) > maxpt) maxpt = (int)segments(i);
        in.segmentlist[i] = (int)segments(i);
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
    sprintf(s,"DQzqpa%.60f",maxarea);
    // opserr<<"discretize : ";
    // opserr<<"Start Delaunay Triangulation -- If got segmentation fault, please check inputs \n";
    triangulate(s, &in, &out, &vout);
    // opserr<<"Finish Delaunay Triangulation\n";
    //timer.pause();
    //opserr<<timer;
    freeTri(in);

    // fix, vel and mass
    int numvel = vel.Size();
    if(numvel > ndf) numvel = ndf;
    Vector tvel(ndf);
    for(int j=0; j<numvel; j++) {
        tvel(j) = vel(j);
    }
    int numfix = fix.Size();
    if(numfix > ndf) numfix = ndf;

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

        // initial velocity
        theNode->setTrialVel(tvel);
        theNode->commitState();
        
        // fix
        for(int j=0; j<numfix; j++) {
            if(fix(j) != 0) {
                SP_Constraint* theSP = new SP_Constraint(newtag, j, 0.0, true);
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
            for(int j=0; j<nummass; j++) {
                theMass(j,j) = mass(j);
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

    endnodetag = newtag;

    return 0;
}

// particles
int 
PFEMMesher2D::discretize(int startnodetag, const Vector& points, 
                         int ndf, const Vector& fix, 
                         const Vector& vel, const Vector& mass, 
                         Domain* theDomain, int& endnodetag)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }

    // number of points
    int num = points.Size()/2;

    // no points
    if(num < 1) {
        return 0;
    }

    // fix, vel and mass
    int numvel = vel.Size();
    if(numvel > ndf) numvel = ndf;
    Vector tvel(ndf);
    for(int j=0; j<numvel; j++) {
        tvel(j) = vel(j);
    }
    int numfix = fix.Size();
    if(numfix > ndf) numfix = ndf;

    // read outputs
    int newtag = startnodetag-1;
    for(int i=0; i<num; i++) {
        double x = points(2*i);
        double y = points(2*i+1);

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

        // initial velocity
        theNode->setTrialVel(tvel);
        theNode->commitState();
        
        // fix
        for(int j=0; j<numfix; j++) {
            if(fix(j) != 0) {
                SP_Constraint* theSP = new SP_Constraint(newtag, j, 0.0, true);
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
            for(int j=0; j<nummass; j++) {
                theMass(j,j) = mass(j);
            }
            if(theNode->setMass(theMass) < 0) {
                opserr<<"WARNING: failed to set mass of node "<<newtag;
                opserr<<" -- PFEMMesher2D::discretize\n";
                return -1;
            }
        }
    }
    
    endnodetag = newtag;

    return 0;
}



// rectangle
int 
PFEMMesher2D::discretize(int startnodetag, double x1, double y1, double hx, double hy, double angle,
                         int nx, int ny, int ndf, const Vector& fix, const Vector& vel, 
                         const Vector& mass, const Vector& boundary,
                         Domain* theDomain, int& endnodetag)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }
    if(boundary.Size() < 4) {
        opserr<<"WARNING: boundary vector must have 4 entries";
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
    int start=0, end=ny;
    if(!boundary(0)) start++;
    if(!boundary(1)) end--;
    Vector linebound(2);
    linebound(0) = boundary(2);
    linebound(1) = boundary(3);
    for(int i=start; i<=end; i++) {
        double lx = i*lhx+x1;
        double ly = i*lhy+y1;
        int res = discretize(++newtag, lx, ly, hx, angle, nx, ndf, fix, vel, 
                             mass, linebound, theDomain, newtag);
        if(res < 0) return -1;
    }

    endnodetag = newtag;
    return 0;
}


// line
int 
PFEMMesher2D::discretize(int startnodetag, double x1, double y1, double h, double angle, 
                         int num, int ndf, const Vector& fix, const Vector& vel, const Vector& mass, 
                         const Vector& boundary, Domain* theDomain, int& endnodetag)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }
    if(boundary.Size() < 2) {
        opserr<<"WARNING: boundary vector must have 2 entries";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }

    if(num < 0) return 0;
    double s = sin(angle*PI/180.0);
    double c = cos(angle*PI/180.0);

    double hx = c*h;
    double hy = s*h;

    int numvel = vel.Size();
    if(numvel > ndf) numvel = ndf;
    Vector tvel(ndf);
    for(int j=0; j<numvel; j++) {
        tvel(j) = vel(j);
    }
    int numfix = fix.Size();
    if(numfix > ndf) numfix = ndf;

    int newtag = startnodetag-1;
    int start = 0, end = num;
    if(!boundary(0)) start = 1;
    if(!boundary(1)) end = num-1;
    
    for(int i=start; i<=end; i++) {
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

        // initial velocity
        theNode->setTrialVel(tvel);
        theNode->commitState();

        // fix
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

    endnodetag = newtag;
    return 0;
}

// triangle
int
PFEMMesher2D:: discretize(int startnode, double x1, double y1, double x2, double y2, 
                          double x3, double y3, int n1, int n2, int ndf, const Vector& fix,
                          const Vector& vel, const Vector& mass, const Vector& boundary,
                          Domain* theDomain, int& endnode)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }
    if(boundary.Size() < 3) {
        opserr<<"WARNING: boundary vector must have 3 entries";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }

    if(n1<=0 || n2<=0) return 0;

    double h1 = 1.0/n1;
    double h2 = 1.0/n2;

    int numvel = vel.Size();
    if(numvel > ndf) numvel = ndf;
    Vector tvel(ndf);
    for(int j=0; j<numvel; j++) {
        tvel(j) = vel(j);
    }

    endnode = startnode-1;
    for(int i=0; i<=n1; i++) {
        if(!boundary(0) && i==0) continue;
        double L1 = i*h1;
        for(int j=0; j<=n2; j++) {
            if(!boundary(1) && j==0) continue;
            double L2 = j*h2;
            double L3 = 1.0-L1-L2;
            if(L3<-1e-6) continue;
            if(!boundary(2) && fabs(L3)<1e-6) continue;
            double x = L1*x1+L2*x2+L3*x3;
            double y = L1*y1+L2*y2+L3*y3;

            // define node
            Node* theNode = new Node(++endnode, ndf, x, y);
            if (theNode == 0) {
                opserr << "WARNING ran out of memory creating node\n";
                opserr << "node: " << endnode << "\n";
                opserr << "PFEMMesher2D::discretize\n";
                return -1;
            }
            if(theDomain->addNode(theNode) == false) {
                opserr << "WARNING failed to add node to the domain\n";
                opserr << "node: " << endnode << "\n";
                opserr << "PFEMMesher2D::discretize\n";
                delete theNode; // otherwise memory leak
                return -1;
            }

            // initial velocity
            theNode->setTrialVel(tvel);
            theNode->commitState();

            // fix
            int numfix = fix.Size();
            if(numfix > ndf) numfix = ndf;
            for(int i=0; i<numfix; i++) {
                if(fix(i) != 0) {
                    SP_Constraint* theSP = new SP_Constraint(endnode, i, 0.0, true);
                    if(theSP == 0) {
                        opserr << "WARNING ran out of memory creating SP_Constraint\n";
                        opserr << "node: " << endnode << "\n";
                        opserr << "PFEMMesher2D::discretize\n";
                        return -1;
                    }
                    if(theDomain->addSP_Constraint(theSP) == false) {
                        opserr << "WARNING failed to add SP_Constraint to the domain\n";
                        opserr << "node: " << endnode << "\n";
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
                    opserr<<"WARNING: failed to set mass of node "<<endnode;
                    opserr<<" -- PFEMMesher2D::discretize\n";
                    return -1;
                }
            }
        }
    }
    
    return 0;
}


// circle
int 
PFEMMesher2D::discretize(int startnode, double xc, double yc, double r1, double r2,
                         int nc, int nr, int ndf, const Vector& fix, const Vector& vel,
                         const Vector& mass, const Vector& boundary, const Vector& angle,
                         Domain* theDomain, int& endnode)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }
    if(boundary.Size() < 2) {
        opserr<<"WARNING: boundary vector must have 2 entries";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }
    if(angle.Size() < 2) {
        opserr<<"WARNING: angle vector must have 2 entries";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }

    if(nc<=0 || nr<0) return 0;
    double angle1 = angle(0);
    double angle2 = angle(1);
    if(angle1 < 0) angle1 = 0;
    if(angle2 > 360) angle2 = 360;
    angle1 = angle1*PI/180.0;
    angle2 = angle2*PI/180.0;

    double hc = (angle2-angle1)/nc;
    double hr = 0.0;
    if(nr>0) hr = (r2-r1)/nr;

    
    int numvel = vel.Size();
    if(numvel > ndf) numvel = ndf;
    Vector tvel(ndf);
    for(int j=0; j<numvel; j++) {
        tvel(j) = vel(j);
    }

    endnode = startnode-1;
    int start = 0, end = nr;
    if(!boundary(0)) start++;
    if(!boundary(1)) end--;
    for(int i=start; i<=end; i++) {
        double ri = hr*i+r1;
        for(int j=0; j<=nc; j++) {
            double angle = angle1+j*hc;
            double x = ri*cos(angle)+xc;
            double y = ri*sin(angle)+yc;

            // define node
            Node* theNode = new Node(++endnode, ndf, x, y);
            if (theNode == 0) {
                opserr << "WARNING ran out of memory creating node\n";
                opserr << "node: " << endnode << "\n";
                opserr << "PFEMMesher2D::discretize\n";
                return -1;
            }
            if(theDomain->addNode(theNode) == false) {
                opserr << "WARNING failed to add node to the domain\n";
                opserr << "node: " << endnode << "\n";
                opserr << "PFEMMesher2D::discretize\n";
                delete theNode; // otherwise memory leak
                return -1;
            }

            // initial velocity
            theNode->setTrialVel(tvel);
            theNode->commitState();

            // fix
            int numfix = fix.Size();
            if(numfix > ndf) numfix = ndf;
            for(int i=0; i<numfix; i++) {
                if(fix(i) != 0) {
                    SP_Constraint* theSP = new SP_Constraint(endnode, i, 0.0, true);
                    if(theSP == 0) {
                        opserr << "WARNING ran out of memory creating SP_Constraint\n";
                        opserr << "node: " << endnode << "\n";
                        opserr << "PFEMMesher2D::discretize\n";
                        return -1;
                    }
                    if(theDomain->addSP_Constraint(theSP) == false) {
                        opserr << "WARNING failed to add SP_Constraint to the domain\n";
                        opserr << "node: " << endnode << "\n";
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
                    opserr<<"WARNING: failed to set mass of node "<<endnode;
                    opserr<<" -- PFEMMesher2D::discretize\n";
                    return -1;
                }
            }
        }
    }
    
    return 0;
}


// frame
int 
PFEMMesher2D::discretize(int startnode, double x, double y,
                         double hcol, double hbeam, int ndf,
                         const Vector& spans, const Vector& heights,
                         const Vector& fix, const Vector& colvel, 
                         const Vector& beamvel, const Vector& colmass, 
                         const Vector& beammass, const Vector& roofmass,
                         Domain* theDomain, ID& nodelist)

{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }

    // quick return
    int ncol = spans.Size()+1;
    int nbeam = heights.Size();
    if(ncol<2 || nbeam<1 || hcol<=0 || hbeam<=0 || ndf <1) {
        return 0;
    }
    for(int i=0; i<ncol-1;i++) {
        if(spans(i) <= 0) {
            opserr<<"WARNING: "<<i<<"th span <= 0 ";
            opserr<<"-- PFEMMesher2D::discretize";
            return -1;
        }
    }
    for(int i=0; i<nbeam;i++) {
        if(heights(i) <= 0) {
            opserr<<"WARNING: "<<i<<"th height <= 0 ";
            opserr<<"-- PFEMMesher2D::discretize";
            return -1;
        }
    }

    // nodelist size
    nodelist.resize(ncol*(nbeam+1)+(ncol-1)*nbeam*2);
    nodelist(0) = startnode;

    // boundary and fix list
    Vector nobound(2), nofix(ndf);

    // columns
    int ndindex = 0;
    int endnode = startnode;
    Vector pos(2);
    pos(0) = x;
    pos(1) = y;
    for(int i=0; i<ncol; i++) {
        pos(1) = y;
        for(int j=0; j<nbeam; j++) {

            // start node
            int res = 0;
            if(j==0) {
                res = discretize(nodelist(ndindex),pos,ndf,fix,
                                 colvel,colmass,theDomain,endnode);
            } else {
                res = discretize(nodelist(ndindex),pos,ndf,nofix,
                                 colvel,colmass,theDomain,endnode);
            }
            if(res < 0) return res;

            // mid nodes
            int num = floor(heights(j)/hcol+0.5);
            res = discretize(nodelist(ndindex)+1,pos(0),pos(1),hcol,90.0,
                             num,ndf,nofix,colvel,colmass,
                             nobound,theDomain,endnode);
            if(res < 0) return res;
            ndindex++;
            nodelist(ndindex) = endnode+1;
            pos(1) += heights(j);

            // end node
            if(j == nbeam-1) {
                res = discretize(nodelist(ndindex),pos,ndf,nofix,
                                 colvel,colmass,theDomain,endnode);
                if(res < 0) return res;
                ndindex++;
                nodelist(ndindex) = endnode+1;
            }
        }
        if(i<ncol-1) pos(0) += spans(i);
    }

    // beams
    pos(0) = x;
    pos(1) = y;
    for(int j=0; j<nbeam; j++) {
        pos(0) = x;
        pos(1) += heights(j);        
        for(int i=0; i<ncol-1; i++) {
            int num = floor(spans(i)/hbeam+0.5);
            
            // mid nodes
            int res;
            if(j<nbeam-1) {
                res = discretize(nodelist(ndindex),pos(0),pos(1),hbeam,0.0,
                                 num,ndf,nofix,beamvel,beammass,
                                 nobound,theDomain,endnode);
            } else {
                res = discretize(nodelist(ndindex),pos(0),pos(1),hbeam,0.0,
                                 num,ndf,nofix,beamvel,roofmass,
                                 nobound,theDomain,endnode);
            }
            if(res < 0) return res;
            ndindex++;
            nodelist(ndindex) = endnode;
            if(ndindex<nodelist.Size()-1) {
                ndindex++;
                nodelist(ndindex) = endnode+1;
            }
            pos(0) += spans(i);
        }
    }

    return 0;
}

// linear elements
int
PFEMMesher2D::doTriangulation(int starteletag, double alpha, const ID& groups, 
                              const ID& addgroups,Domain* theDomain, 
                              double rho, double mu, double b1, double b2, 
                              double thk, double kappa, int type, int& endele)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }

    //Timer timer;
    //timer.start();

    // do triangulation
    ID eles;
    int res = 0;
    res = doTriangulation(alpha,groups,addgroups,theDomain,eles);
    if(res < 0) {
        opserr<<"WARNING: failed to do triangulation --";
        opserr<<"PFEMMesher2D::soTriangulation\n";
        return res;
    }

    // add PFEM elements
    int numeles = eles.Size()/3;
    if(numeles == 0) return 0;
    int etag = starteletag-1;
    for(int i=0; i<numeles; i++) {
        Element* theEle = 0;
        if(type == 1) {
            theEle = new PFEMElement2D(++etag, eles(3*i), eles(3*i+1), eles(3*i+2),rho, mu, b1, b2, thk);
        } else if(type == 3) {
            theEle = new PFEMElement2DCompressible(++etag, eles(3*i), eles(3*i+1), eles(3*i+2),rho, mu, b1, b2, thk, kappa);
        } else if(type == 4) {
            theEle = new PFEMElement2DBubble(++etag, eles(3*i), eles(3*i+1), eles(3*i+2),rho, mu, b1, b2, thk, kappa);
        }
        
        if(theEle == 0) {
            opserr<<"WARNING: no enough memory -- ";
            opserr<<" -- PFEMMesher2D::doTriangulation\n";
            return -1;
        }
        if(theDomain->addElement(theEle) == false) {
            opserr<<"WARNING: failed to add element to domain -- ";
            opserr<<" -- PFEMMesher2D::doTriangulation\n";
            delete theEle;
            return -1;
        }
    }

    // identify
    identify(b2,theDomain);
    //timer.pause();
    //opserr<<"meshing :"<<timer.getCPU()<<"\n";

    endele = etag;
    return res;
    
}

// solid elements
int
PFEMMesher2D::doTriangulation(int starteletag, double alpha, const ID& groups, 
                              const ID& addgroups,Domain* theDomain, 
                              double t, const char* type, int matTag,
                              double p, double rho, double b1, double b2,
                              int& endele)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }
    //Timer timer;
    // do triangulation
    ID eles;
    int res = doTriangulation(alpha,groups,addgroups,theDomain,eles);
    if(res < 0) {
        opserr<<"WARNING: failed to do triangulation --";
        opserr<<"PFEMMesher2D::doTriangulation\n";
        return res;
    }

    NDMaterial *theMaterial = OPS_getNDMaterial(matTag);
    if(theMaterial == 0) {
        opserr << "WARNING:  Material " << matTag << "not found\n";
        opserr<<"PFEMMesher2D::doTriangulation\n";
        return -1;
    }

    // add Tri31 elements
    //timer.start();
    int numeles = eles.Size()/3;
    if(numeles == 0) return 0;
    int etag = starteletag-1;
    for(int i=0; i<numeles; i++) {
        Tri31* theEle = new Tri31(++etag, eles(3*i), eles(3*i+1), eles(3*i+2),
                                  *theMaterial, type, t, p, rho, b1, b2);

        if(theEle == 0) {
            opserr<<"WARNING: no enough memory -- ";
            opserr<<" -- PFEMMesher2D::doTriangulation\n";
            return -1;
        }
        if(theDomain->addElement(theEle) == false) {
            opserr<<"WARNING: failed to add element to domain -- ";
            opserr<<" -- PFEMMesher2D::doTriangulation\n";
            delete theEle;
            return -1;
        }
    }
    //timer.pause();
    //opserr<<"create PFEM elements :"<<timer;

    endele = etag;

    return res;
    
}

// core function
int
PFEMMesher2D::doTriangulation(double alpha, const ID& groups, const ID& addgroups,
                              Domain* theDomain, ID& eles)
{
    //Timer theTimer;
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }

    // fluidnodes
    std::map<int,int> fluidNodes;
    getNodes(groups,fluidNodes,theDomain);
    getNodes(addgroups,fluidNodes,theDomain);

    // structs
    triangulateio in, out, vout;
    initializeTri(in);
    initializeTri(out);
    initializeTri(vout);

    // number of nodes
    int numnodes = fluidNodes.size();
    in.numberofpoints = numnodes;
    if(in.numberofpoints < 3) {
        return 0;
    }

    // get space for pointlist
    in.pointlist = (double*) calloc(in.numberofpoints*2, sizeof(double));
    if(in.pointlist == NULL) {
        opserr<<"WARNING: no enough memory -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }
    //opserr<<"in.numberofpoints = "<<in.numberofpoints<<"\n";
    // copy nodes to pointlist, use current coordinates
    int numpoints = 0;
    ID p2nd(0, in.numberofpoints);
    for(std::map<int,int>::iterator it=fluidNodes.begin(); it!=fluidNodes.end(); it++) {

        // pointer to node
        int ndtag = it->first;
        // int type = it->second;
        // bool haveit = false;
        // for(int i=0; i<groups.Size()+addgroups.Size(); i++) {
        //     if((i<groups.Size() && type==groups(i)) ||
        //        (i>=groups.Size() && type==addgroups(i-groups.Size()))) {
        //         haveit = true;
        //         break;
        //     }
        // }
        // if(!haveit) continue;

        Node* node = theDomain->getNode(ndtag);
        if(node == 0) continue;

        // get coordinates
        const Vector& coord = node->getCrds();
        if(coord.Size() < 2) {
            opserr<<"WARNING: 2d is required -- PFEMMesher2D::doTriangulation\n";
            return -1;
        }
        const Vector& disp = node->getTrialDisp();
        if(disp.Size() < 2) {
            opserr<<"WARNING: 2 ndf is required -- PFEMMesher2D::doTriangulation\n";
            return -1;
        }

        // set pointlist
        for(int j=0; j<2; j++) {
            in.pointlist[2*numpoints+j] = coord(j)+disp(j);
        }
        p2nd[numpoints++] = ndtag;

    }
    in.numberofpoints = numpoints;

    // Delaunay Triangulation
    char s[] = "Qzv";
    //theTimer.start();
    // opserr<<"Start Delaunay Triangulation -- If got segmentation fault, please check inputs \n";
    triangulate(s, &in, &out, &vout);
    // opserr<<"Finish Delaunay Triangulation\n";
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
        Vector radius(out.numberoftriangles);
        double avesize = 0.0;
        for(int i=0; i<out.numberoftriangles; i++) {

            // circumcenter of triangle
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
        eles.resize(out.numberoftriangles*out.numberofcorners);
        for(int i=0; i<out.numberoftriangles; i++) {
            if(radius(i) / avesize <= alpha) {
                // pass the test

                // check if all nodes are additional nodes
                int add[3] = {-1,-1,-1};
                for(int j=0; j<3; j++) {
                    int tag = p2nd(out.trianglelist[out.numberofcorners*i+j]);
                    int type = fluidNodes[tag];
                    for(int k=0; k<groups.Size(); k++) {
                        if(type == groups(k)) {
                            add[j] = 0;
                            break;
                        }
                    }
                }
                    
                // add ele
                // if(add[0]==add[1] && add[1]==add[2] && add[2]!=0) continue;
                if(add[0]!=0 && add[1]!=0 && add[2]!=0) continue;


                for(int j=0; j<3; j++) {
                    int tag = p2nd(out.trianglelist[out.numberofcorners*i+j]);
                    eles(num++) = tag;
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
        eles.resize(out.numberoftriangles*out.numberofcorners);

        // copy triangles
        for(int i=0; i<out.numberoftriangles; i++) {
            for(int j=0; j<3; j++) {
                int tag = p2nd(out.trianglelist[out.numberofcorners*i+j]);
                eles(out.numberofcorners*i+j) = tag;
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
PFEMMesher2D::save(const char* filename, Domain* theDomain, int maxelenodes)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::save\n";
        return -1;
    }

    // get nodes
    std::map<int,int> fluidNodes;
    ID regions;
    theDomain->getRegionTags(regions);
    getNodes(regions,fluidNodes,theDomain);

    std::ofstream outfile(std::string(filename).append(".node").c_str());
    outfile<<theDomain->getCurrentTime()<<" ";
    for(int i=0; i<8; i++) outfile<<0<<" ";
    outfile<<"\n";

    // write nodes
    for(std::map<int,int>::iterator it=fluidNodes.begin(); it!=fluidNodes.end(); it++) {
        int tag = it->first;
        int type = it->second;
        Node* node = theDomain->getNode(tag);
        if(node == 0) {
            opserr<<"WARNING: node "<<tag<<" dose not exist -- ";
            opserr<<" -- PFEMMesher2D::save\n";
            return -1;
        }
        const Vector& coord = node->getCrds();
        const Vector& disp = node->getDisp();
        const Vector& vel = node->getVel();
        const Vector& accel = node->getAccel();
        outfile<<tag<<" "<<coord(0)+disp(0)<<" "<<coord(1)+disp(1)<<" ";
        outfile<<type<<" ";
        outfile<<vel(0)<<" "<<vel(1)<<" "<<accel(0)<<" "<<accel(1);
        Pressure_Constraint* thePC = theDomain->getPressure_Constraint(tag);
        if(thePC != 0) {
            outfile<<" "<<thePC->getPressure()<<"\n";
        } else {
            outfile<<" "<<0<<"\n";
        }
    }
    outfile.close();

    // write elements
    outfile.open(std::string(filename).append(".ele").c_str());
    ElementIter& theEles = theDomain->getElements();
    Element* theEle = 0;
    while((theEle = theEles()) != 0) {
        const ID& ntags = theEle->getExternalNodes();
        int num = 0;
        for(int i=0; i<ntags.Size(); i++) {
            if(fluidNodes.find(ntags(i)) != fluidNodes.end()) {
                outfile<<ntags(i)<<" ";
                num++;
            }
        }
        for(int i=num; i<maxelenodes; i++) {
            outfile<<ntags(num-1)<<" ";
        }
        outfile<<"\n";
    }
    outfile.close();

    return 0;
}

int 
PFEMMesher2D::save(const char* filename, const ID& nodeRegs,
                   const ID& eleRegs, Domain* theDomain, int maxelenodes)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::save\n";
        return -1;
    }

    // get nodes
    std::map<int,int> fluidNodes, fluidEles;
    getNodes(nodeRegs,fluidNodes,theDomain);
    getElements(eleRegs,fluidEles,theDomain);

    std::ofstream outfile(std::string(filename).append(".node").c_str());
    outfile<<theDomain->getCurrentTime()<<" ";
    for(int i=0; i<8; i++) outfile<<0<<" ";
    outfile<<"\n";

    // write nodes
    for(std::map<int,int>::iterator it=fluidNodes.begin(); it!=fluidNodes.end(); it++) {
        int tag = it->first;
        int type = it->second;
        Node* node = theDomain->getNode(tag);
        if(node == 0) {
            opserr<<"WARNING: node "<<tag<<" dose not exist -- ";
            opserr<<" -- PFEMMesher2D::save\n";
            return -1;
        }
        const Vector& coord = node->getCrds();
        const Vector& disp = node->getDisp();
        const Vector& vel = node->getVel();
        const Vector& accel = node->getAccel();
        outfile<<tag<<" "<<coord(0)+disp(0)<<" "<<coord(1)+disp(1)<<" ";
        outfile<<type<<" ";
        outfile<<vel(0)<<" "<<vel(1)<<" "<<accel(0)<<" "<<accel(1);
        Pressure_Constraint* thePC = theDomain->getPressure_Constraint(tag);
        if(thePC != 0) {
            outfile<<" "<<thePC->getPressure()<<"\n";
        } else {
            outfile<<" "<<0<<"\n";
        }
    }
    outfile.close();

    // write elements
    outfile.open(std::string(filename).append(".ele").c_str());
    for(std::map<int,int>::iterator it=fluidEles.begin(); it!=fluidEles.end(); it++) {
        int tag = it->first;
        Element* ele = theDomain->getElement(tag);
        if(ele == 0) {
            opserr<<"WARNING: element "<<tag<<" dose not exist -- ";
            opserr<<" -- PFEMMesher2D::save\n";
            return -1;
        }

        const ID& ntags = ele->getExternalNodes();
        int num = 0;
        for(int i=0; i<ntags.Size(); i++) {
            if(fluidNodes.find(ntags(i)) != fluidNodes.end()) {
                outfile<<ntags(i)<<" ";
                num++;
            }
        }
        for(int i=num; i<maxelenodes; i++) {
            outfile<<ntags(num-1)<<" ";
        }
        outfile<<"\n";
    }
    outfile.close();
    
    return 0;
}

Vector
PFEMMesher2D::calculateForces(const ID& boundary, int basenode, 
                              Vector& dragdir, Vector& liftdir, 
                              Domain* theDomain)
{
    // drag, lift and moment
    Vector forces(3);

    // check inputs
    if(dragdir.Size()<2 || liftdir.Size()<2 || boundary.Size()<2) {
        return forces;
    }

    // for each fluid elements
    ElementIter& theEles = theDomain->getElements();
    Element* theEle = 0;
    while((theEle = theEles()) != 0) {
        std::string type(theEle->getClassType());
        if(type == "PFEMElement2D") {

            // nodal position, pressure
            const ID& elenodes = theEle->getExternalNodes();
            Vector xx(3), yy(3), pressures(3);
            bool haveit[3] = {false,false,false};
            for(int i=0; i<3; i++) {
                for(int j=0; j<boundary.Size()/2; j++) {
                    if(elenodes(2*i)>=boundary(2*j) && elenodes(2*i)<=boundary(2*j+1)) {
                        haveit[i] = true;
                    }
                }
                Node* node = theDomain->getNode(elenodes(2*i));
                if(node == 0) {
                    opserr<<"WARNING: nodes "<<elenodes(2*i);
                    opserr<<" does not exist -- PFEMMesher2D::calculateForces\n";
                    forces.Zero();
                    return forces;
                }
                const Vector& coord = node->getCrds();
                if(coord.Size()<2) {
                    opserr<<"WARNING: nodes "<<elenodes(2*i);
                    opserr<<" has wrong ndm < 2 -- PFEMMesher2D::calculateForces\n";
                    forces.Zero();
                    return forces;
                }
                const Vector& disp = node->getDisp();
                if(disp.Size()<2) {
                    opserr<<"WARNING: nodes "<<elenodes(2*i);
                    opserr<<" has wrong ndf < 2 -- PFEMMesher2D::calculateForces\n";
                    forces.Zero();
                    return forces;
                }
                xx(i) = coord(0)+disp(0);
                yy(i) = coord(1)+disp(1);

                Pressure_Constraint* thePC = theDomain->getPressure_Constraint(elenodes(2*i));
                if(thePC != 0) {
                    Node* pnode = thePC->getPressureNode();
                    if(pnode != 0) {
                        const Vector& vel = pnode->getVel();
                        if(vel.Size() > 0) {
                            pressures(i) = vel(0);
                        }
                    }
                }
            }
            
            // drag and lift direction
            dragdir.Normalize();
            liftdir.Normalize();

            for(int i=0; i<3; i++) {
                int j = i+1>2?i-2:i+1;
                int k = i+2>2?i-1:i+2;
                if(haveit[j] && haveit[k]) {
                    Vector dir(2);
                    dir(0) = yy(j)-yy(k);
                    dir(1) = xx(j)-xx(k);
                    Vector dir1(2);
                    dir1(0) = xx(j)-xx(i);
                    dir1(1) = yy(j)-yy(i);
                    if((dir^dir1) < 0) {
                        dir(0) = -dir(0);
                        dir(1) = -dir(1);
                    }

                    // length of side
                    double len = dir.Norm();

                    // pressure direction
                    dir.Normalize();

                    // pressur force
                    double pf = (pressures(j)+pressures(k))*len/2.0;

                    // add pressure to drag, lift
                    forces(0) += pf*(dir^dragdir);
                    forces(1) += pf*(dir^liftdir);
                    
                    // moement : to be added 
                }
            }
        }
    }

    return forces;
}

void 
PFEMMesher2D::removeOutBoundNodes(const ID& groups, double x1, double y1, 
                                  double x2, double y2, Domain* theDomain)
{
    // get nodes
    std::map<int,int> fluidNodes;
    getNodes(groups, fluidNodes, theDomain);

    // check nodes
    ID removelist;
    for(std::map<int,int>::iterator it=fluidNodes.begin(); it!=fluidNodes.end(); it++) {
        int tag = it->first;
        // int type = it->second;
        // bool haveit = false;
        // for(int i=0; i<groups.Size(); i++) {
        //     if(type==groups(i)) {
        //         haveit = true;
        //         break;
        //     }
        // }
        // if(!haveit) continue;

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
        if(x<x1 || x>x2 || y<y1 || y>y2) {
            removelist[removelist.Size()] = tag;
        }
    }

    // remove nodes
    for(int i=0; i<removelist.Size(); i++) {
        int tag = removelist[i];
        Node* theNode = theDomain->removeNode(tag);
        if(theNode != 0) {
            delete theNode;
        }
        Pressure_Constraint* thePC = theDomain->removePressure_Constraint(removelist[i]);
        if(thePC != 0) {
            delete thePC;
        }
        fluidNodes.erase(tag);
    }
}

int
PFEMMesher2D::addPC(const ID& groups, int pndf, int startpnode, Domain* theDomain, int& endpnode)
{
    if(theDomain==0) return -1;

    // get nodes
    std::map<int,int> fluidNodes;
    getNodes(groups, fluidNodes, theDomain);

    endpnode = startpnode-1;
    for(std::map<int,int>::iterator it=fluidNodes.begin(); it!=fluidNodes.end(); it++) {
        int tag = it->first;
        // int type = it->second;
        // bool haveit = false;
        // for(int i=0; i<groups.Size(); i++) {
        //     if(type==groups(i)) {
        //         haveit = true;
        //         break;
        //     }
        // }
        // if(!haveit) continue;
        Node* theNode = theDomain->getNode(tag);
        if(theNode==0) continue;
        Pressure_Constraint* thePC = theDomain->getPressure_Constraint(tag);
        if(thePC == 0) {
            thePC = new Pressure_Constraint(tag, ++endpnode, pndf);
            if(thePC == 0) {
                opserr<<"WARNING: no enough memory for Pressure_Constraint -- ";
                opserr<<"PFEMMesher2D::addPC \n";
                return -1;
            }
            if(theDomain->addPressure_Constraint(thePC) == false) {
                opserr<<"WARNING: failed to add Pressure_Constraint to domain -- ";
                opserr<<"PFEMMesher2D::addPC\n ";
                delete thePC;
                thePC = 0;
                return -1;
            }
        }
    }

    return 0;
}

void
PFEMMesher2D::setNodes(const ID& newnodes, int type, bool series, 
                       int action, Domain* theDomain)
{
    std::set<int> ndset;
    
    // nodes in existing region
    MeshRegion* region = theDomain->getRegion(type);
    if(region != 0) {
        const ID& nodes = region->getNodes();
        for(int i=0; i<nodes.Size(); i++) {
            ndset.insert(nodes(i));
        }
    }
    
    // apply action:0-set,1-add,2-remove
    if(action == 0) ndset.clear();

    if(series) {
        for(int i=0; i<newnodes.Size(); i++) {
            int nd = newnodes(i);
            if(action == 2) {
                ndset.erase(nd);
            } else {
                ndset.insert(nd);
            }
        }
        
    } else {

        for(int i=0; i<newnodes.Size()/2; i++) {
            for(int nd=newnodes(2*i); nd<=newnodes(2*i+1); nd++) {
                if(action == 2) {
                    ndset.erase(nd);
                } else {
                    ndset.insert(nd);
                }
            }
        }

    }

    // copy to an ID
    int size = ndset.size();
    ID currentnodes(size);
    if(size > 0) {
        std::copy(ndset.begin(),ndset.end(),&currentnodes(0));
    }
    
    // set region
    if(region == 0) {
        region = new MeshRegion(type);
        if(region == 0) {
            opserr<<"WARNING: run out of memory in creating MeshRegion ";
            opserr<<"-- PFEMMesher2D::setNodes\n";
            return;
        }
        if(theDomain->addRegion(*region)==-1) {
            opserr<<"WARNING: failed to add MeshRegion to domain ";
            opserr<<"-- PFEMMesher2D::setNodes\n";
            return;            
        }
    }
    region->setNodes(currentnodes);
}

void
PFEMMesher2D::getNodes(const ID& regions, std::map<int,int>& nodes, 
                       Domain* theDomain)
{
    for(int i=0; i<regions.Size(); i++) {
        int type = regions(i);
        MeshRegion* region = theDomain->getRegion(type);
        if(region == 0) continue;
        const ID& nds = region->getNodes();
        for(int j=0; j<nds.Size(); j++) {
            nodes[nds(j)] = type;
        }
    }
}

void
PFEMMesher2D::setElements(const ID& neweles, int regTag, bool series, 
                          int action, Domain* theDomain)
{
    std::set<int> eleset;
    
    // elements in existing region
    MeshRegion* region = theDomain->getRegion(regTag);
    if(region != 0) {
        const ID& eles = region->getElements();
        for(int i=0; i<eles.Size(); i++) {
            eleset.insert(eles(i));
        }
    }
    
    // apply action
    if(action == 0) eleset.clear();

    if(series) {
        for(int i=0; i<neweles.Size(); i++) {
            int nd = neweles(i);
            if(action == 2) {
                eleset.erase(nd);
            } else {
                eleset.insert(nd);
            }
        }
        
    } else {

        for(int i=0; i<neweles.Size()/2; i++) {
            for(int nd=neweles(2*i); nd<=neweles(2*i+1); nd++) {
                if(action == 2) {
                    eleset.erase(nd);
                } else {
                    eleset.insert(nd);
                }
            }
        }

    }

    // copy to an ID
    int size = eleset.size();
    ID currenteles(size);
    if(size > 0) {
        std::copy(eleset.begin(),eleset.end(),&currenteles(0));
    }
    
    // set region
    if(region == 0) {
        region = new MeshRegion(regTag);
        if(region == 0) {
            opserr<<"WARNING: run out of memory in creating MeshRegion ";
            opserr<<"-- PFEMMesher2D::setNodes\n";
            return;
        }
        if(theDomain->addRegion(*region)==-1) {
            opserr<<"WARNING: failed to add MeshRegion to domain ";
            opserr<<"-- PFEMMesher2D::setNodes\n";
            return;            
        }
    }
    region->setElements(currenteles);
}

void
PFEMMesher2D::getElements(const ID& regions, std::map<int,int>& elements, 
                          Domain* theDomain)
{
    for(int i=0; i<regions.Size(); i++) {
        int type = regions(i);
        MeshRegion* region = theDomain->getRegion(type);
        if(region == 0) continue;
        const ID& eles = region->getElements();
        for(int j=0; j<eles.Size(); j++) {
            elements[eles(j)] = type;
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

int
PFEMMesher2D::updateNode(int tag, int dof, double value, 
                         int type, Domain* theDomain)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::updateNode\n";
        return -1;
    }
    Node* theNode = theDomain->getNode(tag);
    if(theNode == 0) {
        opserr<<"WARNNG: node "<<tag<<" does not exist ";
        opserr<<" -- PFEMMesher2D::updateNode\n";
        return -1;
    }
    dof--;
    if(dof >= theNode->getNumberDOF()) {
        opserr<<"WARNNG: dof "<<dof<<" exceed the ndf ";
        opserr<<" -- PFEMMesher2D::updateNode\n";
        return -1;
    }
    if(type==1) {   // disp
        theNode->setTrialDisp(value, dof);
    } else if(type==2) {  // vel
        Vector vel = theNode->getVel();
        vel(dof) = value;
        theNode->setTrialVel(vel);
    } else if(type==3) {  // accel
        Vector accel = theNode->getAccel();
        accel(dof) = value;
        theNode->setTrialAccel(accel);        
    }

    theNode->commitState();

    return 0;
}

void
PFEMMesher2D::identify(double g, Domain* theDomain)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::identify\n";
    }

    // disconnect pcs from structural elements
    Pressure_ConstraintIter& thePCs = theDomain->getPCs();
    Pressure_Constraint* thePC = 0;
    while((thePC = thePCs()) != 0) {
        thePC->disconnect();
    }

    // to connect pc to all elements
    Element *elePtr;
    ElementIter &theElemIter = theDomain->getElements();    
    while((elePtr = theElemIter()) != 0) {

        const ID& nodeids = elePtr->getExternalNodes();
        for(int i=0; i<nodeids.Size(); i++) {
            Pressure_Constraint* thePC = theDomain->getPressure_Constraint(nodeids(i));
            if(thePC != 0) {
                thePC->connect(elePtr->getTag(), false);
            }
        }
    }

    // get all SPs and MPs
    SP_ConstraintIter& theSPs = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint* theSP = 0;
    std::map<int,ID> sps;
    while((theSP = theSPs()) != 0) {
        sps[theSP->getNodeTag()].insert(theSP->getDOF_Number());
    }

    // get number of parameters
    int numParams = theDomain->getNumParameters();

    // move the isolated nodes
    Pressure_ConstraintIter& thePCs1 = theDomain->getPCs();
    while((thePC = thePCs1()) != 0) {

        if(thePC->isIsolated()) {
            int tag = thePC->getTag();
            Node* node = theDomain->getNode(tag);
            if(node == 0) {
                opserr<<"WARNING: node "<<tag;
                opserr<<" does not exist -- PFEMMesher2D::identify\n";
            }
            std::map<int,ID>::iterator it = sps.find(tag);

            // responses
            const Vector& disp = node->getDisp();
            const Vector& vel = node->getVel();
            const Vector& accel = node->getAccel();
            Vector ndisp = disp;
            Vector nvel = vel;
            Vector naccel = accel;
            naccel.Zero();
            naccel(1) = g;
            ndisp.addVector(1.0, vel, ops_Dt);
            ndisp.addVector(1.0, accel, 0.5*ops_Dt*ops_Dt);
            nvel.addVector(1.0, accel, ops_Dt);

            if(it != sps.end()) {
                const ID& dofs = it->second;
                for(int i=0; i<dofs.Size(); i++) {
                    ndisp(dofs(i)) = disp(dofs(i));
                    nvel(dofs(i)) = vel(dofs(i));
                    naccel(dofs(i)) = accel(dofs(i));
                }
            }


            node->setTrialDisp(ndisp);
            node->setTrialVel(nvel);
            node->setTrialAccel(naccel);
            node->commitState();

            // sensitivity
            for(int grad=0; grad<numParams; grad++) {
                Vector sensdisp(disp.Size());
                Vector sensvel(vel.Size());
                Vector sensaccel(accel.Size());

                for(int dof=0; dof<disp.Size(); dof++) {
                    bool fixed = false;
                    if(it != sps.end()) {
                        const ID& dofs = it->second;
                        for(int i=0; i<dofs.Size(); i++) {
                            if(dof == dofs(i)) {
                                fixed = true;
                                break;
                            }
                        }
                    }
                    if(!fixed) {
                        sensvel(dof) = node->getVelSensitivity(dof,grad);
                        sensdisp(dof) = node->getDispSensitivity(dof,grad)+ops_Dt*sensvel(dof);
                    }
                }

                node->saveDispSensitivity(sensdisp,grad,numParams);
                node->saveVelSensitivity(sensvel,grad,numParams);
                node->saveAccelSensitivity(sensaccel,grad,numParams);
            }
        }
    }
}

int
PFEMMesher2D::findNodeTag(Domain* theDomain) 
{
    int nodetag = 0;
    NodeIter& nodes = theDomain->getNodes();
    Node* theNode = 0;
    while((theNode = nodes()) != 0) {
        nodetag = theNode->getTag();
    }
    return nodetag+1;
}

int
PFEMMesher2D::findEleTag(Domain* theDomain) 
{
    int eletag = 0;
    ElementIter& elements = theDomain->getElements();
    Element* theEle = 0;
    while((theEle = elements()) != 0) {
        eletag = theEle->getTag();
    }
    return eletag+1;
}

void
PFEMMesher2D::removeElements(int regTag, Domain* theDomain)
{
    MeshRegion* eleReg = theDomain->getRegion(regTag);
    if(eleReg != 0) {
        const ID& regEles = eleReg->getElements();
        for(int i=0; i<regEles.Size(); i++) {
            Element* ele = theDomain->removeElement(regEles(i));
            if(ele != 0) delete ele;
        }
        eleReg->setElements(ID());
    }
}
