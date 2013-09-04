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
#include <PFEMElement2DCompressible.h>
#include <PFEMElement2DBubble.h>
#include <NDMaterial.h>
#include <Tri31.h>
#include <set>
#include <fstream>
#include <iostream>
#include <Timer.h>
#include <algorithm>
#include <elementAPI.h>

PFEMMesher2D::PFEMMesher2D()
    :PI(3.1415926535897932384626433), bound(4), frameBase(2), Lspan(), Height(), avesize(0.0),
     fluidNodes()
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
    sprintf(s,"Qzqpa%.60f",maxarea);
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
PFEMMesher2D::discretize(int startnode, char type, int n,
                         int nth, int nthfloor, int ndf,
                         const Vector& fix, const Vector& vel, 
                         const Vector& mass, Domain* theDomain, int& endnode)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::discretize\n";
        return -1;
    }

    if(n<0) return 0;
    if(nth<0 || nthfloor<0) return 0;
    if(nth>Lspan.Size()) return 0;
    if(nthfloor>Height.Size()) return 0;    

    double h = 0.0;
    double x = frameBase(0), y = frameBase(1), angle = 0.0;
    int num = 0.0;
    if(type=='b' || type=='B') {
        num = n-2;
        h = Lspan(nth-1)/n;
        for(int i=1; i<nth; i++) {
            x += Lspan(i-1);
        }
        x += h;
        for(int j=0; j<nthfloor; j++) {
            y += Height(j);
        }

    } else {
        num = n;
        h = Height[nthfloor-1]/n;
        for(int i=0; i<nth; i++) {
            x += Lspan(i);
        }
        for(int j=1; j<nthfloor; j++) {
            y += Height(j-1);
        }
        if(nthfloor>1) {
            y+=h;
            num--;
        }
        angle = 90.0;
    }
    Vector linebound(2);
    return this->discretize(startnode, x, y, h, angle, num, ndf, fix, vel, mass, 
                            linebound, theDomain, endnode);

    
}

int
PFEMMesher2D::doTriangulation(int starteletag, double alpha, const ID& groups, 
                              const ID& addgroups,Domain* theDomain, 
                              double rho, double mu, double b1, double b2, 
                              double thk, double kappa, int type)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }

    //Timer timer;
    // do triangulation
    ID eles;
    int res = 0;
    res = doTriangulation(alpha,groups,addgroups,theDomain, eles);
    if(res < 0) {
        opserr<<"WARNING: failed to do triangulation --";
        opserr<<"PFEMMesher2D::soTriangulation\n";
        return res;
    }

    // add PFEM elements
    //timer.start();
    int numeles = eles.Size()/3;
    if(numeles == 0) return 0;
    int etag = starteletag-1;
    for(int i=0; i<numeles; i++) {
        Element* theEle = 0;
        if(type == 1) {
            theEle = new PFEMElement2D(++etag, eles(3*i), eles(3*i+1), eles(3*i+2),
                                                      rho, mu, b1, b2, thk);
        } else if(type == 3) {
            theEle = new PFEMElement2DCompressible(++etag, eles(3*i), eles(3*i+1), eles(3*i+2),
                                                   rho, mu, b1, b2, thk, kappa);
        } else if(type == 4) {
            theEle = new PFEMElement2DBubble(++etag, eles(3*i), eles(3*i+1), eles(3*i+2),
                                             rho, mu, b1, b2, thk, kappa);
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
    //timer.pause();
    //opserr<<"create PFEM elements :"<<timer;
    return etag;
    
}

int
PFEMMesher2D::doTriangulation(int starteletag, double alpha, const ID& groups, 
                              const ID& addgroups,Domain* theDomain, 
                              double t, const char* type, int matTag,
                              double p, double rho, double b1, double b2)
{
    avesize = 0.0;
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::doTriangulation\n";
        return -1;
    }
    //Timer timer;
    // do triangulation
    ID eles;
    int res = doTriangulation(alpha,groups,addgroups,theDomain, eles);
    if(res < 0) {
        opserr<<"WARNING: failed to do triangulation --";
        opserr<<"PFEMMesher2D::doTriangulation\n";
        return res;
    }

    NDMaterial *theMaterial = OPS_GetNDMaterial(matTag);
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
    return etag;
    
}

int
PFEMMesher2D::doTriangulation(double alpha, const ID& groups, const ID& addgroups,
                              Domain* theDomain, ID& eles, bool o2)
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
        int type = it->second;
        bool haveit = false;
        for(int i=0; i<groups.Size()+addgroups.Size(); i++) {
            if((i<groups.Size() && type==groups(i)) ||
               (i>=groups.Size() && type==addgroups(i-groups.Size()))) {
                haveit = true;
                break;
            }
        }
        if(!haveit) continue;

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
    char s1[] = "Qzv";
    char s2[] = "Qzvo2";
    char* s = 0;
    if(o2) {
        s = s2;
    } else {
        s = s1;
    }
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
                if(o2) {
                    for(int j=3; j<out.numberofcorners; j++) {
                        int tag = out.trianglelist[out.numberofcorners*i+j];
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
        eles.resize(out.numberoftriangles*out.numberofcorners);

        // copy triangles
        for(int i=0; i<out.numberoftriangles; i++) {
            for(int j=0; j<3; j++) {
                int tag = p2nd(out.trianglelist[out.numberofcorners*i+j]);
                eles(out.numberofcorners*i+j) = tag;
            }
            if(o2) {
                for(int j=3; j<out.numberofcorners; j++) {
                    int tag = out.trianglelist[out.numberofcorners*i+j];
                    eles(out.numberofcorners*i+j) = tag;
                }
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
PFEMMesher2D::save(const char* filename, const ID& snodes, Domain* theDomain)
{
    if(theDomain == 0) {
        opserr<<"WARNING: null domain";
        opserr<<" -- PFEMMesher2D::save\n";
        return -1;
    }

    // write nodes
    std::ofstream file2(std::string(filename).append(".node").c_str());
    file2<<theDomain->getCurrentTime();
    for(int i=0; i<8; i++) {
        file2 << " " << 0;
    }
    file2<<"\n";

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
            file2<<ntag<<" "<<coord(0)+disp(0)<<" "<<coord(1)+disp(1)<<" ";
            int type = 1;
            for(int i=0; i<snodes.Size()/2; i++) {
                if(ntag>=snodes(2*i) && ntag<=snodes(2*i+1)) {
                    type = i+2;
                    break;
                }
            }
            file2<<type<<" ";
            file2<<vel(0)<<" "<<vel(1)<<" "<<accel(0)<<" "<<accel(1);
            file2<<" "<<pvel(2)<<"\n";
        }
    }

    // nodes
    int type = snodes.Size()/2+2;
    for(int i=0; i<snodes.Size()/2; i++) {
        for(int tag=snodes(2*i); tag<=snodes(2*i+1); tag++) {
            Node* theNode = theDomain->getNode(tag);
            if(theNode != 0) {
                std::map<int,int>::iterator it = nodeIndex.find(tag);
                if(it != nodeIndex.end()) continue;
                nodeIndex[tag] = index++;
                const Vector& coord = theNode->getCrds();
                const Vector& disp = theNode->getDisp();
                const Vector& vel = theNode->getVel();
                const Vector& accel = theNode->getAccel();
                file2<<tag<<" "<<coord(0)+disp(0)<<" "<<coord(1)+disp(1)<<" "<<type<<" ";
                file2<<vel(0)<<" "<<vel(1)<<" "<<accel(0)<<" "<<accel(1);
                file2<<" "<<0<<"\n";
            }
        }
    }
    file2.close();

    // max number of nodes in elements
    int maxnum = 0;
    ElementIter& theEles = theDomain->getElements();
    Element* theEle = 0;
    while((theEle = theEles()) != 0) {
        int num = 0;
        const ID& ntags = theEle->getExternalNodes();
        for(int i=0; i<ntags.Size(); i++) {
            std::map<int,int>::iterator it = nodeIndex.find(ntags(i));
            if(it == nodeIndex.end()) continue;
            num++;
        }
        if(num > maxnum) maxnum = num;
    }

    // write elements
    std::ofstream file1(std::string(filename).append(".ele").c_str());
    ElementIter& theEles1 = theDomain->getElements();
    while((theEle = theEles1()) != 0) {
        int num = 0;
        const ID& ntags = theEle->getExternalNodes();
        for(int i=0; i<ntags.Size(); i++) {
            std::map<int,int>::iterator it = nodeIndex.find(ntags(i));
            if(it == nodeIndex.end()) continue;
            num++;
            file1<<ntags(i)<<" ";
        }
        for(int i=num; i<maxnum; i++) {
            file1<<"nan ";
        }
        for(int i=0; i<ntags.Size(); i++) {
            std::map<int,int>::iterator it = nodeIndex.find(ntags(i));
            if(it == nodeIndex.end()) continue;
            file1<<it->second<<" ";
        }
        for(int i=num; i<maxnum; i++) {
            file1<<"nan ";
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
PFEMMesher2D::setFrame(double x1, double y1, const Vector& span, 
                       const Vector& height)
{
    frameBase(0) = x1;
    frameBase(1) = y1;
    Lspan = span;
    Height = height;
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
                    Node* pnode = theDomain->getNode(thePC->getPressureNode());
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
PFEMMesher2D::removeOutBoundNodes(const ID& nodes, Domain* theDomain)
{
    // check nodes
    ID removelist;
    for(int i=0; i<nodes.Size()/2; i++) {
        for(int tag=nodes(2*i); tag<=nodes(2*i+1); tag++) {
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
                removelist[removelist.Size()] = tag;
            }
        }
    }

    // remove nodes
    for(int i=0; i<removelist.Size(); i++) {
        Node* theNode = theDomain->removeNode(removelist[i]);
        if(theNode != 0) {
            delete theNode;
        }
        Pressure_Constraint* thePC = theDomain->removePressure_Constraint(removelist[i]);
        if(thePC != 0) {
            delete thePC;
        }
    }
}

int
PFEMMesher2D::addPC(const ID& nodes, int pndf, int startpnode, Domain* theDomain, int& endpnode)
{
    if(theDomain==0) return -1;

    endpnode = startpnode-1;
    for(int i=0; i<nodes.Size()/2; i++) {
        for(int j=nodes(2*i); j<=nodes(2*i+1); j++) {
            Node* theNode = theDomain->getNode(j);
            if(theNode==0) continue;
            Pressure_Constraint* thePC = theDomain->getPressure_Constraint(j);
            if(thePC == 0) {
                thePC = new Pressure_Constraint(j, ++endpnode, pndf);
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
    }

    return 0;
}

void
PFEMMesher2D::setNodes(const ID& newnodes, int type, bool series, int action)
{
    if(action == 0) fluidNodes.clear();

    if(series) {
        for(int i=0; i<newnodes.Size(); i++) {
            int nd = newnodes(i);
            if(action == 2) {
                fluidNodes.erase(nd);
            } else {
                fluidNodes[nd] = type;
            }
        }
        
    } else {

        for(int i=0; i<newnodes.Size()/2; i++) {
            for(int nd=newnodes(2*i); nd<=newnodes(2*i+1); nd++) {
                if(action == 2) {
                    fluidNodes.erase(nd);
                } else {
                    fluidNodes[nd] = type;
                }
            }
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
