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
// $Source: SRC/element/PFEMElement/PFEMMesher2D.h,v $
                                                                        
// Written: Minjie Zhu
// Created: Oct. 31
//
// Description: The class PFEMMesher2D encapsulates all necessary
// routines for PFEM analysis for fluid-structure interaction
// The triangulation uses the 
// Triangle program written by Jonathan Richard Shewchuk
// http://www.cs.cmu.edu/~quake/triangle.html

//
// What: "@(#) PFEMMesher2D.h, revA"

#ifndef PFEMMesher2D_H
#define PFEMMesher2D_H

#define REAL double
#define VOID void

extern "C" {
#include <triangle.h>
}
#include <string>
#include <Vector.h>
#include <ID.h>

class Domain;


class PFEMMesher2D 
{

public:
    

    // constructor
    PFEMMesher2D();

    // destructor
    ~PFEMMesher2D();

    // discretize domain
    int discretize(int startnodetag, const Vector& points, const Vector& segments, 
                   const Vector& holes,
                   double maxarea, int ndf, const Vector& fix, const Vector& vel,
                   const Vector& mass, Domain* theDomain, int& endnodetag); // PSPG

    int discretize(int startnodetag, double x1, double y1, double hx, double hy, double angle,
                   int nx, int ny, int ndf, const Vector& fix, const Vector& vel,
                   const Vector& mass, const Vector& boundary, 
                   Domain* theDomain, int& endnodetag); // rectangle

    int discretize(int startnodetag, double x1, double y1, double h, double angle, 
                   int num, int ndf, const Vector& fix, const Vector& vel,
                   const Vector& mass, const Vector& boundary,
                   Domain* theDomain, int& endnodetag); // line

    int discretize(int startnode, double x1, double y1, double x2, double y2, double x3, double y3,
                   int ni, int nj, int ndf, const Vector& fix, const Vector& vel,
                   const Vector& mass, const Vector& boundary,
                   Domain* theDomain, int& endnode);    // triangle

    int discretize(int startnode, double xc, double yc, double r1, double r2,
                   int nc, int nr, int ndf, const Vector& fix, const Vector& vel,
                   const Vector& mass, const Vector& boundary,
                   Domain* theDomain, int& endnode);    // circle
    int discretize(int startnode, char type, int n,
                   int nth, int nthfloor,  int ndf,
                   const Vector& fix, const Vector& vel, 
                   const Vector& mass, Domain* theDomain, int& endnode);   // frame

    int addPC(const ID& nodes, int startpnode, Domain* theDomain, int& endpnode); // add PC for nodes

    // triangulation
    int doTriangulation(const ID& nodes, double alpha, 
                        const ID& addnodes, Domain* theDomain, ID& eles);
    int doTriangulation(int startele, const ID& nodes, double alpha, 
                        const ID& addnodes, Domain* theDomain,
                        double rho, double mu, double b1, double b2, double thk);
    int doTriangulation(int startele, const ID& nodes, double alpha, 
                        const ID& addnodes, Domain* theDomain,
                        double t, const char* type, int matTag,
                        double p, double rho, double b1, double b2);

    // save
    int save(const char* name, const ID& snode, Domain* theDomain);

    // set boundary of fluids
    void setBoundary(double x1, double y1, double x2, double y2);
    void removeOutBoundNodes(const ID& nodes, Domain* theDomain);

    // set frame
    void setFrame(double x1, double y1, const Vector& span, 
                  const Vector& height);

    // calculate lift, drag, overturning moment from pressure
    Vector calculateForces(const ID& boundary, int basenode, 
                           Vector& dragdir, Vector& liftdir, 
                           Domain* theDomain);

    double geth() {return avesize;}
    void setNodes(const ID& nodes, bool fluid, bool append = false);

private:

    // initialize triangulateio
    void initializeTri(triangulateio& tri);

    // free triangulateio
    void freeTri(triangulateio& tri);
    void freeTriOut(triangulateio& tri);
    
    double PI;
    Vector bound;
    Vector frameBase;
    Vector Lspan;
    Vector Height;
    double avesize;
    ID fluidNodes;
    ID structureNodes;
};




#endif
