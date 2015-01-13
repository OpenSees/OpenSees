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
#include <map>
#include <vector>

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

    int discretize(int startnodetag, const Vector& points, 
                   int ndf, const Vector& fix, const Vector& vel,
                   const Vector& mass, Domain* theDomain, int& endnodetag); // particles

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
                   const Vector& mass, const Vector& boundary, const Vector& angle,
                   Domain* theDomain, int& endnode);    // circle
    int discretize(int startnode, double x, double y,
                   double hcol, double hbeam, int ndf,
                   const Vector& spans, const Vector& heights,
                   const Vector& fix, const Vector& colvel, 
                   const Vector& beamvel, const Vector& colmass, 
                   const Vector& beammass, const Vector& roofmass,
                   Domain* theDomain, ID& nodelist);   // frame

    int addPC(const ID& nodes, int pndf, int startpnode, Domain* theDomain, int& endpnode); // add PC for nodes

    // triangulation
    // core function
    int doTriangulation(double alpha, const ID& groups, const ID& addgroups,
                        Domain* theDomain, ID& eles);  
    // linear elements
    int doTriangulation(int startele, double alpha, const ID& groups, 
                        const ID& addgroups,Domain* theDomain,
                        double rho, double mu, double b1, double b2, 
                        double thk, double kappa, int type, int& endele);
    // solid element
    int doTriangulation(int startele, double alpha, const ID& groups, 
                        const ID& addgroups, Domain* theDomain,
                        double t, const char* type, int matTag,
                        double p, double rho, double b1, double b2, 
                        int& endele);

    // Crouzeix-Raviart element
    int doTriangulation(int newNodeRegTag, int eleRegTag,
                        double alpha, const ID& groups, const ID& addgroups,
                        Domain* theDomain, double rho, double mu,
                        double b1, double b2, double thk, double kappa);

    // save
    int save(const char* name, Domain* theDomain, int maxelenodes=3);
    int save(const char* name, const ID& nodeRegs,
             const ID& eleRegs, Domain* theDomain, int maxelenodes=3);

    // set boundary of fluids
    void removeOutBoundNodes(const ID& groups, double x1, double y1,
                             double x2, double y2, Domain* theDomain);

    // calculate lift, drag, overturning moment from pressure
    Vector calculateForces(const ID& boundary, int basenode, 
                           Vector& dragdir, Vector& liftdir, 
                           Domain* theDomain);

    // set region states
    int updateNode(int tag, int dof, double vale, 
                   int type, Domain* theDomain);

    void setNodes(const ID& nodes, int region, bool series, 
                  int action, Domain* theDomain);
    void getNodes(const ID& regions, std::map<int,int>& nodes, 
                  Domain* theDomain);
    void setElements(const ID& elements, int region, bool series, 
                     int action, Domain* theDomain);
    void getElements(const ID& regions, std::map<int,int>& elements, 
                     Domain* theDomain);
    void removeElements(int regTag, Domain* theDomain);

    // identify interface
    void identify(double g, Domain* theDomain);

    // find a node tag
    int findNodeTag(Domain* theDomain);
    int findEleTag(Domain* theDomain);

private:

    // initialize triangulateio
    void initializeTri(triangulateio& tri);

    // free triangulateio
    void freeTri(triangulateio& tri);
    void freeTriOut(triangulateio& tri);
    
    // PI
    static double PI;
};




#endif
