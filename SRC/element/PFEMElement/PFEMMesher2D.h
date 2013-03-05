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

class ID;
class Domain;


class PFEMMesher2D 
{

public:

    // constructor
    PFEMMesher2D();

    // destructor
    ~PFEMMesher2D();

    // discretize domain
    int discretize(int startnodetag, const Vector& points, const ID& segments, const Vector& holes,
                   double maxarea, int ndf, const ID& fix, 
                   const Vector& mass, Domain* theDomain);

    int discretize(int startnodetag, double x1, double y1, double hx, double hy, double slope,
                   int nx, int ny, int ndf, const ID& fix, 
                   const Vector& mass, Domain* theDomain);

    int discretize(int startnodetag, double x1, double y1, double h, double slope, 
                   int num, int ndf, const ID& fix,
                   const Vector& mass, Domain* theDomain);

    // triangulation
    int doTriangulation(int startnodetag, int endnodetag, double alpha, 
                        int startanodetag, int endanodetag, Domain* theDomain, ID& eles);
    int doTriangulation(int starteletag, int startnodetag, int endnodetag, double alpha, 
                        int startanodetag, int endanodetag, Domain* theDomain,
                        double rho, double mu, double b1, double b2);

    // save
    int save(const char* name, int startsnode, int endsnode, Domain* theDomain);

    // set boundary of fluids
    void setBoundary(double x1, double y1, double x2, double y2);
    void removeOutBoundNodes(int startnode, int endnode, Domain* theDomain);

private:

    // initialize triangulateio
    void initializeTri(triangulateio& tri);

    // free triangulateio
    void freeTri(triangulateio& tri);
    void freeTriOut(triangulateio& tri);
    
    double PI;
    Vector bound;
};




#endif
