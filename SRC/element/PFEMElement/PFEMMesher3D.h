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

class ID;
class Domain;


class PFEMMesher3D 
{

public:

    // constructor
    PFEMMesher3D();

    // destructor
    ~PFEMMesher3D();

    // discretize domain
    // int discretize(const Vector& points, const ID& segments, 
    //                       double maxarea, int ndf, const ID& fix, 
    //                       const Vector& mass, Domain* theDomain, 
    //                       ID& nodes);

    // int discretize(double x1, double y1, double hx, double hy, double slope,
    //                       int nx, int ny, int ndf, const ID& fix, 
    //                       const Vector& mass, Domain* theDomain, 
    //                       ID& nodes);

    // int discretize(double x1, double y1, double h, double slope, 
    //                       int num, int ndf, const ID& fix, 
    //                       const Vector& mass, Domain* theDomain, 
    //                       ID& nodes);

    // triangulation
    int doTriangulation(int startnode, int endnode, double alpha, 
                        int startanode, int endanode, Domain* theDomain, 
                        ID& eles);
    // int doTriangulation(const ID& nodes, double alpha, 
    //                            const ID& addNodes, Domain* theDomain, 
    //                            ID& eles, double rho, double mu,
    //                            double b1, double b2);

    // save
    // int save(const char* name, const ID& nodes, Domain* theDomain);

    // set boundary of fluids
    // void setBoundary(double x1, double y1, double x2, double y2);
    // void removeOutBoundNodes(const ID& nodes, ID& nodes2, Domain* theDomain);

private:

    double PI;
    Vector bound;

};




#endif
