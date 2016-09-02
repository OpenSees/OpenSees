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
// $Date: 2015-10-19 $
                                                                        
// Written: Minjie Zhu
// Created: Oct. 19
//
// Description: The class TriMesh is for meshing triangles
//

#ifndef TriMesh_h
#define TriMesh_h

class Domain;
class ID;

class TriMesh
{
public:
    TriMesh(Domain& theDomain, int ndf);
    ~TriMesh();

    virtual int mesh(int rtag, double size, const ID& nodes,const ID& bound);
    virtual int mesh(double alpha, const ID& rtagsfree, const ID& rtagsfix);
    
private:
    Domain* theDomain;
    int ndf;
};

#endif
