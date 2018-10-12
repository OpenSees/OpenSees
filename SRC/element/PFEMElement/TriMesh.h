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

#include "Mesh.h"

class TriMesh : public Mesh
{
public:
    explicit TriMesh(int tag);
    ~TriMesh();

    virtual const ID& getLineTags() const {return ltags;}
    virtual void setLineTags(const ID& tags) {ltags=tags;}

    int mesh();

    // remesh all
    static int remesh(double alpha);

private:
    ID ltags;
};

#endif
