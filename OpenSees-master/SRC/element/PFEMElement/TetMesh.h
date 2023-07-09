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
//
// Description: The class TetMesh is for meshing tetrahedron
//

#ifndef TetMesh_h
#define TetMesh_h

#include "Mesh.h"

class TetMesh : public Mesh
{
public:
    explicit TetMesh(int tag);
    ~TetMesh();

    virtual const ID& getMeshTags() const {return mtags;}
    virtual void setMeshTags(const ID& tags) {mtags=tags;}

    int mesh();

    // remesh all
    static int remesh(double alpha);

private:
    ID mtags;
};

#endif
