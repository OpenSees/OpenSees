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
// Description: The class LineMesh is for meshing lines.
//

#ifndef LineMesh_h
#define LineMesh_h

class Domain;
class ID;

class LineMesh
{
public:
    LineMesh(Domain& theDomain, int ndm, int ndf);
    ~LineMesh();

    virtual int mesh(int rtag, double size, const ID& nodes,const ID& bound);
    
private:
    virtual int mesh(double size, int nd1, int nd2, ID& nodes, ID& elenodes);

    Domain* theDomain;
    int ndm, ndf;
};

#endif
