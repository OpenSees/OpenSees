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

// $Revision: 1.3 $
// $Date: 2005-12-15 00:30:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CrdTransf3d.h,v $


// File: ~/CrdTransf/CrdTransf3d.h
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the class definition for 
// CrdTransf3d. CrdTransf3d provides the abstraction of spatial 
// coordinate transformation for a 3d frame. 

//
// What: "@(#) CrdTransf3d.h, revA"

#ifndef CrdTransf3d_h
#define CrdTransf3d_h
#include <CrdTransf.h>

class Vector;

class CrdTransf3d: public CrdTransf
{
public:
    CrdTransf3d (int tag, int classTag);
    virtual ~CrdTransf3d();
    
    virtual CrdTransf3d *getCopy(void) = 0;
    virtual int  getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis) = 0; 
    
protected:
    
private:
};

#endif
