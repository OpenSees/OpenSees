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
                                                                        
// $Revision: 1.1 $
// $Date: 2006-11-08 20:06:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/utility/StringContainer.h,v $
                                                                        
// Written: Jade Cohen
// Created: 19/05
//
// Description: This file contains the method for GPYS2d.
// GPYS2d is used to build general polynomial yield surface for 2d
// frame elements.
//

#ifndef GPYS2d_h
#define GPYS2d_h

#include <Matrix.h>
#include <Vector.h>

void GPYS2d(const Matrix &GPYSC, const Vector &xyref, const Vector &ScVec, double &f, Vector &g, Matrix &h);

#endif
