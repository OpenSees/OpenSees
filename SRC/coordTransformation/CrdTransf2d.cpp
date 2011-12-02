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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CrdTransf2d.cpp,v $
                                                                        
                                                                        
// File: ~/crdTransf/CrdTransf2d.C
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the implementation for the CrdTransf2d class.
// CrdTransf2d provides the abstraction of a 2d frame 
// coordinate transformation. It is an abstract base class and 
// thus no objects of  it's type can be instatiated. 

#include <CrdTransf2d.h>


// constructor:
CrdTransf2d::CrdTransf2d(int tag, int classTag): CrdTransf(tag, classTag)
{
}


// destructor:
CrdTransf2d::~CrdTransf2d()
{
}

