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
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/CrdTransf.cpp,v $


// File: ~/crdTransf/CrdTransf.C
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
//
// Description: This file contains the implementation for the CrdTransf class.
// CrdTransf provides the abstraction of a frame 
// coordinate transformation. It is an abstract base class and 
// thus no objects of  it's type can be instatiated. 

#include <CrdTransf.h>
#include <Vector.h>


// constructor:
CrdTransf::CrdTransf(int tag, int classTag):TaggedObject(tag), MovableObject(classTag)
{
}

// destructor:
CrdTransf::~CrdTransf()
{
}

const Vector &
CrdTransf::getBasicDisplSensitivity(int gradNumber)
{
    opserr << "WARNING CrdTransf::getBasicDisplSensitivity() - this method "
        << " should not be called." << endln;
    
    static Vector dummy(1);
    return dummy;
}

const Vector &
CrdTransf::getGlobalResistingForceShapeSensitivity(const Vector &pb, const Vector &p0)
{
    opserr << "ERROR CrdTransf::getGlobalResistingForceSensitivity() - has not been"
        << " implemented yet for the chosen transformation." << endln;
    
    static Vector dummy(1);
    return dummy;
}

const Vector &
CrdTransf::getBasicTrialDispShapeSensitivity(void)
{
    opserr << "ERROR CrdTransf::getBasicTrialDispShapeSensitivity() - has not been"
        << " implemented yet for the chosen transformation." << endln;
    
    static Vector dummy(1);
    return dummy;
}