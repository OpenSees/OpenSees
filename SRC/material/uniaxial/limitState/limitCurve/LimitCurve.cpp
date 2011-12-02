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
// $Date: 2006-02-07 23:15:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/limitState/limitCurve/LimitCurve.cpp,v $

// Written: kje 
// Created: 08/01
// Revision: A
//
// Description: This file contains the class implementation for 
// LimitCurve.
//
// What: "@(#) LimitCurve.C, revA"

#include <LimitCurve.h>
#include <string.h>
#include <Information.h>


LimitCurve::LimitCurve(int tag, int clasTag)
:TaggedObject(tag), MovableObject(clasTag)
{

}

LimitCurve::~LimitCurve()
{
  // does nothing
}




int
LimitCurve::setParameter(const char **argv, int argc, Information &info)
{
    return -1;
}

int
LimitCurve::updateParameter(int parameterID, Information &info)
{
    return -1;
}
