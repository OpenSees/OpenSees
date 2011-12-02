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
// $Date: 2000-09-15 08:23:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/groundMotion/GroundMotionIntegrator.cpp,v $
                                                                        
                                                                        
// File: ~/domain/groundMotion/GroundMotionIntegrator.C
// 
// Written: MHS
// Created: 10/99
// Revision: A
//
// Description: This file contains the class definition for 
// GroundMotionIntegrator.
//
// What: "@(#) GroundMotionIntegrator.C, revA"

#include <GroundMotionIntegrator.h>

GroundMotionIntegrator::GroundMotionIntegrator (int classTag, double del)
:MovableObject(classTag), delta(del)
{

}

GroundMotionIntegrator::~GroundMotionIntegrator()
{

}
