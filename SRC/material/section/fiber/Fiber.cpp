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
                                                                        
// $Revision: 1.4 $
// $Date: 2006-08-04 18:32:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/fiber/Fiber.cpp,v $
                                                                        
                                                                        
// File: ~/fiber/Fiber.C
//
// Written: Remo Magalhaes de Souza
// Created: 10/98
// Revision: 
//
// Description: This file contains the implementation for the
// Fiber class. Fiber provides the abstraction of a section fiber.
// 
// What: "@(#) Fiber.C, revA"


#include <Fiber.h>

// constructor:

Fiber::Fiber(int tag, int classTag):
TaggedObject(tag), MovableObject(classTag)
{

}

// destructor:
Fiber::~Fiber()
{

}

Response*
Fiber::setResponse(const char **argv, int argc, Information &info, OPS_Stream &s)
{
	return 0;
}

int
Fiber::getResponse(int responseID, Information &info)
{
	return -1;
}

