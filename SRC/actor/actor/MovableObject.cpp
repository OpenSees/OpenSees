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
                                                                        
// $Revision: 1.2 $
// $Date: 2006-08-03 23:12:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/actor/actor/MovableObject.cpp,v $
                                                                        
                                                                        
// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of MovableObject.
//
// What: "@(#) MovableObject.C, revA"

#include <MovableObject.h>

MovableObject::MovableObject(int cTag, int dTag)
:classTag(cTag), dbTag(dTag)
{

}


MovableObject::MovableObject(int theTag)
:classTag(theTag), dbTag(0)
{

}

MovableObject::~MovableObject()
{

}




int 
MovableObject::getClassTag(void) const
{
    return classTag;
}

const char *
MovableObject::getClassType(void) const
{
    return "UnknownMovableObject";
}


int 
MovableObject::getDbTag(void) const
{
    return dbTag;
}

void
MovableObject::setDbTag(int newTag)
{
  dbTag = newTag;
}
