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
// $Date: 2003-02-14 23:02:10 $
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/TaggedObject.cpp,v $
                                                                        
                                                                        
// File: ~/taggedt/TaggedObject.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class implementation for TaggedObject.
// A TaggedObject is an object with an integer identifier. It is used as
// a base class by DomainComponent, Graph and other classes in the framework.
//
// What: "@(#) TaggedObject.C, revA"

#include "TaggedObject.h"

TaggedObject::TaggedObject(int tag)
  :theTag(tag)
{

}

TaggedObject::~TaggedObject()
{
    // does nothing
}
 
void
TaggedObject::setTag(int newTag) 
{
    theTag = newTag;
}

OPS_Stream &operator<<(OPS_Stream &s, TaggedObject &m)
{
    m.Print(s);
    return s;
}
