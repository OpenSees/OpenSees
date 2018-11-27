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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Load.cpp,v $
                                                                        
                                                                        
// File: ~/domain/load/Load.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementaion for the Load class.
//
// What: "@(#) Load.C, revA"


#include <Load.h>


Load::Load(int tag, int clasTag)
:DomainComponent(tag, clasTag), loadPatternTag(-1)
{
    // does nothing
}


Load::~Load()
{
    // does nothing
}

void
Load::setLoadPatternTag(int tag)
{
  loadPatternTag = tag;
}

int
Load::getLoadPatternTag(void) const
{
  return loadPatternTag;
}
