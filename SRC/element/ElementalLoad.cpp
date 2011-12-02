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
// $Source: /usr/local/cvs/OpenSees/SRC/element/ElementalLoad.cpp,v $
                                                                        
                                                                        
// File: ~/element/ElementalLoad.h
//
// Written: fmk 11/95
// Revised:
//
// Purpose: This file contains the methods for class ElementalLoad.
//
// The interface:

#include "ElementalLoad.h"

ElementalLoad::ElementalLoad(int elementTag, int tag, int cTag)
:Load(tag, cTag), theElementTag(elementTag)
{
}

// provided for the FEM_Object broker; the tag and elementTag need
// to be supplied in recvSelf();
ElementalLoad::ElementalLoad(int cTag)
:Load(0, cTag), theElementTag(0)
{
}

ElementalLoad::~ElementalLoad()
{
    ;
}


int 
ElementalLoad::getElementTag(void) const
{
    return theElementTag;
}


