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
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/repres/section/SectionRepres.cpp,v $
                                                                        
                                                                        
// File: SectionRepres.C
//
// Written by Remo M. de Souza
// November 1998

#include <SectionRepres.h>
#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>

static MapOfTaggedObjects theSectionRepresObjects;

bool OPS_addSectionRepres(SectionRepres *newComponent)
{
    return theSectionRepresObjects.addComponent(newComponent);
}

SectionRepres *OPS_getSectionRepres(int tag)
{
    TaggedObject *theResult = theSectionRepresObjects.getComponentPtr(tag);
    if(theResult == 0) {
	return 0;
    }
    SectionRepres *theRep = (SectionRepres *)theResult;

    return theRep;
}

void OPS_clearAllSectionRepres(void)
{
    theSectionRepresObjects.clearAll();
}

SectionRepres::SectionRepres(int tag):
                 TaggedObject(tag)
{
   
}
   
SectionRepres::~SectionRepres(void)
{

}













