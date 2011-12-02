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
// $Date: 2005-11-28 21:34:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/ConstraintHandler.cpp,v $
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the implementation of ConstraintHandler.
//
// What: "@(#) ConstraintHandler.h, revA"

#include <ConstraintHandler.h>
#include <Domain.h>
#include <AnalysisModel.h>
#include <Integrator.h>
#include <FE_EleIter.h>
#include <FE_Element.h>

ConstraintHandler::ConstraintHandler(int clasTag)
:MovableObject(clasTag),
 theDomainPtr(0),theAnalysisModelPtr(0),theIntegratorPtr(0)
{
}


ConstraintHandler::~ConstraintHandler()
{
    
}

int
ConstraintHandler::doneNumberingDOF(void)
{
  // iterate through the FE_Element getting them to set their IDs
  FE_EleIter &theEle = theAnalysisModelPtr->getFEs();
  FE_Element *elePtr;
  while ((elePtr = theEle()) != 0)
    elePtr->setID();
  return 0;
}

void 
ConstraintHandler::setLinks(Domain &theDomain, 
			    AnalysisModel &theModel,
			    Integrator &theIntegrator)
{
    theDomainPtr = &theDomain;
    theAnalysisModelPtr = &theModel;
    theIntegratorPtr = &theIntegrator;
}
	

int
ConstraintHandler::update(void)
{
  return 0;
}

int
ConstraintHandler::applyLoad(void)
{
  return 0;
}


Domain *
ConstraintHandler::getDomainPtr(void) const
{
    return theDomainPtr;
}

AnalysisModel *
ConstraintHandler::getAnalysisModelPtr(void) const
{
    return theAnalysisModelPtr;
}

Integrator *
ConstraintHandler::getIntegratorPtr(void) const
{
    return theIntegratorPtr;
}



