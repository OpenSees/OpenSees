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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/single/SingleDomAllSP_Iter.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/domain/domain/SingleDomAllSP_Iter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// SingleDomAllSP_Iter. SingleDomAllSP_Iter is a class for iterating through the 
// elements of a domain. 

#include "SingleDomAllSP_Iter.h"

#include <Domain.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <SP_Constraint.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>

// SingleDomAllSP_Iter(SingleDomAllain &theDomain):
//	constructor that takes the model, just the basic iter

SingleDomAllSP_Iter::SingleDomAllSP_Iter(Domain &domain)
  :theDomain(&domain), doneDomainSPs(false)
{

}

SingleDomAllSP_Iter::~SingleDomAllSP_Iter()
{
}    


void
SingleDomAllSP_Iter::reset(void)
{
  theDomainSPs = &(theDomain->getSPs());
  theLoadPatterns = &(theDomain->getLoadPatterns());
  currentLoadPattern = (*theLoadPatterns)();
  if (currentLoadPattern != 0) {
      theLoadPatternSPs = &(currentLoadPattern->getSPs());
  }

  doneDomainSPs = false;
}


SP_Constraint *
SingleDomAllSP_Iter::operator()(void)
{
  SP_Constraint *theRes = 0;

  if (doneDomainSPs == false) {
    theRes = (*theDomainSPs)();
    if (theRes != 0)
      return theRes;
    else
      doneDomainSPs = true;
  }

  while (currentLoadPattern != 0) {
    theRes = (*theLoadPatternSPs)();
    if (theRes == 0) {
      currentLoadPattern = (*theLoadPatterns)();
      if (currentLoadPattern != 0)
	theLoadPatternSPs = &(currentLoadPattern->getSPs());
    } else
	return theRes;
  }

  return 0;
}
    
    


    
    
