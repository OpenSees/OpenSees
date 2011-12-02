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
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/PerformanceFunctionCoefficientIter.cpp,v $
                                                                        
#include <PerformanceFunctionCoefficientIter.h>
#include <PerformanceFunctionCoeff.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


// SingleDomNodIter(SingleDomain &theDomain):
//	constructor that takes the model, just the basic iter

PerformanceFunctionCoefficientIter::PerformanceFunctionCoefficientIter
(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}

PerformanceFunctionCoefficientIter::~PerformanceFunctionCoefficientIter()
{
}    


void
PerformanceFunctionCoefficientIter::reset(void)
{
    myIter.reset();
}


PerformanceFunctionCoeff *
PerformanceFunctionCoefficientIter::operator()(void)
{
    // check if we still have Nodes in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
	return 0;
    else {
	PerformanceFunctionCoeff *result = 
		(PerformanceFunctionCoeff *)theComponent;
	return result;
    }
}


    
    
