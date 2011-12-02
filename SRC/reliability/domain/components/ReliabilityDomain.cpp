/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.4 $
// $Date: 2003-03-04 00:44:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ReliabilityDomain.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ReliabilityDomain.h>

#include <CorrelationCoefficient.h>
#include <RandomVariable.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>
#include <ParameterPositioner.h>
#include <ArrayOfTaggedObjects.h>
#include <ModulatingFunction.h>
#include <Filter.h>
#include <Spectrum.h>


ReliabilityDomain::ReliabilityDomain()
{
	theRandomVariablesPtr = new ArrayOfTaggedObjects (256);
	theCorrelationCoefficientsPtr = new ArrayOfTaggedObjects (256);
	theLimitStateFunctionsPtr = new ArrayOfTaggedObjects (256);
	theRandomVariablePositionersPtr = new ArrayOfTaggedObjects (256);
	theParameterPositionersPtr = new ArrayOfTaggedObjects (256);
	theModulatingFunctionsPtr = new ArrayOfTaggedObjects (256);
	theFiltersPtr = new ArrayOfTaggedObjects (256);
	theSpectraPtr = new ArrayOfTaggedObjects (256);
	tagOfActiveLimitStateFunction = 1;
}

ReliabilityDomain::~ReliabilityDomain()
{
	if (!theRandomVariablesPtr)
		delete theRandomVariablesPtr;
	if (!theCorrelationCoefficientsPtr)
		delete theCorrelationCoefficientsPtr;
	if (!theLimitStateFunctionsPtr)
		delete theLimitStateFunctionsPtr;
	if (!theRandomVariablePositionersPtr)
		delete theRandomVariablePositionersPtr;
	if (!theParameterPositionersPtr)
		delete theParameterPositionersPtr;
	if (!theModulatingFunctionsPtr)
		delete theModulatingFunctionsPtr;
	if (!theSpectraPtr)
		delete theSpectraPtr;
	if (!theFiltersPtr)
		delete theFiltersPtr;
}


bool
ReliabilityDomain::addRandomVariable(RandomVariable *theRandomVariable)
{
	bool result = theRandomVariablesPtr->addComponent(theRandomVariable);
	return result;
}

bool
ReliabilityDomain::addCorrelationCoefficient(CorrelationCoefficient *theCorrelationCoefficient)
{
	bool result = theCorrelationCoefficientsPtr->addComponent(theCorrelationCoefficient);
	return result;
}

bool
ReliabilityDomain::addLimitStateFunction(LimitStateFunction *theLimitStateFunction)
{
	bool result = theLimitStateFunctionsPtr->addComponent(theLimitStateFunction);
	return result;
}

bool
ReliabilityDomain::addRandomVariablePositioner(RandomVariablePositioner *theRandomVariablePositioner)
{
	bool result = theRandomVariablePositionersPtr->addComponent(theRandomVariablePositioner);
	return result;
}

bool
ReliabilityDomain::addParameterPositioner(ParameterPositioner *theParameterPositioner)
{
	bool result = theParameterPositionersPtr->addComponent(theParameterPositioner);
	return result;
}

bool
ReliabilityDomain::addModulatingFunction(ModulatingFunction *theModulatingFunction)
{
	bool result = theModulatingFunctionsPtr->addComponent(theModulatingFunction);
	return result;
}

bool
ReliabilityDomain::addSpectrum(Spectrum *theSpectrum)
{
	bool result = theSpectraPtr->addComponent(theSpectrum);
	return result;
}

bool
ReliabilityDomain::addFilter(Filter *theFilter)
{
	bool result = theFiltersPtr->addComponent(theFilter);
	return result;
}




RandomVariable *
ReliabilityDomain::getRandomVariablePtr(int tag)
{
	TaggedObject *theComponent = theRandomVariablesPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	RandomVariable *result = (RandomVariable *) theComponent;
	return result;
}


CorrelationCoefficient * 
ReliabilityDomain::getCorrelationCoefficientPtr(int tag)
{
	TaggedObject *theComponent = theCorrelationCoefficientsPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	CorrelationCoefficient *result = (CorrelationCoefficient *) theComponent;
	return result;
}


LimitStateFunction *
ReliabilityDomain::getLimitStateFunctionPtr(int tag)
{
	TaggedObject *theComponent = theLimitStateFunctionsPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	LimitStateFunction *result = (LimitStateFunction *) theComponent;
	return result;
}


RandomVariablePositioner *
ReliabilityDomain::getRandomVariablePositionerPtr(int tag)
{
	TaggedObject *theComponent = theRandomVariablePositionersPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	RandomVariablePositioner *result = (RandomVariablePositioner *) theComponent;
	return result;
}

ParameterPositioner *
ReliabilityDomain::getParameterPositionerPtr(int tag)
{
	TaggedObject *theComponent = theParameterPositionersPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	ParameterPositioner *result = (ParameterPositioner *) theComponent;
	return result;
}


ModulatingFunction *
ReliabilityDomain::getModulatingFunction(int tag)
{
	TaggedObject *theComponent = theModulatingFunctionsPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	ModulatingFunction *result = (ModulatingFunction *) theComponent;
	return result;
}


Spectrum *
ReliabilityDomain::getSpectrum(int tag)
{
	TaggedObject *theComponent = theSpectraPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	Spectrum *result = (Spectrum *) theComponent;
	return result;
}


Filter *
ReliabilityDomain::getFilter(int tag)
{
	TaggedObject *theComponent = theFiltersPtr->getComponentPtr(tag);
//	if ( *theComponent == 0 )
//		return 0;
	Filter *result = (Filter *) theComponent;
	return result;
}


int
ReliabilityDomain::removeRandomVariablePositioner(int tag)
{
	theRandomVariablePositionersPtr->removeComponent(tag);

	return 0;
}

int
ReliabilityDomain::removePerformanceFunction(int tag)
{
	theLimitStateFunctionsPtr->removeComponent(tag);

	return 0;
}


int
ReliabilityDomain::getNumberOfRandomVariables()
{
	return theRandomVariablesPtr->getNumComponents();
}
int
ReliabilityDomain::getNumberOfCorrelationCoefficients()
{
	return theCorrelationCoefficientsPtr->getNumComponents();
}
int
ReliabilityDomain::getNumberOfLimitStateFunctions()
{
	return theLimitStateFunctionsPtr->getNumComponents();
}
int
ReliabilityDomain::getNumberOfRandomVariablePositioners()
{
	return theRandomVariablePositionersPtr->getNumComponents();
}
int
ReliabilityDomain::getNumberOfParameterPositioners()
{
	return theParameterPositionersPtr->getNumComponents();
}
int
ReliabilityDomain::getNumberOfModulatingFunctions()
{
	return theModulatingFunctionsPtr->getNumComponents();
}
int
ReliabilityDomain::getNumberOfFilters()
{
	return theFiltersPtr->getNumComponents();
}
int
ReliabilityDomain::getNumberOfSpectra()
{
	return theSpectraPtr->getNumComponents();
}





int
ReliabilityDomain::getTagOfActiveLimitStateFunction()
{
	return tagOfActiveLimitStateFunction;
}
void
ReliabilityDomain::setTagOfActiveLimitStateFunction(int passedTag)
{
	tagOfActiveLimitStateFunction = passedTag;
}

