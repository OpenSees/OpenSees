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
                                                                        
// $Revision: 1.11 $
// $Date: 2007-10-26 17:37:36 $
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

#include <RandomVariableIter.h>
#include <RandomVariablePositionerIter.h>
#include <ParameterPositionerIter.h>
#include <LimitStateFunctionIter.h>
#include <CorrelationCoefficientIter.h>

ReliabilityDomain::ReliabilityDomain():
  numRandomVariables(0), numRandomVariablePositioners(0),
  numParameterPositioners(0), numLimitStateFunctions(0),
  numCorrelationCoefficients(0)
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

	theRVIter = new RandomVariableIter(theRandomVariablesPtr);
	theRVPosIter = new RandomVariablePositionerIter(theRandomVariablePositionersPtr);
	theParamPosIter = new ParameterPositionerIter(theParameterPositionersPtr);
	theLSFIter = new LimitStateFunctionIter(theLimitStateFunctionsPtr);
	theCCIter = new CorrelationCoefficientIter(theCorrelationCoefficientsPtr);
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

	if (theRVIter != 0)
	  delete theRVIter;
	if (theRVPosIter != 0)
	  delete theRVPosIter;
	if (theParamPosIter != 0)
	  delete theParamPosIter;
	if (theLSFIter != 0)
	  delete theLSFIter;
	if (theCCIter != 0)
	  delete theCCIter;
}


bool
ReliabilityDomain::addRandomVariable(RandomVariable *theRandomVariable)
{
  theRandomVariable->setIndex(numRandomVariables);
  numRandomVariables++;

  bool result = theRandomVariablesPtr->addComponent(theRandomVariable);
  return result;
}

bool
ReliabilityDomain::addCorrelationCoefficient(CorrelationCoefficient *theCorrelationCoefficient)
{
  theCorrelationCoefficient->setIndex(numCorrelationCoefficients);
  numCorrelationCoefficients++;

  bool result = theCorrelationCoefficientsPtr->addComponent(theCorrelationCoefficient);
  return result;
}

bool
ReliabilityDomain::addLimitStateFunction(LimitStateFunction *theLimitStateFunction)
{
  theLimitStateFunction->setIndex(numLimitStateFunctions);
  numLimitStateFunctions++;

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



RandomVariableIter &
ReliabilityDomain::getRandomVariables(void)
{
  theRVIter->reset();
  return *theRVIter;
}

RandomVariablePositionerIter &
ReliabilityDomain::getRandomVariablePositioners(void)
{
  theRVPosIter->reset();
  return *theRVPosIter;
}

ParameterPositionerIter &
ReliabilityDomain::getParameterPositioners(void)
{
  theParamPosIter->reset();
  return *theParamPosIter;
}

LimitStateFunctionIter &
ReliabilityDomain::getLimitStateFunctions(void)
{
  theLSFIter->reset();
  return *theLSFIter;
}

CorrelationCoefficientIter &
ReliabilityDomain::getCorrelationCoefficients(void)
{
  theCCIter->reset();
  return *theCCIter;
}

RandomVariable *
ReliabilityDomain::getRandomVariablePtr(int tag)
{
	TaggedObject *theComponent = theRandomVariablesPtr->getComponentPtr(tag);
	if ( theComponent == 0 )
		return 0;
	RandomVariable *result = (RandomVariable *) theComponent;
	return result;
}


CorrelationCoefficient * 
ReliabilityDomain::getCorrelationCoefficientPtr(int tag)
{
	TaggedObject *theComponent = theCorrelationCoefficientsPtr->getComponentPtr(tag);
	if ( theComponent == 0 )
		return 0;
	CorrelationCoefficient *result = (CorrelationCoefficient *) theComponent;
	return result;
}


LimitStateFunction *
ReliabilityDomain::getLimitStateFunctionPtr(int tag)
{
  TaggedObject *theComponent = theLimitStateFunctionsPtr->getComponentPtr(tag);
  if (theComponent == 0)
    return 0;
  LimitStateFunction *result = (LimitStateFunction *) theComponent;
  return result;
}


RandomVariablePositioner *
ReliabilityDomain::getRandomVariablePositionerPtr(int tag)
{
	TaggedObject *theComponent = theRandomVariablePositionersPtr->getComponentPtr(tag);
	if ( theComponent == 0 )
		return 0;
	RandomVariablePositioner *result = (RandomVariablePositioner *) theComponent;
	return result;
}

ParameterPositioner *
ReliabilityDomain::getParameterPositionerPtr(int tag)
{
	TaggedObject *theComponent = theParameterPositionersPtr->getComponentPtr(tag);
	if ( theComponent == 0 )
		return 0;
	ParameterPositioner *result = (ParameterPositioner *) theComponent;
	return result;
}


ModulatingFunction *
ReliabilityDomain::getModulatingFunction(int tag)
{
	TaggedObject *theComponent = theModulatingFunctionsPtr->getComponentPtr(tag);
	if ( theComponent == 0 )
		return 0;
	ModulatingFunction *result = (ModulatingFunction *) theComponent;
	return result;
}


Spectrum *
ReliabilityDomain::getSpectrum(int tag)
{
	TaggedObject *theComponent = theSpectraPtr->getComponentPtr(tag);
	if ( theComponent == 0 )
		return 0;
	Spectrum *result = (Spectrum *) theComponent;
	return result;
}


Filter *
ReliabilityDomain::getFilter(int tag)
{
	TaggedObject *theComponent = theFiltersPtr->getComponentPtr(tag);
	if ( theComponent == 0 )
		return 0;
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
ReliabilityDomain::removeParameterPositioner(int tag)
{
	theParameterPositionersPtr->removeComponent(tag);

	return 0;
}

int
ReliabilityDomain::removeAllParameterPositioners(void)
{
	theParameterPositionersPtr->clearAll();

	return 0;
}

int
ReliabilityDomain::removeRandomVariable(int tag)
{
  RandomVariable *theRV = (RandomVariable*) theRandomVariablesPtr->getComponentPtr(tag);
  int indexToRemove = theRV->getIndex();

  // shift indices down by one
  for (int i = indexToRemove; i < numRandomVariables-1; i++) {
    // find rv with index equal to i+1
    theRVIter->reset();
    while ((theRV = (*theRVIter)()) != 0) {
      if (theRV->getIndex() == i+1)
	break;
    }
    // now set its index to i
    theRV->setIndex(i);
  }

  // Now remove the component
  theRandomVariablesPtr->removeComponent(tag);
  numRandomVariables--;

  return 0;
}


int
ReliabilityDomain::removeCorrelationCoefficient(int tag)
{
  CorrelationCoefficient *theCC = (CorrelationCoefficient*) theCorrelationCoefficientsPtr->getComponentPtr(tag);
  int indexToRemove = theCC->getIndex();

  // shift indices down by one
  for (int i = indexToRemove; i < numCorrelationCoefficients-1; i++) {
    // find rv with index equal to i+1
    theCCIter->reset();
    while ((theCC = (*theCCIter)()) != 0) {
      if (theCC->getIndex() == i+1)
	break;
    }
    // now set its index to i
    theCC->setIndex(i);
  }

  // Now remove the component
  theCorrelationCoefficientsPtr->removeComponent(tag);
  numCorrelationCoefficients--;

  return 0;
}

int
ReliabilityDomain::removePerformanceFunction(int tag)
{
  LimitStateFunction *theLSF = (LimitStateFunction*) theLimitStateFunctionsPtr->getComponentPtr(tag);
  int indexToRemove = theLSF->getIndex();

  // shift indices down by one
  for (int i = indexToRemove; i < numLimitStateFunctions-1; i++) {
    // find rv with index equal to i+1
    theLSFIter->reset();
    while ((theLSF = (*theLSFIter)()) != 0) {
      if (theLSF->getIndex() == i+1)
	break;
    }
    // now set its index to i
    theLSF->setIndex(i);
  }

  // Now remove the component
  theLimitStateFunctionsPtr->removeComponent(tag);
  numLimitStateFunctions--;

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

