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
                                                                        
// $Revision: 1.15 $
// $Date: 2008-04-10 16:25:09 $
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
  numRandomVariables(0), numLimitStateFunctions(0),
  rvIndex(0), lsfIndex(0)
{
	theRandomVariablesPtr = new ArrayOfTaggedObjects (256);
	theCorrelationCoefficientsPtr = new ArrayOfTaggedObjects (256);
	theLimitStateFunctionsPtr = new ArrayOfTaggedObjects (256);
	theRandomVariablePositionersPtr = new ArrayOfTaggedObjects (256);
	theParameterPositionersPtr = new ArrayOfTaggedObjects (256);
	theModulatingFunctionsPtr = new ArrayOfTaggedObjects (256);
	theFiltersPtr = new ArrayOfTaggedObjects (256);
	theSpectraPtr = new ArrayOfTaggedObjects (256);

	theDesignVariablesPtr = new ArrayOfTaggedObjects (256);
	theDesignVariablePositionersPtr = new ArrayOfTaggedObjects (256);
	theObjectiveFunctionsPtr = new ArrayOfTaggedObjects (256);
	theConstraintFunctionsPtr = new ArrayOfTaggedObjects (256);

	tagOfActiveLimitStateFunction = 1;

	theRVIter = new RandomVariableIter(theRandomVariablesPtr);
	theRVPosIter = new RandomVariablePositionerIter(theRandomVariablePositionersPtr);
	theParamPosIter = new ParameterPositionerIter(theParameterPositionersPtr);
	theLSFIter = new LimitStateFunctionIter(theLimitStateFunctionsPtr);
	theCCIter = new CorrelationCoefficientIter(theCorrelationCoefficientsPtr);

	rvIndex = new int[rvSize_init];
	rvSize = rvSize_init;

	lsfIndex = new int[lsfSize_init];
	lsfSize = lsfSize_init;
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

	if (!theDesignVariablesPtr)
		delete theDesignVariablesPtr;
	if (!theDesignVariablePositionersPtr)
		delete theDesignVariablePositionersPtr;
	if (!theConstraintFunctionsPtr)
		delete theConstraintFunctionsPtr;
	if (!theObjectiveFunctionsPtr)
		delete theObjectiveFunctionsPtr;

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

	if (rvIndex != 0)
	  delete [] rvIndex;
	if (lsfIndex != 0)
	  delete [] lsfIndex;
}


bool
ReliabilityDomain::addRandomVariable(RandomVariable *theRandomVariable)
{
  bool result = theRandomVariablesPtr->addComponent(theRandomVariable);

  if (result == true) {

    // Array is full
    if (numRandomVariables == rvSize) {

      // Increase size and allocate new array
      rvSize += rvSize_grow;
      int *tmp_rvIndex = new int[rvSize];

      // Copy values from old array to new
      for (int i = 0; i < numRandomVariables; i++)
	tmp_rvIndex[i] = rvIndex[i];

      // Get rid of old array
      delete [] rvIndex;

      // Set pointer to new array
      rvIndex = tmp_rvIndex;
    }

    // Add to index
    rvIndex[numRandomVariables] = theRandomVariable->getTag();
    numRandomVariables++;
  }

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

  if (result == true) {

    // Array is full
    if (numLimitStateFunctions == lsfSize) {

      // Increase size and allocate new array
      lsfSize += lsfSize_grow;
      int *tmp_lsfIndex = new int[lsfSize];

      // Copy values from old array to new
      for (int i = 0; i < numLimitStateFunctions; i++)
	tmp_lsfIndex[i] = lsfIndex[i];

      // Get rid of old array
      delete [] lsfIndex;

      // Set pointer to new array
      lsfIndex = tmp_lsfIndex;
    }

    // Add to index
    lsfIndex[numLimitStateFunctions] = theLimitStateFunction->getTag();
    numLimitStateFunctions++;
  }

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



//- Quan
bool
ReliabilityDomain::addDesignVariable(DesignVariable *theDesignVariable)
{
	bool result = theDesignVariablesPtr->addComponent(theDesignVariable);
	return result;
}


bool
ReliabilityDomain::addDesignVariablePositioner(DesignVariablePositioner *theDesignVariablePositioner)
{
	bool result = theDesignVariablePositionersPtr->addComponent(theDesignVariablePositioner);
	return result;
}


bool
ReliabilityDomain::addConstraintFunction(ConstraintFunction *theConstraintFunction)
{
	bool result = theConstraintFunctionsPtr->addComponent(theConstraintFunction);
	return result;
}


bool
ReliabilityDomain::addObjectiveFunction(ObjectiveFunction *theObjectiveFunction)
{
	bool result = theObjectiveFunctionsPtr->addComponent(theObjectiveFunction);
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

RandomVariable *
ReliabilityDomain::getRandomVariablePtrFromIndex(int index)
{
  if (index >= 0 && index < numRandomVariables)
    return this->getRandomVariablePtr(rvIndex[index]);

  else {
    opserr << "ReliabilityDomain::getRandomVariablePtrFromIndex -- index " << index << " out of bounds 0 ... " << numRandomVariables-1 << endln;
    return 0;
  }

}

int
ReliabilityDomain::getRandomVariableIndex(int tag)
{
  int index;

  // Find index of RV with specified tag
  for (index = 0; index < numRandomVariables; index++) {
    if (rvIndex[index] == tag)
      break;
  }

  if (index == numRandomVariables) {
    opserr << "ReliabilityDomain::getRandomVariableIndex -- rv with tag " << tag << " not found" << endln;
    return -1;
  }

  return index;
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

LimitStateFunction *
ReliabilityDomain::getLimitStateFunctionPtrFromIndex(int index)
{
  if (index >= 0 && index < numLimitStateFunctions)
    return this->getLimitStateFunctionPtr(lsfIndex[index]);

  else {
    opserr << "ReliabilityDomain::getLimitStateFunctionPtrFromIndex -- index " << index << " out of bounds 0 ... " << numLimitStateFunctions-1 << endln;
    return 0;
  }

}

int
ReliabilityDomain::getLimitStateFunctionIndex(int tag)
{
  int index;

  // Find index of LSF with specified tag
  for (index = 0; index < numLimitStateFunctions; index++) {
    if (lsfIndex[index] == tag)
      break;
  }

  if (index == numLimitStateFunctions) {
    opserr << "ReliabilityDomain::getLimitStateFunctionIndex -- lsf with tag " << tag << " not found" << endln;
    return -1;
  }

  return index;
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

//- Quan
DesignVariable *
ReliabilityDomain::getDesignVariablePtr(int tag)
{
	TaggedObject *theComponent = theDesignVariablesPtr->getComponentPtr(tag);
	
	DesignVariable *result = (DesignVariable *) theComponent;
	return result;
	
}

DesignVariablePositioner *
ReliabilityDomain::getDesignVariablePositionerPtr(int tag)
{
	TaggedObject *theComponent = theDesignVariablePositionersPtr->getComponentPtr(tag);
	DesignVariablePositioner *result = (DesignVariablePositioner *) theComponent;
	return result;
}

ConstraintFunction *
ReliabilityDomain::getConstraintFunctionPtr(int tag)
{
	TaggedObject *theComponent = theConstraintFunctionsPtr->getComponentPtr(tag);
	ConstraintFunction *result = (ConstraintFunction *) theComponent;
	return result;
}

ObjectiveFunction *
ReliabilityDomain::getObjectiveFunctionPtr(int tag)
{
	TaggedObject *theComponent = theObjectiveFunctionsPtr->getComponentPtr(tag);
	ObjectiveFunction *result = (ObjectiveFunction *) theComponent;
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
  
  if (theRV != 0) {

    // Find where RV is located
    int index;
    for (index = 0; index < numRandomVariables; index++) {
      if (rvIndex[index] == tag)
	break;
    }
    
    // Shift indices down by one
    for (int i = index; i < numRandomVariables-1; i++)
      rvIndex[i] = rvIndex[i+1];
    
    // Now remove the component
    theRandomVariablesPtr->removeComponent(tag);
    numRandomVariables--;
  }

  return 0;
}


int
ReliabilityDomain::removeCorrelationCoefficient(int tag)
{
  theCorrelationCoefficientsPtr->removeComponent(tag);

  return 0;
}

int
ReliabilityDomain::removeLimitStateFunction(int tag)
{
  LimitStateFunction *theLSF = (LimitStateFunction*) theLimitStateFunctionsPtr->getComponentPtr(tag);
  
  if (theLSF != 0) {

    // Find where LSF is located
    int index;
    for (index = 0; index < numLimitStateFunctions; index++) {
      if (lsfIndex[index] == tag)
	break;
    }
    
    // Shift indices down by one
    for (int i = index; i < numLimitStateFunctions-1; i++)
      lsfIndex[i] = lsfIndex[i+1];
    
    // Now remove the component
    theLimitStateFunctionsPtr->removeComponent(tag);
    numLimitStateFunctions--;
  }

  return 0;
}


//-Quan

int
ReliabilityDomain::removeDesignVariable(int tag)
{
	theDesignVariablesPtr->removeComponent(tag);

	return 0;
}

int
ReliabilityDomain::removeDesignVariablePositioner(int tag)
{
	theDesignVariablePositionersPtr->removeComponent(tag);

	return 0;
}

int
ReliabilityDomain::removeConstraintFunction(int tag)
{
	theConstraintFunctionsPtr->removeComponent(tag);

	return 0;
}

int
ReliabilityDomain::removeObjectiveFunction(int tag)
{
	theObjectiveFunctionsPtr->removeComponent(tag);

	return 0;
}







int
ReliabilityDomain::getNumberOfRandomVariables()
{
	return theRandomVariablesPtr->getNumComponents();
}


// Quan --

int ReliabilityDomain::setNumberOfRandomVariables( int pNum){
  opserr << "ReliabilityDomain::setNumberOfRandomVariables() - SHOULD NOT BE CALLED\n";
  return -1;
};

// --Quan
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

// -Quan

int
ReliabilityDomain::getNumberOfDesignVariables()
{
	return theDesignVariablesPtr->getNumComponents();
}

int
ReliabilityDomain::getNumberOfDesignVariablePositioners()
{
	return theDesignVariablePositionersPtr->getNumComponents();
}

int
ReliabilityDomain::getNumberOfConstraintFunctions()
{
	return theConstraintFunctionsPtr->getNumComponents();
}

int
ReliabilityDomain::getNumberOfObjectiveFunctions()
{
	return theObjectiveFunctionsPtr->getNumComponents();
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

void
ReliabilityDomain::Print(OPS_Stream &s, int flag) 
{
  s << "Current Reliability Domain Information\n";

  s << theRandomVariablesPtr->getNumComponents() << " random variables\n";
  if (flag == 1)
    theRandomVariablesPtr->Print(s, flag);

  s << theCorrelationCoefficientsPtr->getNumComponents() << " correlation coefficients\n";
  if (flag == 1)
    theCorrelationCoefficientsPtr->Print(s, flag);

  s << theLimitStateFunctionsPtr->getNumComponents() << " limit state functions\n";
  if (flag == 1)
    theLimitStateFunctionsPtr->Print(s, flag);
}
