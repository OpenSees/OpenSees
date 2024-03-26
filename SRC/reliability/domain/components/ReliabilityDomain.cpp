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
                                                                        
// $Revision: 1.21 $
// $Date: 2008-10-22 16:42:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ReliabilityDomain.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ReliabilityDomain.h>
#include <Domain.h>
#include <Vector.h>

#include <CorrelationCoefficient.h>
#include <Cutset.h>
#include <RandomVariable.h>
#include <LimitStateFunction.h>
#include <Parameter.h>
#include <ArrayOfTaggedObjects.h>
#include <MapOfTaggedObjects.h>
#include <ModulatingFunction.h>
#include <Filter.h>
#include <Spectrum.h>

#include <RandomVariableIter.h>
#include <LimitStateFunctionIter.h>
#include <CorrelationCoefficientIter.h>
#include <CutsetIter.h>
#include <ModulatingFunctionIter.h>
#include <FilterIter.h>
#include <SpectrumIter.h>

ReliabilityDomain::ReliabilityDomain(Domain *passedDomain):
    theOpenSeesDomain(passedDomain),
    numRandomVariables(0), numLimitStateFunctions(0), numCutsets(0)
{

	theRandomVariablesPtr = new ArrayOfTaggedObjects (256);
	theCorrelationCoefficientsPtr = new ArrayOfTaggedObjects (256);
	theLimitStateFunctionsPtr = new ArrayOfTaggedObjects (256);
	theCutsetsPtr = new ArrayOfTaggedObjects (256);
	theModulatingFunctionsPtr = new ArrayOfTaggedObjects (256);
	theFiltersPtr = new ArrayOfTaggedObjects (256);
	theSpectraPtr = new ArrayOfTaggedObjects (256);

	theDesignVariablesPtr = new ArrayOfTaggedObjects (256);
	theDesignVariablePositionersPtr = new ArrayOfTaggedObjects (256);
	theObjectiveFunctionsPtr = new ArrayOfTaggedObjects (256);
	theConstraintFunctionsPtr = new ArrayOfTaggedObjects (256);

	tagOfActiveLimitStateFunction = 1;

	theRVIter = new RandomVariableIter(theRandomVariablesPtr);
	theLSFIter = new LimitStateFunctionIter(theLimitStateFunctionsPtr);
	theCutIter = new CutsetIter(theCutsetsPtr);
	theCCIter = new CorrelationCoefficientIter(theCorrelationCoefficientsPtr);
	theMFIter = new ModulatingFunctionIter(theModulatingFunctionsPtr);
	theFilterIter = new FilterIter(theFiltersPtr);
	theSpectrumIter = new SpectrumIter(theSpectraPtr);

	rvIndex = new int[rvSize_init];
	rvSize = rvSize_init;

	lsfIndex = new int[lsfSize_init];
	lsfSize = lsfSize_init;
	
	cutsetIndex = new int[cutsetSize_init];
	cutsetSize = cutsetSize_init;
}

void
ReliabilityDomain::clearAll()
{
  if (theRandomVariablesPtr != 0) {
    theRandomVariablesPtr->clearAll();
    numRandomVariables = 0;
  }
  if (theLimitStateFunctionsPtr != 0) {
    theLimitStateFunctionsPtr->clearAll();
    numLimitStateFunctions = 0;
  }
  if (theCutsetsPtr != 0) {
    theCutsetsPtr->clearAll();
    numCutsets = 0;
  }
 
  if (theCorrelationCoefficientsPtr != 0) 
    theCorrelationCoefficientsPtr->clearAll();
  if (theModulatingFunctionsPtr != 0)
    theModulatingFunctionsPtr->clearAll();
  if (theSpectraPtr != 0)
    theSpectraPtr->clearAll();
  if (theFiltersPtr != 0)
    theFiltersPtr->clearAll();

  if (theDesignVariablesPtr != 0) 
    theDesignVariablesPtr->clearAll();
  if (theDesignVariablePositionersPtr != 0) 
    theDesignVariablePositionersPtr->clearAll();
  if (theConstraintFunctionsPtr != 0) 
    theConstraintFunctionsPtr->clearAll();
  if (theObjectiveFunctionsPtr != 0) 
    theObjectiveFunctionsPtr->clearAll();
}

ReliabilityDomain::~ReliabilityDomain()
{
  if (theRandomVariablesPtr != 0) {
    theRandomVariablesPtr->clearAll();
    delete theRandomVariablesPtr;
  }
  if (theLimitStateFunctionsPtr != 0) {
    theLimitStateFunctionsPtr->clearAll();
    delete theLimitStateFunctionsPtr;
  }
  if (theCutsetsPtr != 0) {
    theCutsetsPtr->clearAll();
    delete theCutsetsPtr;
  }
  if (theCorrelationCoefficientsPtr != 0) {
    theCorrelationCoefficientsPtr->clearAll();
    delete theCorrelationCoefficientsPtr;
  }
  if (theModulatingFunctionsPtr != 0) {
    theModulatingFunctionsPtr->clearAll();
    delete theModulatingFunctionsPtr;
  }
  if (theSpectraPtr != 0) {
    theSpectraPtr->clearAll();
    delete theSpectraPtr;
  }
  if (theFiltersPtr != 0) {
    theFiltersPtr->clearAll();
    delete theFiltersPtr;
  }

  if (theDesignVariablesPtr != 0) {
    theDesignVariablesPtr->clearAll();
    delete theDesignVariablesPtr;
  }
  if (theDesignVariablePositionersPtr != 0) {
    theDesignVariablePositionersPtr->clearAll();
    delete theDesignVariablePositionersPtr;
  }
  if (theConstraintFunctionsPtr != 0) {
    theConstraintFunctionsPtr->clearAll();
    delete theConstraintFunctionsPtr;
  }
  if (theObjectiveFunctionsPtr != 0) {
    theObjectiveFunctionsPtr->clearAll();
    delete theObjectiveFunctionsPtr;
  }

  if (theRVIter != 0)
    delete theRVIter;
  if (theLSFIter != 0)
    delete theLSFIter;
  if (theCutIter != 0)
    delete theCutIter;
  if (theCCIter != 0)
    delete theCCIter;
  if (theMFIter != 0)
    delete theMFIter;
  if (theFilterIter != 0)
    delete theFilterIter;
  if (theSpectrumIter != 0)
    delete theSpectrumIter;

  if (rvIndex != 0)
    delete [] rvIndex;
  if (lsfIndex != 0)
    delete [] lsfIndex;
  if (cutsetIndex != 0)
    delete [] cutsetIndex;
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
ReliabilityDomain::addCutset(Cutset *theCutset)
{
  bool result = theCutsetsPtr->addComponent(theCutset);

  if (result == true) {

    // Array is full
    if (numCutsets == cutsetSize) {

      // Increase size and allocate new array
      cutsetSize += cutsetSize_grow;
      int *tmp_cutsetIndex = new int[cutsetSize];

      // Copy values from old array to new
      for (int i = 0; i < numCutsets; i++)
	    tmp_cutsetIndex[i] = cutsetIndex[i];

      // Get rid of old array
      delete [] cutsetIndex;

      // Set pointer to new array
      cutsetIndex = tmp_cutsetIndex;
    }

    // Add to index
    cutsetIndex[numCutsets] = theCutset->getTag();
    numCutsets++;
  }

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

CutsetIter &
ReliabilityDomain::getCutsets(void)
{
  theCutIter->reset();
  return *theCutIter;
}

ModulatingFunctionIter &
ReliabilityDomain::getModulatingFunctions(void)
{
  theMFIter->reset();
  return *theMFIter;
}

FilterIter &
ReliabilityDomain::getFilters(void)
{
  theFilterIter->reset();
  return *theFilterIter;
}

SpectrumIter &
ReliabilityDomain::getSpectra(void)
{
  theSpectrumIter->reset();
  return *theSpectrumIter;
}

int 
ReliabilityDomain::getRandomVariableIndexFromParameterIndex(int param_index)
{
    // note this map is not guaranteed (there are parameters that are not RVs)
    int numberOfParameters = theOpenSeesDomain->getNumParameters();
    int result = -1;
    
    if (param_index >= 0 && param_index < numberOfParameters) {
        // map from parameters to random variables 
        Parameter *theParam = theOpenSeesDomain->getParameterFromIndex(param_index);
        if ( strcmp( theParam->getType(), "RandomVariable" ) == 0 )
            result = getRandomVariableIndex( theParam->getPointerTag() );
    }
    else {
        opserr << "ReliabilityDomain::getRandomVariableIndexFromParameterIndex -- index " << param_index 
            << " out of bounds 0 ... " << numberOfParameters-1 << endln;
        return -1;
    }
    
    return result;
}

int 
ReliabilityDomain::getParameterIndexFromRandomVariableIndex(int rv_index)
{
    // note this map IS guaranteed because there is a parameter for every RV
    int numberOfParameters = theOpenSeesDomain->getNumParameters();
    int result;
    int *RVmap = new int[numRandomVariables];
    
    for (int j = 0; j < numberOfParameters; j++) {
        Parameter *theParam = theOpenSeesDomain->getParameterFromIndex(j);
        if ( strcmp( theParam->getType(), "RandomVariable" ) == 0 ) {
            result = theParam->getPointerTag();
            RVmap[ getRandomVariableIndex(result) ] = j;
        }
    }
    
    if (rv_index >= 0 && rv_index < numRandomVariables) {
        // map from random variables to parameters
        result = RVmap[rv_index];
    }
    else {
        opserr << "ReliabilityDomain::getParameterIndexFromRandomVariableIndex -- index " << rv_index 
            << " out of bounds 0 ... " << numRandomVariables-1 << endln;
        return -1;
    }
    delete [] RVmap;
    return result;
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


Cutset * 
ReliabilityDomain::getCutsetPtr(int tag)
{
	TaggedObject *theComponent = theCutsetsPtr->getComponentPtr(tag);
	if ( theComponent == 0 )
		return 0;
	Cutset *result = (Cutset *) theComponent;
	return result;
}

Cutset *
ReliabilityDomain::getCutsetPtrFromIndex(int index)
{
  if (index >= 0 && index < numCutsets)
    return this->getCutsetPtr(cutsetIndex[index]);

  else {
    opserr << "ReliabilityDomain::getCutsetPtrFromIndex -- index " << index << " out of bounds 0 ... " << numCutsets-1 << endln;
    return 0;
  }

}

int
ReliabilityDomain::getCutsetIndex(int tag)
{
  int index;

  // Find index of cutset with specified tag
  for (index = 0; index < numCutsets; index++) {
    if (cutsetIndex[index] == tag)
      break;
  }

  if (index == numCutsets) {
    opserr << "ReliabilityDomain::getCutsetIndex -- cutset with tag " << tag << " not found" << endln;
    return -1;
  }

  return index;
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
    for (int i = index; i < numRandomVariables-1; i++) {
      rvIndex[i] = rvIndex[i+1];
    }
    
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


int
ReliabilityDomain::removeCutset(int tag)
{
  Cutset *theCutset = (Cutset*) theCutsetsPtr->getComponentPtr(tag);
  
  if (theCutset != 0) {

    // Find where cutset is located
    int index;
    for (index = 0; index < numCutsets; index++) {
      if (cutsetIndex[index] == tag)
	  break;
    }
    
    // Shift indices down by one
    for (int i = index; i < numCutsets-1; i++)
      cutsetIndex[i] = cutsetIndex[i+1];
    
    // Now remove the component
    theCutsetsPtr->removeComponent(tag);
    numCutsets--;
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
ReliabilityDomain::getNumberOfCutsets()
{
	return theCutsetsPtr->getNumComponents();
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
