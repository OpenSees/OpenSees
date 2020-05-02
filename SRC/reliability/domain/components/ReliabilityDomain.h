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
                                                                        
// $Revision: 1.19 $
// $Date: 2008-10-22 16:42:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ReliabilityDomain.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef ReliabilityDomain_h
#define ReliabilityDomain_h

#include <RandomVariable.h>
#include <CorrelationCoefficient.h>
#include <Cutset.h>
#include <LimitStateFunction.h>
#include <ModulatingFunction.h>
#include <Filter.h>
#include <Spectrum.h>

//Quan
#include <DesignVariable.h>
#include <DesignVariablePositioner.h>
#include <ConstraintFunction.h>
#include <ObjectiveFunction.h>

class Vector;

class TaggedObjectStorage;
class ConstraintFunction;
class ObjectiveFunction;
class DesignVariable;
class DesignVariablePositioner;

class RandomVariableIter;
class LimitStateFunctionIter;
class CutsetIter;
class CorrelationCoefficientIter;
class ModulatingFunctionIter;
class FilterIter;
class SpectrumIter;

class ReliabilityDomain
{

public:
	ReliabilityDomain(Domain *theOpenSeesDomain);
	virtual ~ReliabilityDomain();

	// Member functions to add components to the domain
	virtual bool addRandomVariable(RandomVariable *theRandomVariable);
	virtual bool addCorrelationCoefficient(CorrelationCoefficient *theCorrelationCoefficient);
	virtual bool addLimitStateFunction(LimitStateFunction *theLimitStateFunction);
	virtual bool addCutset(Cutset *theCutset);
	virtual bool addModulatingFunction(ModulatingFunction *theModulatingFunction);
	virtual bool addFilter(Filter *theFilter);
	virtual bool addSpectrum(Spectrum *theSpectrum);

	virtual bool addDesignVariable(DesignVariable *theDesignVariable);
	virtual bool addDesignVariablePositioner(DesignVariablePositioner *theDesignVariablePositioner);
	virtual bool addObjectiveFunction(ObjectiveFunction *theObjectiveFunction);
	virtual bool addConstraintFunction(ConstraintFunction *theConstraintFunction);

	// Member functions to get components from the domain
	RandomVariable *getRandomVariablePtr(int tag);
	// Following two methods to map index to tag and vice versa
	RandomVariable *getRandomVariablePtrFromIndex(int index);
	int getRandomVariableIndex(int tag);
    
    // for parameter to RV map
    int getRandomVariableIndexFromParameterIndex(int param_index);
    int getParameterIndexFromRandomVariableIndex(int rv_index);

	LimitStateFunction *getLimitStateFunctionPtr(int tag);
	// Following two methods to map index to tag and vice versa
	LimitStateFunction *getLimitStateFunctionPtrFromIndex(int index);
	int getLimitStateFunctionIndex(int tag);

	Cutset *getCutsetPtr(int tag);
	// Following two methods to map index to tag and vice versa
	Cutset *getCutsetPtrFromIndex(int index);
	int getCutsetIndex(int tag);

	CorrelationCoefficient *getCorrelationCoefficientPtr(int tag);
	ModulatingFunction *getModulatingFunction(int tag);
	Filter *getFilter(int tag);
	Spectrum *getSpectrum(int tag);

	DesignVariable *getDesignVariablePtr(int tag);
	DesignVariablePositioner *getDesignVariablePositionerPtr(int tag);
	ObjectiveFunction *getObjectiveFunctionPtr(int tag);
	ConstraintFunction *getConstraintFunctionPtr(int tag);
	int setNumberOfRandomVariables( int num);

	// Member functions giving information about the domain
	int getNumberOfRandomVariables(void);
	int getNumberOfCorrelationCoefficients(void);
	int getNumberOfLimitStateFunctions(void);
	int getNumberOfCutsets(void);
	int getNumberOfModulatingFunctions(void);
	int getNumberOfFilters(void);
	int getNumberOfSpectra(void);

	int getNumberOfDesignVariables(void);
	int getNumberOfDesignVariablePositioners(void);
	int getNumberOfObjectiveFunctions(void);
	int getNumberOfConstraintFunctions(void);

	// Member functions to set/get active limit-state function
	int getTagOfActiveLimitStateFunction(void); // Probably not need anymore -- MHS
	void setTagOfActiveLimitStateFunction(int tag); // Probably not need anymore -- MHS

	// Member functions to remove single components from the domain
	void clearAll(void);
	int removeRandomVariable(int tag);
	int removeCorrelationCoefficient(int tag);
	int removeLimitStateFunction(int tag);
	int removeCutset(int tag);
	int removeModulatingFunction(int tag);
	int removeFilter(int tag);
	int removeSpectrum(int tag);

	int removeDesignVariable(int tag);
	int removeDesignVariablePositioner(int tag);
	int removeObjectiveFunction(int tag);
	int removeConstraintFunction(int tag);

	int removeAllParameterPositioners(void);

	RandomVariableIter &getRandomVariables(void);
	LimitStateFunctionIter &getLimitStateFunctions(void);
	CorrelationCoefficientIter &getCorrelationCoefficients(void);
	CutsetIter &getCutsets(void);
	ModulatingFunctionIter &getModulatingFunctions(void);
	FilterIter &getFilters(void);
	SpectrumIter &getSpectra(void);
	
	virtual void Print(OPS_Stream &s, int flag =0);

protected:

private:
	TaggedObjectStorage *theRandomVariablesPtr;
	TaggedObjectStorage *theCorrelationCoefficientsPtr;
	TaggedObjectStorage *theLimitStateFunctionsPtr;
	TaggedObjectStorage *theCutsetsPtr;
	TaggedObjectStorage *theModulatingFunctionsPtr;
	TaggedObjectStorage *theFiltersPtr;
	TaggedObjectStorage *theSpectraPtr;
	int tagOfActiveLimitStateFunction; // Probably not need anymore -- MHS

	TaggedObjectStorage *theDesignVariablesPtr;
	TaggedObjectStorage *theObjectiveFunctionsPtr;
	TaggedObjectStorage *theConstraintFunctionsPtr;
	TaggedObjectStorage *theDesignVariablePositionersPtr;

	RandomVariableIter *theRVIter;
	LimitStateFunctionIter *theLSFIter;
	CutsetIter *theCutIter;
	CorrelationCoefficientIter *theCCIter;
	ModulatingFunctionIter *theMFIter;
	FilterIter *theFilterIter;
	SpectrumIter *theSpectrumIter;
    
    // store RV to parameter map in the reliability domain
    Domain *theOpenSeesDomain;

	// Integer array: index[i] = tag of component i
	// Should put these in another class eventually because
	// we may need to do the same thing for lsf, cc, etc. -- MHS
	int *rvIndex;
	enum {rvSize_init = 100};
	enum {rvSize_grow = 20};
	int rvSize;
	int numRandomVariables;

	int *lsfIndex;
	enum {lsfSize_init = 10};
	enum {lsfSize_grow = 2};
	int lsfSize;
	int numLimitStateFunctions;
	
	int *cutsetIndex;
	enum {cutsetSize_init = 10};
	enum {cutsetSize_grow = 2};
	int cutsetSize;
	int numCutsets;
};

#endif

