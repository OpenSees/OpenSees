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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ReliabilityDomain.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef ReliabilityDomain_h
#define ReliabilityDomain_h

#include <RandomVariable.h>
#include <CorrelationCoefficient.h>
#include <LimitStateFunction.h>
#include <RandomVariablePositioner.h>
#include <ParameterPositioner.h>
#include <ModulatingFunction.h>
#include <Filter.h>
#include <Spectrum.h>

class TaggedObjectStorage;

class ReliabilityDomain
{

public:
	ReliabilityDomain();
	virtual ~ReliabilityDomain();

	// Member functions to add components to the domain
	virtual bool addRandomVariable(RandomVariable *theRandomVariable);
	virtual bool addCorrelationCoefficient(CorrelationCoefficient *theCorrelationCoefficient);
	virtual bool addLimitStateFunction(LimitStateFunction *theLimitStateFunction);
	virtual bool addRandomVariablePositioner(RandomVariablePositioner *theRandomVariablePositioner);
	virtual bool addParameterPositioner(ParameterPositioner *theParameterPositioner);
	virtual bool addModulatingFunction(ModulatingFunction *theModulatingFunction);
	virtual bool addFilter(Filter *theFilter);
	virtual bool addSpectrum(Spectrum *theSpectrum);

	// Member functions to get components from the domain
	RandomVariable *getRandomVariablePtr(int tag);
	CorrelationCoefficient *getCorrelationCoefficientPtr(int tag);
	LimitStateFunction *getLimitStateFunctionPtr(int tag);
	RandomVariablePositioner *getRandomVariablePositionerPtr(int tag);
	ParameterPositioner *getParameterPositionerPtr(int tag);
	ModulatingFunction *getModulatingFunction(int tag);
	Filter *getFilter(int tag);
	Spectrum *getSpectrum(int tag);

	// Member functions giving information about the domain
	int getNumberOfRandomVariables(void);
	int getNumberOfCorrelationCoefficients(void);
	int getNumberOfLimitStateFunctions(void);
	int getNumberOfRandomVariablePositioners(void);
	int getNumberOfParameterPositioners(void);
	int getNumberOfModulatingFunctions(void);
	int getNumberOfFilters(void);
	int getNumberOfSpectra(void);

	// Member functions to set/get active limit-state function
	int getTagOfActiveLimitStateFunction(void);
	void setTagOfActiveLimitStateFunction(int tag);

	// Member functions to remove single components from the domain
	int removeRandomVariablePositioner(int tag);
	int removePerformanceFunction(int tag);

protected:

private:
	TaggedObjectStorage *theRandomVariablesPtr;
	TaggedObjectStorage *theCorrelationCoefficientsPtr;
	TaggedObjectStorage *theLimitStateFunctionsPtr;
	TaggedObjectStorage *theRandomVariablePositionersPtr;
	TaggedObjectStorage *theParameterPositionersPtr;
	TaggedObjectStorage *theModulatingFunctionsPtr;
	TaggedObjectStorage *theFiltersPtr;
	TaggedObjectStorage *theSpectraPtr;
	int tagOfActiveLimitStateFunction;

};

#endif

