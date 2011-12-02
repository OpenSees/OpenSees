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
**   Optimization module developed by:                                **
**   Quan Gu  (qgu@ucsd.edu)                                          **
**   Joel Conte (jpconte@ucsd.edu)                                    **
**   Philip Gill (pgill@ucsd.edu)                                     **
** ****************************************************************** */


//
// Written by Quan Gu (qgu@ucsd.edu)    March 2010
//

#ifndef OptimizationDomain_h
#define OptimizationDomain_h

#include <Parameter.h>


class Vector;

class TaggedObjectStorage;
class ConstraintFunction;
class ObjectiveFunction;
class DesignVariable;
class DesignVariablePositioner;


class OptimizationDomain
{

public:
	OptimizationDomain();
	virtual ~OptimizationDomain();

	// Member functions to add components to the domain

	virtual bool addDesignVariable(DesignVariable *theDesignVariable);
	virtual bool addDesignVariablePositioner(DesignVariablePositioner *theDesignVariablePositioner);
	virtual bool addObjectiveFunction(ObjectiveFunction *theObjectiveFunction);
	virtual bool addConstraintFunction(ConstraintFunction *theConstraintFunction);

	// Member functions to get components from the domain

	DesignVariable *getDesignVariablePtr(int tag);
	DesignVariablePositioner *getDesignVariablePositionerPtr(int tag);
	ObjectiveFunction *getObjectiveFunctionPtr(int tag);
	ConstraintFunction *getConstraintFunctionPtr(int tag);


	// Member functions giving information about the domain

	int getNumberOfDesignVariables(void);
	int getNumberOfDesignVariablePositioners(void);
	int getNumberOfObjectiveFunctions(void);
	int getNumberOfConstraintFunctions(void);

	// Member functions to remove single components from the domain

	int removeDesignVariable(int tag);
	int removeDesignVariablePositioner(int tag);
	int removeObjectiveFunction(int tag);
	int removeConstraintFunction(int tag);

	
	virtual void Print(OPS_Stream &s, int flag =0);

protected:

private:

	TaggedObjectStorage *theDesignVariablesPtr;
	TaggedObjectStorage *theObjectiveFunctionsPtr;
	TaggedObjectStorage *theConstraintFunctionsPtr;
	TaggedObjectStorage *theDesignVariablePositionersPtr;


};

#endif

