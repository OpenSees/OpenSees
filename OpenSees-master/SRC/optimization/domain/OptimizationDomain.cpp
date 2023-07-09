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

#include <OptimizationDomain.h>
#include <ArrayOfTaggedObjects.h>

/*#include <DesignVariable.h>
#include <DesignVariablePositioner.h>
#include <ConstraintFunction.h>
#include <ObjectiveFunction.h>
*/
OptimizationDomain::OptimizationDomain()
{

	theDesignVariablesPtr = new ArrayOfTaggedObjects (256);
	theDesignVariablePositionersPtr = new ArrayOfTaggedObjects (256);
	theObjectiveFunctionsPtr = new ArrayOfTaggedObjects (256);
	theConstraintFunctionsPtr = new ArrayOfTaggedObjects (256);

}

OptimizationDomain::~OptimizationDomain()
{

	
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

}



bool
OptimizationDomain::addDesignVariable(DesignVariable *theDesignVariable)
{
	bool result = theDesignVariablesPtr->addComponent((TaggedObject *)theDesignVariable);
	return result;
}


bool
OptimizationDomain::addDesignVariablePositioner(DesignVariablePositioner *theDesignVariablePositioner)
{
	bool result = theDesignVariablePositionersPtr->addComponent((TaggedObject *)theDesignVariablePositioner);
	return result;
}


bool
OptimizationDomain::addConstraintFunction(ConstraintFunction *theConstraintFunction)
{
	bool result = theConstraintFunctionsPtr->addComponent((TaggedObject *)theConstraintFunction);
	return result;
}


bool
OptimizationDomain::addObjectiveFunction(ObjectiveFunction *theObjectiveFunction)
{
	bool result = theObjectiveFunctionsPtr->addComponent((TaggedObject *)theObjectiveFunction);
	return result;
}



DesignVariable *
OptimizationDomain::getDesignVariablePtr(int tag)
{
	TaggedObject *theComponent = theDesignVariablesPtr->getComponentPtr(tag);
	
	DesignVariable *result = (DesignVariable *) theComponent;
	return result;
	
}

DesignVariablePositioner *
OptimizationDomain::getDesignVariablePositionerPtr(int tag)
{
	TaggedObject *theComponent = theDesignVariablePositionersPtr->getComponentPtr(tag);
	DesignVariablePositioner *result = (DesignVariablePositioner *) theComponent;
	return result;
}

ConstraintFunction *
OptimizationDomain::getConstraintFunctionPtr(int tag)
{
	TaggedObject *theComponent = theConstraintFunctionsPtr->getComponentPtr(tag);
	ConstraintFunction *result = (ConstraintFunction *) theComponent;
	return result;
}

ObjectiveFunction *
OptimizationDomain::getObjectiveFunctionPtr(int tag)
{
	TaggedObject *theComponent = theObjectiveFunctionsPtr->getComponentPtr(tag);
	ObjectiveFunction *result = (ObjectiveFunction *) theComponent;
	return result;
}



int
OptimizationDomain::removeDesignVariable(int tag)
{
	theDesignVariablesPtr->removeComponent(tag);

	return 0;
}

int
OptimizationDomain::removeDesignVariablePositioner(int tag)
{
	theDesignVariablePositionersPtr->removeComponent(tag);

	return 0;
}

int
OptimizationDomain::removeConstraintFunction(int tag)
{
	theConstraintFunctionsPtr->removeComponent(tag);

	return 0;
}

int
OptimizationDomain::removeObjectiveFunction(int tag)
{
	theObjectiveFunctionsPtr->removeComponent(tag);

	return 0;
}



int
OptimizationDomain::getNumberOfDesignVariables()
{
	return theDesignVariablesPtr->getNumComponents();
}

int
OptimizationDomain::getNumberOfDesignVariablePositioners()
{
	return theDesignVariablePositionersPtr->getNumComponents();
}

int
OptimizationDomain::getNumberOfConstraintFunctions()
{
	return theConstraintFunctionsPtr->getNumComponents();
}

int
OptimizationDomain::getNumberOfObjectiveFunctions()
{
	return theObjectiveFunctionsPtr->getNumComponents();
}

void
OptimizationDomain::Print(OPS_Stream &s, int flag) 
{
  s << "Current Optimization Domain Information\n";

  s << theDesignVariablesPtr->getNumComponents() << " random variables\n";
  if (flag == 1)
    theDesignVariablesPtr->Print(s, flag);

   
}
