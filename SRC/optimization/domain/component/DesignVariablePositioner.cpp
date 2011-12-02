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

#include "DesignVariable.h"
#include "DesignVariablePositioner.h"
#include <classTags.h>


#define DESIGN_VARIABLE_POSITIONER 100078734



DesignVariablePositioner::DesignVariablePositioner(int passedTag, OptimizationDomain * theReliablityDomain, 
												   int passedDVnumber, DomainComponent *passedObject, const char **argv, int argc)
:OptimizationDomainComponent(passedTag, DESIGN_VARIABLE_POSITIONER),
theParameter(passedTag,  passedObject, argv, argc ), dVNumber(passedDVnumber)
{
	DesignVariable * theDesignVariable = theReliablityDomain->getDesignVariablePtr(dVNumber);
	theDesignVariable->addDVPositioner(passedTag);


	
	/* -----------NEED TO DO!!! Quan  2010 ----------

	dVNumber = passedDVnumber;
	theObject = passedObject;
	theInfo.theInt = passedDVnumber; // Used by random process discretizer

	DesignVariable * theDesignVariable = theOptimizationDomain->getDesignVariablePtr(dVNumber);
	theDesignVariable->addDVPositioner(passedTag);



	if (theObject) 
		parameterID = theObject->setParameter (argv, argc, theInfo);






	if (parameterID < 0)
		opserr << "DesignVariablePositioner::DesignVariablePositioner "<< this->getTag() <<" -- unable to set parameter" << endln;
*/

}


DesignVariablePositioner::DesignVariablePositioner(int passedTag, OptimizationDomain  * theOptimizationDomain, int passedDVnumber, Parameter &param)
:OptimizationDomainComponent(passedTag, DESIGN_VARIABLE_POSITIONER), dVNumber(passedDVnumber),
theParameter(param)
{
	DesignVariable * theDesignVariable = theOptimizationDomain->getDesignVariablePtr(dVNumber);
	theDesignVariable->addDVPositioner(passedTag);
}

DesignVariablePositioner::~DesignVariablePositioner()
{
}


int
DesignVariablePositioner::update (double newValue)
{
	return theParameter.update(newValue);
}




void
DesignVariablePositioner::Print(OPS_Stream &s, int flag)  
{
}


int 
DesignVariablePositioner::getDVNumber(void)
{
	return dVNumber;
}


int 
DesignVariablePositioner::setNewTag(int newTag)
{
	this->setTag(newTag);
	return 0;
}


int 
DesignVariablePositioner::setDVNumber(int newDvNumber)
{
	dVNumber = newDvNumber;
	return 0;
}
