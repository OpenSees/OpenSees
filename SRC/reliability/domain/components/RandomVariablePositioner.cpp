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
                                                                        
// $Revision: 1.6 $
// $Date: 2003-10-27 23:04:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariablePositioner.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <RandomVariablePositioner.h>
#include <classTags.h>


RandomVariablePositioner::RandomVariablePositioner (int passedTag,
		int passedRVnumber,
		DomainComponent *object,
		const char **argv, int argc)
:ReliabilityDomainComponent(passedTag, RANDOM_VARIABLE_POSITIONER)
{
	rvNumber = passedRVnumber;
	theObject = object;
	theInfo.theInt = rvNumber; // Used by random process discretizer

	if (theObject) 
		parameterID = theObject->setParameter (argv, argc, theInfo);

	if (parameterID < 0)
		opserr << "RandomVariablePositioner::RandomVariablePositioner "<< this->getTag() <<" -- unable to set parameter" << endln;
}


RandomVariablePositioner::~RandomVariablePositioner()
{
}


int
RandomVariablePositioner::update (double newValue)
{
	theInfo.theDouble = newValue;

	if (parameterID >= 0)
		return theObject->updateParameter (parameterID, theInfo);
	else
		return -1;
}

int
RandomVariablePositioner::activate(bool active)
{
	if (active) {
		theObject->activateParameter(parameterID);
	}
	else {
		theObject->activateParameter(0);
	}
	return 0;
}


void
RandomVariablePositioner::Print(OPS_Stream &s, int flag)  
{
}


int 
RandomVariablePositioner::getRvNumber(void)
{
	return rvNumber;
}


int 
RandomVariablePositioner::setNewTag(int newTag)
{
	this->setTag(newTag);
	return 0;
}


int 
RandomVariablePositioner::setRvNumber(int newRvNumber)
{
	rvNumber = newRvNumber;
	return 0;
}
