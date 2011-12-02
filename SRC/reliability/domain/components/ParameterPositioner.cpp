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
                                                                        
// $Revision: 1.1 $
// $Date: 2003-03-04 00:44:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ParameterPositioner.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ParameterPositioner.h>
#include <classTags.h>


ParameterPositioner::ParameterPositioner (int passedTag,
		DomainComponent *object,
		const char **argv, int argc)
:ReliabilityDomainComponent(passedTag, RANDOM_VARIABLE_POSITIONER)
{
	tag = passedTag;
	theObject = object;

	if (theObject) 
		parameterID = theObject->setParameter (argv, argc, theInfo);

	if (parameterID < 0)
		opserr << "ParameterPositioner::ParameterPositioner "<< tag <<" -- unable to set parameter" << endln;
}


ParameterPositioner::~ParameterPositioner()
{
}


int
ParameterPositioner::update (double newValue)
{
	theInfo.theDouble = newValue;

	if (parameterID >= 0)
		return theObject->updateParameter (parameterID, theInfo);
	else
		return -1;
}

int
ParameterPositioner::activate(bool active)
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
ParameterPositioner::Print(OPS_Stream &s, int flag)  
{
}


