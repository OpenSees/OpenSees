// DesignVariablePositioner.cpp: implementation of the DesignVariablePositioner class.
//
//////////////////////////////////////////////////////////////////////

#include "DesignVariablePositioner.h"
#include <classTags.h>
#define DESIGN_VARIABLE_POSITIONER 100078




//////////////////////////////////////////////////////////////////////
DesignVariablePositioner::DesignVariablePositioner(int passedTag, ReliabilityDomain * theReliablityDomain, 
						   int passedDVnumber, DomainComponent *passedObject, const char **argv, int argc)
:ReliabilityDomainComponent(passedTag, DESIGN_VARIABLE_POSITIONER)
{
  dVNumber = passedDVnumber;
  theObject = passedObject;
  //  theInfo.theInt = passedDVnumber; // Used by random process discretizer
  
  DesignVariable * theDesignVariable = theReliablityDomain->getDesignVariablePtr(dVNumber);
  theDesignVariable->addDVPositioner(passedTag);
  
  int res = 0;
  if (theObject) 
    res = theObject->setParameter (argv, argc, theParameter);
  
  if (res < 0)
    opserr << "DesignVariablePositioner::DesignVariablePositioner "<< this->getTag() <<" -- unable to set parameter" << endln;
}


DesignVariablePositioner::~DesignVariablePositioner()
{

}


int
DesignVariablePositioner::update (double newValue)
{
  int result = theParameter.update(newValue);
  
  return result;
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
