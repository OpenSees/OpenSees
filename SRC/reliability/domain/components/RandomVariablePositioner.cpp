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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-10-31 21:36:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariablePositioner.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <RandomVariablePositioner.h>
#include <classTags.h>

RandomVariablePositioner::RandomVariablePositioner(int passedTag,
		int passedRVindex,
		DomainComponent *object,
		const char **argv, int argc)
  :ReliabilityDomainComponent(passedTag, RANDOM_VARIABLE_POSITIONER),
   rvIndex(passedRVindex), theParameter(passedTag, object, argv, argc)
{

}

RandomVariablePositioner::RandomVariablePositioner(int passedTag,
		int passedRVindex, Parameter &param)
  :ReliabilityDomainComponent(passedTag, RANDOM_VARIABLE_POSITIONER),
   rvIndex(passedRVindex), theParameter(param)
{

}

RandomVariablePositioner::~RandomVariablePositioner()
{
  // Nothing to do
}

int
RandomVariablePositioner::update(double newValue)
{
  return theParameter.update(newValue);
}

int
RandomVariablePositioner::activate(bool active)
{
  return theParameter.activate(active);
}

void
RandomVariablePositioner::Print(OPS_Stream &s, int flag)  
{
  s << "RandomVariablePositoner, tag = " << this->getTag() << '\n';
  s << "\trvIndex = " << rvIndex << '\n';
  theParameter.Print(s, flag);
}

int 
RandomVariablePositioner::getRvIndex(void)
{
  return rvIndex;
}

int 
RandomVariablePositioner::setNewTag(int newTag)
{
  this->setTag(newTag);

  return 0;
}

int 
RandomVariablePositioner::setRvIndex(int newRvIndex)
{
  rvIndex = newRvIndex;
  
  return 0;
}
