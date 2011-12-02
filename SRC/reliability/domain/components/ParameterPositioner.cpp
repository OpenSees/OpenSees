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
// $Date: 2008-05-22 20:05:28 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/ParameterPositioner.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#include <ParameterPositioner.h>
#include <classTags.h>
#include <DomainComponent.h>
#include <Parameter.h>

ParameterPositioner::ParameterPositioner (int passedTag,
					  DomainComponent *parentObject,
					  const char **argv, int argc)
  :ReliabilityDomainComponent(passedTag, PARAMETER_POSITIONER),
   theParameter(passedTag, parentObject, argv, argc), gradNumber(passedTag)
{

}

ParameterPositioner::ParameterPositioner (int passedTag, Parameter &param)
  :ReliabilityDomainComponent(passedTag, PARAMETER_POSITIONER),
   theParameter(param), gradNumber(passedTag)
{

}

ParameterPositioner::~ParameterPositioner()
{

}

int
ParameterPositioner::update(int newValue)
{
  return theParameter.update(newValue);
}

int
ParameterPositioner::update(double newValue)
{
  return theParameter.update(newValue);
}

int
ParameterPositioner::activate(bool active)
{
  return theParameter.activate(active);
}

void
ParameterPositioner::Print(OPS_Stream &s, int flag)  
{
  s << "ParameterPositioner, tag = " << this->getTag() << endln;
  s << "\tgrad number " << gradNumber << endln;
}
