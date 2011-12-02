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
                                                                        
// $Revision: 1.5 $
// $Date: 2007-02-17 21:27:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/domain/components/RandomVariable.cpp,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#include <RandomVariable.h>
#include <classTags.h>
#include <OPS_Globals.h>

RandomVariable::RandomVariable(int tag, int classTag, double startVal)
  :ReliabilityDomainComponent(tag,classTag), rvNumber(tag), startValue(startVal)
{

}

RandomVariable::~RandomVariable()
{
}



int 
RandomVariable::setNewTag(int newTag)
{
  TaggedObject::setTag(newTag);
  
  return 0;
}

double
RandomVariable::getParameter1()
{
  opserr << "No parameter 1 in r.v. #" << this->getTag() << endln;
  return 0.0;
}

double
RandomVariable::getParameter2()
{
  opserr << "No parameter 2 in r.v. #" << this->getTag() << endln;
  return 0.0;
}

double
RandomVariable::getParameter3()
{
  opserr << "No parameter 3 in r.v. #" << this->getTag() << endln;
  return 0.0;
}

double
RandomVariable::getParameter4()
{
  opserr << "No parameter 4 in r.v. #" << this->getTag() << endln;
  return 0.0;  
}
