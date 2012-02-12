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
** ****************************************************************** */
                                                                        
// $Revision: 1.9 $
// $Date: 2008-08-26 15:43:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/RVParameter.cpp,v $

#include <classTags.h>
#include <RVParameter.h>
#include <RandomVariable.h>

RVParameter::RVParameter(int passedTag, RandomVariable *theRV, Parameter *theParam)
  :Parameter(passedTag,1976), myRV(theRV), myParam(theParam), currentValue(0.0)
{
  if (myRV != 0)
    currentValue = myRV->getCurrentValue();
}


RVParameter::~RVParameter()
{
  if (myParam != 0)
    delete myParam;
}

int
RVParameter::update(int newValue)
{
  currentValue = newValue;

  return 0;
}

int
RVParameter::update(double newValue)
{
  currentValue = newValue;

  myRV->setCurrentValue(newValue);

  if (myParam != 0)
    myParam->setValue(newValue);

  return 0;
}

int
RVParameter::activate(bool active)
{
  return 0;
}

double
RVParameter::getValue(void)
{
  //return currentValue;
  return myRV->getCurrentValue();
}

void
RVParameter::setValue(double newValue)
{
  currentValue = newValue;

  myRV->setCurrentValue(newValue);

  if (myParam != 0)
    myParam->setValue(newValue);
}

double
RVParameter::getPerturbation(void)
{
  return 0.001*myRV->getStdv();
}

void
RVParameter::Print(OPS_Stream &s, int flag)  
{
  s << "RVParameter, tag = " << this->getTag() << endln;
  myRV->Print(s, flag);
  if (myParam != 0) {
    myParam->Print(s, flag);
  }
}

void
RVParameter::setDomain(Domain *theDomain)
{
  return;
}

int 
RVParameter::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int 
RVParameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}
