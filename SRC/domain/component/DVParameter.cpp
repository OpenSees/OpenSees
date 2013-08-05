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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/DVParameter.cpp,v $

#include <classTags.h>
#include <DVParameter.h>
#include <DesignVariable.h>

DVParameter::DVParameter(int passedTag, DesignVariable *theDV, Parameter *theParam)
  :Parameter(passedTag,1976), myDV(theDV), myParam(theParam), currentValue(0.0)
{
  if (myDV != 0)
    currentValue = myDV->getValue();
}


DVParameter::~DVParameter()
{
  if (myParam != 0)
    delete myParam;
}

int
DVParameter::update(int newValue)
{
  currentValue = newValue;

  return 0;
}

int
DVParameter::update(double newValue)
{
  currentValue = newValue;

  myDV->setValue(newValue);

  if (myParam != 0) {
    myParam->setValue(newValue);
    myParam->update(newValue);
  }

  return 0;
}

int
DVParameter::activate(bool active)
{
  if (myParam != 0)
    return myParam->activate(active);
  else
    return 0;
}

double
DVParameter::getValue(void)
{
  //return currentValue;
  return myDV->getValue();
}

void
DVParameter::setValue(double newValue)
{
  currentValue = newValue;

  myDV->setValue(newValue);

  if (myParam != 0) {
    myParam->setValue(newValue);
    myParam->update(newValue);
  }
}

bool
DVParameter::isImplicit(void)
{
  return false;
}

double
DVParameter::getSensitivity(int index)
{
  return Parameter::getSensitivity(index);
}

double
DVParameter::getPerturbation(void)
{
  //return 0.001*myDV->getStdv();

  // Can probably do something better here
  return 0.001*(myDV->getUpperBound() - myDV->getLowerBound()); 
}

int 
DVParameter::getPointerTag(void) 
{
    return myDV->getTag();
}

void
DVParameter::Print(OPS_Stream &s, int flag)  
{
  s << "DVParameter, tag = " << this->getTag() << endln;
  myDV->Print(s, flag);
  if (myParam != 0) {
    myParam->Print(s, flag);
  }
}

void
DVParameter::setDomain(Domain *theDomain)
{
  if (myParam != 0)
    myParam->setDomain(theDomain);
  return;
}

int 
DVParameter::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int 
DVParameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}
