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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/LoadFactorParameter.cpp,v $

#include <classTags.h>
#include <LoadFactorParameter.h>
#include <LoadPattern.h>

LoadFactorParameter::LoadFactorParameter(int passedTag, LoadPattern *thePattern)
  :Parameter(passedTag,1976), myPattern(thePattern), currentValue(0.0)
{

}


LoadFactorParameter::~LoadFactorParameter()
{

}

int
LoadFactorParameter::update(int newValue)
{
  currentValue = newValue;

  return 0;
}

int
LoadFactorParameter::update(double newValue)
{
  // ignore passed value, get straight from node
  // use setValue to assign a value to parameter directly
  //const Vector *u = myNode->getResponse(myType);
  //currentValue = (*u)(myDOF-1);
  
  currentValue = myPattern->getLoadFactor();

  return 0;
}

/*
int
LoadFactorParameter::activate(bool active)
{
  return 0;
}
*/

double
LoadFactorParameter::getValue(void)
{
  return currentValue;
}

void
LoadFactorParameter::setValue(double newValue)
{
  currentValue = newValue;
}

bool
LoadFactorParameter::isImplicit(void)
{
  return true;
}

double
LoadFactorParameter::getSensitivity(int index)
{
  double dudh = myPattern->getLoadFactorSensitivity(index);

  return dudh;
}

double
LoadFactorParameter::getPerturbation(void)
{
    return 0.005;
}

int 
LoadFactorParameter::getPointerTag(void) 
{
  return myPattern->getTag();
}

void
LoadFactorParameter::Print(OPS_Stream &s, int flag)  
{
  s << "LoadFactorParameter, tag = " << this->getTag() << endln;
  myPattern->Print(s, flag);
}

void
LoadFactorParameter::setDomain(Domain *theDomain)
{
  return;
}

int 
LoadFactorParameter::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int 
LoadFactorParameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}
