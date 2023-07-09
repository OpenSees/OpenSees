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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/NodeResponseParameter.cpp,v $

#include <classTags.h>
#include <NodeResponseParameter.h>
#include <Node.h>

NodeResponseParameter::NodeResponseParameter(int passedTag, Node *theNode, 
					     NodeResponseType theType,
					     int theDOF)
  :Parameter(passedTag,1976), myNode(theNode), myType(theType), myDOF(theDOF), currentValue(0.0)
{

}


NodeResponseParameter::~NodeResponseParameter()
{

}

int
NodeResponseParameter::update(int newValue)
{
  currentValue = newValue;

  return 0;
}

int
NodeResponseParameter::update(double newValue)
{
    // ignore passed value, get straight from node
    // use setValue to assign a value to parameter directly
    const Vector *u = myNode->getResponse(myType);
    
    currentValue = (*u)(myDOF-1);

    return 0;
}

/*
int
NodeResponseParameter::activate(bool active)
{
  return 0;
}
*/

double
NodeResponseParameter::getValue(void)
{
  return currentValue;
}

void
NodeResponseParameter::setValue(double newValue)
{
  currentValue = newValue;
}

bool
NodeResponseParameter::isImplicit(void)
{
  return true;
}

double
NodeResponseParameter::getSensitivity(int index)
{
  double dudh = myNode->getDispSensitivity(myDOF, index);
  //opserr << "NRP " << this->getTag() << " index = " << index << " dudh = " << dudh << endln;
  return dudh;
}

double
NodeResponseParameter::getPerturbation(void)
{
    return 0.005;
}

int 
NodeResponseParameter::getPointerTag(void) 
{
    return myNode->getTag();
}

void
NodeResponseParameter::Print(OPS_Stream &s, int flag)  
{
  s << "NodeResponseParameter, tag = " << this->getTag() << endln;
  myNode->Print(s, flag);
}

void
NodeResponseParameter::setDomain(Domain *theDomain)
{
  return;
}

int 
NodeResponseParameter::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int 
NodeResponseParameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return 0;
}
